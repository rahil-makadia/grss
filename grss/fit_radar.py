import numpy as np
from astropy.time import Time
from json import loads
from requests import request

from .fit_utilities import *

__all__ = [ 'get_radar_obs_array',
]

def get_radar_raw_data(tdes):
    """Get json data from JPL small-body radar API entry for desired small body

    Args:
        tdes (str/int): IMPORTANT: must be the designation of the small body, not the name.

    Returns:
        raw_data: JSON output of small body information query from JPL small-body radar API
    """
    url = f"https://ssd-api.jpl.nasa.gov/sb_radar.api?des={tdes}"
    
    r = request("GET",url)

    return loads(r.text)

def get_radar_obs_array(tdes):
    """Get radar observations from JPL small-body radar API entry for desired small body

    Args:
        tdes (str/int): IMPORTANT: must be the designation of the small body, not the name.

    Returns:
        obs_array_radar (np.ndarray): Radar observation data for the given body
        observer_codes_radar (tuple): Observer locations for each observation in obs_array_radar
    """
    # map to MPC code for radar observatories if it exists (otherwise use JPL code)
    radar_observer_map = {  '-1': '251', # Arecibo (300-m, 1963 to 2020)
                            '-2': '254', # Haystack (37-m)
                            '-9': '256', # Green Bank Telescope (100-m, GBT)
                            '-13': '252', # DSS-13 (34-m BWG, R&D)
                            '-14': '253', # DSS-14 (70-m)
                            '-25': '257', # Goldstone DSS 25
                            '-35': '-35', # DSS-35 (34-m BWG)
                            '-36': '-36', # DSS-36 (34-m BWG)
                            '-38': '255', # Evpatoria (70-m)
                            '-43': '-43', # DSS-43 (70-m)
                            '-47': '-47', # DSS-47 (ATCA ref. W196)
                            '-73': '259',} # Tromso (32-m, EISCAT)
    raw_data = get_radar_raw_data(tdes)
    if isinstance(raw_data, str):
        return None, None
    num_obs = int(raw_data['count'])
    data = raw_data['data']
    obs_array_radar = np.zeros((num_obs, 5))
    observer_codes_radar = []
    rows_to_delete = []
    for i in range(num_obs):
        obs = data[num_obs-i-1]
        date = Time(obs[1], format='iso', scale='utc')
        obs_val = float(obs[2])
        obs_sigma = float(obs[3])
        delay = obs[4] == 'us'
        doppler = obs[4] == 'Hz'
        freq = float(obs[5])*1e6 # MHz -> Hz
        tx = obs[6]
        rx = obs[7]

        obs_array_radar[i,0] = date.mjd
        if delay:
            obs_array_radar[i,1] = obs_val
            obs_array_radar[i,2] = np.nan
            obs_array_radar[i,3] = obs_sigma
            obs_array_radar[i,4] = np.nan
            observer_codes_radar.append((radar_observer_map[tx], radar_observer_map[rx]))
        elif doppler:
            # skip doppler observations for now
            rows_to_delete.append(i)
            continue
            obs_array_radar[i,1] = np.nan
            obs_array_radar[i,2] = obs_val
            obs_array_radar[i,3] = np.nan
            obs_array_radar[i,4] = obs_sigma
            observer_codes_radar.append(((radar_observer_map[rx], radar_observer_map[tx]), freq))
        else:
            raise ValueError("Observation type not recognized")
    obs_array_radar = np.delete(obs_array_radar, rows_to_delete, axis=0)
    return obs_array_radar, tuple(observer_codes_radar)
