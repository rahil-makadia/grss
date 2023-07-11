"""Radar observation handling for the GRSS orbit determination code"""
from json import loads
from requests import request
from astropy.time import Time
import numpy as np

__all__ = [ 'get_radar_obs_array',
]

def get_radar_raw_data(tdes):
    """
    Get JSON data from JPL small-body radar API entry for desired small body

    Parameters
    ----------
    tdes : str/int
        IMPORTANT: must be the designation of the small body, not the name.

    Returns
    -------
    raw_data : dict
        JSON output of small body information query from JPL small-body radar API
    """
    url = f"https://ssd-api.jpl.nasa.gov/sb_radar.api?des={tdes}"
    req = request("GET", url, timeout=30)
    return loads(req.text)

def get_radar_obs_array(tdes, t_min_tdb=None, t_max_tdb=None, verbose=False):
    # sourcery skip: low-code-quality
    """
    Get radar observations from JPL small-body radar API entry for a desired small body
    and assemble the radar observations for the given body into an array for an orbit fit.

    Parameters
    ----------
    tdes : str/int
        IMPORTANT: must be the designation of the small body, not the name.
    t_min_tdb : float, optional
        Minimum time (MJD TDB) for observations to be included, by default None
    t_max_tdb : float, optional
        Maximum time (MJD TDB) for observations to be included, by default None
    verbose : bool, optional
        Flag to print out information about the observations, by default False

    Returns
    -------
    obs_array_radar : array
        Radar observation data for the given body
    observer_codes_radar : tuple
        Observer locations for each observation in obs_array_radar

    Raises
    ------
    ValueError
        If the observation type is not recognized
    """
    if t_min_tdb is None:
        t_min_utc = -np.inf
    else:
        t_min_utc = Time(t_min_tdb, format='mjd', scale='tdb').utc.mjd
    if t_max_tdb is None:
        t_max_utc = np.inf
    else:
        t_max_utc = Time(t_max_tdb, format='mjd', scale='tdb').utc.mjd
    # map to MPC code for radar observatories if it exists (otherwise use JPL code)
    use_mpc_mapping = False
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
    obs_array_radar = np.zeros((num_obs, 6))
    observer_codes_radar = []
    rows_to_delete = []
    for i in range(num_obs):
        obs = data[num_obs-i-1]
        date = Time(obs[1], format='iso', scale='utc')
        if date.utc.mjd < t_min_utc or date.utc.mjd > t_max_utc:
            rows_to_delete.append(i)
            continue
        obs_val = float(obs[2])
        obs_sigma = float(obs[3])
        delay = obs[4] == 'us'
        doppler = obs[4] == 'Hz'
        freq = float(obs[5])*1e6 # MHz -> Hz
        # transmitter and receiver codes, use radar_observer_map if you
        # want to use MPC station info (less accurate in my experience)
        tx_code = radar_observer_map[obs[6]] if use_mpc_mapping else obs[6]
        rx_code = radar_observer_map[obs[7]] if use_mpc_mapping else obs[7]
        bounce_point = obs[8]
        bounce_point_int = 0 if bounce_point == 'C' else 1

        obs_array_radar[i,0] = date.utc.mjd
        obs_array_radar[i,5] = np.nan
        if delay:
            obs_array_radar[i,1] = obs_val
            obs_array_radar[i,2] = np.nan
            obs_array_radar[i,3] = obs_sigma
            obs_array_radar[i,4] = np.nan
            observer_codes_radar.append(((tx_code, rx_code), bounce_point_int))
        elif doppler:
            obs_array_radar[i,1] = np.nan
            obs_array_radar[i,2] = obs_val
            obs_array_radar[i,3] = np.nan
            obs_array_radar[i,4] = obs_sigma
            observer_codes_radar.append(((tx_code, rx_code), bounce_point_int, freq))
        else:
            raise ValueError("Observation type not recognized")
    if verbose:
        print(f"Deleted {len(rows_to_delete)} radar observations outside of time range")
    obs_array_radar = np.delete(obs_array_radar, rows_to_delete, axis=0)
    return obs_array_radar, tuple(observer_codes_radar)
