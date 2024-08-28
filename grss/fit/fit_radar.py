"""Radar observation handling for the GRSS orbit determination code"""
from json import loads
import requests
from astropy.time import Time
import numpy as np

__all__ = [ 'add_radar_obs',
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
    req = requests.request("GET", url, timeout=30)
    return loads(req.text)

def add_radar_obs(obs_df, t_min_tdb=None, t_max_tdb=None, verbose=False):
    # sourcery skip: low-code-quality
    """
    Get radar observations from JPL small-body radar API entry for a desired small body
    and assemble the radar observations for the given body into an array for an orbit fit.

    Parameters
    ----------
    obs_df : pandas DataFrame
        Optical observation data for the given body
    t_min_tdb : float, optional
        Minimum time (MJD TDB) for observations to be included, by default None
    t_max_tdb : float, optional
        Maximum time (MJD TDB) for observations to be included, by default None
    verbose : bool, optional
        Flag to print out information about the observations, by default False

    Returns
    -------
    obs_df : pandas DataFrame
        Radar+Optical observation data for the given body

    Raises
    ------
    ValueError
        If the observation type is not recognized
    """
    perm_id = obs_df.iloc[-1]['permID']
    prov_id = obs_df.iloc[-1]['provID']
    if t_min_tdb is None:
        t_min_utc = -np.inf
    else:
        t_min_utc = Time(t_min_tdb, format='mjd', scale='tdb').utc.mjd
    if t_max_tdb is None:
        t_max_utc = np.inf
    else:
        t_max_utc = Time(t_max_tdb, format='mjd', scale='tdb').utc.mjd
    body_id = perm_id if isinstance(perm_id, str) else prov_id
    raw_data = get_radar_raw_data(body_id)
    if isinstance(raw_data, str) or ('code' in raw_data and raw_data['code'] == '400'):
        return obs_df
    num_obs = int(raw_data['count'])
    if verbose:
        print(f"Read in {num_obs} radar observations from JPL radar API.")
    data = raw_data['data']
    time_range_count = 0
    for i in range(num_obs):
        obs = data[num_obs-i-1]
        date = Time(obs[1], format='iso', scale='utc')
        if date.utc.mjd < t_min_utc or date.utc.mjd > t_max_utc:
            time_range_count += 1
            continue
        obs_val = float(obs[2])
        obs_sigma = float(obs[3])
        delay = obs[4] == 'us'
        doppler = obs[4] == 'Hz'
        freq = float(obs[5])
        # transmitter and receiver codes, use radar_observer_map if you
        # want to use MPC station info (less accurate in my experience)
        rx_code = obs[6]
        tx_code = obs[7]
        bounce_point_int = 1 if obs[8] == 'C' else 0
        idx = len(obs_df)
        obs_df.loc[idx,'permID'] = perm_id
        obs_df.loc[idx,'provID'] = prov_id
        obs_df.loc[idx,'obsTime'] = f'{date.utc.isot}Z'
        obs_df.loc[idx,'obsTimeMJD'] = date.utc.mjd
        obs_df.loc[idx,'obsTimeMJDTDB'] = date.tdb.mjd
        obs_df.loc[idx,'mode'] = 'RAD'
        obs_df.loc[idx,'rcv'] = rx_code
        obs_df.loc[idx,'trx'] = tx_code
        obs_df.loc[idx,'delay'] = obs_val/1.0e6 if delay else np.nan
        obs_df.loc[idx,'rmsDelay'] = obs_sigma if delay else np.nan
        obs_df.loc[idx,'doppler'] = obs_val if doppler else np.nan
        obs_df.loc[idx,'rmsDoppler'] = obs_sigma if doppler else np.nan
        obs_df.loc[idx,'com'] = bounce_point_int
        obs_df.loc[idx,'frq'] = freq
        obs_df.loc[idx,'sigDelay'] = obs_sigma if delay else np.nan
        obs_df.loc[idx,'sigDoppler'] = obs_sigma if doppler else np.nan
    if verbose:
        print(f"\tFiltered to {num_obs-time_range_count} observations that satisfy the "
                "time range constraints.")
    return obs_df
