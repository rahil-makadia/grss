import os
file_path = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.append(f'{file_path}/../extern/')
import astrocat_debiasing.debias as ad

from astropy.time import Time
from astroquery.mpc import MPC
import numpy as np
import pandas as pd

from .fit_utilities import *

__all__ = [ 'get_optical_obs_array',
]

def get_optical_data(body_id, optical_obs_file, t_min_tdb=None, t_max_tdb=None, verbose=False):
    """Get observation data for a given body.

    Args:
        body_id (int or str): Target id, numbers are interpreted as asteroids, append 'P' for comets, start with comet type and a '/' for comet designations
        t_min_tdb (float or None, optional): Minimum date to cap the optical observation epochs. Defaults to None.
        t_max_tdb (float or None, optional): Maximum date to cap the optical observation epochs. Defaults to None.
        verbose (bool, optional): Toggle for whether to print out information about the observations. Defaults to False.

    Returns:
        obs_array_optical (np.ndarray): Optical observation data for the given body
        star_catalog_codes (tuple): Star star_catalog_codes used for reducing the observations in obs_array_optical
        observer_codes_optical (tuple): Observer locations for each observation in obs_array_optical
    """
    if t_min_tdb is None:
        t_min_utc = -np.inf
    else:
        t_min_utc = Time(t_min_tdb, format='mjd', scale='tdb').utc.mjd
    if t_max_tdb is None:
        t_max_utc = np.inf
    else:
        t_max_utc = Time(t_max_tdb, format='mjd', scale='tdb').utc.mjd
    if optical_obs_file is None:
        obs_raw = MPC.get_observations(body_id, get_mpcformat=True, cache=False)['obs']
    else:
        with open(optical_obs_file, 'r') as f:
            obs_raw = f.readlines()
            obs_raw = [x[:-1] for x in obs_raw] # remove \n at the end of each line
    obs_array_optical = np.zeros((len(obs_raw), 5))
    star_catalog_codes = []
    observer_codes_optical = []
    observation_type_codes = []
    observer_program_codes = []
    space_observatory_codes = ['245', '249', '250', '258', '274', 'C49', 'C50', 'C51', 'C52', 'C53', 'C54', 'C55', 'C56', 'C57', 'C59',]
    non_geocentric_occultation_codes = ['275']
    unsupported_code_counter = 0
    unsupported_type_counter = 0
    out_of_range_counter = 0
    skip_counter = 0
    for i in range(len(obs_raw)):
        data = obs_raw[i]
        obs_type = 'P' if data[14] == ' ' else data[14]
        date = data[15:32]
        obs_time_mjd = Time(date[:-7].replace(' ','-'), format='iso', scale='utc').utc.mjd+float(date[-7:])
        obs_code = data[77:80]
        disallow_code = obs_code in space_observatory_codes + non_geocentric_occultation_codes
        disallow_type = obs_type in ['R', 'r', 'V', 'v']
        disallow_time = obs_time_mjd < t_min_utc or obs_time_mjd > t_max_utc
        if not disallow_code and not disallow_type and not disallow_time:
            observer_program = data[13]
            ra = get_ra_from_hms(data[32:44])
            dec = get_dec_from_dms(data[44:56])
            obs_array_optical[i, 0] = obs_time_mjd
            obs_array_optical[i, 1] = ra
            obs_array_optical[i, 2] = dec
            star_catalog = data[71]
            star_catalog_codes.append(star_catalog)
            observer_codes_optical.append(obs_code)
            observation_type_codes.append(obs_type)
            observer_program_codes.append(observer_program)
        else:
            if disallow_code:
                increment = 1
                unsupported_code_counter += increment
            elif disallow_type:
                increment = 0.5
                unsupported_type_counter += increment
            else: # disallow_time
                increment = 1
                out_of_range_counter += increment
            skip_counter += increment
            obs_array_optical[i, :] = np.nan
            star_catalog_codes.append(np.nan)
            observer_codes_optical.append(np.nan)
            observation_type_codes.append(np.nan)
            observer_program_codes.append(np.nan)
    if verbose:
        print(f"Skipped {int(skip_counter)} observations \n\t {unsupported_code_counter} of which were non-geocentric occultations or space-based observations, \n\t {int(unsupported_type_counter)} were either roving or radar observations (radar is handled separately), \n\t {out_of_range_counter} of which were outside the specified time range.")
    star_catalog_codes = tuple(np.array(star_catalog_codes)[~np.isnan(obs_array_optical[:, 0])])
    observer_codes_optical = tuple(np.array(observer_codes_optical)[~np.isnan(obs_array_optical[:, 0])])
    observation_type_codes = tuple(np.array(observation_type_codes)[~np.isnan(obs_array_optical[:, 0])])
    observer_program_codes = tuple(np.array(observer_program_codes)[~np.isnan(obs_array_optical[:, 0])])
    obs_array_optical = obs_array_optical[~np.isnan(obs_array_optical[:, 0])]
    return obs_array_optical, star_catalog_codes, observer_codes_optical, observation_type_codes, observer_program_codes

def apply_debiasing_scheme(obs_array_optical, star_catalog_codes, observer_codes_optical, lowres=False, verbose=False):
    """Apply the Eggl et al. (2020) debiasing scheme to the observation array from https://doi.org/10.1016/j.icarus.2019.113596.

    Args:
        obs_array_optical (np.ndarray): Optical observation data for the given body
        star_catalog_codes (tuple): Star star_catalog_codes used for reducing the observations in obs_array_optical
        observer_codes_optical (tuple): Observer locations for each observation in obs_array_optical
        lowres (bool, optional): Toggle for whether to use the low resolution debiasing scheme. Defaults to False.
        verbose (bool, optional): Toggle for whether to print information about the debiasing scheme. Defaults to False.

    Returns:
        obs_array_optical (np.ndarray): Optical observation data for the given body with debiasing applied
        star_catalog_codes (tuple): Star star_catalog_codes used for reducing the observations in obs_array_optical with debiasing applied
        observer_codes_optical (tuple): Observer locations for each observation in obs_array_optical with debiasing applied
    """
    star_catalog_codes = list(star_catalog_codes)
    observer_codes_optical = list(observer_codes_optical)
    biased_catalogs = ['a', 'b', 'c', 'd', 'e', 'g', 'i', 'j', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 't', 'u', 'v', 'w', 'L', 'N', 'Q', 'R', 'S', 'U', 'W'] # MPC star_catalog_codes codes, https://www.minorplanetcenter.net/iau/info/CatalogueCodes.html
    unbiased_catalogs = [   'g', # Tycho-2
                            'l', # ACT
                            'U', # Gaia DR1
                            'V', # Gaia DR2
                            'W', # Gaia DR3
                            'X', # Gaia EDR3
                            'Y', # UCAC-5
                        ]
    columns = []
    for cat in biased_catalogs:
        columns.extend([f"{cat}_ra", f"{cat}_dec", f"{cat}_pm_ra", f'{cat}_pm_dec'])
    if lowres:
        biasfile = 'bias_lowres.dat'
        nside = 64
    else:
        biasfile = 'bias.dat'
        nside = 256
    biasdf = pd.read_csv(biasfile,sep='\s+',skiprows=23,names=columns)
    unbias_counter = 0
    debias_counter = 0
    for i, row in enumerate(obs_array_optical):
        obs_time_jd = Time(row[0], format='mjd', scale='utc').tt.jd
        ra = row[1]/3600*np.pi/180
        dec = row[2]/3600*np.pi/180
        star_catalog = star_catalog_codes[i]
        if star_catalog in unbiased_catalogs:
            unbias_counter += 1
            continue
        elif star_catalog in biased_catalogs:
            debias_counter += 1
            # apply debiasing from Eggl et al. 2020, https://doi.org/10.1016/j.icarus.2019.113596
            ra_new, dec_new = ad.debiasRADec(ra, dec, obs_time_jd, star_catalog, biasdf, nside=nside)
            obs_array_optical[i, 1] = ra_new*180/np.pi*3600
            obs_array_optical[i, 2] = dec_new*180/np.pi*3600
            margin = 1
            if (verbose and (np.rad2deg(ra_new - ra) > margin or np.rad2deg(dec_new - dec) > margin)):
                print(f"Debiased observation {i} from observatory {observer_codes_optical[i]} using star_catalog_codes {star_catalog}. RA bias = {np.rad2deg(ra_new-ra)*3600:0.4f} arcsec, DEC bias = {np.rad2deg(dec_new-dec)*3600:0.4f} arcsec")
        # else:
            # skip_counter += 1
            # skip_catalogs.append(star_catalog)
            # obs_array_optical[i, :] = np.nan
            # star_catalog_codes[i] = np.nan
            # observer_codes_optical[i] = np.nan
            # if verbose:
            #     print(f"Skipping observation {i} because it is not in a biased MPC catalog (catalog '{star_catalog}') and debiasing is turned on.")
    if verbose:
        print(f"No debiasing needed for {unbias_counter} observations. Debiased {debias_counter} observations. No biasing information for {len(obs_array_optical)-unbias_counter-debias_counter} observations.")
    star_catalog_codes = tuple(np.array(star_catalog_codes)[~np.isnan(obs_array_optical[:, 0])])
    observer_codes_optical = tuple(np.array(observer_codes_optical)[~np.isnan(obs_array_optical[:, 0])])
    obs_array_optical = obs_array_optical[~np.isnan(obs_array_optical[:, 0])]
    return obs_array_optical, star_catalog_codes, observer_codes_optical

def apply_weights_old(obs_array_optical, star_catalog_codes, observer_codes_optical, observation_type_codes, observer_program_codes, verbose=False):
    """Apply the weights from Vereš et al. (2017) to the observation array from https://doi.org/10.1016/j.icarus.2017.05.021.

    Args:
        obs_array_optical (np.ndarray): Optical observation data for the given body
        star_catalog_codes (tuple): Star star_catalog_codes used for reducing the observations in obs_array_optical
        observer_codes_optical (tuple): Observer locations for each observation in obs_array_optical
        verbose (bool, optional): Toggle for whether to print information about the applied weights. Defaults to False.

    Returns:
        obs_array_optical (np.ndarray): Optical observation data for the given body with weights applied
    """
    # weighting scheme from vereš et al. 2017, tables 2-4, http://dx.doi.org/10.1016/j.icarus.2017.05.021
    default_weight_counter = 0
    for i in range(len(obs_array_optical)):
        star_catalog = star_catalog_codes[i]
        obs_code = observer_codes_optical[i]
        obs_date_mjd = obs_array_optical[i, 0]
        # table 2
        if obs_code == '703':
            obs_array_optical[i, 3:5] = 0.8 if obs_date_mjd >= Time('2014-01-01', format='iso', scale='utc').mjd else 1.0
        elif obs_code == '691':
            obs_array_optical[i, 3:5] = 0.5 if obs_date_mjd >= Time('2003-01-01', format='iso', scale='utc').mjd else 0.6
        elif obs_code == '644':
            obs_array_optical[i, 3:5] = 0.4 if obs_date_mjd >= Time('2003-09-01', format='iso', scale='utc').mjd else 0.6
        # table 3
        elif obs_code in ['704', 'C51', 'J75']:
            obs_array_optical[i, 3:5] = 1.0
        elif obs_code == 'G96':
            obs_array_optical[i, 3:5] = 0.5
        elif obs_code == 'F51':
            obs_array_optical[i, 3:5] = 0.2
        elif obs_code in ['G45', '608']:
            obs_array_optical[i, 3:5] = 0.6
        elif obs_code == '699':
            obs_array_optical[i, 3:5] = 0.8
        elif obs_code in ['D29', 'E12']:
            obs_array_optical[i, 3:5] = 0.75
        # table 4
        elif obs_code in ['645', '673', 'H01']:
            obs_array_optical[i, 3:5] = 0.3
        elif obs_code in ['689', '950']:
            obs_array_optical[i, 3:5] = 0.5
        elif obs_code in ['J04', 'K92', 'K93', 'Q63', 'Q64', 'V37', 'W84', 'W85', 'W86', 'W87', 'K91', 'E10', 'F65']:
            obs_array_optical[i, 3:5] = 0.4
        elif obs_code == 'G83' and star_catalog in ['U', 'V', 'W', 'X']:
            obs_array_optical[i, 3:5] = 0.2
        elif obs_code == 'G83' and star_catalog in ['q', 't']:
            obs_array_optical[i, 3:5] = 0.3
        elif obs_code == 'Y28' and star_catalog in ['t', 'U', 'V', 'W', 'X']:
            obs_array_optical[i, 3:5] = 0.3
        elif obs_code == '568' and star_catalog in ['o', 's']:
            obs_array_optical[i, 3:5] = 0.5
        elif obs_code == '568' and star_catalog in ['U', 'V', 'W', 'X']:
            obs_array_optical[i, 3:5] = 0.1
        elif obs_code == '568' and star_catalog == 't':
            obs_array_optical[i, 3:5] = 0.2
        elif obs_code in ['T09', 'T12', 'T14'] and star_catalog in ['U', 'V', 'W', 'X']:
            obs_array_optical[i, 3:5] = 0.1
        elif obs_code == '309' and star_catalog == ['q', 't']:
            obs_array_optical[i, 3:5] = 0.3
        elif obs_code == '309' and star_catalog in ['U', 'V', 'W', 'X']:
            obs_array_optical[i, 3:5] = 0.2
        else:
            obs_array_optical[i, 3:5] = 1.0
            default_weight_counter += 1
    if verbose:
        print(f"Applied default weight of 1 arcsec to {default_weight_counter} observations")
        # find where the weights are 0
        zero_weight_indices = np.where(obs_array_optical[:, 3] == 0)[0]
        if len(zero_weight_indices) > 0:
            raise ValueError(f"apply_weights: Found {len(zero_weight_indices)} observations with zero weight. Check the star_catalog_codes and observer location codes and make sure proper weights are applied.")
    return obs_array_optical

def apply_weights(obs_array_optical, star_catalog_codes, observer_codes_optical, observation_type_codes, observer_program_codes, verbose=False):
    """Apply the weights from Vereš et al. (2017) to the observation array from https://doi.org/10.1016/j.icarus.2017.05.021.

    Args:
        obs_array_optical (np.ndarray): Optical observation data for the given body
        star_catalog_codes (tuple): Star star_catalog_codes used for reducing the observations in obs_array_optical
        observer_codes_optical (tuple): Observer locations for each observation in obs_array_optical
        verbose (bool, optional): Toggle for whether to print information about the applied weights. Defaults to False.

    Returns:
        obs_array_optical (np.ndarray): Optical observation data for the given body with weights applied
    """
    # weighting scheme from vereš et al. 2017, obswgt_veres2017.inp
    default_weight_counter = 0
    for i in range(len(obs_array_optical)):
        star_catalog = star_catalog_codes[i]
        obs_code = observer_codes_optical[i]
        obs_date_mjd = obs_array_optical[i, 0]
        obs_type = observation_type_codes[i]
        obs_program = observer_program_codes[i]
        # obs_array_optical[i, 3:5] = 1.0
        # continue
        if obs_type in ['P', 'A', 'N', ' ']:
            if obs_date_mjd <= Time('1890-01-01', format='iso', scale='utc').utc.mjd:
                obs_array_optical[i, 3:5] = 10.0
            elif obs_date_mjd > Time('1890-01-01', format='iso', scale='utc').utc.mjd and obs_date_mjd <= Time('1950-01-01', format='iso', scale='utc').utc.mjd:
                obs_array_optical[i, 3:5] = 5.0
            elif obs_date_mjd > Time('1950-01-01', format='iso', scale='utc').utc.mjd:
                obs_array_optical[i, 3:5] = 2.5
        elif obs_type in ['E', 'H']: # occultation or hipparcos
            obs_array_optical[i, 3:5] = 0.2
        elif obs_type  in ['T']: # transit circle
            obs_array_optical[i, 3:5] = 0.5
        elif obs_type in ['e']: # encoder
            obs_array_optical[i, 3:5] = 0.75
        elif obs_type in ['M']: # micrometer
            obs_array_optical[i, 3:5] = 2.0
        elif obs_type in ['S']: # satellite
            if obs_code in '275': # non-geocentric occultation
                obs_array_optical[i, 3:5] = 0.2
            else:
                obs_array_optical[i, 3:5] = 1.5
        elif obs_type in ['c', 'C', 'V', 'n', 'B']:
            if obs_type in ['c', 'C']:
                if obs_code in ['F51', 'F52']:
                    obs_array_optical[i, 3:5] = 0.2
                    if obs_program == 'Z':
                        obs_array_optical[i, 3:5] = 1.0
                elif obs_code in ['G96']:
                    obs_array_optical[i, 3:5] = 0.5
                    if obs_program == 'Z':
                        obs_array_optical[i, 3:5] = 1.0
                elif obs_code in ['703']:
                    if obs_date_mjd < Time('2014-01-01', format='iso', scale='utc').utc.mjd:
                        obs_array_optical[i, 3:5] = 1.0
                    else:
                        obs_array_optical[i, 3:5] = 0.8
                    if obs_program == 'Z':
                        obs_array_optical[i, 3:5] = 1.0
                elif obs_code in ['E12']:
                    obs_array_optical[i, 3:5] = 0.75
                elif obs_code in ['704']:
                    obs_array_optical[i, 3:5] = 1.0
                elif obs_code in ['691', '291']: # 291 as well?
                    if obs_date_mjd < Time('2003-01-01', format='iso', scale='utc').utc.mjd:
                        obs_array_optical[i, 3:5] = 0.6
                    else:
                        obs_array_optical[i, 3:5] = 0.5
                elif obs_code in ['644']:
                    if obs_date_mjd < Time('2003-09-01', format='iso', scale='utc').utc.mjd:
                        obs_array_optical[i, 3:5] = 0.6
                    else:
                        obs_array_optical[i, 3:5] = 0.4
                elif obs_code in ['699']:
                    obs_array_optical[i, 3:5] = 0.8
                elif obs_code in ['G45']:
                    obs_array_optical[i, 3:5] = 0.6
                elif obs_code in ['D29']:
                    obs_array_optical[i, 3:5] = 0.75
                elif obs_code in ['T07', 'M22', 'W68']:
                    obs_array_optical[i, 3:5] = 0.8
                elif obs_code in ['T05', 'T08']:
                    obs_array_optical[i, 3:5] = 1.0
                elif obs_code in ['568']:
                    obs_array_optical[i, 3:5] = 0.5
                    if obs_program == '_' and star_catalog in ['t', 'q']:
                        obs_array_optical[i, 3:5] = 0.25
                    elif obs_program == '_' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.15
                    elif obs_program == '^' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.2
                    elif obs_program == '2' and star_catalog in ['o', 's']:
                        obs_array_optical[i, 3:5] = 0.5
                    elif obs_program == '2' and star_catalog in ['t']:
                        obs_array_optical[i, 3:5] = 0.2
                    elif obs_program == '2' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.1
                    elif obs_program == '0':
                        obs_array_optical[i, 3:5] = 1.0
                elif obs_code in ['T09']:
                    obs_array_optical[i, 3:5] = 0.5
                    if obs_program == '0':
                        obs_array_optical[i, 3:5] = 0.1
                    elif obs_program == '1':
                        obs_array_optical[i, 3:5] = 1.0
                    elif obs_program == '2':
                        obs_array_optical[i, 3:5] = 0.5
                    elif obs_program == '3':
                        obs_array_optical[i, 3:5] = 0.5
                    if obs_date_mjd < Time('2023-10-06', format='iso', scale='utc').utc.mjd:
                        obs_array_optical[i, 3:5] = 0.1
                elif obs_code in ['T10', 'T11', 'T13', 'T16', 'T17']:
                    obs_array_optical[i, 3:5] = 0.5
                elif obs_code in ['T12']:
                    obs_array_optical[i, 3:5] = 0.5
                    if obs_program == '0':
                        obs_array_optical[i, 3:5] = 0.1
                    elif obs_program == '1':
                        obs_array_optical[i, 3:5] = 0.2
                    elif obs_program == '2':
                        obs_array_optical[i, 3:5] = 0.5
                    elif obs_program == '3':
                        obs_array_optical[i, 3:5] = 0.5
                    if obs_date_mjd < Time('2022-10-06', format='iso', scale='utc').utc.mjd:
                        obs_array_optical[i, 3:5] = 0.1
                elif obs_code in ['T14']:
                    obs_array_optical[i, 3:5] = 0.5
                    if obs_program == '0':
                        obs_array_optical[i, 3:5] = 0.1
                    elif obs_program == '1':
                        obs_array_optical[i, 3:5] = 0.5
                    elif obs_program == '2':
                        obs_array_optical[i, 3:5] = 0.5
                    elif obs_program == '3':
                        obs_array_optical[i, 3:5] = 0.2
                    elif obs_program == '4':
                        obs_array_optical[i, 3:5] = 0.5
                    elif obs_program == '5':
                        obs_array_optical[i, 3:5] = 1.0
                    elif obs_program == '6':
                        obs_array_optical[i, 3:5] = 0.5
                    elif obs_program == '7':
                        obs_array_optical[i, 3:5] = 0.1
                    if obs_date_mjd < Time('2022-10-06', format='iso', scale='utc').utc.mjd:
                        obs_array_optical[i, 3:5] = 0.1
                elif obs_code in ['T15']:
                    obs_array_optical[i, 3:5] = 0.5
                    if obs_program  == '2':
                        obs_array_optical[i, 3:5] = 1.0
                elif obs_code in ['H01']:
                    obs_array_optical[i, 3:5] = 1.0
                    if star_catalog in ['L', 't', 'U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.3
                elif obs_code in ['673']:
                    obs_array_optical[i, 3:5] = 1.0
                    if obs_program == '1':
                        obs_array_optical[i, 3:5] = 0.3
                elif obs_code in ['645']:
                    obs_array_optical[i, 3:5] = 0.3
                elif obs_code in ['689']:
                    obs_array_optical[i, 3:5] = 0.5
                elif obs_code in ['J04', 'Z84']:
                    obs_array_optical[i, 3:5] = 1.0
                    if star_catalog in ['U', 'V', 'W', 'X', 't', 'L', 'q', 'r', 'u', 'e']:
                        obs_array_optical[i, 3:5] = 0.5
                elif obs_code in ['608']:
                    obs_array_optical[i, 3:5] = 1.0
                    obs_array_optical[i, 3:5] = 0.6
                elif obs_code in ['W84']:
                    obs_array_optical[i, 3:5] = 0.5
                    if obs_program == '&':
                        obs_array_optical[i, 3:5] = 1.0
                elif obs_code in ['950']:
                    obs_array_optical[i, 3:5] = 1.0
                    if star_catalog in ['U', 'V', 'W', 'X', 't', 'L', 'q', 'r', 'u', 'e']:
                        obs_array_optical[i, 3:5] = 0.5
                    if obs_program == '@':
                        obs_array_optical[i, 3:5] = 1.0
                elif obs_code in ['E10']:
                    obs_array_optical[i, 3:5] = 0.4
                    if obs_program == '(' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.2
                elif obs_code in ['F65']:
                    obs_array_optical[i, 3:5] = 0.4
                    if obs_program == '8' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.2
                    elif obs_program == '9' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.2
                elif obs_code in ['K91', 'K92', 'K93', 'Q63', 'Q64', 'V37', 'V39', 'W85', 'W86', 'W87', 'Z24', 'Z31']:
                    obs_array_optical[i, 3:5] = 0.4
                    if obs_program == '0':
                        obs_array_optical[i, 3:5] = 0.3
                elif obs_code in ['Y28']:
                    obs_array_optical[i, 3:5] = 1.0
                    if star_catalog in ['t', 'U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.3
                elif obs_code in ['309']:
                    obs_array_optical[i, 3:5] = 1.0
                    if star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.2
                    elif star_catalog in ['t', 'q']:
                        obs_array_optical[i, 3:5] = 0.3
                    if obs_program in ['&', '%']:
                        if star_catalog in ['U', 'V', 'W', 'X']:
                            obs_array_optical[i, 3:5] = 0.1
                        elif star_catalog in ['t', 'q']:
                            obs_array_optical[i, 3:5] = 0.3
                elif obs_code in ['G83']:
                    obs_array_optical[i, 3:5] = 1.0
                    if obs_program == '2':
                        if star_catalog in ['U', 'V', 'W', 'X']:
                            obs_array_optical[i, 3:5] = 0.2
                        elif star_catalog in ['t', 'q']:
                            obs_array_optical[i, 3:5] = 0.3
                elif obs_code in ['C65']:
                    obs_array_optical[i, 3:5] = 1.0
                    if obs_program == '3' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.2
                elif obs_code in ['094']:
                    obs_array_optical[i, 3:5] = 1.0
                    if obs_program == '9' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.5
                elif obs_code in ['Z18']:
                    obs_array_optical[i, 3:5] = 1.0
                    if obs_program == '1' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.2
                elif obs_code in ['181']:
                    obs_array_optical[i, 3:5] = 1.0
                    if obs_program == '1' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.5
                elif obs_code in ['G37']:
                    obs_array_optical[i, 3:5] = 1.0
                    if obs_program == '5' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.2
                elif obs_code in ['J04']:
                    obs_array_optical[i, 3:5] = 1.0
                    if obs_program == '2' and star_catalog in ['t', 'q', 'U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.3
                    elif obs_program == '#':
                        obs_array_optical[i, 3:5] = 1.0
                elif obs_code in ['Z84']:
                    obs_array_optical[i, 3:5] = 1.0
                    if obs_program == '1' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.3
                elif obs_code in ['D20']:
                    obs_array_optical[i, 3:5] = 1.0
                    if obs_program == '1' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.5
                elif obs_code in ['Z18']:
                    obs_array_optical[i, 3:5] = 1.0
                    if obs_program == '1' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.1
                elif obs_code in ['N50']:
                    obs_array_optical[i, 3:5] = 1.0
                    if obs_program == '2' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.3
                elif obs_code in ['I41']:
                    obs_array_optical[i, 3:5] = 1.0
                    if obs_program == 'Z':
                        obs_array_optical[i, 3:5] = 1.0
                else:
                    default_weight_counter += 1
                    obs_array_optical[i, 3:5] = 1.0
            if star_catalog == ' ':
                obs_array_optical[i, 3:5] = 1.5
            if obs_type in ['V']:
                obs_array_optical[i, 3:5] = 2.0
        else:
            raise ValueError(f"apply_weights: Observation type {obs_type} not recognized.")
    if verbose:
        print(f"Applied default weight of 1 arcsec to {default_weight_counter} observations")
        # find where the weights are 0
        zero_weight_indices = np.where(obs_array_optical[:, 3] == 0)[0]
        if len(zero_weight_indices) > 0:
            raise ValueError(f"apply_weights: Found {len(zero_weight_indices)} observations with zero weight. Check the star_catalog_codes and observer location codes and make sure proper weights are applied.")
    return obs_array_optical

def deweight_obs(obs_array_optical, star_catalog_codes, observer_codes_optical, eff_obs_per_night=5, verbose=False):
    """Deweight observations that occur on the same night at the same observatory.

    Args:
        obs_array_optical (np.ndarray): Optical observation data for the given body
        star_catalog_codes (tuple): Star star_catalog_codes used for processing the observations in obs_array_optical
        observer_codes_optical (tuple): Observer locations for each observation in obs_array_optical
        eff_obs_per_night (int, optional): Effective number of observations per night from the same observatory. Defaults to 5.
        verbose (bool, optional): Toggle for whether to print information about the deweighting scheme. Defaults to False.

    Returns:
        obs_array_deweighted: Deweighted Optical observation data for the given body
        catalog_deweighted: Star star_catalog_codes used for processing the observations in obs_array_deweighted
        observer_loc_deweighted: Observer locations for each observation in obs_array_deweighted
    """
    obs_array_optical = obs_array_optical[obs_array_optical[:,0].argsort()]
    obs_array_deweighted = obs_array_optical[0,None].copy()
    catalog_deweighted = [star_catalog_codes[0]]
    observer_loc_deweighted = [observer_codes_optical[0]]
    night_count = 1
    deweight_count = 0
    for i in range(1, len(obs_array_optical)):
        curr_jd_day = np.floor(obs_array_optical[i,0]+2400000.5)
        prev_jd_day = np.floor(obs_array_optical[i-1,0]+2400000.5)
        curr_observatory = observer_codes_optical[i]
        prev_observatory = observer_codes_optical[i-1]
        obs_array_deweighted = np.vstack((obs_array_deweighted, obs_array_optical[i]))
        catalog_deweighted.append(star_catalog_codes[i])
        observer_loc_deweighted.append(observer_codes_optical[i])
        if curr_jd_day == prev_jd_day and curr_observatory == prev_observatory:
            night_count += 1
            if night_count >= eff_obs_per_night:
                deweight_count += 1
                factor = night_count**0.5/eff_obs_per_night**0.5
                obs_array_deweighted[i-night_count+1:i+1,3:5] = factor * obs_array_optical[i-night_count+1:i+1,3:5]
        else:
            night_count = 1
    if verbose:
        print(f"Deweighted {deweight_count} observations as part of deweighting scheme.")
    return obs_array_deweighted, tuple(catalog_deweighted), tuple(observer_loc_deweighted)

def eliminate_obs(obs_array_optical, star_catalog_codes, observer_codes_optical, max_obs_per_night=5, verbose=False):
    """Limit the number of observations that occur on the same night at the same observatory.

    Args:
        obs_array_optical (np.ndarray): Optical observation data for the given body
        star_catalog_codes (tuple): Star star_catalog_codes used for processing the observations in obs_array_optical
        observer_codes_optical (tuple): Observer locations for each observation in obs_array_optical
        max_obs_per_night (int, optional): Maximum number of observations per night from the same observatory. Defaults to 5.
        verbose (bool, optional): Toggle for whether to print information about the deweighting scheme. Defaults to False.

    Returns:
        obs_array_eliminated: Eliminated Optical observation data for the given body
        catalog_eliminated: Star star_catalog_codes used for processing the observations in obs_array_eliminated
        observer_loc_eliminated: Observer locations for each observation in obs_array_eliminated
    """
    obs_array_optical = obs_array_optical[obs_array_optical[:,0].argsort()]
    obs_array_eliminated = obs_array_optical[0,None].copy()
    catalog_eliminated = [star_catalog_codes[0]]
    observer_loc_eliminated = [observer_codes_optical[0]]
    night_count = 1
    for i in range(1, len(obs_array_optical)):
        curr_jd_day = np.floor(obs_array_optical[i,0]+2400000.5)
        prev_jd_day = np.floor(obs_array_optical[i-1,0]+2400000.5)
        curr_observatory = observer_codes_optical[i]
        prev_observatory = observer_codes_optical[i-1]
        if curr_jd_day == prev_jd_day and curr_observatory == prev_observatory:
            night_count += 1
            if night_count <= max_obs_per_night:
                obs_array_eliminated = np.vstack((obs_array_eliminated, obs_array_optical[i]))
                catalog_eliminated.append(star_catalog_codes[i])
                observer_loc_eliminated.append(observer_codes_optical[i])
        else:
            obs_array_eliminated = np.vstack((obs_array_eliminated, obs_array_optical[i]))
            catalog_eliminated.append(star_catalog_codes[i])
            observer_loc_eliminated.append(observer_codes_optical[i])
            night_count = 1
    if verbose:
        print(f"Removed {len(obs_array_optical)-len(obs_array_eliminated)} observations as part of elimination scheme.")
    return obs_array_eliminated, tuple(catalog_eliminated), tuple(observer_loc_eliminated)

def get_optical_obs_array(body_id, optical_obs_file=None, t_min_tdb=None, t_max_tdb=None, debias=True, debias_lowres=False, deweight=True, eliminate=False, max_obs_per_night=5, verbose=False, old_weights=False):
    """Assemble the optical observations for a given body into an array for an orbit fit.

    Args:
        body_id (int or str): Target id, numbers are interpreted as asteroids, append 'P' for comets, start with comet type and a '/' for comet designations
        t_min_tdb (float or None, optional): Minimum date to cap the optical observation epochs. Defaults to None.
        t_max_tdb (float or None, optional): Maximum date to cap the optical observation epochs. Defaults to None.
        debias (bool, optional): Toggle for whether to use the debiasing scheme. Defaults to True.
        debias_lowres (bool, optional): Toggle for whether to use the low resolution debiasing scheme. Defaults to False.
        deweight (bool, optional): Toggle for whether to deweight observations. Defaults to True.
        eliminate (bool, optional): Toggle for whether to limit the number of observations per night. Defaults to False.
        max_obs_per_night (int, optional): Maximum number of observations to allow per night. Defaults to 5.
        verbose (bool, optional): Toggle for whether to print out information about the optical observation array. Defaults to False.

    Raises:
        ValueError: If deweight and eliminate are both True.
        ValueError: If deweight and eliminate are both False.
        ValueError: If debias and debias_lowres are both True.
        ValueError: If debias and debias_lowres are both False.

    Returns:
        obs_array_optical (np.ndarray): Optical observation data for the given body
        observer_codes_optical (tuple): Observer locations for each observation in obs_array_optical
    """
    if eliminate and deweight:
        raise ValueError('Cannot deweight and eliminate observations at the same time.')
    elif not eliminate and not deweight:
        print('WARNING: No deweighting or elimination scheme applied for observations during the same night.')
        # raise ValueError('Must deweight or eliminate observations.')
    if debias and debias_lowres:
        raise ValueError('Cannot debias and debias_lowres at the same time.')
    elif not debias and not debias_lowres:
        raise ValueError('Must debias or debias_lowres.')
    
    obs_array_optical, star_catalog_codes, observer_codes_optical, observation_type_codes, observer_program_codes = get_optical_data(body_id, optical_obs_file, t_min_tdb, t_max_tdb, verbose)
    if debias:
        obs_array_optical, star_catalog_codes, observer_codes_optical = apply_debiasing_scheme(obs_array_optical, star_catalog_codes, observer_codes_optical, debias_lowres, verbose)
    if old_weights:
        print('WARNING: Using old weighting scheme.')
        obs_array_optical = apply_weights_old(obs_array_optical, star_catalog_codes, observer_codes_optical, observation_type_codes, observer_program_codes, verbose)
    else:
        obs_array_optical = apply_weights(obs_array_optical, star_catalog_codes, observer_codes_optical, observation_type_codes, observer_program_codes, verbose)
    if deweight:
        obs_array_optical, star_catalog_codes, observer_codes_optical = deweight_obs(obs_array_optical, star_catalog_codes, observer_codes_optical, max_obs_per_night, verbose)
    if eliminate:
        obs_array_optical, star_catalog_codes, observer_codes_optical = eliminate_obs(obs_array_optical, star_catalog_codes, observer_codes_optical, max_obs_per_night, verbose)
    
    return obs_array_optical, observer_codes_optical
