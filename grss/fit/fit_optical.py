"""Optical observation handling for the GRSS orbit determination code"""
from astropy.time import Time
from astroquery.mpc import MPC
from astroquery.gaia import Gaia
import healpy as hp
import numpy as np
import pandas as pd
import spiceypy as spice

from .fit_utils import get_ra_from_hms, get_dec_from_dms, radec2icrf, icrf2radec, rec2lat
from .. import utils

__all__ = [ 'get_ades_optical_obs_array',
            'get_gaia_optical_obs_array',
            'get_mpc_optical_obs_array',
]

def get_space_based_observer_code(obs_time_mjd_tdb, first_80, last_80):
    """
    Get the full observer code for a space-based observation from the Minor Planet Center.

    Parameters
    ----------
    obs_time_mjd_tdb : float
        Observation time (MJD TDB)
    first_80 : str
        First 80 characters of the space-based observation
    last_80 : str
        Last 80 characters of the space-based observation

    Returns
    -------
    full_obs_code : tuple
        Full observer code for the space-based observation
    """
    code = first_80[77:80]
    if last_80[32] == '1':
        fac = 1.0
    elif last_80[32] == '2':
        fac = 1.495978707e8
    else:
        raise ValueError(f"Space-based observation units code not recognized. Given: {last_80[32]}")
    pos_x = float(last_80[34:45].replace(' ', ''))*fac
    pos_y = float(last_80[46:57].replace(' ', ''))*fac
    pos_z = float(last_80[58:69].replace(' ', ''))*fac
    lon, lat, rho = rec2lat(obs_time_mjd_tdb, pos_x, pos_y, pos_z)
    return (code, lon, lat, rho*1e3)

def get_optical_data(body_id, de_kernel_path, optical_obs_file=None, t_min_tdb=None,
                        t_max_tdb=None, verbose=False):
    # sourcery skip: low-code-quality
    """
    Get optical observation data from the Minor Planet Center for a desired small body
    and assemble the optical observations for the given body into an array for an orbit fit.

    Parameters
    ----------
    body_id : str/int
        Target id, numbers are interpreted as asteroids,
        append 'P' for comets, start with comet type and a '/' for comet designations
    de_kernel_path : str
        Path to the DE kernel file
    optical_obs_file : str, optional
        Filepath to the optical observations file, by default None
    t_min_tdb : float, optional
        Minimum time (MJD TDB) for observations to be included, by default None
    t_max_tdb : float, optional
        Maximum time (MJD TDB) for observations to be included, by default None
    verbose : bool, optional
        Flag to print out information about the observations, by default False

    Returns
    -------
    obs_array_optical : array
        Optical observation data for the given body
    star_catalog_codes : tuple
        Star catalog codes used for reducing each observation in obs_array_optical
    observer_codes_optical : tuple
        Observer locations for each observation in obs_array_optical
    observation_type_codes : tuple
        Observation types for each observation in obs_array_optical
    observer_program_codes : tuple
        Observatory programs for each observation in obs_array_optical
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
        # pylint: disable=no-member
        obs_raw = MPC.get_observations(body_id, get_mpcformat=True, cache=False)['obs']
    else:
        with open(optical_obs_file, 'r', encoding='utf-8') as file:
            obs_raw = file.readlines()
            obs_raw = [line.rstrip() for line in obs_raw] # remove \n at the end of each line
    obs_array_optical = np.zeros((len(obs_raw), 6))
    star_catalog_codes = []
    observer_codes_optical = []
    observation_type_codes = []
    observer_program_codes = []
    space_observatory_codes = [ '245', '249', '250', '258', '274', 'C49', 'C50', 'C51',
                                'C52', 'C53', 'C54', 'C55', 'C56', 'C57', 'C59',]
    non_geocentric_occultation_codes = ['275']
    unsupported_code_counter = 0
    unsupported_type_counter = 0
    out_of_range_counter = 0
    skip_counter = 0
    spice.furnsh(de_kernel_path)
    for i, data in enumerate(obs_raw):
        obs_type = 'P' if data[14] == ' ' else data[14]
        date = data[15:32]
        date_main = date[:-7].replace(' ','-')
        date_float = float(date[-7:])
        obs_time_mjd = Time(date_main, format='iso', scale='utc').utc.mjd+date_float
        obs_code = data[77:80]
        disallow_code = obs_code in non_geocentric_occultation_codes
        disallow_type = obs_type in ['R', 'r', 'V', 'v', 'O', 'X', 'x']
        disallow_line2 = obs_type in ['s']
        disallow_time = obs_time_mjd < t_min_utc or obs_time_mjd > t_max_utc
        if not disallow_code and not disallow_type and not disallow_time and not disallow_line2:
            observer_program = data[13]
            r_asc = get_ra_from_hms(data[32:44])
            dec = get_dec_from_dms(data[44:56])
            obs_array_optical[i, 0] = obs_time_mjd
            obs_array_optical[i, 1] = r_asc
            obs_array_optical[i, 2] = dec
            star_catalog = data[71]
            star_catalog_codes.append(star_catalog)
            if obs_code in space_observatory_codes:
                obs_time_mjd_tdb = Time(obs_time_mjd, format='mjd', scale='utc').tdb.mjd
                first_80 = data
                last_80 = data[80:] if optical_obs_file is None else obs_raw[i+1][:80]
                full_obs_code = get_space_based_observer_code(obs_time_mjd_tdb, first_80, last_80)
                observer_codes_optical.append(full_obs_code)
            else:
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
            elif disallow_time:
                increment = 1
                out_of_range_counter += increment
            else: # disallow_line2
                increment = 0
            skip_counter += increment
            obs_array_optical[i, :] = np.nan
            star_catalog_codes.append(np.nan)
            observer_codes_optical.append(np.nan)
            observation_type_codes.append(np.nan)
            observer_program_codes.append(np.nan)
    spice.kclear()
    if verbose:
        print(f"Skipped {int(skip_counter)} observations \n\t {unsupported_code_counter} of",
                f"which were non-geocentric occultations, \n\t {int(unsupported_type_counter)}",
                "were either roving or radar observations (radar is handled separately),",
                f"\n\t {out_of_range_counter} of which were outside the specified time range.")
    non_nan_idx = ~np.isnan(obs_array_optical[:, 0])
    star_catalog_codes = tuple(np.array(star_catalog_codes)[non_nan_idx])
    observer_codes_optical = tuple(np.array(observer_codes_optical, dtype=object)[non_nan_idx])
    observation_type_codes = tuple(np.array(observation_type_codes)[non_nan_idx])
    observer_program_codes = tuple(np.array(observer_program_codes)[non_nan_idx])
    obs_array_optical = obs_array_optical[non_nan_idx]
    return (obs_array_optical, star_catalog_codes,
            observer_codes_optical, observation_type_codes, observer_program_codes)

def get_ades_optical_obs_array(psv_obs_file, de_kernel_path, occultation_obs=False):
    """
    Assemble the optical observations for a given body into an array for an orbit fit,
    from the ADES PSV observation file.

    Parameters
    ----------
    psv_obs_file : str
        Path to the ADES PSV observation file
    de_kernel_path : str
        Path to the DE kernel file
    occultation_obs : bool, optional
        Flag for whether file is made of occultation measurement data, by default False

    Returns
    -------
    obs_array_optical : array
        Optical observation data for the given body
    observer_codes_optical : tuple
        Observer locations for each observation in obs_array_optical
    """
    if occultation_obs and de_kernel_path is None:
        raise ValueError(("Must specify path to DE kernel file to calculate"
                            "body-fixed coordinates for occultation observations."))
    obs_df = pd.read_csv(psv_obs_file, sep='|')
    obs_df.columns = obs_df.columns.str.strip()
    obs_array_optical = np.zeros((len(obs_df), 6))
    obs_times = Time(obs_df.obsTime.to_numpy(dtype=str),
                                    format='isot', scale='utc')
    obs_array_optical[:, 0] = obs_times.utc.mjd
    if occultation_obs:
        ra_star = obs_df.raStar.to_numpy(dtype=float)
        delta_ra = obs_df.deltaRA.to_numpy(dtype=float)
        obs_array_optical[:, 1] = ra_star+delta_ra
        dec_star = obs_df.decStar.to_numpy(dtype=float)
        delta_dec = obs_df.deltaDec.to_numpy(dtype=float)
        obs_array_optical[:, 2] = dec_star+delta_dec
    else:
        obs_array_optical[:, 1] = obs_df.ra.to_numpy(dtype=float)
        obs_array_optical[:, 2] = obs_df.dec.to_numpy(dtype=float)
    obs_array_optical[:, 1:3] *= 3600.0
    obs_array_optical[:, 3] = obs_df.rmsRA.to_numpy(dtype=float)
    obs_array_optical[:, 4] = obs_df.rmsDec.to_numpy(dtype=float)
    obs_array_optical[:, 5] = obs_df.rmsCorr.to_numpy(dtype=float)
    observer_codes_optical = []
    spice.furnsh(de_kernel_path)
    for i, obs in obs_df.iterrows():
        if str(obs['stn']) in {'275', 'S/C'}:
            pos_x = float(obs['pos1'])
            pos_y = float(obs['pos2'])
            pos_z = float(obs['pos3'])
            lon, lat, rho = rec2lat(obs_times[i].tdb.mjd, pos_x, pos_y, pos_z)
            observer_codes_optical.append((str(obs['stn']), lon, lat, 1e3*rho))
        else:
            observer_codes_optical.append(obs['stn'])
    spice.kclear()
    return obs_array_optical, tuple(observer_codes_optical)

def _get_gaia_query_results(body_id, release, verbose):
    """
    Submit a Gaia archive query for a given body ID from a specific
    Gaia data release.

    Parameters
    ----------
    body_id : str/int
        Target id, numbers are interpreted as asteroids,
        append 'P' for comets, start with comet type and a '/' for comet designations
    release : str, optional
        Gaia data release version database name (should use 'gaiadr3' if unsure)
    verbose : bool, optional
        Flag to print out information about the observations

    Returns
    -------
    res : astropy.table.Table
        Query results
    """
    table = 'sso_observation'
    if body_id.isdigit():
        match_str = f"WHERE ({release}.{table}.number_mp={body_id})"
    else:
        match_str = f"WHERE ({release}.{table}.denomination='{body_id.lower()}')"
    query = (
        "SELECT transit_id,denomination,number_mp,epoch_utc,epoch_err,"
        + "ra,dec,"
        + "ra_error_systematic,dec_error_systematic,ra_dec_correlation_systematic,"
        + "ra_error_random,dec_error_random,ra_dec_correlation_random,"
        + "x_gaia_geocentric,y_gaia_geocentric,z_gaia_geocentric,"
        + "vx_gaia_geocentric,vy_gaia_geocentric,vz_gaia_geocentric"
        + f" FROM {release}.{table} {match_str} ORDER BY epoch_utc ASC"
    )
    job = Gaia.launch_job_async(query, dump_to_file=False,background=True)
    res = job.get_results()
    res.sort('epoch_utc')
    if verbose:
        print(f"Found {len(res)} observations from {release}")
    return res

def get_gaia_optical_obs_array(body_id, de_kernel_path, t_min_tdb=None,
                                t_max_tdb=None, gaia_dr='gaiadr3', verbose=False):
    """
    Assemble the optical observations for a given body from Gaia DR3.

    Parameters
    ----------
    body_id : str/int
        Target id, numbers are interpreted as asteroids,
        append 'P' for comets, start with comet type and a '/' for comet designations
    de_kernel_path : str
        Path to the DE kernel file
    t_min_tdb : float, optional
        Minimum time (MJD TDB) for observations to be included, by default None
    t_max_tdb : float, optional
        Maximum time (MJD TDB) for observations to be included, by default None
    gaia_dr : str, optional
        Gaia data release version database name, by default 'gaiadr3'
    verbose : bool, optional
        Flag to print out information about the observations, by default False

    Returns
    -------
    obs_array_optical : array
        Optical observation data for the given body
    observer_codes_optical : tuple
        Observer locations for each observation in obs_array_optical
    """
    if t_min_tdb is None:
        t_min_tdb = -np.inf
    if t_max_tdb is None:
        t_max_tdb = np.inf
    # get gaia query results
    res = _get_gaia_query_results(body_id, release=gaia_dr, verbose=verbose)
    spice.furnsh(de_kernel_path)
    au2km = 149597870.7
    obs_array = np.nan*np.ones((len(res), 6))
    observer_codes = []
    curr_transit = int(-1e6)
    for i, data in enumerate(res):
        if curr_transit != data['transit_id']:
            curr_transit = data['transit_id']
            transit_count = 1
            while i+transit_count < len(res) and res[i+transit_count]['transit_id'] == curr_transit:
                transit_count += 1
        obs_time = Time(data['epoch_utc'] + 55197.0, format='mjd', scale='utc')
        if obs_time.tdb.mjd < t_min_tdb or obs_time.tdb.mjd > t_max_tdb:
            continue
        cosdec = np.cos(data['dec']*np.pi/180)
        obs_array[i, 0] = obs_time.utc.mjd
        obs_array[i, 1] = data['ra']*3600
        obs_array[i, 2] = data['dec']*3600
        ra_sig_sys = data['ra_error_systematic']/cosdec/1000*transit_count
        ra_sig_rand = data['ra_error_random']/cosdec/1000
        dec_sig_sys = data['dec_error_systematic']/1000*transit_count
        dec_sig_rand = data['dec_error_random']/1000
        corr_sys = data['ra_dec_correlation_systematic']
        corr_rand = data['ra_dec_correlation_random']
        cov_sys = np.array([[ra_sig_sys**2, corr_sys*ra_sig_sys*dec_sig_sys],
                            [corr_sys*ra_sig_sys*dec_sig_sys, dec_sig_sys**2]])
        cov_rand = np.array([[ra_sig_rand**2, corr_rand*ra_sig_rand*dec_sig_rand],
                            [corr_rand*ra_sig_rand*dec_sig_rand, dec_sig_rand**2]])
        cov = cov_sys + cov_rand
        ra_sig = cov[0,0]**0.5
        dec_sig = cov[1,1]**0.5
        corr = cov[0,1]/ra_sig/dec_sig
        obs_array[i, 3] = ra_sig
        obs_array[i, 4] = dec_sig
        obs_array[i, 5] = corr
        pos_x = data['x_gaia_geocentric']*au2km
        pos_y = data['y_gaia_geocentric']*au2km
        pos_z = data['z_gaia_geocentric']*au2km
        lon, lat, rho = rec2lat(obs_time.tdb.mjd, pos_x, pos_y, pos_z)
        observer_codes.append(('258', lon, lat, rho*1e3))
    spice.kclear()
    non_nan_idx = ~np.isnan(obs_array[:, 0])
    obs_array = obs_array[non_nan_idx, :]
    if verbose:
        print(f"\t Added {len(obs_array)} of those observations.")
    return obs_array, tuple(observer_codes)

def debias_obs(r_asc, dec, epoch, catalog, biasdf, nside=256):
    """
    Calculate the bias on an optical observation given the epoch and the star catalog.

    Parameters
    ----------
    r_asc : float
        Right Ascension in radians.
    dec : float
        Declination in radians.
    epoch : float
        Epoch of observation in JD TT.
    catalog : str
        Star catalog MPC code.
    biasdf : pandas DataFrame
        Bias dataframe of MPC star catalogs.
    nside : int, optional
        Healpix nside parameter, by default 256

    Returns
    -------
    ra_deb : float
        Debiased right ascension in radians.
    dec_deb : float
        Debiased declination in radians.
    """
    j2000_jd = 2451545.0
    # arcseconds to radians
    as2rad = 1/3600*np.pi/180
    # find pixel from RADEC
    idx = hp.ang2pix(nside, np.rad2deg(r_asc), np.rad2deg(dec), nest=False, lonlat=True)
    # find catalog data in pandas Data Frame
    colnames = [f'{catalog}_ra', f'{catalog}_dec', f'{catalog}_pm_ra', f'{catalog}_pm_dec']
    bias = biasdf[colnames].iloc[idx]
    # time from epoch in Julian years
    dt_jy = (epoch-j2000_jd)/365.25
    # bias correction
    ddec = (bias[colnames[1]]+dt_jy*bias[colnames[3]]/1000)*as2rad
    dec_deb = dec-ddec
    dra = (bias[colnames[0]]+dt_jy*bias[colnames[2]]/1000)*as2rad / np.cos(dec)
    ra_deb = r_asc-dra
    # Quadrant correction
    xyz = radec2icrf(ra_deb, dec_deb, deg=False)
    ra_deb, dec_deb = icrf2radec(xyz[0], xyz[1], xyz[2], deg=False)
    return ra_deb, dec_deb

def apply_debiasing_scheme(obs_array_optical, star_catalog_codes, observer_codes_optical,
                            lowres=False, verbose=False):
    # sourcery skip: low-code-quality
    """
    Apply the Eggl et al. (2020) debiasing scheme to the observation array.
    (paper from https://doi.org/10.1016/j.icarus.2019.113596)

    Parameters
    ----------
    obs_array_optical : array
        Optical observation data for the given body
    star_catalog_codes : tuple
        Star catalog codes used for reducing each observation in obs_array_optical
    observer_codes_optical : tuple
        Observer locations for each observation in obs_array_optical
    lowres : bool, optional
        Flag to debias the observations using the low resolution scheme, by default False
    verbose : bool, optional
        Flag to print out information about the observations, by default False

    Returns
    -------
    obs_array_optical : array
        Optical observation data for the given body with debiasing applied
    star_catalog_codes : tuple
        Star catalog codes used for reducing each observation in obs_array_optical
    observer_codes_optical : tuple
        Observer locations for each observation in obs_array_optical
    """
    star_catalog_codes = list(star_catalog_codes)
    observer_codes_optical = list(observer_codes_optical)
    # MPC star_catalog_codes codes, https://www.minorplanetcenter.net/iau/info/CatalogueCodes.html
    biased_catalogs = ['a', 'b', 'c', 'd', 'e', 'g', 'i', 'j', 'l', 'm', 'n', 'o', 'p',
                        'q', 'r', 't', 'u', 'v', 'w', 'L', 'N', 'Q', 'R', 'S', 'U', 'Y']
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
    debias_lowres_path = utils.grss_path+'/debias/lowres_data'
    debias_hires_path = utils.grss_path+'/debias/hires_data'
    if lowres:
        biasfile = f'{debias_lowres_path}/bias.dat'
        nside = 64
    else:
        biasfile = f'{debias_hires_path}/bias.dat'
        nside = 256
    biasdf = pd.read_csv(biasfile,sep=r'\s+',skiprows=23,names=columns)
    unbias_counter = 0
    debias_counter = 0
    for i, row in enumerate(obs_array_optical):
        obs_time_jd = Time(row[0], format='mjd', scale='utc').tt.jd
        r_asc = row[1]/3600*np.pi/180
        dec = row[2]/3600*np.pi/180
        star_catalog = star_catalog_codes[i]
        if star_catalog in unbiased_catalogs:
            unbias_counter += 1
            continue
        elif star_catalog in biased_catalogs:
            debias_counter += 1
            # apply debiasing from Eggl et al. 2020, https://doi.org/10.1016/j.icarus.2019.113596
            ra_new, dec_new = debias_obs(r_asc, dec, obs_time_jd,
                                star_catalog, biasdf, nside=nside)
            obs_array_optical[i, 1] = ra_new*180/np.pi*3600
            obs_array_optical[i, 2] = dec_new*180/np.pi*3600
            margin = 1.5
            if (verbose and
                (np.rad2deg(abs(r_asc - ra_new))*3600 >= margin or
                np.rad2deg(abs(dec - dec_new))*3600 >= margin)):
                print(f"Debiased observation {i+1} at MJD {row[0]} UTC from observatory",
                        f"{observer_codes_optical[i]} using catalog {star_catalog}.",
                        f"RA bias = {np.rad2deg(r_asc-ra_new)*3600:0.4f} arcsec,",
                        f"DEC bias = {np.rad2deg(dec-dec_new)*3600:0.4f} arcsec")
    if verbose:
        no_bias_counter = len(obs_array_optical)-unbias_counter-debias_counter
        print(f"No debiasing needed for {unbias_counter} observations. Debiased {debias_counter}",
                f"observations. No biasing information for {no_bias_counter} observations.")
    non_nan_idx = ~np.isnan(obs_array_optical[:, 0])
    star_catalog_codes = tuple(np.array(star_catalog_codes)[non_nan_idx])
    observer_codes_optical = tuple(np.array(observer_codes_optical, dtype=object)[non_nan_idx])
    obs_array_optical = obs_array_optical[non_nan_idx]
    return obs_array_optical, star_catalog_codes, observer_codes_optical

def apply_weights(obs_array_optical, star_catalog_codes, observer_codes_optical,
                    observation_type_codes, observer_program_codes, verbose=False):
    """
    Apply the weights from Vereš et al. (2017) to the observation array.
    (paper from https://doi.org/10.1016/j.icarus.2017.05.021)

    Parameters
    ----------
    obs_array_optical : array
        Optical observation data for the given body
    star_catalog_codes : tuple
        Star catalog codes used for reducing each observation in obs_array_optical
    observer_codes_optical : tuple
        Observer locations for each observation in obs_array_optical
    observation_type_codes : tuple
        Observation types for each observation in obs_array_optical
    observer_program_codes : tuple
        Observatory programs for each observation in obs_array_optical
    verbose : bool, optional
        Flag to print out information about the observations, by default False

    Returns
    -------
    obs_array_optical : array
        Optical observation data for the given body with weights applied
    

    Raises
    ------
    ValueError
        _description_
    ValueError
        _description_
    """
    # weighting scheme from vereš et al. 2017, obswgt_veres2017.inp
    default_weight_counter = 0
    for i in range(len(obs_array_optical)):
        star_catalog = star_catalog_codes[i]
        obs_code = observer_codes_optical[i]
        if isinstance(obs_code, tuple):
            obs_code = obs_code[0]
        obs_date_mjd = obs_array_optical[i, 0]
        obs_type = observation_type_codes[i]
        obs_program = observer_program_codes[i]
        if obs_type in ['P', 'A', 'N', ' ']:
            if obs_date_mjd <= Time('1890-01-01', format='iso', scale='utc').utc.mjd:
                obs_array_optical[i, 3:5] = 10.0
            elif (  obs_date_mjd > Time('1890-01-01', format='iso', scale='utc').utc.mjd and
                    obs_date_mjd <= Time('1950-01-01', format='iso', scale='utc').utc.mjd):
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
            if obs_code in ['275']: # non-geocentric occultation
                obs_array_optical[i, 3:5] = 0.2
            elif obs_code in ['250']: # HST
                obs_array_optical[i, 3:5] = 1.3
            elif obs_code in ['249', 'C49', 'C50']: # SOHO, STEREO-A, STEREO-B
                obs_array_optical[i, 3:5] = 60.0
            elif obs_code in ['C51']: # WISE
                obs_array_optical[i, 3:5] = 1.0
            elif obs_code in ['C53']: # NEOSSat
                obs_array_optical[i, 3:5] = 3.0
            elif obs_code in ['C57']: # TESS
                obs_array_optical[i, 3:5] = 5.0
            elif obs_code in ['258']: # Gaia
                obs_array_optical[i, 3:5] = 0.5
            else:
                obs_array_optical[i, 3:5] = 1.5
        elif obs_type in ['c', 'C', 'V', 'n', 'B', 'X', 'x']:
            obs_array_optical[i, 3:5] = 1.0
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
                elif obs_code in ['691', '291']:
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
                elif obs_code in ['T05', 'T08', 'T07', 'M22', 'W68']:
                    obs_array_optical[i, 3:5] = 0.8
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
                    if star_catalog in ['L', 't', 'U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.3
                elif obs_code in ['673']:
                    if obs_program == '1':
                        obs_array_optical[i, 3:5] = 0.3
                elif obs_code in ['645']:
                    obs_array_optical[i, 3:5] = 0.3
                elif obs_code in ['689']:
                    obs_array_optical[i, 3:5] = 0.5
                elif obs_code in ['J04', 'Z84']:
                    if star_catalog in ['U', 'V', 'W', 'X', 't', 'L', 'q', 'r', 'u', 'e']:
                        obs_array_optical[i, 3:5] = 0.5
                    if obs_code in ['J04']:
                        if obs_program == '2' and star_catalog in ['t', 'q', 'U', 'V', 'W', 'X']:
                            obs_array_optical[i, 3:5] = 0.3
                        elif obs_program == '#':
                            obs_array_optical[i, 3:5] = 1.0
                    elif obs_code in ['Z84']:
                        if obs_program == '1' and star_catalog in ['U', 'V', 'W', 'X']:
                            obs_array_optical[i, 3:5] = 0.3
                elif obs_code in ['608']:
                    obs_array_optical[i, 3:5] = 0.6
                elif obs_code in ['W84']:
                    obs_array_optical[i, 3:5] = 0.5
                    if obs_program == '&':
                        obs_array_optical[i, 3:5] = 1.0
                elif obs_code in ['950']:
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
                elif obs_code in ['K91', 'K92', 'K93', 'Q63', 'Q64', 'V37', 'V39',
                                    'W85', 'W86', 'W87', 'Z24', 'Z31']:
                    obs_array_optical[i, 3:5] = 0.4
                    if obs_program == '0':
                        obs_array_optical[i, 3:5] = 0.3
                elif obs_code in ['Y28']:
                    if star_catalog in ['t', 'U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.3
                elif obs_code in ['309']:
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
                    if obs_program == '2':
                        if star_catalog in ['U', 'V', 'W', 'X']:
                            obs_array_optical[i, 3:5] = 0.2
                        elif star_catalog in ['t', 'q']:
                            obs_array_optical[i, 3:5] = 0.3
                elif obs_code in ['C65']:
                    if obs_program == '3' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.2
                elif obs_code in ['094']:
                    if obs_program == '9' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.5
                elif obs_code in ['Z18']:
                    if obs_program == '1' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.2
                elif obs_code in ['181']:
                    if obs_program == '1' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.5
                elif obs_code in ['G37']:
                    if obs_program == '5' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.2
                elif obs_code in ['D20']:
                    if obs_program == '1' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.5
                elif obs_code in ['Z18']:
                    if obs_program == '1' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.1
                elif obs_code in ['N50']:
                    if obs_program == '2' and star_catalog in ['U', 'V', 'W', 'X']:
                        obs_array_optical[i, 3:5] = 0.3
                elif obs_code in ['I41']:
                    if obs_program == 'Z':
                        obs_array_optical[i, 3:5] = 1.0
                else:
                    default_weight_counter += 1
                    obs_array_optical[i, 3:5] = 1.0
            if star_catalog == ' ':
                obs_array_optical[i, 3:5] = 1.5
            if obs_code in ['270']:
                obs_array_optical[i, 3:5] = 2.0
        else:
            raise ValueError(f"apply_weights: Observation type {obs_type} not recognized.")
    if verbose:
        print(f"Applied default weight of 1 arcsec to {default_weight_counter} CCD observations")
        # find where the weights are 0
        zero_weight_indices_ra = np.where(obs_array_optical[:, 3] == 0)[0]
        zero_weight_indices_dec = np.where(obs_array_optical[:, 4] == 0)[0]
        zero_weight_indices = np.concatenate((zero_weight_indices_ra, zero_weight_indices_dec))
        if len(zero_weight_indices) > 0:
            raise ValueError(f"apply_weights: Found {len(zero_weight_indices)} observations",
                                    "with zero weight. Check the star_catalog_codes and observer",
                                    "location codes and make sure proper weights are applied.")
    return obs_array_optical

def deweight_obs(obs_array_optical, star_catalog_codes, observer_codes_optical,
                    eff_obs_per_night=5, verbose=False):
    """
    Deweight observations that occur on the same night at the same observatory.

    Parameters
    ----------
    obs_array_optical : array
        Optical observation data for the given body
    star_catalog_codes : tuple
        Star catalog codes used for reducing each observation in obs_array_optical
    observer_codes_optical : tuple
        Observer locations for each observation in obs_array_optical
    eff_obs_per_night : int, optional
        Effective number of observations per night, by default 5
    verbose : bool, optional
        Flag to print out information about the observations, by default False

    Returns
    -------
    obs_array_deweighted : array
        Deweighted Optical observation data for the given body
    catalog_deweighted : tuple
        Star catalog codes used for reducing each observation in obs_array_deweighted
    observer_loc_deweighted : tuple
        Observer locations for each observation in obs_array_deweighted
    """
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
                start = i-night_count+1
                end = i+1
                obs_array_deweighted[start:end,3:5] = factor * obs_array_optical[start:end,3:5]
        else:
            night_count = 1
    if verbose:
        print(f"Deweighted {deweight_count} observations as part of deweighting scheme.")
    return obs_array_deweighted, tuple(catalog_deweighted), tuple(observer_loc_deweighted)

def eliminate_obs(obs_array_optical, star_catalog_codes, observer_codes_optical,
                    max_obs_per_night=5, verbose=False):
    """
    Limit the number of observations that occur on the same night at the same observatory.

    Parameters
    ----------
    obs_array_optical : array
        Optical observation data for the given body
    star_catalog_codes : tuple
        Star catalog codes used for reducing each observation in obs_array_optical
    observer_codes_optical : tuple
        Observer locations for each observation in obs_array_optical
    max_obs_per_night : int, optional
        Maximum number of observations per night, by default 5
    verbose : bool, optional
        Flag to print out information about the observations, by default False

    Returns
    -------
    obs_array_eliminated : array
        Eliminated Optical observation data for the given body
    catalog_eliminated : tuple
        Star catalog codes used for reducing each observation in obs_array_eliminated
    observer_loc_eliminated : tuple
        Observer locations for each observation in obs_array_eliminated
    """
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
        print(f"Removed {len(obs_array_optical)-len(obs_array_eliminated)} observations",
                    "as part of elimination scheme.")
    return obs_array_eliminated, tuple(catalog_eliminated), tuple(observer_loc_eliminated)

def get_mpc_optical_obs_array(body_id, de_kernel_path, optical_obs_file=None, t_min_tdb=None,
                            t_max_tdb=None, debias_hires=True, debias_lowres=False, deweight=True,
                            eliminate=False, num_obs_per_night=5, verbose=False):
    """
    Assemble the optical observations for a given body into an array for an orbit fit.

    Parameters
    ----------
    body_id : str/int
        Target id, numbers are interpreted as asteroids,
        append 'P' for comets, start with comet type and a '/' for comet designations
    de_kernel_path : str
        Path to the JPL DE kernel
    optical_obs_file : str, optional
        Filepath to the optical observations file, by default None
    t_min_tdb : float, optional
        Minimum time (MJD TDB) for observations to be included, by default None
    t_max_tdb : float, optional
        Maximum time (MJD TDB) for observations to be included, by default None
    debias_hires : bool, optional
        Flag to debias the observations using the high resolution scheme, by default True
    debias_lowres : bool, optional
        Flag to debias the observations using the low resolution scheme, by default False
    deweight : bool, optional
        Flag to deweight the observations, by default True
    eliminate : bool, optional
        Flag to eliminate too many observations from the same night, by default False
    num_obs_per_night : int, optional
        Effective/Maximum number of effective observations per night, by default 5
    verbose : bool, optional
        Flag to print out information about the observations, by default False

    Returns
    -------
    obs_array_optical : array
        Optical observation data for the given body
    observer_codes_optical : tuple
        Observer locations for each observation in obs_array_optical

    Raises
    ------
    ValueError
        If deweight and eliminate are both True.
    ValueError
        If deweight and eliminate are both False.
    ValueError
        If debias_hires and debias_lowres are both True.
    ValueError
        If debias_hires and debias_lowres are both False.
    """
    if eliminate and deweight:
        raise ValueError('Cannot deweight and eliminate observations at the same time.')
    if not eliminate and not deweight:
        print("WARNING: No deweighting or elimination scheme applied",
                "for observations during the same night.")
    if debias_hires and debias_lowres:
        raise ValueError('Cannot debias_hires and debias_lowres at the same time.')
    if not debias_hires and not debias_lowres:
        print("WARNING: No debiasing scheme applied for observations.")
    (obs_array_optical, star_catalog_codes,
        observer_codes_optical, observation_type_codes,
        observer_program_codes) = get_optical_data(body_id, de_kernel_path, optical_obs_file,
                                    t_min_tdb, t_max_tdb, verbose)
    if debias_hires or debias_lowres:
        (obs_array_optical, star_catalog_codes,
            observer_codes_optical) = apply_debiasing_scheme(obs_array_optical, star_catalog_codes,
                                        observer_codes_optical, debias_lowres, verbose)
    obs_array_optical = apply_weights(obs_array_optical, star_catalog_codes,
                            observer_codes_optical, observation_type_codes,
                            observer_program_codes, verbose)
    if deweight:
        (obs_array_optical, star_catalog_codes,
            observer_codes_optical) = deweight_obs(obs_array_optical, star_catalog_codes,
                                        observer_codes_optical, num_obs_per_night, verbose)
    if eliminate:
        (obs_array_optical, star_catalog_codes,
            observer_codes_optical) = eliminate_obs(obs_array_optical, star_catalog_codes,
                                        observer_codes_optical, num_obs_per_night, verbose)
    return obs_array_optical, observer_codes_optical
