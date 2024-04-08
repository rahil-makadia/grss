"""Optical observation handling for the GRSS orbit determination code"""
from astropy.time import Time
import requests
from astroquery.gaia import Gaia
import healpy as hp
import numpy as np
import pandas as pd

from .fit_utils import radec2icrf, icrf2radec
from .. import utils

__all__ = [ 'get_ades_optical_obs_array',
            'get_gaia_optical_obs_array',
            'get_mpc_optical_obs_array',
]

def get_mpc_raw_data(tdes):
    """
    Get observation data from MPC observations API for desired small body

    Parameters
    ----------
    tdes : str/int
        IMPORTANT: must be the designation of the small body, not the name.

    Returns
    -------
    raw_data : dict
        JSON output of small body information query from JPL small-body radar API
    """
    response = requests.get("https://data.minorplanetcenter.net/api/get-obs",
                            json={"desigs": [f"{tdes}"], "output_format":["XML"]},
                            timeout=60)
    if response.ok:
        obs_data = response.json()[0]["XML"]
    else:
        print("Error getting MPC XML data: ", response.status_code, response.content)
    return obs_data

def _ades_mode_check(df):
    """
    Check the mode values in the ADES data frame.

    Parameters
    ----------
    df : pandas DataFrame
        ADES data frame

    Raises
    ------
    ValueError
        If the mode values are invalid
    """
    # from https://www.minorplanetcenter.net/iau/info/ADESFieldValues.html
    valid_modes = ['CCD', 'CMO', 'VID', 'PHO', 'ENC', 'PMT', 'MIC', 'MER', 'TDI', 'OCC', 'UNK']
    if not df['mode'].isin(valid_modes).all():
        raise ValueError(f"Invalid mode in the data: {df['mode'].unique()}.\n"
                        f"Acceptable modes are {valid_modes}.")
    if df['mode'].isin(['UNK']).any():
        print("\tWARNING: At least one unknown obs mode in the data.")
    return None

def _ades_astCat_check(df):
    """
    Check the astCat values in the ADES data frame.

    Parameters
    ----------
    df : pandas DataFrame
        ADES data frame

    Raises
    ------
    ValueError
        If the astCat values are invalid
    """
    # from https://www.minorplanetcenter.net/iau/info/ADESFieldValues.html
    valid_cats = [
        'Gaia_Int', 'PS1_DR1', 'PS1_DR2', 'ATLAS2',
        'Gaia3', 'Gaia3E', 'Gaia2', 'Gaia1', 'Gaia_2016',
        'NOMAD', 'PPMXL', 'UBSC', 'UCAC5', 'UCAC4', 'URAT1', '2MASS'
    ]
    deprecated_cats = [
        'UCAC3', 'UCAC2', 'UCAC1',
        'USNOB1', 'USNOA2', 'USNOSA2', 'USNOA1', 'USNOSA1',
        'Tyc2', 'Tyc1', 'Hip2', 'Hip1', 'ACT',
        'GSCACT', 'GSC2.3', 'GSC2.2', 'GSC1.2', 'GSC1.1', 'GSC1.0', 'GSC',
        'SDSS8', 'SDSS7', 'CMC15', 'CMC14', 'SSTRC4', 'SSTRC1',
        'MPOSC3', 'PPM', 'AC', 'SAO1984', 'SAO', 'AGK3', 'FK4',
        'ACRS', 'LickGas', 'Ida93', 'Perth70', 'COSMOS',
        'Yale', 'ZZCAT', 'IHW', 'GZ', 'UNK']
    if df['astCat'].isin(deprecated_cats).any():
        print("\tWARNING: At least one deprecated star catalog in the data.")
    if not df['astCat'].isin(valid_cats+deprecated_cats).all():
        raise ValueError(f"At least one invalid astCat in the data: {df['astCat'].unique()}.\n"
                        f"Acceptable astCat values are {valid_cats+deprecated_cats}.")

def get_optical_obs_data(body_id, optical_obs_file=None, t_min_tdb=None,
                        t_max_tdb=None, verbose=False):
    # sourcery skip: low-code-quality
    """
    Get optical observation data from the Minor Planet Center for a desired small body
    and assemble the optical observations for the given body into an array for an orbit fit.

    Parameters
    ----------
    body_id : str
        Small body designation (any value accepted by the MPC)
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
    obs_df : pandas DataFrame
        Optical observation data for the given body
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
        obs_data = get_mpc_raw_data(body_id)
    else:
        obs_data = optical_obs_file
    # from table 16 of ADES_Description.pdf from the ADES-Master git repo
    # comments correspond to future GRSS capabilities.
    columns_to_keep = {
        'permID':'str', 'provID': 'str', 'mode': 'str', 'stn': 'str', 'prog': 'str',
        'obsTime': 'str', #'rmsTime': 'Float64',
        'ra': 'Float64', 'dec': 'Float64',
        'rmsRA': 'Float64', 'rmsDec': 'Float64', 'rmsCorr': 'Float64', 'astCat': 'str',
        'trx': 'str', 'rcv': 'str', 'delay': 'Float64', 'doppler': 'Float64',
        'rmsDelay': 'Float64', 'rmsDoppler': 'Float64', 'com': 'Int64', 'freq': 'Float64',
        'raStar': 'Float64', 'decStar': 'Float64', 'deltaRA': 'Float64', 'deltaDec': 'Float64',
        'sys': 'str', 'ctr': 'Int64', 'pos1': 'Float64', 'pos2': 'Float64', 'pos3': 'Float64',
    }
    residual_columns_to_add = {
        'resRA': 'Float64', 'resDec': 'Float64', 'selAst': 'str',
        'sigRA': 'Float64', 'sigDec': 'Float64', 'sigCorr': 'Float64', 'sigTime': 'Float64',
        'biasRA': 'Float64', 'biasDec': 'Float64', 'biasTime': 'Float64',
        'resDelay': 'Float64', 'resDoppler': 'Float64', 'selDelay': 'str', 'selDoppler': 'str',
        'sigDelay': 'Float64', 'sigDoppler': 'Float64',
    }
    utility_columns_to_add = {
        'obsTimeMJD': 'Float64', 'cosDec': 'Float64',
    }
    # do the list addition instead of pipe operator to keep the order
    columns_in_df = (
        list(columns_to_keep.keys()) +
        list(residual_columns_to_add.keys()) +
        list(utility_columns_to_add.keys())
    )
    types = columns_to_keep | residual_columns_to_add | utility_columns_to_add
    # read in the data and add extra columns if not present
    obs_df = (
        pd.read_xml(obs_data)
        .reindex(columns=columns_in_df, fill_value=np.nan) # add columns if they are not present
        .astype(types)
        # compute the obsTimeMJD time from obsTime
        .assign(obsTimeMJD=lambda x:Time(x['obsTime'].to_list(),format='isot',scale='utc').utc.mjd)
    )
    if verbose:
        print(f"Read in {len(obs_df)} observations from the MPC.")
    _ades_mode_check(obs_df)
    _ades_astCat_check(obs_df)
    # filter the data based on the time range
    obs_df.query(f"{t_min_utc} <= obsTimeMJD <= {t_max_utc}", inplace=True)
    # filter the data based on the observation station
    stn_to_remove = ('247', '270') # roving observatories
    obs_df.query(f"stn not in {stn_to_remove}", inplace=True)
    # for all indices with OCC mode compute ra and dec as raStar+deltaRA and decStar+deltaDec
    occ_idx = obs_df.query("mode == 'OCC'").index
    obs_df.loc[occ_idx, 'dec'] = obs_df.loc[occ_idx, 'decStar']+obs_df.loc[occ_idx, 'deltaDec']
    obs_df['cosDec'] = np.cos(obs_df['dec']*np.pi/180)
    obs_df.loc[occ_idx, 'ra'] = (obs_df.loc[occ_idx, 'raStar']+
                                    obs_df.loc[occ_idx, 'deltaRA']/obs_df.loc[occ_idx, 'cosDec'])
    if verbose:
        print(f"\tFiltered to {len(obs_df)} observations that satisfy the "
                "time range constraints and are from accepted observatories.")
    return obs_df

def get_ades_optical_obs_array(psv_obs_file, occultation_obs=False):
    """
    Assemble the optical observations for a given body into an array for an orbit fit,
    from the ADES PSV observation file.

    Parameters
    ----------
    psv_obs_file : str
        Path to the ADES PSV observation file
    occultation_obs : bool, optional
        Flag for whether file is made of occultation measurement data, by default False

    Returns
    -------
    obs_array_optical : array
        Optical observation data for the given body
    observer_codes_optical : tuple
        Observer locations for each observation in obs_array_optical
    """
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
    for _, obs in obs_df.iterrows():
        if str(obs['stn']) in {'275', 'S/C'}:
            pos_x = float(obs['pos1'])*1e3
            pos_y = float(obs['pos2'])*1e3
            pos_z = float(obs['pos3'])*1e3
            observer_codes_optical.append((str(obs['stn']), pos_x, pos_y, pos_z, 0, 0, 0))
        else:
            observer_codes_optical.append(obs['stn'])
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
    # change the column names to be lowercase
    for col_name in res.colnames:
        res.rename_column(col_name, col_name.lower())
    res.sort('epoch_utc')
    if verbose:
        print(f"Found {len(res)} observations from {release}")
    return res

def get_gaia_optical_obs_array(body_id, t_min_tdb=None,
                                t_max_tdb=None, gaia_dr='gaiadr3', verbose=False):
    """
    Assemble the optical observations for a given body from Gaia DR3.

    Parameters
    ----------
    body_id : str/int
        Target id, numbers are interpreted as asteroids,
        append 'P' for comets, start with comet type and a '/' for comet designations
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
    au2m = 149597870700
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
        ra_sig_sys = data['ra_error_systematic']/cosdec/1000
        ra_sig_rand = data['ra_error_random']/cosdec/1000
        dec_sig_sys = data['dec_error_systematic']/1000
        dec_sig_rand = data['dec_error_random']/1000
        corr_sys = data['ra_dec_correlation_systematic']
        corr_rand = data['ra_dec_correlation_random']
        cov_sys = np.array([[ra_sig_sys**2, corr_sys*ra_sig_sys*dec_sig_sys],
                            [corr_sys*ra_sig_sys*dec_sig_sys, dec_sig_sys**2]])
        cov_rand = np.array([[ra_sig_rand**2, corr_rand*ra_sig_rand*dec_sig_rand],
                            [corr_rand*ra_sig_rand*dec_sig_rand, dec_sig_rand**2]])
        cov = transit_count*cov_sys + cov_rand
        ra_sig = cov[0,0]**0.5
        dec_sig = cov[1,1]**0.5
        corr = cov[0,1]/ra_sig/dec_sig
        obs_array[i, 3] = ra_sig
        obs_array[i, 4] = dec_sig
        obs_array[i, 5] = corr
        pos_x = data['x_gaia_geocentric']*au2m
        pos_y = data['y_gaia_geocentric']*au2m
        pos_z = data['z_gaia_geocentric']*au2m
        observer_codes.append(('258', pos_x, pos_y, pos_z, 0, 0, 0))
    non_nan_idx = ~np.isnan(obs_array[:, 0])
    obs_array = obs_array[non_nan_idx, :]
    if verbose:
        print(f"\t Added {len(obs_array)} of those observations.")
    return obs_array, tuple(observer_codes)

def debias_obs(r_asc_vals, dec_vals, epoch_vals, catalog, biasdf, nside):
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
    nside : int
        Healpix nside parameter.

    Returns
    -------
    ra_deb : float
        Debiased right ascension in radians.
    dec_deb : float
        Debiased declination in radians.
    """
    j2000_jd = 2451545.0
    # find pixel from RADEC
    idx_vals = hp.ang2pix(nside, r_asc_vals, dec_vals, nest=False, lonlat=True)
    # find catalog data in pandas Data Frame
    colnames = [f'{catalog}_ra', f'{catalog}_dec', f'{catalog}_pm_ra', f'{catalog}_pm_dec']
    bias_vals = biasdf[colnames].iloc[idx_vals].values
    # time from epoch in Julian years
    dt_jy_vals = (epoch_vals-j2000_jd)/365.25
    # compute biases
    ddec_vals = (bias_vals[:, 1]+dt_jy_vals*bias_vals[:, 3]/1000)
    dra_cos_dec_vals = (bias_vals[:, 0]+dt_jy_vals*bias_vals[:, 2]/1000)
    return dra_cos_dec_vals, ddec_vals

def apply_debiasing_scheme(obs_df, lowres, verbose):
    # sourcery skip: low-code-quality
    """
    Apply the Eggl et al. (2020) debiasing scheme to the observation array.
    (paper from https://doi.org/10.1016/j.icarus.2019.113596)

    Parameters
    ----------
    obs_df : pandas DataFrame
        Optical observation data for the given body
    lowres : bool
        Flag to debias the observations using the low resolution scheme
    verbose : bool
        Flag to print out information about the observations

    Returns
    -------
    obs_df : pandas DataFrame
        Optical observation data for the given body with debiasing applied
    """
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
    # from ADES-Master/Python/bin/packUtil.py
    ades_catalog_map = {
        ' ': 'UNK',
        'a': 'USNOA1',
        'b': 'USNOSA1',
        'c': 'USNOA2',
        'd': 'USNOSA2',
        'e': 'UCAC1',
        'f': 'Tyc1',
        'g': 'Tyc2',
        'h': 'GSC1.0',
        'i': 'GSC1.1',
        'j': 'GSC1.2',
        'k': 'GSC2.2',
        'l': 'ACT',
        'm': 'GSCACT',
        'n': 'SDSS8',
        'o': 'USNOB1',
        'p': 'PPM',
        'q': 'UCAC4',
        'r': 'UCAC2',
        's': 'USNOB2',  # USNOB2 missing on ADES web page
        't': 'PPMXL',
        'u': 'UCAC3',
        'v': 'NOMAD',
        'w': 'CMC14',
        'x': 'Hip2',
        'y': 'Hip1',
        'z': 'GSC',
        'A': 'AC',
        'B': 'SAO1984',
        'C': 'SAO',
        'D': 'AGK3',
        'E': 'FK4',
        'F': 'ACRS',
        'G': 'LickGas',
        'H': 'Ida93',
        'I': 'Perth70',
        'J': 'COSMOS',
        'K': 'Yale',
        'L': '2MASS',
        'M': 'GSC2.3',
        'N': 'SDSS7',
        'O': 'SSTRC1',
        'P': 'MPOSC3',
        'Q': 'CMC15',
        'R': 'SSTRC4',
        'S': 'URAT1',
        'T': 'URAT2',  # URAT2 missing on ADES web page
        'U': 'Gaia1',
        'V': 'Gaia2',
        'W': 'Gaia3',
        'X': 'Gaia3E',  
        'Y': 'UCAC5',  # UCAC5 mission on ADES web page
        'Z': 'ATLAS2', 
        '0': 'IHW', 
        '1': 'PS1_DR1', 
        '2': 'PS1_DR2', 
        '3': 'Gaia_Int', 
        '4': 'GZ', 
        '5': 'UBSC', 
        '6': 'Gaia_2016', 
    }
    # flip the dictionary to get the catalog codes
    ades_catalog_map = {v: k for k, v in ades_catalog_map.items()}
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

    if verbose:
        print("Applying Eggl et al. (2020) debiasing scheme to the observations.")
    df_groups = obs_df.groupby('astCat')
    for cat, group in df_groups:
        cat_code = ades_catalog_map[cat]
        if cat_code in unbiased_catalogs:
            obs_df.loc[group.index, 'biasRA'] = 0.0
            obs_df.loc[group.index, 'biasDec'] = 0.0
            unbias_counter += len(group)
            continue
        if cat_code not in biased_catalogs:
            print(f"\tUnknown star catalog: {cat}")
            obs_df.loc[group.index, 'biasRA'] = 0.0
            obs_df.loc[group.index, 'biasDec'] = 0.0
            continue
        ra_vals = group['ra'].values
        dec_vals = group['dec'].values
        epoch_vals = Time(group['obsTimeMJD'].to_list(), format='mjd', scale='utc').tt.jd
        ra_cos_dec_biases, dec_biases = debias_obs(ra_vals, dec_vals, epoch_vals,
                                                    cat_code, biasdf, nside)
        obs_df.loc[group.index, 'biasRA'] = ra_cos_dec_biases
        obs_df.loc[group.index, 'biasDec'] = dec_biases
        debias_counter += len(group)
    if verbose:
        no_bias_counter = len(obs_df)-unbias_counter-debias_counter
        print(f"\tNo debiasing needed for {unbias_counter} observations. Debiased {debias_counter}",
                f"observations. No biasing information for {no_bias_counter} observations.")
    return obs_df

def apply_weighting_scheme(obs_df, verbose):
    """
    Apply the weights from Vereš et al. (2017) to the observation array.
    (paper from https://doi.org/10.1016/j.icarus.2017.05.021)

    Parameters
    ----------
    obs_df : pandas DataFrame
        Optical observation data for the given body
    verbose : bool
        Flag to print out information about the observations

    Returns
    -------
    obs_df : pandas DataFrame
        Optical observation data for the given body with weights applied

    Raises
    ------
    ValueError
        If the observation type is unknown
    ValueError
        If the weight is zero after applying the weighting scheme
    """
    # weighting scheme from vereš et al. 2017, obswgt_veres2017.inp
    # valid_modes = ['CCD', 'CMO', 'VID', 'PHO', 'ENC', 'PMT', 'MIC', 'MER', 'TDI', 'OCC', 'UNK']
    if verbose:
        print("Applying Vereš et al. (2017) weighting scheme to the observations.")
    # all non-nan correlations are used, others are 0
    obs_df['sigCorr'] = 0.0
    subset_corr = obs_df.query("rmsCorr == rmsCorr")
    if verbose:
        print(f"\tUsing {len(subset_corr)} observations with provided RA/Dec correlation values.")
    obs_df.loc[subset_corr.index, 'sigCorr'] = subset_corr['rmsCorr'].values
    cols = ['sigRA', 'sigDec']
    mode_grp = obs_df.groupby('mode')
    for mode, group in mode_grp:
        if mode in {'PHO', 'PHA', 'NOR', 'UNK'}:
            subset_1890 = group.query("obsTimeMJD <= 11368.0") # until 1890-01-01
            subset_1950 = group.query("11368.0 < obsTimeMJD <= 33282.0") # until 1950-01-01
            subset_present = group.query("33282.0 < obsTimeMJD") # 1950-01-01 to present
            obs_df.loc[subset_1890.index, cols] = 10.0
            obs_df.loc[subset_1950.index, cols] = 5.0
            obs_df.loc[subset_present.index, cols] = 2.5
        elif mode in {'CCD', 'VID', 'CMO'}: # CCD, Video, CMOS
            obs_df.loc[group.index, cols] = 1.0
            subset_unk_cat = group.query("astCat == 'UNK'") # unknown star catalog
            obs_df.loc[subset_unk_cat.index, cols] = 1.5
            # if rmsRA and rmsDec are not nan, use them
            subset_rms = group.query("rmsRA == rmsRA and rmsDec == rmsDec")
            obs_df.loc[subset_rms.index, cols] = subset_rms[['rmsRA', 'rmsDec']].values
            obs_df.loc[subset_rms.index, 'sigRA'] /= obs_df.loc[subset_rms.index, 'cosDec']
            if verbose:
                print(f"\tUsing {len(subset_rms)} {mode} observations with provided "
                        "RA/Dec RMS values.")
        elif mode in {'OCC'}: # Occultations
            obs_df.loc[group.index, cols] = 0.05
            # if rmsRA and rmsDec are not nan, use them
            subset_occ_rms = group.query("rmsRA == rmsRA and rmsDec == rmsDec")
            obs_df.loc[subset_occ_rms.index, cols] = subset_occ_rms[['rmsRA', 'rmsDec']].values
            obs_df.loc[subset_occ_rms.index, 'sigRA'] /= obs_df.loc[subset_occ_rms.index, 'cosDec']
            if verbose:
                print(f"\tUsing {len(subset_occ_rms)} {mode} observations with provided "
                        "RA/Dec RMS values.")
        elif mode in {'PMT'}: # Hipparcos
            obs_df.loc[group.index, cols] = 0.2
        elif mode in {'MER'}: # Transit circle
            obs_df.loc[group.index, cols] = 0.5
        elif mode in {'ENC'}: # Encoder
            obs_df.loc[group.index, cols] = 0.75
        elif mode in {'MIC'}: # Micrometer
            obs_df.loc[group.index, cols] = 2.0
        elif mode in {'TDI'}: # TDI
            obs_df.loc[group.index, cols] = 1.5
        else:
            raise ValueError(f"Unknown mode: {mode}")
    # check that none of the weights are zero or nan
    if (obs_df[cols] <= 0).any().any():
        raise ValueError("Some weights are zero after applying the weighting scheme.")
    if (obs_df[cols].isna()).any().any():
        raise ValueError("Some weights are NaN after applying the weighting scheme.")
    return obs_df

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

def get_mpc_optical_obs_array(body_id, optical_obs_file=None, t_min_tdb=None,
                            t_max_tdb=None, debias_lowres=True, deweight=True,
                            eliminate=False, num_obs_per_night=5, verbose=False):
    """
    Assemble the optical observations for a given body into an array for an orbit fit.

    Parameters
    ----------
    body_id : str/int
        Target id, numbers are interpreted as asteroids,
        append 'P' for comets, start with comet type and a '/' for comet designations
    optical_obs_file : str, optional
        Filepath to the optical observations file, by default None
    t_min_tdb : float, optional
        Minimum time (MJD TDB) for observations to be included, by default None
    t_max_tdb : float, optional
        Maximum time (MJD TDB) for observations to be included, by default None
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
    """
    if eliminate and deweight:
        raise ValueError('Cannot deweight and eliminate observations at the same time.')
    if not eliminate and not deweight:
        print("WARNING: No deweighting or elimination scheme applied",
                "for observations during the same night.")
    obs_df = get_optical_obs_data(body_id, optical_obs_file,
                                    t_min_tdb, t_max_tdb, verbose)
    obs_df = apply_debiasing_scheme(obs_df, debias_lowres, verbose)
    obs_df = apply_weighting_scheme(obs_df, verbose)
    if deweight:
        (obs_array_optical, star_catalog_codes,
            observer_codes_optical) = deweight_obs(obs_array_optical, star_catalog_codes,
                                        observer_codes_optical, num_obs_per_night, verbose)
    elif eliminate:
        (obs_array_optical, star_catalog_codes,
            observer_codes_optical) = eliminate_obs(obs_array_optical, star_catalog_codes,
                                        observer_codes_optical, num_obs_per_night, verbose)
    return obs_array_optical, observer_codes_optical
