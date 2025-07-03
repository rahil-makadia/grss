"""Optical observation handling for the GRSS orbit determination code"""
import os
from io import StringIO
from astropy.time import Time
import requests
from astroquery.gaia import Gaia
import healpy as hp
import numpy as np
import pandas as pd

from ..utils import grss_path
from .fit_ades import (
    ades_column_types,
    ades_catalog_map,
    unpack_letters,
    pack_letters,
    prog_codes,
    special_codes,
)

__all__ = [ 'add_psv_obs',
            'add_gaia_obs',
            'get_optical_obs',
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
                            timeout=30)
    if response.ok:
        obs_data = response.json()[0]["XML"]
    else:
        print(f"Error {response.status_code}:")
        print(response.content.decode(encoding='utf-8'))
        raise ValueError("Error getting MPC XML data.")
    return StringIO(obs_data)

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
    return None

def _ades_ast_cat_check(df):
    """
    Check the astCat values in the ADES data frame.

    Parameters
    ----------
    df : pandas DataFrame
        ADES data frame

    Returns
    -------
    df : pandas DataFrame
        ADES data frame with invalid astCat values removed

    Raises
    ------
    ValueError
        If the astCat values are invalid
    """
    # from https://www.minorplanetcenter.net/iau/info/ADESFieldValues.html
    valid_cats = list(ades_catalog_map.keys())
    deprecated_cats = [
        'UCAC3', 'UCAC2', 'UCAC1',
        'USNOB1', 'USNOA2', 'USNOSA2', 'USNOA1', 'USNOSA1',
        'Tyc2', 'Tyc1', 'Hip2', 'Hip1', 'ACT',
        'GSCACT', 'GSC2.3', 'GSC2.2', 'GSC1.2', 'GSC1.1', 'GSC1.0', 'GSC',
        'SDSS8', 'SDSS7', 'CMC15', 'CMC14', 'SSTRC4', 'SSTRC1',
        'MPOSC3', 'PPM', 'AC', 'SAO1984', 'SAO', 'AGK3', 'FK4',
        'ACRS', 'LickGas', 'Ida93', 'Perth70', 'COSMOS',
        'Yale', 'ZZCAT', 'IHW', 'GZ', 'UNK']
    df_cats = df['astCat']
    if not df_cats.isin(valid_cats).all():
        invalid_cats = np.setdiff1d(df_cats.unique(), valid_cats)
        print("\tWARNING: At least one unrecognized astCat in the data. "
                f"Unrecognized values are {invalid_cats}. "
                "Force deleting corresponding observations and setting catalog to UNK...")
        delete_idx = df[df['astCat'].isin(invalid_cats)].index
        df.loc[delete_idx, 'selAst'] = 'd'
        df.loc[delete_idx, 'astCat'] = 'UNK'
    return df

def create_optical_obs_df(body_id, optical_obs_file=None, t_min_tdb=None,
                        t_max_tdb=None, verbose=False):
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
        if not os.path.exists(obs_data):
            raise FileNotFoundError(f"File {obs_data} does not exist.")
    # read in the data and add extra columns if not present
    obs_df = pd.read_xml(obs_data, dtype=ades_column_types)
    obs_times = Time(obs_df['obsTime'].to_list(), format='isot', scale='utc')
    obs_df['obsTimeMJD'] = obs_times.utc.mjd
    obs_df['obsTimeMJDTDB'] = obs_times.tdb.mjd
    if 'deprecated' in obs_df:
        # drop rows with deprecated discovery observations
        obs_df.query("deprecated != 'x' and deprecated != 'X'", inplace=True)
    # add columns if they are not present
    str_cols = ['trx', 'rcv', 'sys', 'selAst']
    for col in str_cols:
        if col not in obs_df:
            obs_df[col] = str(np.nan)
    for col in ades_column_types:
        if col not in obs_df:
            obs_df[col] = np.nan
    # remove any columns that are not in ades_column_types
    obs_df = obs_df[ades_column_types.keys()]
    if verbose:
        source = "MPC" if optical_obs_file is None else "file"
        print(f"Read in {len(obs_df)} observations from the {source}.")
    _ades_mode_check(obs_df)
    obs_df = _ades_ast_cat_check(obs_df)
    # filter the data based on the time range
    obs_df.query(f"{t_min_utc} <= obsTimeMJD <= {t_max_utc}", inplace=True)
    # reindex the data frame
    obs_df.reset_index(drop=True, inplace=True)
    # for all indices with OCC mode compute ra and dec as raStar+deltaRA and decStar+deltaDec
    occ_idx = obs_df.query("mode == 'OCC'").index
    obs_df.loc[occ_idx, 'dec'] = obs_df.loc[occ_idx, 'decStar']
    obs_df.loc[occ_idx, 'dec'] += obs_df.loc[occ_idx, 'deltaDec']/3600
    obs_df['cosDec'] = np.cos(obs_df['dec']*np.pi/180)
    obs_df.loc[occ_idx, 'ra'] = obs_df.loc[occ_idx, 'raStar']
    obs_df.loc[occ_idx, 'ra'] += obs_df.loc[occ_idx, 'deltaRA']/3600/obs_df.loc[occ_idx, 'cosDec']
    if verbose:
        print(f"\tFiltered to {len(obs_df)} observations that satisfy the "
                "time range and accepted observatory constraints.")
    return obs_df

def add_psv_obs(psv_obs_file, obs_df, t_min_tdb=None, t_max_tdb=None, verbose=False):
    """
    Assemble the optical observations for a given body into an array for an orbit fit,
    from the ADES PSV observation file.

    Parameters
    ----------
    psv_obs_file : str
        Path to the ADES PSV observation file
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
        Optical observation data for the given body with the PSV observations added
    """
    if t_min_tdb is None:
        t_min_tdb = -np.inf
    if t_max_tdb is None:
        t_max_tdb = np.inf
    # read file and skip header rows starting with '#' or '!'
    skip_count = 0
    with open(psv_obs_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('#') or line.startswith('!'):
                skip_count += 1
            else:
                break
    psv_df = pd.read_csv(psv_obs_file, sep='|', skipinitialspace=True, skiprows=skip_count)
    psv_df.columns = psv_df.columns.str.strip()
    psv_df = psv_df[[col for col in psv_df.columns if col in ades_column_types]]
    psv_column_types = {col: ades_column_types[col] for col in psv_df.columns}
    psv_df = psv_df.astype(psv_column_types)
    psv_df['permID'] = obs_df.iloc[-1]['permID']
    psv_df['provID'] = obs_df.iloc[-1]['provID']
    # strip every entry in a column of leading and trailing whitespace
    psv_df = psv_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    occ_idx = psv_df.query("mode == 'OCC'").index
    if len(occ_idx) > 0:
        psv_df.loc[occ_idx,'dec'] = psv_df.loc[occ_idx,'decStar']
        psv_df.loc[occ_idx,'dec'] += psv_df.loc[occ_idx,'deltaDec']/3600
        psv_df.loc[occ_idx,'cosDec'] = np.cos(psv_df.loc[occ_idx,'dec']*np.pi/180)
        psv_df.loc[occ_idx,'ra'] = psv_df.loc[occ_idx,'raStar']
        psv_df.loc[occ_idx,'ra'] += psv_df.loc[occ_idx,'deltaRA']/3600/psv_df.loc[occ_idx,'cosDec']
    # if biasRA and biasDec are not present, set them to 0
    if 'biasRA' not in psv_df:
        psv_df['biasRA'] = 0.0
    if 'biasDec' not in psv_df:
        psv_df['biasDec'] = 0.0
    psv_df['cosDec'] = np.cos(psv_df['dec']*np.pi/180)
    psv_df['sigRA'] = psv_df['rmsRA']
    psv_df['sigDec'] = psv_df['rmsDec']
    if 'rmsCorr' not in psv_df:
        psv_df['sigCorr'] = 0.0
    if 'rmsTime' not in psv_df:
        psv_df['sigTime'] = 1.0
    times = Time(psv_df['obsTime'].to_list(), format='isot', scale='utc')
    psv_df['obsTimeMJD'] = times.utc.mjd
    psv_df['obsTimeMJDTDB'] = times.tdb.mjd
    add_counter = 0
    if verbose:
        print(f"Read in {len(psv_df)} observations from the ADES PSV file.")
    for i, row in psv_df.iterrows():
        time = times[i].tdb.mjd
        if time < t_min_tdb or time > t_max_tdb:
            continue
        idx = len(obs_df)
        for field in psv_df.columns:
            if field in obs_df.columns:
                obs_df.loc[idx, field] = row[field]
        add_counter += 1
    if verbose:
        print(f"\tFiltered to {add_counter} observations that satisfy the time range constraints.")
    return obs_df

def _get_gaia_query_results(body_id, release):
    """
    Submit a Gaia archive query for a given body ID from a specific
    Gaia data release.

    Parameters
    ----------
    body_id : str/int
        Target id, numbers are interpreted as asteroids,
        append 'P' for comets, start with comet type and a '/' for comet designations
    release : str
        Gaia data release version database name ('gaiadr3', 'gaiafpr', etc.)

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
        + "vx_gaia_geocentric,vy_gaia_geocentric,vz_gaia_geocentric,"
        + "astrometric_outcome_ccd, astrometric_outcome_transit"
        + f" FROM {release}.{table} {match_str} ORDER BY epoch_utc ASC"
    )
    job = Gaia.launch_job_async(query, dump_to_file=False,background=True)
    res = job.get_results()
    # change the column names to be lowercase
    for col_name in res.colnames:
        res.rename_column(col_name, col_name.lower())
    res.sort('epoch_utc')
    return res

def add_gaia_obs(obs_df, t_min_tdb=None, t_max_tdb=None, gaia_dr='gaiafpr', verbose=False):
    """
    Assemble the optical observations for a given body from Gaia FPR.

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
        Gaia data release version database name, by default 'gaiafpr'
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
    perm_id = obs_df.iloc[-1]['permID']
    prov_id = obs_df.iloc[-1]['provID']
    body_id = perm_id if isinstance(perm_id, str) else prov_id
    res = _get_gaia_query_results(body_id, release=gaia_dr)
    if verbose:
        print(f"Read in {len(res)} Gaia observations from {gaia_dr}")
    sys = 'ICRF_AU'
    ctr = 500
    curr_transit = int(-1e6)
    gaia_add_counter = 0
    for i, data in enumerate(res):
        # # print(data["astrometric_outcome_ccd"], data["astrometric_outcome_transit"])
        # # print(type(data["astrometric_outcome_ccd"]), type(data["astrometric_outcome_transit"]))
        # if data['astrometric_outcome_ccd'] != 1 or data['astrometric_outcome_transit'] != 1:
        #     continue
        if curr_transit != data['transit_id']:
            curr_transit = data['transit_id']
            transit_count = 1
            while i+transit_count < len(res) and res[i+transit_count]['transit_id'] == curr_transit:
                transit_count += 1
        obs_time = Time(data['epoch_utc'] + 55197.0, format='mjd', scale='utc')
        if obs_time.tdb.mjd < t_min_tdb or obs_time.tdb.mjd > t_max_tdb:
            continue
        idx = len(obs_df)
        gaia_add_counter += 1
        obs_df.loc[idx, 'permID'] = perm_id
        obs_df.loc[idx, 'provID'] = prov_id
        obs_df.loc[idx, 'obsTime'] = f'{obs_time.utc.iso}Z'
        obs_df.loc[idx, 'obsTimeMJD'] = data['epoch_utc'] + 55197.0
        obs_df.loc[idx, 'obsTimeMJDTDB'] = obs_time.tdb.mjd
        obs_df.loc[idx, 'mode'] = 'CCD'
        obs_df.loc[idx, 'stn'] = '258'
        obs_df.loc[idx, 'ra'] = data['ra']
        obs_df.loc[idx, 'dec'] = data['dec']
        cosdec = np.cos(data['dec']*np.pi/180)
        obs_df.loc[idx, 'cosDec'] = cosdec
        ra_sig_sys = data['ra_error_systematic']/1000
        ra_sig_rand = data['ra_error_random']/1000
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
        obs_df.loc[idx, 'sigRA'] = ra_sig
        obs_df.loc[idx, 'sigDec'] = dec_sig
        obs_df.loc[idx, 'sigCorr'] = corr
        # obs_df.loc[idx, 'sigTime'] = data['epoch_err']*86400.0
        obs_df.loc[idx, 'biasRA'] = 0.0
        obs_df.loc[idx, 'biasDec'] = 0.0
        obs_df.loc[idx, 'ctr'] = ctr
        obs_df.loc[idx, 'sys'] = sys
        tcb_tdb_fac = 1 - 1.550519768e-8
        obs_df.loc[idx, 'pos1'] = data['x_gaia_geocentric']*tcb_tdb_fac
        obs_df.loc[idx, 'pos2'] = data['y_gaia_geocentric']*tcb_tdb_fac
        obs_df.loc[idx, 'pos3'] = data['z_gaia_geocentric']*tcb_tdb_fac
        # obs_df.loc[idx, 'vel1'] = data['vx_gaia_geocentric']*tcb_tdb_fac
        # obs_df.loc[idx, 'vel2'] = data['vy_gaia_geocentric']*tcb_tdb_fac
        # obs_df.loc[idx, 'vel3'] = data['vz_gaia_geocentric']*tcb_tdb_fac
        # # add some position uncertainty (testing)
        # au2km = 149597870.7
        # pos_sig = 1.0/au2km
        # obs_df.loc[idx, 'posCov11'] = pos_sig**2
        # obs_df.loc[idx, 'posCov12'] = 0.0
        # obs_df.loc[idx, 'posCov13'] = 0.0
        # obs_df.loc[idx, 'posCov22'] = pos_sig**2
        # obs_df.loc[idx, 'posCov23'] = 0.0
        # obs_df.loc[idx, 'posCov33'] = pos_sig**2
    if verbose:
        print(f"\tFiltered to {gaia_add_counter} observations that",
                "satisfy the time range constraints.")
    return obs_df

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
    bias_vals = biasdf.loc[idx_vals, colnames].values
    # time from epoch in Julian years
    dt_jy_vals = (epoch_vals-j2000_jd)/365.25
    # compute biases
    ddec_vals = (bias_vals[:, 1]+dt_jy_vals*bias_vals[:, 3]/1000)
    dra_cos_dec_vals = (bias_vals[:, 0]+dt_jy_vals*bias_vals[:, 2]/1000)
    return dra_cos_dec_vals, ddec_vals

def apply_debiasing_scheme(obs_df, lowres, verbose):
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
    # MPC catalog codes, https://www.minorplanetcenter.net/iau/info/CatalogueCodes.html
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
    debias_lowres_path = grss_path+'/debias/lowres_data'
    debias_hires_path = grss_path+'/debias/hires_data'
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
            if verbose:
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
        print(f"\tNo debiasing needed for {unbias_counter} observations.")
        print(f"\tDebiased {debias_counter} observations.")
        print(f"\tNo bias information for {no_bias_counter} observations.")
    return obs_df

# from ADES-Master/Python/bin/packUtil.py
def get_packed_prog_id(code):
    """
    Get the packed program ID from the ADES data.

    Parameters
    ----------
    code : str
        Full ADES Program code ID

    Returns
    -------
    packed : str
        Packed program ID
    """
    # should be radix-52 with 0-1 first
    if len(code) != 2 or code[0] not in "01":
        raise RuntimeError ("Illegal packed prog ID " + code + " in xml")
    try:
        index = unpack_letters[code[1]]
        index += 62*unpack_letters[code[0]]
        packed = prog_codes[index] if index <= 93 else ' '
    except Exception as exc:
        raise RuntimeError ("Illegal packed prog ID " + code + " in xml") from exc
    return packed

# from ADES-Master/Python/bin/packUtil.py
def get_unpacked_prog_id(code):
    """
    Get the unpacked program ID for the ADES data.

    Parameters
    ----------
    code : str
        Packed program code ID

    Returns
    -------
    unpacked : str
        Full ADES Program code ID
    """
    # should be radix-52 with 0-1 first
    try:
        index = prog_codes.index(code)
    except Exception as exc:
        raise RuntimeError("Illegal program code " + code) from exc
    first = '1' if index > 61 else '0'
    second = pack_letters[index%62]
    return first + second

def apply_station_weight_rules(group, obs_df, cols, verbose):
    """
    Apply the station-specific weight rules to the observation data.

    Parameters
    ----------
    group : pandas DataFrame
        Group of CCD observations in the observation data
    obs_df : pandas DataFrame
        Optical observation data for the given body
    cols : list
        Columns to apply the weight rules to (sigRA, sigDec)
    verbose : bool
        Flag to print out information about the observations
    """
    # sourcery skip: low-code-quality
    default_weight_counter = 0
    all_times = obs_df.loc[group.index, 'obsTimeMJD'].values
    all_stns = obs_df.loc[group.index, 'stn'].values
    all_cats = obs_df.loc[group.index, 'astCat'].values
    all_progs = obs_df.loc[group.index, 'prog'].values
    ccd_weight_arr = np.ones((len(group), 2))*np.nan
    for i in range(len(group)):
        time = all_times[i]
        stn = all_stns[i]
        cat = all_cats[i]
        prog = all_progs[i]
        # star catalog and program code need to be switched back to old 80-column format
        cat = ades_catalog_map[cat]
        if isinstance(prog, str):
            prog = get_packed_prog_id(prog)
        if stn in {'F51', 'F52'}:
            ccd_weight_arr[i,:] = 0.2
            if prog in {'Z'}:
                ccd_weight_arr[i,:] = 1.0
        elif stn in {'G96'}:
            ccd_weight_arr[i,:] = 0.5
            if prog in {'Z'}:
                ccd_weight_arr[i,:] = 1.0
        elif stn in {'703'}:
            # until 2014-01-01
            ccd_weight_arr[i,:] = 1.0 if time <= 56658.0 else 0.8
            if prog in {'Z'}:
                ccd_weight_arr[i,:] = 1.0
        elif stn in {'E12'}:
            ccd_weight_arr[i,:] = 0.75
        elif stn in {'704'}:
            ccd_weight_arr[i,:] = 1.0
        elif stn in {'691', '291'}:
            # until 2003-01-01
            ccd_weight_arr[i,:] = 0.6 if time <= 56240.0 else 0.5
        elif stn in {'644'}:
            # until 2003-09-01
            ccd_weight_arr[i,:] = 0.6 if time <= 52883.0 else 0.4
        elif stn in {'699'}:
            ccd_weight_arr[i,:] = 0.8
        elif stn in {'G45'}:
            ccd_weight_arr[i,:] = 0.6
        elif stn in {'D29'}:
            ccd_weight_arr[i,:] = 0.75
        elif stn in {'T05', 'T08', 'T07', 'M22', 'W68'}:
            ccd_weight_arr[i,:] = 0.8
        elif stn in {'568'}:
            ccd_weight_arr[i,:] = 0.5
            if prog in {'_'}:
                if cat in {'t', 'q'}:
                    ccd_weight_arr[i,:] = 0.25
                elif cat in {'U', 'V', 'W', 'X'}:
                    ccd_weight_arr[i,:] = 0.15
            elif prog in {'^'} and cat in {'U', 'V', 'W', 'X'}:
                ccd_weight_arr[i,:] = 0.2
            elif prog in {'2'}:
                if cat in {'o', 's'}:
                    ccd_weight_arr[i,:] = 0.5
                elif cat in {'t'}:
                    ccd_weight_arr[i,:] = 0.2
                elif cat in {'U', 'V', 'W', 'X'}:
                    ccd_weight_arr[i,:] = 0.1
            elif prog in {'0'}:
                ccd_weight_arr[i,:] = 1.0
        elif stn in {'T09'}:
            ccd_weight_arr[i,:] = 0.5
            if prog in {'0'}:
                ccd_weight_arr[i,:] = 0.1
            elif prog in {'1', '4'}:
                ccd_weight_arr[i,:] = 1.0
            elif prog in {'2', '3'}:
                ccd_weight_arr[i,:] = 0.5
            if time <= 59858.0: # until 2022-10-06
                ccd_weight_arr[i,:] = 0.1
        elif stn in {'T10', 'T11', 'T13', 'T16', 'T17'}:
            ccd_weight_arr[i,:] = 0.5
        elif stn in {'T12'}:
            ccd_weight_arr[i,:] = 0.5
            if prog in {'0'}:
                ccd_weight_arr[i,:] = 0.1
            elif prog in {'1'}:
                ccd_weight_arr[i,:] = 0.2
            elif prog in {'2', '3'}:
                ccd_weight_arr[i,:] = 0.5
            if time <= 59858.0: # until 2022-10-06
                ccd_weight_arr[i,:] = 0.1
        elif stn in {'T14'}:
            ccd_weight_arr[i,:] = 0.5
            if prog in {'0', '7'}:
                ccd_weight_arr[i,:] = 0.1
            elif prog in {'1', '2', '4', '6', '8'}:
                ccd_weight_arr[i,:] = 0.5
            elif prog in {'3'}:
                ccd_weight_arr[i,:] = 0.2
            elif prog in {'5'}:
                ccd_weight_arr[i,:] = 1.0
            if time <= 59858.0: # until 2022-10-06
                ccd_weight_arr[i,:] = 0.1
        elif stn in {'T15'}:
            ccd_weight_arr[i,:] = 0.5
            if prog in {'2'}:
                ccd_weight_arr[i,:] = 1.0
        elif stn in {'H01'} and cat in {'L', 't', 'U', 'V', 'W', 'X'}:
                ccd_weight_arr[i,:] = 0.3
        elif stn in {'673'} and prog in {'1'}:
                ccd_weight_arr[i,:] = 0.3
        elif stn in {'645'}:
            ccd_weight_arr[i,:] = 0.3
        elif stn in {'689'}:
            ccd_weight_arr[i,:] = 0.5
        elif stn in {'J04', 'Z84'}:
            ccd_weight_arr[i,:] = 1.0
            if cat in {'U', 'V', 'W', 'X', 't', 'L', 'q', 'r', 'u', 'e'}:
                ccd_weight_arr[i,:] = 0.5
            if stn in {'J04'}:
                if prog in {'2', '$'} and cat in {'t', 'q', 'U', 'V', 'W', 'X'}:
                    ccd_weight_arr[i,:] = 0.3
                elif prog in {'#'}:
                    ccd_weight_arr[i,:] = 1.0
            elif stn in {'Z84'}:
                if prog in {'1', '7'} and cat in {'U', 'V', 'W', 'X'}:
                    ccd_weight_arr[i,:] = 0.3
        elif stn in {'608'}:
            ccd_weight_arr[i,:] = 0.6
        elif stn in {'W84'}:
            ccd_weight_arr[i,:] = 0.5
            if prog in {'&'}:
                ccd_weight_arr[i,:] = 1.0
        elif stn in {'950'}:
            ccd_weight_arr[i,:] = 1.0
            if cat in {'U', 'V', 'W', 'X', 't', 'L', 'q', 'r', 'u', 'e'}:
                ccd_weight_arr[i,:] = 0.5
            if prog in {'@'}:
                ccd_weight_arr[i,:] = 1.0
        elif stn in {'E10'}:
            ccd_weight_arr[i,:] = 0.4
            if prog in {'(', '.'} and cat in {'U', 'V', 'W', 'X'}:
                ccd_weight_arr[i,:] = 0.2
        elif stn in {'F65'}:
            ccd_weight_arr[i,:] = 0.4
            if prog in {'8'} and cat in {'U', 'V', 'W', 'X'}:
                ccd_weight_arr[i,:] = 0.2
            elif prog in {'9'} and cat in {'U', 'V', 'W', 'X'}:
                ccd_weight_arr[i,:] = 0.2
        elif stn in {'K91', 'K92', 'K93', 'Q63', 'Q64', 'V39',
                    'W85', 'W86', 'W87', 'Z24', 'Z31'}:
            ccd_weight_arr[i,:] = 0.4
            if prog in {'0', '3'}:
                ccd_weight_arr[i,:] = 0.3
        elif stn in {'V37'}:
            ccd_weight_arr[i,:] = 0.4
            if prog in {'0', '6'}:
                ccd_weight_arr[i,:] = 0.3
        elif stn in {'Y28'} and cat in {'t', 'U', 'V', 'W', 'X'}:
            ccd_weight_arr[i,:] = 0.3
        elif stn in {'309'}:
            ccd_weight_arr[i,:] = 1.0
            if cat in {'U', 'V', 'W', 'X'}:
                ccd_weight_arr[i,:] = 0.2
                if prog in {'&', '%'}:
                    ccd_weight_arr[i,:] = 0.1
            elif cat in {'t', 'q'}:
                ccd_weight_arr[i,:] = 0.3
                if prog in {'&', '%'}:
                    ccd_weight_arr[i,:] = 0.3
        elif stn in {'G83'} and prog in {'2'} and cat in {'U', 'V', 'W', 'X'}:
                ccd_weight_arr[i,:] = 0.2
        elif stn in {'C65'} and prog in {'3'} and cat in {'U', 'V', 'W', 'X'}:
                ccd_weight_arr[i,:] = 0.2
        elif stn in {'094'} and prog in {'9'} and cat in {'U', 'V', 'W', 'X'}:
                ccd_weight_arr[i,:] = 0.5
        elif stn in {'Z18'} and prog in {'1'} and cat in {'U', 'V', 'W', 'X'}:
                ccd_weight_arr[i,:] = 0.1
        elif stn in {'181'} and prog in {'1'} and cat in {'U', 'V', 'W', 'X'}:
                ccd_weight_arr[i,:] = 0.5
        elif stn in {'D20'} and prog in {'1', '4'} and cat in {'U', 'V', 'W', 'X'}:
                ccd_weight_arr[i,:] = 0.5
        elif stn in {'G37'} and prog in {'5'} and cat in {'U', 'V', 'W', 'X'}:
                ccd_weight_arr[i,:] = 0.2
        elif stn in {'N50'} and prog in {'2'} and cat in {'U', 'V', 'W', 'X'}:
                ccd_weight_arr[i,:] = 0.3
        elif stn in {'I41'} and prog in {'Z'}:
            ccd_weight_arr[i,:] = 1.0
        elif stn in {'W57'} and prog in {'0', '1'}:
                ccd_weight_arr[i,:] = 0.4
        elif stn in {'M28'}:
            ccd_weight_arr[i,:] = 3.0
        elif stn in {'247', '270'}:
            ccd_weight_arr[i,:] = 2.0
        else:
            ccd_weight_arr[i,:] = 1.0
            default_weight_counter += 1
    # set the values in obs_df to values from ccd_weight_arr
    obs_df.loc[group.index, cols] = ccd_weight_arr
    if verbose:
        print(f"\tUsing {len(group)-default_weight_counter} CCD",
                "observations with station-specific weight rules.")
    return None

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
    # sourcery skip: low-code-quality
    # weighting scheme from vereš et al. 2017, obswgt_veres2017.inp
    if verbose:
        print("Applying Vereš et al. (2017) weighting scheme to the observations.")
    # all correlations are 0 by default
    non_radar = obs_df.query("mode != 'RAD'")
    obs_df.loc[non_radar.index, 'sigCorr'] = 0.0
    cols = ['sigRA', 'sigDec']
    mode_grp = obs_df.groupby('mode')
    default_time_uncert = 1.0
    no_time_uncert = (
        special_codes["gaia"]
        | special_codes["occultation"]
        | special_codes["spacecraft"]
    )
    non_radar_occ_sc = obs_df.query(
        "mode != 'RAD' and stn not in @no_time_uncert and sigTime != sigTime"
    )
    obs_df.loc[non_radar_occ_sc.index, 'sigTime'] = default_time_uncert
    for mode, group in mode_grp:
        if mode in {'PHO', 'PHA', 'NOR', 'UNK'}:
            date_1890 = group.query("obsTimeMJD <= 11368.0") # until 1890-01-01
            date_1950 = group.query("11368.0 < obsTimeMJD <= 33282.0") # until 1950-01-01
            date_present = group.query("33282.0 < obsTimeMJD") # 1950-01-01 to present
            obs_df.loc[date_1890.index, cols] = 10.0
            obs_df.loc[date_1950.index, cols] = 5.0
            obs_df.loc[date_present.index, cols] = 2.5
        elif mode in {'CCD', 'VID', 'CMO'}: # CCD, Video, CMOS
            obs_df.loc[group.index, cols] = 1.0
            unknown_cat = group.query("astCat == 'UNK'") # unknown star catalog
            obs_df.loc[unknown_cat.index, cols] = 1.5
            if mode in {'CCD'}:
                apply_station_weight_rules(group, obs_df, cols, verbose)
            # # if rmsRA and rmsDec are not nan, use them
            # ccd_rms = group.query("rmsRA == rmsRA and rmsDec == rmsDec and rmsCorr == rmsCorr")
            # obs_df.loc[ccd_rms.index, cols] = ccd_rms[['rmsRA', 'rmsDec']].values
            # obs_df.loc[ccd_rms.index, 'sigCorr'] = ccd_rms['rmsCorr'].values
            # if verbose:
            #     print(f"\tUsing {len(ccd_rms)} {mode} observations with provided "
            #             "RA/Dec RMS and correlation values.")
        elif mode in {'OCC'}: # Occultations
            obs_df.loc[group.index, cols] = 0.01
            # if rmsRA and rmsDec are not nan, use them
            occ_rms = group.query("rmsRA == rmsRA and rmsDec == rmsDec and rmsCorr == rmsCorr")
            obs_df.loc[occ_rms.index, cols] = occ_rms[['rmsRA', 'rmsDec']].values
            obs_df.loc[occ_rms.index, 'sigCorr'] = occ_rms['rmsCorr'].values
            if verbose:
                print(f"\tUsing {len(occ_rms)} {mode} observations with provided "
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

def deweight_obs(obs_df, eff_obs_per_night, verbose):
    """
    Deweight observations that occur on the same night at the same observatory.

    Parameters
    ----------
    obs_df : pandas DataFrame
        Optical observation data for the given body
    eff_obs_per_night : int
        Effective number of observations per night
    verbose : bool
        Flag to print out information about the observations

    Returns
    -------
    obs_df : pandas DataFrame
        Optical observation data for the given body with deweighting applied
    """
    night_count = 1
    deweight_count = 0
    if verbose:
        print(f"Applying sqrt(N/{eff_obs_per_night}) deweighting scheme.")
    times = obs_df['obsTimeMJD'].values
    stations = obs_df['stn'].values
    weights = obs_df[['sigRA', 'sigDec']].values
    batch_first_time = times[0]
    for i in range(1, len(obs_df)):
        night_match = times[i] - batch_first_time < 8/24
        curr_observatory = stations[i]
        prev_observatory = stations[i-1]
        observatory_match = curr_observatory == prev_observatory
        if night_match and observatory_match:
            night_count += 1
        else:
            if night_count > eff_obs_per_night:
                deweight_count += night_count
                factor = night_count**0.5/eff_obs_per_night**0.5
                weights[i-night_count:i] *= factor
            night_count = 1
            batch_first_time = times[i]
    # edge case where last observation is part of a batch that needs deweighting
    if night_count > eff_obs_per_night:
        deweight_count += night_count
        factor = night_count**0.5/eff_obs_per_night**0.5
        weights[i+1-night_count:i+1] *= factor
    obs_df[['sigRA', 'sigDec']] = weights
    if verbose:
        print(f"\tDeweighted {deweight_count} observations.")
    return obs_df

def eliminate_obs(obs_df, max_obs_per_night, verbose):
    """
    Limit the number of observations that occur on the same night at the same observatory.

    Parameters
    ----------
    obs_df : pandas DataFrame
        Optical observation data for the given body
    max_obs_per_night : int
        Maximum number of observations per night
    verbose : bool
        Flag to print out information about the observations

    Returns
    -------
    obs_df : pandas DataFrame
        Optical observation data for the given body with elimination applied
    """
    night_count = 1
    eliminate_count = 0
    idx_to_drop = []
    if verbose:
        print(f"Applying {max_obs_per_night}-observation per night elimination scheme.")
    times = obs_df['obsTimeMJD'].values
    stations = obs_df['stn'].values
    batch_first_time = times[0]
    for i in range(1, len(obs_df)):
        night_match = times[i] - batch_first_time < 8/24
        curr_observatory = stations[i]
        prev_observatory = stations[i-1]
        observatory_match = curr_observatory == prev_observatory
        if night_match and observatory_match:
            night_count += 1
            if night_count > max_obs_per_night:
                eliminate_count += 1
                idx_to_drop.append(i)
        else:
            night_count = 1
            batch_first_time = times[i]
    obs_df.drop(idx_to_drop, inplace=True)
    obs_df.reset_index(drop=True, inplace=True)
    if verbose:
        print(f"\tEliminated {eliminate_count} observations.")
    return obs_df

def get_optical_obs(body_id, optical_obs_file=None, t_min_tdb=None,
                    t_max_tdb=None, debias_lowres=True, deweight=True,
                    eliminate=False, num_obs_per_night=5, verbose=False,
                    accept_weights=False):
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
    accept_weights : bool, optional
        Flag to accept the weights from the input file, by default False

    Returns
    -------
    obs_df : pandas DataFrame
        Optical observation data for the given body with the desired biasing and weighting

    Raises
    ------
    ValueError
        If deweight and eliminate are both True.
    ValueError
        If deweight and eliminate are both False.
    """
    if eliminate and deweight:
        raise ValueError('Cannot deweight and eliminate observations at the same time.')
    if not eliminate and not deweight and verbose:
        print("WARNING: No deweighting or elimination scheme applied",
                "for observations during the same night.")
    obs_df = create_optical_obs_df(body_id, optical_obs_file,
                                    t_min_tdb, t_max_tdb, verbose)
    if debias_lowres is not None:
        obs_df = apply_debiasing_scheme(obs_df, debias_lowres, verbose)
    else:
        if verbose:
            print("WARNING: No debiasing scheme applied to the observations.",
                    "Setting biases to zero.")
        opt_idx = obs_df.query("mode != 'RAD'").index
        obs_df.loc[opt_idx, ['biasRA', 'biasDec']] = 0.0
    if not accept_weights:
        obs_df = apply_weighting_scheme(obs_df, verbose)
    if deweight:
        obs_df = deweight_obs(obs_df, num_obs_per_night, verbose)
    elif eliminate:
        obs_df = eliminate_obs(obs_df, num_obs_per_night, verbose)
    return obs_df
