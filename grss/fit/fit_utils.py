"""Utilities for the GRSS orbit determination code"""
import json
import requests
import numpy as np
from numba import jit

from ..utils import grss_path
from .fit_ades import special_codes

__all__ = [ 'mjd2et',
            'et2mjd',
            'get_codes_dict',
            'get_observer_info',
            'get_sbdb_info',
            'get_similarity_stats',
]

@jit('float64(float64)', nopython=True, cache=True)
def mjd2et(mjd):
    """
    Converts Modified Julian Date to JPL NAIF SPICE Ephemeris time.

    Parameters
    ----------
    mjd : float
        Modified Julian Date

    Returns
    -------
    et : float
        Ephemeris time (ephemeris seconds beyond epoch J2000)
    """
    return (mjd+2400000.5-2451545.0)*86400

@jit('float64(float64)', nopython=True, cache=True)
def et2mjd(ephem_time):
    """
    Converts JPL NAIF SPICE Ephemeris time to Modified Julian Date.

    Parameters
    ----------
    ephem_time : float
        Ephemeris time (ephemeris seconds beyond epoch J2000)

    Returns
    -------
    mjd : float
        Modified Julian Date
    """
    return ephem_time/86400-2400000.5+2451545.0


@jit('Tuple((float64, float64, float64))(float64, float64, float64)',
        nopython=True, cache=True)
def parallax_to_lat_lon_alt(lon, rho_cos_lat, rho_sin_lat):
    """
    Convert parallax constants to geocentric longitude, latitude, and distance.

    Parameters
    ----------
    lon : float
        degrees, longitude of observer
    rho_cos_lat : float
        rho cos phi', where rho is the distance from the geocenter to the observer
        and phi' is the geocentric latitude of the observer
    rho_sin_lat : float
        rho sin phi', where rho is the distance from the geocenter to the observer
        and phi' is the geocentric latitude of the observer

    Returns
    -------
    lon : float
        radians, longitude of observer
    lat : float
        radians, geocentric latitude of observer
    rho : float
        distance from the geocenter to the observer, m
    """
    # WGS84 ellipsoid parameters
    r_eq = 6378.137/1.495978707e8 # equatorial radius in AU
    lat = np.arctan2(rho_sin_lat, rho_cos_lat)
    rho = np.sqrt(rho_cos_lat**2 + rho_sin_lat**2)*r_eq
    return lon*np.pi/180, lat, rho

def get_mpc_observatory_info():
    """
    Read in the MPC observatory codes and return them as a dictionary.

    Returns
    -------
    codes_dict : dict
        dictionary of observatory codes and their corresponding longitude,
        rho*cos(latitude), and rho*sin(latitude)
    """
    with open(f'{grss_path}/fit/codes.json', 'r', encoding='utf-8') as f:
        mpc_info_dict = json.load(f)
    return mpc_info_dict

def get_radar_codes_dict():
    # sourcery skip: inline-immediately-returned-variable
    """
    Creates a dictionary of radar codes and their corresponding longitude,
    geocentric latitude, and distance from the geocenter

    Returns
    -------
    radar_codes_dict : dict
        dictionary of radar codes and their corresponding longitude,
        geocentric latitude, and distance from the geocenter
    """
    # from JPL radar API
    m2au = 1.0/1.495978707e11
    radar_site_loc = {
        '-13': (4.244737449544063, 0.6120174701461084, 6372125.10051922*m2au),
        '-1': (5.1181310262659565, 0.31816602936555305, 6376488.108593255*m2au),
        '-2': (5.035481521533785, 0.7405709918867184, 6368491.961110464*m2au),
        '-14': (4.2430780044539, 0.6151299903972761, 6371993.266933192*m2au),
        '-47': (2.6101423188745216, -0.5261379198959515, 6372960.260589986*m2au),
        '-9': (4.889717249348727, 0.6675169909800179, 6370740.389143447*m2au),
        '-43': (2.600213638178763, -0.6147210003517548, 6371689.000848973*m2au),
        '-73': (0.33555473306040223, 1.2123123239363724, 6359490.218520948*m2au),
        '-38': (0.5792216160079022, 0.7853411225947086, 6367487.942395174*m2au),
        '-36': (2.6001661111179013, -0.6145934819355463, 6371688.2123778965*m2au),
        '-35': (2.6002169281244027, -0.6146055679274739, 6371697.3584917635*m2au),
        '-25': (4.24332540487537, 0.6135924747227121, 6371982.53410041*m2au),
    }
    return radar_site_loc

def get_codes_dict():
    """
    Creates a dictionary of MPC codes and their corresponding longitude,
    geocentric latitude, and distance from the geocenter

    Returns
    -------
    codes_dict : dict
        dictionary of MPC codes and their corresponding longitude,
        geocentric latitude, and distance from the geocenter
    """
    codes_dict = {}
    mpc_info_dict = get_mpc_observatory_info()
    for code, info in mpc_info_dict.items():
        lon, lat, rho = parallax_to_lat_lon_alt(info[0], info[1], info[2])
        codes_dict[code] = (lon, lat, rho)
    radar_codes_dict = get_radar_codes_dict()
    for code, info in radar_codes_dict.items():
        codes_dict[code] = info
    return codes_dict

def get_observer_info(obs_df):
    """
    Creates a list of observer information for each observer code for
    passing into a libgrss PropSimulation object

    Parameters
    ----------
    obs_df : pandas DataFrame
        Observation data conatining observer information

    Returns
    -------
    observer_info : list
        Information about each observer (Spice ID for parent body, longitude,
        geocentric latitude, and altitude). Shape of each entry changes on type
        of observation: 4 (optical), 8 (radar), or 9 (doppler).
    """
    stn_codes_dict = get_codes_dict()
    observer_info = []
    gaia = special_codes['gaia']
    occultation = special_codes['occultation']
    spacecraft = special_codes['spacecraft']
    roving = special_codes['roving']
    conv_to_au = {
        'ICRF_AU': (1,1),
        'ICRF_KM': (1/1.495978707e8,86400/1.495978707e8),
    }
    m2au = 1.0/1.495978707e11
    fields = ['stn', 'mode', 'sys', 'pos1', 'pos2', 'pos3',
              'vel1', 'vel2', 'vel3', 'trx', 'rcv', 'com', 'frq', 'doppler']
    for obs_info in zip(*[obs_df[field] for field in fields]):
        stn, mode, sys, pos1, pos2, pos3, vel1, vel2, vel3, trx, rcv, com, freq, doppler = obs_info
        info_list = []
        if stn in gaia|occultation|spacecraft:
            c = conv_to_au[sys]
            c_pos = c[0]
            if stn in gaia:
                c_vel = c[1]
                info_list.extend( (500, pos1*c_pos, pos2*c_pos, pos3*c_pos,
                                        vel1*c_vel, vel2*c_vel, vel3*c_vel) )
            else:
                info_list.extend((500, pos1*c_pos, pos2*c_pos, pos3*c_pos, 0, 0, 0))
        elif mode == 'RAD':
            rx_lon, rx_lat, rx_rho = stn_codes_dict[rcv]
            tx_lon, tx_lat, tx_rho = stn_codes_dict[trx]
            info_list.extend((
                399, rx_lon, rx_lat, rx_rho,
                399, tx_lon, tx_lat, tx_rho,
                com
            ))
            if np.isfinite(doppler):
                info_list.append(freq)
        elif stn in roving:
            lon = pos1*np.pi/180
            if sys == 'WGS84':
                lat = pos2*np.pi/180
                rho = (pos3+6378137.0)*m2au
            elif sys == 'ITRF':
                lat = np.arctan2(pos2, pos3)
                rho = np.sqrt(pos2**2 + pos3**2)*1.0e3*m2au
            else:
                raise ValueError(f"Invalid system {sys}")
            info_list.extend((399, lon, lat, rho))
        elif mode.startswith('SIM'):
            lon = pos1*np.pi/180
            lat = np.arctan2(pos2, pos3)
            rho = np.sqrt(pos2**2 + pos3**2)*1.0e3*m2au
            if mode in {'SIM_RAD_DEL', 'SIM_RAD_DOP'}:
                info_list.extend((
                    399, lon, lat, rho,
                    399, lon, lat, rho,
                    com
                ))
                if np.isfinite(doppler):
                    info_list.append(freq)
            else:
                if sys == 'WGS84':
                    lat = pos2*np.pi/180
                    rho = (pos3+6378137.0)*m2au
                elif sys == 'ITRF':
                    lat = np.arctan2(pos2, pos3)
                    rho = np.sqrt(pos2**2 + pos3**2)*1.0e3*m2au
                else:
                    raise ValueError(f"Invalid system {sys}")
                info_list.extend((399, lon, lat, rho))
        else:
            lon, lat, rho = stn_codes_dict[stn]
            info_list.extend((399, lon, lat, rho))
        observer_info.append(info_list)
    return observer_info

def get_sbdb_raw_data(tdes):
    """
    Get json data from JPL SBDB entry for desired small body

    Parameters
    ----------
    tdes : str
        IMPORTANT: must be the designation of the small body, not the name.

    Returns
    -------
    raw_data : dict
        JSON output of small body information query from SBDB
    """
    url = f"https://ssd-api.jpl.nasa.gov/sbdb.api?sstr={tdes}&cov=mat&phys-par=true&full-prec=true"
    req = requests.request("GET", url, timeout=30)
    return json.loads(req.text)

def get_sbdb_elems(tdes, cov_elems=True):
    """
    Get a set of desired elements for small body from SBDB API

    Parameters
    ----------
    tdes : str
        IMPORTANT: must be the designation of the small body, not the name.
    cov_elems : bool, optional
        Boolean for whether to extract element set corresponding to the covariance
        information (since the nominal set on the webpage might have a different epoch
        than the full covariance info), by default True

    Returns
    -------
    elements : dict
        Dictionary containing desired cometary elements for the small body
    """
    response = requests.get("https://data.minorplanetcenter.net/api/query-identifier",
                            data=tdes, timeout=30)
    if response.ok:
        tdes = response.json()['unpacked_primary_provisional_designation']
    else:
        print("get_radar_obs_array: ERROR. ", response.status_code, response.content)
        raise ValueError("Failed to get JPL orbit data")
    raw_data = get_sbdb_raw_data(tdes)
    # check that covariance exists
    cov_exists = raw_data['orbit']['covariance'] is not None
    if cov_elems and cov_exists:
        # epoch of orbital elements at reference time [JD -> MJD]
        epoch_mjd = float(raw_data['orbit']['covariance']['epoch']) - 2400000.5
        # cometary elements at epoch_mjd
        if epoch_mjd == float(raw_data['orbit']['epoch']) - 2400000.5:
            elem = raw_data['orbit']['elements']
        else:
            elem = raw_data['orbit']['covariance']['elements']
    else:
        # epoch of orbital elements at reference time [JD -> MJD]
        epoch_mjd = float(raw_data['orbit']['epoch']) - 2400000.5
        # cometary elements at epoch_mjd
        elem = raw_data['orbit']['elements']
    hdr = []
    val = []
    for ele in elem:
        if ele['value'] is not None:
            hdr.append(ele['name'])
            val.append(float(ele['value']))
    full_elements_dict = dict(zip(hdr, val))
    # cometary elements
    # eccentricity, perihelion distance, time of periapse passage (JD),
    # longitude of the ascending node, argument of perihelion, inclination
    keys = ['e', 'q', 'tp', 'om', 'w', 'i']
    elements = {'t': epoch_mjd}
    # add every element key to elements dictionary
    for key in keys:
        elements[key] = full_elements_dict[key]
        if key == 'tp':
            elements[key] = full_elements_dict[key] - 2400000.5
        if key in {'om', 'w', 'i'}:
            elements[key] = full_elements_dict[key]*np.pi/180
    return elements

def get_sbdb_info(tdes):
    """
    Return small body orbit information from JPL SBDB API

    Parameters
    ----------
    tdes : str
        IMPORTANT: must be the designation of the small body, not the name.

    Returns
    -------
    elements : dict
        Dictionary containing desired cometary elements for the small body
    cov_mat : array
        Covariance matrix corresponding to cometary elements
    nongrav_params : dict
        Dictionary containing information about nongravitational acceleration
        constants for target body
    """
    # get json info from SBDB entry for small body
    raw_data = get_sbdb_raw_data(tdes)
    # cometary elements corresponding to the covariance on JPL SBDB
    elements = get_sbdb_elems(tdes, cov_elems=True)
    # covariance matrix for cometary orbital elements
    cov_exists = raw_data['orbit']['covariance'] is not None
    cov_keys = raw_data['orbit']['covariance']['labels'] if cov_exists else []
    if cov_exists:
        cov_mat = (np.array(raw_data['orbit']['covariance']['data'])).astype(float)
        # convert covariance matrix angle blocks from degrees to radians
        cov_mat[3:6, :] *= np.pi/180
        cov_mat[:, 3:6] *= np.pi/180
    else:
        print(f"WARNING: No covariance data available for {tdes}")
        cov_mat = np.zeros((6, 6))
    # nongravitatinoal constants for target body
    nongrav_data = raw_data['orbit']['model_pars']
    hdr = [param['name'] for param in nongrav_data]
    val = [float(param['value']) for param in nongrav_data]
    nongrav_data_dict = dict(zip(hdr, val))
    nongrav_keys = ['A1', 'A2', 'A3', 'ALN', 'NK', 'NM', 'NN', 'R0']
    # from https://ssd.jpl.nasa.gov/horizons/manual.html
    nongrav_default_vals = [0.0, 0.0, 0.0, 0.1112620426, 4.6142, 2.15, 5.093, 2.808]
    nongrav_default_dict = dict(zip(nongrav_keys, nongrav_default_vals))
    nongrav_key_map = { 'A1': 'a1', 'A2': 'a2', 'A3': 'a3', 'ALN': 'alpha',
                        'NK': 'k', 'NM': 'm', 'NN': 'n', 'R0': 'r0_au'}
    nongrav_params = {}
    for key in nongrav_keys:
        try:
            nongrav_params[nongrav_key_map[key]] = nongrav_data_dict[key]
        except KeyError:
            nongrav_params[nongrav_key_map[key]] = nongrav_default_dict[key]
        if key in cov_keys:
            elements[nongrav_key_map[key]] = nongrav_data_dict[key]
    abs_mag = 0
    albedo = 0.125
    for par in raw_data['phys_par']:
        if par['name'] == 'H':
            abs_mag = float(par['value'])
        if par['name'] == 'albedo':
            albedo = float(par['value'])
    body_radius = 1329*10**(-0.2*abs_mag)/(2*albedo**0.5)*1e3
    if abs_mag == 0:
        body_radius = 0
    nongrav_params['radius'] = body_radius
    return [elements, cov_mat, nongrav_params]

def get_similarity_stats(sol_1, cov_1, sol_2, cov_2):
    """
    Get similarity statistics between two solutions. This includes the
    Mahalanobis distance of both solutions from the other, Bhattacharyya
    distance, and Bhattacharyya coefficient.

    Parameters
    ----------
    sol_1 : vector
        Mean solution 1
    cov_1 : array
        Covariance matrix 1
    sol_2 : vector
        Mean solution 2
    cov_2 : _type_
        Covariance matrix 2

    Returns
    -------
    maha_dist_1 : float
        Mahalanobis distance of solution 2 from solution 1
    maha_dist_2 : float
        Mahalanobis distance of solution 1 from solution 2
    bhattacharya : float
        Bhattacharyya distance between solutions
    bhatt_coeff : float
        Bhattacharyya coefficient of solutions
    """
    maha_dist_1 = np.sqrt((sol_1-sol_2) @ np.linalg.inv(cov_1) @ (sol_1-sol_2).T)
    maha_dist_2 = np.sqrt((sol_1-sol_2) @ np.linalg.inv(cov_2) @ (sol_1-sol_2).T)
    big_cov = (cov_1+cov_2)/2
    det_1 = np.linalg.det(cov_1)
    det_2 = np.linalg.det(cov_2)
    det_big = np.linalg.det(big_cov)
    term_1 = 1/8 * (sol_1-sol_2) @ big_cov @ (sol_1-sol_2).T
    # simplify the natural log in term_2 because sometimes the determinant
    # is too small and the product of the two determinants is beyond the lower
    # limit of the float64 type
    # term_2 = 1/2 * np.log(det_big/np.sqrt(det_1*det_2))
    # term_2 = 1/2 * (np.log(det_big) - np.log(np.sqrt(det_1*det_2)))
    # term_2 = 1/2 * (np.log(det_big) - 1/2 * np.log(det_1*det_2))
    term_2 = 1/2 * np.log(det_big) - 1/4 * (np.log(det_1) + np.log(det_2))
    bhattacharya = term_1 + term_2
    bhatt_coeff = np.exp(-bhattacharya)
    return maha_dist_1, maha_dist_2, bhattacharya, bhatt_coeff
