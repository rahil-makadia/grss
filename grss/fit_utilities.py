import numpy as np
from numba import jit
from astroquery.mpc import MPC
from json import loads
from requests import request

__all__ = [ 'get_ra_from_hms',
            'get_dec_from_dms',
            'mjd2et',
            'et2mjd',
            'get_radec',
            'parallax_constants_to_lat_lon_alt',
            'get_codes_dict',
            'get_observer_info',
            'get_sbdb_info',
]

def get_ra_from_hms(ra_hms):
    """Convert right ascension from HH:MM:SS.SSS to radians.

    Args:
        ra_hms (str): Right ascension in HH:MM:SS.SSS format

    Returns:
        ra: Right ascension in radians
    """
    # convert right ascension from HH:MM:SS.SSS to radians
    ra_split = ra_hms.split(' ')
    ra = 15*(float(ra_split[0])+float(ra_split[1])/60+float(ra_split[2])/3600)
    return ra*3600 # np.deg2rad(ra)

def get_dec_from_dms(dec_dms):
    """Convert declination from deg MM:SS.SSS to radians.

    Args:
        dec_dms (str): Declination in deg MM:arcSS.SSS format

    Returns:
        dec: Declination in radians
    """
    # convert declination from deg MM:SS.SSS to radians
    dec_split = dec_dms.replace('+','').replace('-','').split(' ')
    dec = (float(dec_split[0])+float(dec_split[1])/60+float(dec_split[2])/3600)
    if dec_dms[0] == '-':
        dec *= -1
    return dec*3600 # np.deg2rad(dec)

@jit('float64(float64)', nopython=True, cache=True)
def mjd2et(MJD):
    """
    Converts modified Julian Date to JPL NAIF SPICE Ephemeris time.
    Only valid for TDB timescales.

    Parameters:
    -----------
    MJD ... Modified Julian Day

    Returns:
    --------
    ET  ... Ephemeris time (ephemeris seconds beyond epoch J2000)
    """
    return (MJD+2400000.5-2451545.0)*86400

@jit('float64(float64)', nopython=True, cache=True)
def et2mjd(et):
    """
    Converts JPL NAIF SPICE Ephemeris time to Modified Julian Date.
    Only valid for TDB timescales.

    Parameters:
    -----------
    ET  ... Ephemeris time (ephemeris seconds beyond epoch J2000)

    Returns:
    --------
    MJD ... Modified Julian Day
    """
    return et/86400-2400000.5+2451545.0

@jit('Tuple((float64, float64))(float64[:])', nopython=True, cache=True)
def get_radec(state):
    """Convert a cartesian state into right ascension and declination

    Args:
        state (vector): 6-element cartesian state vector

    Returns:
        ra (float): right ascension in radians
        dec (float): declination in radians
    """
    pos = state[:3]
    # vel = state[3:6]
    dist = np.linalg.norm(pos) # calculate distance
    ra = np.arctan2(pos[1], pos[0])
    if ra < 0:
        ra = ra + 2*np.pi
    # end if
    dec = np.arcsin(pos[2]/dist) # calculate declination
    return ra*180/np.pi*3600, dec*180/np.pi*3600 # ra, dec

@jit('Tuple((float64, float64, float64, float64, float64))(float64, float64, float64)', nopython=True, cache=False)
def parallax_constants_to_lat_lon_alt(lon, rho_cos_lat, rho_sin_lat):
    """Convert parallax constants to geodetic latitude, longitude, and altitude.

    Args:
        lon (float): degrees, longitude of observer
        rho_cos_lat (float): rho cos phi', where rho is the distance from the geocenter to the observer and phi' is the geocentric latitude of the observer
        rho_sin_lat (float): rho sin phi', where rho is the distance from the geocenter to the observer and phi' is the geocentric latitude of the observer

    Returns:
        lon (float): radians, longitude of observer
        lat (float): radians, geocentric latitude of observer
        geodet_lat (float): radians, geodetic latitude of observer
        alt (float): km, altitude of observer
    """
    # WGS84 ellipsoid parameters
    re = 6378137.0 # km
    f = 1/298.257223563
    rp = re*(1 - f)
    e_sq = f*(2 - f)
    
    lat = np.arctan2(rho_sin_lat, rho_cos_lat)
    rho = np.sqrt(rho_cos_lat**2 + rho_sin_lat**2)*re
    
    geodet_lat = np.arctan( np.tan(lat) / (1 - e_sq) ) # from Hedgley, D. R., Jr. "An Exact Transformation from Geocentric to Geodetic Coordinates for Nonzero Altitudes." NASA TR R-458, March, 1976. https://ntrs.nasa.gov/citations/19760012748
    rellip = np.sqrt(  ( (re**2*np.cos(geodet_lat))**2 + (rp**2*np.sin(geodet_lat))**2 ) / ( (re*np.cos(geodet_lat))**2 + (rp*np.sin(geodet_lat))**2 )  )
    alt = rho - rellip # km -> m
    return lon*np.pi/180, lat, geodet_lat, alt, rho

def get_codes_dict():
    """Creates a dictionary of MPC codes and their corresponding longitude, geodetic latitude, and altitude

    Args:
        codes (astropy Table): table of MPC codes and their corresponding longitude and parallax constants

    Returns:
        codes_dict (dict): dictionary of MPC codes and their corresponding longitude, geodetic latitude, and altitude
    """
    codes = MPC.get_observatory_codes()
    codes_dict = {}
    for i in range(len(codes)):
        if np.isnan(codes[i]['Longitude'])==False:
            code = codes[i]['Code']
            lon, rho_cos_lat, rho_sin_lat = float(codes[i]['Longitude']), float(codes[i]['cos']), float(codes[i]['sin'])
            # print(code, lon, rho_cos_lat, rho_sin_lat)
            lon, lat, geodet_lat, alt, rho = parallax_constants_to_lat_lon_alt(lon, rho_cos_lat, rho_sin_lat)
            codes_dict[code] = (lon, lat, rho)
    # from https://www.minorplanetcenter.net/iau/lists/ObsCodes.html
    extra_observatories = { 'Y98': (358.68692, 0.617953, +0.783613),
                            'U94': (246.30250, 0.792006, +0.608864),
                            'Y66': (343.49053, 0.881484, +0.471429),
                            'X06': (289.14649, 0.862374, -0.505110)}
    for code, info in extra_observatories.items():
        lon, lat, geodet_lat, alt, rho = parallax_constants_to_lat_lon_alt(info[0], info[1], info[2])
        codes_dict[code] = (lon, lat, rho)
    return codes_dict

def get_observer_info(observer_codes):
    """Creates a list of observer information for each observer code for passing into a grss propSimulation object

    Args:
        observer_codes (list): Observer locations for each observation. Shape of each entry changes on type of observation: 1 (optical), 2 (radar), or 3 (doppler).

    Returns:
        observer_info (list): Information about each observer (Spice ID for parent body, longitude, geodetic latitude, and altitude). Shape of each entry changes on type of observation: 4 (optical), 8 (radar), or 9 (doppler).
    """
    codes_dict = get_codes_dict()
    observer_info = []
    for code in observer_codes:
        info_list = [399]
        # check if code is a tuple - corresponds to a radar observation
        if isinstance(code, tuple) and len(code) in {2, 3}:
            receiver_info = codes_dict[code[0]]
            info_list.extend((receiver_info[0], receiver_info[1], receiver_info[2]))
            transmitter_info = codes_dict[code[1]]
            info_list.extend((399, transmitter_info[0], transmitter_info[1], transmitter_info[2]))
            if len(code) == 3: # add tranmission frequency if dopppler observation
                info_list.append(code[2])
        else:
            info = codes_dict[code]
            info_list.extend((info[0], info[1], info[2]))
        observer_info.append(info_list)
    return observer_info

def get_sbdb_raw_data(tdes):
    """Get json data from JPL SBDB entry for desired small body
    Args:
        tdes (str): IMPORTANT: must be the designation of the small body, not the name.
    Returns:
        raw_data (dict): JSON output of small body information query from SBDB
    """
    url = f"https://ssd-api.jpl.nasa.gov/sbdb.api?sstr={tdes}&cov=mat&phys-par=true&full-prec=true"
    r = request("GET",url)
    return loads(r.text)

def get_sbdb_elems(tdes, cov_elems=True):
    """Get a set of desired elements for small body from SBDB API
    Args:
        tdes (str): IMPORTANT: must be the designation of the small body, not the name.
        cov_elems (bool, optional): Boolean for whether to extract element set corresponding to the covariance information (since the nominal set on the webpage might have a different epoch than the full covariance info). Defaults to True.
    Returns:
        elements (dict): Dictionary containing desired cometary elements for the small body
    """
    raw_data = get_sbdb_raw_data(tdes)
    if cov_elems:
        # epoch of orbital elements at reference time [JD -> MJD]
        epoch_mjd = float(raw_data['orbit']['covariance']['epoch']) - 2400000.5
        # cometary elements at epoch_mjd
        if epoch_mjd == float(raw_data['orbit']['epoch']):
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
    for i in range(len(elem)):
        hdr.append(elem[i]['name'])
        val.append(float(elem[i]['value']))
    full_elements_dict = dict(zip(hdr, val))
    # cometary elements
    # eccentricity, perihelion distance, time of periapse passage (JD), longitude of the ascending node, argument of perihelion, inclination
    keys = ['e', 'q', 'tp', 'om', 'w', 'i']
    values = [full_elements_dict[key] for key in keys]
    elements = {'t': epoch_mjd}
    # add every element key to elements dictionary
    for key in keys:
        elements[key] = full_elements_dict[key]
        if key == 'tp':
            elements[key] = full_elements_dict[key] - 2400000.5
        # if key in {'om', 'w', 'i'}:
        #     elements[key] = full_elements_dict[key]*np.pi/180
    return elements

def get_sbdb_info(tdes):
    """Return small body information from JPL SBDB API
    Args:
        tdes (str): IMPORTANT: must be the designation of the small body, not the name.
    Returns:
        elements (dict): Dictionary containing desired cometary elements for the small body
        cov_mat (np.ndarray): Covariance matrix corresponding to cometary elements
        nongrav_params (dict): Dictionary containing information about nongravitational acceleration constants for target body
    """
    # get json info from SBDB entry for small body
    raw_data = get_sbdb_raw_data(tdes)
    # cometary elements corresponding to the covariance on JPL SBDB
    elements = get_sbdb_elems(tdes, cov_elems=True)
    # covariance matrix for cometary orbital elements
    cov_keys = raw_data['orbit']['covariance']['labels']
    cov_mat = (np.array(raw_data['orbit']['covariance']['data'])).astype(float)
    # nongravitatinoal constants for target body
    nongrav_data = raw_data['orbit']['model_pars']
    hdr = [param['name'] for param in nongrav_data]
    val = [float(param['value']) for param in nongrav_data]
    nongrav_data_dict = dict(zip(hdr, val))
    nongrav_keys = ['A1', 'A2', 'A3', 'ALN', 'NK', 'NM', 'NN', 'R0']
    nongrav_default_vals = [0.0, 0.0, 0.0, 0.1112620426, 4.6142, 2.15, 5.093, 2.808] # from https://ssd.jpl.nasa.gov/horizons/manual.html
    nongrav_default_dict = dict(zip(nongrav_keys, nongrav_default_vals))
    nongrav_key_map = {'A1': 'a1', 'A2': 'a2', 'A3': 'a3', 'ALN': 'alpha', 'NK': 'k', 'NM': 'm', 'NN': 'n', 'R0': 'r0_au'}
    nongrav_params = {}
    for key in nongrav_keys:
        try:
            nongrav_params[nongrav_key_map[key]] = nongrav_data_dict[key]
        except KeyError:
            nongrav_params[nongrav_key_map[key]] = nongrav_default_dict[key]
        if key in cov_keys:
            elements[nongrav_key_map[key]] = nongrav_data_dict[key]
    return [elements, cov_mat, nongrav_params]
