"""Utilities for the GRSS orbit determination code"""
import os
from json import loads
from requests import request
import numpy as np
from numba import jit

__all__ = [ 'get_ra_from_hms',
            'get_dec_from_dms',
            'mjd2et',
            'et2mjd',
            'get_radec',
            'parallax_constants_to_lat_lon_alt',
            'get_codes_dict',
            'get_observer_info',
            'get_sbdb_info',
            'get_similarity_stats',
]

def get_ra_from_hms(ra_hms):
    """
    Convert right ascension from HH:MM:SS.SSS to arcseconds.

    Parameters
    ----------
    ra_hms : str
        Right ascension in HH:MM:SS.SSS format

    Returns
    -------
    r_asc : float
        Right ascension in arcseconds
    """
    # convert right ascension from HH:MM:SS.SSS to arcseconds
    ra_split = ra_hms.split(' ')
    hour = float(ra_split[0])
    mins = 0.0 if ra_split[1] == '' else float(ra_split[1])
    secs = 0.0 if ra_split[2] == '' else float(ra_split[2])
    r_asc = 15*(hour+mins/60+secs/3600)
    return r_asc*3600

def get_dec_from_dms(dec_dms):
    """
    Convert declination from deg MM:SS.SSS to arcseconds.

    Parameters
    ----------
    dec_dms : str
        Declination in deg MM:arcSS.SSS format

    Returns
    -------
    dec : float
        Declination in arcseconds
    """
    # convert declination from deg MM:SS.SSS to arcseconds
    dec_split = dec_dms.replace('+','').replace('-','').split(' ')
    deg = float(dec_split[0])
    mins = 0.0 if dec_split[1] == '' else float(dec_split[1])
    secs = 0.0 if dec_split[2] == '' else float(dec_split[2])
    dec = deg+mins/60+secs/3600
    if dec_dms[0] == '-':
        dec *= -1
    return dec*3600

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

@jit('Tuple((float64, float64))(float64[:])', nopython=True, cache=True)
def get_radec(state):
    """
    Convert a cartesian state into right ascension and declination

    Parameters
    ----------
    state : vector
        6-element cartesian state vector

    Returns
    -------
    r_asc : float
        right ascension in arcseconds
    dec : float
        declination in arcseconds
    """
    pos = state[:3]
    # vel = state[3:6]
    dist = np.linalg.norm(pos) # calculate distance
    r_asc = np.arctan2(pos[1], pos[0])
    if r_asc < 0:
        r_asc = r_asc + 2*np.pi
    # end if
    dec = np.arcsin(pos[2]/dist) # calculate declination
    return r_asc*180/np.pi*3600, dec*180/np.pi*3600

@jit('Tuple((float64, float64, float64, float64, float64))(float64, float64, float64)',
        nopython=True, cache=False)
def parallax_constants_to_lat_lon_alt(lon, rho_cos_lat, rho_sin_lat):
    """
    Convert parallax constants to geodetic latitude, longitude, and altitude.

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
    geodet_lat : float
        radians, geodetic latitude of observer
    alt : float
        altitude of observer, m
    rho : float
        distance from the geocenter to the observer, m
    """
    # WGS84 ellipsoid parameters
    r_eq = 6378137.0 # m
    flatness = 1/298.257223563
    r_po = r_eq*(1 - flatness)
    e_sq = flatness*(2 - flatness)
    lat = np.arctan2(rho_sin_lat, rho_cos_lat)
    rho = np.sqrt(rho_cos_lat**2 + rho_sin_lat**2)*r_eq
    # from Hedgley, D. R., Jr., NASA TR R-458, 1976. https://ntrs.nasa.gov/citations/19760012748
    geodet_lat = np.arctan( np.tan(lat) / (1 - e_sq) )
    num = (r_eq**2*np.cos(geodet_lat))**2 + (r_po**2*np.sin(geodet_lat))**2
    den = (r_eq*np.cos(geodet_lat))**2 + (r_po*np.sin(geodet_lat))**2
    rellip = np.sqrt( num / den )
    alt = rho - rellip
    return lon*np.pi/180, lat, geodet_lat, alt, rho

def get_mpc_observatory_info():
    """
    Get the MPC observatory codes and return them as a dictionary.

    Returns
    -------
    codes_dict : dict
        dictionary of observatory codes and their corresponding longitude,
        rho*cos(latitude), and rho*sin(latitude)
    """
    file_dir = os.path.dirname(os.path.abspath(__file__))
    with open(f"{file_dir}/obs_codes.txt", "r", encoding='utf-8') as obs_code_file:
        data = obs_code_file.read()
    data = data.split('\n')[2:]
    mpc_info_dict = {}
    for line in data:
        lon_deg = line[4:13]
        # if lon_deg only has spaces, continue
        if lon_deg.strip() == '':
            continue
        code = str(line[:3])
        mpc_info_dict[code] = float(lon_deg), float(line[13:21]), float(line[21:30])
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
    # from obscodr.dat
    radar_codes_dict = {'-1' : (5.118130907583567, 0.3181656968486174, 6376485.333881727),
                        '-2' : (5.035480840855376, 0.7405706484785946, 6368498.652158532),
                        '-5' : (4.404902108756085, 0.5916728612960379, 6373576.72604664),
                        '-7' : (5.100642708478584, -0.399524038599547, 6379967.245894235),
                        '-9' : (4.889717923045817, 0.6675170779877466, 6370739.0833090255),
                        '-12': (4.244544747741349, 0.6129362736340168, 6371991.805482853),
                        '-13': (4.244736733959069, 0.6120180410509634, 6372122.4980728),
                        '-14': (4.243078671169674, 0.6151300479071157, 6371998.818089744),
                        '-15': (4.24311881374247, 0.6150585185761809, 6371969.864633745),
                        '-24': (4.2433352345697175, 0.6136333053130102, 6371971.305963156),
                        '-25': (4.243324762594206, 0.613591562161615, 6371985.244902656),
                        '-26': (4.243366650496253, 0.6135596331039699, 6371988.34181429),
                        '-27': (4.245049147895176, 0.6118624490546375, 6372102.910857716),
                        '-28': (4.24500900532238, 0.6118630223297832, 6372123.749433432),
                        '-30': (2.454981871184476, 0.6243606496656824, 6370885.734886958),
                        '-31': (2.41488642227841, 0.6274360123894954, 6372276.380152331),
                        '-34': (2.600226426206192, -0.6146518945371586, 6371696.895259182),
                        '-35': (2.6002176995599324, -0.6146063109997237, 6371693.310483042),
                        '-36': (2.600165339682372, -0.6145926584351252, 6371686.845710653),
                        '-37': (0.6623995750216518, 0.9719174569868173, 6363624.565513317),
                        '-38': (0.5792206735301062, 0.7853416799255691, 6367485.25229262),
                        '-39': (0.6624327362774398, 1.0069184644478415, 6363037.644793144),
                        '-43': (2.600214208901428, -0.6147209467847035, 6371696.186663469),
                        '-45': (2.6001513770483564, -0.6146513303468355, 6371676.046576824),
                        '-47': (2.6101416416867713, -0.5261377306704356, 6372958.625370879),
                        '-49': (2.5876862355306125, -0.5728673795561365, 6372243.0334759895),
                        '-51': (4.920559636257809, 0.798717284494754, 6367362.434760925),
                        '-53': (6.208949472775259, 0.702271825398319, 6370038.345634378),
                        '-54': (6.208937255470495, 0.7022465885003641, 6370019.897420738),
                        '-55': (6.2089634354092755, 0.702223144300857, 6370011.116026489),
                        '-56': (6.208972162055535, 0.7022516910577987, 6370027.803387983),
                        '-58': (0.2032774507822284, 0.7736762901596285, 6367689.985118462),
                        '-59': (0.19335106519443584, 0.8323411726917304, 6367084.059003408),
                        '-62': (0.12445593830121165, 0.8801430084714452, 6365696.804677316),
                        '-63': (6.209043720554868, 0.7023436744200388, 6370053.302865111),
                        '-65': (6.208996596665064, 0.7022748594030155, 6370015.505758979),
                        '-71': (0.12014176364296197, 0.878526528349089, 6365864.205220217),
                        '-72': (4.163102448855539, 0.709077390570306, 6370066.157584699),
                        '-73': (0.3355547330604022, 1.2123119930334905, 6359465.28889241),
                        '-74': (2.3354215814351083, -0.5531874093138018, 6372379.120672884),
                        '-75': (2.5732802878846512, -0.7437500522399563, 6368332.037514396),
                        '-76': (2.5732802878846512, -0.7437500522399563, 6368332.037514396),
                        '-77': (2.3065276556683423, -0.24930848030288796, 6377014.525085852),
                        '-78': (2.0131604976883715, -0.5041219522325624, 6373379.289639647),
                        '-79': (2.6022806787357897, -0.542750689803241, 6373270.145477656),
                        '-80': (3.5699749918455415, 0.34346577655008154, 6379419.660235015),
                        '-81': (4.218860482469002, 0.6465820287430366, 6371559.134510055),
                        '-82': (4.335181441126667, 0.5547313111571804, 6374103.163770516),
                        '-83': (4.3961492825573325, 0.5955447712959684, 6373752.421762615),
                        '-84': (4.4288497714226995, 0.6212127637775329, 6372836.487310473),
                        '-85': (4.469006306852585, 0.5317428681811908, 6374241.915461778),
                        '-86': (4.194317662527457, 0.8367083415513576, 6366576.525540495),
                        '-87': (4.684915752629047, 0.7257130595739807, 6368924.758078582),
                        '-88': (5.0267821198634355, 0.7459841618649351, 6368570.241645971),
                        '-89': (5.155988844388577, 0.3046527239435524, 6376222.04666057),
                        '-91': (0.12553455177894413, 0.8942928579278642, 6365319.001972028),
                        '-92': (0.4832012101987129, -0.4492288614371409, 6375505.90823571),
                        '-93': (0.16135568934687575, 0.685984303210341, 6370122.728403139),
                        '-96': (2.923065213935339, 0.1629031860614893, 6377644.65762988),
                        '-97': (5.484002533349886, -1.5706763615148356, 6359749.6043233145)
                        }
    return radar_codes_dict

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
        lon, lat, _, _, rho = parallax_constants_to_lat_lon_alt(info[0], info[1], info[2])
        codes_dict[code] = (lon, lat, rho)
    radar_codes_dict = get_radar_codes_dict()
    for code, info in radar_codes_dict.items():
        codes_dict[code] = info
    return codes_dict

def get_observer_info(observer_codes):
    """
    Creates a list of observer information for each observer code for
    passing into a grss propSimulation object

    Parameters
    ----------
    observer_codes : list
        Observer locations for each observation. Shape of each entry changes
        on type of observation: 1 (optical), 2 (radar), or 3 (doppler).

    Returns
    -------
    observer_info : list
        Information about each observer (Spice ID for parent body, longitude,
        geodetic latitude, and altitude). Shape of each entry changes on type
        of observation: 4 (optical), 8 (radar), or 9 (doppler).
    """
    codes_dict = get_codes_dict()
    observer_info = []
    body_id = 399 # default to Earth
    for code in observer_codes:
        info_list = []
        # check if code is a tuple - corresponds to a radar observation
        if isinstance(code, tuple) and len(code) in {2, 3}:
            tx_code = code[0][0]
            rx_code = code[0][1]
            bounce_point = code[1]
            receiver_info = codes_dict[rx_code]
            transmitter_info = codes_dict[tx_code]
            info_list.extend((body_id,receiver_info[0],receiver_info[1],receiver_info[2],
                                body_id,transmitter_info[0],transmitter_info[1],transmitter_info[2],
                                bounce_point))
            if len(code) == 3: # add tranmission frequency if dopppler observation
                freq = code[2]
                info_list.append(freq)
        # for geocentric occultations code is a tuple but needs to be decomposed
        elif isinstance(code, tuple) and code[0] in {'275', 'S/C'}:
            info_list.extend((body_id, code[1], code[2], code[3]))
            # info_list.extend((500, code[1], code[2], code[3], code[4], code[5], code[6]))
        else:
            info = codes_dict[code]
            info_list.extend((body_id, info[0], info[1], info[2]))
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
    req = request("GET", url, timeout=30)
    return loads(req.text)

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
    raw_data = get_sbdb_raw_data(tdes)
    if cov_elems:
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
        # if key in {'om', 'w', 'i'}:
        #     elements[key] = full_elements_dict[key]*np.pi/180
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
    cov_keys = raw_data['orbit']['covariance']['labels']
    cov_mat = (np.array(raw_data['orbit']['covariance']['data'])).astype(float)
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
    return [elements, cov_mat, nongrav_params]

def icrf2radec(pos_x, pos_y, pos_z, deg=False):
    """
    Convert ICRF unit vector to Right Ascension and Declination.

    Parameters
    ----------
    pos_x : float or array
        ICRF unit vector x-component.
    pos_y : float or array
        ICRF unit vector y-component.
    pos_z : float or array
        ICRF unit vector z-component.
    deg : bool, optional
        Flag for whether angles are in degrees, by default False

    Returns
    -------
    r_asc : float or array
        Right Ascension in degrees or radians.
    dec : float or array
        Declination in degrees or radians.
    """
    pos = np.array([pos_x, pos_y, pos_z])
    dist = np.linalg.norm(pos,axis=0) if (pos.ndim>1) else np.linalg.norm(pos)
    phi = np.arctan2(pos_y/dist,pos_x/dist)
    delta = np.arcsin(pos_z/dist)
    if deg:
        r_asc = np.mod(np.rad2deg(phi)+360,360)
        dec = np.rad2deg(delta)
    else:
        r_asc = np.mod(phi+2*np.pi,2*np.pi)
        dec = delta
    return r_asc, dec

def radec2icrf(r_asc, dec, deg=False):
    """
    Convert Right Ascension and Declination to ICRF unit vector.

    Parameters
    ----------
    r_asc : float or array
        Right Ascension in degrees or radians.
    dec : float or array
        Declination in degrees or radians.
    deg : bool, optional
        Flag for whether angles are in degrees, by default False

    Returns
    -------
    unit_vector : array
        ICRF unit vector.
    """
    if deg:
        alpha = np.deg2rad(r_asc)
        delta = np.deg2rad(dec)
    else:
        alpha = np.array(r_asc)
        delta = np.array(dec)
    cosd = np.cos(delta)
    pos_x = cosd*np.cos(alpha)
    pos_y = cosd*np.sin(alpha)
    pos_z = np.sin(delta)
    return np.array([pos_x, pos_y, pos_z])

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
    term_2 = 1/2 * np.log(det_big/np.sqrt(det_1*det_2))
    bhattacharya = term_1 + term_2
    bhatt_coeff = np.exp(-bhattacharya)
    return maha_dist_1, maha_dist_2, bhattacharya, bhatt_coeff
