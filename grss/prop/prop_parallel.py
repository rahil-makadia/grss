"""Parallel computing utilities for the GRSS orbit propagation code"""
import os
import time
import numpy as np
import pandas as pd
# pylint: disable=no-name-in-module
from .. import libgrss

__all__ = [
    'parallel_propagate',
    'cluster_ca_or_impacts',
    'reconstruct_all_log_files',
]

def _handle_one_cloned_sim(sol, ref_nongrav):
    """
    Parse a solution dictionary into an IntegBody object.

    Parameters
    ----------
    sol : dict
        Dictionary of solution parameters.
    ref_nongrav : libgrss.NongravParameters
        Reference nongravitational parameters to use for propagating the orbits.

    Returns
    -------
    body : libgrss.IntegBody
        IntegBody object with the solution parameters.
    """
    mass = sol.get('mass', 0.0)
    radius = sol.get('radius', 0.0)
    a1 = sol['a1'] if ref_nongrav.a1Est else ref_nongrav.a1
    a2 = sol['a2'] if ref_nongrav.a2Est else ref_nongrav.a2
    a3 = sol['a3'] if ref_nongrav.a3Est else ref_nongrav.a3
    alpha = ref_nongrav.alpha
    k = ref_nongrav.k
    m = ref_nongrav.m
    n = ref_nongrav.n
    r0_au = ref_nongrav.r0_au
    ng_list = [a1, a2, a3, alpha, k, m, n, r0_au]
    cometary = False
    cart = False
    if all(key in sol for key in ['e', 'q', 'tp', 'om', 'w', 'i']):
        e = sol['e']
        q = sol['q']
        tp = sol['tp']
        om = sol['om']
        w = sol['w']
        i = sol['i']
        com = [e, q, tp, om, w, i]
        cometary = True
    elif all(key in sol for key in ['x', 'y', 'z', 'vx', 'vy', 'vz']):
        pos = [sol['x'], sol['y'], sol['z']]
        vel = [sol['vx'], sol['vy'], sol['vz']]
        cart = True
    else:
        # rause error and print dictionary
        raise ValueError(f'Invalid solution dictionary:\n{sol}')
    if cometary:
        full_list = [sol['t'], mass, radius] + com + ng_list
    elif cart:
        full_list = [sol['t'], mass, radius] + pos + vel + ng_list
    return full_list

def parallel_propagate(ref_sol, ref_nongrav, ref_sim, clones, reconstruct=False):
    """
    Propagate multiple simulations in parallel using a reference simulation.

    Parameters
    ----------
    ref_sol : dict
        Reference solution to use for propagating the reference simulation.
    ref_nongrav : libgrss.NongravParameters
        Reference nongravitational parameters to use for propagating the orbits.
    ref_sim : libgrss.PropSimulation
        Reference simulation to use for propagating the orbits.
    clones : dict
        Dictionary of orbit solutions to propagate in parallel.
    """
    if all(key in ref_sol for key in ['e', 'q', 'tp', 'om', 'w', 'i']):
        is_cometary = True
    elif all(key in ref_sol for key in ['x', 'y', 'z', 'vx', 'vy', 'vz']):
        is_cometary = False
    else:
        raise ValueError('Invalid reference solution dictionary')

    # this does not need parallelism, only takes ~2.5s to do 1e6 clones!
    all_info = [_handle_one_cloned_sim(sol, ref_nongrav) for sol in clones]

    # save directory is ref_sim.name with spaces replaced by underscores
    save_dir = './logdir_'+ref_sim.name.replace(' ', '_')
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    else:
        # clear the directory if it already exists
        for file in os.listdir(save_dir):
            os.remove(os.path.join(save_dir, file))
    start_time = time.time()
    libgrss.propSim_parallel_omp(ref_sim, is_cometary, all_info)
    end_time = time.time()
    duration = end_time - start_time
    mm = int(duration / 60)
    ss = duration % 60
    print(f'Parallel propagation took {mm:02d} minute(s) and {ss:.6f} seconds')

    ca_list, impact_list = None, None
    if reconstruct:
        ca_list, impact_list = reconstruct_all_log_files(save_dir)
    else:
        print(('No reconstruction requested. '
                f'Log files can be reconstructed later from the "{save_dir}" directory'))
    return ca_list, impact_list

def reconstruct_all_log_files(log_dir):
    """
    Reconstruct CloseApproachParameters and ImpactParameters objects
    from all log files in a directory.

    Parameters
    ----------
    log_dir : str
        Path to the directory containing the log files.

    Returns
    -------
    ca_list : list of libgrss.CloseApproachParameters
        List of close approaches in all log files.
    impact_list : list of libgrss.ImpactParameters
        List of impacts in all log files.
    """
    start_time = time.time()
    ca_list = []
    impact_list = []
    files = os.listdir(log_dir)
    files.sort()
    for file in files:
        if file.endswith('.log'):
            log_file = os.path.join(log_dir, file)
            ca_list_, impact_list_ = _reconstruct_one_log_file(log_file)
            ca_list.append(ca_list_)
            impact_list.append(impact_list_)
    end_time = time.time()
    duration = end_time - start_time
    mm = int(duration / 60)
    ss = duration % 60
    print(f'Reconstruction took {mm:02d} minute(s) and {ss:.6f} seconds')
    return ca_list, impact_list

def _reconstruct_one_dataframe(lines, typ):
    """
    Reconstruct a pandas dataframe from a log file for close approaches or impacts.

    Parameters
    ----------
    lines : list of str
        List of lines from the log file.
    typ : str
        Type of the log file. Must be either 'ca' or 'impact'.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame of close approaches or impacts.

    Raises
    ------
    ValueError
        If typ is not 'ca' or 'impact'.
    """
    if typ not in ['ca', 'impact']:
        raise ValueError('typ must be either "ca" or "impact"')
    elif typ == 'ca':
        start_key = '$$CA_START'
        end_key = '$$CA_END'
    elif typ == 'impact':
        start_key = '$$IMPACT_START'
        end_key = '$$IMPACT_END'
    # find the index of the line that starts with "$$IMPACT_START" and "$$IMPACT_END"
    start = [i for i, line in enumerate(lines) if line.startswith(start_key)]
    end = [i for i, line in enumerate(lines) if line.startswith(end_key)]
    # if start and end keys do not exist, return None
    if not start or not end:
        return None
    start = start[0]
    end = end[0]
    # delimiter is "|"
    delimiter = ' | '
    # first line is the header
    header = lines[start+1].strip().split(delimiter)
    # create pandas dataframe
    df = pd.DataFrame([line.strip().split(delimiter) for line in lines[start+2:end]],
                        columns=header)
    # convert the columns to numeric
    str_cols = [
        'flybyBody', 'centralBody',
    ]
    list_cols = [
        'xRel', 'xRelMap', 'bVec',
        'kizner_dx', 'kizner_dy',
        'opik_dx', 'opik_dy',
        'scaled_dx', 'scaled_dy',
        'mtp_dx', 'mtp_dy',
        'xRelBodyFixed',
    ]
    bool_cols = [
        'impact',
    ]
    int_cols = [
        'flybyBodyIdx',
        'centralBodyIdx',
        'centralBodySpiceId',
    ]
    for col in df.columns:
        if col in str_cols:
            df[col] = df[col].astype(str)
        elif col in list_cols:
            df[col] = df[col].apply(lambda x: [float(num) for num in
                                                x.replace('[','').replace(']','').split(',')
                                                if num])
        elif col in bool_cols:
            df[col] = df[col].apply(lambda x: x == 'true')
        elif col in int_cols:
            df[col] = df[col].astype(int)
        else:
            df[col] = df[col].astype(float)
    return df

def _reconstruct_ca_params(row):
    """
    Reconstruct a CloseApproachParameters object from a row in a dataframe.

    Parameters
    ----------
    row : pandas.Series
        Row in a dataframe.

    Returns
    -------
    ca : libgrss.CloseApproachParameters
        CloseApproachParameters object.
    """
    ca = libgrss.CloseApproachParameters()
    ca.t = row['t']
    ca.xRel = row['xRel']
    ca.tMap = row['tMap']
    ca.xRelMap = row['xRelMap']
    ca.dist = row['dist']
    ca.vel = row['vel']
    ca.vInf = row['vInf']
    ca.flybyBody = row['flybyBody']
    ca.flybyBodyIdx = row['flybyBodyIdx']
    ca.centralBody = row['centralBody']
    ca.centralBodyIdx = row['centralBodyIdx']
    ca.centralBodySpiceId = row['centralBodySpiceId']
    ca.impact = row['impact']
    ca.tPeri = row['tPeri']
    ca.tLin = row['tLin']
    ca.bVec = row['bVec']
    ca.bMag = row['bMag']
    ca.gravFocusFactor = row['gravFocusFactor']
    ca.kizner.x = row['kizner_x']
    ca.kizner.y = row['kizner_y']
    ca.kizner.z = row['kizner_z']
    ca.kizner.dx = row['kizner_dx']
    ca.kizner.dy = row['kizner_dy']
    ca.opik.x = row['opik_x']
    ca.opik.y = row['opik_y']
    ca.opik.z = row['opik_z']
    ca.opik.dx = row['opik_dx']
    ca.opik.dy = row['opik_dy']
    ca.scaled.x = row['scaled_x']
    ca.scaled.y = row['scaled_y']
    ca.scaled.z = row['scaled_z']
    ca.scaled.dx = row['scaled_dx']
    ca.scaled.dy = row['scaled_dy']
    ca.mtp.x = row['mtp_x']
    ca.mtp.y = row['mtp_y']
    ca.mtp.z = row['mtp_z']
    ca.mtp.dx = row['mtp_dx']
    ca.mtp.dy = row['mtp_dy']
    return ca

def _reconstruct_impact_params(row):
    """
    Reconstruct an ImpactParameters object from a row in a dataframe.

    Parameters
    ----------
    row : pandas.Series
        Row in a dataframe.

    Returns
    -------
    impact : libgrss.ImpactParameters
        ImpactParameters object.
    """
    impact = libgrss.ImpactParameters()
    impact.t = row['t']
    impact.xRel = row['xRel']
    impact.dist = row['dist']
    impact.vel = row['vel']
    impact.vInf = row['vInf']
    impact.flybyBody = row['flybyBody']
    impact.flybyBodyIdx = row['flybyBodyIdx']
    impact.centralBody = row['centralBody']
    impact.centralBodyIdx = row['centralBodyIdx']
    impact.centralBodySpiceId = row['centralBodySpiceId']
    impact.impact = row['impact']
    impact.tPeri = row['tPeri']
    impact.tLin = row['tLin']
    impact.bVec = row['bVec']
    impact.bMag = row['bMag']
    impact.gravFocusFactor = row['gravFocusFactor']
    impact.kizner.x = row['kizner_x']
    impact.kizner.y = row['kizner_y']
    impact.kizner.z = row['kizner_z']
    impact.kizner.dx = row['kizner_dx']
    impact.kizner.dy = row['kizner_dy']
    impact.opik.x = row['opik_x']
    impact.opik.y = row['opik_y']
    impact.opik.z = row['opik_z']
    impact.opik.dx = row['opik_dx']
    impact.opik.dy = row['opik_dy']
    impact.scaled.x = row['scaled_x']
    impact.scaled.y = row['scaled_y']
    impact.scaled.z = row['scaled_z']
    impact.scaled.dx = row['scaled_dx']
    impact.scaled.dy = row['scaled_dy']
    impact.mtp.x = row['mtp_x']
    impact.mtp.y = row['mtp_y']
    impact.mtp.z = row['mtp_z']
    impact.mtp.dx = row['mtp_dx']
    impact.mtp.dy = row['mtp_dy']
    impact.xRelBodyFixed = row['xRelBodyFixed']
    impact.lon = row['lon']
    impact.lat = row['lat']
    impact.alt = row['alt']
    return impact

def _reconstruct_one_log_file(log_file):
    """
    Reconstruct CloseApproachParameters and ImpactParameters objects from a log file.

    Parameters
    ----------
    log_file : str
        Path to the log file.

    Returns
    -------
    ca_list : list of libgrss.CloseApproachParameters or None
        List of close approaches in the log file.
    impact_list : list of libgrss.ImpactParameters or None
        List of impacts in the log file.
    """
    # first read the log file
    with open(log_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    ca_df = _reconstruct_one_dataframe(lines, 'ca')
    impact_df = _reconstruct_one_dataframe(lines, 'impact')
    ca_list = []
    impact_list = []
    if ca_df is not None:
        for i in range(len(ca_df)):
            ca = _reconstruct_ca_params(ca_df.iloc[i])
            ca_list.append(ca)
    if impact_df is not None:
        for i in range(len(impact_df)):
            impact = _reconstruct_impact_params(impact_df.iloc[i])
            impact_list.append(impact)
    return ca_list, impact_list

def cluster_ca_or_impacts(full_list, max_duration=45, central_body=399):
    """
    Cluster a list of close approaches by time and uniqueness.

    Parameters
    ----------
    full_list : list of libgrss.CloseApproachParameters objects
        List of close approaches to cluster.
    max_duration : float
        Maximum duration (in days) between close approaches in a cluster.
    central_body : int
        SPICE ID of the central body.

    Returns
    -------
    all_clusters : tuple of list of libgrss.CloseApproachParameters objects
        A tuple of close approach clusters (each cluster is a list of
        close approaches).
    """
    all_clusters = []
    full_list = [ca_or_impact for ca_or_impact in full_list
                    if ca_or_impact.centralBodySpiceId == central_body]
    if not full_list:
        return tuple(all_clusters)
    times = [ca_or_impact.t for ca_or_impact in full_list]
    bodies = [ca_or_impact.flybyBody for ca_or_impact in full_list]
    idx_list = list(range(len(full_list)))

    df = pd.DataFrame({'time': times, 'body': bodies, 'idx': idx_list})
    df = df.sort_values(by=['time'])

    # create a new cluster if one of two conditions is met:
    # 1. the time between the current close approach and the last close approach
    #       is greater than max_duration
    # 2. the current close approach body is already in the current cluster
    cluster = [full_list[df.iloc[0]['idx']]]
    cluster_bodies = [df.iloc[0]['body']]
    for i in range(1, len(df)):
        time_condition = df.iloc[i]['time'] - df.iloc[i-1]['time'] > max_duration
        body_condition = df.iloc[i]['body'] in cluster_bodies
        if time_condition or body_condition:
            all_clusters.append(cluster)
            cluster = [full_list[df.iloc[i]['idx']]]
            cluster_bodies = [df.iloc[i]['body']]
        else:
            cluster.append(full_list[df.iloc[i]['idx']])
            cluster_bodies.append(df.iloc[i]['body'])
    all_clusters.append(cluster)
    return tuple(all_clusters)
