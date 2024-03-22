"""Parallel computing utilities for the GRSS orbit propagation code"""
import numpy as np
import pandas as pd
# pylint: disable=no-name-in-module
from .. import libgrss

__all__ = [ 'cluster_ca_or_impacts',
]

def _handle_one_cloned_sim(sol):
    """
    Parse a solution dictionary into an IntegBody object.

    Parameters
    ----------
    sol : dict
        Dictionary of solution parameters.

    Returns
    -------
    body : libgrss.IntegBody
        IntegBody object with the solution parameters.
    """
    mass = sol.get('mass', 0.0)
    radius = sol.get('radius', 0.0)
    a1 = sol.get('a1', 0.0)
    a2 = sol.get('a2', 0.0)
    a3 = sol.get('a3', 0.0)
    alpha = sol.get('alpha', 0.1112620426)
    k = sol.get('k', 4.6142)
    m = sol.get('m', 2.15)
    n = sol.get('n', 5.093)
    r0_au = sol.get('r0_au', 2.808)
    ng_list = [a1, a2, a3, alpha, k, m, n, r0_au]
    cometary = False
    cart = False
    if 'e' in sol and 'q' in sol and 'tp' in sol and 'om' in sol and 'w' in sol and 'i' in sol:
        e = sol['e']
        q = sol['q']
        tp = sol['tp']
        om = sol['om']*np.pi/180.0
        w = sol['w']*np.pi/180.0
        i = sol['i']*np.pi/180.0
        com = [e, q, tp, om, w, i]
        cometary = True
    elif 'x' in sol and 'y' in sol and 'z' in sol and 'vx' in sol and 'vy' in sol and 'vz' in sol:
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

def parallel_propagate(ref_sol, ref_sim, clones):
    """
    Propagate multiple simulations in parallel using a reference simulation.

    Parameters
    ----------
    ref_sol : dict
        Reference solution to use for propagating the reference simulation.
    ref_sim : libgrss.PropSimulation
        Reference simulation to use for propagating the orbits.
    clones : dict
        Dictionary of orbit solutions to propagate in parallel.
    """
    if ('e' in ref_sol and 'q' in ref_sol and 'tp' in ref_sol
        and 'om' in ref_sol and 'w' in ref_sol and 'i' in ref_sol):
        is_cometary = True
    elif ('x' in ref_sol and 'y' in ref_sol and 'z' in ref_sol
            and 'vx' in ref_sol and 'vy' in ref_sol and 'vz' in ref_sol):
        is_cometary = False
    else:
        raise ValueError('Invalid reference solution dictionary')

    # this does not need parallelism, only takes ~2.5s to do 1e6 clones!
    all_info = [_handle_one_cloned_sim(sol) for sol in clones]
    ref_sim.parallelMode = True

    return libgrss.propSim_parallel_omp(ref_sim, is_cometary, all_info)

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
