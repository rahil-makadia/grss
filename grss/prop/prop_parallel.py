"""Parallel computing utilities for the GRSS orbit propagation code"""
import pandas as pd

__all__ = [ 'cluster_approaches',
]

def cluster_approaches(ca_list, max_duration=45):
    """
    Cluster a list of close approaches by time and uniqueness.

    Parameters
    ----------
    ca_list : list of prop_simulation.CloseApproachParameters objects
        List of close approaches to cluster.
    max_duration : float
        Maximum duration (in days) between close approaches in a cluster.

    Returns
    -------
    all_clusters : tuple of list of prop_simulation.CloseApproachParameters objects
        A tuple of close approach clusters (each cluster is a list of
        close approaches).
    """
    all_clusters = []
    times = [ca.t for ca in ca_list]
    bodies = [ca.flybyBody for ca in ca_list]
    idx_list = list(range(len(ca_list)))

    df = pd.DataFrame({'time': times, 'body': bodies, 'idx': idx_list})
    df = df.sort_values(by=['time'])

    # create a new cluster if one of two conditions is met:
    # 1. the time between the current close approach and the last close approach
    #       is greater than max_duration
    # 2. the current close approach body is already in the current cluster
    cluster = [ca_list[df.iloc[0]['idx']]]
    cluster_bodies = [df.iloc[0]['body']]
    for i in range(1, len(df)):
        time_condition = df.iloc[i]['time'] - df.iloc[i-1]['time'] > max_duration
        body_condition = df.iloc[i]['body'] in cluster_bodies
        if time_condition or body_condition:
            all_clusters.append(cluster)
            cluster = [ca_list[df.iloc[i]['idx']]]
            cluster_bodies = [df.iloc[i]['body']]
        else:
            cluster.append(ca_list[df.iloc[i]['idx']])
            cluster_bodies.append(df.iloc[i]['body'])
    all_clusters.append(cluster)
    return tuple(all_clusters)
