"""Utilities for the GRSS orbit propagation code"""
import numpy as np
import matplotlib.pyplot as plt

__all__ = [ 'plot_solar_system'
]

def _add_planet_xy(axis, dist, color, alpha, lstyle):
    """
    Add a planet to the x-y plane plot.

    Parameters
    ----------
    axis : matplotlib axis
        axis to plot on
    dist : float
        semi-major axis of the planet's orbit
    color : str
        color of the planet's orbit
    alpha : float
        transparency of the planet's orbit
    lstyle : str
        linestyle of the planet's orbit

    Returns
    -------
    None : NoneType
        None
    """
    circle = plt.Circle((0, 0), dist, color=color, fill=False, alpha=alpha, linestyle=lstyle)
    axis.add_artist(circle)
    return None

def _add_planet_xz(axis, dist, obliq, color, alpha, lstyle):
    """
    Add a planet to the x-z plane plot.

    Parameters
    ----------
    axis : matplotlib axis
        axis to plot on
    dist : float
        semi-major axis of the planet's orbit
    obliq : float
        obliquity of the planet's orbit
    color : str
        color of the planet's orbit
    alpha : float
        transparency of the planet's orbit
    lstyle : str
        linestyle of the planet's orbit

    Returns
    -------
    None : NoneType
        None
    """
    sin = dist*np.sin(obliq*np.pi/180)
    cos = dist*np.cos(obliq*np.pi/180)
    axis.plot([-cos, cos], [-sin, sin], color=color, alpha=alpha, linestyle=lstyle)
    return None

def plot_solar_system(axis, xy_plane, alpha=0.25, lstyle='--'):
    """
    Plots the solar system in the x-y or x-z plane. For visual purposes only.

    Parameters
    ----------
    axis : matplotlib axis
        axis to plot on
    xy_plane : bool
        True to plot in x-y plane, False to plot in x-z plane
    alpha : float, optional
        transparency of the planet's orbit, by default 0.25
    lstyle : str, optional
        linestyle of the planet's orbit, by default '--'

    Returns
    -------
    None : NoneType
        None
    """
    # draw the sun
    axis.scatter(0, 0, s=50, c='yellow', marker='o', edgecolors='gray', alpha=alpha)
    # draw the orbit of the planets as circles for x-y plane
    num_planets = 9
    dists = [0.387, 0.723, 1.0, 1.524, 5.203, 9.539, 19.18, 30.06, 39.53]
    colors = ['gray', 'orange', 'blue', 'red', 'brown', 'orange', 'cyan', 'blue', 'gray']
    if xy_plane:
        for i in range(num_planets):
            _add_planet_xy(axis, dists[i], colors[i], alpha, lstyle)
    # draw the orbit of the planets as inclined lines for x-z plane
    else:
        obliqs = [7.25, 3.39, 0.0, 1.85, 1.3, 2.49, 0.77, 1.77, 0.0]
        for i in range(num_planets):
            _add_planet_xz(axis, dists[i], obliqs[i]+23.44, colors[i], alpha, lstyle)
    return None
