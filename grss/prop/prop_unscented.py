"""Unscented transformations for the GRSS orbit propagation code"""
import numpy as np

def generate_sigma_points(sol, cov):
    """
    Generate sigma points for propagating uncertainties using the unscented
    transformation.

    Parameters
    ----------
    sol : dict
        solution dictionary from the SBDB API
    cov : numpy.ndarray
        covariance matrix of the solution

    Returns
    -------
    sigma_points : list
        list of dictionaries containing each of the sigma points
    """
    assert len(sol)-1 == cov.shape[0]
    assert cov.shape[0] == cov.shape[1]
    dim = cov.shape[0]
    sqrt_cov = np.linalg.cholesky(cov)
    sol_array = np.array([sol[key] for key in sol if key != 't'])
    sigma_points = [list(sol_array)]
    fac = dim**0.5
    for i in range(dim):
        plus = sol_array + fac*sqrt_cov.T[i]
        minus = sol_array - fac*sqrt_cov.T[i]
        sigma_points.append(list(plus))
        sigma_points.append(list(minus))
    # convert sigma points to dict
    for i in range(2*dim+1):
        sigma_points[i] = {key: sigma_points[i][j-1] for j, key in enumerate(sol) if key != 't'}
    return sigma_points
