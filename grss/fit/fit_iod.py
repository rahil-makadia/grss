"""Angles-only preliminary orbit determination, independent of ephemeris kernels."""
from dataclasses import dataclass

import numpy as np

__all__ = [
    'GaussIODSolution',
    'gauss_iod',
    '_stumpff',
    '_propagate',
    '_residuals',
    '_refine',
    'stumpff',
    'propagate_two_body',
    'gauss_angular_residuals',
    'refine_gauss_candidate',
    'universal_kepler_equation',
]


@dataclass
class GaussIODSolution:
    """One physical Gauss candidate at the middle observation epoch.

    ``state`` is central-body-relative position followed by velocity, in the
    input frame and units. ``ranges`` contains the three positive slant ranges.
    ``root`` is the preliminary middle radius before refinement.
    ``residual`` is the maximum angular residual in radians after refinement.
    """
    epoch: float
    state: np.ndarray
    ranges: np.ndarray
    root: float
    residual: float
    iterations: int

    def to_fit_dict(self, central_body_state):
        """Return a Cartesian ``FitSimulation`` initial-solution dictionary.

        The input observations must use MJD TDB, AU, days and J2000 equatorial
        axes. Supply the central body's six-component barycentric state in
        AU and AU/day at ``epoch``; this explicitly translates the origin.
        This method does not rotate frames or convert units.
        """
        origin = np.asarray(central_body_state, dtype=float)
        if origin.shape != (6,) or not np.all(np.isfinite(origin)):
            raise ValueError('central_body_state must be a finite six-component state.')
        state = self.state + origin
        return dict(t=self.epoch, **dict(zip(('x', 'y', 'z', 'vx', 'vy', 'vz'), state)))


def _stumpff(z):
    """Evaluate universal-variable Stumpff functions without small-z cancellation."""
    if abs(z) < 1e-4:
        c, s = 0.5, 1.0 / 6.0
        tc, ts = c, s
        for k in range(1, 8):
            tc *= -z / ((2*k + 1) * (2*k + 2))
            ts *= -z / ((2*k + 2) * (2*k + 3))
            c += tc
            s += ts
        return c, s
    if z > 0:
        q = np.sqrt(z)
        return (1 - np.cos(q)) / z, (q - np.sin(q)) / q**3
    q = np.sqrt(-z)
    if q > 700:
        raise ValueError('Universal anomaly overflow.')
    return (np.cosh(q) - 1) / (-z), (np.sinh(q) - q) / q**3


def universal_kepler_equation(anomaly, radius, radial_dot, alpha, dt, mu):
    """Return the universal Kepler equation's value and anomaly derivative.

    ``anomaly`` has units sqrt(length), ``radius`` is the initial distance,
    ``radial_dot`` is dot(initial position, initial velocity)/sqrt(mu), and
    ``alpha`` is 2/radius - speed**2/mu (inverse semi-major axis). Units must
    agree with the time interval ``dt`` and gravitational parameter ``mu``.
    """
    c, s = _stumpff(alpha*anomaly*anomaly)
    value = (radial_dot*anomaly**2*c + (1-alpha*radius)*anomaly**3*s
             + radius*anomaly - np.sqrt(mu)*dt)
    derivative = (radial_dot*anomaly*(1-alpha*anomaly**2*s)
                  + (1-alpha*radius)*anomaly**2*c + radius)
    return value, derivative


def _propagate(state, dt, mu):
    """Return a three-component position after two-body propagation.

    ``state`` contains initial central-body-relative position and velocity (six
    components). ``dt`` is the signed interval from that epoch and ``mu`` is the
    positive gravitational parameter. Use consistent length/time units. Supports
    elliptic, near-parabolic and hyperbolic motion. Does not modify the input.
    Raises ValueError if the universal Kepler solve fails.
    """
    r, v = state[:3], state[3:]
    radius = np.linalg.norm(r)
    if radius == 0 or not np.isfinite(radius):
        raise ValueError('Invalid central-body distance.')
    if dt == 0:
        return r.copy()
    sqrt_mu = np.sqrt(mu)
    rv = np.dot(r, v) / sqrt_mu
    alpha = 2 / radius - np.dot(v, v) / mu

    bound = max(sqrt_mu*abs(dt)/radius, np.sqrt(radius)*1e-6)
    lo, hi = (-bound, 0.0) if dt < 0 else (0.0, bound)
    for _ in range(100):
        if (universal_kepler_equation(lo, radius, rv, alpha, dt, mu)[0] <= 0
                <= universal_kepler_equation(hi, radius, rv, alpha, dt, mu)[0]):
            break
        lo, hi = (lo*2, hi) if dt < 0 else (lo, hi*2)
    else:
        raise ValueError('Could not bracket the universal anomaly.')
    x = (lo + hi) / 2
    for _ in range(100):
        value, deriv = universal_kepler_equation(x, radius, rv, alpha, dt, mu)
        if abs(value) <= 1e-13 * sqrt_mu*abs(dt):
            break
        if value < 0:
            lo = x
        else:
            hi = x
        trial = x - value/deriv if deriv > 0 else np.nan
        x = trial if lo < trial < hi else (lo+hi)/2
    else:
        raise ValueError('Universal anomaly did not converge.')
    c, s = _stumpff(alpha*x*x)
    return (1-x*x*c/radius)*r + (dt-x**3*s/sqrt_mu)*v


def _residuals(state, dt, observers, los, mu):
    """Return (cross-product residuals, slant ranges, maximum angular error).

    ``state`` is the six-component middle-epoch state; ``dt`` has three signed
    time offsets. ``observers`` and unit ``los`` have shape (3, 3) in matching
    inertial axes; observer positions are central-body-relative. ``mu`` uses
    the same length/time units. The flattened residual has nine dimensionless
    components, ranges have input length units, and angular error is in radians.
    """
    positions = np.array([_propagate(state, t, mu) for t in dt])
    relative = positions - observers
    ranges = np.linalg.norm(relative, axis=1)
    if np.any(ranges == 0) or not np.all(np.isfinite(relative)):
        raise ValueError('Invalid observer distance.')
    directions = relative / ranges[:, None]
    # Three cross-product components per direction avoid RA wrap/pole singularities.
    residuals = np.cross(los, directions).ravel()
    angles = np.arctan2(np.linalg.norm(np.cross(los, directions), axis=1),
                        np.sum(los*directions, axis=1))
    return residuals, ranges, np.max(angles)


def _refine(state, dt, observers, los, mu, max_iter, tol):
    """Refine a seed with two-body differential correction.

    Inputs follow ``gauss_angular_residuals``. ``max_iter`` is the nonnegative
    correction limit and ``tol`` the maximum accepted angular error in radians.
    Returns (state, three slant ranges, maximum angular error, iterations).
    Raises ValueError when convergence fails. Does not apply perturbations,
    observational corrections or weighting.
    """
    length = np.linalg.norm(state[:3])
    scale = np.array([length]*3 + [np.sqrt(mu/length)]*3)
    for iteration in range(max_iter + 1):
        residuals, ranges, angle = _residuals(state, dt, observers, los, mu)
        if angle <= tol:
            return state, ranges, angle, iteration
        if iteration == max_iter:
            break
        jac = np.empty((9, 6))
        step = 1e-5
        for j in range(6):
            delta = np.zeros(6)
            delta[j] = step*scale[j]
            plus = _residuals(state+delta, dt, observers, los, mu)[0]
            minus = _residuals(state-delta, dt, observers, los, mu)[0]
            jac[:, j] = (plus-minus)/(2*step)
        correction, _, rank, _ = np.linalg.lstsq(jac, -residuals, rcond=None)
        if rank < 6:
            break
        for damping in 0.5**np.arange(16):
            trial = state + damping*scale*correction
            try:
                trial_angle = _residuals(trial, dt, observers, los, mu)[2]
            except (ValueError, FloatingPointError):
                continue
            if trial_angle < angle:
                state = trial
                break
        else:
            break
    raise ValueError('Gauss candidate failed to converge to the angular tolerance.')


def gauss_iod(times, ra, dec, observer_positions, mu, *, max_iter=50, tol=1e-10):
    """Determine all surviving preliminary orbits from three optical directions.

    Parameters
    ----------
    times : array-like, shape (3,)
        Strictly increasing epochs in one uniform dynamical time scale.
        Only epoch differences enter the calculation (e.g. MJD TDB in days).
    ra, dec : array-like, shape (3,)
        Geometric right ascensions and declinations in radians, in the same
        inertial frame as the observer positions. RA may wrap; declinations
        must be in [-pi/2, pi/2].
    observer_positions : array-like, shape (3, 3)
        Central-body-to-observer positions, one row per epoch. For heliocentric
        IOD use Sun-relative observer positions, not barycentric positions.
    mu : float
        Positive central-body gravitational parameter, using the same length
        and time units as positions and times (e.g. AU**3/day**2).
    max_iter : int, optional
        Maximum differential-correction iterations for each positive root.
    tol : float, optional
        Maximum angular residual in radians for accepting a refined candidate.

    Returns
    -------
    list of GaussIODSolution
        Physical, converged candidates in increasing preliminary-root order.
        Different roots converging to the same state are deduplicated. Three
        directions can admit multiple orbits; additional observations must
        resolve ambiguity. No candidate is automatically chosen.

    Raises
    ------
    ValueError
        For invalid inputs, singular line-of-sight geometry, or no physical
        candidate converging within the requested tolerance and iteration limit.

    Notes
    -----
    Uses the classical eighth-degree Gauss polynomial followed by two-body
    differential correction. No light-time, aberration, gravitational
    deflection, measurement weighting, or perturbing-body corrections are
    applied. This is a preliminary orbit, not a covariance or a final fit.
    Extremely short arcs and nearly coplanar directions are ill-conditioned.

    References
    ----------
    Curtis, H. D., Orbital Mechanics for Engineering Students, 4th ed.,
    section 5.10 (Gauss's method of preliminary orbit determination).
    """
    times, ra, dec, observers = [np.asarray(a, dtype=float) for a in
                                (times, ra, dec, observer_positions)]
    if times.shape != (3,) or ra.shape != (3,) or dec.shape != (3,):
        raise ValueError('times, ra and dec must each have shape (3,).')
    if observers.shape != (3, 3):
        raise ValueError('observer_positions must have shape (3, 3).')
    if not all(np.all(np.isfinite(a)) for a in (times, ra, dec, observers)):
        raise ValueError('All inputs must be finite.')
    if not np.all(np.diff(times) > 0):
        raise ValueError('Observation times must be strictly increasing.')
    if np.any(np.abs(dec) > np.pi/2):
        raise ValueError('Declinations must be between -pi/2 and pi/2 radians.')
    if np.ndim(mu) != 0 or not np.isfinite(mu) or mu <= 0:
        raise ValueError('mu must be a finite positive scalar.')
    if (isinstance(max_iter, (bool, np.bool_))
            or not isinstance(max_iter, (int, np.integer)) or max_iter < 0):
        raise ValueError('max_iter must be a nonnegative integer.')
    if np.ndim(tol) != 0 or not np.isfinite(tol) or not 0 < tol < np.pi/2:
        raise ValueError('tol must be positive and less than pi/2 radians.')
    los = np.column_stack((np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)))
    if np.linalg.cond(los) > 1e12:
        raise ValueError('Line-of-sight geometry is singular or nearly coplanar.')
    # Normalize length and time before polynomial root finding. This makes the
    # same calculation work in AU/day and km/s without ill-scaled coefficients.
    length = np.max(np.linalg.norm(observers, axis=1))
    if length == 0:
        raise ValueError('Observer positions cannot all be at the central body.')
    time_scale = np.sqrt(length**3/mu)
    obs = observers/length
    dt = (times-times[1])/time_scale
    t1, t3 = dt[0], dt[2]
    span = t3-t1
    p = np.array([np.cross(los[1], los[2]), np.cross(los[0], los[2]),
                  np.cross(los[0], los[1])])
    d0 = np.dot(los[0], p[0])
    d = obs @ p.T
    a = (-d[0, 1]*t3/span + d[1, 1] + d[2, 1]*t1/span)/d0
    b = (d[0, 1]*(t3*t3-span*span)*t3/span
         + d[2, 1]*(span*span-t1*t1)*t1/span)/(6*d0)
    e = np.dot(obs[1], los[1])
    coefficients = [1, 0, -(a*a+2*a*e+np.dot(obs[1], obs[1])), 0, 0,
                    -2*b*(a+e), 0, 0, -b*b]
    roots = sorted(z.real for z in np.roots(coefficients)
                   if z.real > 0 and abs(z.imag) <= 1e-9*abs(z.real))
    solutions = []
    for root in roots:
        try:
            with np.errstate(over='raise', invalid='raise', divide='raise'):
                f1, f3 = 1-t1*t1/(2*root**3), 1-t3*t3/(2*root**3)
                g1, g3 = t1-t1**3/(6*root**3), t3-t3**3/(6*root**3)
                denom = f1*g3-f3*g1
                c1, c3 = g3/denom, -g1/denom
                matrix = np.column_stack((c1*los[0], -los[1], c3*los[2]))
                ranges = np.linalg.solve(matrix, obs[1]-c1*obs[0]-c3*obs[2])
                if np.any(ranges <= 0):
                    continue
                positions = obs + ranges[:, None]*los
                velocity = (-f3*positions[0]+f1*positions[2])/denom
                seed = np.concatenate((positions[1], velocity))
                state, ranges, residual, iterations = _refine(
                    seed, dt, obs, los, 1.0, max_iter, tol)
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            continue
        state *= np.array([length]*3 + [length/time_scale]*3)
        if any(np.linalg.norm((state-sol.state) / np.array(
                [length]*3 + [length/time_scale]*3)) < 1e-7 for sol in solutions):
            continue
        solutions.append(GaussIODSolution(float(times[1]), state, ranges*length,
                                          root*length, float(residual), iterations))
    if not solutions:
        raise ValueError('No physical Gauss candidate converged; check geometry, '
                         'units, arc length, tolerance and iteration limit.')
    return solutions


# Public helper names; the historical private spellings remain callable.
stumpff = _stumpff
propagate_two_body = _propagate
gauss_angular_residuals = _residuals
refine_gauss_candidate = _refine
