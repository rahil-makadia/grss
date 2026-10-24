GRSS Orbit Determination Module (grss.fit)
==========================================
The orbit determination code within GRSS is completely on the Python side of things, but it heavily relies on the C++ propagator binding. Given a small body orbit that needs to be fitted to a set of observations, the module uses the batch least squares algorithm to solve for an updated nominal orbit.

Gauss initial orbit determination
---------------------------------

``gauss_iod`` supplies preliminary Cartesian orbits from three geometric optical
directions and the corresponding central-body-relative observer positions.
It uses the classical eighth-degree Gauss polynomial and refines each physical
candidate with two-body differential correction. The numerical calculation
requires only NumPy; it does not query ephemerides or download data. Importing
the parent ``grss`` package still performs its usual initialization.

For heliocentric asteroid IOD, supply times in MJD TDB, right ascensions and
declinations in **radians**, Sun-to-observer positions in AU, and the solar
gravitational parameter in AU cubed per day squared. Directions and observer
positions must share the same inertial axes. In particular, use J2000 equatorial
axes when preparing a GRSS Cartesian initial solution::

    import numpy as np
    from grss.fit import gauss_iod

    # Select exactly three chronological, geometric observations.
    # ra_deg, dec_deg: arrays of length 3, in degrees.
    # mjd_tdb: array of length 3; observer_helio_au: array of shape (3, 3).
    # mu_sun: solar gravitational parameter in AU**3/day**2.
    candidates = gauss_iod(
        mjd_tdb, np.deg2rad(ra_deg), np.deg2rad(dec_deg),
        observer_helio_au, mu_sun,
    )

Each ``GaussIODSolution`` contains ``epoch``, a six-component ``state``
(position followed by velocity), three positive ``ranges``, the preliminary
polynomial ``root``, the maximum angular ``residual`` in radians, and the number
of refinement ``iterations``. Solutions are sorted by preliminary radius.
Multiple roots may describe different valid orbits; use additional observations
to choose a candidate. Three angular directions alone do not resolve that
ambiguity. Invalid or singular inputs, or failure of every candidate to converge,
raise ``ValueError``. ``max_iter`` and ``tol`` control refinement.

The returned state is **relative to the central body**. GRSS's Cartesian
``FitSimulation`` input is barycentric, so explicitly add the central body's
barycentric state at the middle observation epoch. After selecting a candidate
using additional observations::

    # selected_candidate is one of candidates, chosen by the caller.
    # sun_barycentric_state_at_epoch: six values in AU and AU/day, J2000.
    x_init = selected_candidate.to_fit_dict(sun_barycentric_state_at_epoch)
    # x_init has t, x, y, z, vx, vy, vz for FitSimulation.
    # Supply a separately justified covariance when constructing the fitter.

The adapter translates the origin; it does not change axes or units. For an
Earth-centered calculation, supply Earth-relative observer positions and the
Earth gravitational parameter in consistent units. Convert km/s results to
AU/day before constructing a GRSS initial solution.

This is a preliminary two-body calculation. It does not apply light-time,
aberration, gravitational deflection, astrometric debiasing, or perturbing-body
corrections, and it does not estimate a covariance. Observational corrections
must be handled by the caller. Very short arcs and nearly coplanar lines of
sight can be ill-conditioned; use well-separated observations and refine the
chosen preliminary orbit with the existing full dynamical fitter.

Offline numerical regression tests can be run from the repository root with::

    python -m unittest discover -s tests/python/iod -v

Batch least-squares fitting
--------------------------

The optical measurements are acquired using the `Minor Planet Center API <https://minorplanetcenter.net/mpcops/documentation/observations-api/>`_, and the radar measurements are acquired using the `JPL Small-body Radar API <https://ssd-api.jpl.nasa.gov/doc/sb_radar.html>`_. The optical astrometry is preprocessed to account for the following :

#. Star catalog biases [#]_
#. Measurement weighting [#]_
#. Gaia astrometry handling [#]_

Once the optical astrometry has been processed, the radar astrometry has been acquired, and the initial orbit is provided by the user, the least squares filter can be run. As of now, the filter can fit the nominal state, the nongravitational acceleration parameters, and any impulsive maneuver events. Currently the partial derivatives in the normal matrix are calculated using 1\ :sup:`st`-order central differences by default, but analytical partial derivatives are also a choice offered to the user. The filter also implements an outlier rejection scheme [#]_ to make sure any spurious measurements do not contaminate the fit.

.. rubric:: References
.. [#] Eggl, S., Farnocchia, D., Chamberlin, A.B., and Chesley, S.R., "Star catalog position and proper motion corrections in asteroid astrometry II: The Gaia era", Icarus, Volume 339, Pages 1-17, 2020. https://doi.org/10.1016/j.icarus.2019.113596.
.. [#] Vereš, P., Farnocchia, D., Chesley, S.R., and Chamberlin, A.B., "Statistical analysis of astrometric errors for the most productive asteroid surveys", Volume 296, Pages 139-149, 2017. https://doi.org/10.1016/j.icarus.2017.05.021.
.. [#] Fuentes-Muñoz, O., Farnocchia, D., Naidu, S.P., and Park, R.S., "Asteroid Orbit Determination Using Gaia FPR: Statistical Analysis", AJ, Volume 167, Pages 290-299, 2024. https://doi.org/10.3847/1538-3881/ad4291.
.. [#] Carpino, M., Milani, A., and Chesley, S.R., "Error statistics of asteroid optical astrometric observations", Icarus, Volume 166, Pages 248-270, 2003. https://doi.org/10.1016/S0019-1035(03)00051-4.
