GRSS Propagator Module (grss.prop)
==================================
The propagator within GRSS is based on the IAS15 algorithm [#]_, which is a 15th-order integrator based on Gauss-Radau quadrature. The algorithm is adaptive, meaning that the step size is adjusted to ensure that the error is below a certain threshold. However, the integrator can also be run in fixed-step mode, where a user-defined fixed time step is used until the end of the integration.

The force model includes the following effects:

#. Newtonian point-mass gravity
#. Einstein-Infeld-Hoffmann (EIH) formulation of point-mass relativistic gravity [#]_ [#]_
#. J\ :sub:`2` zonal harmonic effects from the Sun and the Earth
#. A\ :sub:`1`, A\ :sub:`2`, A\ :sub:`3` nongravitational  acceleration model [#]_

In addition to any integrated bodies defined by the user, the default force model includes effects from the Sun, the planets, the Moon, Pluto, and the Big 16 main-belt asteroids. The states of these bodies are read in from either the JPL DE431 or DE441 ephemerides, based on the user's choice. The propagator can also interpolate any of the integrated bodies' states to any time within the integration (or within a margin around the integration time bounds). This interpolation is done directly via polynomial coefficients computed by IAS15 during propagation.

The propagator module also has the ability to calculate apparent states of any integrated body. This is done through two methods. The first is the simpler of the two, and uses a simple approximation of the light-time between an observer and the integrated body. The second method is more accurate, and uses an iterative method to calculate the light-time between an observer and the integrated body. The second method is slightly more computationally expensive since it requires a few extra interpolations.

This ability to calculate apparent states is used to calculate optical observables (Right Ascension and Declination). Furthermore, the propagator suite can also calculate range observables (delay and doppler). The delay measurement is trivially returned from the light-time calculation. The doppler measurable, however, is calculated using a dedicated doppler measurement model. This entire codebase is written in C++ and a Python binding is using `pybind11 <https://pybind11.readthedocs.io/en/stable/>`_.

Current and previous versions of the integrator are a hybrid amalgamation of the work done in the following sources:

#. `An efficient integrator that uses Gauss-Radau spacings <https://doi.org/10.1017/S0252921100083913>`_
#. `IAS15: a fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine precision over a billion orbits <https://doi.org/10.1093/mnras/stu2164>`_
#. `Moving Planets Around: An Introduction to N-Body Simulations Applied to Exoplanetary Systems <https://mitpress.mit.edu/9780262539340/moving-planets-around/>`_
#. `A new timestep criterion for N-body simulations <https://doi.org/10.21105/astro.2401.02849>`_

.. rubric:: References
.. [#] Rein, H. and Spiegel D.S., "IAS15: a fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine precision over a billion orbits", Monthly Notices of the Royal Astronomical Society, Volume 446, Issue 2, Pages 1424-1437, 2015. https://doi.org/10.1093/mnras/stu2164.
.. [#] Einstein, A., Infeld, L., and Hoffmann, B., “The Gravitational Equations and the Problem of Motion”, Annals of Mathematics, Volume 39, Issue 1, Pages 65-100, 1938. https://doi.org/10.2307/1968714.
.. [#] Moyer, T.D., "Mathematical formulation of the Double Precision Orbit Determination Program (DPODP)", Jet Propulsion Laboratory Technical Report 32-1527, 1971. https://ntrs.nasa.gov/citations/19710017134.
.. [#] Marsden, B.G., Sekanina, Z., and Yeomans, D.K., “Comets and nongravitational forces. V”, The Astronomical Journal, Volume 78, Issue 2, Pages 211-225, 1973. https://doi.org/10.1086/111402.
.. [#] Everhart, E., “An efficient integrator that uses Gauss-Radau spacings”, Dynamics of Comets: Their Origin and Evolution, Proceedings of IAU Colloquium 83, Pages 185-202, 1985. https://doi.org/10.1017/S0252921100083913.
.. [#] Roa, J., Hamers, A.S., Cai, M.X., and Leigh, N.W.C., "Moving Planets Around: An Introduction to N-Body Simulations Applied to Exoplanetary Systems", MIT Press, 2020. https://mitpress.mit.edu/9780262539340/moving-planets-around/.
.. [#] Pham, D., Rein, H., and Spiegel, D.S., “A new timestep criterion for N-body simulations”, The Open Journal of Astrophysics 7, 2024. https://doi.org/10.21105/astro.2401.02849.
