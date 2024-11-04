---
title: 'Gauss-Radau Small-body Simulator (GRSS): An Open-Source Library for Planetary Defense'
tags:
  - Python
  - C++
  - Asteroids
  - Comets
  - Orbit Determination
  - Orbit Propagation
authors:
  - name: Rahil Makadia
    corresponding: true
    orcid: 0000-0001-9265-2230
    affiliation: 1
  - name: Davide Farnocchia
    orcid: 0000-0003-0774-884X
    affiliation: 2
  - name: Steven R. Chesley
    orcid: 0000-0003-3240-6497
    affiliation: 2
  - name: Siegfried Eggl
    orcid: 0000-0002-1398-6302
    affiliation: 1
affiliations:
  - name: Department of Aerospace Engineering, University of Illinois at Urbana-Champaign, Urbana, IL 61801, USA
    index: 1
  - name: Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA 91109, USA
    index: 2
date: XX August 2024
bibliography: paper.bib
aas-doi: 10.3847/PSJ/xxxxx
aas-journal: Planetary Science Journal
---

# Statement of Need

Modeling the motion of small solar system bodies is of utmost importance when looking at the problem in the context of planetary defense. Reliably computing the orbit of an asteroid or a comet from various observations, and then predicting its trajectory in the future is critical in assessing the associated Earth-impact hazard. The NASA Center for Near-Earth Object Studies at the Jet Propulsion Laboratory (JPL) has developed a suite of state-of-the-art tools over the course of decades for this specific purpose. However, these tools are not publicly available. With the expected increase in the number of Near-Earth Object observations as well as discoveries when new observatories such as the Vera C. Rubin Observatory come online [@Schwamb2023], there is a need for a publicly available library that can reliably perform both orbit determination and propagation for small bodies in the solar system. Such a publicly available library will enable community efforts in planetary defense research.

# Summary

In this paper, we present ``GRSS``, the Gauss-Radau Small-body Simulator, an open-source library for orbit determination and propagation of small bodies in the solar system. ``GRSS`` is an open-source software library with a C++11 foundation and a Python binding. The propagator is based on the ``IAS15`` algorithm, a 15<sup>th</sup> order integrator based on Gauss-Radau quadrature [@Rein2014]. Only the particles of interest are integrated within ``GRSS`` to reduce computational cost. The states for the planets and 16 largest main-belt asteroids are computed using memory-mapped JPL digital ephemeris kernels [@Park2021] as done in the ``ASSIST`` orbit propagator library [@Holman2023]. In addition to the propagator, the C++ portion of the library also has the ability to predict impacts and calculate close encounter circumstances using various formulations of the B-plane [@Kizner1961; @Opik1976; @Chodas1999; @Milani1999; @Farnocchia2019].

The C++ functionality is exposed to Python through a binding generated using ``pybind11``[^1]. The Python layer of ``GRSS`` uses the propagator as the foundation to compute the orbits of small bodies from a given set of optical and/or radar astrometry from the Minor Planet Center[^2], the JPL Small Body Radar Astrometry database[^3], and the Gaia Focused Product Release solar system observations database[^4]. Additionally, the orbit determination modules also have the ability to fit especially demanding measurements such as stellar occultations. These capabilities of the ``GRSS`` library have already been used to study the the heliocentric changes in the orbit of the (65803) Didymos binary asteroid system as a result of the DART impact [@Makadia2022; @Makadia2024] and for analyzing the impact locations of two asteroids, 2024 BX<sub>1</sub> and 2024 RW<sub>1</sub> [@Makadia2024b].

``GRSS`` will continue to be developed in the future, with anticipated contributions including the ability to perform mission studies for future asteroid deflections. ``GRSS`` is publicly available to the community through the Python Package Index and the source code is available on GitHub[^5]. Therefore, ``GRSS`` is a reliable and efficient tool that the community has access to for studying the dynamics of small bodies in the solar system.

[^1]: <https://github.com/pybind/pybind11>
[^2]: <https://minorplanetcenter.net/>
[^3]: <https://ssd.jpl.nasa.gov/sb/radar.html>
[^4]: <https://www.cosmos.esa.int/web/gaia/fpr#SSOs>
[^5]: <https://github.com/rahil-makadia/grss>

# Acknowledgements

R.M. acknowledges funding from a NASA Space Technology Graduate Research Opportunities award, NASA contract No. 80NSSC22K1173. The work of S.R.C. and D.F. was carried out at the Jet Propulsion Laboratory, California Institute of Technology, under a contract with the National Aeronautics and Space Administration (No. 80NM0018D0004). This work has made use of data from the European Space Agency mission Gaia, processed by the Gaia Data Processing and Analysis Consortium. This library has also extensively used data and services provided by the International Astronomical Union Minor Planet Center (MPC). Data from the MPC's database is made freely available to the public. Funding for this data and the MPC's operations comes from a NASA Planetary Defense Coordination Office grant (80NSSC22M0024), administered via a University of Maryland - Smithsonian Astrophysical Observatory subaward (106075-Z6415201). The MPC's computing equipment is funded in part by the above award, and in part by funding from the Tamkin Foundation.

# References
