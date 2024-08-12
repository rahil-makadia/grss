---
title: 'GRSS: An open-source small-body science tool for planetary defense'
tags:
  - Python
  - C++
  - Asteroids
  - Comets
  - Orbit Determination
  - Orbit Propagator
authors:
  - name: Rahil Makadia
    corresponding: true
    orcid: 0000-0001-9265-2230
    affiliation: 1
  - name: Steven R. Chesley
    orcid: 0000-0003-3240-6497
    affiliation: 2
  - name: Davide Farnocchia
    orcid: 0000-0003-0774-884X
    affiliation: 2
  - name: Siegfried Eggl
    orcid: 0000-0002-1398-6302
    affiliation: 1
affiliations:
  - name: Department of Aerospace Engineering, University of Illinois at Urbana-Champaign, Urbana, IL 61801, USA
    index: 1
  - name: Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA 91109, USA
    index: 2
date: XX April 2024
bibliography: joss_paper.bib
aas-doi: 10.3847/PSJ/xxxxx
aas-journal: Planetary Science Journal
---

# Statement of Need

Understanding the motion of small solar system bodies is of utmost importance when looking at the problem through the lens of planetary defense. The ability to compute the orbit of an asteroid or a comet from various observations and then predicting the body's motion in the future is critical in understanding its impact hazard. The NASA Center for Near-Earth Object Studies (CNEOS) has developed a suite of state-of-the-art tools over the course of decades for this specific purpose. However, these tools are not publicly available. With the expected increase in the number of Near-Earth Object (NEO) observations as well as discoveries when new observatories such as the Vera C. Rubin Observatory come online [@Schwamb2023], there is a need for a publicly available library that can reliably perform both orbit propagation and determination for asteroids and comets.

# Summary

In this paper, we present ``GRSS``, the Gauss-Radau Small-body Simulator, an open-source library for orbit determination and propagation of small bodies in the solar system. ``GRSS`` is an open-source, MIT licensed software library with a C++ foundation and a Python binding. The propagator is based on the ``IAS15`` algorithm, a 15<sup>th</sup> order integrator based on Gauss-Radau quadrature [@Rein2014]. Only the particles of interest are integrated within ``GRSS`` to reduce computational cost. The states for the planets and Big16 main-belt asteroids are computed using memory-mapped SPICE digital ephemeris kernels as done by @Holman2023 in the ``ASSIST`` orbit propagator library. In addition to the propagator, the C++ portion of the library also has the ability to detect impacts and calculate close encounter circumstances using various formulations of the B-plane [@Kizner1961; @Opik1976; @Chodas1999; @Milani1999].

The C++ functionality is then exposed to Python through a binding generated using ``pybind11`` [@pybind11]. The Python layer of ``GRSS`` uses the propagator as the foundation to compute the orbits of small bodies from a given set of optical and/or radar astrometry from the Minor Planet Center, the JPL Small Body Radar Astrometry database, and the Gaia Focused Product Release (FPR) solar system observations database. Additionally, the orbit determination modules also have the ability to fit especially demanding measurements such as stellar occultations. These capabilities of the ``GRSS`` library have already been used to study the the heliocentric changes in the orbit of the (65803) Didymos binary asteroid system as a result of the DART impact [@Makadia2022; @Makadia2024].

``GRSS`` will continue to be developed in the future, with anticipated contributions including the ability to perform impact monitoring and kinetic impact deflection studies. ``GRSS`` is available to the community through the Python Package Index (PyPI) and the source code is available on GitHub. This availability will allow the research community to have access to a reliable and efficient tool for studying the dynamics of small bodies in the solar system.

# Acknowledgements

R.M. acknowledges funding from a NASA Space Technology Graduate Research Opportunities (NSTGRO) award, NASA contract No. 80NSSC22K1173. The work of S.R.C. and D.F. was carried out at the Jet Propulsion Laboratory, California Institute of Technology, under a contract with the National Aeronautics and Space Administration (No. 80NM0018D0004). This work has made use of data from the European Space Agency (ESA) mission Gaia, processed by the Gaia Data Processing and Analysis Consortium. This research has also extensively used data and services provided by the International Astronomical Union Minor Planet Center. Data from the MPC's database is made freely available to the public. Funding for this data and the MPC's operations comes from a NASA PDCO grant (80NSSC22M0024), administered via a University of Maryland - SAO subaward (106075-Z6415201). The MPC's computing equipment is funded in part by the above award, and in part by funding from the Tamkin Foundation.

# References
