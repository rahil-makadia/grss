---
title: 'GRSS: An open-source small-body science tool for planetary defense'
tags:
  - Python
  - C++
  - astronomy
  - asteroids
  - comets
  - orbit-determination
  - orbit-propagator
authors:
  - name: Rahil Makadia
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
date: 18 September 2023
bibliography: joss_paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/PSJ/xxxxx
aas-journal: Planetary Science Journal
---

# Statement of Need

Understanding the motion of small solar system bodies is of utmost importance when looking at the problem through the lens of planetary defense. The ability to compute the orbit of an asteroid or a comet from various observation sources and then predicting the body's motion in the future is critical in determining whether the Earth will remain safe from impacts in the future. The NASA Center for Near-Earth Object Studies (CNEOS) has developed a suite of state-of-the-art tools over the course of decades for this exact purpose. However, these tools are not publicly available. With the expected increase in the number of Near-Earth Object (NEO) observations as well as discoveries when new observatories such as the Rubin Observatory come online, there is a need for a publicly available library that can reliably perform both orbit determination and propagation for asteroids and comets.

# Summary

In this paper, we present ``GRSS``, the Gauss-Radau Small-body Simulator, an open-source library for orbit determination and propagation of small bodies in the solar system. ``GRSS`` is an open-source, GPL-3.0 licensed software library with a C++ foundation and a Python binding. The propagator is based on the ``IAS15`` algorithm [@Rein2014], a 15<sup>th</sup> order integrator based on Gauss-Radau quadrature. Only the particles of interest are integrated within ``GRSS`` to reduce computational cost. The states for the planets and Big16 main-belt asteroids are computed using memory-mapped SPICE digital ephemeris kernels as done by @Holman2023 in the ``ASSIST`` propagator library. In addition to the propagator, the C++ portion of the library also has the ability to detect impacts and calculate close encounter circumstances using various formulations of the B-plane [@Kizner1961; @Opik1976; @Chodas1999; @Milani1999].

The C++ functionality is then exposed to Python through a binding generated using ``pybind11`` [@pybind11]. The Python layer of ``GRSS`` uses the propagator as the foundation to compute the orbits of small bodies from a given set of optical and/or radar astrometry from the Minor Planet Center, the JPL Small Body Radar Astrometry database, and the Gaia DR3 solar system observations database. Additionally, the orbit determination modules also gave the ability to fit especially demanding measurements such as stellar occultations. These capabilities of the ``GRSS`` library have already been used to study the measurability of the heliocentric changes in the orbit of the (65803) Didymos binary asteroid system as a result of the DART impact[@Makadia2023].

**_Heartwarming conclusion that makes everyone who wants to propagate/determine small body orbits want to use GRSS goes below._**

# Acknowledgements

R.M. acknowledges support from a NASA Space Technology Graduate Research Opportunities (NSTGRO) Fellowship, contract No. 80NSSC22K1173. The work of S.R.C. and D.F. was carried out at the Jet Propulsion Laboratory, California Institute of Technology, under contract No. 80NM0018D0004 with NASA.

# References
