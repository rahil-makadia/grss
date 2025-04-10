# GRSS

[![PyPI](https://img.shields.io/pypi/v/grss?label=PyPI&color=[41,171,71])](https://pypi.python.org/pypi/grss/)
[![Build docs)](https://github.com/rahil-makadia/grss/actions/workflows/docs.yml/badge.svg)](https://github.com/rahil-makadia/grss/actions/workflows/docs.yml)
[![Python tests)](https://github.com/rahil-makadia/grss/actions/workflows/python_tests.yml/badge.svg)](https://github.com/rahil-makadia/grss/actions/workflows/python_tests.yml)
[![C++ tests)](https://github.com/rahil-makadia/grss/actions/workflows/cpp_tests.yml/badge.svg)](https://github.com/rahil-makadia/grss/actions/workflows/cpp_tests.yml)
[![PSJ](https://img.shields.io/badge/PSJ-10.3847/PSJ/adbc88-[41,171,71])](https://doi.org/10.3847/PSJ/adbc88)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.07861/status.svg)](https://doi.org/10.21105/joss.07861)
[![GPLv3](https://img.shields.io/badge/license-GPLv3-[41,171,71])](https://github.com/rahil-makadia/grss/blob/main/LICENSE)

**GRSS** (pronounced "grass"), the *Gauss-Radau Small-body Simulator* is a Python package with a C++ binding for propagating and fitting the orbits of small bodies in the solar system, such as asteroids and comets.

If you use GRSS in your research, please cite at least one of the following:

* Makadia et al. (2025). Gauss-Radau Small-body Simulator (GRSS): An Open-Source Library for Planetary Defense. [The Planetary Science Journal, 6(4), 85.](https://doi.org/10.3847/PSJ/adbc88)
* Makadia et al. (2025). Gauss-Radau Small-body Simulator (GRSS): An Open-Source Library for Planetary Defense. [Journal of Open Source Software, 10(107), 7861.](https://doi.org/10.21105/joss.07861)

## Getting Started

There are currently two different ways to install the GRSS library.

* Using the [Python Package Index (PyPI)](https://pypi.org/project/grss/)
* Using the source code from the [GitHub repository](https://www.github.com/rahil-makadia/grss)

### Install via PyPI

The GRSS library is available on PyPI and can be installed using the following command:

``` console
    pip install grss
```

If this installation fails (i.e., you get an error when importing GRSS), you can try installing it without using the binary wheel on PyPI by using the following command:

``` console
    pip install grss --no-binary grss
```

NOTE: The GRSS library is currently not pip-installable on Intel-based Macs. To use the library on an Intel-based Mac, please install the library using the source code from the GitHub repository (see below for instructions).

### Install via source code (Python)

The source code for the GRSS library is available on GitHub and can be downloaded using the following command:

```console
    git clone https://www.github.com/rahil-makadia/grss
```

Once the source code has been downloaded, the library can be installed using the following command:

```console
    cd grss
    source initialize.sh
    python3 -m pip install .
```

### Install via source code (C++, reduced functionality)

The source code for the GRSS library is available on GitHub and can be downloaded using the following command:

```console
    git clone https://www.github.com/rahil-makadia/grss
```

Once the source code has been downloaded, the library can be installed using the following command:

```console
    cd grss
    source initialize.sh
    source build_cpp.sh
```

You will need to have CMake installed on your system to build the C++ library. Once the build script has completed, you can use the resulting static/shared library from the `build` directory in your C++ projects.

Keep in mind the C++ library only contains support for propagating orbits and calculating observables. If you want to use the orbit fitting functionality, you will need to install the full Python library.

## Usage

Once the GRSS library has been installed, it can be imported into a Python script using the following command:

``` python
   import grss
```

The first time the library is imported, it will download some data files such as NAIF SPICE kernels and the data needed to debias optical astrometry. This should should take a few minutes. Once these files are available to the library, you are ready to use GRSS to its full potential!

Check out the [examples](https://rahil-makadia.github.io/grss/examples.html) on the [GRSS website](https://rahil-makadia.github.io/grss/) to get started.

## Contributing

If you would like to contribute to the GRSS library or raise an issue, please see the [contributing file](https://github.com/rahil-makadia/grss/blob/main/contributing.md).

## Acknowledgements

GRSS Development Team:

* Rahil Makadia
* Steven R. Chesley
* Siegfried Eggl
* Davide Farnocchia

The GRSS library was developed by Rahil Makadia as part of his PhD dissertation at the University of Illinois at Urbana-Champaign. This work was supported by a NASA Space Technologies Graduate Research Opportunities (NSTGRO) Fellowship, Grant #80NSSC22K1173. Rahil would like to thank his advisor, Dr. Siegfried Eggl as well as his collaborators, Dr. Steven R. Chesley, and Dr. Davide Farnocchia for their guidance and support.
