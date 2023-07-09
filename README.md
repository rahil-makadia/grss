# grss

[![PyPi Version](https://img.shields.io/pypi/v/grss.svg)](https://pypi.python.org/pypi/grss/)
[![Build Sphinx docs)](https://github.com/rahil-makadia/grss/actions/workflows/docs.yml/badge.svg)](https://github.com/rahil-makadia/grss/actions/workflows/docs.yml)
[![Python tests)](https://github.com/rahil-makadia/grss/actions/workflows/python_tests.yml/badge.svg)](https://github.com/rahil-makadia/grss/actions/workflows/python_tests.yml)
[![C++ tests)](https://github.com/rahil-makadia/grss/actions/workflows/cpp_tests.yml/badge.svg)](https://github.com/rahil-makadia/grss/actions/workflows/cpp_tests.yml)
[![GPL](https://img.shields.io/badge/license-GPL-green.svg?style=flat)](https://github.com/rahil-makadia/grss/blob/main/LICENSE)

**GRSS** (pronounced "grass"), the *Gauss-Radau Small-body Simulator* is a Python package with a C++ binding for propagating and fitting the orbits of small bodies in the solar system, such as asteroids and comets.

## Getting Started

There are currently two different ways to install the GRSS library.

* Using the [Python Package Index (PyPI)](https://pypi.org/project/grss/)
* Using the source code from the [GitHub repository](https://www.github.com/rahil-makadia/grss)

### Install via PyPI

The GRSS library is available on PyPI and can be installed using the following command:

``` console
    pip install grss
```

### Install via source code

The source code for the GRSS library is available on GitHub and can be downloaded using the following command:

```console
    git clone https://www.github.com/rahil-makadia/grss
```

Once the source code has been downloaded, the library can be installed using the following command:

```console
    python setup.py install
```

## Usage

Once the GRSS library has been installed, it can be imported into a Python script using the following command:

``` python
   import grss
```

The first time the library is imported, it will download some data files such as NAIF SPICE kernels and the data needed to debias optical astrometry. This should not take more than a couple minutes, and if the download was completed, the following message will be printed:

```console
   YYYY-MM-DD HH:MM:SS URL:url-of-downloaded-file [filesize] -> path/to/downloaded/file [1]
```

Once these files are available to the library, you are ready to use GRSS to its full potential!

Check out the [examples](https://rahil-makadia.github.io/grss/examples.html) on the [GRSS website](https://rahil-makadia.github.io/grss/) to get started.

## Acknowledgements

GRSS Development Team:

* Rahil Makadia
* Steven R. Chesley
* Siegfried Eggl
* Davide Farnocchia

The GRSS library was developed by Rahil Makadia as part of his PhD dissertation at the University of Illinois at Urbana-Champaign. This work was supported by a NASA Space Technologies Graduate Research Opportunities (NSTGRO) Fellowship, Grant #80NSSC22K1173. The author would like to thank his advisor, Dr. Siegfried Eggl as well as his collaborators, Dr. Steven R. Chesley, and Dr. Davide Farnocchia for their guidance and support.
