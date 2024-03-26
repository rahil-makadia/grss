Getting Started
===============
There are currently two different ways to install the GRSS library.

* Using the `Python Package Index (PyPI) <https://pypi.org/project/grss/>`_
* Using the source code from the `GRSS GitHub repository <https://www.github.com/rahil-makadia/grss>`_

----------------
Install via PyPI
----------------
The GRSS library is available on PyPI and can be installed using the following command:

.. code-block:: console

   pip install grss

If this installation fails (i.e., you get an error when importing GRSS), you can try installing it without using the binary wheel on PyPI by using the following command:

.. code-block:: console

   pip install grss --no-binary grss

NOTE: The GRSS library is currently not pip-installable on Intel-based Macs. To use the library on an Intel-based Mac, please install the library using the source code from the GitHub repository (see below for instructions).

--------------------------------
Install via source code (Python)
--------------------------------
The source code for the GRSS library is available on GitHub and can be downloaded using the following command:

.. code-block:: console

   git clone https://www.github.com/rahil-makadia/grss

Once the source code has been downloaded, the library can be installed using the following commands:

.. code-block:: console

   cd grss
   source initialize.sh
   python3 -m pip install .

----------------------------------------------------
Install via source code (C++, reduced functionality)
----------------------------------------------------
The source code for the GRSS library is available on GitHub and can be downloaded using the following command:

.. code-block:: console

   git clone https://www.github.com/rahil-makadia/grss

Once the source code has been downloaded, the library can be installed using the following commands:

.. code-block:: console

   cd grss
   source initialize.sh
   source build_cpp.sh

You will need to have CMake installed on your system to build the C++ library. Once the build script has completed, you can use the resulting static/shared library from the `build` directory in your C++ projects.

Keep in mind the C++ library only contains support for propagating orbits and calculating observables. If you want to use the orbit fitting functionality, you will need to install the full Python library.

Usage
=====
Once the GRSS library has been installed, it can be imported into a Python script using the following command:

.. code-block:: python

   import grss

The first time the library is imported, it will download some data files such as NAIF SPICE kernels and the data needed to debias optical astrometry. This should should take a few minutes, and if the download was completed, the following message will be printed for each file:

.. code-block:: console

   YYYY-MM-DD HH:MM:SS URL:url-of-downloaded-file [filesize] -> path/to/downloaded/file [1]

Once these files are available to the library, you are ready to use GRSS to its full potential!
