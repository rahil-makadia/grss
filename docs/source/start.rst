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

-----------------------
Install via source code
-----------------------
The source code for the GRSS library is available on GitHub and can be downloaded using the following command:

.. code-block:: console

   git clone https://www.github.com/rahil-makadia/grss

Once the source code has been downloaded, the library can be installed using the following commands:

.. code-block:: console

   source initialize.sh
   python -m pip install .

-----
Usage
-----
Once the GRSS library has been installed, it can be imported into a Python script using the following command:

.. code-block:: python

   import grss

The first time the library is imported, it will download some data files such as NAIF SPICE kernels and the data needed to debias optical astrometry. This should not take more than a couple minutes, and if the download was completed, the following message will be printed:

.. code-block:: console

   YYYY-MM-DD HH:MM:SS URL:url-of-downloaded-file [filesize] -> path/to/downloaded/file [1]

Once these files are available to the library, you are ready to use GRSS to its full potential!
