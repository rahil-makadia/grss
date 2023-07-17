"""General utilities for the GRSS package"""
import os
import sys
import datetime
from requests import request

__all__ = [ 'default_kernel_path',
            'grss_path',
            'grss_kernel_path',
            'initialize'
]

grss_path = os.path.dirname(os.path.abspath(__file__))
grss_kernel_path = f'{grss_path}/kernels'

sys.path.append(grss_kernel_path)

def default_kernel_path(kernel_version=0):
    """
    Returns the paths to the default SPICE kernels shipped with GRSS.

    Parameters
    ----------
    kernel_version : int, optional
        SPICE kernel version. Optional, defaults to 0 (DE441)

    Returns
    -------
    path : str
        Full path to the SPICE kernel

    Raises
    ------
    ValueError
        When the kernel version is not known
    """
    if kernel_version == 0:
        return f'{grss_kernel_path}/planets_big16_de441_1950_2350.tm'
    file = f'{grss_kernel_path}/planets_big16_de{kernel_version}_1950_2350.tm'
    if os.path.isfile(file):
        return file
    raise ValueError(f'Unknown kernel version: {kernel_version}')

def _download_obs_codes_file():
    """
    Download a text file of the observatory code info from the Minor Planet
    Center.
    """
    url = "https://www.minorplanetcenter.net/iau/lists/ObsCodes.html"
    req = request("GET", url, timeout=30)
    text = req.content.replace(b'<pre>',b'').replace(b'</pre>',b'')

    with open(f'{grss_path}/obs_codes.txt', "wb") as file:
        now = datetime.datetime.now()
        file.write((f'{now.year}, {now.month}, {now.day}, '
                    f'{now.hour}, {now.minute}, {now.second}\n'.encode()))
        file.write(text[1:-2])
    return None

def get_obs_codes_file():
    # sourcery skip: extract-method
    """
    Download a text file of the observatory code info from the Minor Planet
    Center. If it already exists, update it if the last download was more than 1
    week ago.
    """
    # check if file exists
    if os.path.isfile(f'{grss_path}/obs_codes.txt'):
        with open(f'{grss_path}/obs_codes.txt', 'r', encoding='utf-8') as file:
            now = datetime.datetime.now()
            data = file.read()
            last_updated = data.split('\n')[0]
            last_updated = datetime.datetime(*[int(i) for i in last_updated.split(',')])
            diff_weeks = (now - last_updated).seconds / (60 * 60 * 24 * 7)
            if diff_weeks >= 1.0:
                _download_obs_codes_file()
    else:
        _download_obs_codes_file()
    return None

def initialize():
    """
    GRSS library initialization function. This function will download the
    necessary data files for the GRSS library to work that aren't included
    in the package distribution (e.g. SPICE kernels, debiasing data).

    Returns
    -------
    None : NoneType
        None
    """
    # get the file directory
    cwd = os.path.dirname(os.path.abspath(__file__))
    # run the get_debiasing_data.py script
    os.system(f'python {cwd}/debias/get_debiasing_data.py')
    # run the get_kernels.py script
    os.system(f'python {cwd}/kernels/get_kernels.py')
    # get the observatory codes file
    get_obs_codes_file()
    return None
