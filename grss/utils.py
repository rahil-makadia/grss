"""General utilities for the GRSS package"""
import os
import sys

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
    return None
