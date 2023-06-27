import os
import sys
from . import prop
from . import fit

grss_path = os.path.dirname(os.path.abspath(__file__))
grss_kernel_path = f'{grss_path}/../kernels'
sys.path.append(grss_kernel_path)

def default_kernel_path(kernel_version=0):
    if kernel_version == 0:
        return f'{grss_kernel_path}/planets_big16_de441_1950_2350.tm'
    file = f'{grss_kernel_path}/planets_big16_de{kernel_version}_1950_2350.tm'
    if os.path.isfile(file):
        return file
    raise ValueError(f'Unknown kernel version: {kernel_version}')
