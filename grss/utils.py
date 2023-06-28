import os
import sys

__all__ = [ 'default_kernel_path',
            'grss_path',
            'grss_kernel_path',
]

grss_path = os.path.dirname(os.path.abspath(__file__))
grss_kernel_path = f'{grss_path}/kernels'

sys.path.append(grss_kernel_path)

def default_kernel_path(kernel_version=0):
    if kernel_version == 0:
        return f'{grss_kernel_path}/planets_big16_de441_1950_2350.tm'
    file = f'{grss_kernel_path}/planets_big16_de{kernel_version}_1950_2350.tm'
    if os.path.isfile(file):
        return file
    raise ValueError(f'Unknown kernel version: {kernel_version}')

def initialize():
    # get the file directory
    cwd = os.path.dirname(os.path.abspath(__file__))
    # run the get_debiasing_data.py script
    os.system(f'python {cwd}/debias/get_debiasing_data.py')
    # run the get_kernels.py script
    os.system(f'python {cwd}/kernels/update_kernels.py')
    return None
