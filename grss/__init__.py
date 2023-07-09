"""GRSS: the Gauss-Radau Small-body Simulator"""
from . import prop
from . import fit
from . import utils

__all__ = ['prop', 'fit', 'utils']

# get version from version.txt
with open(f'{utils.grss_path}/version.txt', 'r', encoding='utf-8') as ver_file:
    __version__ = ver_file.read().strip()

utils.initialize()
