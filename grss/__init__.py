"""GRSS: the Gauss-Radau Small-body Simulator"""
from . import prop
from . import fit
from . import utils

__all__ = ['prop', 'fit', 'utils']

# get version from version.txt
with open(f'{utils.grss_path}/version.txt', 'r', encoding='utf-8') as f:
    __version__ = f.read().strip()

utils.initialize()
