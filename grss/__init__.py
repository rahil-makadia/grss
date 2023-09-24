"""GRSS: the Gauss-Radau Small-body Simulator"""
from . import fit
from . import prop
from . import utils

# get version from version.txt
__version__ = open(f'{utils.grss_path}/version.txt', 'r', encoding='utf-8').read().strip()

utils.initialize()
