"""GRSS: the Gauss-Radau Small-body Simulator"""
from . import fit
from . import prop
from . import utils

# get version from version.txt
with open(f'{utils.grss_path}/version.txt', 'r', encoding='utf-8') as f:
    __version__ = f.read().strip()

utils.initialize()
