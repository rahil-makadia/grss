"""GRSS: the Gauss-Radau Small-body Simulator"""
from . import utils
# get version from version.txt
with open(f'{utils.grss_path}/version.txt', 'r', encoding='utf-8') as _version_file:
    __version__ = _version_file.read().strip()
utils.initialize()

from . import fit
from . import prop
from . import kernels
from . import debias
from . import libgrss

__all__ = ['fit', 'prop', 'utils', 'kernels', 'debias', 'libgrss', '__version__']
