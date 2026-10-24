"""GRSS Python orbit determination subpackage"""
import warnings
from . import fit_ades, fit_iod, fit_optical, fit_radar, fit_simulation, fit_utils
from .fit_ades import *
from .fit_iod import *
from .fit_optical import *
from .fit_radar import *
from .fit_simulation import *
from .fit_utils import *
__all__ = list(dict.fromkeys(
    name for module in (fit_ades, fit_iod, fit_optical, fit_radar, fit_simulation, fit_utils)
    for name in module.__all__
))
warnings.filterwarnings(action='ignore', module='erfa') # for the astropy time warning
