"""Plugin docstring.

"""
__version__ = '0.1'
__author__  = 'Eugene DePrince'

# Load Python modules
from pymodule import *

# Load C++ plugin
import os
import PsiMod
plugdir = os.path.split(os.path.abspath(__file__))[0]
sofile = plugdir + '/' + os.path.split(plugdir)[1] + '.so'
PsiMod.plugin_load(sofile)

