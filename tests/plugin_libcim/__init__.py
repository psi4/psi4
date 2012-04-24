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
#   helper localcc.so file (rel path already had to be hand-coded in Makefile)
sofile = os.path.split(plugdir)[0] + '/plugin_ccsd_serial/plugin_ccsd_serial.so'
PsiMod.plugin_load(sofile)
sofile = os.path.split(plugdir)[0] + '/plugin_cepa/plugin_cepa.so'
PsiMod.plugin_load(sofile)
#   libcim.so file
sofile = plugdir + '/' + os.path.split(plugdir)[1] + '.so'
PsiMod.plugin_load(sofile)

