#
#@BEGIN LICENSE
#
# PSI4: an ab initio quantum chemistry software package
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#@END LICENSE
#

"""Module to facilitate quantum chemical computations on chemical
databases. Contains Molecule class and physical constants from psi4 suite.

"""
__version__ = '0.3'
__author__ = 'Lori A. Burns'

# Load Generic Python Modules
try:
    from collections import defaultdict, OrderedDict
except ImportError:
    from collections import defaultdict
    from oldpymodules import OrderedDict

# Load Python modules
from molecule import Molecule
from dbproc import *
from options import *
from qcformat import *
import cfour
import jajo
from orient import OrientMols
#from origdbwrap import Database, FourDatabases, ThreeDatabases  # , ReactionDatum, Reagent, Reaction
from dbwrap import Database, DB4 #DatabaseWrapper  #ReactionDatum, Reagent, Reaction
from libmintspointgrp import SymmetryOperation, PointGroup
from libmintsbasisset import BasisSet

# Load items that are useful to access from an input file
from psiutil import *
from physconst import *
