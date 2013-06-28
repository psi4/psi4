import sys
import os

# Add the directory containing 'grendel' to the sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir))

# Add the directory containing the 'grendel_tests' package to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from grendel.util.units import IncompatibleUnitsError, UnknownUnitError



