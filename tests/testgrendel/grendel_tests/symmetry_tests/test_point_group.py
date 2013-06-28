import unittest
import sys
import os

try:
    profile
except NameError:
    def profile(f):
        return f

# Add the directory containing 'grendel' to the sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir))

# Add the directory containing the 'grendel_tests' package to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from grendel.chemistry import SampleMolecules
from grendel_tests import allow_generators, long_test, skip

@allow_generators
class PointGroupTest(unittest.TestCase):

    molecules = {
        'water'   : 'C_2v',
        'methane' : 'T_d',
        'benzene' : 'D_6h',
        'ethylene': 'D_2h'
    }

    @classmethod
    def symmetry_test_generator(cls):
        def _generated(self, molname, expected_pg_name):
            self.assertEqual(str(SampleMolecules[molname].point_group), expected_pg_name)
        for name in cls.molecules:
            _generated.__name__ = "test_" + name + "_point_group"
            if name in SampleMolecules:
                yield profile(long_test(_generated)), name, cls.molecules[name]
            else:
                yield skip("Sample geometry for " + name + " not available.")(_generated), name, cls.molecules[name]

