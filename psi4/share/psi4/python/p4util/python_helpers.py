import os
import subprocess
import re

# Relative hack for now
import sys, inspect
path_dir = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../")))
sys.path.append(path_dir)
import qcdb

import psi4

@staticmethod
def pybuild_basis(mol, key, target, fitrole='BASIS', other=None, puream=-1):
    basis = qcdb.BasisSet.pyconstruct(mol.create_psi4_string_from_molecule(),
                                      key, target, fitrole, other)
    psibasis = psi4.BasisSet.construct_from_pydict(mol, basis, puream)
    return psibasis

psi4.BasisSet.build = pybuild_basis
