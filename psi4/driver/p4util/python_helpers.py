import os
import subprocess
import re

from psi4.driver import qcdb
from psi4 import psi4core

@staticmethod
def pybuild_basis(mol, key, target, fitrole='BASIS', other=None, puream=-1):
    basis = qcdb.BasisSet.pyconstruct(mol.create_psi4_string_from_molecule(),
                                      key, target, fitrole, other)
    psibasis = psi4core.BasisSet.construct_from_pydict(mol, basis, puream)
    return psibasis

psi4core.BasisSet.build = pybuild_basis
