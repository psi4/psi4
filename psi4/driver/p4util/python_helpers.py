import os
import subprocess
import re

from psi4.driver import qcdb
from psi4 import core

## Python basis helps

@staticmethod
def pybuild_basis(mol, key=None, target=None, fitrole='BASIS', other=None, puream=-1):
    if key is None:
        key = "ORBITAL"
    if target is None:
        target = psi4.get_global_option("BASIS")

    basis = qcdb.BasisSet.pyconstruct(mol.create_psi4_string_from_molecule(),
                                      key, target, fitrole, other)
    psibasis = core.BasisSet.construct_from_pydict(mol, basis, puream)
    return psibasis

core.BasisSet.build = pybuild_basis

## Python wavefunction helps

@staticmethod
def pybuild_wavefunction(mol, basis=None):
    if basis is None:
        basis = core.BasisSet.build(mol)
    elif isinstance(basis, (str, unicode)):
        basis = core.BasisSet.build(mol, "ORBITAL", basis)


    return core.Wavefunction(mol, basis)

def delete(self):
    print('Clearing cdict')
    self.cdict.clear()

core.Wavefunction.build = pybuild_wavefunction
core.Wavefunction.__del__ = delete
core.Wavefunction.__exit__ = delete

## Python JK helps

@staticmethod
def pybuild_JK(basis, aux=None):
    if aux is None:
        aux = core.BasisSet.build(basis.molecule(), "DF_BASIS_SCF",
                                      core.get_option("SCF", "DF_BASIS_SCF"),
                                      "JKFIT", core.get_global_option('BASIS'),
                                      basis.has_puream())

    return core.JK.build_JK(basis, aux)

core.JK.build = pybuild_JK

## Python other helps

core.Molecule.run_dftd3 = qcdb.interface_dftd3.run_dftd3

