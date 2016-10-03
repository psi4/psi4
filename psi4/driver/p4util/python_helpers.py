import os
import subprocess
import re

from psi4.driver import qcdb
from psi4 import psi4core

## Python basis helps

@staticmethod
def pybuild_basis(mol, key=None, target=None, fitrole='BASIS', other=None, puream=-1):
    if key is None:
        key = "ORBITAL"
    if target is None:
        target = psi4.get_global_option("BASIS")

    basis = qcdb.BasisSet.pyconstruct(mol.create_psi4_string_from_molecule(),
                                      key, target, fitrole, other)
    psibasis = psi4core.BasisSet.construct_from_pydict(mol, basis, puream)
    return psibasis

psi4core.BasisSet.build = pybuild_basis

## Python wavefunction helps

@staticmethod
def pybuild_wavefunction(mol, basis=None):
    if basis is None:
        basis = psi4core.BasisSet.build(mol)
    elif isinstance(basis, (str, unicode)):
        basis = psi4core.BasisSet.build(mol, "ORBITAL", basis)


    return psi4core.Wavefunction(mol, basis)

def delete(self):
    print('Clearing cdict')
    self.cdict.clear()

psi4core.Wavefunction.build = pybuild_wavefunction
psi4core.Wavefunction.__del__ = delete
psi4core.Wavefunction.__exit__ = delete

## Python JK helps

@staticmethod
def pybuild_JK(basis, aux=None):
    if aux is None:
        aux = psi4core.BasisSet.build(basis.molecule(), "DF_BASIS_SCF",
                                      psi4core.get_option("SCF", "DF_BASIS_SCF"),
                                      "JKFIT", psi4core.get_global_option('BASIS'),
                                      basis.has_puream())

    return psi4core.JK.build_JK(basis, aux)

psi4core.JK.build = pybuild_JK

## Python other helps

psi4core.Molecule.run_dftd3 = qcdb.interface_dftd3.run_dftd3

