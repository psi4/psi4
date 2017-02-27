#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
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
# @END LICENSE
#

import os
import subprocess
import re
import sys

from . import optproc
from psi4.driver import qcdb
from psi4 import core

## Python basis helps

@staticmethod
def pybuild_basis(mol, key=None, target=None, fitrole='BASIS', other=None, puream=-1, return_atomlist=False):
    horde = qcdb.libmintsbasisset.basishorde

    if key == 'ORBITAL':
        key = 'BASIS'

    if horde and key:
        tmp = horde.get(core.get_global_option(key), None)
        if tmp:
            target = tmp
        elif target:
            pass
        elif tmp is None:
            target = None
    elif target:
        pass
    elif key is None:
        target = core.get_global_option("BASIS")
    else:
        target = core.get_global_option(key)

    basisdict = qcdb.BasisSet.pyconstruct(mol.create_psi4_string_from_molecule(),
                                          key, target, fitrole, other, return_atomlist=return_atomlist)

    if return_atomlist:
        atom_basis_list = []
        for atbs in basisdict:
            atommol = core.Molecule.create_molecule_from_string(atbs['molecule'])
            lmbs = core.BasisSet.construct_from_pydict(atommol, atbs, puream)
            atom_basis_list.append(lmbs)
            #lmbs.print_detail_out()
        return atom_basis_list

    psibasis = core.BasisSet.construct_from_pydict(mol, basisdict, puream)
    return psibasis

core.BasisSet.build = pybuild_basis

## Python wavefunction helps

@staticmethod
def pybuild_wavefunction(mol, basis=None):
    if basis is None:
        basis = core.BasisSet.build(mol)
    elif (sys.version_info[0] == 2) and isinstance(basis, (str, unicode)):
        basis = core.BasisSet.build(mol, "ORBITAL", basis)
    elif (sys.version_info[0] > 2) and isinstance(basis, str):
        basis = core.BasisSet.build(mol, "ORBITAL", basis)


    return core.Wavefunction(mol, basis)

core.Wavefunction.build = pybuild_wavefunction

## Python JK helps

@staticmethod
def pybuild_JK(orbital_basis, aux=None, jk_type=None):
    """
    Constructs a Psi4 JK object from an input basis.

    Parameters
    ----------
    orbital_basis : :py:class:`~psi4.core.BasisSet`
        Orbital basis to use in the JK object.
    aux : :py:class:`~psi4.core.BasisSet`
        Optional auxiliary basis set for density-fitted tensors. Defaults
        to the DF_BASIS_SCF if set, otherwise the correspond JKFIT basis
        to the passed in orbital_basis.
    type : str
        Type of JK object to build (DF, Direct, PK, etc). Defaults to the
        current global SCF_TYPE option.

    Returns
    -------
    :py:class:`~psi4.core.JK`
        Uninitialized JK object.

    Example
    -------

    jk = psi4.core.JK.build(bas)
    jk.set_memory(int(5e8)) # 4GB of memory
    jk.initialize()

    ...

    jk.C_left_add(matirx)
    jk.compute()
    jk.C_clear()

    ...

    """

    optstash = optproc.OptionsState(["SCF_TYPE"])

    if jk_type is not None:
        core.set_global_option("SCF_TYPE", jk_type)

    if aux is None:
        if core.get_option("SCF", "SCF_TYPE") == "DF":
            aux = core.BasisSet.build(orbital_basis.molecule(), "DF_BASIS_SCF",
                                      core.get_option("SCF", "DF_BASIS_SCF"),
                                      "JKFIT", core.get_global_option('BASIS'),
                                      orbital_basis.has_puream())
        else:
            aux = core.BasisSet.zero_ao_basis_set()


    jk = core.JK.build_JK(orbital_basis, aux)

    optstash.restore()
    return jk

core.JK.build = pybuild_JK

## Python other helps

core.Molecule.run_dftd3 = qcdb.interface_dftd3.run_dftd3
core.Molecule.run_gcp = qcdb.interface_gcp.run_gcp


def set_options(options_dict):
    """
    Sets Psi4 global options from an input dictionary.
    """

    for k, v, in options_dict.items():
        core.set_global_option(k.upper(), v)


def set_module_options(module, options_dict):
    """
    Sets Psi4 module options from a module specification and input dictionary.
    """

    for k, v, in options_dict.items():
        core.set_local_option(module.upper(), k.upper(), v)


def pcm_helper(block):
    """Passes multiline string *block* to PCMSolver parser."""

    with open('pcmsolver.inp', 'w') as handle:
        handle.write(block)
    import pcmsolver
    pcmsolver.parse_pcm_input('pcmsolver.inp')
