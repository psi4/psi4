#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2021 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#
"""Module with property-related helper functions."""

import numpy as np
import qcelemental as qcel

import psi4
from psi4 import core


def free_atom_volumes(wfn, **kwargs):
    """ 

    """

    # the level of theory
    module = wfn.module()
    # if we're doing dft, grab the functional
    if module.lower() == 'scf':
        module = wfn.functional().name()
    # list of reference number of unpaired electrons.
    # Note that this is not the same as the 
    # total spin of the ground state atom
    reference_S = [ 0,
                    1,                                                                                           0,
                    1, 0,                                                                         1, 2, 3, 2, 1, 0,
                    1, 0,                                                                         1, 2, 3, 2, 1, 0,
                    1, 0,                                           1, 2, 3, 6, 5, 4, 3, 2, 1, 0, 1, 2, 3, 2, 1, 0,
                    1, 0,                                           1, 2, 5, 6, 5, 4, 3, 0, 1, 0, 1, 2, 3, 2, 1, 0,
                    1, 0, 1, 0, 3, 4, 5, 6, 7, 8, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0, 1, 2, 3, 2, 1, 0 ]





    # A dictionary to map elements to volumes
    # Volumes computed using MBIS with UHF 
    # and the same basis set as the molecular computation
    volumes = {}

    # the parent molecule
    mol = wfn.molecule()
    basis = mol.basis_on_atom(0)

    # Get unique atoms by input symbol
    symbols = {}
    natom = mol.natom()
    for atom in range(natom):
        symbols[mol.symbol(atom)] = int(round(mol.Z(atom)))

    core.print_out(f"  Running {len(symbols)} free-atom UHF computations")

    for a_sym, a_z in symbols.items():

        geom = f"""
0 {int(1+reference_S[a_z])} 
{a_sym} 0.0 0.0 0.0
symmetry c1
"""
    
        # make sure we do UHF/UKS if we're not a singlet
        if reference_S[a_z] != 0:
            psi4.set_options({"reference":"UHF"}) 

        # Set the molecule, here just an atom
        molrec = qcel.molparse.from_string(geom, enable_qm=True, 
        missing_enabled_return_qm='minimal', enable_efp=True, missing_enabled_return_efp='none')
        a_mol = core.Molecule.from_dict(molrec['qm'])
        a_mol.update_geometry() 
        psi4.molutil.activate(a_mol)

        method = module+"/"+basis

        # Get the atomic wfn
        at_e, at_wfn = psi4.energy(method, return_wfn=True)

        # Now, re-run mbis for the atomic density, grabbing only the volume 
        psi4.oeprop(at_wfn, 'MBIS_CHARGES', title=a_sym + " " + method , free_atom=True) 
        vw = at_wfn.array_variable('MBIS RADIAL MOMENTS <R^3>')
        vw = vw.get(0,0) 

        # set the atomic widths as wfn variables
        wfn.set_variable("MBIS FREE ATOM " + a_sym + " VOLUME", vw)


    # reset mol to original
    mol.update_geometry()

    

