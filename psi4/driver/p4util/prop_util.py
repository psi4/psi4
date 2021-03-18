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

import qcelemental as qcel

import psi4
from psi4 import core
from . import optproc


def free_atom_volumes(wfn, **kwargs):
    """ 

    """

    # print level
    print_level = core.get_global_option("PRINT")

    # the level of theory
    current_en = wfn.scalar_variable('CURRENT ENERGY')
    test = [k for k,v in wfn.scalar_variables().items() if abs(v-current_en) <= 1e-12]     
    theory = ""
    for var in test:
        if 'TOTAL ENERGY' in var:
            var = var.split()
            if var[0] == 'SCF':
                continue
            elif var[0] == 'DFT':
                theory =  wfn.functional().name()
            else:
                theory = var[0]

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


    # the parent molecule and reference type
    mol = wfn.molecule()

    # Get unique atoms by input symbol,
    # Be to handle different basis sets
    unq_atoms = set()
    for atom in range(mol.natom()):
        symbol = mol.symbol(atom)
        Z = int(mol.Z(atom))
        basis = mol.basis_on_atom(atom)
        unq_atoms.add((symbol,Z,basis))

    core.print_out(f"  Running {len(unq_atoms)} free-atom UHF computations")

    optstash = optproc.OptionsState(['REFERENCE'])
    for a_sym, a_z, basis in unq_atoms:

        geom = f"""
0 {int(1+reference_S[a_z])} 
{a_sym} 0.0 0.0 0.0
symmetry c1
"""
        
        # make sure we do UHF/UKS if we're not a singlet
        if reference_S[a_z] != 0:
            core.set_global_option("REFERENCE", "UHF") 
        else:
            core.set_global_option("REFERENCE", "RHF") 

        # Set the molecule, here just an atom
        a_mol = core.Molecule.from_arrays(geom=[0,0,0], elem=[a_sym], molecular_charge=0, 
                                          molecular_multiplicity=int(1+reference_S[a_z]))
        a_mol.update_geometry() 
        psi4.molutil.activate(a_mol)

        method = theory + "/" + basis

        # Supress printing
        if print_level <= 1:
            core.be_quiet()

        # Get the atomic wfn
        at_e, at_wfn = psi4.energy(method, return_wfn=True)

        # Now, re-run mbis for the atomic density, grabbing only the volume 
        psi4.oeprop(at_wfn, 'MBIS_CHARGES', title=a_sym + " " + method , free_atom=True) 

        if print_level <= 1:
            core.reopen_outfile()

        vw = at_wfn.array_variable('MBIS RADIAL MOMENTS <R^3>')
        vw = vw.get(0,0) 

        # set the atomic widths as wfn variables
        wfn.set_variable("MBIS FREE ATOM " + a_sym.upper() + " VOLUME", vw)


    # reset mol and reference to original
    optstash.restore()
    mol.update_geometry()
    psi4.molutil.activate(mol) 

