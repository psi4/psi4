#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2025 The Psi4 Developers.
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

import psi4

from .exceptions import ValidationError
from . import optproc

__all__ = ['free_atom_volumes']


def free_atom_volumes(wfn: psi4.core.Wavefunction, **kwargs):
    """ 
    Computes free-atom volumes using MBIS density partitioning.
    The free-atom volumes are computed for all unique (inc. basis set)
    atoms in a molecule and stored as wavefunction variables, :psivar:`MBIS FREE ATOM n VOLUME`.
    Free-atom densities are computed at the same level of theory as the molecule, 
    and we use unrestricted references as needed in computing the ground-state. 

    The free-atom volumes are used to compute volume ratios in routine MBIS computations

    Parameters
    ----------
    wfn
        The wave function associated with the molecule, method, and basis for 
        atomic computations
    """

    # If we're already a free atom, break to avoid recursion
    # We don't ever need volume ratios for free atoms since they
    # are by definition 1.0
    natom = wfn.molecule().natom()
    if natom == 1:
        return 0 
    
    # We need to know the level of theory of the system to compute the free atoms
    # This isn't stored or might have been preloaded so we search the psi4 variables for
    # the best match
    current_en = wfn.scalar_variable('CURRENT ENERGY')
    total_keys = [
        k for k in wfn.scalar_variables().keys() if
        ('TOTAL ENERGY' in k and 'SCF' not in k)
    ]
    total_energy_diffs = sorted([
        [abs(wfn.scalar_variable(k) - current_en), k] for k in total_keys],
        key=lambda x: x[0]
    )
    if len(total_energy_diffs) == 0 or total_energy_diffs[0][0] > 1e-8:
        raise ValidationError(
            "No valid 'method TOTAL ENERGY' found in wavefunction scalar variables. Needed for MBIS free atoms."
        )

    theory = total_energy_diffs[0][1].split()[0]
    if theory == 'DFT':
        theory = wfn.functional().name()

    # list of reference number of unpaired electrons.
    # Note that this is not the same as the
    # total spin of the ground state atom
    reference_S = [
        0, 1, 0, 1, 0, 1, 2, 3, 2, 1, 0, 1, 0, 1, 2, 3, 2, 1, 0, 1, 0, 1, 2, 3, 6, 5, 4, 3, 2, 1, 0, 1, 2, 3, 2, 1, 0,
        1, 0, 1, 2, 5, 6, 5, 4, 3, 0, 1, 0, 1, 2, 3, 2, 1, 0, 1, 0, 1, 0, 3, 4, 5, 6, 7, 8, 5, 4, 3, 2, 1, 0, 1, 2, 3,
        4, 5, 4, 3, 2, 1, 0, 1, 2, 3, 2, 1, 0
    ]

    # the parent molecule and reference type
    mol = wfn.molecule()

    # Get unique atoms by input symbol,
    # Be to handle different basis sets
    unq_atoms = set()
    for atom in range(mol.natom()):
        symbol = mol.symbol(atom)
        Z = int(mol.Z(atom))
        basis = mol.basis_on_atom(atom)
        unq_atoms.add((symbol, Z, basis))

    psi4.core.print_out(f"  Running {len(unq_atoms)} free-atom UHF computations")

    optstash = optproc.OptionsState(["SCF", 'REFERENCE'])
    for a_sym, a_z, basis in unq_atoms:

        # make sure we do UHF/UKS if we're not a singlet
        if reference_S[a_z] != 0:
            psi4.core.set_local_option("SCF", "REFERENCE", "UHF")
        else:
            psi4.core.set_local_option("SCF", "REFERENCE", "RHF")

        # Set the molecule, here just an atom
        a_mol = psi4.core.Molecule.from_arrays(geom=[0, 0, 0],
                                          elem=[a_sym],
                                          molecular_charge=0,
                                          molecular_multiplicity=int(1 + reference_S[a_z]))
        a_mol.update_geometry()
        psi4.molutil.activate(a_mol)

        method = theory + "/" + basis

        # Get the atomic wfn
        at_e, at_wfn = psi4.energy(method, return_wfn=True)

        # Now, re-run mbis for the atomic density, grabbing only the volume
        psi4.oeprop(at_wfn, 'MBIS_CHARGES', title=a_sym + " " + method, free_atom=True)

        vw = at_wfn.array_variable('MBIS RADIAL MOMENTS <R^3>')  # P::e OEPROP
        vw = vw.get(0, 0)

        # set the atomic widths as wfn variables
        wfn.set_variable("MBIS FREE ATOM " + a_sym.upper() + " VOLUME", vw)
        # set_variable("MBIS FREE ATOM n VOLUME")  # P::e OEPROP
        
        psi4.core.clean()
        psi4.core.clean_variables()

    # reset mol and reference to original
    optstash.restore()
    mol.update_geometry()
    psi4.molutil.activate(mol)
