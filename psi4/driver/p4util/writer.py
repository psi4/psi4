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
"""Module for writing input files to external codes."""

from psi4 import core
from psi4.driver import constants


def _write_nbo(self, name):
    basisset = self.basisset()
    mints = core.MintsHelper(basisset)
    mol = self.molecule()

    # Populate header and coordinates.
    NBO_file = f" $GENNBO NATOMS = {mol.natom()} NBAS = {basisset.nbf()} BODM "
    if self.nalpha() != self.nbeta():
        NBO_file += f" OPEN"
    NBO_file += " $END\n $NBO       $END\n $COORD\n"
    NBO_file += " GENNBO expects one comment line here. So, here's a comment line.\n"
    for atom in range(mol.natom()):
        NBO_file += f"{mol.true_atomic_number(atom):2d}  {int(mol.Z(atom)):2d}  {constants.bohr2angstroms * mol.x(atom):20.12f} {constants.bohr2angstroms * mol.y(atom):20.12f} {constants.bohr2angstroms * mol.z(atom):20.12f}\n"
    NBO_file += " $END\n"

    # Populate basis function information.
    pure_order = [
        [1],  # s
        [103, 101, 102],  # p
        [255, 252, 253, 254, 251],  # d: z2 xz yz x2-y2 xy
        [351, 352, 353, 354, 355, 356, 357],  # f
        [451, 452, 453, 454, 455, 456, 457, 458, 459],  #g
        [551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561]  #h
    ]
    # For historical reasons, the code loops over shells first and then basis functions within the shell.
    # This turns out to not give the same ordering as looping over basis functions directly.
    NBO_file += " $BASIS\n"
    center_string = ""
    label_string = ""
    count = 0
    for i in range(basisset.nshell()):
        shell = basisset.shell(i)
        am = shell.am
        for j in range(shell.nfunction):
            if not (count % 10):
                center_string += "\n  CENTER =" if not i else "\n          "
                label_string += "\n   LABEL =" if not i else "\n          "
            center_string += f" {shell.ncenter + 1:6d}"
            if basisset.has_puream():
                label = pure_order[am][j]
            else:
                label = 100 * am + j + 1
            label_string += f" {label:6d}"
            count += 1
    NBO_file += center_string + label_string + "\n $END\n"

    # Populate contraction information.Start with exponents.
    NBO_file += f" $CONTRACT\n  NSHELL = {basisset.nshell():6d}\n    NEXP = {basisset.nprimitive():6d}\n"
    function_nums = ""
    prim_nums = ""
    prim_indices = ""
    exponents = ""
    prim_index = 0
    coefficients = []  # [(AM int, coefficient), ...]
    for i in range(basisset.nshell()):
        if not (i % 10):
            function_nums += "\n   NCOMP =" if not i else "\n          "
            prim_nums += "\n   NPRIM =" if not i else "\n          "
            prim_indices += "\n    NPTR =" if not i else "\n          "
        shell = basisset.shell(i)
        nprim = shell.nprimitive
        function_nums += f" {shell.nfunction:6d}"
        prim_nums += f" {nprim:6d}"
        prim_indices += f" {prim_index + 1:6d}"
        for j in range(nprim):
            if not (prim_index % 4):
                exponents += "\n     EXP =" if not prim_index else "\n          "
            exponents += f"{shell.exp(j):15.6E}"
            prim_index += 1
            coefficients.append((shell.am, shell.coef(j)))
    NBO_file += function_nums + prim_nums + prim_indices + exponents

    # opulate contracton coefficients.Because some basis sets(Poples with S and P) use the same
    # coefficients for multiple angular momenta, we must supply coefficients for all primitives, for all
    # angular momenta.This leads to many zero elements.
    am_labels = ["S", "P", "D", "F", "G", "H"]
    for current_nbo_section_am in range(basisset.max_am() + 1):
        for i, (shell_am, coefficient) in enumerate(coefficients):
            if not (i % 4):
                NBO_file += f"\n      C{am_labels[current_nbo_section_am]} =" if not i else "\n          "
            if shell_am != current_nbo_section_am:
                coefficient = 0
            NBO_file += f"{coefficient:15.6E}"
    NBO_file += "\n $END"

    # That finishes most of the basis information.Next is the overlap.It would be great if we could just dump Psi's AO
    # overlap matrix, but we can 't. Per CCA guidelines, Psi' s Cartesian d and higher AM AOs aren't normalized to 1.
    # While NBO can "fix" this itself, it changes other AO quantities to match and gets the Fock matrix wrong.
    # Let's normalize ourselves instead.
    ao_overlap = mints.ao_overlap().np
    nbf = ao_overlap.shape[0]
    ao_normalizer = ao_overlap.diagonal()**(-1 / 2)

    def normalize(matrix, normalizer):
        return ((matrix * normalizer).T * normalizer).T

    normalized_ao_overlap = normalize(ao_overlap, ao_normalizer)

    def write_ao_quantity(*args):
        string = ""
        count = 0
        for quantity in args:
            for i in range(nbf):
                for j in range(nbf):
                    if not (count % 5):
                        string += "\n  "
                    string += f"{quantity[i][j]:15.6E}"
                    count += 1
        return string

    NBO_file += "\n $OVERLAP"
    NBO_file += write_ao_quantity(normalized_ao_overlap)
    NBO_file += "\n $END"

    normalized_alpha_density = normalize(self.Da_subset("AO"), 1 / ao_normalizer)
    normalized_beta_density = normalize(self.Db_subset("AO"), 1 / ao_normalizer)
    normalized_alpha_fock = normalize(self.Fa_subset("AO"), ao_normalizer)

    NBO_file += "\n $DENSITY"
    if self.same_a_b_dens():
        density = normalized_alpha_density + normalized_beta_density
        NBO_file += write_ao_quantity(density)
    else:
        NBO_file += write_ao_quantity(normalized_alpha_density, normalized_beta_density)
    NBO_file += "\n $END"

    NBO_file += "\n $FOCK"
    if not self.same_a_b_dens():
        normalized_beta_fock = normalize(self.Fb_subset("AO"), ao_normalizer)
        NBO_file += write_ao_quantity(normalized_alpha_fock, normalized_beta_fock)
    else:
        NBO_file += write_ao_quantity(normalized_alpha_fock)
    NBO_file += "\n $END"

    # The last step is to write the MO coefficients.
    NBO_file += "\n $LCAOMO"

    def write_C_matrix(C, count):
        # The C coefficients supplied the missing multiplication by the ao_normalizer in the overlap matrix before.
        # For NBO, we need that multiplication gone.
        C = (C.np.T / ao_normalizer).T
        string = ""
        for i in range(self.nmo()):
            for mu in range(nbf):
                count += 1
                if (count % 5 == 1):
                    string += ("\n  ")
                string += f"{C[mu][i]:15.6E}"
        # Pad linear dependencies
        for i in range((nbf - self.nmo()) * nbf):
            count += 1
            if (count % 5 == 1):
                string += ("\n  ")
            string += f"{0:15.6E}"
        return count, string

    count, alpha_LCAOMO = write_C_matrix(self.Ca_subset("AO", "ALL"), 0)
    NBO_file += alpha_LCAOMO
    if not self.same_a_b_orbs():
        NBO_file += write_C_matrix(self.Cb_subset("AO", "ALL"), count)[1]
    NBO_file += "\n $END\n"

    #Now time to write !
    with open(name, 'w') as f:
        f.write(NBO_file)

core.Wavefunction.write_nbo = _write_nbo
