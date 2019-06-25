#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
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

import numpy as np
import cppe

from psi4 import core


def psi4mol_to_cppemol(psi4mol):
    mol = cppe.Molecule()
    geom = psi4mol.geometry().np
    for i, c in enumerate(geom):
        a = cppe.Atom((int(psi4mol.Z(i))), *c)
        mol.append(a)
    return mol


class CppeInterface:
    def __init__(self, molecule, options, basisset):
        # setup the initial CppeState
        self.molecule = molecule
        self.options = options
        self.basisset = basisset
        self.mints = core.MintsHelper(self.basisset)
        self.cppe_state = cppe.CppeState(
            self.options, psi4mol_to_cppemol(self.molecule)
        )
        # obtain coordinates of polarizable sites
        self.polarizable_coords = np.array([
            site.position for site in self.cppe_state.get_potentials()
            if site.is_polarizable
        ])
        self.V_es = None

    def get_pe_contribution(self, density_matrix, elec_only=False):
        # build electrostatics operator
        if self.V_es is None and not elec_only:
            self.build_electrostatics_operator()
        e_electrostatic = np.sum(density_matrix * self.V_es)
        self.cppe_state.energies["Electrostatic"]["Electronic"] = e_electrostatic

        # obtain expectation values of elec. field at polarizable sites
        # elec_fields = integral_library.electric_field_value(
        #     self.polarizable_coords, density_matrix
        # )
        # solve induced moments
        # self.cppe_state.update_induced_moments(elec_fields)
        # induced_moments = self.cppe_state.induced_moments

        # build induction operator
        # V_ind = np.zeros_like(self.V_es)
        # for coord, ind_mom in zip(self.polarizable_coords, induced_moments):
        #     field_int = integral_library.electric_field_integral(site=coord)
        #     V_ind += -1.0 * sum(ind_mom[i] * field_int[i] for i in range(3))
        E_pe = self.cppe_state.total_energy
        V_pe = self.V_es # + V_ind
        # only take electronic contributions into account
        # if elec_only:
        #     V_pe = V_ind
        #     E_pe = self.cppe_state.energies["Polarization"]["Electronic"]
        return E_pe, core.Matrix.from_array(V_pe)

    def build_electrostatics_operator(self):
        n_bas = self.basisset.nbf()
        self.V_es = np.zeros((n_bas, n_bas))
        for site in self.cppe_state.get_potentials():
            prefactors = []
            for multipole in site.multipoles:
                prefactors.extend(
                    cppe.prefactors(multipole.k) * multipole.values
                )
            integrals = self.mints.ao_multipole_potential(
                site.position, max_k=multipole.k
            )
            self.V_es += sum(
                pref * intv.np for pref, intv in zip(prefactors, integrals)
            )
