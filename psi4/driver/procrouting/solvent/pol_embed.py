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

        def callback(output):
            print(output)
        self.cppe_state = cppe.CppeState(
            self.options, psi4mol_to_cppemol(self.molecule)
        )
        self.cppe_state.calculate_static_energies_and_fields()
        # obtain coordinates of polarizable sites
        self._enable_induction = False
        if self.cppe_state.get_polarizable_site_number():
            self._enable_induction = True
            coords = np.array([
                site.position for site in self.cppe_state.get_potentials()
                if site.is_polarizable
            ])
            self.polarizable_coords = core.Matrix.from_array(coords)
        self.V_es = None

    def get_pe_contribution(self, density_matrix, elec_only=False):
        # build electrostatics operator
        if self.V_es is None and not elec_only:
            self.build_electrostatics_operator()
        e_electrostatic = np.sum(density_matrix.np * self.V_es)
        self.cppe_state.energies["Electrostatic"]["Electronic"] = e_electrostatic

        n_bas = self.basisset.nbf()
        V_pe = np.zeros((n_bas, n_bas))
        if self._enable_induction:
            # obtain expectation values of elec. field at polarizable sites
            elec_fields = self.mints.electric_field_value(
                self.polarizable_coords, density_matrix
            ).np
            # solve induced moments
            self.cppe_state.update_induced_moments(
                elec_fields.flatten(), elec_only
            )
            induced_moments = np.array(
                self.cppe_state.get_induced_moments()
            ).reshape(self.polarizable_coords.shape)

            # build induction operator
            V_ind = self.mints.induction_operator(
                self.polarizable_coords,
                core.Matrix.from_array(induced_moments)
            ).np
            V_pe += V_ind
        # only take electronic contributions into account
        if elec_only:
            E_pe = self.cppe_state.energies["Polarization"]["Electronic"]
        else:
            V_pe += self.V_es
            E_pe = self.cppe_state.total_energy
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
