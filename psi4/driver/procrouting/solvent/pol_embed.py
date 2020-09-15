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
from qcelemental import constants
from pkg_resources import parse_version

from psi4 import core
from psi4.driver.p4util.exceptions import ValidationError


def get_pe_options():
    if core.get_option('SCF', 'PCM'):
        raise ValidationError("""Error: 3-layer QM/PE/PCM not implemented.\n""")
    rmin = core.get_option('PE', 'BORDER_RMIN')
    if core.get_option('PE', 'BORDER_RMIN_UNIT').upper() == "AA":
        rmin *= 1.0 / constants.bohr2angstroms
    pol_embed_options = {
        "potfile": core.get_option('PE', 'POTFILE'),
        "iso_pol": core.get_option('PE', 'ISOTROPIC_POL'),
        "induced_thresh": core.get_option('PE', 'INDUCED_CONVERGENCE'),
        "maxiter": core.get_option('PE', 'MAXITER'),
        # tree options
        "summation_induced_fields": core.get_option('PE', 'SUMMATION_FIELDS').lower(),
        "tree_expansion_order": core.get_option('PE', 'TREE_EXPANSION_ORDER'),
        "theta": core.get_option('PE', 'TREE_THETA'),
        # damping options
        "damp_induced": core.get_option('PE', 'DAMP_INDUCED'),
        "damping_factor_induced": core.get_option('PE', 'DAMPING_FACTOR_INDUCED'),
        "damp_multipole": core.get_option('PE', 'DAMP_MULTIPOLE'),
        "damping_factor_multipole": core.get_option('PE', 'DAMPING_FACTOR_MULTIPOLE'),
        "pe_border": core.get_option('PE', 'BORDER'),
        "border_type": core.get_option('PE', 'BORDER_TYPE').lower(),
        "border_rmin": rmin,
        "border_nredist": core.get_option('PE', 'BORDER_N_REDIST'),
        "border_redist_order": core.get_option('PE', 'BORDER_REDIST_ORDER'),
        "border_redist_pol": core.get_option('PE', 'BORDER_REDIST_POL'),
    }
    return pol_embed_options


def psi4mol_to_cppemol(psi4mol):
    mol = cppe.Molecule()
    geom = psi4mol.geometry().np
    for i, c in enumerate(geom):
        a = cppe.Atom((int(psi4mol.Z(i))), *c)
        mol.append(a)
    return mol


class CppeInterface:
    def __init__(self, molecule, options, basisset):
        # verify that the minimal version is used if CPPE is provided
        # from outside the Psi4 ecosystem
        min_version = "0.2.0"
        if parse_version(cppe.__version__) < parse_version(min_version):
            raise ModuleNotFoundError("CPPE version {} is required at least. "
                                      "Version {}"
                                      " was found.".format(min_version,
                                                           cppe.__version__))
        # setup the initial CppeState
        self.molecule = molecule
        self.options = options
        self.basisset = basisset
        self.mints = core.MintsHelper(self.basisset)

        def callback(output):
            core.print_out("{}\n".format(output))

        self.cppe_state = cppe.CppeState(self.options, psi4mol_to_cppemol(self.molecule), callback)
        core.print_out("CPPE Options:\n")
        for k in cppe.valid_option_keys:
            core.print_out(f"{k} = {self.cppe_state.options[k]}\n")
        core.print_out("-------------------------\n\n")
        self.cppe_state.calculate_static_energies_and_fields()
        # obtain coordinates of polarizable sites
        self._enable_induction = False
        if self.cppe_state.get_polarizable_site_number():
            self._enable_induction = True
            coords = self.cppe_state.positions_polarizable
            self.polarizable_coords = core.Matrix.from_array(coords)
        self.V_es = None

    def get_pe_contribution(self, density_matrix, elec_only=False):
        # build electrostatics operator
        if self.V_es is None and not elec_only:
            self.build_electrostatics_operator()

        n_bas = self.basisset.nbf()
        V_pe = np.zeros((n_bas, n_bas))
        if self._enable_induction:
            # obtain expectation values of elec. field at polarizable sites
            elec_fields = self.mints.electric_field_value(self.polarizable_coords, density_matrix).np
            # solve induced moments
            self.cppe_state.update_induced_moments(elec_fields.flatten(), elec_only)
            induced_moments = np.array(self.cppe_state.get_induced_moments()).reshape(self.polarizable_coords.shape)

            # build induction operator
            V_ind = self.mints.induction_operator(self.polarizable_coords, core.Matrix.from_array(induced_moments)).np
            V_pe += V_ind
        # only take electronic contributions into account
        if elec_only:
            E_pe = self.cppe_state.energies["Polarization"]["Electronic"]
        else:
            e_el = np.sum(density_matrix.np * self.V_es)
            self.cppe_state.energies["Electrostatic"]["Electronic"] = e_el
            V_pe += self.V_es
            E_pe = self.cppe_state.total_energy
        return E_pe, core.Matrix.from_array(V_pe)

    def build_electrostatics_operator(self):
        n_bas = self.basisset.nbf()
        self.V_es = np.zeros((n_bas, n_bas))
        for site in self.cppe_state.potentials:
            prefactors = []
            for multipole in site.multipoles:
                prefactors.extend(cppe.prefactors(multipole.k) * multipole.values)
            integrals = self.mints.ao_multipole_potential(site.position, max_k=multipole.k)
            self.V_es += sum(pref * intv.np for pref, intv in zip(prefactors, integrals))
