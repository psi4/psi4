#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2023 The Psi4 Developers.
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

from tempfile import NamedTemporaryFile

import numpy as np
import cppe
from qcelemental.util import parse_version

from ...constants import constants
from psi4 import core
from ...qcdb import libmintsbasisset
from ...p4util.exceptions import ValidationError


def get_pe_options():
    if core.get_option('SCF', 'PCM') or core.get_option('SCF', 'DDX'):
        raise ValidationError("""Error: 3-layer QM/PE/continuum solvation not implemented.\n""")
    rmin = core.get_option('PE', 'BORDER_RMIN')
    if core.get_option('PE', 'BORDER_RMIN_UNIT').upper() == "AA":
        rmin *= 1.0 / constants.bohr2angstroms

    # potfile option can be filename or contents
    potfile_keyword = core.get_option('PE', 'POTFILE')
    if "@COORDINATES" in potfile_keyword:
        fl = NamedTemporaryFile(mode="w+t", delete=False)
        fl.write(potfile_keyword)
        fl.close()
        potfile_name = fl.name
    else:
        potfile_name = potfile_keyword

    pol_embed_options = {
        "potfile": potfile_name,
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
        # PE(ECP)
        "pe_ecp": core.get_option('PE', 'PE_ECP'),
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
        min_version = "0.3.1"
        if parse_version(cppe.__version__) < parse_version(min_version):
            raise ModuleNotFoundError("CPPE version {} is required at least. "
                                      "Version {}"
                                      " was found.".format(min_version,
                                                           cppe.__version__))
        # setup the initial CppeState
        self.molecule = molecule
        self.options = options
        self.pe_ecp = self.options.pop("pe_ecp", False)
        self.basisset = basisset
        self.mints = core.MintsHelper(self.basisset)

        def callback(output):
            core.print_out(f"{output}\n")

        self.cppe_state = cppe.CppeState(self.options, psi4mol_to_cppemol(self.molecule), callback)
        core.print_out("CPPE Options:\n")
        core.print_out(f"PE(ECP) repulsive potentials = {self.pe_ecp}\n")
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
        if self.pe_ecp:
            self._setup_pe_ecp()

    def _setup_pe_ecp(self):
        mol_clone = self.molecule.clone()
        geom, _, elems, _, _ = mol_clone.to_arrays()
        n_qmatoms = len(elems)
        geom = geom.tolist()
        elems = elems.tolist()
        for p in self.cppe_state.potentials:
            if p.element == "X":
                continue
            elems.append(f"{p.element}_pe")
            geom.append([p.x, p.y, p.z])
        qmmm_mol = core.Molecule.from_arrays(
            geom=geom, elbl=elems, units="Bohr",
            fix_com=True, fix_orientation=True, fix_symmetry="c1"
        )
        n_qmmmatoms = len(elems)

        def __basisspec_pe_ecp(mol, role):
            global_basis = core.get_global_option("BASIS")
            for i in range(n_qmatoms):
                mol.set_basis_by_number(i, global_basis, role=role)
            for i in range(n_qmatoms, n_qmmmatoms):
                mol.set_basis_by_number(i, "pe_ecp", role=role)
            return {}

        libmintsbasisset.basishorde["PE_ECP_BASIS"] = __basisspec_pe_ecp
        self.pe_ecp_basis = core.BasisSet.build(qmmm_mol, "BASIS", "PE_ECP_BASIS")
        ecp_mints = core.MintsHelper(self.pe_ecp_basis)
        self.V_pe_ecp = ecp_mints.ao_ecp().np


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
            E_ecp = 0.0
            if self.pe_ecp:
                E_ecp = np.sum(density_matrix.np * self.V_pe_ecp)
                V_pe += self.V_pe_ecp
            E_pe = self.cppe_state.total_energy + E_ecp
        return E_pe, core.Matrix.from_array(V_pe)

    def build_electrostatics_operator(self):
        n_bas = self.basisset.nbf()
        self.V_es = np.zeros((n_bas, n_bas))
        # TODO: run this in one go...
        for site in self.cppe_state.potentials:
            prefactors = []
            for multipole in site.multipoles:
                prefactors.extend(cppe.prefactors(multipole.k) * multipole.values)
            integrals = self.mints.ao_multipole_potential(order=multipole.k, origin=site.position)
            self.V_es += sum(pref * intv.np for pref, intv in zip(prefactors, integrals))
