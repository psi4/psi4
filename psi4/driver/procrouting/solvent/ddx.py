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
import numpy as np

from qcelemental import constants

from psi4 import core
from psi4.driver.p4util.exceptions import ValidationError

import pyddx
import pyddx.data


def get_ddx_options(molecule):
    if core.get_option("SCF", "PCM"):
        raise ValidationError("Error: DDX and PCM are exclusive.")

    if molecule.units() == "Angstrom":
        convfac = 1 / constants.bohr2angstroms
    else:
        convfac = 1.0

    if core.has_option_changed("DDX", "RADII"):
        radii = convfac * np.fromiter(core.get_option("DDX", "RADII"), float)
        if not radii.shape[0] == molecule.natom():
            raise ValidationError("Length of vector of custom cavity radii does "
                                  "not agree with number of atoms.")
    else:
        radii_set = core.get_option("DDX", "RADII_SET").lower()
        radii_lookup = getattr(pyddx.data, "radius_" + radii_set)

        # Use user-specified or ddx-recommended radius offset/scaling
        if core.has_option_changed("DDX", "RADII_SCALING"):
            radii_scaling = core.get_option("DDX", "RADII_SCALING")
        else:
            radii_scaling = getattr(pyddx.data, "radius_" + radii_set + "_scaling")

        radii = np.empty((molecule.natom()))
        for i in range(molecule.natom()):
            element = molecule.label(i).lower()
            if element not in radii_lookup:
                raise ValidationError(f"Element {element} not known in radii_set {radii_set}. "
                                      "Please choose a different radii_set or provide radii "
                                      "manually using the 'RADII' option.")
            radii[i] = radii_scaling * radii_lookup[element]

    solvent = core.get_option("DDX", "SOLVENT").lower()
    if core.has_option_changed("DDX", "SOLVENT_EPSILON"):
        solvent_epsilon = core.get_option("DDX", "SOLVENT_EPSILON")
    elif solvent == "":
        raise ValidationError("Required option 'SOLVENT' is missing.")
    elif solvent not in pyddx.data.solvent_epsilon:
        raise ValidationError("Unknown solvent {solvent}.")
    else:
        solvent_epsilon = pyddx.data.solvent_epsilon[solvent]

    model_options = {
        "model": core.get_option("DDX", "MODEL").lower(),
        "sphere_charges": np.array([molecule.Z(i) for i in range(molecule.natom())]),
        "sphere_centres": molecule.geometry().np.T,
        "sphere_radii": radii,
        "solvent_epsilon": solvent_epsilon,
        "eta": core.get_option("DDX", "ETA"),
        "lmax": core.get_option("DDX", "LMAX"),
        "n_lebedev": core.get_option("DDX", "N_LEBEDEV"),
        "maxiter": 100,       # TODO Configurable
        "jacobi_n_diis": 20,  # TODO Configurable
    }
    solver_options = {
        "tol": 1e-8,  # TODO Configurable!
    }
    return {"model": model_options, "solver": solver_options }


def _print_cavity(charges, centres, radii, unit="Angstrom"):
    if unit == "Angstrom":
        convfac = constants.bohr2angstroms
    elif unit == "Bohr":
        convfac = 1.0
    core.print_out(f"\n    Cavity sphere setup (in {unit}):\n\n")
    core.print_out(("    {0:^6s}   {1:^10s}   {2:^10s}   {3:^10s}   {4:^10s}"
                   "\n").format("Charge", "X", "Y", "Z", "Radius"))
    core.print_out("    " + "-" * 6 + ("   " + "-" * 10) * 4 + "\n")
    for i in range(len(charges)):
        core.print_out(f"    {charges[i]:6.2f}")
        for j in range(3):
            core.print_out(f"   {convfac * centres[j, i]:10.6f}")
        core.print_out(f"   {convfac * radii[i]:10.6f}\n")


class DdxInterface:
    def __init__(self, molecule, options, basisset):
        self.basisset = basisset
        self.mints = core.MintsHelper(self.basisset)
        self.op_solver = options["solver"]

        # Setup the model
        try:
            self.model = pyddx.Model(**options["model"])
        except ValidationError as e:
            raise ValidationError(str(e))
        self.cavity = core.Matrix.from_array(self.model.cavity.T)

        # Print summary of options
        core.print_out("  ==> DDX options <==\n")
        for k in sorted(list(self.model.input_parameters)):
            if k not in ("sphere_charges", "sphere_centres", "sphere_radii"):
                core.print_out(f"    {k:<15s} = {self.model.input_parameters[k]}\n")
        for k in sorted(list(self.op_solver.keys())):
            core.print_out(f"    {k:<15s} = {self.op_solver[k]}\n")
        _print_cavity(self.model.sphere_charges, self.model.sphere_centres,
                      self.model.sphere_radii, "Angstrom")
        _print_cavity(self.model.sphere_charges, self.model.sphere_centres,
                      self.model.sphere_radii, "Bohr")
        core.print_out("\n")

        # DFT grid for numerical integration
        # TODO Make configurable!
        int_opts = {
            # "DFT_BLOCK_MAX_POINTS": 0,
            # "DFT_BLOCK_MIN_POINTS": 0,
            # "DFT_SPHERICAL_POINTS": 0,
            # "DFT_RADIAL_POINTS": 0,
            # "PRINT": 0,
            # "DEBUG": 0,
            # "BENCH": 0,
        }
        string_opts = {
            # "DFT_RADIAL_SCHEME": "",
            # "DFT_PRUNING_SCHEME": "",
            # "DFT_NUCLEAR_SCHEME": "",
            # "DFT_GRID_NAME": "",
            "DFT_BLOCK_SCHEME": "ATOMIC",
        }
        self.dftgrid = core.DFTGrid.build(molecule, basisset, int_opts, string_opts)
        self.numints = core.NumIntHelper(self.dftgrid)

        # Build the scaled ylms
        self.scaled_ylms = []
        for block in self.dftgrid.blocks():
            mtx = np.zeros((block.npoints(), self.model.n_basis))
            for (i, (x, y, z)) in enumerate(zip(block.x().np, block.y().np, block.z().np)):
                # Notice mtx[i, :] is contiguous in (row-major) memory
                self.model.scaled_ylm([x, y, z], block.parent_atom(), out=mtx[i, :])
            self.scaled_ylms.append(core.Matrix.from_array(mtx.T))

        self.state = self.model.initial_guess()

    def get_solvation_contributions(self, density_matrix, elec_only=False):
        assert not elec_only  # Not yet implemented

        # Initialise with nuclear contributions:
        nuclear = self.model.solute_nuclear_contribution()
        phi, psi = nuclear["phi"], nuclear["psi"]

        # Add electronic contributions
        psi += self.numints.dd_density_integral(self.scaled_ylms, density_matrix).np.T
        dummy_charges = core.Vector.from_array(np.ones(self.model.n_cav))
        coords = core.Matrix.from_array(self.model.cavity.T)
        phi += self.mints.electrostatic_potential_value(dummy_charges, coords, density_matrix).np

        # Solve problem and adjoint problem
        self.state = self.model.solve(self.state, phi, **self.op_solver)
        self.state = self.model.adjoint_solve(self.state, psi, **self.op_solver)

        # Compute solvation energy
        f_epsilon = 1.0
        if self.model.model == "cosmo":
            epsilon = self.model.solvent_epsilon
            fepsilon = (epsilon - 1) / epsilon
        E_ddx = 0.5 * fepsilon * np.sum(self.state.x * psi)

        # Fock-matrix contributions
        eta = [core.Vector.from_array(ylm.np.T @ self.state.x[:, block.parent_atom()])
               for (block, ylm) in zip(self.dftgrid.blocks(), self.scaled_ylms)]
        V_ddx = self.numints.potential_integral(eta)

        # This hack is needed because ExternalPotential() assumes the same
        # units as the molecule.
        convfac = 1.0
        if self.basisset.molecule().units() == "Angstrom":
            convfac = constants.bohr2angstroms

        extern = core.ExternalPotential()
        for xi_nj, pos_nj in zip(self.state.xi, self.model.cavity.T):
            pos_nj = convfac * pos_nj
            extern.addCharge(xi_nj, pos_nj[0], pos_nj[1], pos_nj[2])
        V_ddx.add(extern.computePotentialMatrix(self.basisset))

        V_ddx.scale(-0.5 * fepsilon)  # Scale total potential
        return E_ddx, V_ddx
