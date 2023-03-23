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
import numpy as np

from qcelemental import constants
from pkg_resources import parse_version

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

    if core.has_option_changed("DDX", "DDX_RADII"):
        radii = convfac * np.fromiter(core.get_option("DDX", "DDX_RADII"), float)
        if not radii.shape[0] == molecule.natom():
            raise ValidationError("Length of vector of custom cavity radii does "
                                  "not agree with number of atoms.")
    else:
        radii_set = core.get_option("DDX", "DDX_RADII_SET").lower()
        radii_lookup = getattr(pyddx.data, "radius_" + radii_set)

        # Use user-specified or ddx-recommended radius offset/scaling
        if core.has_option_changed("DDX", "DDX_RADII_SCALING"):
            radii_scaling = core.get_option("DDX", "DDX_RADII_SCALING")
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

    solvent = core.get_option("DDX", "DDX_SOLVENT").lower()
    if core.has_option_changed("DDX", "DDX_SOLVENT_EPSILON"):
        solvent_epsilon = core.get_option("DDX", "DDX_SOLVENT_EPSILON")
    elif solvent == "":
        raise ValidationError("Required option 'DDX_SOLVENT' is missing.")
    elif solvent not in pyddx.data.solvent_epsilon:
        raise ValidationError("Unknown solvent {solvent}.")
    else:
        solvent_epsilon = pyddx.data.solvent_epsilon[solvent]

    if core.has_option_changed("DDX", "DDX_SOLVENT_EPSILON_OPTICAL"):
        solvent_epsilon_optical = core.get_option("DDX", "DDX_SOLVENT_EPSILON_OPTICAL")
    else:
        solvent_epsilon_optical = pyddx.data.solvent_epsilon_optical.get(solvent, None)

    # TODO Deal with unit conversion (A^{-1} -> bohr^{-1})
    solvent_kappa = 0.0
    if core.get_option("DDX", "DDX_MODEL").lower() == "lpb":
        solvent_kappa = core.get_option("DDX", "DDX_SOLVENT_KAPPA")
        if solvent_kappa <= 0:
            raise ValidationError("DDX_SOLVENT_KAPPA is required for LPB and should be a "
                                  "positive quantity")

    fmm_multipole_lmax = core.get_option("DDX", "DDX_LMAX")
    if core.has_option_changed("DDX", "DDX_FMM_MULTIPOLE_LMAX"):
        fmm_multipole_lmax = core.get_option("DDX", "DDX_FMM_MULTIPOLE_LMAX")

    model_options = {
        "model": core.get_option("DDX", "DDX_MODEL").lower(),
        "solvent_kappa": solvent_kappa,
        "sphere_centres": molecule.geometry().np.T,
        "sphere_radii": radii,
        "lmax": core.get_option("DDX", "DDX_LMAX"),
        "n_lebedev": core.get_option("DDX", "DDX_N_LEBEDEV"),
        "maxiter": core.get_option("DDX", "DDX_MAXITER"),
        "jacobi_n_diis": core.get_option("DDX", "DDX_DIIS_MAX_VECS"),
        "n_proc": core.get_num_threads(),
        "incore": core.get_option("DDX", "DDX_INCORE"),
        "enable_fmm": core.get_option("DDX", "DDX_FMM"),
        "fmm_local_lmax": core.get_option("DDX", "DDX_FMM_LOCAL_LMAX"),
        "fmm_multipole_lmax": fmm_multipole_lmax,
        "logfile": core.get_option("DDX", "DDX_LOGFILE"),
        "eta": core.get_option("DDX", "DDX_ETA"),
        "shift": core.get_option("DDX", "DDX_SHIFT"),
    }
    dielectric_options = {
        "solvent_epsilon": solvent_epsilon,
        "solvent_epsilon_optical": solvent_epsilon_optical,
    }
    solver_options = {
        "tol": core.get_option("DDX", "DDX_SOLVATION_CONVERGENCE"),
    }
    grid_options = {
        # This overrides some DFT grid parameters just for the integrals
        # needed for DDX. The defaults are safe to avoid people from falling
        # into traps if they change their DFT grid setup.
        "DFT_SPHERICAL_POINTS": core.get_option("DDX", "DDX_SOLUTE_SPHERICAL_POINTS"),
        "DFT_RADIAL_POINTS": core.get_option("DDX", "DDX_SOLUTE_RADIAL_POINTS"),
        "DFT_NUCLEAR_SCHEME": "BECKE",  # Treutler and others might work here,
        "DFT_RADIAL_SCHEME": "BECKE",   # but this is so far untested with Psi4
        "DFT_PRUNING_SCHEME": "ROBUST",
        "DFT_BLOCK_SCHEME": "ATOMIC",
    }
    return {"model": model_options, "dielectric": dielectric_options,
            "solver": solver_options, "grid": grid_options}


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
        # verify that the minimal version is used if pyddx is provided
        # from outside the Psi4 ecosystem
        min_version = "0.2.0"
        if parse_version(pyddx.__version__) < parse_version(min_version):
            raise ModuleNotFoundError("pyddx version {} is required at least. "
                                      "Version {}"
                                      " was found.".format(min_version,
                                                           pyddx.__version__))

        self.basisset = basisset
        self.mints = core.MintsHelper(self.basisset)
        self.op_solver = options["solver"]
        op_grid = options["grid"]
        op_dielectric = options["dielectric"]

        # Setup the model
        try:
            self.model = pyddx.Model(**options["model"],
                                     solvent_epsilon=op_dielectric["solvent_epsilon"])
            e_optical = op_dielectric["solvent_epsilon_optical"]
            if e_optical:
                self.model_optical = pyddx.Model(**options["model"], solvent_epsilon=e_optical)
            else:
                self.model_optical = None
        except ValidationError as e:
            raise ValidationError(str(e))
        self.cavity = core.Matrix.from_array(self.model.cavity.T)
        self.sphere_charges = np.array([molecule.Z(i) for i in range(molecule.natom())])

        # Print summary of options
        core.print_out(pyddx.banner())
        core.print_out("\n")
        for k in sorted(list(self.model.input_parameters)):
            if k not in ("sphere_centres", "sphere_radii", "solvent_epsilon"):
                core.print_out(f"    {k:<23s} = {self.model.input_parameters[k]}\n")
        for k in sorted(list(op_dielectric.keys())):
            core.print_out(f"    {k:<23s} = {op_dielectric[k]}\n")
        for k in sorted(list(self.op_solver.keys())):
            core.print_out(f"    {k:<23s} = {self.op_solver[k]}\n")
        core.print_out(f"\n    DDX numerical integration setup:\n\n")
        for k in sorted(list(options["grid"].keys())):
            core.print_out(f"    {k.lower():<20s} = {op_grid[k]}\n")

        _print_cavity(self.sphere_charges, self.model.sphere_centres,
                      self.model.sphere_radii, "Angstrom")
        _print_cavity(self.sphere_charges, self.model.sphere_centres,
                      self.model.sphere_radii, "Bohr")
        core.print_out("\n")

        grid_int_opts = {k: v for (k, v) in op_grid.items() if isinstance(v, int)}
        grid_string_opts = {k: v for (k, v) in op_grid.items() if isinstance(v, str)}
        self.dftgrid = core.DFTGrid.build(molecule, basisset, grid_int_opts, grid_string_opts)
        self.numints = core.NumIntHelper(self.dftgrid)

        # Build the scaled ylms
        self.scaled_ylms = []
        for block in self.dftgrid.blocks():
            mtx = np.zeros((block.npoints(), self.model.n_basis))
            for (i, (x, y, z)) in enumerate(zip(block.x().np, block.y().np, block.z().np)):
                # Notice mtx[i, :] is contiguous in (row-major) memory
                self.model.scaled_ylm([x, y, z], block.parent_atom(), out=mtx[i, :])
            self.scaled_ylms.append(core.Matrix.from_array(mtx.T))

        # Build the nuclear contributions
        solute_multipoles = self.sphere_charges.reshape(1, -1) / np.sqrt(4 * np.pi)
        self.nuclear = self.model.multipole_electrostatics(solute_multipoles)
        self.nuclear["psi"] = self.model.multipole_psi(solute_multipoles)

    def get_solvation_contributions(self, density_matrix, state=None,
                                    elec_only=False, nonequilibrium=False):
        # TODO elec_only=True has not yet been tested. Will be properly integrated in a follow-up
        # TODO nonequilibrium=True has not yet been tested. Will be properly integrated in a follow-up
        #
        # Compute electronic contributions
        psi = self.numints.dd_density_integral(self.scaled_ylms, density_matrix).np.T
        dummy_charges = core.Vector.from_array(np.ones(self.model.n_cav))
        coords = core.Matrix.from_array(self.model.cavity.T)
        phi = self.mints.electrostatic_potential_value(dummy_charges, coords, density_matrix).np

        elec_field = None
        derivative_order = self.model.required_phi_derivative_order(compute_forces=False)
        if derivative_order > 0:
            elec_field = self.mints.electric_field_value(coords, density_matrix).np.T
            # print()
            # print()

            # delta = np.random.rand(3)
            # eps = 1e-4

            # coords_plus = core.Matrix.from_array(self.model.cavity.T + eps * delta)
            # coords_minus = core.Matrix.from_array(self.model.cavity.T - eps * delta)
            # phi_minus = self.mints.electrostatic_potential_value(dummy_charges, coords_minus, density_matrix).np
            # phi_plus  = self.mints.electrostatic_potential_value(dummy_charges, coords_plus, density_matrix).np
            # ref = (phi_plus - phi_minus) / 2 / eps

            # print((delta @ elec_field)[0:5])
            # print((ref)[0:5])

            # print(np.max(np.abs(delta @ elec_field - ref)))
            # print()
            # print()

        # elec_only = True
        if not elec_only:
            psi += self.nuclear["psi"]
            phi += self.nuclear["phi"]
            if elec_field is not None:
                elec_field += self.nuclear["e"]

        # psi = self.nuclear["psi"]
        # phi = self.nuclear["phi"]
        # elec_field = self.nuclear["e"]

        if elec_field is not None:
        #     elec_field = -elec_field
            print(elec_field[:, 0:10].T)

        if not nonequilibrium:
            model = self.model  # Use standard dielectric constant
        elif self.model_optical is None:  # Non-equilibrium not available
            raise ValueError("Non-equilibrium solvation not available, likely because no optical dielectric "
                             "constant is tabulated for the chosen solvent. Pleise provide the optical "
                             "dielectric constant manually using the option DDX_SOLVENT_EPSILON_OPTICAL "
                             "or change to a supported solvent (one of " +
                             ", ".join(pyddx.data.solvent_epsilon_optical) + ").")
        else:
            model = self.model_optical  # Use optical dielectric constant
        print("    Using epsilon: ", model.solvent_epsilon)

        # Update solvation problem with current phi, elec_field and psi
        if state is None:
            state = pyddx.State(model, psi, phi, elec_field)
            state.fill_guess(**self.op_solver)
            state.fill_guess_adjoint(**self.op_solver)
        else:
            state.update_problem(psi, phi, elec_field)

        # Solve problem and adjoint problem
        state.solve(**self.op_solver)
        state.solve_adjoint(**self.op_solver)

        print("    DDX Iters: ", state.x_n_iter, " ", state.s_n_iter)

        # Compute solvation energy
        fepsilon = 1.0
        if model.model == "cosmo":
            epsilon = model.solvent_epsilon
            fepsilon = (epsilon - 1) / epsilon
        E_ddx = fepsilon * state.energy()

        # Fock-matrix contributions
        eta = [core.Vector.from_array(ylm.np.T @ state.x[:, block.parent_atom()])
               for (block, ylm) in zip(self.dftgrid.blocks(), self.scaled_ylms)]
        V_ddx = self.numints.potential_integral(eta)

        extern = core.ExternalPotential()
        for xi_nj, pos_nj in zip(state.xi, model.cavity.T):
            extern.addCharge(xi_nj, pos_nj[0], pos_nj[1], pos_nj[2])
        V_ddx.add(extern.computePotentialMatrix(self.basisset))

        V_ddx.scale(-0.5 * fepsilon)  # Scale total potential
        return E_ddx, V_ddx, state
