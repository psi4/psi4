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

"""Support for using Psi4 as an MDI engine.
For details regarding MDI, see https://molssi.github.io/MDI_Library/html/index.html.
"""

__all__ = [
    "mdi_init",
    "mdi_run",
    "MDIEngine",
]

import numpy as np
import qcelemental as qcel
import psi4

_have_mdi = False
try:
    from mdi import (
        MDI_DOUBLE,
        MDI_INT,
        MDI_Accept_Communicator,
        MDI_Init,
        MDI_MPI_get_world_comm,
        MDI_Recv,
        MDI_Recv_Command,
        MDI_Register_Command,
        MDI_Register_Node,
        MDI_Send,
    )
    _have_mdi = True
except ImportError:
    pass

try:
    from mpi4py import MPI
    use_mpi4py = True
except ImportError:
    use_mpi4py = False

class MDIEngine():
    def __init__(self, scf_method: str, **kwargs):
        self.scf_method = scf_method
        self.kwargs = kwargs

        input_molecule = kwargs.pop('molecule', psi4.core.get_active_molecule())
        _ini_cart = getattr(input_molecule, "_initial_cartesian", None)
        self.molecule = input_molecule.clone()
        if _ini_cart:
            self.molecule._initial_cartesian = _ini_cart
        psi4.core.set_active_molecule(self.molecule)

        # Cache variables
        self.computed_this_step = False
        self.cached_energy = 0.0
        self.cached_qm_forces = None
        self.cached_lattice_forces = None

        self.nlattice = 0
        self.clattice = []
        self.lattice = []
        self.lattice_field = psi4.core.ExternalPotential()
        self.set_lattice = False

        self.mpi_world = None
        self.world_rank = 0
        if use_mpi4py:
            self.mpi_world = MDI_MPI_get_world_comm()
            self.world_rank = self.mpi_world.Get_rank()
            if self.mpi_world.Get_size() != 1:
                MPI.COMM_WORLD.Abort()

        self.comm = MDI_Accept_Communicator()
        self.molecule.reset_point_group('c1')
        self.molecule.fix_orientation(True)
        self.molecule.fix_com(True)
        self.molecule.reinterpret_coordentry(False)
        self.molecule.update_geometry()

        self.stop_listening = False

        self.commands = {
            "<NATOMS": self.send_natoms,
            "<COORDS": self.send_coords,
            "<CHARGES": self.send_charges,
            "<ELEMENTS": self.send_elements,
            "<MASSES": self.send_masses,
            "<ENERGY": self.send_energy,
            "<FORCES": self.send_forces,
            "<LATTICE_FORCES": self.send_lattice_forces,
            ">COORDS": self.recv_coords,
            ">NLATTICE": self.recv_nlattice,
            ">CLATTICE": self.recv_clattice,
            ">LATTICE": self.recv_lattice,
            ">MASSES": self.recv_masses,
            "SCF": self.run_scf,
            "<DIMENSIONS": self.send_dimensions,
            "<TOTCHARGE": self.send_total_charge,
            ">TOTCHARGE": self.recv_total_charge,
            "<ELEC_MULT": self.send_multiplicity,
            ">ELEC_MULT": self.recv_multiplicity,
            "EXIT": self.exit
        }

        MDI_Register_Node("@DEFAULT")
        for command in self.commands.keys():
            MDI_Register_Command("@DEFAULT", command)

    def reset_cache(self):
        """ Reset the computation cache whenever geometry or lattice changes """
        self.computed_this_step = False

    def compute_all_results(self):
        """ Perform the heavy lifting (SCF + Gradient) only once per MDI step """
        if not self.computed_this_step:
            # gradient() triggers energy() internally
            # return_wfn=True is required to access the ExternalPotential back-reaction
            grad_matrix, wfn = psi4.driver.gradient(self.scf_method, return_wfn=True, **self.kwargs)
            
            # QM Forces (Nuclear Gradients * -1.0)
            self.cached_qm_forces = grad_matrix.np.ravel() * -1.0
            self.cached_energy = wfn.energy()
            
            # Lattice Forces (Back-reaction on point charges)
            if self.set_lattice and wfn.get_external_potential():
                # ExternalPotential.get_v_gradient returns dE/dR_ext
                # Force is -dE/dR_ext
                lat_grad = wfn.get_external_potential().get_v_gradient()
                self.cached_lattice_forces = lat_grad.np.ravel() * -1.0
            else:
                self.cached_lattice_forces = np.zeros(3 * self.nlattice)
            
            self.computed_this_step = True

    def send_energy(self):
        self.compute_all_results()
        MDI_Send(self.cached_energy, 1, MDI_DOUBLE, self.comm)
        return self.cached_energy

    def send_forces(self):
        self.compute_all_results()
        MDI_Send(self.cached_qm_forces, len(self.cached_qm_forces), MDI_DOUBLE, self.comm)
        return self.cached_qm_forces

    def send_lattice_forces(self):
        self.compute_all_results()
        MDI_Send(self.cached_lattice_forces, len(self.cached_lattice_forces), MDI_DOUBLE, self.comm)
        return self.cached_lattice_forces

    def recv_coords(self, coords=None):
        natom = self.molecule.natom()
        if coords is None:
            coords = MDI_Recv(3 * natom, MDI_DOUBLE, self.comm)
        matrix = psi4.core.Matrix.from_array(np.array(coords).reshape(-1, 3))
        self.molecule.set_geometry(matrix)
        self.molecule._initial_cartesian = matrix
        self.reset_cache()

    def recv_nlattice(self, nlattice=None):
        if nlattice is None:
            self.nlattice = MDI_Recv(1, MDI_INT, self.comm)
        else:
            self.nlattice = nlattice
        self.clattice = [0.0 for ilat in range(3 * self.nlattice)]
        self.lattice = [0.0 for ilat in range(self.nlattice)]
        self.reset_cache()
        self.set_lattice_field()

    def recv_clattice(self, clattice=None):
        if clattice is None:
            self.clattice = MDI_Recv(3 * self.nlattice, MDI_DOUBLE, self.comm)
        else:
            self.clattice = clattice
        self.reset_cache()
        self.set_lattice_field()

    def recv_lattice(self, lattice=None):
        if lattice is None:
            self.lattice = MDI_Recv(self.nlattice, MDI_DOUBLE, self.comm)
        else:
            self.lattice = lattice
        self.reset_cache()
        self.set_lattice_field()

    def set_lattice_field(self):
        arr = []
        for ilat in range(self.nlattice):
            arr.append(self.lattice[ilat])
            arr.append(self.clattice[3 * ilat + 0])
            arr.append(self.clattice[3 * ilat + 1])
            arr.append(self.clattice[3 * ilat + 2])
        self.kwargs["external_potentials"] = np.array(arr).reshape((-1, 4))
        self.set_lattice = True

    # [Remaining helper methods: send_natoms, run_scf, etc. same as original]
    def send_natoms(self):
        natom = self.molecule.natom()
        MDI_Send(natom, 1, MDI_INT, self.comm)
        return natom

    def send_coords(self):
        coords = self.molecule.geometry().np.ravel()
        MDI_Send(coords, len(coords), MDI_DOUBLE, self.comm)
        return coords

    def send_charges(self):
        natom = self.molecule.natom()
        charges = [self.molecule.charge(iatom) for iatom in range(natom)]
        MDI_Send(charges, natom, MDI_DOUBLE, self.comm)
        return charges

    def send_masses(self):
        natom = self.molecule.natom()
        molecule_dict = self.molecule.to_dict()
        masses = molecule_dict['mass']
        MDI_Send(masses, natom, MDI_DOUBLE, self.comm)
        return masses

    def send_elements(self):
        natom = self.molecule.natom()
        elements = [self.molecule.true_atomic_number(iatom) for iatom in range(natom)]
        MDI_Send(elements, natom, MDI_INT, self.comm)
        return elements

    def run_scf(self):
        self.compute_all_results()

    def send_dimensions(self):
        dimensions = [1, 1, 1]
        MDI_Send(dimensions, 3, MDI_INT, self.comm)
        return dimensions

    def send_total_charge(self):
        charge = self.molecule.molecular_charge()
        MDI_Send(charge, 1, MDI_DOUBLE, self.comm)
        return charge

    def recv_total_charge(self, charge=None):
        if charge is None:
            charge = MDI_Recv(1, mdi.MDI_DOUBLE, self.comm)
        self.molecule.set_molecular_charge(int(round(charge)))

    def send_multiplicity(self):
        multiplicity = self.molecule.multiplicity()
        MDI_Send(multiplicity, 1, MDI_INT, self.comm)
        return multiplicity

    def recv_multiplicity(self, multiplicity=None):
        if multiplicity is None:
            multiplicity = MDI_Recv(1, MDI_INT, self.comm)
        self.molecule.set_multiplicity(multiplicity)

    def exit(self):
        self.stop_listening = True
        if self.set_lattice:
            self.kwargs.pop("external_potentials", None)

    def listen_for_commands(self):
        while not self.stop_listening:
            if self.world_rank == 0:
                command = MDI_Recv_Command(self.comm)
            else:
                command = None
            if use_mpi4py:
                command = self.mpi_world.bcast(command, root=0)
            if self.world_rank == 0:
                psi4.core.print_out('\nMDI command received: ' + str(command) + ' \n')

            found_command = False
            for supported_command in self.commands:
                if not found_command and command == supported_command:
                    self.commands[supported_command]()
                    found_command = True
            if not found_command:
                raise Exception('Unrecognized command: ' + str(command))

def mdi_init(mdi_arguments):
    MDI_Init(mdi_arguments)

def mdi_run(scf_method: str, **kwargs):
    engine = MDIEngine(scf_method, **kwargs)
    engine.listen_for_commands()
