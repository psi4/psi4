#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2024 The Psi4 Developers.
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
        """ Initialize an MDIEngine object for communication with MDI (MolSSI driver interface)

        Parameters
        ----------
        scf_method
            Method (SCF or post-SCF) used when calculating energies or gradients.
        molecule
            The target molecule, if not the last molecule defined.
        kwargs
            Any additional arguments to pass to :func:`psi4.driver.energy` or
            :func:`psi4.driver.gradient` computation.

        """

        # Method used when the SCF command is received
        self.scf_method = scf_method

        # Additional arguments for energy, gradient, or optimization calculations
        self.kwargs = kwargs

        # Molecule all MDI operations are performed on
        input_molecule = kwargs.pop('molecule', psi4.core.get_active_molecule())
        _ini_cart = getattr(input_molecule, "_initial_cartesian", None)
        self.molecule = input_molecule.clone()
        if _ini_cart:
            self.molecule._initial_cartesian = _ini_cart
        psi4.core.set_active_molecule(self.molecule)

        # Most recent SCF energy
        self.energy = 0.0

        # Variables used when MDI sets a lattice of point charges
        self.nlattice = 0  # number of lattice point charges
        self.clattice = []  # list of lattice coordinates
        self.lattice = []  # list of lattice charges
        self.lattice_field = psi4.core.ExternalPotential()  # Psi4 chargefield

        # MPI variables
        self.mpi_world = None
        self.world_rank = 0

        # Flag for if a lattice of point charges has been set
        self.set_lattice = False

        # Get correct intra-code MPI communicator
        if use_mpi4py:
            self.mpi_world = MDI_MPI_get_world_comm()
            self.world_rank = self.mpi_world.Get_rank()

            # Psi4 does not currently support multiple MPI ranks
            if self.mpi_world.Get_size() != 1:
                MPI.COMM_WORLD.Abort()

        # Accept a communicator to the driver code
        self.comm = MDI_Accept_Communicator()

        # Ensure that the molecule is using c1 symmetry
        self.molecule.reset_point_group('c1')
        self.molecule.fix_orientation(True)
        self.molecule.fix_com(True)
        self.molecule.reinterpret_coordentry(False)
        self.molecule.update_geometry()

        # Flag to stop listening for MDI commands
        self.stop_listening = False

        # Dictionary of all supported MDI commands
        self.commands = {
            "<NATOMS": self.send_natoms,
            "<COORDS": self.send_coords,
            "<CHARGES": self.send_charges,
            "<ELEMENTS": self.send_elements,
            "<MASSES": self.send_masses,
            "<ENERGY": self.send_energy,
            "<FORCES": self.send_forces,
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

        # Register all the supported commands
        MDI_Register_Node("@DEFAULT")
        for command in self.commands.keys():
            MDI_Register_Command("@DEFAULT", command)

    def length_conversion(self):
        """ Obtain the conversion factor between the geometry specification units and bohr

        :returns: *unit_conv* Conversion factor between the geometry specification units and bohr
        """
        unit_name = self.molecule.units()
        if unit_name == "Angstrom":
            # beware if MDI and psi4 choose different sets of constants
            unit_conv = psi4.driver.constants.bohr2angstroms
        elif unit_name == "Bohr":
            unit_conv = 1.0
        else:
            raise Exception('Unrecognized unit type: ' + str(unit_name))

        return unit_conv

    # Respond to the <NATOMS command
    def send_natoms(self):
        """ Send the number of atoms through MDI
        """
        natom = self.molecule.natom()
        MDI_Send(natom, 1, MDI_INT, self.comm)
        return natom

    # Respond to the <COORDS command
    def send_coords(self):
        """ Send the nuclear coordinates through MDI
        """
        coords = self.molecule.geometry().np.ravel()
        MDI_Send(coords, len(coords), MDI_DOUBLE, self.comm)
        return coords

    # Respond to the <CHARGES command
    def send_charges(self):
        """ Send the nuclear charges through MDI

        :returns: *charges* Atomic charges
        """
        natom = self.molecule.natom()
        charges = [self.molecule.charge(iatom) for iatom in range(natom)]
        MDI_Send(charges, natom, MDI_DOUBLE, self.comm)
        return charges

    # Respond to the <MASSES command
    def send_masses(self):
        """ Send the nuclear masses through MDI

        :returns: *masses* Atomic masses
        """
        natom = self.molecule.natom()
        molecule_dict = self.molecule.to_dict()
        masses = molecule_dict['mass']
        MDI_Send(masses, natom, MDI_DOUBLE, self.comm)
        return masses

    # Respond to the <ELEMENTS command
    def send_elements(self):
        """ Send the atomic number of each nucleus through MDI

        :returns: *elements* Element of each atom
        """
        natom = self.molecule.natom()
        elements = [self.molecule.true_atomic_number(iatom) for iatom in range(natom)]
        MDI_Send(elements, natom, MDI_INT, self.comm)
        return elements

    # Respond to the <ENERGY command
    def send_energy(self):
        """ Send the total energy through MDI

        :returns: *energy* Energy of the system
        """
        self.run_scf()
        MDI_Send(self.energy, 1, MDI_DOUBLE, self.comm)
        return self.energy

    # Respond to the <FORCES command
    def send_forces(self):
        """ Send the nuclear forces through MDI

        :returns: *forces* Atomic forces
        """
        force_matrix = psi4.driver.gradient(self.scf_method, **self.kwargs)
        forces = force_matrix.np.ravel()
        MDI_Send(forces, len(forces), MDI_DOUBLE, self.comm)
        return forces

    # Respond to the >CHARGES command
    def recv_charges(self, charges=None):
        """ Receive a set of nuclear charges through MDI and assign them to the atoms in the current molecule

        Arguments:
            charges: New nuclear charges. If None, receive through MDI.
        """
        natom = self.molecule.natom()
        if charges is None:
            charges = MDI_Recv(natom, MDI_DOUBLE, self.comm)

        # Assign the charge of all atoms, taking care to avoid ghost atoms
        jatom = 0
        for iatom in range(natom):
            while self.molecule.fZ(jatom) == 0 and jatom < self.molecule.nallatom():
                jatom = jatom + 1
                if jatom >= self.molecule.nallatom():
                    raise Exception('Unexpected number of ghost atoms when receiving masses')
            self.molecule.set_nuclear_charge(iatom, charges[jatom])
            jatom = jatom + 1


    # Respond to the >COORDS command
    def recv_coords(self, coords=None):
        """ Receive a set of nuclear coordinates through MDI and assign them to the atoms in the current molecule

        Arguments:
            coords: New nuclear coordinates. If None, receive through MDI.
        """
        natom = self.molecule.natom()
        if coords is None:
            coords = MDI_Recv(3 * natom, MDI_DOUBLE, self.comm)
        matrix = psi4.core.Matrix.from_array(np.array(coords).reshape(-1, 3))
        self.molecule.set_geometry(matrix)
        self.molecule._initial_cartesian = matrix

    # Respond to the >MASSES command
    def recv_masses(self, masses=None):
        """ Receive a set of nuclear masses through MDI and assign them to the atoms in the current molecule

        Arguments:
            masses: New nuclear masses. If None, receive through MDI.
        """
        natom = self.molecule.natom()
        if masses is None:
            masses = MDI_Recv(natom, MDI_DOUBLE, self.comm)

        # Assign the mass of all atoms, taking care to avoid ghost atoms
        jatom = 0
        for iatom in range(natom):
            while self.molecule.fZ(jatom) == 0 and jatom < self.molecule.nallatom():
                jatom = jatom + 1
                if jatom >= self.molecule.nallatom():
                    raise Exception('Unexpected number of ghost atoms when receiving masses')
            self.molecule.set_mass(iatom, masses[jatom])
            jatom = jatom + 1

    # Set a lattice of point charges
    def set_lattice_field(self):
        """ Set a field of lattice point charges using information received through MDI
        """
        arr = []
        for ilat in range(self.nlattice):
            arr.append(self.lattice[ilat])
            arr.append(self.clattice[3 * ilat + 0])
            arr.append(self.clattice[3 * ilat + 1])
            arr.append(self.clattice[3 * ilat + 2])
        self.kwargs["external_potentials"] = np.array(arr).reshape((-1, 4))
        self.set_lattice = True

    # Respond to the >NLATTICE command
    def recv_nlattice(self, nlattice=None):
        """ Receive the number of lattice point charges through MDI

        Arguments:
            nlattice: New number of point charges. If None, receive through MDI.
        """
        if nlattice is None:
            self.nlattice = MDI_Recv(1, MDI_INT, self.comm)
        else:
            self.nlattice = nlattice
        self.clattice = [0.0 for ilat in range(3 * self.nlattice)]
        self.lattice = [0.0 for ilat in range(self.nlattice)]
        self.set_lattice_field()

    # Respond to the >CLATTICE command
    def recv_clattice(self, clattice=None):
        """ Receive the coordinates of a set of lattice point charges through MDI

        Arguments:
            clattice: New coordinates of the lattice of point charges. If None, receive through MDI.
        """
        if clattice is None:
            self.clattice = MDI_Recv(3 * self.nlattice, MDI_DOUBLE, self.comm)
        else:
            self.clattice = clattice
        self.set_lattice_field()

    # Respond to the >LATTICE command
    def recv_lattice(self, lattice=None):
        """ Receive the charges of a set of lattice point charges through MDI

        Arguments:
            lattice: New charges of the lattice of point charges. If None, receive through MDI.
        """
        if lattice is None:
            self.lattice = MDI_Recv(self.nlattice, MDI_DOUBLE, self.comm)
        else:
            self.lattice = lattice
        self.set_lattice_field()

    # Respond to the SCF command
    def run_scf(self):
        """ Run an energy calculation
        """
        self.energy = psi4.energy(self.scf_method, **self.kwargs)

    # Respond to the <DIMENSIONS command
    def send_dimensions(self):
        """ Send the dimensionality of the system through MDI

        :returns: *dimensions* Dimensionality of the system
        """
        dimensions = [1, 1, 1]
        MDI_Send(dimensions, 3, MDI_INT, self.comm)
        return dimensions

    # Respond to the <TOTCHARGE command
    def send_total_charge(self):
        """ Send the total system charge through MDI

        :returns: *charge* Total charge of the system
        """
        charge = self.molecule.molecular_charge()
        MDI_Send(charge, 1, MDI_DOUBLE, self.comm)
        return charge

    # Respond to the >TOTCHARGE command
    def recv_total_charge(self, charge=None):
        """ Receive the total system charge through MDI

        Arguments:
            charge: New charge of the system. If None, receive through MDI.
        """
        if charge is None:
            charge = MDI_Recv(1, MDI_DOUBLE, self.comm)
        self.molecule.set_molecular_charge(int(round(charge)))

    # Respond to the <ELEC_MULT command
    def send_multiplicity(self):
        """ Send the electronic multiplicity through MDI

        :returns: *multiplicity* Multiplicity of the system
        """
        multiplicity = self.molecule.multiplicity()
        MDI_Send(multiplicity, 1, MDI_INT, self.comm)
        return multiplicity

    # Respond to the >ELEC_MULT command
    def recv_multiplicity(self, multiplicity=None):
        """ Receive the electronic multiplicity through MDI

        Arguments:
            multiplicity: New multiplicity of the system. If None, receive through MDI.
        """
        if multiplicity is None:
            multiplicity = MDI_Recv(1, MDI_INT, self.comm)
        self.molecule.set_multiplicity(multiplicity)

    # Respond to the EXIT command
    def exit(self):
        """ Stop listening for MDI commands
        """
        self.stop_listening = True

        # If a lattice of point charges was set, unset it now
        if self.set_lattice:
            self.kwargs.pop("external_potentials", None)
            

    # Enter server mode, listening for commands from the driver
    def listen_for_commands(self):
        """ Receive commands through MDI and respond to them as defined by the MDI Standard
        """

        while not self.stop_listening:
            if self.world_rank == 0:
                command = MDI_Recv_Command(self.comm)
            else:
                command = None
            if use_mpi4py:
                command = self.mpi_world.bcast(command, root=0)
            if self.world_rank == 0:
                psi4.core.print_out('\nMDI command received: ' + str(command) + ' \n')

            # Search for this command in self.commands
            found_command = False
            for supported_command in self.commands:
                if not found_command and command == supported_command:
                    # Run the function corresponding to this command
                    self.commands[supported_command]()
                    found_command = True
            if not found_command:
                raise Exception('Unrecognized command: ' + str(command))


def mdi_init(mdi_arguments):
    """ Initialize the MDI Library

    Parameters
    ----------
    mdi_arguments
        MDI configuration options

    """
    MDI_Init(mdi_arguments)


def mdi_run(scf_method: str, **kwargs):
    """ Begin functioning as an MDI (MolSSI driver interface) engine

    Parameters
    ----------
    scf_method
        Method (SCF or post-SCF) used when calculating energies or gradients.
    molecule
        The target molecule, if not the last molecule defined.
    kwargs
        Any additional arguments to pass to :func:`psi4.driver.energy` or
        :func:`psi4.driver.gradient` computation.

    """
    engine = MDIEngine(scf_method, **kwargs)
    engine.listen_for_commands()
