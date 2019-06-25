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
"""Support for using Psi4 as an MDI engine.
For details regarding MDI, see https://molssi.github.io/MDI_Library/html/index.html.

"""
import psi4
from psi4.driver import p4util
from psi4.driver import driver
from psi4.driver.p4util.exceptions import *

from MDI_Library.mdi import MDI_Init, MDI_Get_Intra_Code_MPI_Comm, MDI_Accept_Communicator
from MDI_Library.mdi import MDI_Send, MDI_Recv, MDI_Recv_Command, MDI_Conversion_Factor
from MDI_Library.mdi import MDI_INT, MDI_DOUBLE, MDI_CHAR

try:
    from mpi4py import MPI
    use_mpi4py = True
except ImportError:
    use_mpi4py = False

class MDIEngine():

    def __init__(self, scf_method):
        """ Initialize an MDIEngine object for communication with MDI

        Arguments:
           scf_method: Method used when calculating energies or gradients
        """

        # Method used when the SCF command is received
        self.scf_method = scf_method

        self.molecule = psi4.core.get_active_molecule()
        self.energy = 0.0

        # lattice variables
        self.nlattice = 0 # number of lattice point charges
        self.clattice = [] # list of lattice coordinates
        self.lattice = [] # list of lattice charges
        self.lattice_field = psi4.QMMM() # Psi4 chargefield

        # MPI variables
        self.mpi_world = None
        self.world_rank = 0

        # Get correct intra-code MPI communicator
        if use_mpi4py:
            self.mpi_world = MDI_Get_Intra_Code_MPI_Comm()
            self.world_rank = self.mpi_world.Get_rank()

        # Accept a communicator to the driver code
        self.comm = MDI_Accept_Communicator()

        # Ensure that the molecule is using c1 symmetry
        self.molecule.reset_point_group('c1')
        self.molecule.fix_orientation(True)
        self.molecule.fix_com(True)
        self.molecule.update_geometry()

    def length_conversion(self):
        """ Obtain the conversion factor between the geometry specification units and bohr
        """
        unit_name = self.molecule.units()
        unit_conv = 1.0
        if unit_name == "Angstrom":
            unit_conv = MDI_Conversion_Factor("bohr","angstrom")
        elif unit_name == "Bohr":
            unit_conv = 1.0
        else:
            raise Exception('Unrecognized unit type: ' + str(unit_name))

        return unit_conv

    # respond to <NATOMS command
    def send_natoms(self):
        """ Send the number of atoms through MDI
        """
        natom = self.molecule.natom()
        MDI_Send(natom, 1, MDI_INT, self.comm)
        print("Natoms: " + str(natom))

    # respond to <COORDS command
    def send_coords(self):
        """ Send the nuclear coordinates through MDI
        """
        natom = self.molecule.natom()
        coords = [ self.molecule.xyz(iatom)[icoord] for iatom in range(natom) for icoord in range(3) ]
        MDI_Send(coords, 3*natom, MDI_DOUBLE, self.comm)
        print("Coords: " + str(coords))

    # respond to <CHARGES command
    def send_charges(self):
        """ Send the nuclear charges through MDI
        """
        natom = self.molecule.natom()
        charges = [ self.molecule.charge(iatom) for iatom in range(natom) ]
        MDI_Send(charges, natom, MDI_DOUBLE, self.comm)
        print("Charges: " + str(charges))

    # respond to <MASSES command
    def send_masses(self):
        """ Send the nuclear masses through MDI
        """
        natom = self.molecule.natom()
        masses = [ self.molecule.mass(iatom) for iatom in range(natom) ]
        MDI_Send(masses, natom, MDI_DOUBLE, self.comm)
        print("Masses: " + str(masses))

    # respond to <ELEMENTS command
    def send_elements(self):
        """ Send the atomic number of each nucleus through MDI
        """
        natom = self.molecule.natom()
        elements = [ self.molecule.true_atomic_number(iatom) for iatom in range(natom) ]
        MDI_Send(elements, natom, MDI_INT, self.comm)
        print("Elements: " + str(elements))

    # respond to <ENERGY command
    def send_energy(self):
        """ Send the total energy through MDI
        """
        MDI_Send(self.energy, 1, MDI_DOUBLE, self.comm)
        print("Energy: " + str(self.energy))

    # respond to <FORCES command
    def send_forces(self):
        """ Send the nuclear forces through MDI
        """
        natom = self.molecule.natom()
        force_matrix = psi4.driver.gradient(self.scf_method)
        unit_conv = 1.0
        forces = [ force_matrix.get(iatom,icoord)*unit_conv for iatom in range(natom) for icoord in range(3) ]
        MDI_Send(forces, 3*natom, MDI_DOUBLE, self.comm)
        print("Forces: " + str(forces))

    # respond to >CHARGES command
    def recv_charges(self):
        """ Receive a set of nuclear charges through MDI and assign them to the atoms in the current molecule
        """
        natom = self.molecule.natom()
        charges = MDI_Recv(natom, MDI_DOUBLE, self.comm)
        for iatom in range(natom):
            self.molecule.set_nuclear_charge(iatom, charges[iatom])

    # respond to >COORDS command
    def recv_coords(self):
        """ Receive a set of nuclear coordinates through MDI and assign them to the atoms in the current molecule
        """
        natom = self.molecule.natom()
        coords = MDI_Recv(3*natom, MDI_DOUBLE, self.comm)
        matrix = psi4.core.Matrix(natom, 3)
        for iatom in range(natom):
            matrix.set(iatom,0,coords[3*iatom+0])
            matrix.set(iatom,1,coords[3*iatom+1])
            matrix.set(iatom,2,coords[3*iatom+2])
        self.molecule.set_geometry(matrix)

    # respond to the >MASSES command
    def recv_masses(self):
        """ Receive a set of nuclear masses through MDI and assign them to the atoms in the current molecule
        """
        natom = self.molecule.natom()
        masses = MDI_Recv(natom, MDI_DOUBLE, self.comm)

        # assign the mass of all atoms, taking care to avoid ghost atoms
        jatom = 0
        for iatom in range(natom):
            while self.molecule.fcharge(jatom) == 0.0 and jatom < self.molecule.nallatom():
                jatom = jatom + 1
                if jatom >= self.molecule.nallatom():
                    raise Exception('Unexpected number of ghost atoms when receiving masses: ')
            self.molecule.set_mass(iatom, masses[jatom])
            jatom = jatom + 1

    # set a lattice of point charges
    def set_lattice_field(self):
        """ Set a field of lattice point charges using information received through MDI
        """
        self.lattice_field = psi4.QMMM()
        for ilat in range(self.nlattice):
            latx = self.clattice[3*ilat+0]
            laty = self.clattice[3*ilat+1]
            latz = self.clattice[3*ilat+2]
            self.lattice_field.extern.addCharge(self.lattice[ilat], latx, laty, latz)
        psi4.core.set_global_option_python('EXTERN', self.lattice_field.extern)

    # respond to >NLATTICE command
    def recv_nlattice(self):
        """ Receive the number of lattice point charges through MDI
        """
        self.nlattice = MDI_Recv(1, MDI_INT, self.comm)
        self.clattice = [ 0.0 for ilat in range(3*self.nlattice) ]
        self.lattice = [ 0.0 for ilat in range(self.nlattice) ]
        set_lattice_field(self.molecule)

    # respond to >CLATTICE command
    def recv_clattice(self):
        """ Receive the coordinates of a set of lattice point charges through MDI
        """
        self.clattice = MDI_Recv(3*self.nlattice, MDI_DOUBLE, self.comm)
        set_lattice_field(self.molecule)

    # respond to >LATTICE command
    def recv_lattice(self):
        """ Receive the charges of a set of lattice point charges through MDI
        """
        self.lattice = MDI_Recv(self.nlattice, MDI_DOUBLE, self.comm)
        set_lattice_field(self.molecule)

    # respond to SCF command
    def run_scf(self):
        """ Run an energy calculation
        """
        self.energy = psi4.energy(self.scf_method)

    # enter server mode, listening for commands from the driver
    def listen_for_commands(self):
        """ Receive commands through MDI and respond to them as defined by the MDI Standard
        """

        while True:
            if self.world_rank == 0:
                command = MDI_Recv_Command(self.comm)
            else:
                command = None
            if use_mpi4py:
                command = mpi_world.bcast(command, root=0)
            if self.world_rank == 0:
                core.print_out('\nMDI command received: ' + str(command) + ' \n')

            # respond to commands
            if command == "<NATOMS":
                self.send_natoms()
            elif command == "<COORDS":
                self.send_coords()
            elif command == "<CHARGES":
                self.send_charges()
            elif command == "<ELEMENTS":
                self.send_elements()
            elif command == "<MASSES":
                self.send_masses()
            elif command == "<ENERGY":
                self.send_energy()
            elif command == "<FORCES":
                self.send_forces()
            elif command == ">COORDS":
                self.recv_coords()
            elif command == ">NLATTICE":
                self.recv_nlattice()
            elif command == ">CLATTICE":
                self.recv_clattice()
            elif command == ">LATTICE":
                self.recv_lattice()
            elif command == ">MASSES":
                self.recv_masses()
            elif command == "SCF":
                self.run_scf()
            elif command == "EXIT":
                break
            else:
                raise Exception('Unrecognized command: ' + str(command))




def mdi_init(mdi_arguments):
    """ Initialize the MDI Library

    Arguments:
        mdi_arguments: MDI configuration options
    """
    mpi_world = None
    if use_mpi4py:
        mpi_world = MPI.COMM_WORLD
    MDI_Init(mdi_arguments, mpi_world)

def mdi(scf_method):
    """ Begin functioning as an MDI engine

    Arguments:
        scf_method: Method used when calculating energies or gradients
    """
    engine = MDIEngine(scf_method)
    engine.listen_for_commands()
