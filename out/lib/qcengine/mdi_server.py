"""Support for using QCEngine as an MDI engine.
For details regarding MDI, see https://molssi.github.io/MDI_Library/html/index.html.

"""
from typing import Any, Dict, List, Optional

import numpy as np
import qcelemental as qcel
from qcelemental.util import which_import

from .compute import compute

try:
    from mdi import (
        MDI_CHAR,
        MDI_COMMAND_LENGTH,
        MDI_DOUBLE,
        MDI_INT,
        MDI_MAJOR_VERSION,
        MDI_Accept_Communicator,
        MDI_Init,
        MDI_MPI_get_world_comm,
        MDI_Recv,
        MDI_Recv_Command,
        MDI_Register_Command,
        MDI_Register_Node,
        MDI_Send,
    )

    use_mdi = True
except ImportError:
    use_mdi = False

try:
    from mpi4py import MPI

    use_mpi4py = True
except ImportError:
    use_mpi4py = False


class MDIServer:
    def __init__(
        self,
        mdi_options: str,
        program: str,
        molecule,
        model,
        keywords,
        raise_error: bool = False,
        local_options: Optional[Dict[str, Any]] = None,
    ):
        """Initialize an MDIServer object for communication with MDI

        Parameters
        ----------
        mdi_options: str
            Options used during MDI initialization.
        program : str
            The program to execute the input with.
        molecule
            The initial state of the molecule.
        model
            The simulation model to use.
        keywords
            Program-specific keywords.
        raise_error : bool, optional
            Determines if compute should raise an error or not.
        local_options : Optional[Dict[str, Any]], optional
            A dictionary of local configuration options
        """

        if not use_mdi:
            raise Exception("Trying to run as an MDI engine, but the MDI Library was not found")

        if MDI_MAJOR_VERSION < 1:
            raise Exception("QCEngine requires version 1.0.0 or higher of the MDI Library")

        # Confirm that the MDI library has been located
        which_import("mdi", raise_error=True, raise_msg="Please install via 'conda install pymdi -c conda-forge'")

        # Initialize MDI
        MDI_Init(mdi_options)

        # Input variables
        self.molecule = molecule
        self.model = model
        self.keywords = keywords
        self.program = program
        self.raise_error = raise_error
        self.local_options = local_options

        # The MDI interface does not currently support multiple fragments
        if len(self.molecule.fragments) != 1:
            raise Exception("The MDI interface does not support multiple fragments")

        # Molecule charge and multiplicity
        self.total_charge = self.molecule.molecular_charge
        self.multiplicity = self.molecule.molecular_multiplicity

        # Flag to track whether the latest molecule specification has been validated
        self.molecule_validated = True

        # Output of most recent compute call
        self.compute_return = None

        # MPI variables
        self.mpi_world = None
        self.world_rank = 0

        # Get correct intra-code MPI communicator
        if use_mpi4py:
            self.mpi_world = MDI_MPI_get_world_comm()
            self.world_rank = self.mpi_world.Get_rank()

            # QCEngine does not currently support multiple MPI ranks
            if self.mpi_world.Get_size() != 1:
                MPI.COMM_WORLD.Abort()

        # Flag to stop listening for MDI commands
        self.stop_listening = False

        # Dictionary of all supported MDI commands
        self.commands = {
            "<@": self.send_node,
            "<NATOMS": self.send_natoms,
            "<COORDS": self.send_coords,
            "<ENERGY": self.send_energy,
            "<FORCES": self.send_forces,
            ">COORDS": self.recv_coords,
            "SCF": self.run_energy,
            "<ELEMENTS": self.send_elements,
            ">ELEMENTS": self.recv_elements,
            "<MASSES": self.send_masses,
            ">MASSES": self.recv_masses,
            "<TOTCHARGE": self.send_total_charge,
            ">TOTCHARGE": self.recv_total_charge,
            "<ELEC_MULT": self.send_multiplicity,
            ">ELEC_MULT": self.recv_multiplicity,
            "EXIT": self.stop,
        }

        # Register the @DEFAULT node
        MDI_Register_Node("@DEFAULT")

        # Register all supported commands
        for c in self.commands.keys():
            MDI_Register_Command("@DEFAULT", c)

        # Set the current node
        self.current_node = "@DEFAULT"

        # Accept a communicator to the driver code
        self.comm = MDI_Accept_Communicator()

    def update_molecule(self, key: str, value):
        """Update the molecule

        Parameters
        ----------
        key : str
            Key of the molecular element to update
        value
            Update value
        """
        if key == "molecular_charge" or key == "molecular_multiplicity":
            # In order to validate correctly, the charges and multiplicities must be set simultaneously
            try:
                self.molecule = qcel.models.Molecule(
                    **{
                        **self.molecule.dict(),
                        **{
                            "molecular_charge": self.total_charge,
                            "fragment_charges": [self.total_charge],
                            "molecular_multiplicity": self.multiplicity,
                            "fragment_multiplicities": [self.multiplicity],
                        },
                    }
                )
                self.molecule_validated = True
            except qcel.exceptions.ValidationError:
                # The molecule didn't validate correctly, but a future >TOTCHARGE or >ELEC_MULT command might fix it
                self.molecule_validated = False
        else:
            try:
                self.molecule = qcel.models.Molecule(**{**self.molecule.dict(), **{key: value}})
                self.molecule_validated = True
            except qcel.exceptions.ValidationError:
                if self.molecule_validated:
                    # This update caused the validation error
                    raise Exception("MDI command caused a validation error")

    # Respond to the <@ command
    def send_node(self) -> str:
        """Send the name of the current node through MDI

        Returns
        -------
        node : str
            Name of the current node
        """
        node = self.current_node
        MDI_Send(node, MDI_COMMAND_LENGTH, MDI_CHAR, self.comm)
        return node

    # Respond to the <NATOMS command
    def send_natoms(self) -> int:
        """Send the number of atoms through MDI

        Returns
        -------
        natom : int
            Number of atoms
        """
        natom = len(self.molecule.geometry)
        MDI_Send(natom, 1, MDI_INT, self.comm)
        return natom

    # Respond to the <COORDS command
    def send_coords(self) -> np.ndarray:
        """Send the nuclear coordinates through MDI

        Returns
        -------
        coords : np.ndarray
            Nuclear coordinates
        """
        natom = len(self.molecule.geometry)

        coords = np.reshape(self.molecule.geometry, (3 * natom))
        MDI_Send(coords, 3 * natom, MDI_DOUBLE, self.comm)

        return coords

    # Respond to the >COORDS command
    def recv_coords(self, coords: Optional[np.ndarray] = None) -> None:
        """Receive a set of nuclear coordinates through MDI and assign them to the atoms in the current molecule

        Parameters
        ----------
        coords : np.ndarray, optional
            New nuclear coordinates. If None, receive through MDI.
        """
        natom = len(self.molecule.geometry)
        if coords is None:
            coords = np.zeros(3 * natom)
            MDI_Recv(3 * natom, MDI_DOUBLE, self.comm, buf=coords)
        new_geometry = np.reshape(coords, (-1, 3))
        self.molecule = qcel.models.Molecule(**{**self.molecule.dict(), **{"geometry": new_geometry}})

    # Respond to the <ENERGY command
    def send_energy(self) -> float:
        """Send the total energy through MDI

        Returns
        -------
        energy : float
            Energy of the system
        """
        # Ensure that the molecule currently passes validation
        if not self.molecule_validated:
            raise Exception("MDI attempting to compute energy on an unvalidated molecule")
        self.run_energy()
        energy = self.compute_return.return_result
        MDI_Send(energy, 1, MDI_DOUBLE, self.comm)
        return energy

    # Respond to the <FORCES command
    def send_forces(self) -> np.ndarray:
        """Send the nuclear forces through MDI

        Returns
        -------
        forces : np.ndarray
            Forces on the nuclei
        """
        # Ensure that the molecule currently passes validation
        if not self.molecule_validated:
            raise Exception("MDI attempting to compute gradients on an unvalidated molecule")

        input = qcel.models.AtomicInput(
            molecule=self.molecule, driver="gradient", model=self.model, keywords=self.keywords
        )
        self.compute_return = compute(
            input_data=input, program=self.program, raise_error=self.raise_error, local_options=self.local_options
        )

        forces = self.compute_return.return_result
        MDI_Send(forces, len(forces), MDI_DOUBLE, self.comm)
        return forces

    # Respond to the SCF command
    def run_energy(self) -> None:
        """Run an energy calculation"""
        input = qcel.models.AtomicInput(
            molecule=self.molecule, driver="energy", model=self.model, keywords=self.keywords
        )
        self.compute_return = compute(
            input_data=input, program=self.program, raise_error=self.raise_error, local_options=self.local_options
        )

    # Respond to the <ELEMENTS command
    def send_elements(self):
        """Send the atomic number of each nucleus through MDI

        Returns
        -------
        elements : :obj:`list` of :obj:`int`
            Element of each atom
        """
        natom = len(self.molecule.geometry)
        elements = [qcel.periodictable.to_atomic_number(self.molecule.symbols[iatom]) for iatom in range(natom)]
        MDI_Send(elements, natom, MDI_INT, self.comm)
        return elements

    # Respond to the >ELEMENTS command
    def recv_elements(self, elements: Optional[List[int]] = None):
        """Receive a set of atomic numbers through MDI and assign them to the atoms in the current molecule

        Parameters
        ----------
        elements : :obj:`list` of :obj:`int`, optional
            New element numbers. If None, receive through MDI.
        """
        natom = len(self.molecule.geometry)
        if elements is None:
            elements = MDI_Recv(natom, MDI_INT, self.comm)

        for iatom in range(natom):
            self.molecule.symbols[iatom] = qcel.periodictable.to_symbol(elements[iatom])

        return elements

    # Respond to the <MASSES command
    def send_masses(self) -> np.ndarray:
        """Send the nuclear masses through MDI

        Returns
        -------
        masses : np.ndarray
            Atomic masses
        """
        natom = len(self.molecule.geometry)
        masses = self.molecule.masses
        MDI_Send(masses, natom, MDI_DOUBLE, self.comm)
        return masses

    # Respond to the >MASSES command
    def recv_masses(self, masses: Optional[List[float]] = None) -> None:
        """Receive a set of nuclear masses through MDI and assign them to the atoms in the current molecule

        Parameters
        ----------
        masses : :obj:`list` of :obj:`float`, optional
            New nuclear masses. If None, receive through MDI.
        """
        natom = len(self.molecule.geometry)
        if masses is None:
            masses = MDI_Recv(natom, MDI_DOUBLE, self.comm)
        self.update_molecule("masses", masses)

    # Respond to the <TOTCHARGE command
    def send_total_charge(self) -> float:
        """Send the total system charge through MDI

        Returns
        -------
        charge : float
            Total charge of the system
        """
        charge = self.molecule.molecular_charge
        MDI_Send(charge, 1, MDI_DOUBLE, self.comm)
        return charge

    # Respond to the >TOTCHARGE command
    def recv_total_charge(self, charge: Optional[float] = None) -> None:
        """Receive the total system charge through MDI

        Parameters
        ----------
        charge : float, optional
            New charge of the system. If None, receive through MDI.
        """
        if charge is None:
            charge = MDI_Recv(1, MDI_DOUBLE, self.comm)
        self.total_charge = charge

        # Allow a validation error here, because a future >ELEC_MULT command might resolve it
        try:
            self.update_molecule("molecular_charge", self.total_charge)
        except qcel.exceptions.ValidationError:
            pass

    # Respond to the <ELEC_MULT command
    def send_multiplicity(self) -> int:
        """Send the electronic multiplicity through MDI

        Returns
        -------
        multiplicity : int
            Multiplicity of the system
        """
        multiplicity = self.molecule.molecular_multiplicity
        MDI_Send(multiplicity, 1, MDI_INT, self.comm)
        return multiplicity

    # Respond to the >ELEC_MULT command
    def recv_multiplicity(self, multiplicity: Optional[int] = None) -> None:
        """Receive the electronic multiplicity through MDI

        Parameters
        ----------
        multiplicity : int, optional
            New multiplicity of the system. If None, receive through MDI.
        """
        if multiplicity is None:
            multiplicity = MDI_Recv(1, MDI_INT, self.comm)
        self.multiplicity = multiplicity

        # Allow a validation error here, because a future >TOTCHARGE command might resolve it
        try:
            self.update_molecule("molecular_multiplicity", self.multiplicity)
        except qcel.exceptions.ValidationError:
            pass

    # Respond to the EXIT command
    def stop(self) -> None:
        """Stop listening for MDI commands"""
        self.stop_listening = True

    # Enter server mode, listening for commands from the driver
    def start(self) -> None:
        """Receive commands through MDI and respond to them as defined by the MDI Standard"""

        while not self.stop_listening:
            if self.world_rank == 0:
                command = MDI_Recv_Command(self.comm)
            else:
                command = None
            if use_mpi4py:
                command = self.mpi_world.bcast(command, root=0)
            if self.world_rank == 0:
                print("MDI command received: " + str(command))

            # Search for this command in self.commands
            found_command = False
            for supported_command in self.commands:
                if not found_command and command == supported_command:
                    # Run the function corresponding to this command
                    self.commands[supported_command]()
                    found_command = True
            if not found_command:
                raise Exception("Unrecognized command: " + str(command))
