"""
Calls the Molpro executable.
"""

import string
from typing import Any, Dict, List, Optional, Set, Tuple
from xml.etree import ElementTree as ET

from qcelemental.models import AtomicResult
from qcelemental.util import parse_version, safe_version, which

from ..exceptions import InputError, UnknownError
from ..util import execute
from .model import ProgramHarness


class MolproHarness(ProgramHarness):
    _defaults: Dict[str, Any] = {
        "name": "Molpro",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    # Set of implemented dft functionals in Molpro according to dfunc.registry (version 2019.2)
    # fmt: off
    _dft_functionals: Set[str] = {
        "B86MGC", "B86R", "B86", "B88C", "B88", "B95", "B97DF", "B97RDF", "BR", "BRUEG", "BW", "CS1", "CS2", "DIRAC",
        "ECERFPBE", "ECERF", "EXACT", "EXERFPBE", "EXERF", "G96", "HCTH120", "HCTH147", "HCTH93", "HJSWPBEX", "LTA",
        "LYP", "M052XC", "M052XX", "M05C", "M05X", "M062XC", "M062XX", "M06C", "M06HFC", "M06HFX", "M06LC", "M06LX",
        "M06X", "M12C", "MK00B", "MK00", "P86", "PBEC", "PBESOLC", "PBESOLX", "PBEXREV", "PBEX", "PW86", "PW91C",
        "PW91X", "PW92C", "STEST", "TFKE", "TH1", "TH2", "TH3", "TH4", "THGFCFO", "THGFCO", "THGFC", "THGFL", "TPSSC",
        "TPSSX", "VSXC", "VW", "VWN3", "VWN5", "XC-M05-2X", "XC-M05", "XC-M06-2X", "XC-M06-HF", "XC-M06-L", "XC-M06",
        "XC-M08-HX", "XC-M08-SO", "XC-M11-L", "XC-SOGGA11", "XC-SOGGA11-X", "XC-SOGGA", "FRMTST", "LHF", "TLHF",
        "LXBECKE", "ELP", "NULL", "YTEST", "TREF2", "TREF", "TTEST", "GLE", "GREEN", "SRB88", "SRLYP", "LB94", "EI",
        "SAOP", "USER", "INE", "ECERF2", "ECERFINTER", "ECERFLOCAL2", "ECERFLOCAL", "EXERFLOCAL", "FC", "FCFO", "FCO",
        "FL", "XC-M11", "PBEXANAL", "PBECANAL", "PBESOLCANAL", "PBESOLXANAL", "EXSRLDA", "EXSRLPBE", "ECSRLPBE",
        "ECSRLLPBE", "ECSQRTLPBE", "ECMUDIVLPBE", "EXERFPHS", "ECLERFMUPBE", "ECERFERFCPBE", "ECSQRTLDA", "REVPBEX",
        "B", "B-LYP", "BLYP", "B-P", "BP86", "B-VWN", "B3LYP", "B3LYP3", "B3LYP5", "B88X", "B97", "B97R", "BECKE",
        "BH-LYP", "CS", "D", "HFB", "HFS", "LDA", "LSDAC", "LSDC", "KYP88", "MM05", "MM05-2X", "MM06", "MM06-2X",
        "MM06-L", "MM06-HF", "PBE", "PBE0", "PBE0MOL", "PBEREV", "PW91", "S", "S-VWN", "SLATER", "VS99", "VWN",
        "VWN80", 'M05', "M05-2X", "M06", "M06-2X", "M06-L", "M06-HF", "M08-HX", "M08-SO", "M11-L", "TPSS", "TPSSH",
        "M12HFC", "HJSWPBE", "HJSWPBEH", "TCSWPBE", "PBESOL"
    }
    # fmt: on

    # Different open-shell scenarios:
    # - Restricted references
    #       RHF-RMP2 (like RHF-UMP2 using Molpro's naming convention for CC methods)
    #       RHF-UCCSD
    #       RHF-UCCSD(T)
    #       RHF-RCCSD
    #       RHF-RCCSD(T)
    # - Unrestricted references (Only supported up to UMP2, no CC support)
    #       UHF-UMP2
    # NOTE: Unrestricted SCF methods must be specified by using keyword reference
    _hf_methods: Set[str] = {"HF", "RHF"}
    _restricted_post_hf_methods: Set[str] = {"MP2", "CCSD", "CCSD(T)"}  # RMP2, RCCSD, RCCSD(T)}
    # TODO Add keyword to specify unrestricted for WF method
    # _unrestricted_post_hf_methods: Set[str] = {"UMP2", "UCCSD", "UCCSD(T)"}
    _post_hf_methods: Set[str] = {*_restricted_post_hf_methods}

    class Config(ProgramHarness.Config):
        pass

    def found(self, raise_error: bool = False) -> bool:
        return which(
            "molpro", return_bool=True, raise_error=raise_error, raise_msg="Please install via https://www.molpro.net/"
        )

    # TODO Consider changing this to use molpro --version instead of performing a full execute
    def get_version(self) -> str:
        self.found(raise_error=True)

        name_space = {"molpro_uri": "http://www.molpro.net/schema/molpro-output"}
        which_prog = which("molpro")
        if which_prog not in self.version_cache:
            success, output = execute(
                [which_prog, "version.inp", "-d", ".", "-W", "."],
                infiles={"version.inp": ""},
                outfiles=["version.out", "version.xml"],
            )

            if success:
                tree = ET.ElementTree(ET.fromstring(output["outfiles"]["version.xml"]))
                root = tree.getroot()
                version_tree = root.find("molpro_uri:job/molpro_uri:platform/molpro_uri:version", name_space)
                if version_tree is None:
                    # some older schema
                    name_space = {"molpro_uri": "http://www.molpro.net/schema/molpro2006"}
                    version_tree = root.find("molpro_uri:job//molpro_uri:version", name_space)
                if version_tree is None:
                    return safe_version("0.0.0")
                year = version_tree.attrib["major"]
                minor = version_tree.attrib["minor"]
                molpro_version = year + "." + minor
                self.version_cache[which_prog] = safe_version(molpro_version)

        return self.version_cache[which_prog]

    def compute(self, input_data: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        """
        Run Molpro
        """
        # Check if Molpro executable is found
        self.found(raise_error=True)

        # Check Molpro version
        if parse_version(self.get_version()) < parse_version("2018.1"):
            raise TypeError("Molpro version '{}' not supported".format(self.get_version()))

        # Setup the job
        job_inputs = self.build_input(input_data, config)

        # Run Molpro
        exe_success, proc = self.execute(job_inputs)

        # Determine whether the calculation succeeded
        if exe_success:
            # If execution succeeded, collect results
            result = self.parse_output(proc["outfiles"], input_data)
            return result
        else:
            # Return UnknownError for error propagation
            return UnknownError(proc["stderr"])

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_infiles: Optional[Dict[str, str]] = None,
        extra_outfiles: Optional[List[str]] = None,
        as_binary: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        scratch_messy: bool = False,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:
        """
        For option documentation go look at qcengine/util.execute
        """

        # Collect all input files and update with extra_infiles
        infiles = inputs["infiles"]
        if extra_infiles is not None:
            infiles.update(extra_infiles)

        # Collect all output files and update with extra_outfiles
        outfiles = ["dispatch.out", "dispatch.xml"]
        if extra_outfiles is not None:
            outfiles.extend(extra_outfiles)

        # Replace commands with extra_commands if present
        commands = inputs["commands"]
        if extra_commands is not None:
            commands = extra_commands

        # Run the Molpro program
        exe_success, proc = execute(
            commands,
            infiles=infiles,
            outfiles=outfiles,
            as_binary=as_binary,
            scratch_name=scratch_name,
            scratch_directory=inputs["scratch_directory"],
            scratch_messy=scratch_messy,
            timeout=timeout,
        )
        return exe_success, proc

    def build_input(
        self, input_model: "AtomicInput", config: "TaskConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:
        if template is None:
            input_file = []

            # Resolving keywords
            caseless_keywords = {k.lower(): v for k, v in input_model.keywords.items()}
            unrestricted = False
            # Following Molpro default that ROHF is done for open-shell calculations unless unrestricted is specified
            if "reference" in caseless_keywords and caseless_keywords["reference"] == "unrestricted":
                unrestricted = True

            # Memory is in megawords per core for Molpro
            memory_mw_core = int(config.memory * (1024 ** 3) / 8e6 / config.ncores)
            input_file.append("memory,{},M".format(memory_mw_core))
            input_file.append("")

            # Write the geom
            xyz_block = input_model.molecule.to_string(dtype="molpro", units="Bohr")
            input_file.append(xyz_block)

            # Write the basis set
            input_file.append("basis={")
            input_file.append(f"default,{input_model.model.basis}")
            input_file.append("}")
            input_file.append("")

            # Determine what SCF type (restricted vs. unrestricted)
            hf_type = "RHF"
            dft_type = "RKS"
            if unrestricted:
                hf_type = "UHF"
                dft_type = "UKS"

            # Write energy call
            energy_call = []
            # If post-hf method is called then make sure to write a HF call first
            if input_model.model.method.upper() in self._post_hf_methods:  # post SCF case
                energy_call.append(f"{{{hf_type}}}")
                energy_call.append("")
                energy_call.append(f"{{{input_model.model.method}}}")
            # If DFT call make sure to write {rks,method}
            elif input_model.model.method.upper() in self._dft_functionals:  # DFT case
                energy_call.append(f"{{{dft_type},{input_model.model.method}}}")
            elif input_model.model.method.upper() in self._hf_methods:  # HF case
                energy_call.append(f"{{{hf_type}}}")
            else:
                raise InputError(f"Method {input_model.model.method} not implemented for Molpro.")

            # Write appropriate driver call
            if input_model.driver == "energy":
                input_file.extend(energy_call)
            elif input_model.driver == "gradient":
                input_file.extend(energy_call)
                input_file.append("")
                input_file.append("{force}")
            else:
                raise InputError(f"Driver {input_model.driver} not implemented for Molpro.")

            input_file = "\n".join(input_file)
        else:
            # Some of the potential different template options
            # (A) ordinary build_input (need to define a base template)
            # (B) user wants to add stuff after normal template (A)
            # (C) user knows their domain language (doesn't use any QCSchema quantities)

            # # Build dictionary for substitute
            # sub_dict = {
            #     "method": input_model.model.method,
            #     "basis": input_model.model.basis,
            #     "charge": input_model.molecule.molecular_charge
            # }

            # Perform substitution to create input file
            str_template = string.Template(template)
            input_file = str_template.substitute()

        return {
            "commands": ["molpro", "dispatch.mol", "-d", ".", "-W", ".", "-n", str(config.ncores)],
            "infiles": {"dispatch.mol": input_file},
            "scratch_directory": config.scratch_directory,
            "input_result": input_model.copy(deep=True),
        }

    def parse_output(self, outfiles: Dict[str, str], input_model: "AtomicInput") -> "AtomicResult":
        tree = ET.ElementTree(ET.fromstring(outfiles["dispatch.xml"]))
        root = tree.getroot()
        # print(root.tag)

        # TODO Read information from molecule tag
        #      - cml:molecule, cml:atomArray (?)
        #      - basisSet
        #       - Be aware of symmetry. Might only be able to support if symmetry,nosym
        #      - orbitals
        #       - Be aware of symmetry. Might only be able to support if symmetry,nosym
        # NOTE: Spherical basis set ordering in Molpro (with no symmetry)
        # S -->  0
        # P --> +1, -1, 0
        # D -->  0, -2, +1, +2, -1
        # F --> +1, -1, 0, +3, -2, -3, +2
        # G -->  0, -2, +1, +4, -1, +2, -4, +3, -3
        # H --> +1, -1, +2, +3, -4, -3, +4, -5, 0, +5, -2
        # I --> +6, -2, +5, +4, -5, +2, -6, +3, -4, 0, -3, -1, +1
        properties = {}
        extras = {}
        name_space = {"molpro_uri": "http://www.molpro.net/schema/molpro-output"}

        # Molpro commands map
        molpro_map = {
            "Energy": {
                "HF": "scf_total_energy",
                "RHF": "scf_total_energy",
                "UHF": "scf_total_energy",
                "KS": "scf_total_energy",
                "RKS": "scf_total_energy",
                "UKS": "scf_total_energy",
            },
            "total energy": {
                "MP2": "mp2_total_energy",
                "CCSD": "ccsd_total_energy",
                "CCSD(T)": "ccsd_prt_pr_total_energy",
            },
            "correlation energy": {
                "MP2": "mp2_correlation_energy",
                "CCSD": "ccsd_correlation_energy",
                "CCSD(T)": "ccsd_prt_pr_correlation_energy",  # Need both CCSD(T) and Total
                "Total": "ccsd_prt_pr_correlation_energy",  # Total corresponds to CCSD(T) correlation energy
            },
            "singlet pair energy": {"MP2": "mp2_singlet_pair_energy", "CCSD": "ccsd_singlet_pair_energy"},
            "triplet pair energy": {"MP2": "mp2_triplet_pair_energy", "CCSD": "ccsd_triplet_pair_energy"},
            "Dipole moment": {
                "HF": "scf_dipole_moment",
                "RHF": "scf_dipole_moment",
                "UHF": "scf_dipole_moment",
                "KS": "scf_dipole_moment",
                "RKS": "scf_dipole_moment",
                "UKS": "scf_dipole_moment",
                "MP2": "mp2_dipole_moment",
                "CCSD": "ccsd_dipole_moment",
                "CCSD(T)": "ccsd_prt_pr_dipole_moment",
            },
        }

        # Started adding basic support for local correlation methods in Molpro
        molpro_extras_map = {
            "total energy": {
                "LMP2": "local_mp2_total_energy",
                "LCCSD": "local_ccsd_total_energy",
                "LCCSD(T0)": "local_ccsd_prt0_pr_total_energy",
                "LCCSD(T)": "local_ccsd_prt_pr_total_energy",
            },
            "correlation energy": {"LMP2": "local_mp2_correlation_energy", "LCCSD": "local_ccsd_correlation_energy"},
            "singlet pair energy": {"LMP2": "local_mp2_singlet_pair_energy", "LCCSD": "local_ccsd_singlet_pair_energy"},
            "triplet pair energy": {"LMP2": "local_mp2_triplet_pair_energy", "LCCSD": "local_ccsd_triplet_pair_energy"},
            "singles energy": {"LCCSD": "local_ccsd_singles_energy"},
            # "strong pair energy": {
            #     "LCCSD": "local_ccsd_strong_pair_energy"
            # },
            # "weak pair energy": {
            #     "LMP2": "local_mp2_weak_pair_energy"
            # }
        }

        # Molpro variables map used for quantities not found in the command map
        molpro_variable_map = {
            "_ENUC": "nuclear_repulsion_energy",
            "_DFTFUN": "scf_xc_energy",
            "_NELEC": ["calcinfo_nalpha", "calcinfo_nbeta"]
            # "_EMP2_SCS": "scs_mp2_total_energy"
        }

        # Process data in molpro_map by looping through each jobstep
        # The jobstep tag in Molpro contains output from commands (e.g. {HF}, {force})
        for jobstep in root.findall("molpro_uri:job/molpro_uri:jobstep", name_space):
            command = jobstep.attrib["command"]
            if "FORCE" in command:  # Grab gradient
                for child in jobstep.findall("molpro_uri:gradient", name_space):
                    # Stores gradient as a single list where the ordering is [1x, 1y, 1z, 2x, 2y, 2z, ...]
                    properties["gradient"] = [float(x) for x in child.text.split()]
            else:  # Grab energies and dipole moment
                for child in jobstep.findall("molpro_uri:property", name_space):
                    property_name = child.attrib["name"]
                    property_method = child.attrib["method"]
                    value = child.attrib["value"]
                    if property_name in molpro_map and property_method in molpro_map[property_name]:
                        if property_name == "Dipole moment":
                            properties[molpro_map[property_name][property_method]] = [float(x) for x in value.split()]
                        else:
                            properties[molpro_map[property_name][property_method]] = float(value)
                    elif property_name in molpro_extras_map and property_method in molpro_extras_map[property_name]:
                        extras[molpro_extras_map[property_name][property_method]] = float(value)

        # Convert triplet and singlet pair correlation energies to opposite-spin and same-spin correlation energies
        if "mp2_singlet_pair_energy" in properties and "mp2_triplet_pair_energy" in properties:
            properties["mp2_same_spin_correlation_energy"] = (2.0 / 3.0) * properties["mp2_triplet_pair_energy"]
            properties["mp2_opposite_spin_correlation_energy"] = (1.0 / 3.0) * properties[
                "mp2_triplet_pair_energy"
            ] + properties["mp2_singlet_pair_energy"]
            del properties["mp2_singlet_pair_energy"]
            del properties["mp2_triplet_pair_energy"]

        if "ccsd_singlet_pair_energy" in properties and "ccsd_triplet_pair_energy" in properties:
            properties["ccsd_same_spin_correlation_energy"] = (2.0 / 3.0) * properties["ccsd_triplet_pair_energy"]
            properties["ccsd_opposite_spin_correlation_energy"] = (1.0 / 3.0) * properties[
                "ccsd_triplet_pair_energy"
            ] + properties["ccsd_singlet_pair_energy"]
            del properties["ccsd_singlet_pair_energy"]
            del properties["ccsd_triplet_pair_energy"]

        # Process data in molpro_variable_map
        # Note: For the DFT case molecule_method is the name of the functional plus R or U in front
        molecule = root.find("molpro_uri:job/molpro_uri:molecule", name_space)
        molecule_method = molecule.attrib["method"]
        molecule_final_energy = float(molecule.attrib["energy"])  # Energy from the molecule tag in case its needed
        # Loop over each variable under the variables tag to grab additional info from molpro_variable_map
        for variable in molecule.findall("molpro_uri:variables/molpro_uri:variable", name_space):
            variable_name = variable.attrib["name"]
            if variable_name in molpro_variable_map:
                if variable_name == "_NELEC":
                    nelec = int(float(variable[0].text))
                    nunpaired = input_model.molecule.molecular_multiplicity - 1
                    nbeta = (nelec - nunpaired) // 2
                    nalpha = nelec - nbeta
                    properties[molpro_variable_map[variable_name][0]] = nalpha
                    properties[molpro_variable_map[variable_name][1]] = nbeta
                else:
                    properties[molpro_variable_map[variable_name]] = float(variable[0].text)

        # Process basis set data
        basis_set = root.find("molpro_uri:job/molpro_uri:molecule/molpro_uri:basisSet", name_space)
        nbasis = int(basis_set.attrib["length"])
        # angular_type = basis_set.attrib['angular']  # cartesian vs spherical
        properties["calcinfo_nbasis"] = nbasis

        # Grab the method from input
        method = input_model.model.method.upper()

        # Determining the final energy
        # Throws an error if the energy isn't found for the method specified from the input_model.
        if method in molpro_map["total energy"].keys() and molpro_map["total energy"][method] in properties:
            final_energy = properties[molpro_map["total energy"][method]]
        elif method in molpro_map["Energy"].keys() and molpro_map["Energy"][method] in properties:
            final_energy = properties[molpro_map["Energy"][method]]
        else:
            # Back up method for determining final energy if not already present in properties
            # Use the total energy from the molecule tag if it matches the input method
            # if input_model.model.method.upper() in molecule_method:
            if method in molecule_method:
                final_energy = molecule_final_energy
                if method in self._post_hf_methods:
                    properties[molpro_map["total energy"][method]] = molecule_final_energy
                    properties[molpro_map["correlation energy"][method]] = (
                        molecule_final_energy - properties["scf_total_energy"]
                    )
                elif method in self._dft_functionals:
                    properties[molpro_map["Energy"]["KS"]] = molecule_final_energy
                elif method in self._hf_methods:
                    properties[molpro_map["Energy"][method]] = molecule_final_energy
            else:
                raise KeyError(f"Could not find {method} total energy")

        # Initialize output_data by copying over input_model.dict()
        output_data = input_model.dict()

        # Determining return_result
        if input_model.driver == "energy":
            output_data["return_result"] = final_energy
        elif input_model.driver == "gradient":
            output_data["return_result"] = properties.pop("gradient")

        # Final output_data assignments needed for the AtomicResult object
        output_data["properties"] = properties
        output_data["extras"].update(extras)
        output_data["schema_name"] = "qcschema_output"
        output_data["stdout"] = outfiles["dispatch.out"]
        output_data["success"] = True

        return AtomicResult(**output_data)
