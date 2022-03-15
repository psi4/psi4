"""
Calls the Psi4 executable.
"""
import os
from typing import Any, Dict, List, Optional, Tuple

from qcelemental.models import AtomicResult
from qcelemental.util import which

from ..exceptions import InputError, UnknownError
from ..util import execute
from .model import ProgramHarness


class MopacHarness(ProgramHarness):

    _defaults = {
        "name": "MOPAC",
        "scratch": True,  # Input/output file
        "thread_safe": True,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    def __init__(self, **kwargs):
        extras = {  # All units taken from within MOPAC
            "bohr_to_angstroms": 0.5291772083,
            "hartree_to_ev": 27.2113834,
            "ev_to_kcalmol": 23.060529,
        }
        extras["au_to_debye"] = 2.99792458e10 * 1.602176462e0 * 1e-10 * extras["bohr_to_angstroms"]
        extras["hartree_to_kcalmol"] = extras["hartree_to_ev"] * extras["ev_to_kcalmol"]

        kwargs["extras"] = extras
        super().__init__(**kwargs)

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which(
            "mopac", return_bool=True, raise_error=raise_error, raise_msg="Please install via http://openmopac.net."
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        # Not really possible to pull at the moment, MolSSI will add a version ability
        return "2016"

    def compute(self, input_model: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        self.found(raise_error=True)

        exec_command = self.build_input(input_model, config)

        success, output = self.execute(exec_command)

        ret = self.parse_output(output["outfiles"], input_model)

        # Determine whether the calculation succeeded
        if success:
            # raise Exception()
            return ret
        else:
            return UnknownError(ret["stderr"])

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_infiles: Optional[Dict[str, str]] = None,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        scratch_messy: bool = False,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:
        """
        For option documentation go look at qcengine/util.execute
        """

        infiles = inputs["infiles"]
        if extra_infiles is not None:
            infiles.update(extra_infiles)

        outfiles = inputs["outfiles"]
        if extra_outfiles is not None:
            outfiles.extend(extra_outfiles)

        command = inputs["command"]
        if extra_commands is not None:
            command.extend(extra_commands)

        exe_success, proc = execute(
            command,
            infiles=infiles,
            outfiles=outfiles,
            scratch_directory=inputs["scratch_directory"],
            scratch_name=scratch_name,
            scratch_messy=scratch_messy,
            timeout=timeout,
            environment=inputs.get("environment", None),
        )
        return exe_success, proc

    def build_input(
        self, input_model: "AtomicInput", config: "TaskConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:

        if template is not None:
            raise KeyError("MOPAC does not currently support input templates.")

        method = input_model.model.method.lower()
        if method not in {
            "mndo",
            "am1",
            "pm3",
            "rm1",
            "mndod",
            "pm6",
            "pm6-d3",
            "pm6-dh+",
            "pm6-dh2",
            "pm6-dh2x",
            "pm6-d3h4",
            "pm6-3dh4x",
            "pm7",
            "pm7-ts",
        }:
            raise InputError(f"MOPAC does not have method: {method.upper()}")

        if input_model.driver not in ["energy", "gradient"]:
            raise InputError(f"Driver {input_model.driver} not implemented for MOPAC.")

        input_file = []

        keywords = {"ITER": 100, "PULAY": True}
        keywords.update({k.upper(): v for k, v in input_model.keywords.items()})

        if keywords["PULAY"]:
            pulay = " PULAY"
        else:
            pulay = " "

        unknown_keywords = keywords.keys() - {"ITER", "PULAY"}
        if unknown_keywords:
            raise InputError(f"Unknown keywords given to MOPAC: {unknown_keywords}")

        # 1SCF says not to compute a geometry optimization, always compute a gradient (free), and dump the aux file
        input_file.append(
            f"{method.upper()} "
            f"CHARGE={input_model.molecule.molecular_charge} "
            f"MS={(input_model.molecule.molecular_multiplicity-1)/2}&"
        )
        input_file.append(f"1SCF GRADIENTS AUX(PRECISION=9, XP, XS, XW) A0{pulay} ITER={keywords['ITER']}&")

        # input_file[-1] = input_file[-1][:-1]
        input_file.append("")

        mol = input_model.molecule

        for x in range(len(mol.symbols)):
            input_file.append(
                f"{mol.symbols[x]}  {mol.geometry[x, 0]:17.12f}  {mol.geometry[x, 1]:17.12f}  {mol.geometry[x, 2]:17.12f}"
            )

        env = os.environ.copy()
        env["MKL_NUM_THREADS"] = str(config.ncores)
        env["OMP_NUM_THREADS"] = str(config.ncores)

        return {
            "command": ["mopac", "dispatch.mop"],
            "infiles": {"dispatch.mop": "\n".join(input_file)},
            "outfiles": {"dispatch.out", "dispatch.aux"},
            "scratch_directory": config.scratch_directory,
            "environment": env,
        }

    def parse_output(self, outfiles: Dict[str, str], input_model: "AtomicInput") -> "AtomicResult":

        keep_keys = {
            "heat_of_formation",
            "energy_electronic",
            "energy_nuclear",
            "gradient_norm",
            "dip_vec",
            "spin_component",
            "total_spin",
            "molecular_weight",
            "molecular_weight",
            "total_energy",
            "gradients",
            "mopac_version",
            "atom_charges",
            "point_group",
        }

        # Convert back to atomic units
        conversions = {
            "KCAL/MOL": 1 / self.extras["hartree_to_kcalmol"],
            "KCAL/MOL/ANGSTROM": self.extras["bohr_to_angstroms"] / self.extras["hartree_to_kcalmol"],
            "EV": 1 / self.extras["hartree_to_ev"],
            "DEBYE": 1 / self.extras["au_to_debye"],
            "AMU": 1,
            None: 1,
        }

        data = {}
        last_key = None

        # Parse the weird structure
        if outfiles["dispatch.aux"] is None:
            error = "An unknown error occured and no results were captured."
            if outfiles["dispatch.out"] is not None:
                error = outfiles["dispatch.out"]
            raise UnknownError(error)

        for line in outfiles["dispatch.aux"].splitlines():
            if ("START" in line) or ("END" in line) or ("#" in line):
                continue

            if "=" in line:

                # Primary split
                key, value = line.split("=", 1)

                # Format key, may have units
                # IONIZATION_POTENTIAL:EV
                # GRADIENTS:KCAL/MOL/ANGSTROM[09]
                key_list = key.split(":", 1)
                if len(key_list) == 1:
                    key, units = key_list[0], None
                else:
                    key, units = key.split(":", 1)

                # Pop off [xx] items
                if units and "[" in units:
                    units, _ = units.split("[", 1)

                if "[" in key:
                    key, _ = key.split("[", 1)

                key = key.strip().lower()
                last_key = key

                # Skip keys that are not useful
                if key not in keep_keys:
                    last_key = None
                    continue

                # 1D+3 -> 1E3 conversion
                cf = conversions[units]

                value = value.strip().replace("D+", "E+").replace("D-", "E-")
                if ("E+" in value) or ("E-" in value):
                    if value.count("E") > 1:
                        value = [float(x) * cf for x in value.split()]
                    else:
                        value = float(value) * cf

                if value == "":
                    value = []

                data[key] = (cf, value)
            else:
                if last_key is None:
                    continue

                cf = data[last_key][0]
                data[last_key][1].extend([float(x) * cf for x in line.split()])

        data = {k: v[1] for k, v in data.items()}
        if ("gradients" not in data) or ("mopac_version" not in data):
            raise UnknownError("Could not correctly parse the MOPAC output file.")

        gradient = data.pop("gradients")

        output = input_model.dict()
        output["provenance"] = {"creator": "mopac", "version": data.pop("mopac_version")}

        output["properties"] = {}
        output["properties"]["return_energy"] = data["heat_of_formation"]

        output["extras"].update(data)

        if input_model.driver == "energy":
            output["return_result"] = data["heat_of_formation"]
        else:
            output["return_result"] = gradient

        output["stdout"] = outfiles["dispatch.out"]
        output["success"] = True

        return AtomicResult(**output)
