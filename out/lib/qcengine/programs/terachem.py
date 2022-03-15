"""
Calls the TeraChem executable.
"""

import re
from typing import Any, Dict, Optional

from qcelemental.models import AtomicResult, FailedOperation
from qcelemental.molparse.regex import DECIMAL, NUMBER
from qcelemental.util import parse_version, safe_version, which

from qcengine import util as uti

from ..exceptions import UnknownError
from ..util import popen
from .model import ProgramHarness


class TeraChemHarness(ProgramHarness):

    _defaults = {
        "name": "TeraChem",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which(
            "terachem",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via http://www.petachem.com/index.html",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which("terachem")
        if which_prog not in self.version_cache:
            with popen([which_prog, "--version"]) as exc:
                exc["proc"].wait(timeout=5)
            mobj = re.search(NUMBER, exc["stdout"], re.VERBOSE)
            version = mobj.group(0)
            self.version_cache[which_prog] = safe_version(version)

        return self.version_cache[which_prog]

    def compute(self, input_data: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        """
        Run TeraChem
        """
        self.found(raise_error=True)

        # Check TeraChem version
        if parse_version(self.get_version()) < parse_version("1.5"):
            raise TypeError("TeraChem version '{}' not understood".format(self.get_version()))

        # Setup the job
        job_inputs = self.build_input(input_data, config)
        # Run terachem
        exe_outputs = self.execute(job_inputs, extra_outfiles=input_data.extras)
        exe_success, proc = exe_outputs
        # Determine whether the calculation succeeded
        output_data = {}
        if not exe_success:
            output_data["success"] = False
            output_data["error"] = {"error_type": "unknown_error", "error_message": proc["stderr"]}
            return FailedOperation(
                success=output_data.pop("success", False), error=output_data.pop("error"), input_data=output_data
            )

        # If execution succeeded, collect results
        result = self.parse_output(proc["outfiles"], input_data)
        return result

    def build_input(
        self, input_model: "AtomicInput", config: "TaskConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:
        # Write the geom xyz file with unit au
        xyz_file = input_model.molecule.to_string(dtype="terachem", units="Bohr")

        # Write input file
        input_file = []
        input_file.append("# molecule definition")
        input_file.append("units bohr")
        input_file.append("charge " + str(int(input_model.molecule.molecular_charge)))
        input_file.append("spinmult " + str(input_model.molecule.molecular_multiplicity))
        input_file.append("coordinates geometry.xyz")

        input_file.append("\n# model")
        input_file.append("basis " + str(input_model.model.basis))
        input_file.append("method " + str(input_model.model.method))

        input_file.append("\n# driver")
        input_file.append("run " + input_model.driver)

        input_file.append("\n# keywords")
        for k, v in input_model.keywords.items():
            input_file.append("{} {}".format(k, v))

        input_file = "\n".join(input_file)

        return {
            "commands": ["terachem", "tc.in"],
            "infiles": {"tc.in": input_file, "geometry.xyz": xyz_file},
            "scratch_directory": config.scratch_directory,
            "input_result": input_model.copy(deep=True),
        }

    def parse_output(self, outfiles: Dict[str, str], input_model: "AtomicInput") -> "AtomicResult":
        output_data = {}
        properties = {}

        # Parse the output file, collect properties and gradient
        output_lines = outfiles["tc.out"].split("\n")
        gradients = []
        natom = 0
        line_final_energy = -1
        line_scf_header = -1
        for idx, line in enumerate(output_lines):
            if "FINAL ENERGY" in line:
                properties["scf_total_energy"] = float(line.strip("\n").split()[2])
                line_final_energy = idx
            elif "Start SCF Iterations" in line:
                line_scf_header = idx
            elif "Total atoms" in line:
                natom = int(line.split()[-1])
            elif "DIPOLE MOMENT" in line:
                newline = line.replace(",", "").replace("}", "").replace("{", "")
                properties["scf_dipole_moment"] = [float(x) for x in newline.split()[2:5]]
            elif "Nuclear repulsion energy" in line:
                properties["nuclear_repulsion_energy"] = float(line.split()[-2])
            elif "Gradient units are Hartree/Bohr" in line:
                # Gradient is stored as (dE/dx1,dE/dy1,dE/dz1,dE/dx2,dE/dy2,...)
                for i in range(idx + 3, idx + 3 + natom):
                    grad = output_lines[i].strip("\n").split()
                    for x in grad:
                        gradients.append(float(x))

        last_scf_line = ""
        for idx in reversed(range(line_scf_header, line_final_energy)):
            mobj = re.search(
                r"^\s*\d+\s+" + DECIMAL + r"\s+" + DECIMAL + r"\s+" + DECIMAL + r"\s+" + DECIMAL,
                output_lines[idx],
                re.VERBOSE,
            )
            if mobj:
                last_scf_line = output_lines[idx]
                break

        if len(last_scf_line) > 0:
            properties["scf_iterations"] = int(last_scf_line.split()[0])
            if "XC Energy" in output_lines:
                properties["scf_xc_energy"] = float(last_scf_line.split()[4])
        else:
            raise UnknownError("SCF iteration lines not found in TeraChem output")

        if len(gradients) > 0:
            output_data["return_result"] = gradients

        # Commented out the properties currently not supported by QCSchema
        # properites["spin_S2"] = 1 # calculated S(S+1)
        #   elif "SPIN S-SQUARED" in line:
        #       properties["spin_S2"] = float(line.strip('\n').split()[2])
        # Parse files in scratch folder
        # properties["atomic_charge"] = []
        # atomic_charge_lines =  open(outfiles["charge.xls"]).readlines()
        # for line in atomic_charge_lines:
        #    properties["atomic_charge"].append(line.strip('\n').split()[-1])

        if "return_result" not in output_data:
            if "scf_total_energy" in properties:
                output_data["return_result"] = properties["scf_total_energy"]
            else:
                raise KeyError("Could not find SCF total energy")

        output_data["properties"] = properties

        output_data["schema_name"] = "qcschema_output"
        output_data["stdout"] = outfiles["tc.out"]
        # TODO Should only return True if TeraChem calculation terminated properly
        output_data["success"] = True

        # return extra files requested by user as extras
        for extra in input_model.extras.keys():
            input_model.extras[extra] = outfiles[extra]

        return AtomicResult(**{**input_model.dict(), **output_data})

    def execute(self, inputs, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None):
        binaries = []
        for filename in extra_outfiles:
            if filename in ["scr/ca0", "scr/cb0", "scr/c0"]:
                binaries.append(filename)
        exe_success, proc = uti.execute(
            inputs["commands"],
            infiles=inputs["infiles"],
            outfiles=extra_outfiles,
            as_binary=binaries,
            scratch_directory=inputs["scratch_directory"],
            timeout=timeout,
        )
        proc["outfiles"]["tc.out"] = proc["stdout"]
        return exe_success, proc
