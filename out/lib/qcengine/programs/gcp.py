"""Compute geometrical counterpoise correction using Kruse & Grimme's GCP executable."""

import os
import pathlib
import pprint
import re
import socket
import sys
from decimal import Decimal
from typing import TYPE_CHECKING, Any, Dict, Optional, Tuple

import numpy as np
import qcelemental as qcel
from qcelemental.models import AtomicResult, FailedOperation, Provenance
from qcelemental.util import safe_version, which

from ..exceptions import InputError, UnknownError
from ..util import execute
from .model import ProgramHarness

if TYPE_CHECKING:
    from qcelemental.models import AtomicInput

    from ..config import TaskConfig


pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)


class GCPHarness(ProgramHarness):

    _defaults = {
        "name": "GCP",
        "scratch": True,
        "thread_safe": True,
        "thread_parallel": False,
        "node_parallel": False,
        "managed_memory": False,
    }
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which(
            "gcp",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install gcp -c psi4`.",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        version = None
        which_prog = which("gcp")
        if which_prog not in self.version_cache:
            # option not (yet) available, instead find in help output
            command = [which_prog, "-version"]
            import subprocess

            proc = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            for ln in proc.stdout.decode("utf-8").splitlines():
                if re.match(".*Version", ln):
                    version = ln.split()[2]

            if version is None:
                raise UnknownError(f"Could not identify GCP version")
            # pending upcoming gcp update
            # self.version_cache[which_prog] = safe_version(
            #     proc.stdout.decode("utf-8").strip())
            self.version_cache[which_prog] = safe_version(version)

        return self.version_cache[which_prog]

    def compute(self, input_model: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        self.found(raise_error=True)

        job_inputs = self.build_input(input_model, config)

        success, dexe = self.execute(job_inputs)

        if success:
            dexe["outfiles"]["stdout"] = dexe["stdout"]
            dexe["outfiles"]["stderr"] = dexe["stderr"]
            output_model = self.parse_output(dexe["outfiles"], input_model)

        else:
            output_model = FailedOperation(
                success=False,
                error={"error_type": "execution_error", "error_message": dexe["stderr"]},
                input_data=input_model.dict(),
            )

        return output_model

    def execute(
        self, inputs: Dict[str, Any], *, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None
    ) -> Tuple[bool, Dict]:

        success, dexe = execute(
            inputs["command"],
            inputs["infiles"],
            inputs["outfiles"],
            scratch_messy=inputs["scratch_messy"],
            scratch_directory=inputs["scratch_directory"],
            blocking_files=inputs["blocking_files"],
        )
        return success, dexe

    def build_input(
        self, input_model: "AtomicInput", config: "TaskConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:

        if (input_model.driver.derivative_int() > 1) or (input_model.driver == "properties"):
            raise InputError(f"Driver {input_model.driver} not implemented for GCP.")

        # live somewhere else?
        available_levels = [
            "HF/MINIS",
            "DFT/MINIS",
            "HF/MINIX",
            "DFT/MINIX",
            "HF/SV",
            "DFT/SV",
            "HF/def2-SV(P)",
            "DFT/def2-SV(P)",
            "HF/def2-SVP",
            "DFT/def2-SVP",
            "HF/DZP",
            "DFT/DZP",
            "HF/def-TZVP",
            "DFT/def-TZVP",
            "HF/def2-TZVP",
            "DFT/def2-TZVP",
            "HF/631Gd",
            "DFT/631Gd",
            "HF/cc-pVDZ",
            "DFT/cc-pVDZ",
            "HF/aug-cc-pVDZ",
            "DFT/aug-cc-pVDZ",
            # DFT only
            "DFT/SVX",  # = def2-SV(P/h,c)
            "DFT/LANL",
            "DFT/pobTZVP",
            # HF only, needs update
            # 'HF/631G',
            # functional specific
            "TPSS/def2-SVP",
            "PW6B95/def2-SVP",
            # 3c specials
            "hf3c",
            "pbeh3c",
            # custom
            # will need gcp code changes if $HOME is to be avoided
            # thus code blocks with FILE below are not used yet.
            # 'file',
        ]
        available_levels = [f.upper() for f in available_levels]
        # temp until actual options object
        method = input_model.model.method.upper()
        if method not in available_levels:
            raise InputError(f"GCP does not have method: {method}")

        # Need 'real' field later and that's only guaranteed for molrec
        molrec = qcel.molparse.from_schema(input_model.molecule.dict())

        executable = self._defaults["name"].lower()
        calldash = {"gcp": "-", "mctc-gcp": "--"}[executable]

        command = [executable, "gcp_geometry.xyz", calldash + "level", method]

        if input_model.driver == "gradient":
            command.append(calldash + "grad")

        infiles = {
            "gcp_geometry.xyz": qcel.molparse.to_string(molrec, dtype="xyz", units="Angstrom", ghost_format=""),
        }
        if method == "FILE":
            infiles[".gcppar"] = input_model.extras["parameters"]

        return {
            "command": command,
            "infiles": infiles,
            "outfiles": ["gcp_gradient"],
            "scratch_messy": config.scratch_messy,
            "scratch_directory": config.scratch_directory,
            "input_result": input_model.copy(deep=True),
            "blocking_files": [os.path.join(pathlib.Path.home(), ".gcppar." + socket.gethostname())],
        }

    def parse_output(self, outfiles: Dict[str, str], input_model: "AtomicInput") -> "AtomicResult":
        stdout = outfiles.pop("stdout")

        # parse energy output (could go further and break into E6, E8, E10 and Cn coeff)
        real = np.array(input_model.molecule.real)
        full_nat = real.shape[0]
        real_nat = np.sum(real)

        for ln in stdout.splitlines():
            if re.match("  Egcp:", ln):
                ene = Decimal(ln.split()[1])
            elif re.match("     normal termination of gCP", ln):
                break
        else:
            if self._defaults["name"] == "GCP" and not ((real_nat == 1) and (input_model.driver == "gradient")):
                raise UnknownError(
                    f"Unsuccessful run. Check input, particularly geometry in [a0]. Model: {input_model.model}"
                )

        # parse gradient output
        if outfiles["gcp_gradient"] is not None:
            srealgrad = outfiles["gcp_gradient"].replace("D", "E")
            realgrad = np.fromstring(srealgrad, count=3 * real_nat, sep=" ").reshape((-1, 3))
        elif real_nat == 1:
            realgrad = np.zeros((1, 3))

        if input_model.driver == "gradient":
            ireal = np.argwhere(real).reshape((-1))
            fullgrad = np.zeros((full_nat, 3))
            try:
                fullgrad[ireal, :] = realgrad
            except NameError as exc:
                raise UnknownError("Unsuccessful gradient collection.") from exc

        qcvkey = input_model.model.method.upper()

        calcinfo = []

        calcinfo.append(qcel.Datum("CURRENT ENERGY", "Eh", ene))
        calcinfo.append(qcel.Datum("GCP CORRECTION ENERGY", "Eh", ene))
        if qcvkey:
            calcinfo.append(qcel.Datum(f"{qcvkey} GCP CORRECTION ENERGY", "Eh", ene))

        if input_model.driver == "gradient":
            calcinfo.append(qcel.Datum("CURRENT GRADIENT", "Eh/a0", fullgrad))
            calcinfo.append(qcel.Datum("GCP CORRECTION GRADIENT", "Eh/a0", fullgrad))
            if qcvkey:
                calcinfo.append(qcel.Datum(f"{qcvkey} GCP CORRECTION GRADIENT", "Eh/a0", fullgrad))

        calcinfo = {info.label: info.data for info in calcinfo}

        # Decimal --> str preserves precision
        calcinfo = {k.upper(): str(v) if isinstance(v, Decimal) else v for k, v in calcinfo.items()}

        retres = calcinfo[f"CURRENT {input_model.driver.upper()}"]
        if isinstance(retres, Decimal):
            retres = float(retres)
        elif isinstance(retres, np.ndarray):
            retres = retres.ravel().tolist()

        output_data = {
            "extras": input_model.extras,
            "properties": {},
            "provenance": Provenance(
                creator="GCP", version=self.get_version(), routine=__name__ + "." + sys._getframe().f_code.co_name
            ),
            "return_result": retres,
            "stdout": stdout,
        }

        output_data["extras"]["qcvars"] = calcinfo

        output_data["success"] = True
        return AtomicResult(**{**input_model.dict(), **output_data})


class MCTCGCPHarness(GCPHarness):

    _defaults = {
        "name": "MCTC-GCP",
        "scratch": True,
        "thread_safe": True,
        "thread_parallel": False,
        "node_parallel": False,
        "managed_memory": False,
    }
    version_cache: Dict[str, str] = {}

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which(
            "mctc-gcp",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install gcp-correction -c conda-forge`.",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which("mctc-gcp")
        if which_prog not in self.version_cache:
            command = [which_prog, "--version"]
            import subprocess

            proc = subprocess.run(command, stdout=subprocess.PIPE)
            self.version_cache[which_prog] = safe_version(proc.stdout.decode("utf-8").strip().split()[-1])

        return self.version_cache[which_prog]
