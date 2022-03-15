"""Compute dispersion correction using Greenwell & Beran's MP2D executable."""

import pprint
import re
import sys
from decimal import Decimal
from typing import Any, Dict, Optional, Tuple

import numpy as np
import qcelemental as qcel
from qcelemental.models import AtomicResult, Provenance
from qcelemental.util import safe_version, which

from ..exceptions import InputError, ResourceError, UnknownError
from ..util import execute
from . import empirical_dispersion_resources
from .model import ProgramHarness

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)


class MP2DHarness(ProgramHarness):

    _defaults = {
        "name": "MP2D",
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
            "mp2d",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install mp2d -c psi4`",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which("mp2d")
        if which_prog not in self.version_cache:
            # Note: anything below v1.1 will return an input error message here. but that's fine as version compare evals to False.
            command = [which_prog, "--version"]
            import subprocess

            proc = subprocess.run(command, stdout=subprocess.PIPE)
            self.version_cache[which_prog] = safe_version(proc.stdout.decode("utf-8").strip())

        return self.version_cache[which_prog]

    def compute(self, input_model: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        from ..testing import is_program_new_enough

        self.found(raise_error=True)

        if not is_program_new_enough("mp2d", "1.1"):
            raise ResourceError(f"MP2D version '{self.get_version()}' too old. Please update to at least '1.1'.")

        job_inputs = self.build_input(input_model, config)

        success, dexe = self.execute(job_inputs)

        if success:
            dexe["outfiles"]["stdout"] = dexe["stdout"]
            dexe["outfiles"]["stderr"] = dexe["stderr"]
            output_model = self.parse_output(dexe["outfiles"], input_model)

        else:
            output_model = input_model
            output_model["error"] = {"error_type": "execution_error", "error_message": dexe["stderr"]}

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
        )
        return success, dexe

    def build_input(
        self, input_model: "AtomicInput", config: "TaskConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:

        # strip engine hint
        mtd = input_model.model.method
        if mtd.startswith("mp2d-"):
            mtd = mtd[5:]

        if input_model.driver.derivative_int() > 1:
            raise InputError(f"Driver {input_model.driver} not implemented for MP2D.")

        # temp until actual options object
        input_model.extras["info"] = empirical_dispersion_resources.from_arrays(
            name_hint=mtd,
            level_hint=input_model.keywords.get("level_hint", None),
            param_tweaks=input_model.keywords.get("params_tweaks", None),
            dashcoeff_supplement=input_model.keywords.get("dashcoeff_supplement", None),
        )

        # Need 'real' field later and that's only guaranteed for molrec
        molrec = qcel.molparse.from_schema(input_model.molecule.dict())
        xyz = qcel.molparse.to_string(molrec, dtype="xyz", units="Angstrom", ghost_format="")
        infiles = {"mp2d_geometry": xyz}
        # jobrec['molecule']['real'] = molrec['real']

        # env = {
        #    'HOME': os.environ.get('HOME'),
        #    'PATH': os.environ.get('PATH'),
        #    #'PATH': os.pathsep.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(os.pathsep) if x != '']) + \
        #    #        os.pathsep + os.environ.get('PATH'),
        #    #'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH'),
        # }

        command = ["mp2d", "mp2d_geometry"]
        command.extend(
            """--TT_a1={a1} --TT_a2={a2} --rcut={rcut} --w={w} --s8={s8}""".format(
                **input_model.extras["info"]["dashparams"]
            ).split()
        )
        if input_model.driver == "gradient":
            command.append("--gradient")

        return {
            "command": command,
            "infiles": infiles,
            "outfiles": ["mp2d_gradient"],
            "scratch_messy": config.scratch_messy,
            "scratch_directory": config.scratch_directory,
            "input_result": input_model.copy(deep=True),
        }

    def parse_output(self, outfiles: Dict[str, str], input_model: "AtomicInput") -> "AtomicResult":
        stdout = outfiles.pop("stdout")

        for fl, contents in outfiles.items():
            if contents is not None:
                # LOG text += f'\n  MP2D scratch file {fl} has been read.\n'
                pass

        # parse energy output (could go further and break into UCHF, CKS)
        real = np.array(input_model.molecule.real)
        full_nat = real.shape[0]
        real_nat = np.sum(real)

        for ln in stdout.splitlines():
            if re.match("   MP2D dispersion correction Eh", ln):
                ene = Decimal(ln.split()[4])
            elif re.match("Atomic Coordinates in Angstroms", ln):
                break
        else:
            if not ((real_nat == 1) and (input_model.driver == "gradient")):
                raise UnknownError("Unknown issue occured.")

        # parse gradient output
        if outfiles["mp2d_gradient"] is not None:
            srealgrad = outfiles["mp2d_gradient"]
            realgrad = np.fromstring(srealgrad, count=3 * real_nat, sep=" ").reshape((-1, 3))

        if input_model.driver == "gradient":
            ireal = np.argwhere(real).reshape((-1))
            fullgrad = np.zeros((full_nat, 3))
            try:
                fullgrad[ireal, :] = realgrad
            except NameError as exc:
                raise UnknownError("Unsuccessful gradient collection.") from exc

        qcvkey = input_model.extras["info"]["fctldash"].upper()

        calcinfo = []
        calcinfo.append(qcel.Datum("CURRENT ENERGY", "Eh", ene))
        calcinfo.append(qcel.Datum("DISPERSION CORRECTION ENERGY", "Eh", ene))
        calcinfo.append(qcel.Datum("2-BODY DISPERSION CORRECTION ENERGY", "Eh", ene))
        if qcvkey:
            calcinfo.append(qcel.Datum(f"{qcvkey} DISPERSION CORRECTION ENERGY", "Eh", ene))

        if input_model.driver == "gradient":
            calcinfo.append(qcel.Datum("CURRENT GRADIENT", "Eh/a0", fullgrad))
            calcinfo.append(qcel.Datum("DISPERSION CORRECTION GRADIENT", "Eh/a0", fullgrad))
            calcinfo.append(qcel.Datum("2-BODY DISPERSION CORRECTION GRADIENT", "Eh/a0", fullgrad))
            if qcvkey:
                calcinfo.append(qcel.Datum(f"{qcvkey} DISPERSION CORRECTION GRADIENT", "Eh/a0", fullgrad))

        # LOGtext += qcel.datum.print_variables({info.label: info for info in calcinfo})
        calcinfo = {info.label: info.data for info in calcinfo}
        # calcinfo = qcel.util.unnp(calcinfo, flat=True)

        # got to even out who needs plump/flat/Decimal/float/ndarray/list
        # Decimal --> str preserves precision
        calcinfo = {k.upper(): str(v) if isinstance(v, Decimal) else v for k, v in calcinfo.items()}

        # jobrec['properties'] = {"return_energy": ene}
        # jobrec["molecule"]["real"] = list(jobrec["molecule"]["real"])

        retres = calcinfo[f"CURRENT {input_model.driver.upper()}"]
        if isinstance(retres, Decimal):
            retres = float(retres)
        elif isinstance(retres, np.ndarray):
            retres = retres.ravel().tolist()

        output_data = {
            "extras": input_model.extras,
            "properties": {},
            "provenance": Provenance(
                creator="MP2D", version=self.get_version(), routine=__name__ + "." + sys._getframe().f_code.co_name
            ),
            "return_result": retres,
            "stdout": stdout,
        }
        output_data["extras"]["local_keywords"] = input_model.extras["info"]
        output_data["extras"]["qcvars"] = calcinfo

        output_data["success"] = True
        return AtomicResult(**{**input_model.dict(), **output_data})
