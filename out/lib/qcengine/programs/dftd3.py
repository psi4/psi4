"""Compute dispersion correction using Grimme's DFTD3 executable."""

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

from ..exceptions import InputError, ResourceError, UnknownError
from ..util import execute
from . import empirical_dispersion_resources
from .model import ProgramHarness

if TYPE_CHECKING:
    from qcelemental.models import AtomicInput

    from ..config import TaskConfig


pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)


class DFTD3Harness(ProgramHarness):

    _defaults = {
        "name": "DFTD3",
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
            "dftd3",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install dftd3 -c psi4`.",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which("dftd3")
        if which_prog not in self.version_cache:
            # Note: anything below v3.2.1 will return the help menu here. but that's fine as version compare evals to False.
            command = [which_prog, "-version"]
            import subprocess

            proc = subprocess.run(command, stdout=subprocess.PIPE)
            self.version_cache[which_prog] = safe_version(proc.stdout.decode("utf-8").strip())

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
            # env=inputs["env"],
            scratch_directory=inputs["scratch_directory"],
            blocking_files=inputs["blocking_files"],
        )
        return success, dexe

    def build_input(
        self, input_model: "AtomicInput", config: "TaskConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:

        # strip engine hint
        mtd = input_model.model.method
        if mtd.startswith("d3-"):
            mtd = mtd[3:]

        if (input_model.driver.derivative_int() > 1) or (input_model.driver == "properties"):
            raise InputError(f"Driver {input_model.driver} not implemented for DFTD3.")

        # temp until actual options object
        input_model.extras["info"] = empirical_dispersion_resources.from_arrays(
            name_hint=mtd,
            level_hint=input_model.keywords.get("level_hint", None),
            param_tweaks=input_model.keywords.get("params_tweaks", None),
            dashcoeff_supplement=input_model.keywords.get("dashcoeff_supplement", None),
        )

        # this is what the dftd3 program needs, not what the job needs
        # * form dftd3_parameters string that governs dispersion calc
        # * form dftd3_geometry string that supplies geometry to dispersion calc
        # * form command and arguments

        # Need 'real' field later and that's only guaranteed for molrec
        molrec = qcel.molparse.from_schema(input_model.molecule.dict())
        # jobrec['molecule']['real'] = molrec['real']

        command = ["dftd3", "dftd3_geometry.xyz"]
        if input_model.driver == "gradient":
            command.append("-grad")
        if input_model.extras["info"]["dashlevel"] == "atmgr":
            command.append("-abc")

        # Append `-anal` for pairwise atomic analysis
        if input_model.keywords.get("pair_resolved", False):
            command.append("-anal")

        infiles = {
            ".dftd3par.local": dftd3_coeff_formatter(
                input_model.extras["info"]["dashlevel"], input_model.extras["info"]["dashparams"]
            ),
            "dftd3_geometry.xyz": qcel.molparse.to_string(molrec, dtype="xyz", units="Angstrom", ghost_format=""),
        }

        return {
            "command": command,
            "infiles": infiles,
            "outfiles": ["dftd3_gradient", "dftd3_abc_gradient"],
            "scratch_messy": config.scratch_messy,
            "scratch_directory": config.scratch_directory,
            "input_result": input_model.copy(deep=True),
            "blocking_files": [os.path.join(pathlib.Path.home(), ".dftd3par." + socket.gethostname())],
        }

    #   Notes
    #   -----
    #   Central to harvesting is the fact (to the planting, not to the DFTD3
    #   program) that 2-body and 3-body are run separately. Because of how
    #   damping functions work (see GH:psi4/psi4#1407), some 2-body damping
    #   schemes can give wrong answers for 3-body. And because 3-body is
    #   set to run with some dummy values, the 2-body values are no good.

    def parse_output(self, outfiles: Dict[str, str], input_model: "AtomicInput") -> "AtomicResult":
        Grimme_h2kcal = 627.509541
        stdout = outfiles.pop("stdout")

        for fl, contents in outfiles.items():
            if contents is not None:
                # LOG text += f'\n  DFTD3 scratch file {fl} has been read.\n'
                pass

        # parse energy output (could go further and break into E6, E8, E10 and Cn coeff)
        real = np.array(input_model.molecule.real)
        full_nat = real.shape[0]
        real_nat = np.sum(real)
        for ln in stdout.splitlines():
            if re.match(" Edisp /kcal,au", ln):
                ene = Decimal(ln.split()[3])
            elif re.match(r" E6\(ABC\) \"   :", ln):  # c. v3.2.0
                raise ResourceError("Cannot process ATM results from DFTD3 prior to v3.2.1.")
            elif re.match(r""" E6\(ABC\) /kcal,au:""", ln):
                atm = Decimal(ln.split()[-1])
            elif re.match(" analysis of pair-wise terms", ln):
                D3pairs = np.zeros((full_nat, full_nat))
                # Iterate over block
                start = stdout.splitlines().index(ln) + 2
                for l in stdout.splitlines()[start:]:
                    data = l.replace("-", " -").split()
                    # print(data)
                    if len(data) == 0:
                        break
                    atom1 = int(data[0]) - 1
                    atom2 = int(data[1]) - 1
                    Edisp = Decimal(data[-1])
                    D3pairs[atom1, atom2] = Edisp / Decimal(Grimme_h2kcal)
                    D3pairs[atom2, atom1] = D3pairs[atom1, atom2]

            elif re.match(" normal termination of dftd3", ln):
                break
        else:
            if not ((real_nat == 1) and (input_model.driver == "gradient")):
                raise UnknownError(
                    f"Unsuccessful run. Check input, particularly geometry in [a0]. Model: {input_model.model}"
                )

        # parse gradient output
        # * DFTD3 crashes on one-atom gradients. Avoid the error (above) and just force the correct result (below).
        if outfiles["dftd3_gradient"] is not None:
            srealgrad = outfiles["dftd3_gradient"].replace("D", "E")
            realgrad = np.fromstring(srealgrad, count=3 * real_nat, sep=" ").reshape((-1, 3))
        elif real_nat == 1:
            realgrad = np.zeros((1, 3))

        if outfiles["dftd3_abc_gradient"] is not None:
            srealgrad = outfiles["dftd3_abc_gradient"].replace("D", "E")
            realgradabc = np.fromstring(srealgrad, count=3 * real_nat, sep=" ").reshape((-1, 3))
        elif real_nat == 1:
            realgradabc = np.zeros((1, 3))

        if input_model.driver == "gradient":
            ireal = np.argwhere(real).reshape((-1))
            fullgrad = np.zeros((full_nat, 3))
            rg = realgradabc if (input_model.extras["info"]["dashlevel"] == "atmgr") else realgrad
            try:
                fullgrad[ireal, :] = rg
            except NameError as exc:
                raise UnknownError("Unsuccessful gradient collection.") from exc

        qcvkey = input_model.extras["info"]["fctldash"].upper()

        calcinfo = []
        if input_model.extras["info"]["dashlevel"] == "atmgr":
            calcinfo.append(qcel.Datum("CURRENT ENERGY", "Eh", atm))
            calcinfo.append(qcel.Datum("DISPERSION CORRECTION ENERGY", "Eh", atm))
            calcinfo.append(qcel.Datum("3-BODY DISPERSION CORRECTION ENERGY", "Eh", atm))
            calcinfo.append(qcel.Datum("AXILROD-TELLER-MUTO 3-BODY DISPERSION CORRECTION ENERGY", "Eh", atm))

            if input_model.driver == "gradient":
                calcinfo.append(qcel.Datum("CURRENT GRADIENT", "Eh/a0", fullgrad))
                calcinfo.append(qcel.Datum("DISPERSION CORRECTION GRADIENT", "Eh/a0", fullgrad))
                calcinfo.append(qcel.Datum("3-BODY DISPERSION CORRECTION GRADIENT", "Eh/a0", fullgrad))
                calcinfo.append(
                    qcel.Datum("AXILROD-TELLER-MUTO 3-BODY DISPERSION CORRECTION GRADIENT", "Eh/a0", fullgrad)
                )

        else:
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
            "properties": {
                "return_energy": calcinfo[f"CURRENT ENERGY"],
            },
            "provenance": Provenance(
                creator="DFTD3", version=self.get_version(), routine=__name__ + "." + sys._getframe().f_code.co_name
            ),
            "return_result": retres,
            "stdout": stdout,
        }
        output_data["extras"]["local_keywords"] = input_model.extras["info"]
        output_data["extras"]["qcvars"] = calcinfo
        if input_model.keywords.get("pair_resolved", False):
            output_data["extras"]["qcvars"]["2-BODY PAIRWISE DISPERSION CORRECTION ANALYSIS"] = D3pairs
        output_data["success"] = True

        return AtomicResult(**{**input_model.dict(), **output_data})


def dftd3_coeff_formatter(dashlvl: str, dashcoeff: Dict) -> str:
    """Return strings for DFTD3 program parameter file.

             s6      rs6      s18     rs8     alpha6      version
             ------------------------------------------------------
    d2:      s6      sr6      s8=0.0  a2=None alpha6      version=2
    d3zero:  s6      sr6      s8      a2=sr8  alpha6      version=3
    d3bj:    s6      a1       s8      a2      alpha6=None version=4
    d3mzero: s6      sr6      s8      beta    alpha6=14.0 version=5
    d3mbj:   s6      a1       s8      a2      alpha6=None version=6
    atmgr:   s6=1.0  sr6=None s8=None a2=None alpha6      version=3 (needs -abc, too)

    Parameters
    ----------
    dashlvl : {'d2', 'd3zero', d3bj', 'd3mzero', 'd3mbj', 'atmgr'}
        Level of dispersion correction.
    dashcoeff : dict
        Dictionary fully specifying non-fixed parameters (table above) for `dashlvl` to drive DFTD3.

    Notes
    -----
    The `atmgr` dashlvl is intended for use only to get the three-body Axilrod-Teller-Muto
    three body dispersion correction. Therefore, dummy parameters are passed for two-body damping
    function, and it will give garbage for two-body component of dispersion correction.

    Returns
    -------
    str
        Suitable for `.dftd3par` file.

    """
    dashformatter = """{:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:6}\n"""

    dashlvl = dashlvl.lower()
    if dashlvl == "d2":
        return dashformatter.format(dashcoeff["s6"], dashcoeff["sr6"], 0.0, 0.0, dashcoeff["alpha6"], 2)
    elif dashlvl == "d3zero":
        return dashformatter.format(
            dashcoeff["s6"], dashcoeff["sr6"], dashcoeff["s8"], dashcoeff["sr8"], dashcoeff["alpha6"], 3
        )
    elif dashlvl == "d3bj":
        return dashformatter.format(dashcoeff["s6"], dashcoeff["a1"], dashcoeff["s8"], dashcoeff["a2"], 0.0, 4)
    elif dashlvl == "d3mzero":
        return dashformatter.format(dashcoeff["s6"], dashcoeff["sr6"], dashcoeff["s8"], dashcoeff["beta"], 14.0, 5)
    elif dashlvl == "d3mbj":
        return dashformatter.format(dashcoeff["s6"], dashcoeff["a1"], dashcoeff["s8"], dashcoeff["a2"], 0.0, 6)
    elif dashlvl == "atmgr":
        # need to set first four parameters to something other than None, otherwise DFTD3 gets mad or a bit wrong
        return dashformatter.format(1.0, 0.0, 0.0, 0.0, dashcoeff["alpha6"], 3)
    else:
        raise InputError(f"""-D correction level {dashlvl} is not available. Choose among {dashcoeff.keys()}.""")
