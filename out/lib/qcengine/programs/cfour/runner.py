"""Compute quantum chemistry using Mainz-Austin-Budapest-Gainesville's CFOUR executable."""

import copy
import pprint
import sys
import traceback
from decimal import Decimal
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np
from qcelemental.models import AtomicInput, AtomicResult, BasisSet, Provenance
from qcelemental.util import safe_version, which

from ...exceptions import InputError, UnknownError
from ...util import execute
from ..model import ProgramHarness
from ..qcvar_identities_resources import build_atomicproperties, build_out
from .germinate import muster_modelchem
from .harvester import harvest
from .keywords import format_keywords

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)


class CFOURHarness(ProgramHarness):
    """

    Notes
    -----
    * Looks for basis set file ``../basis/GENBAS`` from ``xcfour`` executable. If this doesn't work, file an issue.

    """

    _defaults = {
        "name": "CFOUR",
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
            "xcfour", return_bool=True, raise_error=raise_error, raise_msg="Please install via http://cfour.de/"
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which("xcfour")
        if which_prog not in self.version_cache:
            success, output = execute([which_prog, "ZMAT"], {"ZMAT": "\nHe\n\n"})

            if success:
                for line in output["stdout"].splitlines():
                    if "Version" in line:
                        branch = " ".join(line.strip().split()[1:])
                self.version_cache[which_prog] = safe_version(branch)

        return self.version_cache[which_prog]

    def compute(self, input_model: AtomicInput, config: "TaskConfig") -> AtomicResult:
        self.found(raise_error=True)

        job_inputs = self.build_input(input_model, config)
        success, dexe = self.execute(job_inputs)

        if success:
            dexe["outfiles"]["stdout"] = dexe["stdout"]
            dexe["outfiles"]["stderr"] = dexe["stderr"]
            return self.parse_output(dexe["outfiles"], input_model)

    def build_input(
        self, input_model: AtomicInput, config: "TaskConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:
        cfourrec = {
            "infiles": {},
            "scratch_directory": config.scratch_directory,
            "scratch_messy": config.scratch_messy,
        }

        opts = copy.deepcopy(input_model.keywords)

        # Handle memory
        # for cfour, [GiB] --> [QW]
        opts["memory_size"] = int(config.memory * (1024 ** 3) / 8)
        opts["mem_unit"] = "integerwords"

        # Handle molecule
        molcmd, moldata = input_model.molecule.to_string(dtype="cfour", units="Bohr", return_data=True)
        opts.update(moldata["keywords"])

        # Handle calc type and quantum chemical method
        mdcopts = muster_modelchem(input_model.model.method, input_model.driver.derivative_int())
        opts.update(mdcopts)

        # Handle basis set
        if isinstance(input_model.model.basis, BasisSet):
            raise InputError("QCSchema BasisSet for model.basis not implemented. Use string basis name.")
        if input_model.model.basis is None:
            raise InputError("None for model.basis is not useable.")

        # * why, yes, this is highly questionable
        #   * assuming relative file location between xcfour exe and GENBAS file
        #   * reading a multi MB file into the inputs dict
        if all(input_model.molecule.real):
            opts["basis"] = input_model.model.basis
            bascmd = ""
        else:
            # * note not getting per-basis casing like if it passed through format_keywords
            opts["basis"] = "SPECIAL"
            text = [
                (
                    f"""H:6-31G"""
                    if (elem == "H" and input_model.model.basis.upper() == "6-31G*")
                    else f"""{elem.upper()}:{input_model.model.basis.upper()}"""
                )
                for iat, elem in enumerate(input_model.molecule.symbols)
            ]
            text.append("")
            text.append("")
            bascmd = "\n".join(text)

        # Handle conversion from schema (flat key/value) keywords into local format
        optcmd = format_keywords(opts)

        xcfour = which("xcfour")
        genbas = Path(xcfour).parent.parent / "basis" / "GENBAS"
        cfourrec["infiles"]["ZMAT"] = molcmd + optcmd + bascmd
        cfourrec["infiles"]["GENBAS"] = genbas.read_text()
        cfourrec["command"] = [xcfour]

        return cfourrec

    def execute(
        self, inputs: Dict[str, Any], *, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None
    ) -> Tuple[bool, Dict]:

        # llel works b/c util.environ_context sets OMP_NUM_THREADS = config.ncores

        success, dexe = execute(
            inputs["command"],
            inputs["infiles"],
            ["GRD", "FCMFINAL", "DIPOL"],
            # "DIPDER", "POLAR", "POLDER"],
            scratch_messy=inputs["scratch_messy"],
            scratch_directory=inputs["scratch_directory"],
        )
        return success, dexe

    def parse_output(
        self, outfiles: Dict[str, str], input_model: AtomicInput
    ) -> AtomicResult:  # lgtm: [py/similar-function]

        stdout = outfiles.pop("stdout")
        stderr = outfiles.pop("stderr")

        method = input_model.model.method.lower()
        method = method[3:] if method.startswith("c4-") else method

        # c4mol, if it exists, is dinky, just a clue to geometry of cfour results
        try:
            # July 2021: c4mol & vector returns now atin/outfile orientation depending on fix_com,orientation=T/F. previously always atin orientation
            qcvars, c4hess, c4grad, c4mol, version, module, errorTMP = harvest(
                input_model.molecule, method, stdout, **outfiles
            )
        except Exception as e:
            raise UnknownError(
                "STDOUT:\n"
                + stdout
                + "\nSTDERR:\n"
                + stderr
                + "\nTRACEBACK:\n"
                + "".join(traceback.format_exception(*sys.exc_info()))
            )

        if errorTMP != "":
            raise UnknownError("STDOUT:\n" + stdout + "\nSTDERR:\n" + stderr)

        try:
            if c4grad is not None:
                qcvars["CURRENT GRADIENT"] = c4grad
                qcvars[f"{method.upper()} TOTAL GRADIENT"] = c4grad

            if c4hess is not None:
                qcvars[f"{method.upper()} TOTAL HESSIAN"] = c4hess
                qcvars["CURRENT HESSIAN"] = c4hess

            if input_model.driver.upper() == "PROPERTIES":
                retres = qcvars[f"CURRENT ENERGY"]
            else:
                retres = qcvars[f"CURRENT {input_model.driver.upper()}"]
        except KeyError as e:
            raise UnknownError(
                "STDOUT:\n"
                + stdout
                + "\nSTDERR:\n"
                + stderr
                + "\nTRACEBACK:\n"
                + "".join(traceback.format_exception(*sys.exc_info()))
            )

        # TODO: "xalloc(): memory allocation failed!"

        if isinstance(retres, Decimal):
            retres = float(retres)
        elif isinstance(retres, np.ndarray):
            retres = retres.ravel().tolist()

        build_out(qcvars)
        atprop = build_atomicproperties(qcvars)

        provenance = Provenance(creator="CFOUR", version=self.get_version(), routine="xcfour").dict()
        if module is not None:
            provenance["module"] = module

        output_data = {
            "schema_version": 1,
            "molecule": c4mol,  # overwrites with outfile Cartesians in case fix_*=F
            "extras": {"outfiles": outfiles, **input_model.extras},
            "properties": atprop,
            "provenance": provenance,
            "return_result": retres,
            "stderr": stderr,
            "stdout": stdout,
            "success": True,
        }

        # got to even out who needs plump/flat/Decimal/float/ndarray/list
        # Decimal --> str preserves precision
        # * formerly unnp(qcvars, flat=True).items()
        output_data["extras"]["qcvars"] = {
            k.upper(): str(v) if isinstance(v, Decimal) else v for k, v in qcvars.items()
        }

        return AtomicResult(**{**input_model.dict(), **output_data})
