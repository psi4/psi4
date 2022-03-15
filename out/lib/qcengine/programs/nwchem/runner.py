"""
Calls the NWChem executable.
"""
import copy
import logging
import pprint
import re
import sys
import traceback
from decimal import Decimal
from typing import Any, Dict, Optional, Tuple

import numpy as np
from qcelemental.models import AtomicInput, AtomicResult, BasisSet, Provenance
from qcelemental.util import safe_version, which, which_import

from qcengine.config import TaskConfig, get_config
from qcengine.exceptions import UnknownError

from ...exceptions import InputError
from ...util import create_mpi_invocation, execute
from ..model import ErrorCorrectionProgramHarness
from ..qcvar_identities_resources import build_atomicproperties, build_out
from ..util import error_stamp
from .errors import all_errors
from .germinate import muster_modelchem
from .harvester import harvest
from .keywords import format_keywords

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)
logger = logging.getLogger(__name__)


class NWChemHarness(ErrorCorrectionProgramHarness):
    """

    Notes
    -----
    * To use the TCE, specify ``AtomicInput.model.method`` as usual, then also include ``qc_module = True`` in ``AtomicInput.keywords``.

    """

    _defaults = {
        "name": "NWChem",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": False,
        "node_parallel": True,
        "managed_memory": True,
    }
    # ATL: OpenMP only >=6.6 and only for Phi; potential for Mac using MKL and Intel compilers
    version_cache: Dict[str, str] = {}

    class Config(ErrorCorrectionProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        """Whether NWChem harness is ready for operation, with both the QC program and any particular dependencies found.

        Parameters
        ----------
        raise_error: bool
            Passed on to control negative return between False and ModuleNotFoundError raised.

        Returns
        -------
        bool
            If both nwchem and its harness dependency networkx are found, returns True.
            If raise_error is False and nwchem or networkx are missing, returns False.
            If raise_error is True and nwchem or networkx are missing, the error message for the first missing one is raised.

        """
        qc = which(
            "nwchem",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via http://www.nwchem-sw.org/index.php/Download",
        )

        dep = which_import(
            "networkx",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="For NWChem harness, please install via `conda install networkx -c conda-forge`.",
        )

        return qc and dep

    def get_version(self) -> str:
        self.found(raise_error=True)

        # Get the node configuration
        config = get_config()

        # Run NWChem
        which_prog = which("nwchem")
        if config.use_mpiexec:
            command = create_mpi_invocation(which_prog, config)
        else:
            command = [which_prog]
        command.append("v.nw")

        if which_prog not in self.version_cache:
            success, output = execute(command, {"v.nw": ""}, scratch_directory=config.scratch_directory)

            if success:
                for line in output["stdout"].splitlines():
                    if "nwchem branch" in line:
                        branch = line.strip().split()[-1]
                    if "nwchem revision" in line:
                        revision = line.strip().split()[-1]
                self.version_cache[which_prog] = safe_version(branch + "+" + revision)
            else:
                raise UnknownError(output["stderr"])

        return self.version_cache[which_prog]

    def _compute(self, input_model: AtomicInput, config: "TaskConfig") -> AtomicResult:
        """
        Runs NWChem in executable mode
        """
        self.found(raise_error=True)

        job_inputs = self.build_input(input_model, config)
        success, dexe = self.execute(job_inputs)

        if "There is an error in the input file" in dexe["stdout"]:
            raise InputError(dexe["stdout"])
        if "not compiled" in dexe["stdout"]:
            # recoverable with a different compilation with optional modules
            raise InputError(dexe["stdout"])

        if success:
            dexe["outfiles"]["stdout"] = dexe["stdout"]
            dexe["outfiles"]["stderr"] = dexe["stderr"]
            return self.parse_output(dexe["outfiles"], input_model)
        else:
            # Check if any of the errors are known
            for error in all_errors:
                error.detect_error(dexe)
            raise UnknownError(f"STDOUT:\n{dexe['stdout']}\nSTDERR:\n{dexe['stderr']}")

    def build_input(
        self, input_model: AtomicInput, config: TaskConfig, template: Optional[str] = None
    ) -> Dict[str, Any]:
        nwchemrec = {
            "infiles": {},
            "scratch_directory": config.scratch_directory,
            "scratch_messy": config.scratch_messy,
        }

        opts = copy.deepcopy(input_model.keywords)
        opts = {k.lower(): v for k, v in opts.items()}

        # Handle memory
        # * [GiB] --> [QW]
        # * int() rounds down
        # * was [GiB] --> [B] c. v6.6 but fails in v7.0 probably b/c https://github.com/nwchemgit/nwchem/commit/fca382eab477c3e85548457bfceb1fc9be31b47c#diff-7baaf4807cc9b853af14d9127f63db47d706e12f697a98560bc98bb647ef8326
        #   * memory_size = int(config.memory * (1024 ** 3))
        memory_size = int(config.memory * (1024 ** 3) / 8)
        if config.use_mpiexec:  # It is the memory per MPI rank
            memory_size //= config.nnodes * config.ncores // config.cores_per_rank
        opts["memory"] = memory_size

        # Handle molecule
        molcmd, moldata = input_model.molecule.to_string(dtype="nwchem", units="Bohr", return_data=True)
        opts.update(moldata["keywords"])

        if opts.pop("geometry__noautoz", False):
            molcmd = re.sub(r"geometry ([^\n]*)", r"geometry \1 noautoz", molcmd)

        # Handle calc type and quantum chemical method
        mdccmd, mdcopts = muster_modelchem(input_model.model.method, input_model.driver, opts.pop("qc_module", False))
        opts.update(mdcopts)

        # Handle basis set
        if isinstance(input_model.model.basis, BasisSet):
            raise InputError("QCSchema BasisSet for model.basis not implemented. Use string basis name.")

        # * for nwchem, still needs sph and ghost
        for el in set(input_model.molecule.symbols):
            opts[f"basis__{el}"] = f"library {input_model.model.basis}"
            opts[f"basis__bq{el}"] = f"library {el} {input_model.model.basis}"

        # Log the job settings
        logger.debug("JOB_OPTS")
        logger.debug(pp.pformat(opts))

        # Handle conversion from schema (flat key/value) keywords into local format
        optcmd = format_keywords(opts)

        # Combine the molecule description, options and method command together
        nwchemrec["infiles"]["nwchem.nw"] = "echo\n" + molcmd + optcmd + mdccmd

        # For gradient methods, add a Python command to save the gradients in higher precision
        #  Note: The Hessian is already stored in high precision in a file named "*.hess"
        if input_model.driver == "gradient":
            # Get the name of the theory used for computing the gradients
            theory = re.search(r"^task (\w+) ", mdccmd, re.MULTILINE).group(1)
            logger.debug(f"Adding a Python task to retrieve gradients. Theory: {theory}")

            # Create a Python function to get the gradient from NWChem's checkpoint file (rtdb)
            #  and save it to a JSON file. The stdout for NWChem only prints 6 _decimal places_
            #  (not 6 significant figures)
            pycmd = f"""
python
  grad = rtdb_get('{theory}:gradient')
  if ga_nodeid() == 0:
      import json
      with open('nwchem.grad', 'w') as fp:
          json.dump(grad, fp)
end

task python
            """
            nwchemrec["infiles"]["nwchem.nw"] += pycmd

        # Determine the command
        if config.use_mpiexec:
            nwchemrec["command"] = create_mpi_invocation(which("nwchem"), config)
            logger.info(f"Launching with mpiexec: {' '.join(nwchemrec['command'])}")
        else:
            nwchemrec["command"] = [which("nwchem")]

        return nwchemrec

    def execute(
        self, inputs: Dict[str, Any], *, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None
    ) -> Tuple[bool, Dict]:

        success, dexe = execute(
            inputs["command"],
            inputs["infiles"],
            ["nwchem.hess", "nwchem.grad"],
            scratch_messy=inputs["scratch_messy"],
            scratch_exist_ok=True,
            scratch_directory=inputs["scratch_directory"],
        )
        return success, dexe

    def parse_output(
        self, outfiles: Dict[str, str], input_model: "AtomicInput"
    ) -> AtomicResult:  # lgtm: [py/similar-function]

        # Get the stdout from the calculation (required)
        stdout = outfiles.pop("stdout")
        stderr = outfiles.pop("stderr")

        method = input_model.model.method.lower()
        method = method[4:] if method.startswith("nwc-") else method

        # Read the NWChem stdout file and, if needed, the hess or grad files
        # July 2021: nwmol & vector returns now atin/outfile orientation depending on fix_com,orientation=T/F. previously always atin orientation
        try:
            qcvars, nwhess, nwgrad, nwmol, version, module, errorTMP = harvest(
                input_model.molecule, method, stdout, **outfiles
            )
        except Exception as e:
            raise UnknownError(error_stamp(stdout=stdout, stderr=stderr))

        try:
            if nwgrad is not None:
                qcvars[f"{method.upper()} TOTAL GRADIENT"] = nwgrad
                qcvars["CURRENT GRADIENT"] = nwgrad

            if nwhess is not None:
                qcvars[f"{method.upper()} TOTAL HESSIAN"] = nwhess
                qcvars["CURRENT HESSIAN"] = nwhess

            # Normalize the output as a float or list of floats
            if input_model.driver.upper() == "PROPERTIES":
                retres = qcvars[f"CURRENT ENERGY"]
            else:
                retres = qcvars[f"CURRENT {input_model.driver.upper()}"]
        except KeyError as e:
            raise UnknownError(error_stamp(stdout=stdout, stderr=stderr))

        if isinstance(retres, Decimal):
            retres = float(retres)
        elif isinstance(retres, np.ndarray):
            retres = retres.tolist()

        # Get the formatted properties
        build_out(qcvars)
        atprop = build_atomicproperties(qcvars)

        provenance = Provenance(creator="NWChem", version=self.get_version(), routine="nwchem").dict()
        if module is not None:
            provenance["module"] = module

        # Format them inout an output
        output_data = {
            "schema_version": 1,
            "molecule": nwmol,  # overwrites with outfile Cartesians in case fix_*=F
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
