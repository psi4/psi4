"""
Calls the Psi4 executable.
"""
import json
import os
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Dict

from qcelemental.models import AtomicResult, BasisSet
from qcelemental.util import deserialize, parse_version, safe_version, which, which_import

from ..exceptions import InputError, RandomError, ResourceError, UnknownError
from ..util import execute, popen, temporary_directory
from .model import ProgramHarness

if TYPE_CHECKING:
    from qcelemental.models import AtomicInput

    from ..config import TaskConfig


class Psi4Harness(ProgramHarness):

    _defaults = {
        "name": "Psi4",
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
        """Whether Psi4 harness is ready for operation.

        Parameters
        ----------
        raise_error: bool
            Passed on to control negative return between False and ModuleNotFoundError raised.

        Returns
        -------
        bool
            If psi4 (psithon or psiapi) is found, returns True.
            If raise_error is False and psi4 is missing, returns False.
            If raise_error is True and psi4 is missing, the error message is raised.

        """
        psithon = which("psi4", return_bool=True)
        psiapi = which_import("psi4", return_bool=True)

        if psithon and not psiapi:
            with popen([which("psi4"), "--module"]) as exc:
                exc["proc"].wait(timeout=30)
            if "module does not exist" in exc["stderr"]:
                pass  # --module argument only in Psi4 DDD branch
            else:
                sys.path.append(exc["stdout"].split()[-1])

        if psiapi and not psithon:
            psiimport = str(Path(which_import("psi4")).parent.parent)
            env = os.environ.copy()
            env["PYTHONPATH"] = psiimport
            with popen(["python", "-c", "import psi4; print(psi4.executable[:-5])"], popen_kwargs={"env": env}) as exc:
                exc["proc"].wait(timeout=30)
            os.environ["PATH"] += os.pathsep + exc["stdout"].split()[-1]

        if psithon or psiapi:
            return True

        return which(
            "psi4",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install psi4 -c psi4`. Check it's in your PATH with `which psi4`.",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which("psi4")
        if which_prog not in self.version_cache:
            with popen([which_prog, "--version"]) as exc:
                exc["proc"].wait(timeout=30)
            self.version_cache[which_prog] = safe_version(exc["stdout"].split()[-1])

        candidate_version = self.version_cache[which_prog]

        if "undef" in candidate_version:
            raise TypeError(
                """Using custom build without tags. Please pull git tags with `git pull origin master --tags`. If building from source, `git fetch upstream "refs/tags/*:refs/tags/*"`."""
            )

        return candidate_version

    def _handle_errors(self, output_data):
        if "error" in output_data:
            if "error_message" not in output_data["error"]:
                error_message = output_data["error"]
                error_type = "internal_error"
            else:
                error_message = output_data["error"]["error_message"]
                error_type = output_data["error"].get("error_type", "unknown_error")
        else:
            error_message = "Unknown error, error message is not found"
            error_type = "internal_error"

        return error_message, error_type

    def compute(self, input_model: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        """
        Runs Psi4 in API mode
        """
        self.found(raise_error=True)
        pversion = parse_version(self.get_version())

        if pversion < parse_version("1.2"):
            raise ResourceError("Psi4 version '{}' not understood.".format(self.get_version()))

        # Location resolution order config.scratch_dir, $PSI_SCRATCH, /tmp
        parent = config.scratch_directory
        if parent is None:
            parent = os.environ.get("PSI_SCRATCH", None)

        error_type = None
        error_message = None
        compute_success = False

        if isinstance(input_model.model.basis, BasisSet):
            raise InputError("QCSchema BasisSet for model.basis not implemented. Use string basis name.")

        # Basis must not be None for HF3c
        old_basis = input_model.model.basis
        input_model.model.__dict__["basis"] = old_basis or ""

        with temporary_directory(parent=parent, suffix="_psi_scratch", messy=config.scratch_messy) as tmpdir:

            caseless_keywords = {k.lower(): v for k, v in input_model.keywords.items()}
            if (input_model.molecule.molecular_multiplicity != 1) and ("reference" not in caseless_keywords):
                input_model.keywords["reference"] = "uhf"

            # Old-style JSON-based command line
            if pversion < parse_version("1.4a2.dev160"):

                # Setup the job
                input_data = input_model.dict(encoding="json")
                input_data["nthreads"] = config.ncores
                input_data["memory"] = int(config.memory * 1024 * 1024 * 1024 * 0.95)  # Memory in bytes
                input_data["success"] = False
                input_data["return_output"] = True

                if input_data["schema_name"] == "qcschema_input":
                    input_data["schema_name"] = "qc_schema_input"

                # Execute the program
                success, output = execute(
                    [which("psi4"), "--scratch", tmpdir, "--json", "data.json"],
                    {"data.json": json.dumps(input_data)},
                    ["data.json"],
                    scratch_directory=tmpdir,
                )

                output_data = input_data.copy()
                if success:
                    output_data = json.loads(output["outfiles"]["data.json"])
                    if "extras" not in output_data:
                        output_data["extras"] = {}

                    # Check QCVars
                    local_qcvars = output_data.pop("psi4:qcvars", None)
                    if local_qcvars:
                        # Edge case where we might already have qcvars, should not happen
                        if "qcvars" in output_data["extras"]:
                            output_data["extras"]["local_qcvars"] = local_qcvars
                        else:
                            output_data["extras"]["qcvars"] = local_qcvars

                    if output_data.get("success", False) is False:
                        error_message, error_type = self._handle_errors(output_data)
                    else:
                        compute_success = True

                else:
                    error_message = output.get("stderr", "No STDERR output")
                    error_type = "execution_error"

                # Reset the schema if required
                output_data["schema_name"] = "qcschema_output"
                output_data.pop("memory", None)
                output_data.pop("nthreads", None)
                output_data["stdout"] = output_data.pop("raw_output", None)

            else:

                if input_model.extras.get("psiapi", False):
                    import psi4

                    orig_scr = psi4.core.IOManager.shared_object().get_default_path()
                    psi4.core.set_num_threads(config.ncores, quiet=True)
                    psi4.set_memory(f"{config.memory}GiB", quiet=True)
                    # psi4.core.IOManager.shared_object().set_default_path(str(tmpdir))
                    if pversion < parse_version("1.5rc1"):  # adjust to where DDD merged
                        # slightly dangerous in that if `qcng.compute({..., psiapi=True}, "psi4")` called *from psi4
                        #   session*, session could unexpectedly get its own files cleaned away.
                        output_data = psi4.schema_wrapper.run_qcschema(input_model).dict()
                    else:
                        output_data = psi4.schema_wrapper.run_qcschema(input_model, postclean=False).dict()
                    # success here means execution returned. output_data may yet be qcel.models.AtomicResult or qcel.models.FailedOperation
                    success = True
                    if output_data.get("success", False):
                        output_data["extras"]["psiapi_evaluated"] = True
                    psi4.core.IOManager.shared_object().set_default_path(orig_scr)
                else:
                    run_cmd = [
                        which("psi4"),
                        "--scratch",
                        str(tmpdir),
                        "--nthread",
                        str(config.ncores),
                        "--memory",
                        f"{config.memory}GiB",
                        "--qcschema",
                        "data.msgpack",
                    ]
                    if config.scratch_messy:
                        run_cmd.append("--messy")
                    input_files = {"data.msgpack": input_model.serialize("msgpack-ext")}
                    success, output = execute(
                        run_cmd, input_files, ["data.msgpack"], as_binary=["data.msgpack"], scratch_directory=tmpdir
                    )
                    if success:
                        output_data = deserialize(output["outfiles"]["data.msgpack"], "msgpack-ext")
                    else:
                        output_data = input_model.dict()

                if success:
                    if output_data.get("success", False) is False:
                        error_message, error_type = self._handle_errors(output_data)
                    else:
                        compute_success = True
                else:
                    error_message = output.get("stderr", "No STDERR output")
                    error_type = "execution_error"

        # Dispatch errors, PSIO Errors are not recoverable for future runs
        if compute_success is False:

            if "PSIO Error" in error_message:
                if "scratch directory" in error_message:
                    # Psi4 cannot access the folder or file
                    raise ResourceError(error_message)
                else:
                    # Likely a random error, worth retrying
                    raise RandomError(error_message)
            elif ("SIGSEV" in error_message) or ("SIGSEGV" in error_message) or ("segmentation fault" in error_message):
                raise RandomError(error_message)
            elif ("TypeError: set_global_option" in error_message) or (error_type == "ValidationError"):
                raise InputError(error_message)
            elif "RHF reference is only for singlets" in error_message:
                raise InputError(error_message)
            else:
                raise UnknownError(error_message)

        # Reset basis
        output_data["model"]["basis"] = old_basis

        # Move several pieces up a level
        output_data["provenance"]["memory"] = round(config.memory, 3)
        output_data["provenance"]["nthreads"] = config.ncores

        # Delete keys
        output_data.pop("return_output", None)

        return AtomicResult(**output_data)
