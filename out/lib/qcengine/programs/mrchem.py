"""
Calls the MRChem executable.
"""
import copy
import json
import logging
import pprint
import sys
from collections import Counter
from functools import reduce
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Optional, Tuple

from qcelemental.models import AtomicResult
from qcelemental.util import safe_version, which

from ..exceptions import InputError, RandomError, UnknownError
from ..util import create_mpi_invocation, execute, popen, temporary_directory
from .model import ProgramHarness

if TYPE_CHECKING:
    from qcelemental.models import AtomicInput

    from ..config import TaskConfig

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)
logger = logging.getLogger(__name__)


class MRChemHarness(ProgramHarness):

    _defaults = {
        "name": "MRChem",
        "scratch": False,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": True,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        """Whether MRChem harness is ready for operation.

        Parameters
        ----------
        raise_error: bool
            Passed on to control negative return between False and ModuleNotFoundError raised.

        Returns
        -------
        bool
            If mrchem and mrchem.x are found, returns True.
            If raise_error is False and MRChem is missing, returns False.
            If raise_error is True and MRChem is missing, the error message is raised.

        """
        mrchem_x = which(
            "mrchem.x",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via https://mrchem.readthedocs.io",
        )
        mrchem = which(
            "mrchem",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via https://mrchem.readthedocs.io",
        )

        return mrchem and mrchem_x

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which("mrchem.x")
        if which_prog not in self.version_cache:
            with popen([which_prog, "--version"]) as exc:
                exc["proc"].wait(timeout=30)
            self.version_cache[which_prog] = safe_version(exc["stdout"].split()[-1])

        candidate_version = self.version_cache[which_prog]

        return candidate_version

    def compute(self, input_model: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        """
        Runs MRChem in executable mode
        """
        self.found(raise_error=True)

        # Location resolution order config.scratch_dir, /tmp
        parent = config.scratch_directory

        error_message = None
        compute_success = False

        job_input = self.build_input(input_model, config)
        input_data = copy.deepcopy(job_input["mrchem_json"])
        output_data = {
            "keywords": input_data,
            "schema_name": "qcschema_output",
            "schema_version": 1,
            "model": input_model.model,
            "molecule": input_model.molecule,
            "driver": input_model.driver,
        }

        with temporary_directory(parent=parent, suffix="_mrchem_scratch") as tmpdir:
            # create folders
            for d in job_input["folders"]:
                if not Path(d).exists():
                    Path(d).mkdir()

            # Execute the program
            success, output = execute(
                command=job_input["command"] + ["data.json"],
                infiles={"data.json": json.dumps(job_input["mrchem_json"])},
                outfiles=["data.json"],
                scratch_directory=tmpdir,
            )

            if success:
                output_data["stdout"] = output["stdout"]
                # get data from the MRChem JSON output and transfer it to the QCSchema output
                mrchem_json = json.loads(output["outfiles"]["data.json"])
                mrchem_output = mrchem_json["output"]
                output_data["success"] = mrchem_output["success"]
                output_data["provenance"] = mrchem_output["provenance"]
                # update the "routine" under "provenance"
                output_data["provenance"]["routine"] = " ".join(job_input["command"])

                # fill up properties
                output_data["properties"] = extract_properties(mrchem_output)

                # prepare a list of computed response properties
                known_rsp_props = [
                    ("dipole_moment", "vector"),
                    ("quadrupole_moment", "tensor"),
                    ("polarizability", "tensor"),
                    ("magnetizability", "tensor"),
                    ("nmr_shielding", "tensor"),
                ]
                computed_rsp_props = [
                    ("properties", x, y, z)
                    for x, z in known_rsp_props
                    if x in mrchem_output["properties"]
                    for y in mrchem_output["properties"][x].keys()
                ]

                # fill up extras:
                # * under "raw_output" the whole JSON output from MRChem
                # * under "properties" all the properties computed by MRChem
                output_data["extras"] = {
                    "raw_output": mrchem_json,
                    "properties": {
                        f"{ks[1]}": {f"{ks[2]}": _nested_get(mrchem_output, ks)} for ks in computed_rsp_props
                    },
                }

                # fill up return_result
                if input_model.driver == "energy":
                    output_data["return_result"] = mrchem_output["properties"]["scf_energy"]["E_tot"]
                elif input_model.driver == "properties":
                    output_data["return_result"] = {
                        f"{ks[1]}": {f"{ks[2]}": _nested_get(mrchem_output, ks)} for ks in computed_rsp_props
                    }
                else:
                    raise InputError(f"Driver {input_model.driver} not implemented for MRChem.")

                compute_success = mrchem_output["success"]

            else:
                output_data["stderr"] = output["stderr"]
                output_data["error"] = {
                    "error_message": output["stderr"],
                    "error_type": "execution_error",
                }

        # Dispatch errors, PSIO Errors are not recoverable for future runs
        if compute_success is False:

            if ("SIGSEV" in error_message) or ("SIGSEGV" in error_message) or ("segmentation fault" in error_message):
                raise RandomError(error_message)
            else:
                raise UnknownError(error_message)

        return AtomicResult(**output_data)

    def build_input(self, input_model: "AtomicInput", config: "TaskConfig") -> Dict[str, Any]:
        with popen([which("mrchem"), "--module"]) as exc:
            exc["proc"].wait(timeout=30)
        sys.path.append(exc["stdout"].split()[-1])
        from mrchem import translate_input, validate

        mrchemrec = {
            "scratch_directory": config.scratch_directory,
        }

        opts = copy.deepcopy(input_model.keywords)

        # Handle molecule
        _, moldict = input_model.molecule.to_string(
            dtype="mrchem", units=opts.get("world_unit", "bohr"), return_data=True
        )
        opts["Molecule"] = moldict["keywords"]

        if "WaveFunction" in opts.keys():
            opts["WaveFunction"]["method"] = input_model.model.method
        else:
            opts["WaveFunction"] = {"method": input_model.model.method}
        # Log the job settings as constructed from the input model
        logger.debug("JOB_OPTS from InputModel")
        logger.debug(pp.pformat(opts))

        try:
            opts = validate(ir_in=opts)
        except Exception as e:
            raise InputError(f"Failure preparing input to MRChem\n {str(e)}")
        # Log the validated job settings
        logger.debug("JOB_OPTS after validation")
        logger.debug(pp.pformat(opts))
        mrchemrec["folders"] = [
            opts["SCF"]["path_checkpoint"],
            opts["SCF"]["path_orbitals"],
            opts["Response"]["path_checkpoint"],
            opts["Response"]["path_orbitals"],
            opts["Plotter"]["path"],
        ]

        try:
            opts = translate_input(opts)
        except Exception as e:
            raise InputError(f"Failure preparing input to MRChem\n {str(e)}")
        opts["printer"]["file_name"] = "data.inp"
        # Log the final job settings
        logger.debug("JOB_OPTS after translation")
        logger.debug(pp.pformat(opts))

        mrchemrec["mrchem_json"] = {
            "input": opts,
        }

        # Determine the command
        if config.use_mpiexec:
            mrchemrec["command"] = create_mpi_invocation(which("mrchem.x"), config)
            logger.info(f"Launching with mpiexec: {' '.join(mrchemrec['command'])}")
        else:
            mrchemrec["command"] = [which("mrchem.x")]

        return mrchemrec


def extract_properties(mrchem_output: Dict[str, Any]) -> Dict[str, Any]:
    """Translate MRChem output to QCSChema properties.

    Parameters
    ----------

    Returns
    -------
    """

    occs = Counter(mrchem_output["properties"]["orbital_energies"]["spin"])
    properties = {
        "calcinfo_nmo": 2 * occs["p"] + occs["a"] + occs["b"],
        "calcinfo_nalpha": occs["p"] + occs["a"],
        "calcinfo_nbeta": occs["p"] + occs["b"],
        "calcinfo_natom": len(mrchem_output["properties"]["geometry"]),
        "return_energy": mrchem_output["properties"]["scf_energy"]["E_tot"],
        "scf_one_electron_energy": mrchem_output["properties"]["scf_energy"]["E_kin"]
        + mrchem_output["properties"]["scf_energy"]["E_en"]
        + mrchem_output["properties"]["scf_energy"]["E_next"]
        + mrchem_output["properties"]["scf_energy"]["E_eext"],
        "scf_two_electron_energy": mrchem_output["properties"]["scf_energy"]["E_ee"]
        + mrchem_output["properties"]["scf_energy"]["E_x"]
        + mrchem_output["properties"]["scf_energy"]["E_xc"],
        "nuclear_repulsion_energy": mrchem_output["properties"]["scf_energy"]["E_nn"],
        "scf_xc_energy": mrchem_output["properties"]["scf_energy"]["E_xc"],
        "scf_total_energy": mrchem_output["properties"]["scf_energy"]["E_tot"],
        "scf_iterations": len(mrchem_output["scf_calculation"]["scf_solver"]["cycles"]),
        "scf_dipole_moment": mrchem_output["properties"]["dipole_moment"]["dip-1"]["vector"],
    }

    return properties


def _nested_get(d: Dict[str, Any], ks: Tuple[str, ...]) -> Optional[Any]:
    """Get value from a nested dictionary.

    Parameters
    ----------
    d : Dict[str, Any]
    ks : str

    Returns
    -------
    v : Optional[Any]

    Notes
    -----
    Adapted from: https://stackoverflow.com/a/40675868/2528668
    """

    def _func(x: Dict[str, Any], k: str) -> Optional[Dict[str, Any]]:
        return x.get(k, None) if isinstance(x, dict) else None

    return reduce(_func, ks, d)  # type: ignore
