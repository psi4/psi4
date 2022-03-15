import logging
import sys
import traceback
from io import StringIO
from typing import Any, Dict, Union

import numpy as np
from qcelemental.models import OptimizationInput, OptimizationResult
from qcelemental.util import which_import

import qcengine

from ..config import TaskConfig
from .model import ProcedureHarness


class BernyProcedure(ProcedureHarness):
    _defaults = {"name": "Berny", "procedure": "optimization"}

    def found(self, raise_error: bool = False) -> bool:
        return which_import(
            "berny",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `pip install pyberny`.",
        )

    def build_input_model(self, data: Union[Dict[str, Any], "OptimizationInput"]) -> "OptimizationInput":
        return self._build_model(data, OptimizationInput)

    def compute(self, input_data: "OptimizationInput", config: "TaskConfig") -> "OptimizationResult":
        try:
            import berny
        except ModuleNotFoundError:
            raise ModuleNotFoundError("Could not find Berny in the Python path.")

        # Get berny version from the installed package, use setuptools'
        # pkg_resources for python < 3.8
        if sys.version_info >= (3, 8):
            from importlib.metadata import distribution
        else:
            from pkg_resources import get_distribution as distribution
        berny_version = distribution("pyberny").version

        # Berny uses the stdlib logging module and by default uses per-module
        # loggers. For QCEngine, we create one logger per BernyProcedure
        # instance, by using the instance's id(), and send all logging messages
        # to a string stream
        log_stream = StringIO()
        log = logging.getLogger(f"{__name__}.{id(self)}")
        log.addHandler(logging.StreamHandler(log_stream))
        log.setLevel("INFO")

        input_data = input_data.dict()
        geom_qcng = input_data["initial_molecule"]
        comput = {**input_data["input_specification"], "molecule": geom_qcng}
        program = input_data["keywords"].pop("program")
        trajectory = []
        output_data = input_data.copy()
        try:
            # Pyberny uses angstroms for the Cartesian geometry, but atomic
            # units for everything else, including the gradients (hartree/bohr).
            geom_berny = berny.Geometry(geom_qcng["symbols"], geom_qcng["geometry"] / berny.angstrom)
            opt = berny.Berny(geom_berny, logger=log, **input_data["keywords"])
            for geom_berny in opt:
                geom_qcng["geometry"] = np.stack(geom_berny.coords * berny.angstrom)
                ret = qcengine.compute(comput, program)
                trajectory.append(ret.dict())
                opt.send((ret.properties.return_energy, ret.return_result))
        except Exception:
            output_data["success"] = False
            output_data["error"] = {"error_type": "unknown", "error_message": f"Berny error:\n{traceback.format_exc()}"}
        else:
            output_data["success"] = True
        output_data.update(
            {
                "schema_name": "qcschema_optimization_output",
                "final_molecule": trajectory[-1]["molecule"],
                "energies": [r["properties"]["return_energy"] for r in trajectory],
                "trajectory": trajectory,
                "provenance": {"creator": "Berny", "routine": "berny.Berny", "version": berny_version},
                "stdout": log_stream.getvalue(),  # collect logged messages
            }
        )
        if output_data["success"]:
            output_data = OptimizationResult(**output_data)
        return output_data
