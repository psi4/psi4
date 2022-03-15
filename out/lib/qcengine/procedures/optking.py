from typing import Any, Dict, Union

from qcelemental.models import OptimizationInput, OptimizationResult
from qcelemental.util import safe_version, which_import

from .model import ProcedureHarness


class OptKingProcedure(ProcedureHarness):

    _defaults = {"name": "OptKing", "procedure": "optimization"}

    version_cache: Dict[str, str] = {}

    class Config(ProcedureHarness.Config):
        pass

    def found(self, raise_error: bool = False) -> bool:
        return which_import(
            "optking",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install optking -c psi4/label/dev`.",
        )

    def build_input_model(self, data: Union[Dict[str, Any], "OptimizationInput"]) -> "OptimizationInput":
        return self._build_model(data, OptimizationInput)

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which_import("optking")
        if which_prog not in self.version_cache:
            import optking

            self.version_cache[which_prog] = safe_version(optking.__version__)

        return self.version_cache[which_prog]

    def compute(self, input_model: "OptimizationInput", config: "TaskConfig") -> "Optimization":
        if self.found(raise_error=True):
            import optking

        input_data = input_model.dict()

        # Set retries to two if zero while respecting local_config
        local_config = config.dict()
        local_config["retries"] = local_config.get("retries", 2) or 2
        input_data["input_specification"]["extras"]["_qcengine_local_config"] = local_config

        # Run the program
        output_data = optking.optwrapper.optimize_qcengine(input_data)

        output_data["schema_name"] = "qcschema_optimization_output"
        output_data["input_specification"]["extras"].pop("_qcengine_local_config", None)
        if output_data["success"]:
            output_data = OptimizationResult(**output_data)

        return output_data
