from typing import Union, Dict, Any

from qcelemental.models import OptimizationInput, AtomicInput, OptimizationResult, Provenance

from qcengine.config import TaskConfig
from qcengine.exceptions import UnknownError, InputError
from qcengine.procedures.nwchem_opt.harvester import harvest_as_atomic_result
from qcengine.programs.nwchem.runner import NWChemHarness
from qcengine.procedures.model import ProcedureHarness


class NWChemDriverProcedure(ProcedureHarness):
    """Structural relaxation using NWChem's optimizer"""

    _defaults = {"name": "NWChemDriver", "procedure": "optimization"}

    class Config(ProcedureHarness.Config):
        pass

    def found(self, raise_error: bool = False) -> bool:
        nwc_harness = NWChemHarness()
        return nwc_harness.found(raise_error)

    def get_version(self) -> str:
        nwc_harness = NWChemHarness()
        return nwc_harness.get_version()

    def build_input_model(self, data: Union[Dict[str, Any], "OptimizationInput"]) -> OptimizationInput:
        return self._build_model(data, OptimizationInput)

    def compute(self, input_data: OptimizationInput, config: TaskConfig) -> "BaseModel":
        nwc_harness = NWChemHarness()
        self.found(raise_error=True)

        # Unify the keywords from the OptimizationInput and QCInputSpecification
        #  Optimization input will override, but don't tell users this as it seems unnecessary
        keywords = input_data.keywords.copy()
        keywords.update(input_data.input_specification.keywords)
        if keywords.get("program", "nwchem").lower() != "nwchem":
            raise InputError("NWChemDriver procedure only works with NWChem")

        # Make an atomic input
        atomic_input = AtomicInput(
            molecule=input_data.initial_molecule,
            driver="energy",
            keywords=keywords,
            **input_data.input_specification.dict(exclude={"driver", "keywords"}),
        )

        # Build the inputs for the job
        job_inputs = nwc_harness.build_input(atomic_input, config)

        # Replace the last line with a "task {} optimize"
        input_file: str = job_inputs["infiles"]["nwchem.nw"].strip()
        beginning, last_line = input_file.rsplit("\n", 1)
        assert last_line.startswith("task")
        last_line = f"task {last_line.split(' ')[1]} optimize"
        job_inputs["infiles"]["nwchem.nw"] = f"{beginning}\n{last_line}"

        # Run it!
        success, dexe = nwc_harness.execute(job_inputs)

        # Check for common errors
        if "There is an error in the input file" in dexe["stdout"]:
            raise InputError(dexe["stdout"])
        if "not compiled" in dexe["stdout"]:
            # recoverable with a different compilation with optional modules
            raise InputError(dexe["stdout"])

        # Parse it
        if success:
            dexe["outfiles"]["stdout"] = dexe["stdout"]
            dexe["outfiles"]["stderr"] = dexe["stderr"]
            return self.parse_output(dexe["outfiles"], input_data)
        else:
            raise UnknownError(dexe["stdout"])

    def parse_output(self, outfiles: Dict[str, str], input_model: OptimizationInput) -> OptimizationResult:

        # Get the stdout from the calculation (required)
        stdout = outfiles.pop("stdout")
        stderr = outfiles.pop("stderr")

        # Parse out the atomic results from the file
        atomic_results = harvest_as_atomic_result(input_model, stdout)

        # Isolate the converged result
        final_step = atomic_results[-1]

        return OptimizationResult(
            initial_molecule=input_model.initial_molecule,
            input_specification=input_model.input_specification,
            final_molecule=final_step.molecule,
            trajectory=atomic_results,
            energies=[float(r.extras["qcvars"]["CURRENT ENERGY"]) for r in atomic_results],
            stdout=stdout,
            stderr=stderr,
            success=True,
            provenance=Provenance(creator="NWChemRelax", version=self.get_version(), routine="nwchem_opt"),
        )
