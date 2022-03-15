"""Known errors for NWChem"""
import logging
from typing import Any, Dict

from qcelemental.models import AtomicInput

from qcengine.exceptions import SimpleKnownErrorException
from qcengine.programs.nwchem.germinate import xc_functionals

logger = logging.getLogger(__name__)


class GeomBinvrError(SimpleKnownErrorException):
    """GeomBinvr Error"""

    error_name = "geom_binvr"
    description = "Error when generating redundant atomic coordinates. Often fixed by turning autoz off (`noautoz`)"

    @classmethod
    def _detect(cls, outputs: Dict[str, str]) -> bool:
        return "geom_binvr: #indep variables incorrect" in outputs["stdout"]

    def create_keyword_update(self, input_data: AtomicInput) -> Dict[str, Any]:
        return {"geometry__noautoz": True}


class ConvergenceFailedError(SimpleKnownErrorException):
    """Failed to converge on electronic steps. This will increase the limit by 4x"""

    error_name = "convergence_failed"
    description = "The computation failed to converge."

    @classmethod
    def _detect(cls, outputs: Dict[str, str]) -> bool:
        return (
            "This type of error is most commonly" in outputs["stdout"] and "convergence criteria" in outputs["stdout"]
        )

    def create_keyword_update(self, input_data: AtomicInput) -> Dict[str, Any]:
        # Fit the correct keyword we are looking to update is different for different methods
        method = input_data.model.method
        use_tce = input_data.keywords.get("qc_module", False)

        if method == "dft" or method.split()[0] in xc_functionals:
            if "dft__iterations" in input_data.keywords:
                kwd = "dft__iterations"
                cur_iter = input_data.keywords["dft__iterations"]
            elif "dft__maxiter" in input_data.keywords:
                kwd = "dft__maxiter"
                cur_iter = input_data.keywords["dft__maxiter"]
            else:
                kwd = "dft__maxiter"
                cur_iter = 20  # The NWChem default
        elif method in ["scf", "hf", "mp2"]:
            kwd = "scf__maxiter"
            cur_iter = input_data.keywords.get(kwd, 8)
        elif method.startswith("ccsd") and not use_tce:
            kwd = "ccsd__maxiter"
            cur_iter = input_data.keywords.get(kwd, 20)
        elif method.startswith("ccsd") and use_tce:
            kwd = "tce__maxiter"
            cur_iter = input_data.keywords.get(kwd, 100)
        else:
            raise ValueError(f'Method "{method}" is not yet supported')

        return {kwd: cur_iter * 4}


# List of all the known errors
all_errors = [GeomBinvrError, ConvergenceFailedError]
