"""
Harness for the DFT-D4 dispersion correction.
This implementation interfaces with the dftd4 Python-API, which provides
native support for QCSchema.

Therefore, this harness only has to provide a thin wrapper to integrate dftd4.
"""

from typing import Dict

from qcelemental.models import AtomicInput, AtomicResult
from qcelemental.util import safe_version, which_import

from ..config import TaskConfig
from ..exceptions import InputError
from .empirical_dispersion_resources import from_arrays, get_dispersion_aliases
from .model import ProgramHarness


class DFTD4Harness(ProgramHarness):
    """Calculation harness for the DFT-D4 dispersion correction."""

    _defaults = {
        "name": "dftd4",
        "scratch": False,
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
        """Check for the availability of the Python API of dftd4"""

        return which_import(
            "dftd4",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via a dftd4 version with enabled Python API",
        )

    def get_version(self) -> str:
        """Return the currently used version of dftd4"""
        self.found(raise_error=True)

        which_prog = which_import("dftd4")
        if which_prog not in self.version_cache:
            import dftd4

            self.version_cache[which_prog] = safe_version(dftd4.__version__)

        return self.version_cache[which_prog]

    def compute(self, input_model: AtomicInput, config: TaskConfig) -> AtomicResult:
        """
        Actual interface to the dftd4 package. The compute function is just a thin
        wrapper around the native QCSchema interface of the dftd4 Python-API.
        """

        self.found(raise_error=True)

        import dftd4
        from dftd4.qcschema import run_qcschema

        # strip engine hint
        input_data = input_model.dict()
        method = input_model.model.method
        if method.startswith("d4-"):
            method = method[3:]
            input_data["model"]["method"] = method
        qcvkey = method.upper() if method is not None else None

        # send `from_arrays` the dftd4 behavior of functional specification overrides explicit parameters specification
        # * differs from dftd3 harness behavior where parameters extend or override functional
        # * stash the resolved plan in extras or, if errored, leave it for the proper dftd4 api to reject
        param_tweaks = None if method else input_model.keywords.get("params_tweaks", None)
        try:
            planinfo = from_arrays(
                verbose=1,
                name_hint=method,
                level_hint=input_model.keywords.get("level_hint", None),
                param_tweaks=param_tweaks,
                dashcoeff_supplement=input_model.keywords.get("dashcoeff_supplement", None),
            )
        except InputError:
            pass
        else:
            input_data["extras"]["info"] = planinfo

        # strip dispersion level from method
        for alias, d4 in get_dispersion_aliases().items():
            if d4 == "d4bjeeqatm" and method.lower().endswith(alias):
                method = method[: -(len(alias) + 1)]
                input_data["model"]["method"] = method

        # consolidate dispersion level aliases
        level_hint = input_model.keywords.get("level_hint", None)
        if level_hint and get_dispersion_aliases()[level_hint.lower()] == "d4bjeeqatm":
            level_hint = "d4"
            input_data["keywords"]["level_hint"] = level_hint

        input_model = AtomicInput(**input_data)

        # Run the Harness
        output = run_qcschema(input_model)

        if "info" in output.extras:
            qcvkey = output.extras["info"]["fctldash"].upper()

        calcinfo = {}
        energy = output.properties.return_energy
        calcinfo["CURRENT ENERGY"] = energy
        calcinfo["DISPERSION CORRECTION ENERGY"] = energy
        if qcvkey:
            calcinfo[f"{qcvkey} DISPERSION CORRECTION ENERGY"] = energy

        if output.driver == "gradient":
            gradient = output.return_result
            calcinfo["CURRENT GRADIENT"] = gradient
            calcinfo["DISPERSION CORRECTION GRADIENT"] = gradient
            if qcvkey:
                calcinfo[f"{qcvkey} DISPERSION CORRECTION GRADIENT"] = gradient

        if output.keywords.get("pair_resolved", False):
            pw2 = output.extras["dftd4"]["additive pairwise energy"]
            pw3 = output.extras["dftd4"]["non-additive pairwise energy"]
            assert abs(pw2.sum() + pw3.sum() - energy) < 1.0e-8, f"{pw2.sum()} + {pw3.sum()} != {energy}"
            calcinfo["2-BODY DISPERSION CORRECTION ENERGY"] = pw2.sum()
            calcinfo["3-BODY DISPERSION CORRECTION ENERGY"] = pw3.sum()
            calcinfo["2-BODY PAIRWISE DISPERSION CORRECTION ANALYSIS"] = pw2
            calcinfo["3-BODY PAIRWISE DISPERSION CORRECTION ANALYSIS"] = pw3

        output.extras["qcvars"] = calcinfo

        return output
