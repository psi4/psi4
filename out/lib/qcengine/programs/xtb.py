"""
Harness for the extended tight binding (xtb) package.
This implementation interfaces with the xtb C-API via the xtb-python project,
which is providing Python bindings and native extensions for common computational
chemistry frameworks, like QCEngine.

Therefore, this harness only has to provide a thin wrapper to integrate xtb.
For more information on xtb-python and the actual QCSchema integration,
visit `its documentation <https://xtb-python.readthedocs.io>`_.
"""

from typing import Dict

from qcelemental.models import AtomicInput, AtomicResult
from qcelemental.util import safe_version, which_import

from ..config import TaskConfig
from .model import ProgramHarness


class XTBHarness(ProgramHarness):
    """Calculation harness for the extended tight binding (xtb) package."""

    _defaults = {
        "name": "xtb",
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
        """Check for the availability of the Python API of xtb"""

        return which_import(
            "xtb",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install xtb-python -c conda-forge`.",
        )

    def get_version(self) -> str:
        """Return the currently used version of xtb-python"""
        self.found(raise_error=True)

        which_prog = which_import("xtb")
        if which_prog not in self.version_cache:
            import xtb

            self.version_cache[which_prog] = safe_version(xtb.__version__)

        return self.version_cache[which_prog]

    def compute(self, input_data: AtomicInput, config: TaskConfig) -> AtomicResult:
        """
        Actual interface to the xtb package. The compute function is just a thin
        wrapper around the native QCSchema interface of xtb-python.
        """

        self.found(raise_error=True)

        import xtb
        from xtb.qcschema.harness import run_qcschema

        # Run the Harness
        output = run_qcschema(input_data)

        # Make sure all keys from the initial input spec are sent along
        output.extras.update(input_data.extras)
        return output
