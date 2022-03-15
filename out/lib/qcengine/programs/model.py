import abc
import logging
from typing import Any, Dict, List, Optional, Tuple

from pydantic import BaseModel
from qcelemental.models import AtomicInput, AtomicResult

from qcengine.exceptions import KnownErrorException

logger = logging.getLogger(__name__)


class ProgramHarness(BaseModel, abc.ABC):

    _defaults: Dict[str, Any] = {}
    name: str
    scratch: bool
    thread_safe: bool
    thread_parallel: bool
    node_parallel: bool
    managed_memory: bool
    extras: Optional[Dict[str, Any]]

    class Config:
        allow_mutation: False
        extra: "forbid"

    def __init__(self, **kwargs):
        super().__init__(**{**self._defaults, **kwargs})

    @abc.abstractmethod
    def compute(self, input_data: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        pass

    @staticmethod
    @abc.abstractmethod
    def found(raise_error: bool = False) -> bool:
        """
        Checks if the program can be found.

        Parameters
        ----------
        raise_error : bool, optional
            If True, raises an error if the program cannot be found.

        Returns
        -------
        bool
            Returns True if the program was found, False otherwise.
        """

    ## Utility

    def get_version(self) -> str:
        """Finds program, extracts version, returns normalized version string.

        Returns
        -------
        str
            Return a valid, safe python version string.
        """

    ## Computers

    def build_input(
        self, input_model: "AtomicInput", config: "TaskConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:
        raise ValueError("build_input is not implemented for {}.", self.__class__)

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:
        raise ValueError("execute is not implemented for {}.", self.__class__)

    def parse_output(self, outfiles: Dict[str, str], input_model: "AtomicInput") -> "AtomicResult":
        raise ValueError("parse_output is not implemented for {}.", self.__class__)


class ErrorCorrectionProgramHarness(ProgramHarness, abc.ABC):
    """Base class for Harnesses that include logic to correct common errors

    Classes which implement this Harness must override the :meth:`_compute` method
    rather than :meth:`compute`. The ``compute`` method from this class will make
    calls to ``_compute`` with different actions as it attempts to correct errors in the
    input files.

    The error corrections are defined by first implementing a :class:`KnownErrorException`
    that contains logic to determine if a certain error occurs and returns an appropriate
    update to the keywords of an AtomicInput.

    Then, modify the ``_compute`` method of your Harness to run the ``detect_error`` method
    of each ``KnownErrorException`` for the Harness after it finishes performing the computation.
    The ``detect_error`` method will raise an exception that is caught by the
    ``ErrorCorrectionProgramHarness`` and used to determine if/how to re-run the computation.
    """

    def _compute(self, input_data: AtomicInput, config: "TaskConfig") -> AtomicResult:
        raise NotImplementedError()

    def compute(self, input_data: AtomicInput, config: "TaskConfig") -> AtomicResult:
        # Get the error correction configuration
        error_policy = input_data.protocols.error_correction

        # Create a local copy of the input data
        local_input_data = input_data

        # Run the method and, if it fails, assess if the failure is restartable
        observed_errors = {}  # Errors that have been observed previously
        while True:
            try:
                result = self._compute(local_input_data, config)
                break
            except KnownErrorException as e:
                logger.info(f"Caught a {type(e)} error.")

                # Determine whether this specific type of error is allowed
                correction_allowed = error_policy.allows(e.error_name)
                if not correction_allowed:
                    logger.info(f'Error correction for "{e.error_name}" is not allowed')
                    raise e
                logger.info(f'Error correction for "{e.error_name}" is allowed')

                # Check if it has run before
                # TODO (wardlt): Should we allow errors to be run >1 time?
                previously_run = e.error_name in observed_errors
                if previously_run:
                    logger.info(
                        "Error has been observed before and mitigation did not fix the issue. Raising exception"
                    )
                    raise e

                # Generate and apply the updated keywords
                keyword_updates = e.create_keyword_update(local_input_data)
                new_keywords = local_input_data.keywords.copy()
                new_keywords.update(keyword_updates)
                local_input_data = AtomicInput(**local_input_data.dict(exclude={"keywords"}), keywords=new_keywords)

                # Store the error details and mitigations employed
                observed_errors[e.error_name] = {"details": e.details, "keyword_updates": keyword_updates}

        # Add the errors observed and corrected for, if any
        if len(observed_errors) > 0:
            result.extras["observed_errors"] = observed_errors
        return result
