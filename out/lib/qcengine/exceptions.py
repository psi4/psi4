import traceback
from typing import Any, Dict, Optional

from qcelemental.models import AtomicInput


class QCEngineException(Exception):
    """
    Base QCEngine exception, should never be called explicitly.
    """

    error_type = "base_error"
    header = "QCEngine Base Error"

    def __init__(self, message: str):

        # Call the base class constructor with the parameters it needs
        super().__init__(message)

        # Now for your custom code...
        self.raw_message = message
        self.traceback = traceback.format_exc()

    @property
    def error_message(self) -> str:
        return f"{self.header}: {self.raw_message}"


class UnknownError(QCEngineException):
    """
    Unknown QCEngine error, the type was not able to be specified.
    """

    error_type = "unknown_error"
    header = "QCEngine Unknown Error"


class InputError(QCEngineException):
    """
    Incorrect user parameters, not recoverable. Also may indicate a version issue.
    """

    error_type = "input_error"
    header = "QCEngine Input Error"


class ResourceError(QCEngineException):
    """
    Not enough resources for computation such as not enough memory, cores, or disk was not available.
    Recoverable on different compute hardware.
    """

    error_type = "resource_error"
    header = "QCEngine Resource Error"


class ConvergenceError(QCEngineException):
    """
    Failed iteration convergence error, likely recoverable with tweaked parameters.
    """

    error_type = "convergence_error"
    header = "QCEngine Convergence Error"


class RandomError(QCEngineException):
    """
    Likely recoverable errors such as segmentation faults or disk io problems.
    """

    error_type = "random_error"
    header = "QCEngine Random Error"


class KnownErrorException(QCEngineException):
    """Base class for detecting errors in the output of a ProgramHarness

    Subclasses for this class must define `error_name` and `description`, which provide
    a short name for the error used in the error correction code and a human-readable
    description of the error and how it is corrected.

    Subclasses for this exception must fulfill the `detect_error` class method,
     which detects if an error has occurred and raises an exception if one was detected.
     When raising the exception, the `detect_error` routine can optionally provide details
     about the error that are needed to correction.
    """

    error_type = "known_error"
    header = "QCEngine Known Error"

    error_name: str  # Unique name for this error. Used in autodetection logic
    description: str  # Human-readable description of the error.
    details: Optional[Dict[str, Any]]  # Any details for this error. Used for user feedback and autocorrection logic

    def __init__(self, details: Optional[dict] = None):
        super().__init__(self.description)
        self.details = details

    @classmethod
    def detect_error(cls, outputs: Dict[str, str]):
        """Detect whether this operation throws an error
        Args:
            outputs (dict): Output files needed to check for errors
        Raises:
            KnownErrorException: Error, if one was detected
        """
        raise NotImplementedError()

    def create_keyword_update(self, input_data: AtomicInput) -> Dict[str, Any]:
        """Create an keyword used to the update the dictionary given observed error

        Parameters
        ----------
            input_data
                Input specification used to perform the calculation that failed

        Returns
        -------
            Dictionary used to update the keywords field of ``input_data``
        """
        raise NotImplementedError()


class SimpleKnownErrorException(KnownErrorException):
    """Subclass for errors with simple detection logic.

    Most useful for error types that do not need any additional details to correct."""

    @classmethod
    def detect_error(cls, outputs: Dict[str, str]):
        if cls._detect(outputs):
            raise cls()

    @classmethod
    def _detect(cls, outputs: Dict[str, str]) -> bool:
        """Detect whether an error is present"""
        raise NotImplementedError()
