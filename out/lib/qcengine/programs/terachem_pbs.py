"""
Calls TeraChem in its "server mode" via a protobuf interface.
"""
import logging
import os
from typing import TYPE_CHECKING

from qcelemental.models import AtomicResult
from qcelemental.util import which_import

from .model import ProgramHarness

if TYPE_CHECKING:
    from qcelemental.models import AtomicInput

    from ..config import TaskConfig

logger = logging.getLogger(__name__)


class TeraChemPBSHarness(ProgramHarness):

    _defaults = {
        "name": "terachem_pbs",
        "scratch": False,
        "thread_safe": False,
        "thread_parallel": False,
        "node_parallel": False,
        "managed_memory": True,
    }

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        """Whether TeraChemPBS harness is ready for operation.
        Parameters
        ----------
        raise_error: bool
            Passed on to control negative return between False and ModuleNotFoundError raised.
        Returns
        -------
        bool
            If tcpb package is found and server available, returns True.
            If raise_error is False and tcpb package missing and/or server us unavailable, returns False.
            If raise_error is True and tcpb package missing and/or server us unavailable, the error message is raised.
        """
        tcpb_pkg_available = which_import(
            "tcpb",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="TeraChem protobuf client package (tcpb) not found. Please install tcpb>=0.7.0.",
        )
        if not tcpb_pkg_available:
            return False

        from tcpb.exceptions import ServerError
        from tcpb.tcpb import TCProtobufClient

        try:
            with TCProtobufClient(
                host=os.getenv("TERACHEM_PBS_HOST"), port=int(os.getenv("TERACHEM_PBS_PORT"))
            ) as client:
                return client.is_available()
        except TypeError as e:
            # TERACHEM_PBS_HOST/PORT environment variables unset
            msg = "Environment variables 'TERACHEM_PBS_HOST' and 'TERACHEM_PBS_PORT' must be set!"
            logger.error(msg)
            if raise_error:
                raise ValueError(msg) from e
        except ServerError as e:
            msg = (
                f"Unable to connect to TeraChem server at "
                f"{os.getenv('TERACHEM_PBS_HOST')}:{os.getenv('TERACHEM_PBS_PORT')}"
            )
            logger.error(msg)
            if raise_error:
                raise OSError(msg) from e
        return False

    def get_version(self) -> str:
        """Returns version of TeraChem Protocol Buffer Server"""
        try:
            import tcpb
        except ModuleNotFoundError:
            return None
        else:
            try:
                return tcpb.__version__
            except AttributeError:
                return None

    def compute(self, input_model: "AtomicInput", config: "TaskConfig" = None) -> "AtomicResult":
        """
        Submit AtomicInput to TeraChem running in "server mode"
        """
        self.found(raise_error=True)

        from tcpb.tcpb import TCProtobufClient

        with TCProtobufClient(host=os.getenv("TERACHEM_PBS_HOST"), port=int(os.getenv("TERACHEM_PBS_PORT"))) as client:
            return client.compute(input_model)
