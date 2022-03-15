import abc
from typing import Any, Dict, Union

from pydantic import BaseModel

from ..util import model_wrapper


class ProcedureHarness(BaseModel, abc.ABC):

    name: str
    procedure: str

    class Config:
        allow_mutation: False
        extra: "forbid"

    def __init__(self, **kwargs):
        super().__init__(**{**self._defaults, **kwargs})

    @abc.abstractmethod
    def build_input_model(self, data: Union[Dict[str, Any], "BaseModel"], raise_error: bool = True) -> "BaseModel":
        """
        Build and validate the input model, passes if the data was a normal BaseModel input.

        Parameters
        ----------
        data : Union[Dict[str, Any], 'BaseModel']
            A data blob to construct the model from or the input model itself
        raise_error : bool, optional
            Raise an error or not if the operation failed.

        Returns
        -------
        BaseModel
            The input model for the procedure.
        """

    @abc.abstractmethod
    def compute(self, input_data: "BaseModel", config: "TaskConfig") -> "BaseModel":
        pass

    @abc.abstractmethod
    def found(self, raise_error: bool = False) -> bool:
        """
        Checks if the program can be found.

        Returns
        -------
        bool
            If the proceudre was found or not.
        """

    def _build_model(self, data: Dict[str, Any], model: "BaseModel") -> "BaseModel":
        """
        Quick wrapper around util.model_wrapper for inherited classes
        """

        return model_wrapper(data, model)

    def get_version(self) -> str:
        """Finds procedure, extracts version, returns normalized version string.
        Returns
        -------
        str
            Return a valid, safe python version string.
        """
