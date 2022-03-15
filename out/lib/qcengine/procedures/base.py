"""
Imports the various procedure backends
"""

from typing import Set

from ..exceptions import InputError, ResourceError
from .berny import BernyProcedure
from .geometric import GeometricProcedure
from .nwchem_opt import NWChemDriverProcedure
from .optking import OptKingProcedure
from .torsiondrive import TorsionDriveProcedure

__all__ = ["register_procedure", "get_procedure", "list_all_procedures", "list_available_procedures"]

procedures = {}


def register_procedure(entry_point: "ProcedureHarness") -> None:
    """
    Register a new ProcedureHarness with QCEngine
    """

    name = entry_point.name
    if name.lower() in procedures.keys():
        raise ValueError("{} is already a registered procedure.".format(name))

    procedures[name.lower()] = entry_point


def get_procedure(name: str) -> "ProcedureHarness":
    """
    Returns a procedures executor class
    """

    name = name.lower()

    if name not in procedures:
        raise InputError(f"Procedure {name} is not registered to QCEngine.")

    ret = procedures[name]
    if not ret.found():
        raise ResourceError(f"Procedure {name} is registered with QCEngine, but cannot be found.")

    return ret


def list_all_procedures() -> Set[str]:
    """
    List all procedures registered by QCEngine.
    """
    return set(procedures.keys())


def list_available_procedures() -> Set[str]:
    """
    List all procedures that can be exectued (found) by QCEngine.
    """

    ret = set()
    for k, p in procedures.items():
        if p.found():
            ret.add(k)

    return ret


register_procedure(GeometricProcedure())
register_procedure(OptKingProcedure())
register_procedure(BernyProcedure())
register_procedure(NWChemDriverProcedure())
register_procedure(TorsionDriveProcedure())
