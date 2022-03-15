from typing import Any, Dict

from qcengine.exceptions import InputError


def muster_modelchem(method: str, derint: int) -> Dict[str, Any]:
    """Converts the QC method into GAMESS keywords."""

    method = method.lower()
    opts = {}

    runtyp = {
        0: "energy",
        1: "gradient",
        2: "hessian",
        #'properties': 'prop',
    }[derint]

    opts["contrl__runtyp"] = runtyp

    if method == "gamess":
        pass

    elif method in ["scf", "hf"]:
        pass

        # opts['contrl__mplevl'] = 0
        # opts['contrl__cityp'] = 'none'
        # opts['contrl__cctyp'] = 'none'

    elif method == "mp2":
        opts["contrl__mplevl"] = 2

    elif method == "lccd":
        opts["contrl__cctyp"] = "lccd"

    elif method == "ccd":
        opts["contrl__cctyp"] = "ccd"

    elif method == "ccsd":
        opts["contrl__cctyp"] = "ccsd"

    elif method in ["ccsd+t(ccsd)", "ccsd(t)"]:
        opts["contrl__cctyp"] = "ccsd(t)"

    elif method == "pbe":
        opts["contrl__dfttyp"] = "pbe"

    elif method == "b3lyp":
        opts["contrl__dfttyp"] = "b3lypv1r"

    elif method == "b3lyp5":
        opts["contrl__dfttyp"] = "b3lyp"

    else:
        raise InputError(f"Method not recognized: {method}")

    return opts
