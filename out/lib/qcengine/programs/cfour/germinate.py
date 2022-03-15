from typing import Any, Dict

from qcengine.exceptions import InputError


def muster_modelchem(method: str, derint: int) -> Dict[str, Any]:
    """Converts the QC method into CFOUR keywords."""

    method = method.lower()
    opts = {}

    if derint == 0:
        if method == "cfour":
            pass  # permit clean operation of sandwich mode
        else:
            opts["deriv_level"] = "zero"

    elif derint == 1:
        opts["deriv_level"] = "first"

    elif derint == 2:
        opts["vibration"] = "exact"

    if method == "cfour":
        pass

    elif method in ["scf", "hf"]:
        opts["calc_level"] = "scf"

    elif method == "mp2":
        opts["calc_level"] = "mp2"

    elif method == "mp3":
        opts["calc_level"] = "mp3"

    elif method == "mp4(sdq)":
        opts["calc_level"] = "sdq-mp4"

    elif method == "mp4":
        opts["calc_level"] = "mp4"

    elif method == "cc2":
        opts["calc_level"] = "cc2"

    elif method == "ccsd":
        opts["calc_level"] = "ccsd"

    elif method == "cc3":
        opts["calc_level"] = "cc3"

    elif method == "ccsd(t)":
        # Can't use (T) b/c bug in xsymcor lops it off
        opts["calc_level"] = "ccsd[t]"

    elif method == "ccsdt":
        opts["calc_level"] = "ccsdt"

    else:
        raise InputError(f"Method not recognized: {method}")

    return opts
