import math
import re

import numpy as np
from qcelemental.exceptions import ValidationError
from qcelemental.util import filter_comments


def load_hessian(shess: str, dtype: str) -> np.ndarray:
    """Construct a Hessian array from any recognized string format.

    Parameters
    ----------
    shess
        Multiline string specification of Hessian in a recognized format.
    dtype
        {"fcmfinal", "cfour", "gamess"}
        Hessian format name.

    Returns
    -------
    np.ndarray
        Hessian array in square shape.

    """
    # list o'lines w/o comments or blanks
    shess = filter_comments(shess)
    lhess = list(filter(None, map(str.strip, shess.splitlines())))

    if dtype in ["fcmfinal", "cfour"]:
        nat = int(lhess[0].split()[0])
        ndof = 3 * nat
        datastr = "\n".join(lhess[1:])
        nhess = np.fromstring(datastr, sep=" ")
        nhess = nhess.reshape(ndof, ndof)
    elif dtype == "gamess":
        if "ENERGY" in lhess[0]:
            lhess.pop(0)
        datastr = []
        for ln in lhess:
            numbers = re.findall(r"([-+]?\d+\.\d+[DdEe][-+]\d\d)", ln)
            if numbers:
                datastr.extend(numbers)

        nhess = np.fromstring(" ".join(datastr), sep=" ")
        ndof = int(math.sqrt(len(nhess)))
        nhess = nhess.reshape((ndof, ndof))
    else:
        raise ValidationError("Unknown dtype: {}".format(dtype))

    return nhess


def hess_to_string(hess, handle, dtype):
    nat = hess.shape[0] // 3
    assert hess.shape == (3 * nat, 3 * nat)

    header = "{:5}{:5}".format(nat, 6 * nat)
    np.savetxt(handle, hess.reshape((-1, 3)), fmt="%20.10f", delimiter="", newline="\n", header=header, comments="")
