"""
Misc information and runtime information.
"""

import re

from . import _version

__all__ = ["get_information", "provenance_stamp"]

versions = _version.get_versions()

__info = {"version": versions["version"], "git_revision": versions["full-revisionid"]}


def get_information(key):
    """
    Obtains a variety of runtime information about QCEngine.
    """
    key = key.lower()
    if key not in __info:
        raise KeyError("Information key '{}' not understood.".format(key))

    return __info[key]


def provenance_stamp(routine):
    """Return dictionary satisfying QCSchema,
    https://github.com/MolSSI/QCSchema/blob/master/qcschema/dev/definitions.py#L23-L41
    with QCEngine's credentials for creator and version. The
    generating routine's name is passed in through `routine`.

    """
    return {"creator": "QCEngine", "version": get_information("version"), "routine": routine}


_yes = re.compile(r"^(yes|true|on|1)", re.IGNORECASE)
_no = re.compile(r"^(no|false|off|0)", re.IGNORECASE)
_der0th = re.compile(r"^(0|none|energy)", re.IGNORECASE)
_der1st = re.compile(r"^(1|first|gradient)", re.IGNORECASE)
_der2nd = re.compile(r"^(2|second|hessian)", re.IGNORECASE)
_der3rd = re.compile(r"^(3|third)", re.IGNORECASE)
_der4th = re.compile(r"^(4|fourth)", re.IGNORECASE)
_der5th = re.compile(r"^(5|fifth)", re.IGNORECASE)
