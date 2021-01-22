import pytest

from qcelemental.util import parse_version, which, which_import
from qcengine.testing import _programs as _programs_qcng

import psi4

__all__ = [
    "hardware_nvidia_gpu",
    "using",
]


def is_psi4_new_enough(version_feature_introduced):
    if not which_import('psi4', return_bool=True):
        return False
    import psi4
    return parse_version(psi4.__version__) >= parse_version(version_feature_introduced)


def is_nvidia_gpu_present():
    try:
        import GPUtil
    except ModuleNotFoundError:
        try:
            import gpu_dfcc
        except ModuleNotFoundError:
            # who knows?
            return False
        else:
            return gpu_dfcc.cudaGetDeviceCount() > 0
    else:
        try:
            ngpu = len(GPUtil.getGPUs())
        except FileNotFoundError:
            # no `nvidia-smi`
            return False
        else:
            return ngpu > 0


# Figure out what is imported
# * only using psi4.addons for C-linked b/c _programs more responsive to env
# * `which` not taking PSIPATH into account
_programs = {
    # non-QC
    "memory_profiler": which_import('memory_profiler', return_bool=True),
    "networkx": which_import("networkx", return_bool=True),

    # QC
    "adcc": which_import("adcc", return_bool=True),
    "ambit": psi4.addons("ambit"),
    "cct3": which_import("cct3", return_bool=True),
    "chemps2": psi4.addons("chemps2"),
    "cppe": which_import("cppe", return_bool=True),  # package pycppe, import cppe
    "dkh": psi4.addons("dkh"),
    "libefp": which_import("pylibefp", return_bool=True),
    "erd": psi4.addons("erd"),
    "fockci": which_import("psi4fockci", return_bool=True),  # package fockci, import psi4fockci
    "forte": which_import("forte", return_bool=True),
    "gdma": psi4.addons("gdma"),
    "gpu_dfcc": which_import("gpu_dfcc", return_bool=True),
    "ipi": which_import("ipi", return_bool=True),
    "mrcc": which("dmrcc", return_bool=True),
    "pcmsolver": psi4.addons("pcmsolver"),
    "psixas": which_import("psixas", return_bool=True),
    "resp": which_import("resp", return_bool=True),
    "simint": psi4.addons("simint"),
    "snsmp2": which_import("snsmp2", return_bool=True),
    "v2rdm_casscf": which_import("v2rdm_casscf", return_bool=True),
}


def has_program(name):
    if name in _programs:
        return _programs[name]
    elif name in _programs_qcng:
        return _programs_qcng[name]
    else:
        raise KeyError(f"Program {name} not registered with Psi4 testing.")


_using_cache = {}


def using(program):

    if program not in _using_cache:
        import_message = f"Not detecting module {program}. Install package if necessary to enable tests."
        skip = pytest.mark.skipif(has_program(program) is False, reason=import_message)
        _using_cache[program] = skip

    return _using_cache[program]


hardware_nvidia_gpu = pytest.mark.skipif(
    True,  #is_nvidia_gpu_present() is False,
    reason='Psi4 not detecting Nvidia GPU via `nvidia-smi`. Install one')
