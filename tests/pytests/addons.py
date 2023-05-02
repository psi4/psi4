import os
import sys
from pathlib import Path
from typing import List

import pytest

from qcelemental.util import parse_version, which, which_import
from qcengine.testing import _programs as _programs_qcng

import psi4

__all__ = [
    "hardware_nvidia_gpu",
    "using",
    "uusing",
    "ctest_labeler",
    "ctest_runner",
]


def is_psi4_new_enough(version_feature_introduced):
    if not which_import('psi4', return_bool=True):
        return False
    import psi4
    return parse_version(psi4.__version__) >= parse_version(version_feature_introduced)


def is_qcfractal_new_enough(version_feature_introduced):
    if not which_import('qcfractal', return_bool=True):
        return False
    import qcfractal
    return parse_version(qcfractal.__version__) >= parse_version(version_feature_introduced)


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
    "ddx": which_import("pyddx", return_bool=True),
    "dkh": psi4.addons("dkh"),
    "ecpint": psi4.addons("ecpint"),
    "libefp": which_import("pylibefp", return_bool=True),
    "erd": psi4.addons("erd"),
    "fockci": which_import("psi4fockci", return_bool=True),  # package fockci, import psi4fockci
    "forte": which_import("forte", return_bool=True),
    "gdma": psi4.addons("gdma"),
    "gpu_dfcc": which_import("gpu_dfcc", return_bool=True),
    "ipi": which_import("ipi", return_bool=True),
    "mrcc": which("dmrcc", return_bool=True),
    "openfermionpsi4": which_import("openfermionpsi4", return_bool=True),
    "pcmsolver": psi4.addons("pcmsolver"),
    "psixas": which_import("psixas", return_bool=True),
    "resp": which_import("resp", return_bool=True),
    "simint": psi4.addons("simint"),
    "snsmp2": which_import("snsmp2", return_bool=True),
    "v2rdm_casscf": which_import("v2rdm_casscf", return_bool=True),
    "qcdb": False,  # capabilities of in-psi and out-of-psi qcdb not aligned
    "qcfractal": which_import("qcfractal", return_bool=True),
    "qcfractal_next": is_qcfractal_new_enough("0.49"),
    "qcportal": which_import("qcportal", return_bool=True),
    "bse": which_import("basis_set_exchange", return_bool=True),
}


def has_program(name):
    # Note for d3/gcp that EmpiricalDispersion is choosing between engines and which can't be globally controlled, so there's no use or truth in separate selectors
    if name == "dftd3":
        return _programs_qcng["s-dftd3"] or _programs_qcng["dftd3"]
    elif name == "gcp":
        return _programs_qcng["mctc-gcp"] or _programs_qcng["gcp"]
    elif name == "classic-dftd3":
        return _programs_qcng["dftd3"]
    elif name in _programs:
        return _programs[name]
    elif name in _programs_qcng:
        return _programs_qcng[name]
    else:
        raise KeyError(f"Program {name} not registered with Psi4 testing.")


_using_cache = {}


def _using(program: str) -> None:

    if program not in _using_cache:
        import_message = f"Not detecting module {program}. Install package if necessary to enable tests."
        skip = pytest.mark.skipif(has_program(program) is False, reason=import_message)
        general = pytest.mark.addon
        particular = getattr(pytest.mark, program)

        all_marks = (skip, general, particular)
        _using_cache[program] = [_compose_decos(all_marks), all_marks]


def _compose_decos(decos):
    # thanks, https://stackoverflow.com/a/45517876
    def composition(func):
        for deco in reversed(decos):
            func = deco(func)
        return func
    return composition


def uusing(program: str):
    """Apply 3 marks: skipif program not detected, label "addon", and label program.
    This is the decorator form for whole test functions: `@mark\n@mark`.

    """
    _using(program)
    return _using_cache[program][0]


def using(program: str) -> List:
    """Apply 3 marks: skipif program not detected, label "addon", and label program.
    This is the inline form for parameterizations: `marks=[]`.
    In combo, do `marks=[*using(), pytest.mark.quick]`

    """
    _using(program)
    return _using_cache[program][1]


def ctest_labeler(labels: str):
    """Apply each label in ``labels`` as PyTest marks. Also adds "psi" and "cli" marks.

    Parameters
    ----------
    labels
        A semicolon-separated list of labels to be applied as PyTest marks.
        These are usually copied from CMakeLists.txt.

    """
    if labels:
        marks = [getattr(pytest.mark, m) for m in labels.split(";")]
    else:
        marks = []
    all_marks = (*marks, pytest.mark.psi, pytest.mark.cli)
    return _compose_decos(all_marks)


hardware_nvidia_gpu = pytest.mark.skipif(
    True,  #is_nvidia_gpu_present() is False,
    reason='Psi4 not detecting Nvidia GPU via `nvidia-smi`. Install one')


def ctest_runner(inputdatloc, *, extra_infiles: List = None, outfiles: List = None, setenv: List = None):
    """Called from a mock PyTest function, this takes a full path ``inputdatloc`` to an ``"input.dat"`` file set up for
    CTest and submits it to the ``psi4`` executable. Any auxiliary files with names listed in ``extra_infiles`` that reside
    alongside ``inputdatloc`` are placed in the Psi4 execution directory. Any strings listed in ``setenv`` are available
    in the test as environment variables set to "1".

    """
    from qcengine.util import execute
    import psi4

    # Pass runtime env through to `execute`
    # * appending Psi4 import path (after all, it worked previous line) since partial/relative paths not robust
    psiimport = Path(psi4.__file__).parent.parent
    env = os.environ.copy()
    if "PYTHONPATH" in env:
        env["PYTHONPATH"] += os.pathsep + str(psiimport)
    else:
        env["PYTHONPATH"] = str(psiimport)

    if setenv:
        for ev in setenv:
            env[ev] = "1"

    ctestdir = Path(inputdatloc).resolve().parent

    if (ctestdir / "input.dat").exists():
        inputdat = "input.dat"
    elif (ctestdir / "input.py").exists():
        inputdat = "input.py"

    infiles = [inputdat]
    if extra_infiles:
        infiles.extend(extra_infiles)
    infiles_with_contents = {(Path(fl).name if Path(fl).is_absolute() else fl): (ctestdir / fl).read_text() for fl in infiles}

    # Note:  The simple `command = ["psi4", "input.dat"]` works fine for Linux and Mac but not for Windows.
    #   L/M/W   ok with `command = [which("psi4"), "input.dat"]` where `which` on Windows finds the psi4.bat file that points to the psi4 python script. -or-
    #   L/M/W   ok with `command = [sys.executable, psi4.executable, "input.dat"]` aka `python /full/path/bin/psi4 input.dat`.
    #   Latter chosen as `psi4.executable` is path computed by `import psi4`, so assured correspondence.
    # Note:  The input.py in json/, python/, and psi4numpy/ are not being treated best.
    #   Properly, as in CTest, it's `command = [sys.executable, "input.py"]`.
    #   Have to either have 3-item `command` or pass PYTHONPATH through env. Since some tests (fsapt) "import psi4" internally, doing both.
    #   The "exe" rerouting is for Windows conda package Scripts/psi4.exe file. It gives an encoding error if Python in command.
    # * Toggle `scratch_messy` to examine scratch contents.
    if Path(psi4.executable).suffix == ".exe":
        command = [psi4.executable, inputdat]
    else:
        command = [sys.executable, psi4.executable, inputdat]
    _, output = execute(command, infiles_with_contents, outfiles, environment=env, scratch_messy=False)

    success = output["proc"].poll() == 0
    assert success, output["stdout"] + output["stderr"]
