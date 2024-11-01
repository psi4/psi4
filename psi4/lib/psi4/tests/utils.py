import os
import pprint
import sys
import shutil

import qcengine as qcng

import psi4
import pytest

pp = pprint.PrettyPrinter(width=120)


__all__ = [
    'a2a',
    'compare',
    'compare_integers',
    'compare_strings',
    'compare_values',
    'compare_arrays',
    'compare_recursive',
    'compare_molrecs',
    'compare_cubes',
    'compare_vectors',
    'compare_matrices',
    'compare_wavefunctions',
    'compare_fcidumps',
    'compare_fchkfiles',
    'run_psi4_cli',
    'tnm',
]

# CODATA ratio 2014 / 2010 Bohr to Angstroms conversion factor
a2a = 0.52917721067 / 0.52917720859


def true_false_decorator(compare_fn, *args, **kwargs):
    """Turns `compare_fn` that raises `psi4.TestComparisonError` on failure into a function that
    returns True on success and False on failure, suitable for assertions in pytest.

    """

    def true_false_wrapper(*args, **kwargs):
        try:
            compare_fn(*args, **kwargs)
        except psi4.TestComparisonError as err:
            return (False, err) if "return_message" in kwargs else False
        else:
            return (True, "") if "return_message" in kwargs else True

    return true_false_wrapper


compare = true_false_decorator(psi4.compare)
compare_integers = true_false_decorator(psi4.compare_integers)
compare_strings = true_false_decorator(psi4.compare_strings)
compare_values = true_false_decorator(psi4.compare_values)
compare_arrays = true_false_decorator(psi4.compare_arrays)

compare_recursive = true_false_decorator(psi4.compare_recursive)
compare_molrecs = true_false_decorator(psi4.compare_molrecs)

compare_cubes = true_false_decorator(psi4.compare_cubes)
compare_vectors = true_false_decorator(psi4.compare_vectors)
compare_matrices = true_false_decorator(psi4.compare_matrices)
compare_fcidumps = true_false_decorator(psi4.compare_fcidumps)
compare_wavefunctions = true_false_decorator(psi4.compare_wavefunctions)
compare_fchkfiles = true_false_decorator(psi4.compare_fchkfiles)

def tnm():
    """Returns the name of the calling function, usually name of test case."""

    return sys._getframe().f_back.f_code.co_name


def run_psi4_cli(inputs, outputs, extra_commands=None, as_binary=None):
    """
    Runs Psi4 from the CLI in a subprocess.
    """

    if extra_commands is None:
        extra_commands = []

    cmds = None
    psidir = os.path.dirname(os.path.abspath(psi4.__file__))

    # Check inplace
    psi_runner = os.path.join(psidir, "run_psi4.py")
    if os.path.isfile(psi_runner):
        cmds = [sys.executable, psi_runner, "--inplace"]

    # Check test install
    if cmds is None:
        binpath = os.path.join(os.path.dirname(os.path.dirname(psidir)), "bin")
        psi_bin = shutil.which('psi4', path=binpath)
        if psi_bin:
            cmds = [psi_bin]

    # Check if Psi4 in path
    if cmds is None:
        pytest.skip("Could not find Psi4 executable.")

    cmds = cmds + extra_commands + list(inputs.keys())

    success, ret = qcng.util.execute(cmds, inputs, outputs, as_binary=as_binary)
    return (success, ret)
