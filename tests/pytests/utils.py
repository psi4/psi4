import os
import pprint
import sys

import qcengine as qcng

import psi4

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
            return False
        else:
            return True

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


def tnm():
    """Returns the name of the calling function, usually name of test case."""

    return sys._getframe().f_back.f_code.co_name


def run_psi4_cli(inputs, outputs, extra_commands=None, as_binary=None):
    """
    Runs Psi4 from the CLI in a subprocess.
    """

    if extra_commands is None:
        extra_commands = []
    psi_runner = os.path.join(os.path.dirname(os.path.abspath(psi4.__file__)), "run_psi4.py")
    cmds = [sys.executable, psi_runner, "--inplace"] + extra_commands + list(inputs.keys())

    success, ret = qcng.util.execute(cmds, inputs, outputs, as_binary=as_binary)
    return (success, ret)
