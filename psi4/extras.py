#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2021 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import os
import atexit
import datetime
from typing import Union

from qcelemental.util import which, which_import

from . import core

# Numpy place holder for files and cleanup
numpy_files = []


def register_numpy_file(filename):
    if not filename.endswith('.npy'): filename += '.npy'
    if filename not in numpy_files:
        numpy_files.append(filename)


def clean_numpy_files():
    for nfile in numpy_files:
        os.unlink(nfile)


atexit.register(clean_numpy_files)


def exit_printing(start_time=None, success=None):
    """Prints the exit time and status.

    Parameters
    ----------
    start_time : datetime.datetime, optional
        starting time from which the elapsed time is computed.
    success : bool
        Provides a success flag, otherwise uses the _success_flag_ global variable

    Returns
    -------
    None

    """
    end_time = datetime.datetime.now()
    core.print_out("\n    Psi4 stopped on: {}".format(end_time.strftime('%A, %d %B %Y %I:%M%p')))
    if start_time is not None:
        run_time = end_time - start_time
        run_time = str(run_time).split('.')
        run_time = run_time[0] + '.' + run_time[1][:2]
        core.print_out("\n    Psi4 wall time for execution: {}\n".format(run_time))

    if success is None:
        success = _success_flag_

    if success:
        core.print_out("\n*** Psi4 exiting successfully. Buy a developer a beer!\n")
    else:
        core.print_out("\n*** Psi4 encountered an error. Buy a developer more coffee!\n")
        core.print_out("*** Resources and help at github.com/psi4/psi4.\n")


_success_flag_ = False

# Working directory
_input_dir_ = os.getcwd()


def get_input_directory():
    return _input_dir_


# Add-Ons
def _CMake_to_Py_boolean(cmakevar):
    if cmakevar.upper() in ["1", "ON", "YES", "TRUE", "Y"]:
        return True
    else:
        return False


def psi4_which(command, *, return_bool: bool = False, raise_error: bool = False,
               raise_msg: str = None) -> Union[bool, None, str]:
    """Test to see if a command is available in Psi4 search path.

    Returns
    -------
    str or None
        By default, returns command path if command found or `None` if not.
        Environment is $PSIPATH:$PATH, less any None values.
    bool
        When `return_bool=True`, returns whether or not found.

    Raises
    ------
    ModuleNotFoundError
        When `raises_error=True` and command not found.

    """
    lenv = (os.pathsep.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(os.pathsep) if x != '']) +
            os.pathsep + os.environ.get('PATH', ''))

    return which(command=command, return_bool=return_bool, raise_error=raise_error, raise_msg=raise_msg, env=lenv)


_addons_ = {
    "ambit": _CMake_to_Py_boolean("@ENABLE_ambit@"),
    "chemps2": _CMake_to_Py_boolean("@ENABLE_CheMPS2@"),
    "dkh": _CMake_to_Py_boolean("@ENABLE_dkh@"),
    "libefp": which_import("pylibefp", return_bool=True),
    "erd": _CMake_to_Py_boolean("@ENABLE_erd@"),
    "gdma": _CMake_to_Py_boolean("@ENABLE_gdma@"),
    "ipi": which_import("ipi", return_bool=True),
    "pcmsolver": _CMake_to_Py_boolean("@ENABLE_PCMSolver@"),
    "cppe": which_import("cppe", return_bool=True),
    "simint": _CMake_to_Py_boolean("@ENABLE_simint@"),
    "dftd3": psi4_which("dftd3", return_bool=True),
    "cfour": psi4_which("xcfour", return_bool=True),
    "mrcc": psi4_which("dmrcc", return_bool=True),
    "gcp": psi4_which("gcp", return_bool=True),
    "v2rdm_casscf": which_import("v2rdm_casscf", return_bool=True),
    "gpu_dfcc": which_import("gpu_dfcc", return_bool=True),
    "forte": which_import("forte", return_bool=True),
    "snsmp2": which_import("snsmp2", return_bool=True),
    "resp": which_import("resp", return_bool=True),
    "psi4fockci": which_import("psi4fockci", return_bool=True),
    "adcc": which_import("adcc", return_bool=True),
    "mdi": which_import("mdi", return_bool=True),
    "cct3": which_import("cct3", return_bool=True),
}


def addons(request=None):
    """Returns boolean of whether Add-On *request* is available to Psi4,
    either compiled in or searchable in $PSIPATH:$PATH, as relevant. If
    *request* not passed, returns list of available Add-Ons.

    """
    if request is None:
        return sorted([k for k, v in _addons_.items() if v])
    return _addons_[request.lower()]


# Testing
def test(extent='full', extras=None):
    """Runs a test suite through pytest.

    Parameters
    ----------
    extent : {'smoke', 'quick', 'full', 'long'}
        All choices are defined, but choices may be redundant in some projects.
        _smoke_ will be minimal "is-working?" test(s).
        _quick_ will be as much coverage as can be got quickly, approx. 1/3 tests.
        _full_ will be the whole test suite, less some exceedingly long outliers.
        _long_ will be the whole test suite.
    extras : list
        Additional arguments to pass to `pytest`.

    Returns
    -------
    int
        Return code from `pytest.main()`. 0 for pass, 1 for fail.

    """
    try:
        import pytest
    except ImportError:
        raise RuntimeError('Testing module `pytest` is not installed. Run `conda install pytest`')
    abs_test_dir = os.path.sep.join([os.path.abspath(os.path.dirname(__file__)), "tests"])

    command = ['-rws', '-v']
    if extent.lower() == 'smoke':
        command.extend(['-m', 'smoke'])
    elif extent.lower() == 'quick':
        command.extend(['-m', 'quick or smoke'])
    elif extent.lower() == 'full':
        command.extend(['-m', 'not long'])
    elif extent.lower() == 'long':
        pass
    if extras is not None:
        command.extend(extras)
    command.extend(['--capture=sys', abs_test_dir])

    retcode = pytest.main(command)
    return retcode

