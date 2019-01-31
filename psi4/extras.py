#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
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
import shutil
import datetime

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

def exit_printing(start_time=None):
    """Prints the exit time and status.

    Parameters
    ----------
    start_time : datetime.datetime, optional
        starting time from which the elapsed time is computed.
 
    Returns
    -------
    None
 
    """
    end_time = datetime.datetime.now()
    core.print_out( "\n    Psi4 stopped on: {}".format(end_time.strftime('%A, %d %B %Y %I:%M%p')))
    if start_time is not None:
        run_time = end_time - start_time
        run_time = str(run_time).split('.')
        run_time = run_time[0] + '.' + run_time[1][:2]
        core.print_out( "\n    Psi4 wall time for execution: {}\n".format(run_time))
    if _success_flag_:
        core.print_out( "\n*** Psi4 exiting successfully. Buy a developer a beer!\n")
    else:
        core.print_out( "\n*** Psi4 encountered an error. Buy a developer more coffee!\n")
        core.print_out( "*** Resources and help at github.com/psi4/psi4.\n")

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


def _psi4_which(command, return_bool=False):
    # environment is $PSIPATH:$PATH, less any None values
    lenv = {'PATH': ':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) +
                    ':' + os.environ.get('PATH')}
    lenv = {k: v for k, v in lenv.items() if v is not None}

    ans = shutil.which(command, mode=os.F_OK | os.X_OK, path=lenv['PATH'])

    if return_bool:
        return bool(ans)
    else:
        return ans


def _plugin_import(plug):
    import sys
    if sys.version_info >= (3, 4):
        from importlib import util
        plug_spec = util.find_spec(plug)
    else:
        import pkgutil
        plug_spec = pkgutil.find_loader(plug)
    if plug_spec is None:
        return False
    else:
        return True


_addons_ = {
    "ambit": _CMake_to_Py_boolean("@ENABLE_ambit@"),
    "chemps2": _CMake_to_Py_boolean("@ENABLE_CheMPS2@"),
    "dkh": _CMake_to_Py_boolean("@ENABLE_dkh@"),
    "libefp": _plugin_import("pylibefp"),
    "erd": _CMake_to_Py_boolean("@ENABLE_erd@"),
    "gdma": _CMake_to_Py_boolean("@ENABLE_gdma@"),
    "pcmsolver": _CMake_to_Py_boolean("@ENABLE_PCMSolver@"),
    "simint": _CMake_to_Py_boolean("@ENABLE_simint@"),
    "dftd3": _psi4_which("dftd3", return_bool=True),
    "cfour": _psi4_which("xcfour", return_bool=True),
    "mrcc": _psi4_which("dmrcc", return_bool=True),
    "gcp": _psi4_which("gcp", return_bool=True),
    "v2rdm_casscf": _plugin_import("v2rdm_casscf"),
    "gpu_dfcc": _plugin_import("gpu_dfcc"),
    "forte": _plugin_import("forte"),
    "snsmp2": _plugin_import("snsmp2"),
    "resp": _plugin_import("resp"),
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
