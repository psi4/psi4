#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import os
import atexit
import subprocess

from . import core

# Numpy place holder for files and cleanup
numpy_files = []
def register_numpy_file(filename):
    if filename not in numpy_files:
        numpy_files.append(filename)

def clean_numpy_files():
    for nfile in numpy_files:
        os.unlink(nfile)

atexit.register(clean_numpy_files)

# Exit printing
def exit_printing():
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


def _psi4_which(command):
    # environment is $PSIPATH:$PATH, less any None values
    lenv = {'PATH': ':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) +
                    ':' + os.environ.get('PATH')}
    lenv = {k: v for k, v in lenv.items() if v is not None}

    # thanks, http://stackoverflow.com/a/11270665
    try:
        from subprocess import DEVNULL  # py33
    except ImportError:
        DEVNULL = open(os.devnull, 'wb')

    try:
        dashout = subprocess.Popen(command, stdout=DEVNULL, stderr=subprocess.STDOUT, env=lenv)
    except OSError as e:
        return False
    else:
        return True


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
    "libefp": _CMake_to_Py_boolean("@ENABLE_libefp@"),
    "erd": _CMake_to_Py_boolean("@ENABLE_erd@"),
    "gdma": _CMake_to_Py_boolean("@ENABLE_gdma@"),
    "pcmsolver": _CMake_to_Py_boolean("@ENABLE_PCMSolver@"),
    "simint": _CMake_to_Py_boolean("@ENABLE_simint@"),
    "dftd3": _psi4_which("dftd3"),
    "cfour": _psi4_which("xcfour"),
    "mrcc": _psi4_which("dmrcc"),
    "gcp": _psi4_which("gcp"),
    "v2rdm_casscf": _plugin_import("v2rdm_casscf"),
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
def test():
    """Runs a smoke test suite through pytest."""

    try:
        import pytest
    except ImportError:
        raise RuntimeError('Testing module `pytest` is not installed. Run `conda install pytest`')
    abs_test_dir = os.path.sep.join([os.path.abspath(os.path.dirname(__file__)), "tests"])
    pytest.main(['-v', abs_test_dir])
