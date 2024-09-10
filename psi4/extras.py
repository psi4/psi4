#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2024 The Psi4 Developers.
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

import atexit
import datetime
import itertools
import os
from pathlib import Path
from typing import List, Optional, Union

from qcelemental.util import which, which_import

from . import core
from .header import print_header as _print_header

# Numpy place holder for files and cleanup
numpy_files = []


def register_numpy_file(filename):
    if not filename.endswith('.npy'): filename += '.npy'
    if filename not in numpy_files:
        numpy_files.append(filename)


def register_scratch_file(filename):
    if filename not in numpy_files:
        numpy_files.append(filename)


def clean_numpy_files():
    for nfile in numpy_files:
        try:
            os.unlink(nfile)
        except OSError:
            pass


atexit.register(clean_numpy_files)


def exit_printing(start_time: datetime.datetime = None, success: bool = None) -> None:
    """Prints the exit time and status.

    Parameters
    ----------
    start_time
         starting time from which the elapsed time is computed.
    success
        Provides a success flag, otherwise uses the ``_success_flag_`` global variable

    Returns
    -------
    None

    """
    end_time = datetime.datetime.now()
    core.print_out("\n    Psi4 stopped on: {}".format(end_time.strftime('%A, %d %B %Y %I:%M%p')))
    if start_time is not None:
        run_time = end_time - start_time
        run_time = str(run_time).split('.')
        # python "helpfully" truncates microseconds if zero. Undo that.
        if len(run_time) == 1: run_time.append("000000")
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
    "adcc": which_import("adcc", return_bool=True),
    "ambit": _CMake_to_Py_boolean("@ENABLE_ambit@"),
    "chemps2": _CMake_to_Py_boolean("@ENABLE_CheMPS2@"),
    "dkh": _CMake_to_Py_boolean("@ENABLE_dkh@"),
    "ecpint": _CMake_to_Py_boolean("@ENABLE_ecpint@"),
    "libefp": which_import("pylibefp", return_bool=True),
    "gdma": which_import("gdma", return_bool=True),  # package pygdma, import gdma
    "ipi": which_import("ipi", return_bool=True),
    "pcmsolver": _CMake_to_Py_boolean("@ENABLE_PCMSolver@"),
    "cppe": which_import("cppe", return_bool=True),
    "ddx": which_import("pyddx", return_bool=True),
    "simint": _CMake_to_Py_boolean("@ENABLE_simint@"),
    "dftd3": which_import("dftd3", return_bool=True),
    "cfour": psi4_which("xcfour", return_bool=True),
    "mrcc": psi4_which("dmrcc", return_bool=True),
    "gcp": psi4_which("mctc-gcp", return_bool=True),
    "v2rdm_casscf": which_import("v2rdm_casscf", return_bool=True),
    "gpu_dfcc": which_import("gpu_dfcc", return_bool=True),
    "forte": which_import("forte", return_bool=True),
    "snsmp2": which_import("snsmp2", return_bool=True),
    "resp": which_import("resp", return_bool=True),
    "psi4fockci": which_import("psi4fockci", return_bool=True),
    "mdi": which_import("mdi", return_bool=True),
    "cct3": which_import("cct3", return_bool=True),
    "dftd4": which_import("dftd4", return_bool=True),
    "mp2d": psi4_which("mp2d", return_bool=True),
    "openfermionpsi4": which_import("openfermionpsi4", return_bool=True),
    "geometric": which_import("geometric", return_bool=True),
    #"optking": which_import("optking", return_bool=True),
    "psixas": which_import("psixas", return_bool=True),
    #"mctc-gcp": psi4_which("mctc-gcp", return_bool=True),
    "bse": which_import("basis_set_exchange", return_bool=True),
    "einsums": _CMake_to_Py_boolean("@ENABLE_Einsums@"),
    "gauxc": _CMake_to_Py_boolean("@ENABLE_gauxc@"),
}


def addons(request: str = None) -> Union[bool, List[str]]:
    """Returns boolean of whether Add-On *request* is available to Psi4,
    either compiled in or searchable in $PSIPATH:$PATH, as relevant. If
    *request* not passed, returns list of available Add-Ons: `['adcc', 'ambit', 'c̶c̶t̶3̶', ...` .

    """
    def strike(text):
        if os.name == "nt":
            # Windows has a probably correctable problem with unicode, but I can't iterate it quickly, so use tilde for strike.
            #   UnicodeEncodeError: 'charmap' codec can't encode character '\u0336' in position 3: character maps to <undefined>
            return "~" + text + "~"
        else:
            return ''.join(itertools.chain.from_iterable(zip(text, itertools.repeat('\u0336'))))

    if request is None:
        return [(k if v else strike(k)) for k, v in sorted(_addons_.items())]
    return _addons_[request.lower()]


# Testing
def test(extent: str = "full", extras: List = None) -> int:
    """Runs a test suite through pytest.

    Parameters
    ----------
    extent
        {'smoke', 'quick', 'full', 'long'}
        All choices are defined, but choices may be redundant in some projects.

          * _smoke_ will be minimal "is-working?" test(s).
          * _quick_ will be as much coverage as can be got quickly, approx. 1/3 tests.
          * _full_ will be the whole test suite, less some exceedingly long outliers.
          * _long_ will be the whole test suite.
    extras
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

    command = ['-rws', '-v', '--color', 'yes']
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


def set_output_file(
    ofile: str,
    append: bool = False,
    *,
    loglevel: int = 20,
    execute: bool = True,
    print_header: Optional[bool] = None,
    inherit_loglevel: bool = False) -> Path:
    """Set the name for output and logging files.

    Parameters
    ----------
    ofile
        Name of ASCII output file including extension. The logging file is set from this string with a ``.log`` extension.
    append
        Do append to the output and logging files rather than (the default) truncating them?
    loglevel
        The criticality level at which to log. 30 for WARN (Python default), 20 for INFO, 10 for DEBUG
    execute
        Do set ``ofile`` via :py:func:`psi4.core.set_output_file` and add the logger, rather than just returning ``ofile`` path.
    print_header
        Whether to write the Psi4 header to the ASCII output file. (Only applicable if ``execute=True``.) By default,
        writes if file is truncated (``append=False``) but not if appended.
    inherit_loglevel
        If true, do not set loglevel even to default value. Instead, allow level to be inherited from existing logger.

    Returns
    -------
    ~pathlib.Path
        ``Path(ofile)``

    Notes
    -----
    This :py:func:`psi4.set_output_file` command calls :py:func:`psi4.core.set_output_file` and should be used in
    preference to it as this additionally sets up logging.

    """
    out = Path(ofile)
    log = out.with_suffix(".log")

    # Get the custom logger
    import logging

    from psi4 import logger
    if not inherit_loglevel:
        logger.setLevel(loglevel)

    # Create formatters
    # * detailed:  example: 2019-11-20:01:13:46,811 DEBUG    [psi4.driver.task_base:156]
    f_format_detailed = logging.Formatter("%(asctime)s,%(msecs)d %(levelname)-8s [%(name)s:%(lineno)d] %(message)s", datefmt="%Y-%m-%d:%H:%M:%S")
    # * light:     example: 2019-11-20:10:45:21 FINDIFREC CLASS INIT DATA
    f_format_light = logging.Formatter("%(asctime)s %(message)s", datefmt="%Y-%m-%d:%H:%M:%S")

    # Create handlers, add formatters to handlers, and add handlers to logger (StreamHandler() also available)
    filemode = "a" if append else "w"
    f_handler = logging.FileHandler(log, filemode)
    f_handler.setLevel(logging.DEBUG)
    f_handler.setFormatter(f_format_detailed)

    if execute:
        core.set_output_file(str(out), append)
        if print_header is True or (print_header is None and not append):
            _print_header()
        # Warning: baseFilename is not part of the documented API for the logging module and could change.
        filenames = [handle.baseFilename for handle in logger.handlers]
        if not f_handler.baseFilename in filenames:
            logger.addHandler(f_handler)
    return out
