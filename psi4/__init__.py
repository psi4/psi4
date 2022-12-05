#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2022 The Psi4 Developers.
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


# Figure out psidatadir: envvar trumps staged/installed
import os
psi4_module_loc = os.path.dirname(os.path.abspath(__file__))
pymod = os.path.normpath(os.path.sep.join(['@PYMOD_INSTALL_LIBDIR@', '@CMAKE_INSTALL_LIBDIR@', 'psi4']))
if pymod.startswith(os.path.sep + os.path.sep):
    pymod = pymod[1:]
pymod_dir_step = os.path.sep.join(['..'] * pymod.count(os.path.sep))
data_dir = os.path.sep.join([psi4_module_loc, pymod_dir_step, '@CMAKE_INSTALL_DATADIR@', 'psi4'])
executable = os.path.abspath(os.path.sep.join([psi4_module_loc, pymod_dir_step, '@CMAKE_INSTALL_BINDIR@', 'psi4']))
from pathlib import Path
if not Path(executable).exists():
    # Win conda recipe moves psi4 executable unknown to CMake
    executable = str((Path(psi4_module_loc) / ".." / ".." / ".." / "Scripts" / "psi4.exe").resolve())

if "PSIDATADIR" in os.environ.keys():
    data_dir = os.path.expanduser(os.environ["PSIDATADIR"])
elif "CMAKE_INSTALL_DATADIR" in data_dir:
    data_dir = os.path.sep.join([os.path.abspath(os.path.dirname(__file__)), "share", "psi4"])

data_dir = os.path.abspath(data_dir)
if not os.path.isdir(data_dir):
    raise KeyError(f"Unable to read the Psi4 Python folder - check the PSIDATADIR environmental variable - current value is {data_dir}")

# Init core
from . import core

from psi4.core import get_num_threads, set_num_threads
core.initialize()

if "PSI_SCRATCH" in os.environ.keys():
    envvar_scratch = os.environ["PSI_SCRATCH"]
    if not os.path.isdir(envvar_scratch):
        raise Exception("Passed in scratch is not a directory (%s)." % envvar_scratch)
    core.IOManager.shared_object().set_default_path(envvar_scratch)

core.set_datadir(data_dir)
del psi4_module_loc, pymod, pymod_dir_step, data_dir

# Cleanup core at exit
import atexit
atexit.register(core.set_legacy_molecule, None)
atexit.register(core.clean_options)
atexit.register(core.clean)
atexit.register(core.finalize)

# Make official plugins accessible in input
from .driver import endorsed_plugins

# Manage threads. Must be after endorsed plugins, honestly.
core.set_num_threads(1, quiet=True)

# Load driver and outfile paraphernalia
from .driver import *
from .header import print_header
from .metadata import __version__, version_formatter

# A few extraneous functions
from .extras import get_input_directory, addons, test, set_output_file
from psi4.core import get_variable  # kill off in 1.4
from psi4.core import variable, set_variable

# Python portions of compiled-in Add-Ons
# * Note that this is a "battening down the hatches" for the many
#   rather than letting PYTHONPATH rule for the few.
import sys
if "@ENABLE_PCMSolver@".upper() in ["1", "ON", "YES", "TRUE", "Y"]:  # PCMSolver
    sys.path.insert(1, "@PCMSolver_PYMOD@")
if "@ENABLE_cppe@".upper() in ["1", "ON", "YES", "TRUE", "Y"]:  # cppe
    sys.path.insert(1, "@cppe_PYMOD@")
if "@ENABLE_ddx@".upper() in ["1", "ON", "YES", "TRUE", "Y"]:  # pyddx
    sys.path.insert(1, "@pyddx_PYMOD@")
if "@ENABLE_libefp@".upper() in ["1", "ON", "YES", "TRUE", "Y"]:  # pylibefp
    sys.path.insert(1, "@pylibefp_PYMOD@")

# Create a custom logger
import logging
root = logging.getLogger()  # create root
root.setLevel("ERROR")
logger = logging.getLogger(__name__)  # create initial psi4 from root to be configured later in extras
