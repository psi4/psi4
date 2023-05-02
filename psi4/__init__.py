#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2023 The Psi4 Developers.
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


# Figure out paths
# * in figuring out psidatadir: envvar trumps staged/installed
# * some full paths are computed here using the prefix, but all outputs are relative to __file__, so relocatability preserved
# * note that all path entities are directories except for "executable" that is a file
import os
from pathlib import Path
psi4_module_loc = Path(__file__).resolve().parent

prefix = Path("@CMAKE_INSTALL_PREFIX@".replace("\\", "/"))
cmake_install_bindir = "@CMAKE_INSTALL_BINDIR@".replace("\\", "/")
cmake_install_datadir = "@CMAKE_INSTALL_DATADIR@".replace("\\", "/")
cmake_install_libdir = "@CMAKE_INSTALL_LIBDIR@".replace("\\", "/")
pymod_install_libdir = "@PYMOD_INSTALL_LIBDIR@".lstrip("/")
full_pymod = (prefix / cmake_install_libdir / pymod_install_libdir / "psi4").resolve()
full_data = prefix / cmake_install_datadir / "psi4"
full_bin = prefix / cmake_install_bindir
rel_data = os.path.relpath(full_data, start=full_pymod)
rel_bin = os.path.relpath(full_bin, start=full_pymod)

executable = psi4_module_loc.joinpath(rel_bin, "psi4")
executable_exe = (Path("/opt/anaconda1anaconda2anaconda3") / "Scripts" / "psi4.exe").resolve(strict=False)
if executable_exe.exists():
    # Win conda-build generates this unbeknownst to CMake
    executable = executable_exe
executable = str(executable.resolve())

data_dir = psi4_module_loc.joinpath(rel_data)
if "PSIDATADIR" in os.environ.keys():
    data_dir = Path(os.path.expanduser(os.environ["PSIDATADIR"]))
elif "CMAKE_INSTALL_DATADIR" in str(data_dir):
    data_dir = Path(os.path.sep.join([os.path.abspath(os.path.dirname(__file__)), "share", "psi4"]))

data_dir = data_dir.resolve(strict=False)
if not data_dir.is_dir():
    raise KeyError(f"Unable to read the Psi4 Python folder - check the PSIDATADIR environmental variable - current value is {str(data_dir)}")
data_dir = str(data_dir)

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
del cmake_install_bindir, cmake_install_datadir, cmake_install_libdir, pymod_install_libdir
del psi4_module_loc, prefix, full_pymod, full_data, full_bin, rel_data, rel_bin, data_dir, executable_exe

# Cleanup core at exit
import atexit
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
logger = logging.getLogger(__name__)  # create initial psi4 from root to be configured later in extras
