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


# Figure out psidatadir: envvar trumps staged/installed
import os
psi4_module_loc = os.path.dirname(os.path.abspath(__file__))
pymod = os.path.normpath(os.path.sep.join(['@PYMOD_INSTALL_LIBDIR@', '@CMAKE_INSTALL_LIBDIR@', 'psi4']))
if pymod.startswith(os.path.sep + os.path.sep):
    pymod = pymod[1:]
pymod_dir_step = os.path.sep.join(['..'] * pymod.count(os.path.sep))
data_dir = os.path.sep.join([psi4_module_loc, pymod_dir_step, '@CMAKE_INSTALL_DATADIR@', 'psi4'])

# from . import config
# data_dir = config.psidatadir

if "PSIDATADIR" in os.environ.keys():
    data_dir = os.path.expanduser(os.environ["PSIDATADIR"])
elif "CMAKE_INSTALL_DATADIR" in data_dir:
    data_dir = os.path.sep.join([os.path.abspath(os.path.dirname(__file__)), "share", "psi4"])

data_dir = os.path.abspath(data_dir)
if not os.path.isdir(data_dir):
    raise KeyError("Unable to read the Psi4 Python folder - check the PSIDATADIR environmental variable"
                    "      Current value of PSIDATADIR is %s" % data_dir)
os.environ["PSIDATADIR"] = data_dir

# Init core
try:
    from . import core
except ImportError as err:
    if 'CXXABI' in str(err):
        raise ImportError("{0}\nLikely cause: GCC >= 4.9 not in [DY]LD_LIBRARY_PATH".format(err))
    else:
        raise ImportError("{0}".format(err))

from psi4.core import set_output_file, get_variable, set_variable, get_num_threads, set_num_threads
core.initialize()
core.efp_init()

if "PSI_SCRATCH" in os.environ.keys():
    envvar_scratch = os.environ["PSI_SCRATCH"]
    if not os.path.isdir(envvar_scratch):
        raise Exception("Passed in scratch is not a directory (%s)." % envvar_scratch)
    core.IOManager.shared_object().set_default_path(envvar_scratch)

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
from .extras import get_input_directory, addons, test

# Python portions of compiled-in Add-Ons
import sys
if "@ENABLE_PCMSolver@".upper() in ["1", "ON", "YES", "TRUE", "Y"]:
    sys.path.insert(1, "@PCMSolver_PYMOD@")

