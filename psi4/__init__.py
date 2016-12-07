#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2016 The Psi4 Developers.
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

# Figure out psidatadir: envvar trumps staged/installed
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

from psi4.core import set_output_file, set_variable  # noqa: F401
core.initialize()
core.efp_init()

# Cleanup core at exit
atexit.register(core.set_legacy_molecule, None)
atexit.register(core.clean)
atexit.register(core.finalize)


from .driver import *  # noqa: F401
from .header import print_header  # noqa: F401
from .metadata import __version__, version_formatter  # noqa: F401

# A few extraneous functions
from .extras import get_input_directory  # noqa: F401
