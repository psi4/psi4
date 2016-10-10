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

import sys
import os

# Figure out psidatadir: envvar trumps staged/installed
psi4_module_loc = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.sep.join([psi4_module_loc, '..', '..', '@CMAKE_INSTALL_DATADIR@', 'psi4'])

if "PSIDATADIR" in os.environ.keys():
    data_dir = os.path.expanduser(os.environ["PSIDATADIR"])

data_dir = os.path.abspath(data_dir)
if not os.path.isdir(data_dir):
    raise KeyError("Unable to read the Psi4 Python folder - check the PSIDATADIR environmental variable"
                    "      Current value of PSIDATADIR is %s" % data_dir)
os.environ["PSIDATADIR"] = data_dir

# Find and import the core
try:
    from . import core
except ImportError as err:
    if 'CXXABI' in str(err):
        raise ImportError("{0}\nLikely cause: GCC >= 4.9 not in [DY]LD_LIBRARY_PATH".format(err))

    # Check if we are running in place
    check_inplace_file = os.path.abspath(os.path.dirname(__file__)) + os.path.sep + "run_psi4.py.in"
    if not os.path.isfile(check_inplace_file):
        raise ImportError("{0}".format(err))

    print("\nRunning psi4 from the source directory, looking for the 'objdir' build directory for core.so ...")
    base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    obj_dir = os.environ["CORE_PSI4_YO"]
    base_path += os.path.sep + obj_dir + os.path.sep + 'stage'
    matches = []
    for root, dirnames, filenames in os.walk(base_path):
        if 'include' in root: continue
        if 'share' in root: continue
        if 'core.so' in filenames:
            matches.append(root)

    if len(matches) == 0:
        raise ImportError("Could not find core.so in basepath: %s" % base_path)

    print("Found core.so at %s\n" % matches[0])
    sys.path.insert(1, matches[0])
    import core

# Init core
core.initialize()
core.efp_init()

# Cleanup core at exit
import atexit
atexit.register(core.set_legacy_molecule, None)
atexit.register(core.clean)
atexit.register(core.finalize)

# Numpy place holder for files and cleanup
numpy_files = []
def register_numpy_file(filename):
    if filename not in numpy_files:
        numpy_files.append(filename)

def clean_numpy_files():
    for nfile in numpy_files:
        os.unlink(nfile)

atexit.register(clean_numpy_files)
