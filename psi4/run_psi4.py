#! /usr/bin/env python

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
import argparse

parser = argparse.ArgumentParser(description="A hybrid C++/Python quantum chemistry module.")
parser.add_argument("-i", "--input", default="input.dat", help="Input file name. Default input.dat.")
parser.add_argument("-o", "--output", help="Redirect output elsewhere.\n"
                                           "Default: when input filename is 'input.dat', then defaults to 'output.dat'. "
                                           "Otherwise, output filename defaults to input filename with any '.in' or 'dat' extension replaced by '.out'\n")

parser.add_argument("-v", "--verbose", action='store_true', help="Print a lot of information.")
parser.add_argument("-V", "--version", action='store_true', help="Print version information.")

# parser.add_argument("-d", "--debug", action='store_true', help="Flush the outfile at every print statement.")
parser.add_argument("-k", "--skip-preprocessor", action='store_true', help="Skips input preprocessing. Expert mode.")
parser.add_argument("-m", "--messy", action='store_true', help="Leave temporary files after the run is completed.")
#parser.add_argument("-r", "--restart", action='store_true', help="Number to be used instead of process id.")

parser.add_argument("-s", "--scratch", help="Psi4 scratch directory to use.")
parser.add_argument("-a", "--append", action='store_true', help="Append results to output file. Default Truncate first")
parser.add_argument("-l", "--psidatadir", help="Specify where to look for the Psi data directory. Overrides PSIDATADIR.")
parser.add_argument("-n", "--nthread", default=1, help="Number of threads to use (overrides OMP_NUM_THREADS)")
parser.add_argument("-p", "--prefix", help="Prefix name for psi files. Default psi")
parser.add_argument("--inplace", action='store_true', help="Runs psi4 from the source directory. "
                                                           "!Warning! expert option.")

# For plugins
parser.add_argument("--new-plugin", action='store_true',
                    help="Creates a new directory with files for writing a "
                    "new plugin. You can specify an additional argument "
                    "that specifies a template to use, for example "
                    "--new-plugin name +mointegrals")
parser.add_argument("--new-plugin-makefile", action='store_true',
                    help="Creates Makefile that can be used to compile"
                    "plugins. The Makefile is placed in the current"
                    "directory.")

# print("Environment Variables\n");
# print("     PSI_SCRATCH           Directory where scratch files are written.")
# print("                           Default: $TMPDIR (or /tmp/ when not set)")
# print("                           This should be a local, not network, disk")

#parser.print_help()
args, unknown = parser.parse_known_args()
args = args.__dict__ # Namespace object seems silly

# Deal with kwargs not implemented
if args["new_plugin"] or args["new_plugin_makefile"]:
    raise Exception("Plugins are not currently implemented")


# Figure out pythonpath
cmake_install_prefix = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + '..'
lib_dir = os.path.sep.join([cmake_install_prefix, "@CMAKE_INSTALL_LIBDIR@", "@PYMOD_INSTALL_LIBDIR@"])

if args["inplace"]:
    if "CMAKE_INSTALL_LIBDIR" not in lib_dir:
        raise ImportError("Cannot run inplace from a installed directory.")

    core_location = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + "core.so"
    if not os.path.isfile(core_location):
        raise ImportError("A compiled Psi4 core.so needs to be symlinked to the repository/psi4 folder")

    lib_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if ("PSIDATADIR" not in os.environ.keys()) and (not args["psidatadir"]):
        data_dir = os.path.sep.join([os.path.abspath(os.path.dirname(__file__)), "share", "psi4"])
        os.environ["PSIDATADIR"] = data_dir


elif "CMAKE_INSTALL_LIBDIR" in lib_dir:
    raise ImportError("Psi4 was not installed correctly!")

# Replace input/output if unknown kwargs
if len(unknown) > 0:
    args["input"] = unknown[0]
if len(unknown) > 1:
    args["output"] = unknown[1]
if len(unknown) > 2:
    raise KeyError("Too many unknown arguments: %s" % str(unknown))

# Figure out output arg
if args["output"] is None:
    if args["input"] == "input.dat":
        args["output"] = "output.dat"
    elif args["input"].endswith(".in"):
        args["output"] = args["input"].replace(".in", ".out")
    elif args["input"].endswith(".dat"):
        args["output"] = args["input"].replace(".dat", ".out")
    else:
        args["output"] = args["input"] + ".dat"

# Transmit any argument psidatadir through environ
if args["psidatadir"] is not None:
    data_dir = os.path.abspath(os.path.expanduser(args["psidatadir"]))
    os.environ["PSIDATADIR"] = data_dir


### Actually import psi4 and apply setup ###

# Import installed psi4
sys.path.insert(1, lib_dir)
import psi4

if args["version"]:
    print(psi4.__version__)
    sys.exit()

if not os.path.isfile(args["input"]):
    raise KeyError("The file %s does not exist." % args["input"])

# Setup outfile
if args["append"] is None:
    args["append"] = False
if args["output"] != "stdout":
    psi4.core.set_output_file(args["output"], args["append"])

# Set a few options
if args["prefix"] is not None:
    psi4.core.set_psi_file_prefix(args["prefix"])
psi4.core.set_nthread(int(args["nthread"]))
psi4.core.set_memory(524288000, True)
psi4.extras._input_dir_ = os.path.dirname(os.path.abspath(args["input"]))
psi4.print_header()

# Read input
with open(args["input"]) as f:
    content = f.read()

if args["scratch"] is not None:
    if not os.path.isdir(args["scratch"]):
        raise Exception("Passed in scratch is not a directory (%s)." % args["scratch"])
    psi4.core.set_environment("PSI_SCRATCH", args["scratch"])

# Preprocess
if not args["skip_preprocessor"]:
    # PSI_SCRATCH must be set before this call!
    content = psi4.process_input(content)

# Handle Verbose
if args["verbose"]:
    psi4.core.print_out('\nParsed Psithon:')
    psi4.core.print_out(content)
    psi4.core.print_out('-' * 75)

# Handle Messy
if args["messy"]:
    import atexit
    for handler in atexit._exithandlers:
        if handler[0] == psi4.core.clean:
            atexit._exithandlers.remove(handler)

# Only register exit if exec was successful
import atexit
atexit.register(psi4.extras.exit_printing)

# Run the program!
exec(content)

#    elif '***HDF5 library version mismatched error***' in str(err):
#        raise ImportError("{0}\nLikely cause: HDF5 used in compilation not prominent enough in RPATH/[DY]LD_LIBRARY_PATH".format(err))

# Celebrate with beer
psi4.extras._success_flag_ = True

