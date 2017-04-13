#!@PYTHON_EXECUTABLE@

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

import sys
import os
import json
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description="Psi4: Open-Source Quantum Chemistry", formatter_class=RawTextHelpFormatter)
parser.add_argument("-i", "--input", default="input.dat",
                    help="Input file name. Default: input.dat.")
parser.add_argument("-o", "--output", help="""\
Redirect output elsewhere.
Default: when input filename is 'input.dat', 'output.dat'.
Otherwise, output filename defaults to input filename with
any '.in' or 'dat' extension replaced by '.out'""")
parser.add_argument("-a", "--append", action='store_true',
                    help="Appends results to output file. Default: Truncate first")
parser.add_argument("-V", "--version", action='store_true',
                    help="Prints version information.")
parser.add_argument("-n", "--nthread", default=1,
                    help="Number of threads to use. Psi4 disregards OMP_NUM_THREADS/MKL_NUM_THREADS.")
parser.add_argument("-s", "--scratch",
                    help="Scratch directory to use. Overrides PSI_SCRATCH.")
parser.add_argument("-m", "--messy", action='store_true',
                    help="Leaves temporary files after the run is completed.")
# parser.add_argument("-d", "--debug", action='store_true', help="Flush the outfile at every print statement.")
# parser.add_argument("-r", "--restart", action='store_true', help="Number to be used instead of process id.")
parser.add_argument("-p", "--prefix",
                    help="Prefix name for psi files. Default psi")
parser.add_argument("--psiapi-path", action='store_true',
                    help="""Generates a bash command to source correct Python """
                         """interpreter and path for ``python -c "import psi4"``""")
parser.add_argument("-v", "--verbose", action='store_true', help="Prints Psithon to Python translation.")
parser.add_argument("--inplace", action='store_true',
                    help="Runs Psi4 from the source directory. !Warning! expert option.")
parser.add_argument("-l", "--psidatadir",
                    help="Specifies where to look for the Psi4 data directory. Overrides PSIDATADIR. !Warning! expert option.")
parser.add_argument("-k", "--skip-preprocessor", action='store_true',
                    help="Skips input preprocessing. !Warning! expert option.")
parser.add_argument("--json", action='store_true',
                    help="Runs a JSON input file. !Warning! experimental option.")
parser.add_argument("-t", "--test", action='store_true',
                    help="Runs smoke tests.")

# For plugins
parser.add_argument("--plugin-name", help="""\
Creates a new directory with files for writing a new plugin.
You can specify an additional argument that specifies a
template to use, for example
>>> psi4 --plugin-name mygreatcode --plugin-template mointegrals""")
parser.add_argument('--plugin-template', default='basic',
                    choices=['aointegrals', 'basic', 'dfmp2', 'mointegrals', 'scf', 'sointegrals', 'wavefunction'],
                    help='Selects new plugin template to use.')
parser.add_argument('--plugin-compile', action='store_true', help="""\
Generates a CMake command for building a plugin against this Psi4 installation.
>>> cd <plugin_directory>
>>> `psi4 --plugin-compile`
>>> make
>>> psi4""")

# print("Environment Variables\n");
# print("     PSI_SCRATCH           Directory where scratch files are written.")
# print("                           Default: $TMPDIR (or /tmp/ when not set)")
# print("                           This should be a local, not network, disk")

# parser.print_help()
args, unknown = parser.parse_known_args()
args = args.__dict__  # Namespace object seems silly

# Figure out pythonpath
cmake_install_prefix = os.path.normpath(os.path.dirname(os.path.abspath(__file__)) + os.path.sep + '..')
lib_dir = os.path.sep.join([cmake_install_prefix, "@CMAKE_INSTALL_LIBDIR@", "@PYMOD_INSTALL_LIBDIR@"])

if args["inplace"]:
    if "CMAKE_INSTALL_LIBDIR" not in lib_dir:
        raise ImportError("Cannot run inplace from a installed directory.")

    core_location = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + "core.so"
    if not os.path.isfile(core_location):
        raise ImportError("A compiled Psi4 core.so needs to be symlinked to the %s folder" % os.path.dirname(__file__))

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

# Plugin compile line
if args['plugin_compile']:
    share_cmake_dir = os.path.sep.join([cmake_install_prefix, 'share', 'cmake', 'psi4'])
    print("""cmake -C {}/psi4PluginCache.cmake -DCMAKE_PREFIX_PATH={} .""".format(share_cmake_dir, cmake_install_prefix))
    sys.exit()

if args['psiapi_path']:
    pyexe_dir = os.path.dirname("@PYTHON_EXECUTABLE@")
    print("""export PATH={}:$PATH\nexport PYTHONPATH={}:$PYTHONPATH""".format(pyexe_dir, lib_dir))
    sys.exit()

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

if args['plugin_name']:
    # This call does not return.
    psi4.plugin.create_plugin(args['plugin_name'], args['plugin_template'])

if args["test"]:
    psi4.test()
    sys.exit()

if not os.path.isfile(args["input"]):
    raise KeyError("The file %s does not exist." % args["input"])
args["input"] = os.path.normpath(args["input"])

# Setup outfile
if args["append"] is None:
    args["append"] = False
if args["output"] != "stdout":
    psi4.core.set_output_file(args["output"], args["append"])

# Set a few options
if args["prefix"] is not None:
    psi4.core.set_psi_file_prefix(args["prefix"])

psi4.core.set_num_threads(int(args["nthread"]), quiet=True)
psi4.core.set_memory_bytes(524288000, True)
psi4.extras._input_dir_ = os.path.dirname(os.path.abspath(args["input"]))
psi4.print_header()

# Prepare scratch for inputparser
if args["scratch"] is not None:
    if not os.path.isdir(args["scratch"]):
        raise Exception("Passed in scratch is not a directory (%s)." % args["scratch"])
    psi4.core.set_environment("PSI_SCRATCH", os.path.abspath(os.path.expanduser(args["scratch"])))

# If this is a json call, compute and stop
if args["json"]:

    with open(args["input"], 'r') as f:
        json_data = json.load(f)

    psi4.extras._success_flag_ = True
    psi4.extras.exit_printing()
    psi4.json_wrapper.run_json(json_data)

    with open(args["input"], 'w') as f:
        json.dump(json_data, f)

    if args["output"] != "stdout":
        os.unlink(args["output"])

    sys.exit()


# Read input
with open(args["input"]) as f:
    content = f.read()

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

    if sys.version_info >= (3, 0):
        atexit.unregister(psi4.core.clean)
    else:
        for handler in atexit._exithandlers:
            if handler[0] == psi4.core.clean:
                atexit._exithandlers.remove(handler)

# Register exit printing, failure GOTO coffee ELSE beer
import atexit
atexit.register(psi4.extras.exit_printing)

# Run the program!
try:
    exec(content)
    psi4.extras._success_flag_ = True

# Capture _any_ python error message
except Exception as exception:
    import traceback
    exc_type, exc_value, exc_traceback = sys.exc_info()
    tb_str = "Traceback (most recent call last):\n"
    tb_str += ''.join(traceback.format_tb(exc_traceback))
    tb_str += '\n'
    tb_str += type(exception).__name__
    tb_str += ': '
    tb_str += str(exception)
    psi4.core.print_out("\n")
    psi4.core.print_out(tb_str)
    psi4.core.print_out("\n")
    if psi4.core.get_output_file() != "stdout":
        print(tb_str)
    sys.exit(1)


#    elif '***HDF5 library version mismatched error***' in str(err):
#        raise ImportError("{0}\nLikely cause: HDF5 used in compilation not prominent enough in RPATH/[DY]LD_LIBRARY_PATH".format(err))

