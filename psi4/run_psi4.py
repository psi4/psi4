#!@Python_EXECUTABLE@

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

import argparse
import atexit
import datetime
import json
import os
import sys
import warnings
from pathlib import Path

# yapf: disable
parser = argparse.ArgumentParser(description="Psi4: Open-Source Quantum Chemistry", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i", "--input", default="input.dat",
                    help="Input file name. Default: input.dat.")
parser.add_argument("-o", "--output", help="""\
Redirect output elsewhere.
Default: when input filename is 'input.dat', 'output.dat'.
Otherwise, output filename defaults to input filename with
'.out' extension""")
parser.add_argument("-a", "--append", action='store_true',
                    help="Appends results to output file. Default: Truncate first")
parser.add_argument("-V", "--version", action='store_true',
                    help="Prints version information.")
parser.add_argument("-n", "--nthread", default=1,
                    help="Number of threads to use. Psi4 disregards OMP_NUM_THREADS/MKL_NUM_THREADS.")
parser.add_argument("--memory", default=524288000,
                    help="The amount of memory to use. Can be specified with units (e.g., '10MB') otherwise bytes is assumed.")
parser.add_argument("-s", "--scratch",
                    help="Scratch directory to use. Overrides PSI_SCRATCH.")
parser.add_argument("-m", "--messy", action='store_true',
                    help="Leaves temporary files after the run is completed.")
# parser.add_argument("-d", "--debug", action='store_true', help="Flush the outfile at every print statement.")
# parser.add_argument("-r", "--restart", action='store_true', help="Number to be used instead of process id.")
# parser.add_argument("-p", "--prefix", help="Prefix name for psi files. Default psi")
parser.add_argument("--psiapi-path", action='store_true',
                    help="""Generates a bash command to source correct Python """
                         """interpreter and path for ``python -c "import psi4"``""")
parser.add_argument("--module", action='store_true',
                    help="""Generates the path to PsiAPI loading. That is, the following file exists: `psi4 --module`/psi4/__init__.py . Also, adding `psi4 --module` to PYTHONPATH allows "import psi4".""")
parser.add_argument("-v", "--verbose", action='store_true', help="Prints Psithon to Python translation.")
parser.add_argument("--inplace", action='store_true',
                    help="Runs Psi4 from the source directory. !Warning! expert option.")
parser.add_argument("-l", "--psidatadir",
                    help="Specifies where to look for the Psi4 data directory. Overrides PSIDATADIR. !Warning! expert option.")
parser.add_argument("-k", "--skip-preprocessor", action='store_true',
                    help="Skips input preprocessing. !Warning! expert option.")
parser.add_argument("--qcschema", "--schema", action='store_true',
                    help="Runs input file as QCSchema. Can either be JSON or MessagePack input. Use `--output` to not overwrite schema input file.")
parser.add_argument("--json", action='store_true',
                    help="Runs a JSON input file. !Warning! depcrated option in 1.4, use --qcschema instead.")
parser.add_argument("-t", "--test", nargs='?', const='smoke', default=None,
                    help="Runs pytest tests (requires pytest installed). If `pytest-xdist` installed, parallel with `--nthread`.")
parser.add_argument("--mdi", default=None,
                    help="Sets MDI configuration options")
parser.add_argument("--loglevel", default=20,
                    help="Sets logging level: WARN=30, INFO=20, DEBUG=10.")
parser.add_argument("--inherit-loglevel", action='store_true',
                    help="Takes no action on logging level, not even setting the default, thereby inheriting from existing logger.")

# For plugins
parser.add_argument("--plugin-name", help="""\
Creates a new directory with files for writing a new plugin.
You can specify an additional argument that specifies a
template to use, for example
>>> psi4 --plugin-name mygreatcode --plugin-template mointegrals""")
parser.add_argument('--plugin-template',
                    choices=['ambit', 'aointegrals', 'basic', 'dfmp2', 'mointegrals', 'scf', 'sointegrals', 'wavefunction'],
                    help='Selects new plugin template to use.')
parser.add_argument('--plugin-compile', action='store_true', help="""\
Generates a CMake command for building a plugin against this Psi4 installation.
>>> cd <plugin_directory>
>>> `psi4 --plugin-compile`
>>> make
>>> psi4""")
# yapf: enable

# print("Environment Variables\n");
# print("     PSI_SCRATCH           Directory where scratch files are written.")
# print("                           Default: $TMPDIR (or /tmp/ when not set)")
# print("                           This should be a local, not network, disk")

# parser.print_help()
args, unknown = parser.parse_known_args()
args = args.__dict__  # Namespace object seems silly

# Figure out paths
# * some full paths are computed here using the prefix, but all outputs are relative to __file__, so relocatability preserved
# * note that all path entities are directories except for "executable" that is a file
executable = Path(__file__).resolve()
psi4_exe_loc = executable.parent

prefix = Path("@CMAKE_INSTALL_PREFIX@".replace("\\", "/"))
cmake_install_bindir = "@CMAKE_INSTALL_BINDIR@".replace("\\", "/")
cmake_install_datadir = "@CMAKE_INSTALL_DATADIR@".replace("\\", "/")
cmake_install_libdir = "@CMAKE_INSTALL_LIBDIR@".replace("\\", "/")
pymod_install_libdir = "@PYMOD_INSTALL_LIBDIR@".lstrip("/")
psi4_install_cmakedir = "@psi4_INSTALL_CMAKEDIR@".replace("\\", "/")
full_pymod = (prefix / cmake_install_libdir / pymod_install_libdir / "psi4").resolve()
full_data = prefix / cmake_install_datadir / "psi4"
full_bin = prefix / cmake_install_bindir
full_cmake = prefix / psi4_install_cmakedir
rel_pymod = os.path.relpath(full_pymod, start=full_bin)
rel_data = os.path.relpath(full_data, start=full_bin)
rel_cmake = os.path.relpath(full_cmake, start=full_bin)

data_dir = psi4_exe_loc.joinpath(rel_data).resolve()
psi4_module_loc = psi4_exe_loc.joinpath(rel_pymod).resolve()
cmake_dir = psi4_exe_loc.joinpath(rel_cmake).resolve()
cmake_install_prefix = os.path.commonpath([data_dir, psi4_module_loc, psi4_exe_loc, cmake_dir])
lib_dir = str(psi4_module_loc.parent)
bin_dir = str(psi4_exe_loc)
share_cmake_dir = str(cmake_dir)

if args["inplace"]:
    # not tested after pathlib adjustments
    if "CMAKE_INSTALL_LIBDIR" not in lib_dir:
        raise ImportError("Cannot run inplace from an installed directory.")

    import sysconfig
    core_location = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + "core" + sysconfig.get_config_var("EXT_SUFFIX")
    if not os.path.isfile(core_location):
        raise ImportError("A compiled Psi4 core{} needs to be symlinked to the {} folder".format(
            sysconfig.get_config_var("EXT_SUFFIX"), os.path.dirname(__file__)))

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
    raise KeyError(f"Too many unknown arguments: {unknown}")

# Figure out output arg
if (args["output"] is None) and (args["qcschema"] is False):
    if args["input"] == "input.dat":
        args["output"] = "output.dat"
    else:
        pinput = Path(args["input"])
        presuffix = pinput.suffix if pinput.suffix in [".out", ".log"] else ""
        args["output"] = str(pinput.with_suffix(presuffix + ".out"))

# Plugin compile line
if args['plugin_compile']:
    plugincachealongside = os.path.isfile(share_cmake_dir + os.path.sep + 'psi4PluginCache.cmake')
    if plugincachealongside:
        print(f"""cmake -C {share_cmake_dir}/psi4PluginCache.cmake -DCMAKE_PREFIX_PATH={cmake_install_prefix} .""")
        sys.exit()
    else:
        print("""Install "psi4-dev" via `conda install psi4-dev -c psi4[/label/dev]`, then reissue command.""")

if args['psiapi_path']:
    pyexe_dir = os.path.dirname("@Python_EXECUTABLE@")
    print(f"""export PATH={pyexe_dir}:$PATH  # python interpreter\nexport PATH={bin_dir}:$PATH  # psi4 executable\nexport PYTHONPATH={lib_dir}:$PYTHONPATH  # psi4 pymodule""")
    # TODO Py not quite right on conda Windows and Psi include %PREFIX$. but maybe not appropriate for Win anyways
    sys.exit()

if args["module"]:
    print(lib_dir)
    sys.exit()

# Transmit any argument psidatadir through environ
if args["psidatadir"] is not None:
    data_dir = os.path.abspath(os.path.expanduser(args["psidatadir"]))
    os.environ["PSIDATADIR"] = data_dir

### Actually import psi4 and apply setup ###

# Arrange for warnings to ignore everything except the message
def custom_formatwarning(msg, *args, **kwargs):
    return str(msg) + '\n'

warnings.formatwarning = custom_formatwarning

# Import installed psi4
sys.path.insert(1, lib_dir)
import psi4  # isort:skip

if args["version"]:
    print(psi4.__version__)
    sys.exit()

# Prevents a poor option combination
if args['plugin_template'] and (not args['plugin_name']):
    raise KeyError("Please specify a '--plugin-name' for your plugin template!")

if args['plugin_name']:

    # Set the flag
    if not args['plugin_template']:
        args['plugin_template'] = 'basic'

    # This call does not return.
    psi4.pluginutil.create_plugin(args['plugin_name'], args['plugin_template'])

    sys.exit()

if args["test"] is not None:
    if args["test"] not in ['smoke', 'quick', 'full', 'long', 'smoke_not_d2ints', 'quick_not_d2ints']:
        raise KeyError("The test category {} does not exist.".format(args["test"]))

    nthread = int(args["nthread"])
    if nthread == 1:
        extras = []
    else:
        extras = ['-n', str(nthread)]
    if args["test"] == "smoke_not_d2ints":
        extras.extend(["-m", "smoke and not d2ints"])
    elif args["test"] == "quick_not_d2ints":
        extras.extend(["-m", "(quick or smoke) and not d2ints"])
    retcode = psi4.test(args["test"], extras=extras)
    sys.exit(retcode)

if not os.path.isfile(args["input"]):
    raise KeyError("The file %s does not exist." % args["input"])
args["input"] = os.path.normpath(args["input"])

# Setup scratch_messy
_clean_functions = [psi4.core.clean, psi4.extras.clean_numpy_files]

# Set a few options (silently)
psi4.core.set_num_threads(int(args["nthread"]), quiet=True)
psi4.set_memory(args["memory"], quiet=True)
psi4.extras._input_dir_ = os.path.dirname(os.path.abspath(args["input"]))

# Setup outfile
if args["append"] is None:
    args["append"] = False
if args["inherit_loglevel"] is None:
    args["inherit_loglevel"] = False
if args["qcschema"] is False:
    psi4.set_output_file(
        args["output"],
        args["append"],
        loglevel=int(args["loglevel"]),
        inherit_loglevel=args["inherit_loglevel"])

start_time = datetime.datetime.now()

# Initialize MDI
if args["mdi"] is not None:
    psi4.mdi_engine.mdi_init(args["mdi"])

# Prepare scratch for inputparser
if args["scratch"] is not None:
    if not os.path.isdir(args["scratch"]):
        raise Exception("Passed in scratch is not a directory (%s)." % args["scratch"])
    psi4.core.IOManager.shared_object().set_default_path(os.path.abspath(os.path.expanduser(args["scratch"])))

# If this is a json or qcschema call, compute and stop
if args["qcschema"]:
    import qcelemental as qcel

    # Handle the reading and deserialization manually
    filename = args["input"]
    if filename.endswith("json"):
        encoding = "json"
        with open(filename, 'r') as handle:
            # No harm in attempting to read json-ext over json
            data = qcel.util.deserialize(handle.read(), "json-ext")
    elif filename.endswith("msgpack"):
        encoding = "msgpack-ext"
        with open(filename, 'rb') as handle:
            data = qcel.util.deserialize(handle.read(), "msgpack-ext")
    else:
        raise Exception("qcschema files must either end in '.json' or '.msgpack'.")

    psi4.extras._success_flag_ = True
    clean = True
    if args["messy"]:
        clean = False
        for func in _clean_functions:
            atexit.unregister(func)

    ret = psi4.schema_wrapper.run_qcschema(data, clean=clean)

    if args["output"] is not None:
        filename = args["output"]
        if filename.endswith("json"):
            encoding = "json"
        elif filename.endswith("msgpack"):
            encoding = "msgpack-ext"
        # Else write with whatever encoding came in

    if encoding == "json":
        with open(filename, 'w') as handle:
            handle.write(ret.serialize(encoding))
    elif encoding == "msgpack-ext":
        with open(filename, 'wb') as handle:
            handle.write(ret.serialize(encoding))

    sys.exit()

if args["json"]:

    with open(args["input"], 'r') as f:
        json_data = json.load(f)

    psi4.extras._success_flag_ = True
    psi4.extras.exit_printing(start_time)
    json_data = psi4.schema_wrapper.run_json(json_data)

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
    for func in _clean_functions:
        atexit.unregister(func)

# Register exit printing, failure GOTO coffee ELSE beer
atexit.register(psi4.extras.exit_printing, start_time=start_time)

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
    tb_str += ''.join(traceback.format_exception_only(type(exception), exception))
    psi4.core.print_out("\n")
    psi4.core.print_out(tb_str)
    psi4.core.print_out("\n\n")

    in_str = "Printing out the relevant lines from the Psithon --> Python processed input file:\n"
    lines = content.splitlines()
    try:
        suspect_lineno = traceback.extract_tb(exc_traceback)[1].lineno - 1  # -1 for 0 indexing
    except IndexError:
        # module error where lineno useless (e.g., `print "asdf"`)
        pass
    else:
        first_line = max(0, suspect_lineno - 5)  # Try to show five lines back...
        last_line = min(len(lines), suspect_lineno + 6)  # Try to show five lines forward
        for lineno in range(first_line, last_line):
            mark = "--> " if lineno == suspect_lineno else "    "
            in_str += mark + lines[lineno] + "\n"
        psi4.core.print_out(in_str)

    # extract exception message and print it in a box for attention.
    ex = ','.join(traceback.format_exception_only(type(exception), exception))
    ex_list = ex.split(":", 1)[-1]
    error = ''.join(ex_list)
    psi4.core.print_out(psi4.driver.p4util.text.message_box(error))
    if psi4.core.get_output_file() != "stdout":
        print(tb_str)
        print(in_str)
        print(psi4.driver.p4util.text.message_box(error))
    sys.exit(1)
