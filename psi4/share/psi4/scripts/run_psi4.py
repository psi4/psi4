#! /usr/bin/env python
import sys
import os
import argparse

parser = argparse.ArgumentParser(description="A hybrid C++/Python quantum chemistry module.")
parser.add_argument("-i", "--input", default="input.dat", help="Input file name. Default input.dat.")
parser.add_argument("-o", "--output", help="Redirect output elsewhere.\n"
                                           "Default filename.out if input is filename"
                                           "filename.out if input is filename.in\n"
                                           "output.dat if input is input.dat\n")

parser.add_argument("-v", "--verbose", action='store_true', help="Print a lot of information.")
parser.add_argument("-V", "--version", action='store_true', help="Print version information.")

parser.add_argument("-d", "--debug", action='store_true', help="Flush the outfile at every print statement.")
parser.add_argument("-k", "--skip-preprocessor", action='store_true', help="Skips input preprocessing. Expert mode.")
parser.add_argument("-m", "--messy", action='store_true', help="Leave temporary files after the run is completed.")
parser.add_argument("-r", "--restart", action='store_true', help="Number to be used instead of process id.")
parser.add_argument("-w", "--wipe", action='store_true', help="Clean out your scratch area.")

parser.add_argument("-s", "--scratch", help="Psi4 scratch directory to use.")
parser.add_argument("-a", "--append", help="Append results to output file. Default Truncate first")
parser.add_argument("-l", "--psidatadir", help="Specify where to look for the Psi data directory. Overrides PSIDATADIR.")
parser.add_argument("-n", "--nthread", default=1, help="Number of threads to use (overrides OMP_NUM_THREADS)")
parser.add_argument("-p", "--prefix", help="Prefix name for psi files. Default psi")

# For plugins
parser.add_argument("--new-plugin", help="Creates a new directory with files for writing a "
                                         "new plugin. You can specify an additional argument "
                                         "that specifies a template to use, for example "
                                         "--new-plugin name +mointegrals")
parser.add_argument("--new-plugin-makefile", help="Creates Makefile that can be used to compile"
                                                  "plugins. The Makefile is placed in the current"
                                                  "directory.")

# print("Environment Variables\n");
# print("     PSI_SCRATCH           Directory where scratch files are written.")
# print("                           Default: $TMPDIR (or /tmp/ when not set)")
# print("                           This should be a local, not network, disk")

#parser.print_help()
args, unknown = parser.parse_known_args()
args = args.__dict__ # Namespace object seems silly

# Insert the python path
try:
    import psi4
except ImportError:
    try:
        new_path = os.path.abspath(__file__ + os.path.sep + '/../../../../../')
        print("Psi4 not found in PYTHONPATH, attempting relative import at %s" % new_path)
        print new_path
        sys.path.insert(1, new_path)
        import psi4
    except ImportError:
        raise ImportError("Could not import Psi4. This is likely due to the fact that Psi4 is not in your PYTHONPATH")
    


# Replace input/output if unknown kwargs
if len(unknown) > 0:
    args["input"] = unknown[0]
elif len(unknown) > 1:
    args["output"] = unknown[1]
elif len(unknown) > 2:
    raise KeyError("Too many unknown arguments: %s" % str(unknown))

# Figure out output arg
if args["output"] is None:
    if args["input"] == "input.dat":
        args["output"] = "output.dat"
    elif args["input"].endswith(".in"):
        args["output"] = args["input"].replace(".in", ".out")
    else:
        args["output"] = args["input"] + ".dat"

if not os.path.isfile(args["input"]):
    raise KeyError("The file %s does not exist." % args["input"])


# Figure out psidata dir
if args["psidatadir"] is not None:
    datadir = os.path.abspath(args["psidatadir"])
elif "PSIDATADIR" in os.environ.keys():
    datadir = os.path.abspath(os.environ["PSIDATADIR"])
else:
    datadir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

os.environ["PSIDATADIR"] = datadir

if not os.path.isdir(datadir):
     raise KeyError("Unable to read the Psi4 Python folder - check the PSIDATADIR environmental variable"
                    "      Current value of PSIDATADIR is %s" % datadir)

# Read input
with open(args["input"]) as f:
    content = f.read()

# Preprocess
if not args["skip_preprocessor"]:
    content = psi4.process_input(content)

# Run the program!
exec(content)


