#!/usr/bin/env python
#
# A simple wrapper to CMake, using configure-like syntax
#
# Andy Simmonett (05/13)
#
import os.path
import sys
try:
    import argparse
except ImportError:
    print(
"""Error: Your Python interpreter (version %s) is older than 2.7.
Since psi4 needs 2.7, too, please consider upgrading. Get the python
development libraries (provides "python-config") while you're at it.""" % (sys.version[:6].strip()))
    sys.exit(1)
import subprocess


def execute(command, die_on_error=True):
    #print("\tExecuting %s" % command)
    failed = 0
    process = subprocess.Popen(command)
    (stdout, stderr) = process.communicate()
    status = process.returncode
    if status:
        print("XXXXXX Error executing %s XXXXXX" % command)
        failed = 1
        if die_on_error:
            sys.exit()
    return failed

# Find the location of this script, which is where CMake will be pointed to
thisscript = os.path.realpath(__file__)
scriptdir = os.path.dirname(thisscript)

parser = argparse.ArgumentParser(description='Configure Psi4 using cmake',
formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#
# Add options
#
blankstring = "No default"  # This looks prettier when the user prints out the help
# The C compiler
parser.add_argument('--with-cc',
                    metavar='= CC',
                    type=str,
                    default=blankstring,
                    help='The C compiler to use.  If not specified the environmental variable $CC will be used. '\
                         + 'If that is not set, CMake will use the first working compiler it finds.')
# The C++ compiler
parser.add_argument('--with-cxx',
                    metavar='= CXX',
                    type=str,
                    default=blankstring,
                    help='The C++ compiler to use.  If not specified the environmental variable $CXX will be used. '\
                         + 'If that is not set, CMake will use the first working compiler it finds.')
# CXXFLAGS flags
parser.add_argument('--with-cxxflags',
                    metavar='= CXXFLAGS',
                    type=str,
                    default=blankstring,
                    help="Any extra flags to pass to the C++ compiler.")
parser.add_argument('--with-cmake',
                    type=str,
                    default='cmake',
                    help="The CMake executable to use.  If not specified, attempts to use cmake in your path")
# Debug symbols
debuggroup = parser.add_mutually_exclusive_group()
debuggroup.add_argument('--with-debug',
                        action="store_true",
                        help='Add debug flags.')
debuggroup.add_argument('--without-debug',
                        action="store_false",
                        help='Do not add debug flags.')
# ERD package
erdgroup = parser.add_mutually_exclusive_group()
erdgroup.add_argument('--with-erd',
                        action="store_true",
                        help='Add support for the ERD integral package.')
erdgroup.add_argument('--without-erd',
                        action="store_false",
                        help='Do not use the ERD integral package.')
#MPI package
mpigroup = parser.add_mutually_exclusive_group()
mpigroup.add_argument('--with-mpi',
                        action="store_true",
                        help='Add support for MPI.')
mpigroup.add_argument('--without-mpi',
                        action="store_false",
                        help='Do not use MPI.')
# The Fortran compiler
parser.add_argument('--with-f77',
                    metavar='= F77',
                    type=str,
                    default=blankstring,
                    help='The Fortran compiler to use.  If not specified the environmental variable $F77 will be used. '\
                         + 'If that is not set, CMake will use the first working compiler it finds.')
# F77FLAGS flags
parser.add_argument('--with-f77flags',
                    metavar='= F77FLAGS',
                    type=str,
                    default=blankstring,
                    help="Any extra flags to pass to the Fortran compiler.")
# F77 symbol
parser.add_argument('--with-f77symbol',
                    metavar='= lcu | lc | uc | ucu | detect',
                    type=str,
                    choices=['lcu', 'lc', 'uc', 'ucu', 'detect'],
                    default='detect',
                    help="The Fortran compiler name mangling convention, used for linking external Fortran libraries, such as BLAS. Values are lcu (lower case with traling underscore), lc (lower case), ucu (upper case with trailing underscore), uc (upper case). If omitted CMake will detect it automatically if a Fortran compiler is present, if not it will use lcu.")
# Lapack
parser.add_argument('--with-lapack-libs',
                    metavar='= LAPACKLIBS',
                    type=str,
                    default=blankstring,
                    help='The flags to be passed to the linker to use BLAS and LAPACK. For modern mkl, with intel compilers, you can set this to -mkl (N.B. not -lmkl).  If omitted, CMake will try to find a working version for you.')
parser.add_argument('--with-lapack-incs',
                    metavar='= LAPACKINCS',
                    type=str,
                    default=blankstring,
                    help='The flags to be passed to the compiler to use BLAS and LAPACK. If omitted, CMake will try to find a working version for you.')

# LD flags
parser.add_argument('--with-ldflags',
                    metavar='= LDFLAGS',
                    type=str,
                    default=blankstring,
                    help="Any extra flags to pass to the linker (usually -Llibdir -llibname type arguments). You shouldn't need this.")
# Libint max A.M.
parser.add_argument('--with-max-am-eri',
                    metavar="= MAX_ANGULAR_MOMENTUM",
                    type=int,
                    default=5,
                    help='The maximum angular momentum level (1=p, 2=d, 3=f, etc.) for the libint and libderiv packages.  Note: A value of N implies a maximum first derivative of N-1, and maximum second derivative of N-2.')
# Plugins
pluginsgroup = parser.add_mutually_exclusive_group()
pluginsgroup.add_argument('--with-plugins',
                    action="store_true",
                    help='Compile with support for pluginss.')
pluginsgroup.add_argument('--without-plugins',
                    action="store_false",
                    help='Compile without support for plugins.')
# Prefix
parser.add_argument('--prefix',
                    metavar='= PREFIX',
                    type=str,
                    default="/usr/local/psi4",
                    help='Installation directory for Psi4')
# Boost
parser.add_argument('--with-boost-incdir',
                    metavar='= BOOST',
                    type=str,
                    default=blankstring,
                    help='The includes directory for boost.  If this is left blank cmake will attempt to find one on your system.  Failing that it will build one for you')
# Boost-libdir
parser.add_argument('--with-boost-libdir',
                    metavar='= BOOSTLIB',
                    type=str,
                    default=blankstring,
                    help='The libraries directory for boost.  If this is left blank cmake will attempt to find one on your system.  Failing that it will build one for you')
# Python
parser.add_argument('--with-python',
                    metavar='= PYTHON',
                    type=str,
                    default=blankstring,
                    help='The Python interpreter (development version) to use.  CMake will detect one automatically, if omitted.')
# Optimization
optgroup = parser.add_mutually_exclusive_group()
optgroup.add_argument('--with-opt',
                    action="store_false",
                    help='Add optimization flags.')
optgroup.add_argument('--without-opt',
                    action="store_true",
                    help='Do not add optimization flags.')

# External plugins
parser.add_argument('--with-external',
                    type=str,
                    nargs='*',
                    choices=['gpu_dfcc', 'dummy_plugin'],
                    help='Any external plugins to download and build.')


def dict_to_list(dictionary):
    l = []
    for k in sorted(dictionary.keys()):
        s = "-D" + k + "="
        if not dictionary[k]:
            s += '""'
        else:
            if isinstance(dictionary[k], list):
                s += " ".join(dictionary[k])
                #print l
            elif isinstance(dictionary[k], str):
                s += dictionary[k]
            elif isinstance(dictionary[k], (int, float)):
                s += str(dictionary[k])
            else:
                raise Exception("Unexpected keyword type: %r" % dictionary[k])
        l.append(s)
    return l


def dict_to_string(dictionary):
    """Converts a dictionary of keywords into a string of arguments to pass to CMake"""
    string = ""
    for k in sorted(dictionary.keys()):
        string += " -D" + k + "="
        if not dictionary[k]:
            string += '""'
        else:
            if isinstance(dictionary[k], list):
                string += '"' + " ".join(dictionary[k]) + '"'
                #print string
            elif isinstance(dictionary[k], str):
                string += dictionary[k]
            elif isinstance(dictionary[k], (int, float)):
                string += str(dictionary[k])
            else:
                raise Exception("Unexpected keyword type: %r" % dictionary[k])
    return string

#
# Convert user options to CMake arguments
#

args = parser.parse_args()
cmakeflags = {}

# PREFIX
cmakeflags['PREFIX'] = args.prefix
# MAX-AM-ERI
cmakeflags['MAX_AM_ERI'] = args.with_max_am_eri
# CXX/F77 FLAGS
cmakeflags['CXXFLAGS'] = ['']
cmakeflags['F77FLAGS'] = ['']
if args.without_opt:
    cmakeflags['CXXFLAGS'].append("-O0")
    cmakeflags['F77FLAGS'].append("-O0")
else:
    cmakeflags['CXXFLAGS'].append("-O2")
    cmakeflags['F77FLAGS'].append("-O2")

if args.with_debug:
    cmakeflags['CXXFLAGS'].append("-g")
    cmakeflags['F77FLAGS'].append("-g")

if args.with_erd:
    cmakeflags['USEERD'] = ["TRUE"]

if args.with_mpi:
    cmakeflags['USEMPI'] = ["TRUE"]

#For some reason making plugins is false if we are making them
if args.with_plugins:
    cmakeflags['CXXFLAGS'].append("-fPIC")
    cmakeflags['F77FLAGS'].append("-fPIC")

if args.with_cxxflags != blankstring:
    cmakeflags['CXXFLAGS'].append(args.with_cxxflags)
if args.with_f77flags != blankstring:
    cmakeflags['F77FLAGS'].append(args.with_f77flags)

cmakeflags['LDFLAGS'] = []
# LDFLAGS
if args.with_ldflags != blankstring:
    cmakeflags['LDFLAGS'] = [args.with_ldflags]
# F77
if args.with_f77 != blankstring:
    cmakeflags['CMAKE_Fortran_COMPILER'] = args.with_f77
# CC
if args.with_cc != blankstring:
    cmakeflags['CMAKE_C_COMPILER'] = args.with_cc
# CXX
if args.with_cxx != blankstring:
    cmakeflags['CMAKE_CXX_COMPILER'] = args.with_cxx
#BOOST
if args.with_boost_incdir != blankstring:
    cmakeflags['Boost_INCLUDE_DIR'] = args.with_boost_incdir
if args.with_boost_libdir != blankstring:
    cmakeflags['BOOST_LIBRARYDIR'] = args.with_boost_libdir
# PYTHON
if args.with_python != blankstring:
    cmakeflags['PYTHON'] = args.with_python
# LAPACK
if args.with_lapack_libs != blankstring:
    cmakeflags['LAPACKLIBS'] = [args.with_lapack_libs]
if args.with_lapack_incs != blankstring:
    cmakeflags['LAPACKINCS'] = [args.with_lapack_incs]
# F77SYMBOL
cmakeflags['F77SYMBOL'] = args.with_f77symbol
# External
if args.with_external:
    for arg in args.with_external:
        cmakeflags['USEEXT_%s' % (arg.upper())] = ["TRUE"]
CMAKE_CALL = args.with_cmake
args = [CMAKE_CALL, scriptdir]
args.extend(dict_to_list(cmakeflags))
execute(args)
