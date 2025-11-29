#!/usr/bin/env python

import os
import re
import sys
import json
import shutil
import argparse
import itertools
import unicodedata
import textwrap as _textwrap
from pathlib import Path
from subprocess import run
from typing import Dict

try:
    from ruamel.yaml import YAML
except ModuleNotFoundError:
    try:
        import yaml
    except ModuleNotFoundError:
        raise ModuleNotFoundError("Python module PyYAML not found. Solve by installing it: "
            "`conda install -c conda-forge ruamel.yaml` or `pip install ruamel.yaml` or "
            "`conda install -c conda-forge pyyaml` or `pip install PyYAML`.")
    else:
        yaml_load = yaml.safe_load
else:
    yaml=YAML(typ='safe')
    yaml_load = yaml.load

codedeps_yaml = Path(__file__).resolve().parent.parent / "codedeps.yaml"
cmake_S = os.path.relpath(codedeps_yaml.parent, start=Path.cwd())


def native_platform() -> str:
    """Return conda platform code.
    This can't distinguish the metal chip for osx, so prefer to use platform from conda info directly.

    """
    import platform
    if sys.platform.startswith("linux"):
        return "linux-64"
    elif sys.platform == "darwin":
        if platform.machine() == "arm64":
            return "osx-arm64"
        elif platform.machine() == "x86_64":
            # note that an Apple Silicon processor running an Intel-compiled Python sorts here
            return "osx-64"
    elif sys.platform.startswith("win"):
        return "win-64"


def conda_list(*, name: str = None, prefix: str = None) -> Dict:
    """Return `conda list` in json format.

    Thanks, https://stackoverflow.com/a/56363822
       SO convinced me that subprocess was better than import (following block)
       Also, the conda.cli only works on base env (where conda pkg installed)

    #import conda.cli.python_api as Conda
    #env_list_json, stderr, rc = Conda.run_command(Conda.Commands.LIST, ["--json"])
    #env_list_dict = json.loads(env_list_json)

    """
    condaexe = "conda.bat" if (os.name == "nt") else "conda"
    if name:
        proc = run([condaexe, "list", "--json", "--name", name], text=True, capture_output=True)
    elif prefix:
        proc = run([condaexe, "list", "--json", "--prefix", prefix], text=True, capture_output=True)
    else:
        proc = run([condaexe, "list", "--json"], text=True, capture_output=True)
    return json.loads(proc.stdout)


def conda_info() -> Dict:
    condaexe = "conda.bat" if (os.name == "nt") else "conda"
    proc = run([condaexe, "info", "--json"], text=True, capture_output=True)
    return json.loads(proc.stdout)


def strike(text: str, tilde: bool = False) -> str:
    if (os.name == "nt") or tilde:
        # Windows has a probably correctable problem with unicode, but I can't iterate it quickly, so use tilde for strike.
        #   UnicodeEncodeError: 'charmap' codec can't encode character '\u0336' in position 3: character maps to <undefined>
        return "~" + text + "~"
    else:
        return ''.join(itertools.chain.from_iterable(zip(text, itertools.repeat('\u0336'))))


def visible_width(chars):
    return sum(not unicodedata.combining(ch) for ch in chars)


def compute_width(head, width):
    return width + len(head) - visible_width(head)


re_pkgline = re.compile("(?P<suppress>//)?(?P<chnl>.*::)?(?P<pkg>[A-Za-z0-9_-]+)(?P<constraint>[=!<> ].*)?")

conda_available = shutil.which("conda")
mamba_available = shutil.which("mamba")
pm_available = (conda_available or mamba_available)

if pm_available:
    conda_info_dict = conda_info()
    conda_platform_native = conda_info_dict["platform"]
    conda_prefix = conda_info_dict["active_prefix"]
    conda_prefix_short = conda_info_dict["active_prefix_name"]
    conda_host = conda_info_dict["env_vars"].get("CONDA_TOOLCHAIN_HOST", None)  # None if no compilers in env
    conda_list_struct = conda_list()
    base_prefix = conda_info_dict["conda_prefix"]  # env with conda cmd
    base_list_struct = conda_list(prefix=base_prefix)
else:
    conda_platform_native = native_platform()
    conda_prefix = "(no prefix)"
    conda_prefix_short = "(no prefix)"
    conda_host = None
    conda_list_struct = {}
    base_list_struct = {}
if conda_prefix:  # None if base env not activated
    conda_prefix = Path(conda_prefix).as_posix()


conda_list_pkgver = {item["name"]: item["version"] for item in conda_list_struct}
conda_lapack_variant = None  # None if no c-f libblas in env
for itm in conda_list_struct:
    if itm["name"] == "libblas":
        conda_lapack_variant = itm["build_string"].split("_")[-1]
    if itm["name"] == "openblas":
        conda_openblas_variant = itm["build_string"].split("_")[0]

# TODO handle conda_host None in base env with no compilers present
# TODO handle conda_lapack_variant None in base env with no lapack present

conda_libmamba_available = False
for itm in base_list_struct:
    if itm["name"] == "conda-libmamba-solver":
        conda_libmamba_available = True
        break

### conda/mamba

solver_choices = ["conda-else-mamba"]
solver_help = [f"""(default: conda-else-mamba)
conda-else-mamba:
    Use conda else mamba to solve environments.
    Can instead adjust on cmdline, so argument mostly for printing."""]
if conda_available:
    solver_choices.append("conda")
    solver_help.append("""conda:
    Use `conda` to solve environments.""")
else:
    solver_choices.append(strike("conda"))
    solver_help.append(f"""{strike('conda')}
    Can't use `conda` to solve environments
    because package (conda) not installed in base env.""")
if mamba_available:
    solver_choices.append("mamba")
    solver_help.append("""mamba:
    Use `mamba` to solve environments. UNTESTED""")
else:
    solver_choices.append(strike("mamba"))
    solver_help.append(f"""{strike('mamba')}
    Can't use `mamba` to solve environments
    because packages (mamba?) not installed in base env.""")
if conda_libmamba_available:
    solver_choices.append("conda-libmamba")
    solver_help.append("""conda-libmamba:
    Use `conda ... --solver=libmamba` to solve environments.""")
else:
    solver_choices.append(strike("conda-libmamba"))
    solver_help.append(f"""{strike('conda-libmamba')}
    Can't use `conda ... --solver=libmamba` to solve environments
    because packages (conda or conda-libmamba-solver) not installed in base env.""")
if conda_available:
    solver_choices.append("conda-classic")
    solver_help.append("""conda-classic:
    Use `conda ... --solver=classic` to solve environments.""")
else:
    solver_choices.append(strike("conda-classic"))
    solver_help.append(f"""{strike('conda-classic')}
    Can't use `conda ... --solver=classic` to solve environments
    Can't use `conda` to solve environments
    because package (conda) not installed in base env.""")


### blas/lapack

lapack_conda_native = "blas-devel" in conda_list_pkgver
if not lapack_conda_native and conda_lapack_variant == "mkl" and "mkl-devel" in conda_list_pkgver:
    lapack_conda_native = True
lapack_default = "conda" if lapack_conda_native else "byo"
lapack_choices = []
lapack_help = []
conda_lapack_platform = {
    "linux-64":  ["mkl",               "openblas", "blis", "netlib"],
    "osx-64":    ["mkl", "accelerate", "openblas", "blis", "netlib"],
    "osx-arm64": [       "accelerate", "openblas",         "netlib"],
    "win-64":    ["mkl",               "openblas", "blis", "netlib"],
}

lapack_help.append(f"""(default: {lapack_default})""")
# byo
lapack_choices.append("byo")
lapack_help.append("""byo:
    Omit blas/lapack from cache to engage self-provided libraries.
    CMake configuration may still use conda libraries if
    present and not contravened.""")

# conda
if lapack_conda_native:
    lapack_choices.append("conda")
    lapack_choices.append(conda_lapack_variant)
    lapack_help.append(f"""conda ({conda_lapack_variant}):
{conda_lapack_variant} (conda):
    Engage conda-provided (blas-devel=*=*{conda_lapack_variant}) blas/lapack libraries.""")
elif conda_lapack_variant:
    lapack_choices.append(strike("conda"))
    lapack_choices.append(strike(conda_lapack_variant))
    lapack_help.append(f"""{strike(f'conda ({conda_lapack_variant})')}:
{strike(f'{conda_lapack_variant} (conda)')}:
    Can't engage conda-provided blas/lapack libraries because
    package (blas-devel=*=*{conda_lapack_variant} etc.) not installed in current env.""")
else:
    conda_lapack_platform_this = conda_lapack_platform[conda_platform_native]
    options_lap = [f"{strike(f'{lap} (conda)')}:" for lap in conda_lapack_platform_this]
    lapack_default_pkg = conda_lapack_platform_this[0]
    lapack_choices.append(strike("conda"))
    for lap in conda_lapack_platform_this:
        lapack_choices.append(strike(lap))
    lapack_help.append(f"""{strike(f'conda ({lapack_default_pkg})')}:\n""" + "\n".join(options_lap) + f"""
    Can't engage conda-provided blas/lapack libraries because
    package (blas-devel=*=*{lapack_default_pkg} etc.) not installed in current env.""")


def lapack_type(arg):
    if arg.lower() == "conda":
        return conda_lapack_variant
    for itm in lapack_choices:
        if arg.lower() == itm.lower():
            return itm


### compilers

byo_icpc = shutil.which("icpc")
byo_icpx = shutil.which("icpx")
conda_native = "cxx-compiler" in conda_list_pkgver
conda_icpx = "dpcpp_linux-64" in conda_list_pkgver

compiler_default = "conda" if conda_native else "byo"
compiler_choices = []
compiler_help = []

if conda_platform_native == "linux-64":
    compiler_help.append(f"""(default: {compiler_default})""")

    # byo
    compiler_choices.append("byo")
    compiler_help.append("""byo:
    Omit compilers from cache to engage self-provided compilers.
    CMake configuration may still use conda compilers if
    present and not contravened.""")

    # conda
    if conda_native:
        compiler_choices.append("conda")
        compiler_choices.append("GNU")
        compiler_help.append(f"""conda (GNU):
GNU (conda):
    Engage conda-provided (cxx-compiler) gcc/g++ compilers.""")
    else:
        compiler_choices.append(strike("conda"))
        compiler_choices.append(strike("GNU"))
        compiler_help.append(f"""{strike('conda (GNU)')}:
{strike('GNU (conda)')}:
    Can't engage conda-provided gcc/g++ compilers because
    package (cxx-compiler) not installed in current env.""")

    # IntelLLVM
    if conda_icpx and conda_native:
        compiler_choices.append("IntelLLVM")
        compiler_help.append(f"""IntelLLVM (conda) atop GNU (conda):
    Engage conda-provided (dpcpp_linux-64) icx/icpx compilers backed
    by conda-provided (cxx-compiler) gcc/g++.""")
    elif not conda_icpx:
        compiler_choices.append(strike("IntelLLVM"))
        compiler_help.append(f"""{strike('IntelLLVM (conda) atop GNU (conda)')}:
    Can't engage conda-provided icx/icpx compilers because
    package (dpcpp_linux-64) not installed in current env.""")
    elif conda_icpx and not conda_native:
        compiler_choices.append(strike("IntelLLVM"))
        compiler_help.append(f"""{strike('IntelLLVM (conda) atop GNU (conda)')}:
    Can't engage conda-provided icx/icpx compilers backed by
    conda-provided gcc/g++ because package (cxx-compiler)
    not installed in current env.""")

    # Intel
    # Intel-multiarch
    if byo_icpc and conda_native:
        compiler_choices.append("Intel")
        compiler_help.append(f"""Intel Classic (byo) atop GNU (conda):
    Engage self-provided icc/icpc compilers backed by conda-
    provided (cxx-compiler) gcc/g++.""")
        compiler_choices.append("Intel-multiarch")
        compiler_help.append(f"""Intel Classic w/multiarch (byo) atop GNU (conda):
    Engage self-provided icc/icpc compilers backed by conda-
    provided (cxx-compiler) gcc/g++ PLUS compile for multiple
    architectures (useful for cluster deployments).""")
    elif not byo_icpc:
        compiler_choices.append(strike("Intel"))
        compiler_help.append(f"""{strike('Intel Classic (byo) atop GNU (conda)')}:
    Can't engage self-provided icc/icpc compilers because
    not in envvar:PATH.""")
        compiler_choices.append(strike("Intel-multiarch"))
        compiler_help.append(f"""{strike('Intel Classic w/multiarch (byo) atop GNU (conda)')}:
    Can't engage self-provided icc/icpc compilers because
    not in envvar:PATH.""")
    elif byo_icpc and not conda_native:
        compiler_choices.append(strike("Intel"))
        compiler_help.append(f"""{strike('Intel Classic (byo) atop GNU (conda)')}:
    Can't engage self-provided icc/icpc compilers backed by
    conda-provided gcc/g++ because package (cxx-compiler)
    not installed in current env.""")
        compiler_choices.append(strike("Intel-multiarch"))
        compiler_help.append(f"""{strike('Intel Classic w/multiarch (byo) atop GNU (conda)')}:
    Can't engage self-provided icc/icpc compilers backed by
    conda-provided gcc/g++ because package (cxx-compiler)
    not installed in current env.""")

elif conda_platform_native in ["osx-64", "osx-arm64"]:
    compiler_help.append(f"""(default: {compiler_default})""")

    # byo
    compiler_choices.append("byo")
    compiler_help.append("""byo:
    Omit compilers from cache to engage self-provided compilers.
    CMake configuration may still use conda compilers if
    present and not contravened.""")

    # conda
    if conda_native:
        compiler_choices.append("conda")
        compiler_choices.append("Clang")
        compiler_help.append(f"""conda (Clang):
Clang (conda):
    Engage conda-provided (cxx-compiler) clang/clang++ compilers.""")
    else:
        compiler_choices.append(strike("conda"))
        compiler_choices.append(strike("Clang"))
        compiler_help.append(f"""{strike('conda (Clang)')}:
{strike('Clang (conda)')}:
    Can't engage conda-provided clang/clang++ compilers because
    package (cxx-compiler) not installed in current env.""")

elif conda_platform_native == "win-64":
    compiler_help.append(f"""(default: {compiler_default})""")

    # byo
    compiler_choices.append("byo")
    compiler_help.append("""byo:
    Omit compilers from cache to engage self-provided compilers.
    CMake configuration may still use conda compilers if
    present and not contravened.""")

    # conda
    if conda_native:
        compiler_choices.append("conda")
        compiler_choices.append("MSVC")
        compiler_help.append(f"""conda (MSVC):
MSVC (conda):
    Engage conda-provided (cxx-compiler) MSVC compilers.""")
    else:
        compiler_choices.append(strike("conda"))
        compiler_choices.append(strike("MSVC"))
        compiler_help.append(f"""{strike('conda (MSVC)')}:
{strike('MSVC (conda)')}:
    Can't engage conda-provided MSVC compilers because
    package (cxx-compiler) not installed in current env.""")

else:
    raise RuntimeError("unexpected A")


def compiler_type(arg):
    if arg.lower() == "conda":
        return {"linux-64": "GNU", "osx-64": "Clang", "osx-arm64": "Clang", "win-64": "MSVC"}[conda_platform_native]
    for itm in compiler_choices:
        if arg.lower() == itm.lower():
            return itm


#   gcc_linux-64              11.3.0              he6f903b_13    conda-forge
#   gfortran_linux-64         11.3.0              h3c55166_13    conda-forge
#   gxx_linux-64              11.3.0              hc203a17_13    conda-forge
#
#   clang_osx-64              12.0.1               h633439f_9    conda-forge
#   clangxx_osx-64            12.0.1               hdb584c0_9    conda-forge
#
#   clang_osx-arm64           15.0.7               h77e971b_3    conda-forge
#   clangxx_osx-arm64         15.0.7               h768a7fd_3    conda-forge
#   gfortran_osx-arm64        12.2.0               h57527a5_1    conda-forge
#
#   vs2019_win-64
#
#   intel::dpcpp_linux-64     2023.0.0            intel_25370    intel
#   dpcpp_linux-64            2023.2.0                  49495    conda-forge


class PreserveWhiteSpaceWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    # thanks, https://stackoverflow.com/a/35925919

    def __add_whitespace(self, idx, iWSpace, text):
        if idx == 0:
            return text
        return (" " * iWSpace) + text

    def _split_lines(self, text, width):
        textRows = text.splitlines()
        for idx,line in enumerate(textRows):
            search = re.search(r'\s*[0-9\-]{0,}\.?\s*', line)
            if line.strip() == "":
                textRows[idx] = " "
            elif search:
                lWSpace = search.end()
                lines = [self.__add_whitespace(i,lWSpace,x) for i,x in enumerate(_textwrap.wrap(line, 70))]  # width))]
                textRows[idx] = lines

        return [item for sublist in textRows for item in sublist]


##### MENU

parser = argparse.ArgumentParser(
    prog="psi4-path-advisor",
    formatter_class=PreserveWhiteSpaceWrapRawTextHelpFormatter,
    description="""Dependency, Build, and Run path advisor for Psi4.
Mediates file https://github.com/psi4/psi4/blob/master/codedeps.yaml
Run env subcommand. Conda env create and activate. Run cmake subcommand. Build.

=========================================
  (A) black-box usage (copy/paste-able)
=========================================
# (1) get code from GitHub
git clone https://github.com/psi4/psi4.git && cd psi4
# (2) generate env spec file from codedeps.yaml. "eval $(...)" creates and activates conda env.
eval $(conda/psi4-path-advisor.py env)
# (3) generate cmake cache file from conda env. "eval $(...)" configures and builds with cmake.
eval $(conda/psi4-path-advisor.py cmake)

# (4) next steps. repeat each login or add to shell's rc. your paths may vary.
eval $(objdir_p4dev/stage/bin/psi4 --psiapi)
export PSI_SCRATCH=/path/to/existing/writable/local-not-network/directory/for/scratch/files

=========================================
  (B) flexible usage
=========================================
# (1) get code from GitHub
git clone https://github.com/psi4/psi4.git && cd psi4

# (2.0) consider dependency options
conda/psi4-path-advisor.py env -h
# (2.1) generate env spec file from codedeps.yaml.
conda/psi4-path-advisor.py env -n p4dev310 --python 3.10 --disable addons --lapack openblas
#> conda env create -n p4dev310 -f /home/psi4/env_p4dev310.yaml && conda activate p4dev310
# (2.2) edit env_p4dev310.yaml to customize software packages.
# (2.3) issue suggested or customized command to create and activate conda env.
conda env create -n p4dev310 -f /home/psi4/env_p4dev310.yaml && conda activate p4dev310

# (3.0) consider compile options
conda/psi4-path-advisor.py cmake -h
# (3.1) generate cmake cache file from conda env.
conda/psi4-path-advisor.py cmake
#> cmake -S. -GNinja -C/home/psi4/cache_p4dev310.cmake -Bobjdir_p4dev310 && cmake --build objdir_p4dev310
# (3.2) edit cache_p4dev310.cmake to customize build configuration.
# (3.3) issue suggested or customized command to configure & build with cmake.
cmake -S. -GNinja -C/home/psi4/cache_p4dev310.cmake -Bobjdir_p4dev310 -DCMAKE_INSTALL_PREFIX=/path/to/install-psi4 && cmake --build objdir_p4dev310

# (4) next steps. repeat each login or add to shell's rc. your paths may vary.
eval $(objdir_p4dev310/stage/bin/psi4 --psiapi)
export PSI_SCRATCH=/path/to/existing/writable/local-not-network/directory/for/scratch/files
""")

parser.add_argument("-v", action="count", default=0,
    help="""Use for more printing (-vv).
Verbosity arg is NOT compatible with bash command substitution.""")

subparsers = parser.add_subparsers(dest="subparser_name",
    help="Script requires a subcommand.")

parser_env = subparsers.add_parser("env",
    aliases=["conda"],
    formatter_class=PreserveWhiteSpaceWrapRawTextHelpFormatter,
    help="Write conda environment file from codedeps file.")
parser_env.add_argument("--lapack",
    choices=["mkl", "openblas", "accelerate", "blis", "netlib"],
    help=f"""(default: mkl if available else accelerate)
Specify a blas/lapack version.
'accelerate' only for osx-* platforms.
'mkl' and 'blis' only for *-64 platforms.""")
parser_env.add_argument("-n", "--name",
    default="p4dev",
    help=f"""Specify environment name.
Can instead edit generated env spec file by hand
-or- rename on the cmdline, so argument mostly for printing.""")
parser_env.add_argument("--python",
    default="",
    help="""Specify a python version.
Can instead accept latest or edit env spec file by hand.""")
parser_env.add_argument("--disable",
    nargs="+",
    choices=["compilers", "lapack", "addons", "test", "docs"],
    help=f"""Categories of dependencies to not include in env file.
Can instead edit generated env spec file by hand.
For example, '--disable compilers lapack addons test docs' for minimal.""")
parser_env.add_argument("--solver",
    default="conda-else-mamba",
    choices=solver_choices,
    help="\n".join(solver_help))
parser_env.add_argument("--platform",
    choices=["linux-64", "osx-64", "osx-arm64", "win-64"],
    default=conda_platform_native,
    help=f"""Conda platform/subdir for env file, if not the computed native
(default: {conda_platform_native}).
Apple Silicon users, check this value! Argument rarely used.""")
parser_env.add_argument("--offline-conda",
    action="store_true",
    help=f"""Use script without conda/mamba available.
Useful for bootstrapping CI runners. Offline arg is NOT compatible with bash command substitution.""")

# add constraints to env file?

parser_cmake = subparsers.add_parser("cache",
    aliases=["cmake"],
    formatter_class=PreserveWhiteSpaceWrapRawTextHelpFormatter,
    help="Write cmake configuration cache from conda environment.")
parser_cmake.add_argument("--compiler",
    choices=compiler_choices,
    default=compiler_default,
    type=compiler_type,
    help="\n".join(compiler_help))
parser_cmake.add_argument("--lapack",
    choices=lapack_choices,
    default=lapack_default,
    type=lapack_type,
    help="\n".join(lapack_help))
parser_cmake.add_argument("--objdir",
    default=f"objdir_{conda_prefix_short}",
    help=f"""Specify a build directory to cmake, if not (default: objdir_{conda_prefix_short}).
Can instead rename on cmdline, so argument mostly for printing.""")
parser_cmake.add_argument("--insist",
    action="store_true",
    help=f"""Set the cache (`INSIST_FIND_PACKAGE_<pkg>=ON`) to prevent cmake from falling back on internal build for packages present in the conda environment.""")

parser_bulletin = subparsers.add_parser("bulletin",
    formatter_class=PreserveWhiteSpaceWrapRawTextHelpFormatter,
    help="(Info only) Read any build messages or advice.")

parser_deploy = subparsers.add_parser("deploy",
    formatter_class=PreserveWhiteSpaceWrapRawTextHelpFormatter,
    help="(Admin only) Apply codedeps info to codebase.")

args = parser.parse_args()

if args.subparser_name in ["conda", "env"]:
    conda_platform = args.platform
else:
    conda_platform = conda_platform_native
conda_unix = "win-64" if conda_platform == "win-64" else "unix"


if args.v > 1:
    print("#######")
    parser.print_help()
    print("#######")
    parser_env.print_help()
    print("#######")
    parser_cmake.print_help()
    print("#######")
    parser_bulletin.print_help()
    print("#######")
    parser_deploy.print_help()
    print("#######")
    print(f"{conda_platform=}  {conda_lapack_variant=}  {args=}")
    print("#######")

#            # pip
#          - pip
#          - pip:
#             - git+https://github.com/i-pi/i-pi.git@master-py3
#""")

with codedeps_yaml.open() as fp:
    ydict = yaml_load(fp)

##### CONDA ENV

if args.subparser_name in ["conda", "env"]:
    if not (pm_available or args.offline_conda):
        raise RuntimeError("usage: this script requires either the conda or mamba command to be in envvar PATH.")

    if args.solver == "conda-else-mamba":
        if conda_available:
            solver = "conda"
        elif mamba_available:
            solver = "mamba"
        elif args.offline_conda:
            solver = "(no pm)"
    else:
        solver = args.solver
    if solver == "conda-libmamba":
        solver = ["conda", "--solver libmamba"]
    elif solver == "conda-classic":
        solver = ["conda", "--solver classic"]
    else:
        solver = [solver, ""]

    stuff = {
        "build": [],
        "non-qc buildtime required": [],
        "qc buildtime required": [],
        "buildtime optional": [],
        "runtime required": [],
        "runtime optional": [],
        "test": [],
        "docs": [],
    }
    notes = {}
    lapack_packages = []  # TODO: incl openmp, too?

    for ddep in ydict["data"]:

        # skips
        use = ddep["use"]
        if "deprecated" in use:
            continue

        conda = ddep["conda"]
        if not conda:
            continue

        if conda_platform == "win-64" and conda.get("skip_win", False):
            continue

        # start collecting
        aux_run = conda.get("aux_run_names", [])
        aux_bld = conda.get("aux_build_names", [])
        if isinstance(aux_run, dict):
            aux_run = aux_run .get(conda_platform, [])
        if isinstance(aux_bld, dict):
            aux_bld = aux_bld.get(conda_platform, [])

        if isinstance(conda["name"], dict):
            primary = conda["name"][conda_platform]
        else:
            primary = conda["name"]

        if primary == "libblas":
            if args.lapack:
                if args.lapack == "accelerate" and not conda_platform.startswith("osx-"):
                    raise RuntimeError("libblas accelerate only available for platforms osx-[64|arm64].")
                if args.lapack in ["mkl", "blis"] and not conda_platform.endswith("-64"):
                    raise RuntimeError("libblas mkl and blis only available for platforms [linux|osx|win]-64.")
                lapack_constraint = f"=*=*{args.lapack}"

                if args.lapack == "openblas" and conda_platform == "linux-64":
                    aux_bld.append("openblas=*=openmp*")
            else:
                lapack_constraint = conda["constraint"][conda_platform]
            primary += lapack_constraint

            lapack_packages.extend([primary, *aux_bld, *aux_run])

        if conda["channel"] != "conda-forge":
            primary = conda["channel"] + "::" + primary

        if cons := conda["constraint"]:
            if isinstance(cons, str):
                primary += cons

        for pkg in aux_run:
            req = "opt'l" if pkg.startswith("//") else "req'd"
            notes[pkg] = f"{req} with {conda['name']}"
            if note := conda.get("aux_run_names_note", {}).get(pkg, None):
                notes[pkg] += f"; {note}"

        for pkg in aux_bld:
            req = "opt'l" if pkg.startswith("//") else "req'd"
            notes[pkg] = f"{req} with {conda['name']}"
            if note := conda.get("aux_build_names_note", {}).get(pkg, None):
                notes[pkg] += f"; {note}"

        # sort into categories
        if use.get("test_required", None) is not None:
            stuff["test"].append(primary)
            stuff["test"].extend(aux_bld)
            stuff["test"].extend(aux_run)
        elif use.get("docs_required", None) is not None:
            stuff["docs"].append(primary)
            stuff["docs"].extend(aux_bld)
            stuff["docs"].extend(aux_run)
        elif use["buildtime"]:
            if use["cms"]:
                if use["required"]:
                    stuff["qc buildtime required"].append(primary)
                    bareprimary = re.match(re_pkgline, primary).group("pkg")
                    for itm in aux_bld:
                        if re.match(re_pkgline, itm).group("pkg").startswith(bareprimary):
                            # pretty much only to sort the alternate libints to the same section
                            stuff["qc buildtime required"].append(itm)
                        else:
                            stuff["non-qc buildtime required"].append(itm)
                    stuff["runtime required"].extend(aux_run)
                else:
                    stuff["buildtime optional"].append(primary)
                    stuff["buildtime optional"].extend(aux_bld)
                    stuff["runtime optional"].extend(aux_run)
            else:
                if primary in ["cmake", "c-compiler", "cxx-compiler", "fortran-compiler", "ninja"] and use["required"]:
                    stuff["build"].append(primary)
                    stuff["build"].extend(aux_bld)
                    if aux_run:
                        raise RuntimeError("unexpected B")
                elif use["required"]:
                    stuff["non-qc buildtime required"].append(primary)
                    stuff["non-qc buildtime required"].extend(aux_bld)
                    stuff["runtime required"].extend(aux_run)
                else:
                    stuff["buildtime optional"].append(primary)
                    stuff["buildtime optional"].extend(aux_bld)
                    stuff["runtime optional"].extend(aux_run)
        else:
            if use["cms"]:
                if use["required"]:
                    raise RuntimeError("unexpected C")
                else:
                    stuff["runtime optional"].append(primary)
                    stuff["runtime optional"].extend(aux_run)
            else:
                stuff["runtime required"].append(primary)
                stuff["runtime required"].extend(aux_run)

    stuff = {k: sorted(v) for k, v in stuff.items()}

    text = []
    text.append(f"name: {args.name}")
    text.append("channels:")
    text.append("  - conda-forge")
    text.append("  - nodefaults")
    text.append("dependencies:")

    if args.disable:
        if "compilers" in args.disable:
            for category in list(stuff):
                for pkg in list(stuff[category]):
                    if "compiler" in pkg or "dpcpp_linux-64" in pkg:
                        stuff[category].remove(pkg)
        if "lapack" in args.disable:
            for category in list(stuff):
                for pkg in list(stuff[category]):
                    if pkg in lapack_packages:
                        stuff[category].remove(pkg)
        if "addons" in args.disable:
            stuff.pop("buildtime optional")
            stuff.pop("runtime optional")
        if "test" in args.disable:
            stuff.pop("test")
        if "docs" in args.disable:
            stuff.pop("docs")

    for category in stuff:
        text.append(f"    # {category}")
        for pkg in sorted(stuff[category], key=lambda x: re.match(re_pkgline, x).group("pkg")):
            commentout, chnl, barepkg, constraint = re.match(re_pkgline, pkg).groups()

            if note := notes.get(pkg, ""):
                if commentout:
                    text.append(f"  #- {pkg[2:]:<24}  # {note}")
                else:
                    text.append(f"  - {pkg:<24}  # {note}")
            else:
                if barepkg == "python" and args.python:
                    text.append(f"  - {pkg}={args.python}")
                elif commentout:
                    text.append(f"  #- {pkg[2:]}")
                else:
                    text.append(f"  - {pkg}")

    text.append("")
    condaenvspec = Path(f"env_{args.name}.yaml").resolve()
    with condaenvspec.open("w") as fp:
        fp.write("\n".join(text))
    if args.v > 0:
        print("\n".join(text))

    subdircmd = "" if conda_platform == conda_platform_native else f"CONDA_SUBDIR={conda_platform} "
    env_create_and_activate_cmd = f"""{subdircmd}{solver[0]} env create -n {args.name} -f {condaenvspec} {solver[1]} && conda activate {args.name}"""
    print(env_create_and_activate_cmd)


##### CMAKE CACHE

elif args.subparser_name in ["cmake", "cache"]:
    if not pm_available:
        raise RuntimeError("usage: this script requires either the conda or mamba command to be in envvar PATH.")

    text = []
    pyotf_merge = {}
    pyotf = "Primarily OTF runtime detected."
    cmake_program_path = []
    cpp_merge = []
    absent_merge = []
    big_args = {
        "-S": cmake_S,
      # "-G"
      # "-C"
        "-B": args.objdir,
    }
    shlib_ext = {
        "linux-64": ".so",
        "osx-64": ".dylib",
        "osx-arm64": ".dylib",
        "win-64": ".dll",
    }

    for ddep in ydict["data"]:

        use = ddep["use"]
        if "deprecated" in use:
            continue
        conda = ddep["conda"]
        if not conda:
            continue

        if isinstance(conda["name"], dict):
            primary = conda["name"][conda_platform]
        else:
            primary = conda["name"]
        aux_run = conda.get("aux_run_names", [])
        aux_bld = conda.get("aux_build_names", [])
        if isinstance(aux_run, dict):
            aux_run = aux_run .get(conda_platform, [])
        if isinstance(aux_bld, dict):
            aux_bld = aux_bld.get(conda_platform, [])

        constraint_delimiter = "=|!|<|>"
        package_set = [primary, *aux_bld, *aux_run]
        plain_package_set = [re.split(constraint_delimiter, pkg)[0] for pkg in package_set]
        plain_package_set = [pkg for pkg in plain_package_set if not pkg.startswith("//")]

        if all(pkg in conda_list_pkgver for pkg in plain_package_set):
            all_found = True
            pkgstr = plain_package_set
            primary_banner = f"""# <<<  {primary:<27} {", ".join(plain_package_set[1:])}"""
        else:
            all_found = False
            pkgstr = [pkg if (pkg in conda_list_pkgver) else strike(pkg, tilde=True) for pkg in plain_package_set]
            primary_banner = f"""# <<<  {pkgstr[0]:<{compute_width(pkgstr[0], 27)}} {", ".join(pkgstr[1:])}"""
            if not use["required"]:
                absent_merge.append(primary_banner)
                continue

        if primary == "libblas" and args.lapack == "openblas" and conda_platform == "linux-64":
            if conda_openblas_variant != "openmp":
                all_found = False
                pkgstr.append(strike("openblas=*=openmp*", tilde=True))
                primary_banner = f"""# <<<  {pkgstr[0]:<{compute_width(pkgstr[0], 27)}} {", ".join(pkgstr[1:])}"""

        if primary == "libblas" and args.lapack == "mkl" and not all_found:
            # since it's been recc for so long, allow mkl-devel pkg instead of blas-devel
            plain_package_set[:] = ["mkl-devel" if (p == "blas-devel") else p for p in plain_package_set]
            all_found = all(pkg in conda_list_pkgver for pkg in plain_package_set)

        if conda["cmake"] and list(conda["cmake"].keys()) == ["CMAKE_PROGRAM_PATH"]:
            if not all_found:
                raise RuntimeError("required and program_path not handled")
            cmake_program_path.append(conda["cmake"]["CMAKE_PROGRAM_PATH"].replace("${CONDA_PREFIX}", conda_prefix))
            cpp_merge.append(f"""# <<<  {primary:<27} {conda.get("cmake_note", "")}""")

        elif conda["cmake"]:
            text.append("")
            text.append(primary_banner)
            if note := conda.get("cmake_note"):
                text.append(f"# {note}")

            if primary in ["c-compiler", "cxx-compiler", "fortran-compiler"]:
                if args.compiler == "byo":
                    text.append("# Bring-your-own (byo) compilers by setting here -or- by passing on cmdline -or- by letting cmake autodetect.")
                dcmake_vars = conda["cmake"].get(f"{args.compiler}_{conda_platform}", conda["cmake"].get(args.compiler))
            elif primary in ["libblas"]:
                if args.lapack == "byo":
                    text.append("# Bring-your-own (byo) blas/lapack libraries by setting here -or- by passing on cmdline -or- by letting cmake autodetect.")
                    text.append("#   Note that mixing lapack implementations is not advised.")
                dcmake_vars = conda["cmake"].get(f"{args.lapack}_{conda_platform}", conda["cmake"].get(args.lapack))
                if dcmake_vars is None:
                    print(f"Libblas info missing in keys '{args.lapack}_{conda_platform}' or '{args.lapack}'. Contact the developers.")
                    sys.exit(4)
            else:
                dcmake_vars = conda["cmake"]

            for k, v in dcmake_vars.items():

                if isinstance(v, dict):
                    v = v.get(conda_platform, v.get(conda_unix))

                if k.endswith("_DIR"):
                    ctyp = "PATH"
                elif isinstance(v, bool):
                    ctyp = "BOOL"
                    v = "ON" if v else "OFF"
                else:
                    ctyp = "STRING"

                if k == "CMAKE_PROGRAM_PATH":
                    raise ValueError("more than one var", dcmake_vars)

                if k.startswith("//"):
                    if args.insist and "INSIST_FIND_PACKAGE" in k and all_found:
                        text.append(f'set({k[2:]:<28} ON CACHE {ctyp} "")')
                    else:
                        text.append(f'# set({k[2:]:<28} {"<placeholder>" if v is None else v} CACHE {ctyp} "")')
                elif not all_found:
                    text.append(f'# set({k:<28} {"<placeholder>"} CACHE {ctyp} "")')
                elif k.startswith("-"):
                    big_args[k] = v.replace("<src>", cmake_S).replace("<bld>", args.objdir)
                else:
                    v = v.replace("${CONDA_PREFIX}", conda_prefix)
                    v = v.replace("${SHLIB_EXT}", shlib_ext[conda_platform])
                    if "${HOST}" in v:
                        v = v.replace("${HOST}", conda_host)
                    # For PATH variables, check if the path exists before adding to cache
                    # This prevents pre-setting *_DIR variables to non-existent paths which breaks find_package
                    if ctyp == "PATH" and not Path(v).exists():
                        text.append(f'# set({k:<28} {v} CACHE {ctyp} "")  # Path does not exist')
                    else:
                        text.append(f'set({k:<30} {v} CACHE {ctyp} "")')
                        if not (ctyp == "BOOL" or Path(v).exists()):
                            pass
                            # print(f"Warning: Active value in cache for {k} does not exist on filesystem: {v}")

        else:
            if note := conda.get("cmake_note"):
                if note == pyotf or (note.startswith(pyotf) and "inconsequential" in note):
                    pyotf_merge[primary] = note
                else:
                    text.append("")
                    text.append(primary_banner)
                    text.append(f"# {note}")
            else:
                text.append("")
                text.append(primary_banner)
                text.append("# NOTE TODO no cmake hook")

    text.append("")
    text.append("# Sections merged because detection hinted through a common list variable")
    text.extend(cpp_merge)
    if cpp_merge:
        text.append(f'set({"CMAKE_PROGRAM_PATH":<30} {";".join(set(cmake_program_path))} CACHE STRING "")')
    else:
        text.append(f'# set({"CMAKE_PROGRAM_PATH":<28} <placeholder> CACHE STRING "")')

    text.append("")
    text.append("# Sections skipped because packages detected on-the-fly, not through CMake")
    for primary, msg in pyotf_merge.items():
        text.append(f"# <<<  {primary:<27} {msg}")

    text.append("")
    text.append("# Sections skipped because packages absent from current environment")
    text.extend(absent_merge)

    objdir_lbl = Path(args.objdir).stem.replace("objdir_", "")
    cmakecache_lbl = conda_prefix_short if (objdir_lbl == conda_prefix_short) else f"{conda_prefix_short}@{objdir_lbl}"
    cmakecache_lbl = f"cache_{cmakecache_lbl}.cmake"
    cmakecache = Path(cmakecache_lbl).resolve()
    big_args["-C"] = rf'"{cmakecache}"'
    big_args_sorted = [f"{k}{big_args[k]}" for k in sorted(big_args.keys(), reverse=True)]

    pretext = f"""# {cmakecache_lbl}
# {'-' * len(cmakecache_lbl)}
#
# This file is autogenerated by psi4-path-advisor from conda env {conda_prefix_short}.
#   It sets some likely values to initialize the CMake cache for dependencies
#   to build your Psi4 source.
#
# Dependency packages are shown as \"present-package\" or \"{strike('missing-package', tilde=True)}\".
#   Psi4 is not assured to build if any required dependencies show as missing.
#   Feel free to edit this file prior to CMake configuration or to override with
#   \"-D\" arguments on the command line. Note that for most setups, no cache file
#   is needed for CMake configuration to detect all the dependencies in an
#   active conda environment.
#
# An example command line usage is below:
#
# > conda activate {conda_prefix_short}
# > cmake \\
"""
    # using tilde-type strike in cache file so searchable

    for itm in big_args_sorted:
        pretext += "#    " + itm
        pretext += "\n" if itm == big_args_sorted[-1] else " \\\n"
    pretext += f"# > cmake --build {big_args['-B']}"

    pretext = pretext.splitlines()
    pretext.extend(text)
    text = pretext

    with cmakecache.open("w") as fp:
        fp.write("\n".join(text))
    if args.v > 0:
        print("\n".join(text))

    cmake_configure_and_build_cmd = "cmake " + " ".join(big_args_sorted) + f" && cmake --build {big_args['-B']}"
    print(cmake_configure_and_build_cmd)


##### BULLETIN

elif args.subparser_name in ["bulletin"]:

    notes = f"""

  * [21 Dec 2023] The Psi4 source can use most any Libint packaged as psi4*::libint2
    in 2023, as conda-forge/label/libint_dev::libint in 2023 or
    conda-forge::libint>=2.8.0 . They may not detect the same with CMake, but at
    runtime, as far as Psi4 uses the library, they're the same. You are encouraged
    to upgrade to v2.8.1 soon, as the master branch will start to require it in
    early 2024. AM5 production builds are available from conda-forge, and an AM7
    is available for Linux from psi4. Henceforth, all packages will be named
    "libint", not "libint2".

  * [15 Nov 2023] If ever you are building from git source and the version is
    undefined (shows up in CMake output, too), the solution is to pull tags:
    `git fetch upstream "refs/tags/*:refs/tags/*"`, then recompile (trivial).

      psi4 --version
      #> undefined
      psi4
      #>  Traceback (most recent call last):
      #>    File "/path/to/psi4/objdir/stage/bin/psi4", line 232, in <module>
      #>      import psi4  # isort:skip
      #>      ^^^^^^^^^^^
      #>    File "/path/to/psi4/objdir/stage/lib/psi4/__init__.py", line 71, in <module>
      #>      from . import core
      #>  ImportError: initialization failed


  * [15 Nov 2023] The error below is known to happen with the
    `subdirectory::package` syntax. It has diagnosed itself correctly, so you may
    as well update libmamba. Alternately, you can add the subdirectory to the
    environment list at the top of the env spec file and delete it and the :: from
    the package line.

      InvalidMatchSpec: Invalid spec 'conda-forge/label/libint_dev::libint==2.7.3dev1': This is a bug in libmamba 1.5.1 when using 'defaults::<spec>' or 'pkgs/main::<spec>'. Consider using '-c defaults' instead.


  * [11 Nov 2023] Around July 2023 miniconda started shipping mamba-ready (though
    not as the default solver). Likewise, miniforge and mambaforge started making
    the conda command available. Around October 2023, miniconda started shipping
    with mamba as the default solver. For an older anaconda/miniconda installation,
    if you want to instruct it to use the mamba solver by default from behind the
    conda hood, issue the following from
    https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community :

      conda update -n base conda
      conda install -n base conda-libmamba-solver
      conda config --set solver libmamba

  * [11 Nov 2023] For Silicon (osx-arm64) owners whose miniconda is set to Intel
    (osx-64) (check with `conda info`):

      conda config --system --set subdir osx-arm64
"""
    print(notes)


##### DEPLOY

elif args.subparser_name in ["deploy"]:
    if not pm_available:
        raise RuntimeError("usage: this script requires either the conda or mamba command to be in envvar PATH.")

    # sugg from docs CONDA_OVERRIDE_LINUX=1 and CONDA_OVERRIDE_GLIBC=2.17

    full_cmake_S = codedeps_yaml.parent
    full_conda_envs = f"{full_cmake_S}/devtools/conda-envs/"
    script = f"""#!/usr/bin/env bash

PNAME=p4dev8  # some env that doesn't exist

# Core
{full_cmake_S}/conda/psi4-path-advisor.py env --platform linux-64 --lapack mkl --disable addons docs
CONDA_SUBDIR=linux-64 conda env create -n $PNAME -f env_p4dev.yaml --dry-run
mv env_p4dev.yaml {full_conda_envs}/linux-64-buildrun.yaml

{full_cmake_S}/conda/psi4-path-advisor.py env --platform osx-64 --lapack mkl --disable addons docs
CONDA_SUBDIR=osx-64 conda env create -n $PNAME -f env_p4dev.yaml --dry-run
mv env_p4dev.yaml {full_conda_envs}/osx-64-buildrun.yaml

{full_cmake_S}/conda/psi4-path-advisor.py env --platform osx-arm64 --lapack accelerate --disable addons docs
CONDA_SUBDIR=osx-arm64 conda env create -n $PNAME -f env_p4dev.yaml --dry-run
mv env_p4dev.yaml {full_conda_envs}/osx-arm64-buildrun.yaml

{full_cmake_S}/conda/psi4-path-advisor.py env --platform win-64 --lapack mkl --disable addons docs
CONDA_SUBDIR=win-64 conda env create -n $PNAME -f env_p4dev.yaml --dry-run
mv env_p4dev.yaml {full_conda_envs}/win-64-buildrun.yaml

# Docs
{full_cmake_S}/conda/psi4-path-advisor.py env --name p4docs --platform linux-64 --disable addons
CONDA_SUBDIR=linux-64 conda env create -n $PNAME -f env_p4docs.yaml --dry-run
mv env_p4docs.yaml {full_conda_envs}/linux-64-docs.yaml

# Eco
{full_cmake_S}/conda/psi4-path-advisor.py env --platform linux-64 --lapack mkl --disable docs
CONDA_SUBDIR=linux-64 conda env create -n $PNAME -f env_p4dev.yaml --dry-run
mv env_p4dev.yaml {full_conda_envs}/linux-64-buildrun-addons.yaml

{full_cmake_S}/conda/psi4-path-advisor.py env --platform osx-64 --lapack mkl --disable docs
CONDA_OVERRIDE_OSX=13 CONDA_SUBDIR=osx-64 conda env create -n $PNAME -f env_p4dev.yaml --dry-run
mv env_p4dev.yaml {full_conda_envs}/osx-64-buildrun-addons.yaml

{full_cmake_S}/conda/psi4-path-advisor.py env --platform osx-arm64 --lapack accelerate --disable docs
CONDA_OVERRIDE_OSX=13 CONDA_SUBDIR=osx-arm64 conda env create -n $PNAME -f env_p4dev.yaml --dry-run
mv env_p4dev.yaml {full_conda_envs}/osx-arm64-buildrun-addons.yaml

{full_cmake_S}/conda/psi4-path-advisor.py env --platform win-64 --lapack mkl --disable docs
CONDA_SUBDIR=win-64 conda env create -n $PNAME -f env_p4dev.yaml --dry-run
mv env_p4dev.yaml {full_conda_envs}/win-64-buildrun-addons.yaml

"""

    with open("deps_deploy_devtools.sh", "w") as fp:
        fp.write(script)
        print("bash deps_deploy_devtools.sh")

    seds = []
    for ddep in ydict["data"]:
        repo = ddep["repository"]
        if repo is None:
            continue

        comment = f"  # {repo['commit_note']}" if ("commit_note" in repo) else ""
        if repo["host"] == "github":
            repo_url = f"https://github.com/{repo['account']}/{repo['name']}"
            if repo.get("githttps"):
                url = f"{repo_url}.git@{repo['commit']}#egg=proj"
            else:
                url = f"{repo_url}/archive/{repo['commit']}.tar.gz{comment}"
        elif repo["host"] == "gitlab":
            repo_url = f"https://gitlab.com/{repo['account']}/{repo['name']}"
            url = f"{repo_url}/-/archive/{repo['commit']}/{repo['name']}-{repo['commit']}.tar.gz{comment}"
        elif repo["host"] == "url":
            repo_url = "/".join(repo["url"].split("/", 3)[:-1])
            url = repo["url"]

        seds.append(f"sed -i 's;{repo_url}.*;{url}  # edit in codedeps;' {full_cmake_S}/external/*/*/CMakeLists.txt")

    for ddep in ydict["data"]:
        cm = ddep["cmake"]
        if cm is None:
            continue

        components = (" COMPONENTS " + " ".join(cm["components"])) if cm.get("components", False) else ""
        constraint = f" {cm['constraint']}" if cm.get("constraint", False) else ""
        fp = cm['name'] + constraint + components + " "

        seds.append(f"""sed -i 's;find_python_module({cm["name"]} .*QUIET;find_python_module({fp}QUIET;' {full_cmake_S}/external/*/*/CMakeLists.txt""")
        seds.append(f"""sed -i 's;find_package({cm["name"]} .*CONFIG;find_package({fp}CONFIG;' {full_cmake_S}/external/*/*/CMakeLists.txt""")
        seds.append(f"""sed -i 's;find_package({cm["name"]} .*CONFIG;find_package({fp}CONFIG;' {full_cmake_S}/psi4/CMakeLists.txt""")

    with open("deps_deploy_external.sh", "w") as fp:
        fp.write("\n".join(seds))
        print("bash deps_deploy_external.sh")


#elif sys.platform == 'darwin':
#    parser.add_argument('--clang', action='store_true',
#                        help="""Engage conda's psi4-dev-provided clang/clang++/gfortran compilers. You must have downloaded this file https://github.com/phracker/MacOSX-SDKs/releases/download/10.13/MacOSX10.9.sdk.tar.xz, unpacked it, and saved it at ~/SDKs/MacOSX10.9.sdk . !Change! this arg invoked XCode AppleClang prior to Jul 2018.""")
#
#    #                    help="""Engage system-provided clang/clang++ compilers and psi4-dev-provided gfortran.""")
#    #parser.add_argument('--gcc', action='store_true',
#    #                    help="""Engage psi4-dev-provided gcc/g++/gfortran compilers.""")
#
## duplicates from `bin/psi4`
#psi4 = os.path.abspath(os.path.dirname(__file__)) + os.path.sep + 'psi4'
#psi4alongside = os.path.isfile(psi4) and os.access(psi4, os.X_OK)
#
#if psi4alongside:
#    parser.add_argument("--psiapi-path", action='store_true',
#                        help="""(Duplicate from `psi4`) Generates a bash command to source correct Python for `python -c "import psi4"`""")
#
#if psi4alongside:
#    from subprocess import call
#
#    if args.psiapi_path:
#        call([psi4, '--psiapi-path'])
#        sys.exit(0)
#
#else:
#    if args.plugin_compile:
#        print("""Install "psi4" via `conda install psi4 -c psi4[/label/dev]`, then reissue command.""")
#
#
#if sys.platform == 'darwin':
#    if args.clang:
#        recc.insert(-1, '-C/opt/anaconda1anaconda2anaconda3/share/cmake/psi4/psi4DepsClangCache.cmake')
#    #    recc.insert(0, 'CONDA_BUILD_SYSROOT=~/SDKs/MacOSX10.9.sdk')
#    #    recc.insert(-1, '-C/opt/anaconda1anaconda2anaconda3/share/cmake/psi4/psi4DepsAppleClangCache.cmake')
#    #if args.gcc:
#    #    recc.insert(-1, '-C/opt/anaconda1anaconda2anaconda3/share/cmake/psi4/psi4DepsGNUCache.cmake')
