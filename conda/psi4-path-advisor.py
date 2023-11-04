#!/usr/bin/env python

import os
import re
import json
import yaml
import shutil
import argparse
import itertools
import unicodedata
import textwrap as _textwrap
from pathlib import Path
from subprocess import run


codedeps_yaml = Path(__file__).parent.parent / "codedeps.yaml"
cmake_S = os.path.relpath(codedeps_yaml.parent, start=Path.cwd())


def conda_list(*, name: str = None, prefix: str = None) -> str:
    # thanks, https://stackoverflow.com/a/56363822
    #   SO convinced me that subprocess was better than import (following block)
    #   Also, the conda.cli only works on base env (where conda pkg installed)

    #import conda.cli.python_api as Conda
    #env_list_json, stderr, rc = Conda.run_command(Conda.Commands.LIST, ["--json"])
    #env_list_dict = json.loads(env_list_json)

    if name:
        proc = run(["conda", "list", "--json", "--name", name], text=True, capture_output=True)
    elif prefix:
        proc = run(["conda", "list", "--json", "--prefix", prefix], text=True, capture_output=True)
    else:
        proc = run(["conda", "list", "--json"], text=True, capture_output=True)
    return json.loads(proc.stdout)


def conda_info():
    proc = run(["conda", "info", "--json"], text=True, capture_output=True)
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
if not conda_available or mamba_available:
    raise RuntimeError("usage: this script requires either the conda or mamba command to be in envvar PATH.")

conda_info_dict = conda_info()
conda_platform_native = conda_info_dict["platform"]
conda_prefix = conda_info_dict["active_prefix"]
conda_prefix_short = conda_info_dict["active_prefix_name"]
conda_host = conda_info_dict["env_vars"].get("CONDA_TOOLCHAIN_HOST", None)  # None if no compilers in env
conda_list_struct = conda_list()
conda_list_pkgver = {item["name"]: item["version"] for item in conda_list_struct}
conda_lapack_variant = None  # None if no c-f libblas in env
for itm in conda_list_struct:
    if itm["name"] == "libblas":
        conda_lapack_variant = itm["build_string"].split("_")[-1]
    if itm["name"] == "openblas":
        conda_openblas_variant = itm["build_string"].split("_")[0]

# TODO handle conda_host None in base env with no compilers present
# TODO handle conda_lapack_variant None in base env with no lapack present

base_prefix = conda_info_dict["conda_prefix"]  # env with conda cmd
base_list_struct = conda_list(prefix=base_prefix)
conda_libmamba_available = False
for itm in base_list_struct:
    if itm["name"] == "conda-libmamba-solver":
        conda_libmamba_available = True
        break

### conda/mamba

solver_choices = ["mamba-as-possible"]
solver_help = [f"""(default: mamba-as-possible)
mamba-as-possible:
    Use conda-libmamba else mamba else conda to solve environments.
    Can instead adjust on cmdline, so argument mostly for printing."""]
if conda_libmamba_available:
    solver_choices.append("conda-libmamba")
    solver_help.append("""conda-libmamba:
    Use `conda ... --solver=libmamba` to solve environments.""")
else:
    solver_choices.append(strike("conda-libmamba"))
    solver_help.append(f"""{strike('conda-libmamba')}
    Can't use `conda ... --solver=libmamba` to solve environments
    because packages (conda or conda-libmamba-solver) not installed in base env.""")
if mamba_available:
    solver_choices.append("mamba")
    solver_help.append("""mamba:
    Use `mamba` to solve environments. UNTESTED""")
else:
    solver_choices.append(strike("mamba"))
    solver_help.append(f"""{strike('mamba')}
    Can't use `mamba` to solve environments
    because packages (mamba?) not installed in base env.""")
if conda_available:
    solver_choices.append("conda")
    solver_help.append("""conda:
    Use `conda` to solve environments.""")
else:
    solver_choices.append(strike("conda"))
    solver_help.append(f"""{strike('conda')}
    Can't use `conda` to solve environments
    because package (conda) not installed in base env.""")


### blas/lapack

lapack_conda_native = "blas-devel" in conda_list_pkgver
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
                lines = [self.__add_whitespace(i,lWSpace,x) for i,x in enumerate(_textwrap.wrap(line, width))]
                textRows[idx] = lines

        return [item for sublist in textRows for item in sublist]


##### MENU

parser = argparse.ArgumentParser(
    prog="psi4-path-advisor",
    formatter_class=PreserveWhiteSpaceWrapRawTextHelpFormatter,
    description="""Dependency, Build, and Run path advisor for Psi4.
Mediates file https://github.com/psi4/psi4/blob/master/codedeps.yaml
Run env subcommand. Conda env create and activate. Run cmake subcommand. Build.

# =======================================
#  (A) black-box usage (copy/paste-able)
# =======================================
# (1) get code from GitHub
git clone https://github.com/psi4/psi4.git && cd psi4
# (2) generate env spec file from codedeps.yaml. "eval $(...)" creates and activates conda env.
eval $(conda/psi4-path-advisor.py env)
# (3) generate cmake cache file from conda env. "eval $(...)" configures and builds with cmake.
eval $(conda/psi4-path-advisor.py cmake)

shows up in p4dev

# =======================================
#  (B) flexible usage
# =======================================
# (1) get code from GitHub
git clone https://github.com/psi4/psi4.git && cd psi4

# (2) generate env spec file from codedeps.yaml.
conda/psi4-path-advisor.py env -n p4dev310 --python 3.10 --disable addons --lapack openblas
#> conda env create -n p4dev310 -f /home/psi4/env_p4dev310.yaml --solver libmamba && conda activate p4dev310
# (3) edit env_p4dev310.yaml to customize software packages.
# (4) issue suggested or customized command to create and activate conda env.
conda env create -n p4dev310 -f /home/psi4/env_p4dev310.yaml --solver libmamba && conda activate p4dev310

# (5) generate cmake cache file from conda env.
./conda/psi4-path-advisor.py cmake
#> cmake -S. -GNinja -C/home/psi4/cache_p4dev310.cmake -Bobjdir_p4dev310 && cmake --build objdir_p4dev310
# (6) edit cache_p4dev310.cmake to customize build configuration.
# (7) issue suggested or customized command to configure and build with cmake.
cmake -S. -GNinja -C/home/psi4/cache_p4dev310.cmake -Bobjdir_p4dev310 -DCMAKE_INSTALL_PREFIX=/path/to/install-psi4 && cmake --build objdir_p4dev310

conda/psi4-path-advisor.py env -n p4dev310 --python 3.10 --disable addons --lapack openblas
> psi4-path-advisor.py
eval $(conda/psi4-path-advisor.py env -n ppa_plain)
#>>> conda create -n p4dev python={3.6} psi4-dev -c psi4[/label/dev]
> # intermediate cmd: conda env create -n p4dev -f env_p4dev.yaml && conda activate p4dev
cmake -S. -GNinja -C/psi/gits/hrw-tag/cache_ppa_plain.cmake -Bobjdir_ppa_plain && cmake --build objdir_ppa_plain
#>>> conda activate p4dev
#>>> psi4-path-advisor
## execute or adapt `cmake` commands above; DepsCache handles python & addons;
##   DepsMKLCache handles math; further psi4-path-advisor options handle compilers.
#>>> cd objdir && make -j`getconf _NPROCESSORS_ONLN`
    #env_create_and_activate_cmd = f"{subdircmd}conda env create -n {args.name} -f {envspec_fn} --solver libmamba  && conda activate {args.name}"
bash for eval. conda available.
#>>> make install""")

parser.add_argument("-v", action="count", default=0,
    help="""Use for more printing (-vv).
Do not use with bash command substitution: eval $(psi4-path-advisor args)""")

subparsers = parser.add_subparsers(dest="subparser_name",
    help="???")

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
    default="mamba-as-possible",
    choices=solver_choices,
    help="\n".join(solver_help))
parser_env.add_argument("--platform",
    choices=["linux-64", "osx-64", "osx-arm64", "win-64"],
    default=conda_platform_native,
    help=f"""Conda platform/subdir for env file, if not the computed native
(default: {conda_platform_native}). Apple Silicon users,
check this value! Argument rarely used.""")
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

args = parser.parse_args()

if args.subparser_name in ["conda", "env"]:
    conda_platform = args.platform
else:
    conda_platform = conda_platform_native


if args.v > 1:
    print("#######")
    parser.print_help()
    print("#######")
    parser_env.print_help()
    print("#######")
    parser_cmake.print_help()
    print("#######")
    print(f"{conda_platform=}  {conda_lapack_variant=}  {args=}")
    print("#######")

#            # pip
#          - pip
#          - pip:
#             - git+https://github.com/i-pi/i-pi.git@master-py3
#""")

with codedeps_yaml.open() as fp:
    ydict = yaml.load(fp, Loader=yaml.FullLoader)

##### ENV

if args.subparser_name in ["conda", "env"]:

    if args.solver == "mamba-as-possible":
        if conda_libmamba_available:
            solver = "conda-libmamba"
        elif mamba_available:
            solver = "mamba"
        elif conda_available:
            solver = "conda"
    else:
        solver = args.solver
    if solver == "conda-libmamba":
        solver = ["conda", "--solver libmamba"]
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


##### CACHE

elif args.subparser_name in ["cmake", "cache"]:

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
            pkgstr = [pkg if (pkg in conda_list_pkgver) else strike(pkg) for pkg in plain_package_set]
            primary_banner = f"""# <<<  {pkgstr[0]:<{compute_width(pkgstr[0], 27)}} {", ".join(pkgstr[1:])}"""
            if not use["required"]:
                absent_merge.append(primary_banner)
                continue

        if primary == "libblas" and args.lapack == "openblas" and conda_platform == "linux-64":
            if conda_openblas_variant != "openmp":
                all_found = False
                pkgstr.append(strike("openblas=*=openmp*"))
                primary_banner = f"""# <<<  {pkgstr[0]:<{compute_width(pkgstr[0], 27)}} {", ".join(pkgstr[1:])}"""

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
            else:
                dcmake_vars = conda["cmake"]

            for k, v in dcmake_vars.items():

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
                    text.append(f'set({k:<30} {v} CACHE {ctyp} "")')

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
    big_args["-C"] = cmakecache
    big_args_sorted = [f"{k}{big_args[k]}" for k in sorted(big_args.keys(), reverse=True)]

    pretext = f"""# {cmakecache_lbl}
# {'-' * len(cmakecache_lbl)}
#
# This file is autogenerated by psi4-path-advisor from conda env {conda_prefix_short}.
#   It sets some likely values to initialize the CMake cache for dependencies
#   to build your Psi4 source.
#
# Dependency packages are shown as \"present-package\" or \"{strike('missing-package')}\".
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
#parser.add_argument('--plugin-compile', action='store_true', help="""\
#(Duplicate from `psi4`) Generates a CMake command for building a plugin against this Psi4 installation.
#>>> cd <plugin_directory>
#>>> `psi4 --plugin-compile`
#>>> make
#>>> psi4""")
#
#args = parser.parse_args()
#
#
#if psi4alongside:
#    from subprocess import call
#
#    if args.psiapi_path:
#        call([psi4, '--psiapi-path'])
#        sys.exit(0)
#
#    if args.plugin_compile:
#        call([psi4, '--plugin-compile'])
#        sys.exit(0)
#
#else:
#    if args.plugin_compile:
#        print("""Install "psi4" via `conda install psi4 -c psi4[/label/dev]`, then reissue command.""")
#
#
#
##advice = {
##    'cmake':  '/opt/anaconda1anaconda2anaconda3/bin/cmake \\',
##    'here':   '    -S. \\',
##    'deps':   '    -C/opt/anaconda1anaconda2anaconda3/share/cmake/psi4/psi4DepsCache.cmake \\',
##    'nooptl': '    -C/opt/anaconda1anaconda2anaconda3/share/cmake/psi4/psi4DepsDisableCache.cmake \\',
##    'Lintel': '    -C/opt/anaconda1anaconda2anaconda3/share/cmake/psi4/psi4DepsIntelCache.cmake \\',
##    'Lgnu':   '    -C/opt/anaconda1anaconda2anaconda3/share/cmake/psi4/psi4DepsGNUCache.cmake \\',
##    'Mclang': '    -C/opt/anaconda1anaconda2anaconda3/share/cmake/psi4/psi4DepsAppleClangCache.cmake \\',
##    'Mgnu':   '    -C/opt/anaconda1anaconda2anaconda3/share/cmake/psi4/psi4DepsGNUCache.cmake \\',
##    'objdir': '    -Bobjdir',
##}
#
#recc = ['/opt/anaconda1anaconda2anaconda3/bin/cmake',
#        '-S.',
#        '-C/opt/anaconda1anaconda2anaconda3/share/cmake/psi4/psi4DepsCache.cmake',
#        '-Bobjdir']
#
#if args.disable_addons:
#    recc.insert(-1, '-C/opt/anaconda1anaconda2anaconda3/share/cmake/psi4/psi4DepsDisableCache.cmake')
#
#if args.mkl:
#    recc.insert(-1, '-C/opt/anaconda1anaconda2anaconda3/share/cmake/psi4/psi4DepsMKLCache.cmake')
#
#if sys.platform == 'darwin':
#    if args.clang:
#        recc.insert(-1, '-C/opt/anaconda1anaconda2anaconda3/share/cmake/psi4/psi4DepsClangCache.cmake')
#    #    recc.insert(0, 'CONDA_BUILD_SYSROOT=~/SDKs/MacOSX10.9.sdk')
#    #    recc.insert(-1, '-C/opt/anaconda1anaconda2anaconda3/share/cmake/psi4/psi4DepsAppleClangCache.cmake')
#    #if args.gcc:
#    #    recc.insert(-1, '-C/opt/anaconda1anaconda2anaconda3/share/cmake/psi4/psi4DepsGNUCache.cmake')
#
#srecc = """    """.join(recc)
#print(srecc)

def native_platform():
    # fn is good but unused. since this can't distinguish the metal chip for osx, we instead use
    #   conda info platform directly.
    import sys
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


# For Silicon owners whose miniconda is set to Intel:
#   conda config --system --set subdir osx-arm64

