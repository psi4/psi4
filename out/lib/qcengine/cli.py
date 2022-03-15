"""
Provides a CLI for QCEngine
"""

import argparse
import json
import os.path
import sys
from typing import Any, Dict

from . import get_program  # run and run-procedure; info
from . import (
    __version__,
    compute,
    compute_procedure,
    get_procedure,
    list_all_procedures,
    list_all_programs,
    list_available_procedures,
    list_available_programs,
)
from .config import global_repr  # info

__all__ = ["main"]

info_choices = frozenset(["version", "programs", "procedures", "config", "all"])


def parse_args():
    parser = argparse.ArgumentParser(description="A CLI for the QCEngine.")
    parser.add_argument("--version", action="version", version=f"{__version__}")

    subparsers = parser.add_subparsers(dest="command")

    info = subparsers.add_parser("info", help="Print information about QCEngine setup, version, and environment.")
    info.add_argument(
        "category", nargs="*", default="all", choices=info_choices, help="The information categories to show."
    )

    run = subparsers.add_parser("run", help="Run a program on a given task. Output is printed as a JSON blob.")
    run.add_argument("program", type=str, help="The program to run.")
    run.add_argument(
        "data",
        type=str,
        help="Data describing the task to run. "
        "One of: (i) A JSON blob, "
        "(ii) A file name, "
        "(iii) '-', indicating data will be read from STDIN.",
    )

    run_procedure = subparsers.add_parser(
        "run-procedure", help="Run a procedure on a given task. " "Output is printed as a JSON blob."
    )
    run_procedure.add_argument("procedure", type=str, help="The procedure to run.")
    run_procedure.add_argument(
        "data",
        type=str,
        help="Data describing the task to run. "
        "One of: (i) A JSON blob, "
        "(ii) A file name, "
        "(iii) '-', indicating data will be read from STDIN.",
    )

    args = vars(parser.parse_args())
    if args["command"] is None:
        parser.print_help(sys.stderr)
        exit(1)

    return args


def info_cli(args):
    def info_version():
        import qcelemental

        print(">>> Version information")
        print(f"QCEngine:    {__version__}")
        print(f"QCElemental: {qcelemental.__version__}")
        print()

    def info_programs():  # lgtm: [py/similar-function]
        print(">>> Program information")
        all_progs = list_all_programs()
        avail_progs = list_available_programs()
        print("Available programs:")
        for prog_name in sorted(avail_progs):
            program = get_program(prog_name)
            version = program.get_version()
            for loc, ver in program.version_cache.items():
                if ver == version:
                    which = loc
                    break
            else:
                which = "???"
            if version is None:
                version = "???"
            print(f"{prog_name + ':':12} v{version:20} {which}")

        print()
        print("Other supported programs:")
        print(" ".join(sorted(all_progs - avail_progs)))
        print()
        print(
            """If you think available programs are missing, query for details: `python -c "import qcengine as qcng; qcng.get_program('<program>')"`"""
        )
        print()

    def info_procedures():  # lgtm: [py/similar-function]
        print(">>> Procedure information")
        all_procs = list_all_procedures()
        avail_procs = list_available_procedures()
        print("Available procedures:")
        for proc_name in sorted(avail_procs):
            version = get_procedure(proc_name).get_version()
            if version is None:
                version = "???"
            print(f"{proc_name} v{version}")

        print()
        print("Other supported procedures:")
        print(" ".join(sorted(all_procs - avail_procs)))
        print()

    # default=["all"] does is not allowed by argparse
    if not isinstance(args["category"], list):
        args["category"] = [args["category"]]
    cat = set(args["category"])

    if "version" in cat or "all" in cat:
        info_version()
    if "programs" in cat or "all" in cat:
        info_programs()
    if "procedures" in cat or "all" in cat:
        info_procedures()
    if "config" in cat or "all" in cat:
        print(">>> Configuration information")
        print()
        print(global_repr())


def data_arg_helper(data_arg: str) -> Dict[str, Any]:
    """
    Converts the data argument of run and run-procedure commands to a dict for compute or compute_procedure

    Parameters
    ----------
    data_arg: str
        Either a data blob or file name or '-' for STDIN

    Returns
    -------
    Dict[str, Any]
        An input for compute or compute_procedure.
    """
    if data_arg == "-":
        return json.load(sys.stdin)
    elif os.path.isfile(data_arg):
        return json.load(open(data_arg))
    else:
        return json.loads(data_arg)


def main(args=None):
    # Grab CLI args if not present
    if args is None:
        args = parse_args()
    command = args.pop("command")
    if command == "info":
        info_cli(args)
    elif command == "run":
        ret = compute(data_arg_helper(args["data"]), args["program"])
        print(ret.json())
    elif command == "run-procedure":
        ret = compute_procedure(data_arg_helper(args["data"]), args["procedure"])
        print(ret.json())
