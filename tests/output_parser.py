#!/usr/bin/env python3

import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Union


PathLike = Union[str, Path]

ITERATION_RE = re.compile(
    r"^\s*@(?P<label>[A-Za-z0-9+\-]+)\s+iter\s+(?P<iteration>\d+):(?P<rest>.*)$"
)
SCF_ACCELERATOR_RE = re.compile(
    r"\b(?:A|E)?DIIS\b(?:/\b(?:A|E)?DIIS\b)*|\bSOSCF\b",
    re.IGNORECASE,
)
OPT_CONVERGENCE_CHECK_RE = re.compile(r"^\s*==>\s+Convergence Check\s+<==\s*$")
OPT_TABLE_HEADER_RE = re.compile(
    r"^\s+Step\s+Total Energy\s+Delta E\s+MAX Force\s+RMS Force\s+MAX Disp\s+RMS Disp\s*~?\s*$",
    re.IGNORECASE,
)
OPT_STEP_RE = re.compile(
    r"^\s+\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+e[-+]\d+\b.*~\s*$",
    re.IGNORECASE,
)


def extract_scf_iterations(output_file: PathLike) -> Dict:
    """Extract SCF and optimization markers from a Psi4 output file.

    This ignores SAD guess lines such as::

        @RHF iter SAD:

    and counts numeric iteration lines such as::

        @RHF iter   1:
        @DF-RHF iter   1:
        @UHF iter   8:

    Numeric iteration ``0`` is counted because Psi4 emits it as a real SCF
    iteration row for some calculation paths.

    Per-iteration SCF accelerator labels are counted from SCF iteration rows.
    Mixed DIIS-family labels such as ``ADIIS/DIIS`` are split and counted under
    both labels. SOSCF is counted when it appears on a numeric SCF iteration
    row.

    Explicit geometry optimization step rows from OptKing ``Convergence Check``
    sections are counted as ``optimization.cycles``. A bare trailing ``~`` is
    not sufficient because other output sections use it for decoration.
    """
    output_path = Path(output_file)
    scf_total_iterations = 0
    scf_iterations_by_label = Counter()
    scf_accelerators = Counter()
    optimization = Counter()
    in_opt_convergence_check = False
    in_opt_step_table = False

    with output_path.open("r", encoding="utf-8", errors="replace") as infile:
        for line in infile:
            if OPT_CONVERGENCE_CHECK_RE.match(line):
                in_opt_convergence_check = True
                in_opt_step_table = False
            elif in_opt_convergence_check and OPT_TABLE_HEADER_RE.match(line):
                in_opt_step_table = True
            elif in_opt_step_table and OPT_STEP_RE.match(line):
                optimization["cycles"] += 1
                in_opt_convergence_check = False
                in_opt_step_table = False

            match = ITERATION_RE.match(line)
            if not match:
                continue

            label = match.group("label")

            scf_total_iterations += 1
            scf_iterations_by_label[label] += 1
            for accelerator_match in SCF_ACCELERATOR_RE.finditer(match.group("rest")):
                for accelerator_label in accelerator_match.group(0).upper().split("/"):
                    scf_accelerators[accelerator_label] += 1

    return {
        "test_name": output_path.parent.name,
        "source_file": output_path.name,
        "scf": {
            "total_iterations": scf_total_iterations,
            "iterations_by_label": dict(sorted(scf_iterations_by_label.items())),
            "accelerators": dict(sorted(scf_accelerators.items())),
        },
        "optimization": {
            "cycles": optimization["cycles"],
        },
    }


def write_scf_iterations_json(output_file: PathLike, json_file: PathLike = None) -> Path:
    """Extract SCF iteration counts and write them to a JSON file."""
    output_path = Path(output_file)
    json_path = Path(json_file) if json_file is not None else output_path.with_name("scf_iterations.json")

    data = extract_scf_iterations(output_path)
    with json_path.open("w", encoding="utf-8") as outfile:
        json.dump(data, outfile, indent=2)
        outfile.write("\n")

    return json_path


def main(argv=None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) not in (1, 2):
        print(
            "Usage: python tests/output_parser.py <output.dat> [output.json]",
            file=sys.stderr,
        )
        return 1

    output_file = Path(args[0])
    json_file = Path(args[1]) if len(args) == 2 else output_file.with_name("scf_iterations.json")

    data = extract_scf_iterations(output_file)
    with json_file.open("w", encoding="utf-8") as outfile:
        json.dump(data, outfile, indent=2)
        outfile.write("\n")

    print(f"Wrote {json_file}")
    print(f"Total numeric SCF iterations: {data['scf']['total_iterations']}")
    print("SCF iterations by label:")
    for label, count in data["scf"]["iterations_by_label"].items():
        print(f"  {label}: {count}")
    if data["scf"]["accelerators"]:
        print("SCF accelerators:")
        for label, count in data["scf"]["accelerators"].items():
            print(f"  {label}: {count}")
    if data["optimization"]["cycles"]:
        print("Optimization:")
        for label, count in data["optimization"].items():
            print(f"  {label}: {count}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
