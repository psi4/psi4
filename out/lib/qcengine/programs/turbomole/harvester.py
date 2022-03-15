import itertools as it
import re
from decimal import Decimal

import numpy as np

from ..util import PreservingDict


def parse_decimal(regex, text, method="search"):
    with_float = re.compile(regex + "([\d\-\.]+)")
    matches = getattr(with_float, method)(text)

    if method == "search":
        # groups = matches.groups()
        # if len(groups) == 1:
        # groups = ("", groups[0])
        matches = [matches.groups()]
    return [(method, Decimal(energy)) for method, energy in matches]


def parse_reference_energy(stdout: str):
    """Parse stdout and return the energy of the reference wavefunction."""
    energy_dict = PreservingDict()

    # Total energy from dscf or ridft
    total_energy_re = re.compile("total energy\s+=\s+([\d\-\.]+)")
    mobj = total_energy_re.search(stdout)
    total_energy = Decimal(mobj[1])

    # Check for DFT, default to HF
    energy_key = "HF TOTAL ENERGY"
    dft_mobj = re.search("density functional", stdout)
    if dft_mobj:
        energy_key = "DFT TOTAL ENERGY"
    energy_dict[energy_key] = total_energy

    # Take into account energies from ricc2 runs. They will be different
    # from the HF energy.
    current_energy = total_energy
    energy_dict["CURRENT ENERGY"] = current_energy

    return energy_dict


def parse_ricc2(stdout: str):
    ricc2_dict = PreservingDict()

    # As CC2 starts from a MP2 guess that is also reported there may be
    # multiple matches for the following regex. Thats why we capture all
    # matches with 'findall'.
    matches = parse_decimal("Final (.+?) energy\s+:\s+", stdout, "findall")
    if len(matches) == 0:
        matches = parse_decimal("E(MP2)\s+:\s+", stdout, "search")

    ricc2_dict["CURRENT ENERGY"] = matches[-1][1]

    return ricc2_dict


def parse_gradient(gradient):
    grad_re = re.compile(
        "\$grad.+"
        "cycle =\s+(?P<cycle>\d+)\s+"
        "(?P<energy_type>.+?) energy =\s+(?P<energy>[\d\.\-]+)\s+"
        "\|dE/dxyz\| =\s+(?P<grad_norm>[\d\.]+)"
        "(?P<coords_gradients>.+)\$end",
        re.DOTALL,
    )
    mobj = grad_re.match(gradient)

    # Commented out the variables that aren't returned so LGTM doesn't
    # complain.
    # energy_type = mobj.group("energy_type")
    # grad_norm = Decimal(mobj.group("grad_norm"))
    # energy = Decimal(mobj.group("energy"))
    coords_grad = mobj.group("coords_gradients")
    # cycle = int(mobj.group("cycle"))

    *_, grad = re.split("[a-z]{1,3}", coords_grad.strip())
    grad = np.array(grad.strip().replace("D", "E").split(), dtype=float)

    return grad


def parse_nprhessian(nprhessian):
    lines = [l.strip() for l in nprhessian.strip().split("\n")]
    assert lines[0] == "$nprhessian"
    assert lines[-1] == "$end"

    hess_lines = [line.split()[2:] for line in lines[1:-1]]
    atom_num = int(lines[-2].split()[0])
    hessian = np.array(list(it.chain(*hess_lines)), dtype=float)
    hessian = hessian.reshape(atom_num, atom_num)

    return hessian


def parse_hessian(hessian, size):
    lines = [l.strip() for l in hessian.strip().split("\n")]
    assert lines[0] == "$hessian"
    assert lines[-1] == "$end"

    hess_lines = list()
    for line in lines[1:-1]:
        line = line.strip()
        if line == "$hessian (projected)":
            break
        hess_lines.append(line.split()[2:])
    else:
        raise Exception("'$hessian (projected)' line was not encountered!")

    hessian = np.array(list(it.chain(*hess_lines)), dtype=float)
    hessian = hessian.reshape(size, size)

    return hessian


def harvest(input_model, stdout, **outfiles):
    qcvars = PreservingDict()

    ref_energy_dict = parse_reference_energy(stdout)
    qcvars.update(ref_energy_dict)

    if "R I C C 2" in stdout:
        ricc2_dict = parse_ricc2(stdout)
        qcvars.update(ricc2_dict)

    gradient = None
    if "gradient" in outfiles:
        gradient = parse_gradient(outfiles["gradient"])

    hessian = None
    if "nprhessian" in outfiles:
        hessian = parse_nprhessian(outfiles["nprhessian"])
    if "hessian" in outfiles:
        size = input_model.geometry.size
        hessian = parse_hessian(outfiles["hessian"], size)

    return qcvars, gradient, hessian
