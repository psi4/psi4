import logging
from decimal import Decimal as Dm
from typing import Any, Dict, List

import numpy as np
from qcelemental.models import AtomicResultProperties

from .util import PreservingDict

logger = logging.getLogger(__name__)


def _difference(args):
    minuend, subtrahend = args
    return minuend - subtrahend


def _product(args):
    multiplicand, multiplier = args
    return multiplicand * multiplier


def _spin_component_scaling_wsing(args):
    os_scale, ss_scale, tot_corl, ss_corl, s_corl = args
    return os_scale * (tot_corl - ss_corl - s_corl) + ss_scale * ss_corl + s_corl


def _dispersion_weighting_wsing(args):
    omega, hb_mtd, dd_mtd, s_corl = args
    return omega * hb_mtd + (1.0 - omega) * dd_mtd + s_corl


def omega(args):
    alpha, beta, ratio = args
    return 0.5 * (1.0 + np.tanh(alpha + beta * ratio))


def _linear(args):
    return sum(c * v for c, v in zip(args))


def _solve_in_turn(args: List, coeff: List) -> List[Dict[str, Any]]:
    """Expands a description of a linear equation of variables `args` and weights `coeff` into all rearrangements of the equation."""

    assert len(args) == len(coeff)
    pv0 = []

    for itgt in range(len(args)):
        non_target_args = args[:]
        non_target_args.pop(itgt)

        non_target_coeff = coeff[:]
        non_target_coeff.pop(itgt)
        solve_by = -1 // coeff[itgt]
        non_target_coeff = [solve_by * c for c in non_target_coeff]

        pv0.append(
            {
                "form": args[itgt],
                "func": lambda vv, cc=non_target_coeff: sum(c * v for c, v in zip(vv, cc)),
                "args": non_target_args,
            }
        )

    return pv0


def qcvar_identities() -> List[Dict[str, Any]]:
    """Define QCVariable identity equations (e.g., method total = method correlation + HF)."""

    pv0 = []

    # MP2
    pv0.extend(_solve_in_turn(args=["MP2 TOTAL ENERGY", "HF TOTAL ENERGY", "MP2 CORRELATION ENERGY"], coeff=[-1, 1, 1]))
    pv0.extend(
        _solve_in_turn(
            args=["MP2 DOUBLES ENERGY", "MP2 SAME-SPIN CORRELATION ENERGY", "MP2 OPPOSITE-SPIN CORRELATION ENERGY"],
            coeff=[-1, 1, 1],
        )
    )
    pv0.extend(
        _solve_in_turn(args=["MP2 CORRELATION ENERGY", "MP2 DOUBLES ENERGY", "MP2 SINGLES ENERGY"], coeff=[-1, 1, 1])
    )

    pv0.extend(
        _solve_in_turn(
            args=["CUSTOM SCS-MP2 TOTAL ENERGY", "CUSTOM SCS-MP2 CORRELATION ENERGY", "HF TOTAL ENERGY"],
            coeff=[-1, 1, 1],
        )
    )

    # SCS-MP2
    pv0.append(
        {
            "form": "SCS-MP2 CORRELATION ENERGY",
            "func": _spin_component_scaling_wsing,
            "args": [
                Dm(6) / Dm(5),
                Dm(1) / Dm(3),
                "MP2 CORRELATION ENERGY",
                "MP2 SAME-SPIN CORRELATION ENERGY",
                "MP2 SINGLES ENERGY",
            ],
        }
    )  # yapf: disable
    pv0.append({"form": "SCS-MP2 TOTAL ENERGY", "func": sum, "args": ["HF TOTAL ENERGY", "SCS-MP2 CORRELATION ENERGY"]})

    # SCS(N)-MP2
    pv0.append(
        {
            "form": "SCS(N)-MP2 CORRELATION ENERGY",
            "func": _spin_component_scaling_wsing,
            "args": [
                Dm(0),
                Dm(1.76),
                "MP2 CORRELATION ENERGY",
                "MP2 SAME-SPIN CORRELATION ENERGY",
                "MP2 SINGLES ENERGY",
            ],
        }
    )  # yapf: disable
    pv0.append(
        {"form": "SCS(N)-MP2 TOTAL ENERGY", "func": sum, "args": ["HF TOTAL ENERGY", "SCS(N)-MP2 CORRELATION ENERGY"]}
    )

    # SCS-MP2-VDW
    # * https://doi.org/10.1080/00268970802641242
    pv0.append(
        {
            "form": "SCS-MP2-VDW CORRELATION ENERGY",
            "func": _spin_component_scaling_wsing,
            "args": [
                Dm(1.28),
                Dm(0.50),
                "MP2 CORRELATION ENERGY",
                "MP2 SAME-SPIN CORRELATION ENERGY",
                "MP2 SINGLES ENERGY",
            ],
        }
    )  # yapf: disable
    pv0.append(
        {"form": "SCS-MP2-VDW TOTAL ENERGY", "func": sum, "args": ["HF TOTAL ENERGY", "SCS-MP2-VDW CORRELATION ENERGY"]}
    )

    # DW-MP2
    # only defined at the (IE) reaction level (like SAPT)
    #    dwmp2['DW-MP2 OMEGA'][mt] = \
    #        rxnm_contract_expand(df.xs('HF TOTAL ENERGY', level='qcvar').xs(mtl[0], level='meta')) / \
    #        rxnm_contract_expand(df.xs('MP2 TOTAL ENERGY', level='qcvar').xs(mtl[1], level='meta'))
    #    df_omega = omega([0.15276, 1.89952, df_omega])

    # pv0.append({
    #    'form': 'DW-MP2 CORRELATION ENERGY',
    #    'func': _dispersion_weighting_wsing,
    #    'args': ['DW_MP2 OMEGA', 'MP2 CORRELATION ENERGY', 'SCS-MP2 CORRELATION ENERGY', 'MP2 SINGLES ENERGY']})
    # pv0.append({
    #    'form': 'DW-MP2 TOTAL ENERGY',
    #    'func': sum,
    #    'args': ['HF TOTAL ENERGY', 'DW-MP2 CORRELATION ENERGY']})

    # MP3
    pv0.extend(_solve_in_turn(args=["MP3 TOTAL ENERGY", "HF TOTAL ENERGY", "MP3 CORRELATION ENERGY"], coeff=[-1, 1, 1]))
    pv0.extend(
        _solve_in_turn(
            args=["MP3 DOUBLES ENERGY", "MP3 SAME-SPIN CORRELATION ENERGY", "MP3 OPPOSITE-SPIN CORRELATION ENERGY"],
            coeff=[-1, 1, 1],
        )
    )
    pv0.extend(
        _solve_in_turn(args=["MP3 CORRELATION ENERGY", "MP3 DOUBLES ENERGY", "MP3 SINGLES ENERGY"], coeff=[-1, 1, 1])
    )

    # MPN
    for mpn in range(4, 20):
        pv0.extend(
            _solve_in_turn(
                args=["MP{} TOTAL ENERGY".format(mpn), "HF TOTAL ENERGY", "MP{} CORRELATION ENERGY".format(mpn)],
                coeff=[-1, 1, 1],
            )
        )

    # CCSD
    pv0.extend(
        _solve_in_turn(args=["CCSD TOTAL ENERGY", "HF TOTAL ENERGY", "CCSD CORRELATION ENERGY"], coeff=[-1, 1, 1])
    )
    pv0.extend(
        _solve_in_turn(
            args=["CCSD DOUBLES ENERGY", "CCSD SAME-SPIN CORRELATION ENERGY", "CCSD OPPOSITE-SPIN CORRELATION ENERGY"],
            coeff=[-1, 1, 1],
        )
    )
    pv0.extend(
        _solve_in_turn(args=["CCSD CORRELATION ENERGY", "CCSD DOUBLES ENERGY", "CCSD SINGLES ENERGY"], coeff=[-1, 1, 1])
    )

    # CCSD(T)
    pv0.extend(
        _solve_in_turn(
            args=["CCSD(T) CORRELATION ENERGY", "CCSD CORRELATION ENERGY", "(T) CORRECTION ENERGY"], coeff=[-1, 1, 1]
        )
    )
    pv0.extend(
        _solve_in_turn(args=["CCSD(T) TOTAL ENERGY", "HF TOTAL ENERGY", "CCSD(T) CORRELATION ENERGY"], coeff=[-1, 1, 1])
    )
    pv0.extend(
        _solve_in_turn(
            args=["CCSD(T) CORRELATION ENERGY", "CCSD CORRELATION ENERGY", "(T) CORRECTION ENERGY"], coeff=[-1, 1, 1]
        )
    )  # duplicate of first so that all minimal combinations covered

    # CCSD[T]
    pv0.extend(
        _solve_in_turn(args=["CCSD[T] TOTAL ENERGY", "HF TOTAL ENERGY", "CCSD[T] CORRELATION ENERGY"], coeff=[-1, 1, 1])
    )
    pv0.extend(
        _solve_in_turn(
            args=["CCSD[T] CORRELATION ENERGY", "CCSD CORRELATION ENERGY", "[T] CORRECTION ENERGY"], coeff=[-1, 1, 1]
        )
    )

    # FCI
    pv0.extend(_solve_in_turn(args=["FCI TOTAL ENERGY", "HF TOTAL ENERGY", "FCI CORRELATION ENERGY"], coeff=[-1, 1, 1]))

    # Generics
    pv0.extend(_solve_in_turn(args=["CC TOTAL ENERGY", "HF TOTAL ENERGY", "CC CORRELATION ENERGY"], coeff=[-1, 1, 1]))
    pv0.extend(_solve_in_turn(args=["CI TOTAL ENERGY", "HF TOTAL ENERGY", "CI CORRELATION ENERGY"], coeff=[-1, 1, 1]))

    # DFT

    #   fctl
    #       TODO want B97 here?
    for fctl in [
        "B3LYP",
        "B3LYP5",
        "WPBE",
        "PBE",
        "CAM-B3LYP",
        "B97",
        "WB97X",
        "SVWN",
        "PW91",
        "BLYP",
        "PW86PBE",
        "FT97",
        "BOP",
        "MPWPW",
        "SOGGA11",
        "BP86",
        "B86BPBE",
        "PW6B95",
        "PBE0",
        "SOGGA11-X",
        "MN15",
    ]:
        pv0.extend(
            _solve_in_turn(
                args=["{} TOTAL ENERGY".format(fctl), "{} FUNCTIONAL TOTAL ENERGY".format(fctl)], coeff=[-1, 1]
            )
        )

    #   fctl + D
    for dash in ["-D2", "-D3", "-D3(BJ)", "-D3M", "-D3M(BJ)"]:
        for fctl in ["B3LYP", "B3LYP5", "PBE", "B97", "BLYP", "BP86", "PBE0", "WPBE"]:
            pv0.extend(
                _solve_in_turn(
                    args=[
                        "{}{} TOTAL ENERGY".format(fctl, dash),
                        "{} FUNCTIONAL TOTAL ENERGY".format("B97-D" if fctl == "B97" else fctl),
                        "{}{} DISPERSION CORRECTION ENERGY".format(fctl, dash),
                    ],
                    coeff=[-1, 1, 1],
                )
            )

    #   fctl + DH
    for fctl in ["B2PLYP", "DSD-PBEP86", "PBE0-2", "B2GPPLYP", "PTPSS", "PWPB95", "DSD-BLYP", "PBE0-DH"]:
        pv0.extend(
            _solve_in_turn(
                args=[
                    "{} TOTAL ENERGY".format(fctl),
                    "{} FUNCTIONAL TOTAL ENERGY".format(fctl),
                    "{} DOUBLE-HYBRID CORRECTION ENERGY".format(fctl),
                ],
                coeff=[-1, 1, 1],
            )
        )

    #   fctl + D + DH
    #   no qcvar for fctl + dh, which would be the more restrictive def
    for dash in ["-D2", "-D3", "-D3(BJ)", "-D3M", "-D3M(BJ)"]:
        for fctl in ["B2PLYP"]:
            pv0.extend(
                _solve_in_turn(
                    args=[
                        "{}{} TOTAL ENERGY".format(fctl, dash),
                        "{} TOTAL ENERGY".format(fctl),
                        "{}{} DISPERSION CORRECTION ENERGY".format(fctl, dash),
                    ],
                    coeff=[-1, 1, 1],
                )
            )

    #   misc.
    pv0.extend(
        _solve_in_turn(
            args=[
                "DLDF+D09 TOTAL ENERGY",
                "DLDF+D09 FUNCTIONAL TOTAL ENERGY",
                "DLDF-DAS2009 DISPERSION CORRECTION ENERGY",
            ],
            coeff=[-1, 1, 1],
        )
    )

    pv0.extend(
        _solve_in_turn(
            args=["WB97X-D TOTAL ENERGY", "WB97X FUNCTIONAL TOTAL ENERGY", "WB97-CHG DISPERSION CORRECTION ENERGY"],
            coeff=[-1, 1, 1],
        )
    )

    # misc.
    pv0.extend(
        _solve_in_turn(
            args=["CURRENT ENERGY", "CURRENT REFERENCE ENERGY", "CURRENT CORRELATION ENERGY"], coeff=[-1, 1, 1]
        )
    )

    return pv0


def build_out(rawvars: Dict[str, Any], verbose: int = 1) -> None:
    """Apply standard QC identities to QCVariables `rawvars` to build more (e.g., correlation from total and HF energies).

    Dictionary `qcvar_identities` has keys with names of QCVariables to be created and values with dictionary of two
    keys: `args`, the QCVariables that contribute to the key and `func`, a functional (or lambda) to combine them.
    This function builds that key QCVariables if all the contributors are available in `rawvars`. Updates internally
    so multiple passes not needed.

    Parameters
    ----------
    verbose : int, optional
        Controls print level. Per-var printing with >=2.

    Returns
    -------
    None
        But input dictionary `rawvars` is updated.

    """
    for action in qcvar_identities():
        pvar = action["form"]
        buildline = """building {} {}""".format(pvar, "." * (50 - len(pvar)))

        data_rich_args = []

        for pv in action["args"]:
            if isinstance(pv, str):
                if pv in rawvars:
                    data_rich_args.append(rawvars[pv])
                else:
                    if verbose >= 2:
                        logger.debug("""{}EMPTY, missing {}""".format(buildline, pv))
                    break
            else:
                data_rich_args.append(pv)
        else:
            result = action["func"](data_rich_args)
            # rawvars[pvar] = result
            # with data coming from file --> variable, looks more precise than it is. hack
            rawvars.__setitem__(pvar, result, 6)
            if verbose >= 1:
                logger.debug("""{}SUCCESS""".format(buildline))

            if pvar == "CURRENT CORRELATION ENERGY" and abs(float(rawvars[pvar])) < 1.0e-16:
                rawvars.pop(pvar)


qcvars_to_atomicproperties = {
    # calcinfo
    "N BASIS FUNCTIONS": "calcinfo_nbasis",
    "N MOLECULAR ORBITALS": "calcinfo_nmo",
    "N ATOMS": "calcinfo_natom",
    "N ALPHA ELECTRONS": "calcinfo_nalpha",
    "N BETA ELECTRONS": "calcinfo_nbeta",
    "NUCLEAR REPULSION ENERGY": "nuclear_repulsion_energy",
    # pre keywords
    "ONE-ELECTRON ENERGY": "scf_one_electron_energy",
    "TWO-ELECTRON ENERGY": "scf_two_electron_energy",
    "CURRENT ENERGY": "return_energy",
    "CURRENT GRADIENT": "return_gradient",
    "CURRENT HESSIAN": "return_hessian",
    # HF keywords
    "HF DIPOLE": "scf_dipole_moment",
    "HF QUADRUPOLE": "scf_quadrupole_moment",
    "HF TOTAL ENERGY": "scf_total_energy",
    "HF ITERATIONS": "scf_iterations",
    "HF TOTAL GRADIENT": "scf_total_gradient",
    "HF TOTAL HESSIAN": "scf_total_hessian",
    # MP2 keywords
    "MP2 SAME-SPIN CORRELATION ENERGY": "mp2_same_spin_correlation_energy",
    "MP2 OPPOSITE-SPIN CORRELATION ENERGY": "mp2_opposite_spin_correlation_energy",
    "MP2 SINGLES ENERGY": "mp2_singles_energy",
    "MP2 DOUBLES ENERGY": "mp2_doubles_energy",
    "MP2 CORRELATION ENERGY": "mp2_correlation_energy",
    "MP2 TOTAL ENERGY": "mp2_total_energy",
    "MP2 DIPOLE": "mp2_dipole_moment",
    # CCSD keywords
    "CCSD SAME-SPIN CORRELATION ENERGY": "ccsd_same_spin_correlation_energy",
    "CCSD OPPOSITE-SPIN CORRELATION ENERGY": "ccsd_opposite_spin_correlation_energy",
    "CCSD SINGLES ENERGY": "ccsd_singles_energy",
    "CCSD DOUBLES ENERGY": "ccsd_doubles_energy",
    "CCSD CORRELATION ENERGY": "ccsd_correlation_energy",
    "CCSD TOTAL ENERGY": "ccsd_total_energy",
    "CCSD DIPOLE": "ccsd_dipole_moment",
    "CCSD ITERATIONS": "ccsd_iterations",
    # CCSD(T) keywords
    "CCSD(T) CORRELATION ENERGY": "ccsd_prt_pr_correlation_energy",
    "CCSD(T) TOTAL ENERGY": "ccsd_prt_pr_total_energy",
    "CCSD(T) DIPOLE": "ccsd_prt_pr_dipole_moment",
    # informal keywords
    # scf_vv10_energy
    "DFT XC ENERGY": "scf_xc_energy",
    "DFT DISPERSION ENERGY": "scf_dispersion_correction_energy",
    # multi-definition
    "SCF TOTAL ENERGY": "scf_total_energy",
}


def build_atomicproperties(qcvars: "PreservingDict") -> AtomicResultProperties:
    """For results extracted from QC output in QCDB terminology, translate to QCSchema terminology.

    Parameters
    ----------
    qcvars : PreservingDict
        Dictionary of calculation information in QCDB QCVariable terminology.

    Returns
    -------
    atprop : AtomicResultProperties
        Object of calculation information in QCSchema AtomicResultProperties terminology.

    """
    atprop = {}
    for pv, dpv in qcvars.items():
        if pv in qcvars_to_atomicproperties:
            atprop[qcvars_to_atomicproperties[pv]] = dpv

    return AtomicResultProperties(**atprop)
