#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2026 The Psi4 Developers.
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
"""
Interaction-induced one-electron properties from SAPT0 and SAPT(DFT).

The property terms are the derivatives of the SAPT energy corrections with respect to
the strength of a field coupled to a one-electron operator, after propSAPT v0.5.0
(https://github.com/btyrcha/propSAPT) and B. Tyrcha, T. Gupta, K. Patkowski and
P. S. Zuchowski, J. Chem. Theory Comput. 21, 4562 (2025).
"""

import time

import numpy as np

from psi4 import core

from ... import p4util
from ...p4util.exceptions import ValidationError
from .. import proc_util
from ..proc import scf_helper
from .sapt_proc import compute_GRAC_shift, sapt_dft
from .sapt_prop_cache import build_prop_cache
from .sapt_prop_operators import build_operators, validate_properties
from .sapt_prop_terms import property_terms

__all__ = ["run_sapt_property"]

# The CPSCF convergence of SAPT properties unless CPHF_R_CONVERGENCE is set: the
# properties are first derivatives, so the CPSCF error enters them linearly.
_PROPERTY_CPHF_R_CONVERGENCE = 1.0e-10


def _term_labels(name: str, second_order: bool, coupled: bool, sinf: bool) -> list:
    """QCVariable labels of the property terms that are computed.

    The first-order polarization term is the derivative of the electrostatic energy,
    but as a property it is a change of the density, hence POL. SAPT0 terms carry the
    intramonomer correlation order, 0, as the second index; SAPT(DFT) terms do not.
    Both exchange-induction variants are suffixed, (S^2) and (S^INF); the dispersion is
    uncoupled, ",U". The totals IND, DISP and TOTAL are named after Psi4's SAPT energy
    classes, whose derivatives they are.

    Returns
    -------
    list
        (propSAPT key, label, whether the term has parts of monomers A and B) of each term.
    """

    first, second = ("10", "20") if name == "sapt0" else ("1", "2")
    labels = [
        ("x1_pol,r", f"POL{first},R", True),
        ("x1_exch,r_S2", f"EXCH{first}(S^2)", True),
        ("x1_exch,r", f"EXCH{first}", True),
    ]
    if not second_order:
        return labels

    # propSAPT's keys of the coupled (",r") and the uncoupled ("") induction amplitudes
    modes = ([(",r", ",R")] if coupled else []) + [("", ",U")]
    for tag, kind in modes:
        labels.append((f"x2_ind{tag}", f"IND{second}{kind}", True))
    for tag, kind in modes:
        labels.append((f"x2_exch-ind{tag}_S2", f"EXCH-IND{second}{kind}(S^2)", False))
        if sinf:
            labels.append((f"x2_exch-ind{tag}", f"EXCH-IND{second}{kind}(S^INF)", False))
    labels.append(("x2_disp", f"DISP{second},U", True))
    labels.append(("x2_exch-disp_S2", f"EXCH-DISP{second},U(S^2)", False))
    labels += [("ind_S2", "IND", False), ("disp_S2", "DISP", False), ("x_induced_S2", "TOTAL", False)]
    return labels


def _print_header(name: str, properties: list, second_order: bool, coupled: bool, sinf: bool):
    core.print_out("\n")
    core.print_out("         ---------------------------------------------------------\n")
    core.print_out("         " + "SAPT Interaction-Induced Properties".center(58) + "\n")
    core.print_out("\n")
    core.print_out("         " + "ported from propSAPT by B. Tyrcha".center(58) + "\n")
    core.print_out("         ---------------------------------------------------------\n")
    core.print_out("\n")
    core.print_out("  ==> Algorithm <==\n\n")
    core.print_out(f"   Method                  {name.upper():>12s}\n")
    core.print_out(f"   Properties              {', '.join(properties)}\n")
    core.print_out(f"   JK Algorithm            {core.get_global_option('SCF_TYPE'):>12s}\n")
    core.print_out(f"   CPHF R Convergence      {core.get_option('SAPT', 'CPHF_R_CONVERGENCE'):12.3e}\n")
    core.print_out(f"   Second Order            {str(second_order):>12s}\n")
    if second_order:
        core.print_out(f"   Coupled Induction       {str(coupled):>12s}\n")
        core.print_out(f"   S^inf Exch-Ind          {str(sinf):>12s}\n")
    core.print_out("\n")
    if not second_order:
        core.print_out(
            "   The second-order properties are computed for SAPT0 only: with Kohn-Sham monomers, their\n"
            "   expressions are not yet the field derivatives of the induction and dispersion energies.\n\n")


def _check_memory(wfn_A: core.Wavefunction, wfn_B: core.Wavefunction, df_basis: core.BasisSet):
    """Prints the memory the second-order terms need, and raises if it is above the memory
    Psi4 is given: the fitted 3-index integrals of both monomers, the dispersion amplitudes
    and their denominators, and one more occ(A) occ(B) vir(A) vir(B) intermediate at a time
    (the derivative of the dispersion amplitudes, or that of the S^inf exchange-induction
    terms)."""

    naux = df_basis.nbf()
    o_A, o_B = wfn_A.nalpha(), wfn_B.nalpha()
    v_A, v_B = wfn_A.nmo() - o_A, wfn_B.nmo() - o_B
    df_doubles = naux * (o_A * o_A + o_A * v_A + v_A * v_A + o_B * o_B + o_B * v_B + v_B * v_B)
    oovv_doubles = o_A * o_B * v_A * v_B
    needed = 8 * (df_doubles + 3 * oovv_doubles)
    available = core.get_memory()

    gib = 1024.0**3
    core.print_out("  ==> Memory <==\n\n")
    core.print_out(f"   Fitted 3-index integrals [GiB]    {8 * df_doubles / gib:10.3f}\n")
    core.print_out(f"   Dispersion amplitudes [GiB]       {16 * oovv_doubles / gib:10.3f}\n")
    core.print_out(f"   Largest intermediate [GiB]        {8 * oovv_doubles / gib:10.3f}\n")
    core.print_out(f"   Available [GiB]                   {available / gib:10.3f}\n\n")

    if needed > available:
        raise ValidationError(
            f"SAPT properties: the second-order terms need at least {needed / gib:.2f} GiB of memory, "
            f"but Psi4 is given {available / gib:.2f} GiB. Raise the memory with set_memory.")


def _print_property(quantity, rows):
    """Prints one property: a line per term, its components in columns."""

    ncol = 6
    core.print_out(f"\n   {quantity.name} [a.u.]\n")
    for start in range(0, len(quantity.labels), ncol):
        labels = quantity.labels[start:start + ncol]
        core.print_out("   " + " " * 22 + "".join(f"{label:>16s}" for label in labels) + "\n")
        for row_name, values in rows:
            core.print_out(f"   {row_name:22s}" +
                           "".join(f"{v:16.10f}" for v in values[quantity.components][start:start + ncol]) + "\n")


def run_sapt_property(name, **kwargs):
    """
    Interaction-induced one-electron properties of a dimer from SAPT0 or SAPT(DFT).

    The monomer wavefunctions are computed in the dimer-centered basis, as in
    :py:func:`~psi4.driver.procrouting.sapt.sapt_proc.run_sapt_dft`. For every property
    linear in the one-particle density, the monomer values and the first-order
    polarization and exchange terms are set as QCVariables
    ``SAPT <TERM> <PROPERTY>``, with the parts of monomers A and B suffixed ``- A`` and
    ``- B``. For SAPT0, so are the second-order terms: induction and exchange-induction
    with the uncoupled induction amplitudes, and with the coupled ones unless
    COUPLED_INDUCTION is false, S^2 and, with DO_IND_EXCH_SINF, S^inf; uncoupled
    dispersion and S^2 exchange-dispersion; and the totals IND, DISP and TOTAL, as the
    SAPT energy classes. The SAPT interaction energy is computed too.
    """

    name = name.lower()
    if name not in ("sapt0", "sapt(dft)"):
        raise ValidationError(f"SAPT properties: method {name} is not available, only SAPT0 and SAPT(DFT).")

    properties = validate_properties(kwargs.pop("properties", ["dipole", "quadrupole"]))

    optstash = p4util.OptionsState(
        ["SCF_TYPE"],
        ["SAVE_JK"],
        ["DFT_GRAC_SHIFT"],
        ["DF_INTS_IO"],
        ["SCF", "REFERENCE"],
        ["SCF", "DFT_GRAC_SHIFT"],
        ["SCF", "SAVE_JK"],
        ["SAPT", "CPHF_R_CONVERGENCE"],
    )

    core.tstart()
    # Alter default algorithm, as SAPT(DFT) does
    if not core.has_global_option_changed("SCF_TYPE"):
        core.set_global_option("SCF_TYPE", "DF")

    core.prepare_options_for_module("SAPT")

    if core.get_option("SCF", "REFERENCE") != "RHF":
        raise ValidationError("SAPT properties currently only support restricted references.")
    if core.get_global_option("FREEZE_CORE") not in ("FALSE", "0"):
        raise ValidationError("SAPT properties are all-electron, set FREEZE_CORE to false.")

    if not core.has_option_changed("SAPT", "CPHF_R_CONVERGENCE"):
        core.set_local_option("SAPT", "CPHF_R_CONVERGENCE", _PROPERTY_CPHF_R_CONVERGENCE)
    cphf_conv = core.get_option("SAPT", "CPHF_R_CONVERGENCE")
    maxiter = core.get_option("SAPT", "MAXITER")

    # The second-order terms are SAPT0 only: with Kohn-Sham monomers, propSAPT v0.5.0's
    # expressions take the field change of the Fock operator with the HF kernel, and miss
    # the finite-field derivatives of the induction and dispersion energies by about 1%.
    second_order = name == "sapt0"
    coupled = core.get_option("SAPT", "COUPLED_INDUCTION")
    sinf = core.get_option("SAPT", "DO_IND_EXCH_SINF")

    # Get the molecule of interest
    molecule = kwargs.pop("molecule", core.get_active_molecule())
    sapt_dimer, monomerA, monomerB = proc_util.prepare_sapt_molecule(molecule, "dimer")
    for label, monomer in (("A", monomerA), ("B", monomerB)):
        if monomer.multiplicity() != 1:
            raise ValidationError(f"SAPT properties: monomer {label} has multiplicity {monomer.multiplicity()}, "
                                  "only closed-shell monomers are supported.")

    # SAPT(DFT) functional and GRAC shifts, as in run_sapt_dft
    functional = core.get_option("SAPT", "SAPT_DFT_FUNCTIONAL") if name == "sapt(dft)" else "HF"
    do_dft = functional != "HF"
    mon_a_shift = mon_b_shift = 0.0
    if do_dft:
        mon_a_shift = core.get_option("SAPT", "SAPT_DFT_GRAC_SHIFT_A")
        mon_b_shift = core.get_option("SAPT", "SAPT_DFT_GRAC_SHIFT_B")
        grac_compute = core.get_option("SAPT", "SAPT_DFT_GRAC_COMPUTE")
        shift_a_set = core.has_option_changed("SAPT", "SAPT_DFT_GRAC_SHIFT_A")
        shift_b_set = core.has_option_changed("SAPT", "SAPT_DFT_GRAC_SHIFT_B")

        if (not (shift_a_set and shift_b_set)) and grac_compute.upper() == "NONE":
            raise ValidationError(
                'SAPT(DFT): User must set both "SAPT_DFT_GRAC_SHIFT_A" and "_B".  Or, to automatically compute the '
                'GRAC shift, set SAPT_DFT_GRAC_COMPUTE to "ITERATIVE" or "SINGLE".')
        if not shift_a_set:
            mon_a_shift = compute_GRAC_shift(sapt_dimer.extract_subsets(1), grac_compute, "Monomer A")
        if not shift_b_set:
            mon_b_shift = compute_GRAC_shift(sapt_dimer.extract_subsets(2), grac_compute, "Monomer B")

        core.set_variable("SAPT DFT GRAC SHIFT A", mon_a_shift)  # P::e SAPT
        core.set_variable("SAPT DFT GRAC SHIFT B", mon_b_shift)  # P::e SAPT

    _print_header(name, properties, second_order, coupled, sinf)

    # Monomer wavefunctions in the dimer-centered basis, sharing one J/K object
    core.IO.set_default_namespace("dimer")
    core.set_global_option("DF_INTS_IO", "SAVE")
    dimer_wfn = core.Wavefunction.build(sapt_dimer, core.get_global_option("BASIS"))

    if do_dft:
        core.set_local_option("SCF", "REFERENCE", "RKS")
    scf_method = functional if do_dft else "scf"

    core.timer_on("SAPT properties: Monomer A SCF")
    core.IO.change_file_namespace(97, "dimer", "monomerA")
    core.IO.set_default_namespace("monomerA")
    core.set_global_option("SAVE_JK", True)
    core.set_global_option("DFT_GRAC_SHIFT", mon_a_shift)
    wfn_A = scf_helper(scf_method, post_scf=False, molecule=monomerA, banner="SAPT properties: Monomer A", **kwargs)
    core.timer_off("SAPT properties: Monomer A SCF")

    core.timer_on("SAPT properties: Monomer B SCF")
    core.IO.change_file_namespace(97, "monomerA", "monomerB")
    core.IO.set_default_namespace("monomerB")
    core.set_global_option("SAVE_JK", True)
    core.set_global_option("DFT_GRAC_SHIFT", mon_b_shift)
    wfn_B = scf_helper(scf_method,
                       post_scf=False,
                       molecule=monomerB,
                       banner="SAPT properties: Monomer B",
                       jk=wfn_A.jk(),
                       **kwargs)
    core.timer_off("SAPT properties: Monomer B SCF")

    core.set_global_option("SAVE_JK", False)
    core.set_global_option("DFT_GRAC_SHIFT", 0.0)
    sapt_jk = wfn_B.jk()
    wfn_A.set_jk(sapt_jk)

    # Fitting basis of the second-order terms, as for the SAPT0 energy
    df_basis = None
    if second_order:
        df_basis = core.BasisSet.build(sapt_dimer, "DF_BASIS_SAPT", core.get_global_option("DF_BASIS_SAPT"), "RIFIT",
                                       core.get_global_option("BASIS"))
        dimer_wfn.set_basisset("DF_BASIS_SAPT", df_basis)
        _check_memory(wfn_A, wfn_B, df_basis)

    # Property terms
    core.timer_on("SAPT properties: Properties")
    tstart = time.time()
    components, quantities = build_operators(properties, wfn_A.basisset(), sapt_dimer)
    cache = build_prop_cache(wfn_A, wfn_B, sapt_jk, df_basis)
    terms = property_terms(cache, [c.matrix for c in components],
                           wfn_A,
                           wfn_B,
                           cphf_conv,
                           maxiter,
                           coupled=coupled,
                           sinf=sinf,
                           second_order=second_order)

    # monomer values with their nuclear parts
    charges_A = np.array([monomerA.Z(atom) for atom in range(monomerA.natom())])
    charges_B = np.array([monomerB.Z(atom) for atom in range(monomerB.natom())])
    rows = [
        ("MONOMER A", terms["x0_A"] + np.array([c.nuclear(charges_A) for c in components])),
        ("MONOMER B", terms["x0_B"] + np.array([c.nuclear(charges_B) for c in components])),
    ]
    for key, label, has_parts in _term_labels(name, second_order, coupled, sinf):
        rows.append((label, terms[key]))
        if has_parts:
            rows.append((f"{label} - A", terms[f"{key}_A"]))
            rows.append((f"{label} - B", terms[f"{key}_B"]))

    core.print_out("\n  ==> Interaction-Induced Properties <==\n")
    for quantity in quantities:
        _print_property(quantity, rows)

        for row_name, values in rows:
            # e.g. "SAPT POL10,R DIPOLE" and "SAPT POL10,R DIPOLE - A"
            if " - " in row_name:
                label, part = row_name.split(" - ")
                key = f"SAPT {label} {quantity.name} - {part}"
            else:
                key = f"SAPT {row_name} {quantity.name}"
            value = quantity.assemble(values[quantity.components])
            core.set_variable(key, value)
            dimer_wfn.set_variable(key, value)

    core.print_out(f"\n   ...properties finished in {time.time() - tstart:.2f} seconds.\n")
    core.timer_off("SAPT properties: Properties")

    # SAPT interaction energy from the same monomers. This sets CURRENT ENERGY, which
    # properties() returns, and frees the J/K object.
    data = {}
    sapt_dft(dimer_wfn, wfn_A, wfn_B, do_dft=do_dft, sapt_jk=sapt_jk, data=data, print_header=False, delta_hf=False)
    for key, value in data.items():
        core.set_variable(key, value)
        dimer_wfn.set_variable(key, value)

    core.tstop()
    optstash.restore()

    return dimer_wfn
