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
from .sapt_prop_terms import first_order_properties

__all__ = ["run_sapt_property"]

# The CPSCF convergence of SAPT properties unless CPHF_R_CONVERGENCE is set: the
# properties are first derivatives, so the CPSCF error enters them linearly.
_PROPERTY_CPHF_R_CONVERGENCE = 1.0e-10


def _term_labels(name: str) -> dict:
    """QCVariable labels of the property terms, keyed as in propSAPT.

    The first-order polarization term is the derivative of the electrostatic energy,
    but as a property it is a change of the density, hence POL. SAPT0 terms carry the
    intramonomer correlation order, 0, as the second index; SAPT(DFT) terms do not.
    """

    order = "10" if name == "sapt0" else "1"
    return {
        "x1_pol,r": f"POL{order},R",
        "x1_exch,r_S2": f"EXCH{order}(S^2)",
        "x1_exch,r": f"EXCH{order}",
    }


def _print_header(name: str, properties: list):
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
    core.print_out("\n")


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
    ``- B``. The SAPT interaction energy is computed too.
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

    _print_header(name, properties)

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

    # Property terms
    core.timer_on("SAPT properties: Properties")
    tstart = time.time()
    components, quantities = build_operators(properties, wfn_A.basisset(), sapt_dimer)
    cache = build_prop_cache(wfn_A, wfn_B, sapt_jk)
    terms = first_order_properties(cache, [c.matrix for c in components], wfn_A, wfn_B, cphf_conv, maxiter)

    # monomer values with their nuclear parts
    charges_A = np.array([monomerA.Z(atom) for atom in range(monomerA.natom())])
    charges_B = np.array([monomerB.Z(atom) for atom in range(monomerB.natom())])
    rows = [
        ("MONOMER A", terms["x0_A"] + np.array([c.nuclear(charges_A) for c in components])),
        ("MONOMER B", terms["x0_B"] + np.array([c.nuclear(charges_B) for c in components])),
    ]
    for key, label in _term_labels(name).items():
        rows.append((label, terms[key]))
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
