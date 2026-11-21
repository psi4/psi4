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
Terms of the interaction-induced SAPT properties.

Ported from propSAPT v0.5.0 (https://github.com/btyrcha/propSAPT). A property term is
the derivative of a SAPT energy correction with respect to the strength of a field
coupled to a one-electron operator, so the monomers relax in the field: each operator
component needs the CPSCF amplitudes of both monomers.
"""

import time
from functools import partial
from typing import Dict, List, Tuple

import numpy as np

from psi4 import core

from ...p4util import solvers
from ...p4util.exceptions import ConvergenceError
from .sapt_prop_exch_disp import exch_disp2_s2_property
from .sapt_prop_exch_ind import exch_ind2_s2_property, u_ra_amplitudes, u_sb_amplitudes
from .sapt_prop_exch_ind_sinf import exch_ind2_sinf_property

__all__ = ["cpscf_batch", "property_terms"]

# the contractions keep propSAPT's index strings
contract = partial(np.einsum, optimize=True)


def cpscf_batch(wfn: core.Wavefunction, eps_occ: np.ndarray, eps_vir: np.ndarray, perturbations: List[np.ndarray],
                conv: float, maxiter: int) -> List[np.ndarray]:
    """
    Solves the CPSCF (CPHF or CPKS) equations of one monomer for several perturbations.

    All perturbations are solved in one conjugate-gradient run: each iteration builds
    the orbital-Hessian products of the unconverged vectors with one J/K call, and each
    vector converges on its own.

    Parameters
    ----------
    wfn
        The monomer wavefunction, whose J/K object builds the Hessian products.
    eps_occ, eps_vir
        Occupied and virtual orbital energies of the monomer.
    perturbations
        Perturbation matrices, (occ, vir) each.
    conv
        Convergence threshold of the relative residual norm.
    maxiter
        Maximum number of iterations.

    Returns
    -------
    list of numpy.ndarray
        The (vir, occ) amplitudes, in the order of `perturbations`.
    """

    nocc, nvir = len(eps_occ), len(eps_vir)
    amplitudes = [np.zeros((nvir, nocc)) for _ in perturbations]

    # a perturbation without an (occ, vir) block has zero amplitudes
    solve = [num for num, pert in enumerate(perturbations) if np.any(pert)]
    if not solve:
        return amplitudes

    rhs = [core.Matrix.from_array(np.ascontiguousarray(perturbations[num])) for num in solve]
    rhs_norms = [x.sum_of_squares() for x in rhs]

    eps_ov = core.Matrix.from_array(eps_occ.reshape(-1, 1) - eps_vir)

    def hessian_vec(x_vec, act_mask):
        active = [x for x, on in zip(x_vec, act_mask) if on]
        products = iter(wfn.cphf_Hx(active) if active else [])
        return [next(products) if on else False for on in act_mask]

    def apply_precon(x_vec, act_mask):
        preconditioned = []
        for x, on in zip(x_vec, act_mask):
            if on:
                p = x.clone()
                p.apply_denominator(eps_ov)
                preconditioned.append(p)
            else:
                preconditioned.append(False)
        return preconditioned

    x_vec, r_vec = solvers.cg_solver(rhs, hessian_vec, apply_precon, printlvl=2, maxiter=maxiter, rcond=conv)

    resid = [(r.sum_of_squares() / norm)**0.5 for r, norm in zip(r_vec, rhs_norms)]
    if max(resid) >= conv:
        raise ConvergenceError(
            "SAPT property CPSCF", maxiter,
            f"Largest relative residual {max(resid):.3e} is above CPHF_R_CONVERGENCE {conv:.1e}; "
            "consider raising MAXITER of the SAPT module.")

    for num, x in zip(solve, x_vec):
        amplitudes[num] = np.array(x).reshape(nocc, nvir).T

    return amplitudes


def _cpscf_round(monomer: str, wfn: core.Wavefunction, eps_occ: np.ndarray, eps_vir: np.ndarray,
                 perturbations: List[np.ndarray], conv: float, maxiter: int, what: str) -> List[np.ndarray]:
    """:py:func:`cpscf_batch` with a header and the timing in the output."""

    core.print_out(f"\n  ==> CPSCF of monomer {monomer}: {len(perturbations)} {what} <==\n\n")
    tstart = time.time()
    amplitudes = cpscf_batch(wfn, eps_occ, eps_vir, perturbations, conv, maxiter)
    core.print_out(f"\n    ...finished in {time.time() - tstart:.2f} seconds.\n")
    return amplitudes


def _ind2_terms(cache: dict, xt_A_ra: np.ndarray, xt_B_sb: np.ndarray, prop_A_aa: np.ndarray, prop_A_rr: np.ndarray,
                prop_B_bb: np.ndarray, prop_B_ss: np.ndarray) -> Tuple[float, float]:
    """
    The second-order induction terms of monomers A and B with the uncoupled induction
    amplitudes, `x2_ind_A` and `x2_ind_B`: the derivative of IND20,U.
    """

    Qaa, Qar, Qbb, Qbs, Qrr, Qss = (cache[key] for key in ("Qaa", "Qar", "Qbb", "Qbs", "Qrr", "Qss"))
    omegaA_bb, omegaA_ss, tA_sb = (cache[key] for key in ("omegaA_bb", "omegaA_ss", "tA_sb"))
    omegaB_aa, omegaB_rr, tB_ra = (cache[key] for key in ("omegaB_aa", "omegaB_rr", "tB_ra"))

    # yapf: disable
    x2_ind_A = (
        + 16 * contract("ra,Qar,Qbs,sb", xt_A_ra, Qar, Qbs, tA_sb)
        - 4 * contract("ra,ac,rc", tB_ra, omegaB_aa, xt_A_ra)
        + 4 * contract("ra,cr,ca", tB_ra, omegaB_rr, xt_A_ra)
        - 2 * contract("ac,rc,ra", prop_A_aa, tB_ra, tB_ra)
        + 2 * contract("cr,ca,ra", prop_A_rr, tB_ra, tB_ra)
        # X_A -> X_A + [W_A, C_A - C_A^\dagger]_1
        + 4 * contract("ra,cd,de,Qer,Qac", xt_A_ra, tB_ra.T, tB_ra, Qar, Qaa)
        + 8 * contract("ra,cd,ec,Qde,Qar", xt_A_ra, tB_ra.T, tB_ra, Qrr, Qar)
        - 4 * contract("ra,cd,ec,Qdr,Qae", xt_A_ra, tB_ra.T, tB_ra, Qrr, Qar)
        - 8 * contract("ra,cd,de,Qar,Qec", xt_A_ra, tB_ra.T, tB_ra, Qar, Qaa)
    )
    x2_ind_B = (
        + 16 * contract("sb,Qar,Qbs,ra", xt_B_sb, Qar, Qbs, tB_ra)
        - 4 * contract("sb,bc,sc", tA_sb, omegaA_bb, xt_B_sb)
        + 4 * contract("sb,cs,cb", tA_sb, omegaA_ss, xt_B_sb)
        - 2 * contract("bc,sc,sb", prop_B_bb, tA_sb, tA_sb)
        + 2 * contract("cs,cb,sb", prop_B_ss, tA_sb, tA_sb)
        # X_B -> X_B + [W_B, C_B - C_B^\dagger]_1
        + 4 * contract("sb,cd,de,Qes,Qbc", xt_B_sb, tA_sb.T, tA_sb, Qbs, Qbb)
        + 8 * contract("sb,cd,ec,Qde,Qbs", xt_B_sb, tA_sb.T, tA_sb, Qss, Qbs)
        - 4 * contract("sb,cd,ec,Qds,Qbe", xt_B_sb, tA_sb.T, tA_sb, Qss, Qbs)
        - 8 * contract("sb,cd,de,Qbs,Qec", xt_B_sb, tA_sb.T, tA_sb, Qbs, Qbb)
    )
    # yapf: enable

    return float(x2_ind_A), float(x2_ind_B)


def _ind2_r_terms(cache: dict, xt_A_ra: np.ndarray, xt_B_sb: np.ndarray, prop_A_aa: np.ndarray, prop_A_rr: np.ndarray,
                  prop_B_bb: np.ndarray, prop_B_ss: np.ndarray, t_ra: np.ndarray,
                  t_sb: np.ndarray) -> Tuple[float, float]:
    """
    The second-order induction terms of monomers A and B with the coupled (CPSCF)
    induction amplitudes `t_ra` and `t_sb`, `x2_ind,r_A` and `x2_ind,r_B`: the
    derivative of IND20,R.
    """

    Qaa, Qar, Qbb, Qbs, Qrr, Qss = (cache[key] for key in ("Qaa", "Qar", "Qbb", "Qbs", "Qrr", "Qss"))
    omegaA_bb, omegaA_ss = (cache[key] for key in ("omegaA_bb", "omegaA_ss"))
    omegaB_aa, omegaB_rr = (cache[key] for key in ("omegaB_aa", "omegaB_rr"))

    # yapf: disable
    x2_ind_r_A = (
        + 16 * contract("ra,Qar,Qbs,sb", xt_A_ra, Qar, Qbs, t_sb)
        - 4 * contract("ra,ac,rc", t_ra, omegaB_aa, xt_A_ra)
        + 4 * contract("ra,cr,ca", t_ra, omegaB_rr, xt_A_ra)
        - 2 * contract("ac,rc,ra", prop_A_aa, t_ra, t_ra)
        + 2 * contract("cr,ca,ra", prop_A_rr, t_ra, t_ra)
        # 2 Re <C^{(1)}_A | W_A C^{[1]}_A C^{(1)}_A > + Re <C^{(1)}_A C^{(1)}_A | W_A C^{[1]}_A >
        + 4 * contract("ra,cr,de,Qad,Qec", xt_A_ra, t_ra.T, t_ra, Qar, Qaa)
        + 4 * contract("ra,cd,de,Qer,Qac", xt_A_ra, t_ra.T, t_ra, Qar, Qaa)
        + 16 * contract("ra,ac,de,Qcr,Qed", xt_A_ra, t_ra.T, t_ra, Qrr, Qar)
        + 8 * contract("ra,cd,ec,Qde,Qar", xt_A_ra, t_ra.T, t_ra, Qrr, Qar)
        - 4 * contract("ra,ac,de,Qcd,Qer", xt_A_ra, t_ra.T, t_ra, Qrr, Qar)
        - 4 * contract("ra,cd,ec,Qdr,Qae", xt_A_ra, t_ra.T, t_ra, Qrr, Qar)
        - 16 * contract("ra,cr,de,Qed,Qac", xt_A_ra, t_ra.T, t_ra, Qar, Qaa)
        - 8 * contract("ra,cd,de,Qar,Qec", xt_A_ra, t_ra.T, t_ra, Qar, Qaa)
        + 4 * contract("ra,cr,de,Qae,Qcd", xt_A_ra, t_ra.T, t_ra, Qaa, Qar)
        - 4 * contract("ra,ac,de,Qec,Qrd", xt_A_ra, t_ra.T, t_ra, Qar, Qrr)
    )

    x2_ind_r_B = (
        + 16 * contract("sb,Qar,Qbs,ra", xt_B_sb, Qar, Qbs, t_ra)
        - 4 * contract("sb,bc,sc", t_sb, omegaA_bb, xt_B_sb)
        + 4 * contract("sb,cs,cb", t_sb, omegaA_ss, xt_B_sb)
        - 2 * contract("bc,sc,sb", prop_B_bb, t_sb, t_sb)
        + 2 * contract("cs,cb,sb", prop_B_ss, t_sb, t_sb)
        # 2 Re <C^{(1)}_B | W_B C^{[1]}_B C^{(1)}_B > + Re <C^{(1)}_B C^{(1)}_B | W_B C^{[1]}_B >
        + 4 * contract("sb,cs,de,Qbd,Qec", xt_B_sb, t_sb.T, t_sb, Qbs, Qbb)
        + 4 * contract("sb,cd,de,Qes,Qbc", xt_B_sb, t_sb.T, t_sb, Qbs, Qbb)
        + 16 * contract("sb,bc,de,Qcs,Qed", xt_B_sb, t_sb.T, t_sb, Qss, Qbs)
        + 8 * contract("sb,cd,ec,Qde,Qbs", xt_B_sb, t_sb.T, t_sb, Qss, Qbs)
        - 4 * contract("sb,bc,de,Qcd,Qes", xt_B_sb, t_sb.T, t_sb, Qss, Qbs)
        - 4 * contract("sb,cd,ec,Qds,Qbe", xt_B_sb, t_sb.T, t_sb, Qss, Qbs)
        - 16 * contract("sb,cs,de,Qed,Qbc", xt_B_sb, t_sb.T, t_sb, Qbs, Qbb)
        - 8 * contract("sb,cd,de,Qbs,Qec", xt_B_sb, t_sb.T, t_sb, Qbs, Qbb)
        + 4 * contract("sb,cs,de,Qbe,Qcd", xt_B_sb, t_sb.T, t_sb, Qbb, Qbs)
        - 4 * contract("sb,bc,de,Qec,Qsd", xt_B_sb, t_sb.T, t_sb, Qbs, Qss)
    )
    # yapf: enable

    return float(x2_ind_r_A), float(x2_ind_r_B)


def _disp2_terms(cache: dict, xt_A_ra: np.ndarray, xt_B_sb: np.ndarray, prop_A_aa: np.ndarray, prop_A_rr: np.ndarray,
                 prop_B_bb: np.ndarray, prop_B_ss: np.ndarray) -> Tuple[float, float]:
    """
    The second-order dispersion terms of monomers A and B with the uncoupled dispersion
    amplitudes, `x2_disp_A` and `x2_disp_B`: the derivative of DISP20,U.
    """

    Qaa, Qar, Qbb, Qbs, Qrr, Qss = (cache[key] for key in ("Qaa", "Qar", "Qbb", "Qbs", "Qrr", "Qss"))
    t_rsab = cache["t_rsab"]

    # yapf: disable
    x2_disp_A = (
        # 2 Re <R(X)|V R(V)>
        8 * contract("rsab,Qcr,Qbs,ca", t_rsab, Qrr, Qbs, xt_A_ra)
        - 8 * contract("rsab,Qac,Qbs,rc", t_rsab, Qaa, Qbs, xt_A_ra)
        # <R(V)|X R(V)>
        - 4 * contract("ac,rscb,rsab", prop_A_aa, t_rsab, t_rsab)
        + 4 * contract("cr,csab,rsab", prop_A_rr, t_rsab, t_rsab)
        # X_A -> X_A + [W_A, C_A - C_A^\dagger]_1
        + 8 * contract("ra,dscb,dseb,Qer,Qac", xt_A_ra, t_rsab, t_rsab, Qar, Qaa)
        + 16 * contract("ra,dscb,escb,Qar,Qde", xt_A_ra, t_rsab, t_rsab, Qar, Qrr)
        - 8 * contract("ra,dscb,escb,Qdr,Qae", xt_A_ra, t_rsab, t_rsab, Qrr, Qar)
        - 16 * contract("ra,dscb,dseb,Qar,Qec", xt_A_ra, t_rsab, t_rsab, Qar, Qaa)
    )
    x2_disp_B = (
        # 2 Re <R(X)|V R(V)>
        8 * contract("rsab,Qar,Qcs,cb", t_rsab, Qar, Qss, xt_B_sb)
        - 8 * contract("rsab,Qar,Qbc,sc", t_rsab, Qar, Qbb, xt_B_sb)
        # <R(V)|X R(V)>
        - 4 * contract("bc,rsac,rsab", prop_B_bb, t_rsab, t_rsab)
        + 4 * contract("cs,rcab,rsab", prop_B_ss, t_rsab, t_rsab)
        # X_B -> X_B + [W_B, C_B - C_B^\dagger]_1
        + 8 * contract("sb,rdac,rdae,Qes,Qbc", xt_B_sb, t_rsab, t_rsab, Qbs, Qbb)
        + 16 * contract("sb,rdac,reac,Qbs,Qde", xt_B_sb, t_rsab, t_rsab, Qbs, Qss)
        - 8 * contract("sb,rdac,reac,Qds,Qbe", xt_B_sb, t_rsab, t_rsab, Qss, Qbs)
        - 16 * contract("sb,rdac,rdae,Qbs,Qec", xt_B_sb, t_rsab, t_rsab, Qbs, Qbb)
    )
    # yapf: enable

    return float(x2_disp_A), float(x2_disp_B)


def property_terms(cache: dict,
                   operators: List[np.ndarray],
                   wfn_A: core.Wavefunction,
                   wfn_B: core.Wavefunction,
                   conv: float,
                   maxiter: int,
                   coupled: bool = True,
                   sinf: bool = False,
                   second_order: bool = True) -> Dict[str, np.ndarray]:
    """
    The monomer, first-order and second-order SAPT property terms of each operator.

    Parameters
    ----------
    cache
        Intermediates from :py:func:`~psi4.driver.procrouting.sapt.sapt_prop_cache.build_prop_cache`,
        with the fitted 3-index integrals if `second_order`.
    operators
        AO matrices of the property components.
    wfn_A, wfn_B
        The monomer wavefunctions.
    conv, maxiter
        Convergence threshold and maximum number of iterations of the CPSCF equations.
    coupled
        Whether to compute the terms with the coupled (CPSCF) induction amplitudes too.
        The ones with the uncoupled amplitudes are always computed.
    sinf
        Whether to compute the S^inf exchange-induction terms.
    second_order
        Whether to compute the second-order terms at all.

    Returns
    -------
    dict
        One value per operator for each term, keyed as in propSAPT:

        - "x0_A", "x0_B": electronic part of the monomer property;
        - "x1_pol,r": first-order polarization, the derivative of ELST10,R;
        - "x1_exch,r_S2": first-order exchange in the S^2 approximation, the derivative
          of EXCH10(S^2);
        - "x1_exch,r": first-order exchange, the derivative of EXCH10;
        - "x2_ind": second-order induction with the uncoupled amplitudes, the
          derivative of IND20,U;
        - "x2_exch-ind_S2", "x2_exch-ind": second-order exchange-induction with the
          uncoupled amplitudes, S^2 and S^inf;
        - "x2_ind,r", "x2_exch-ind,r_S2", "x2_exch-ind,r": the same with the coupled
          amplitudes, if `coupled`;
        - "x2_disp", "x2_exch-disp_S2": second-order dispersion and S^2
          exchange-dispersion, with the uncoupled amplitudes, the derivatives of
          DISP20,U and EXCH-DISP20,U(S^2);
        - "ind_S2", "disp_S2": the induction and the dispersion classes, each term plus
          its S^2 exchange counterpart, as Psi4's SAPT IND and SAPT DISP energies, the
          induction with the coupled amplitudes if `coupled`;
        - "x_induced_S2": the sum of the first-order terms and the two classes, as Psi4's
          SAPT TOTAL energy;

        and the parts of the first-order, the induction and the dispersion terms from
        monomers A and B, keyed with "_A" and "_B" appended. The S^inf
        exchange-induction terms are there only if `sinf`.
    """

    C_a, C_r = cache["C_a"], cache["C_r"]
    C_b, C_s = cache["C_b"], cache["C_s"]
    eps_a, eps_r = cache["eps_a"], cache["eps_r"]
    eps_b, eps_s = cache["eps_b"], cache["eps_s"]

    prop_A_aa = [C_a.T @ op @ C_a for op in operators]
    prop_A_ar = [C_a.T @ op @ C_r for op in operators]
    prop_B_bb = [C_b.T @ op @ C_b for op in operators]
    prop_B_bs = [C_b.T @ op @ C_s for op in operators]

    # relaxed amplitudes of the monomers in the field of each operator, solved together
    # with the coupled induction amplitudes, the response of each monomer to the
    # electrostatic potential of the other
    induction = second_order and coupled
    rhs_A = prop_A_ar + ([cache["omegaB_ar"]] if induction else [])
    rhs_B = prop_B_bs + ([cache["omegaA_bs"]] if induction else [])
    what = "operators and the induction amplitudes" if induction else "operators"
    xt_A = _cpscf_round("A", wfn_A, eps_a, eps_r, rhs_A, conv, maxiter, what)
    xt_B = _cpscf_round("B", wfn_B, eps_b, eps_s, rhs_B, conv, maxiter, what)
    if induction:
        t_ra, t_sb = xt_A.pop(), xt_B.pop()

    terms = {}

    # SCF monomers
    terms["x0_A"] = np.array([2 * contract("aa", p) for p in prop_A_aa])
    terms["x0_B"] = np.array([2 * contract("bb", p) for p in prop_B_bb])

    # First-order polarization
    terms["x1_pol,r_A"] = np.array([4 * contract("ra,ar", xt, cache["omegaB_ar"]) for xt in xt_A])
    terms["x1_pol,r_B"] = np.array([4 * contract("sb,bs", xt, cache["omegaA_bs"]) for xt in xt_B])

    # First-order exchange, S^2 and S^inf
    for tag in ("_S2", ""):
        terms[f"x1_exch,r{tag}_A"] = np.array([
            2 * contract("ra,ar", xt, cache[f"omega_exchB_ar{tag}"]) +
            2 * contract("ra,ra", cache[f"omega_exchB_ra{tag}"], xt) for xt in xt_A
        ])
        terms[f"x1_exch,r{tag}_B"] = np.array([
            2 * contract("sb,bs", xt, cache[f"omega_exchA_bs{tag}"]) +
            2 * contract("sb,sb", cache[f"omega_exchA_sb{tag}"], xt) for xt in xt_B
        ])

    for term in ("x1_pol,r", "x1_exch,r_S2", "x1_exch,r"):
        terms[term] = terms[f"{term}_A"] + terms[f"{term}_B"]

    if not second_order:
        return terms

    prop_A_rr = [C_r.T @ op @ C_r for op in operators]
    prop_B_ss = [C_s.T @ op @ C_s for op in operators]
    components = range(len(operators))

    def per_component(term, *args):
        """`term` of every operator component, its arguments `args` then those of the component."""
        return [
            term(cache, xt_A[n], xt_B[n], *(arg[n] if isinstance(arg, list) else arg for arg in args))
            for n in components
        ]

    # Second-order induction, with the uncoupled and the coupled amplitudes
    ind_tags = [""]
    parts = per_component(_ind2_terms, prop_A_aa, prop_A_rr, prop_B_bb, prop_B_ss)
    terms["x2_ind_A"], terms["x2_ind_B"] = (np.array(part) for part in zip(*parts))
    if coupled:
        ind_tags.append(",r")
        parts = per_component(_ind2_r_terms, prop_A_aa, prop_A_rr, prop_B_bb, prop_B_ss, t_ra, t_sb)
        terms["x2_ind,r_A"], terms["x2_ind,r_B"] = (np.array(part) for part in zip(*parts))
    for tag in ind_tags:
        terms[f"x2_ind{tag}"] = terms[f"x2_ind{tag}_A"] + terms[f"x2_ind{tag}_B"]

    # Second-order exchange-induction, S^2 and S^inf, with the coupled amplitudes and
    # their derivatives from the CPSCF equations, then with the uncoupled ones
    modes = [(True, t_ra, t_sb, ",r")] if coupled else []
    modes.append((False, cache["tB_ra"], cache["tA_sb"], ""))
    for mode_coupled, amps_ra, amps_sb, tag in modes:
        u_ra = per_component(u_ra_amplitudes, prop_A_aa, prop_A_rr, amps_ra, mode_coupled)
        u_sb = per_component(u_sb_amplitudes, prop_B_bb, prop_B_ss, amps_sb, mode_coupled)
        if mode_coupled:
            what = "derivatives of the induction amplitudes"
            u_ra = _cpscf_round("A", wfn_A, eps_a, eps_r, [u.T for u in u_ra], conv, maxiter, what)
            u_sb = _cpscf_round("B", wfn_B, eps_b, eps_s, [u.T for u in u_sb], conv, maxiter, what)

        amplitudes = (u_ra, u_sb, amps_ra, amps_sb)
        terms[f"x2_exch-ind{tag}_S2"] = np.array(per_component(exch_ind2_s2_property, *amplitudes))

        # the S^inf expression holds the induction term as well
        if sinf:
            ind_exch_ind = np.array(per_component(exch_ind2_sinf_property, *amplitudes))
            terms[f"x2_exch-ind{tag}"] = ind_exch_ind - terms[f"x2_ind{tag}"]

    # Second-order dispersion and S^2 exchange-dispersion, uncoupled. The S^inf
    # exchange-dispersion is left out: propSAPT's is not the derivative of its energy.
    parts = per_component(_disp2_terms, prop_A_aa, prop_A_rr, prop_B_bb, prop_B_ss)
    terms["x2_disp_A"], terms["x2_disp_B"] = (np.array(part) for part in zip(*parts))
    terms["x2_disp"] = terms["x2_disp_A"] + terms["x2_disp_B"]
    terms["x2_exch-disp_S2"] = np.array(
        per_component(exch_disp2_s2_property, prop_A_aa, prop_A_rr, prop_B_bb, prop_B_ss))

    # Totals with the S^2 exchange terms at second order, propSAPT's x_induced_S2; the
    # induction with the coupled amplitudes unless only the uncoupled ones were made
    ind = ",r" if coupled else ""
    terms["ind_S2"] = terms[f"x2_ind{ind}"] + terms[f"x2_exch-ind{ind}_S2"]
    terms["disp_S2"] = terms["x2_disp"] + terms["x2_exch-disp_S2"]
    terms["x_induced_S2"] = terms["x1_pol,r"] + terms["x1_exch,r"] + terms["ind_S2"] + terms["disp_S2"]

    return terms
