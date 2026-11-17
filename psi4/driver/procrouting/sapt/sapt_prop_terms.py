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
from typing import Dict, List

import numpy as np

from psi4 import core

from ...p4util import solvers
from ...p4util.exceptions import ConvergenceError

__all__ = ["cpscf_batch", "first_order_properties"]

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


def first_order_properties(cache: dict, operators: List[np.ndarray], wfn_A: core.Wavefunction,
                           wfn_B: core.Wavefunction, conv: float, maxiter: int) -> Dict[str, np.ndarray]:
    """
    The monomer and first-order SAPT property terms of each operator.

    Parameters
    ----------
    cache
        Intermediates from :py:func:`~psi4.driver.procrouting.sapt.sapt_prop_cache.build_prop_cache`.
    operators
        AO matrices of the property components.
    wfn_A, wfn_B
        The monomer wavefunctions.
    conv, maxiter
        Convergence threshold and maximum number of iterations of the CPSCF equations.

    Returns
    -------
    dict
        One value per operator for each term, keyed as in propSAPT:

        - "x0_A", "x0_B": electronic part of the monomer property;
        - "x1_pol,r": first-order polarization, the derivative of ELST10,R;
        - "x1_exch,r_S2": first-order exchange in the S^2 approximation, the derivative
          of EXCH10(S^2);
        - "x1_exch,r": first-order exchange, the derivative of EXCH10;

        and the parts of each first-order term from monomers A and B, keyed with "_A"
        and "_B" appended.
    """

    C_a, C_r = cache["C_a"], cache["C_r"]
    C_b, C_s = cache["C_b"], cache["C_s"]

    prop_A_aa = [C_a.T @ op @ C_a for op in operators]
    prop_A_ar = [C_a.T @ op @ C_r for op in operators]
    prop_B_bb = [C_b.T @ op @ C_b for op in operators]
    prop_B_bs = [C_b.T @ op @ C_s for op in operators]

    # relaxed amplitudes of the monomers in the field of each operator
    core.print_out(f"\n  ==> CPSCF of monomer A for {len(operators)} operators <==\n\n")
    tstart = time.time()
    xt_A = cpscf_batch(wfn_A, cache["eps_a"], cache["eps_r"], prop_A_ar, conv, maxiter)
    core.print_out(f"\n    ...finished in {time.time() - tstart:.2f} seconds.\n")

    core.print_out(f"\n  ==> CPSCF of monomer B for {len(operators)} operators <==\n\n")
    tstart = time.time()
    xt_B = cpscf_batch(wfn_B, cache["eps_b"], cache["eps_s"], prop_B_bs, conv, maxiter)
    core.print_out(f"\n    ...finished in {time.time() - tstart:.2f} seconds.\n")

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

    return terms
