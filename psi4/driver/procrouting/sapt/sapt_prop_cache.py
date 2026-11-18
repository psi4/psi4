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
Intermediates of the interaction-induced SAPT properties.

Ported from propSAPT v0.5.0 (https://github.com/btyrcha/propSAPT), the S^inf
second-quantized formulation of B. Tyrcha, T. Gupta, K. Patkowski and P. S. Zuchowski,
J. Chem. Theory Comput. 21, 4562 (2025). Orbital indices: a, b occupied and r, s
virtual orbitals of monomers A and B.
"""

from functools import partial

import numpy as np

from psi4 import core

from .sapt_sf_terms import _sf_compute_JK

__all__ = ["build_prop_cache"]

# the contractions keep propSAPT's index strings
contract = partial(np.einsum, optimize=True)


def _get_A_bb(s_ba, s_ab):
    """(1 - s_ba s_ab)^-1 through the eigendecomposition of s_ba s_ab."""

    gamma_bb = contract("Ba,ab->Bb", s_ba, s_ab)
    lam, U = np.linalg.eigh(gamma_bb)
    A_bb_prime = np.diag(lam / (1.0 - lam))
    return np.identity(A_bb_prime.shape[0]) + U @ A_bb_prime @ U.T


def _get_B_aa(s_ba, s_ab):
    """(1 - s_ab s_ba)^-1 through the eigendecomposition of s_ab s_ba."""

    delta_aa = contract("ba,Ab->Aa", s_ba, s_ab)
    lam, U = np.linalg.eigh(delta_aa)
    B_aa_prime = np.diag(lam / (1.0 - lam))
    return np.identity(B_aa_prime.shape[0]) + U @ B_aa_prime @ U.T


def _df_blocks(primary: core.BasisSet, df_basis: core.BasisSet, orbitals: dict) -> dict:
    """The fitted 3-index integrals (Q|pq) of the occupied-occupied, occupied-virtual and
    virtual-virtual blocks of each monomer, as (naux, p, q) arrays.

    As propSAPT's DF_BASIS_SAPT tensors: the Coulomb metric to the power -1/2, in the
    dimer-centered basis.
    """

    dfh = core.DFHelper(primary, df_basis)
    dfh.set_memory(int(0.5 * core.get_memory() / 8))
    dfh.set_method("STORE")
    dfh.set_nthreads(core.get_num_threads())
    dfh.set_MO_core(True)
    dfh.initialize()
    dfh.print_header()

    for space, C in orbitals.items():
        dfh.add_space(space, core.Matrix.from_array(C))
    blocks = ("aa", "ar", "rr", "bb", "bs", "ss")
    for block in blocks:
        dfh.add_transformation(f"Q{block}", block[0], block[1], "Qpq")
    dfh.transform()

    naux = df_basis.nbf()
    tensors = {}
    for block in blocks:
        shape = (naux, orbitals[block[0]].shape[1], orbitals[block[1]].shape[1])
        tensors[f"Q{block}"] = np.array(dfh.get_tensor(f"Q{block}")).reshape(shape)
    dfh.clear_all()

    # transposed views
    tensors["Qra"] = tensors["Qar"].transpose(0, 2, 1)
    tensors["Qsb"] = tensors["Qbs"].transpose(0, 2, 1)

    return tensors


def build_prop_cache(wfn_A: core.Wavefunction,
                     wfn_B: core.Wavefunction,
                     jk: core.JK,
                     df_basis: core.BasisSet = None) -> dict:
    """
    The orbitals, potentials and S^inf intermediates of the SAPT properties.

    Parameters
    ----------
    wfn_A, wfn_B
        Closed-shell monomer wavefunctions in the dimer-centered basis.
    jk
        J/K object of that basis. It must compute J and K.
    df_basis
        Fitting basis of the 3-index integrals of the second-order terms (DF_BASIS_SAPT).
        Without it the cache has the first-order intermediates only.

    Returns
    -------
    dict
        NumPy arrays keyed by name: the orbitals ("C_a", "C_r", "C_b", "C_s"), their
        energies ("eps_a", ...), the MO overlaps ("s_ab", ...), the S^inf intermediates
        ("A_bb", "B_aa", "E_ra", ...), the electrostatic potentials of the monomers
        ("omegaA_bs", "omegaB_ar", ...), the first-order exchange potentials, in the S^2
        approximation ("omega_exchA_bs_S2", ...) and in full ("omega_exchA_bs", ...), and
        the uncoupled induction amplitudes ("tB_ra", "tA_sb") with their denominators
        ("eps_ra", "eps_sb"). With `df_basis`, also the fitted 3-index integrals ("Qaa",
        "Qar", "Qrr", "Qbb", "Qbs", "Qss" and the transposes "Qra", "Qsb").
    """

    cache = {}

    # Orbitals and orbital energies
    C_a = np.array(wfn_A.Ca_subset("AO", "OCC"))
    C_r = np.array(wfn_A.Ca_subset("AO", "VIR"))
    C_b = np.array(wfn_B.Ca_subset("AO", "OCC"))
    C_s = np.array(wfn_B.Ca_subset("AO", "VIR"))
    orbitals = {"a": C_a, "r": C_r, "b": C_b, "s": C_s}
    cache.update({f"C_{x}": c for x, c in orbitals.items()})

    cache["eps_a"] = np.array(wfn_A.epsilon_a_subset("AO", "OCC"))
    cache["eps_r"] = np.array(wfn_A.epsilon_a_subset("AO", "VIR"))
    cache["eps_b"] = np.array(wfn_B.epsilon_a_subset("AO", "OCC"))
    cache["eps_s"] = np.array(wfn_B.epsilon_a_subset("AO", "VIR"))

    S = np.array(wfn_A.S())
    cache["S"] = S

    def s(string):
        """MO overlap matrix of the two orbital spaces of `string`."""
        return orbitals[string[0]].T @ S @ orbitals[string[1]]

    def mo(x, string):
        """AO -> MO transformation of `x` into the orbital spaces of `string`."""
        return orbitals[string[0]].T @ x @ orbitals[string[1]]

    # Nuclear potentials and Coulomb matrices of the monomers
    V_A = np.array(core.MintsHelper(wfn_A.basisset()).ao_potential())
    V_B = np.array(core.MintsHelper(wfn_B.basisset()).ao_potential())
    jk.set_do_J(True)
    jk.set_do_K(True)
    (J_A, J_B), _ = _sf_compute_JK(jk, [C_a, C_b], [C_a, C_b])

    # S^inf intermediates
    s_ab = s("ab")
    s_ba = s("ba")
    s_sa = s("sa")
    s_as = s("as")
    s_rb = s("rb")
    s_sr = s("sr")
    s_rs = s("rs")
    s_br = s("br")

    A_bb = _get_A_bb(s_ba, s_ab)
    B_aa = _get_B_aa(s_ba, s_ab)

    H_ab = contract("Bb,aB->ab", A_bb, s_ab)
    H_ba = contract("Ba,bB->ba", s_ba, A_bb)
    I_rb = contract("Bb,rB->rb", A_bb, s_rb)
    I_br = contract("Br,bB->br", s_br, A_bb)
    J_sa = contract("Aa,sA->sa", B_aa, s_sa)
    J_as = contract("As,aA->as", s_as, B_aa)

    G_sr = s_sr + contract("br,Bb,aB,sa->sr", s_br, A_bb, s_ab, s_sa)
    G_rs = s_rs + contract("as,Aa,bA,rb->rs", s_as, B_aa, s_ba, s_rb)

    C_rr = np.identity(C_r.shape[1]) - contract("bR,Bb,rB->rR", s_br, A_bb, s_rb)
    D_ss = np.identity(C_s.shape[1]) - contract("aS,Aa,sA->sS", s_as, B_aa, s_sa)

    E_ra = contract("ba,Bb,rB->ra", s_ba, A_bb, s_rb)
    E_ar = contract("aB,Bb,br->ar", s_ab, A_bb, s_br)
    F_sb = contract("ab,Aa,sA->sb", s_ab, B_aa, s_sa)
    F_bs = contract("bA,Aa,as->bs", s_ba, B_aa, s_as)

    cache.update({
        "s_ab": s_ab,
        "s_ba": s_ba,
        "s_sa": s_sa,
        "s_as": s_as,
        "s_rb": s_rb,
        "s_br": s_br,
        "s_sr": s_sr,
        "s_rs": s_rs,
        "A_bb": A_bb,
        "B_aa": B_aa,
        "H_ab": H_ab,
        "H_ba": H_ba,
        "I_rb": I_rb,
        "I_br": I_br,
        "J_sa": J_sa,
        "J_as": J_as,
        "G_sr": G_sr,
        "G_rs": G_rs,
        "C_rr": C_rr,
        "D_ss": D_ss,
        "E_ra": E_ra,
        "E_ar": E_ar,
        "F_sb": F_sb,
        "F_bs": F_bs,
    })

    # Electrostatic potentials of the monomers
    omegaA = V_A + 2 * J_A
    omegaB = V_B + 2 * J_B

    omegaA_bs = mo(omegaA, "bs")
    omegaB_ar = mo(omegaB, "ar")
    omegaA_bb = mo(omegaA, "bb")
    omegaB_aa = mo(omegaB, "aa")
    omegaA_ss = mo(omegaA, "ss")
    omegaB_rr = mo(omegaB, "rr")

    cache.update({
        "omegaA_bs": omegaA_bs,
        "omegaB_ar": omegaB_ar,
        "omegaA_bb": omegaA_bb,
        "omegaB_aa": omegaB_aa,
        "omegaA_ss": omegaA_ss,
        "omegaB_rr": omegaB_rr,
    })

    # First-order exchange potentials. Every two-electron piece is a J or K build on a
    # generalized density C_left M C_right^T, all eight in one pass of the J/K object;
    # the first four are needed for J only, the last four for K only.
    jk_left, jk_right, jk_tensor = zip(
        (C_a, C_r, E_ra.T),  # J: sum_ar (..|ar) E_ra
        (C_b, C_s, F_sb.T),  # J: sum_bs (..|bs) F_sb
        (C_a, C_r, s_ab @ s_br),  # J: sum_ar (..|ar) (s_ab s_br)_ar
        (C_b, C_s, s_ba @ s_as),  # J: sum_bs (..|bs) (s_ba s_as)_bs
        (C_a, C_s, s_as),  # K: sum_as (.a|.s) s_as
        (C_r, C_b, s_rb),  # K: sum_rb (.r|.b) s_rb
        (C_r, C_b, I_rb),  # K: sum_rb (.r|.b) I_rb
        (C_a, C_s, J_sa.T),  # K: sum_as (.a|.s) J_sa
    )
    J_list, K_list = _sf_compute_JK(jk, list(jk_left), list(jk_right), list(jk_tensor))
    J_E, J_F, J_Sar, J_Sbs = J_list[:4]
    K_Sas, K_Srb, K_I, K_J = K_list[4:]

    cache["omega_exchA_bs_S2"] = (contract("ab,Ba,bs->Bs", s_ab, s_ba, omegaA_bs) -
                                  contract("aS,sa,bs->bS", s_as, s_sa, omegaA_bs) -
                                  contract("rs,ba,ar->bs", s_rs, s_ba, omegaB_ar) - mo(K_Sas, "rb").T @ s_rs +
                                  s_ba @ mo(K_Srb, "as") - 2 * mo(J_Sar, "bs"))

    cache["omega_exchA_bs"] = (contract("Bb,sS,bs->BS", A_bb, D_ss, omegaA_bs) - contract(
        "rs,ba,ar->bs", G_rs, H_ba, omegaB_ar) - 2 * A_bb @ mo(J_E, "bs") @ D_ss - A_bb @ mo(K_J, "rb").T @ G_rs +
                               H_ba @ mo(K_I, "as") @ D_ss + 2 * H_ba @ mo(J_F, "ar") @ G_rs) - omegaA_bs

    cache["omega_exchA_sb_S2"] = (+contract("as,ba,cb->sc", s_as, s_ba, omegaA_bb) -
                                  contract("as,ba,sc->cb", s_as, s_ba, omegaA_ss) -
                                  contract("br,as,ar->sb", s_br, s_as, omegaB_ar) + s_as.T @ mo(K_Srb, "ab") -
                                  mo(K_Sas, "rs").T @ s_br.T - 2 * mo(J_Sar, "bs").T)

    cache["omega_exchA_sb"] = (-2 * mo(J_E, "bs").T + contract("bB,sb->sB", omegaA_bb, F_sb) -
                               contract("Ss,sb->Sb", omegaA_ss, F_sb) + J_sa @ mo(K_I, "ab") - mo(K_J, "rs").T @ I_rb -
                               contract("ar,sa,rb->sb", omegaB_ar, J_sa, I_rb) -
                               contract("BS,Sb,sB->sb", omegaA_bs, F_sb, F_sb) - J_sa @ mo(K_I, "as") @ F_sb -
                               F_sb @ mo(K_J, "rb").T @ I_rb + 2 * J_sa @ mo(J_F, "ar") @ I_rb +
                               2 * F_sb @ mo(J_E, "bs") @ F_sb + 2 * (F_bs @ mo(J_E, "ss")).T - 2 *
                               (mo(J_E, "bb") @ F_bs).T)

    cache["omega_exchB_ar_S2"] = (contract("ba,Ab,ar->Ar", s_ba, s_ab, omegaB_ar) -
                                  contract("bR,rb,ar->aR", s_br, s_rb, omegaB_ar) -
                                  contract("sr,ab,bs->ar", s_sr, s_ab, omegaA_bs) + s_ab @ mo(K_Sas, "rb").T -
                                  mo(K_Srb, "as") @ s_sr - 2 * mo(J_Sbs, "ar"))

    cache["omega_exchB_ar"] = (contract("Aa,rR,ar->AR", B_aa, C_rr, omegaB_ar) - contract(
        "sr,ab,bs->ar", G_sr, H_ab, omegaA_bs) - 2 * B_aa @ mo(J_F, "ar") @ C_rr - B_aa @ mo(K_I, "as") @ G_sr +
                               H_ab @ mo(K_J, "rb").T @ C_rr + 2 * H_ab @ mo(J_E, "bs") @ G_sr) - omegaB_ar

    cache["omega_exchB_ra_S2"] = (+contract("br,ab,ca->rc", s_br, s_ab, omegaB_aa) -
                                  contract("br,ab,rc->ca", s_br, s_ab, omegaB_rr) -
                                  contract("br,as,bs->ra", s_br, s_as, omegaA_bs) + s_br.T @ mo(K_Sas, "ab").T -
                                  mo(K_Srb, "rs") @ s_as.T - 2 * mo(J_Sbs, "ar").T)

    cache["omega_exchB_ra"] = (-2 * mo(J_F, "ar").T + contract("aA,ra->rA", omegaB_aa, E_ra) -
                               contract("Rr,ra->Ra", omegaB_rr, E_ra) + I_rb @ mo(K_J, "ab").T - mo(K_I, "rs") @ J_sa -
                               contract("bs,rb,sa->ra", omegaA_bs, I_rb, J_sa) -
                               contract("AR,Ra,rA->ra", omegaB_ar, E_ra, E_ra) - I_rb @ mo(K_J, "rb").T @ E_ra -
                               E_ra @ mo(K_I, "as") @ J_sa + 2 * I_rb @ mo(J_E, "bs") @ J_sa +
                               2 * E_ra @ mo(J_F, "ar") @ E_ra + 2 * (E_ar @ mo(J_F, "rr")).T - 2 *
                               (mo(J_F, "aa") @ E_ar).T)

    # Uncoupled induction amplitudes
    cache["eps_ra"] = 1 / (-cache["eps_r"].reshape(-1, 1) + cache["eps_a"])
    cache["tB_ra"] = omegaB_ar.T * cache["eps_ra"]
    cache["eps_sb"] = 1 / (-cache["eps_s"].reshape(-1, 1) + cache["eps_b"])
    cache["tA_sb"] = omegaA_bs.T * cache["eps_sb"]

    # Fitted 3-index integrals of the second-order terms
    if df_basis is not None:
        cache.update(_df_blocks(wfn_A.basisset(), df_basis, orbitals))

    return cache
