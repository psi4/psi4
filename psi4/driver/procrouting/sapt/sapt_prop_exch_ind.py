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
Second-order exchange-induction terms of the interaction-induced SAPT properties, S^2.

Ported from propSAPT v0.5.0 (https://github.com/btyrcha/propSAPT),
``prop_sapt/properties_components/exch_ind2.py``. A term is the derivative of the
exchange-induction energy with respect to the strength of a field coupled to a
one-electron operator; the contractions keep propSAPT's index strings. Orbital indices:
a, b occupied and r, s virtual orbitals of monomers A and B.
"""

from functools import partial

import numpy as np

__all__ = ["u_ra_amplitudes", "u_sb_amplitudes", "exch_ind2_s2_property"]

# the contractions keep propSAPT's index strings
contract = partial(np.einsum, optimize=True)


def u_ra_amplitudes(cache: dict, xt_A_ra: np.ndarray, xt_B_sb: np.ndarray, prop_A_aa: np.ndarray,
                    prop_A_rr: np.ndarray, t_ra: np.ndarray, coupled: bool) -> np.ndarray:
    """
    Field derivative of the induction amplitudes of monomer A, or the CPSCF right-hand
    side it solves.

    The monomers relax in the field either way: `xt_A_ra` and `xt_B_sb` are their CPSCF
    responses to the property operator, and the Fock operator changes with them.

    Parameters
    ----------
    cache
        Intermediates from :py:func:`~psi4.driver.procrouting.sapt.sapt_prop_cache.build_prop_cache`,
        with the fitted 3-index integrals.
    xt_A_ra, xt_B_sb
        Responses of monomers A and B to the property operator.
    prop_A_aa, prop_A_rr
        The property operator in the occupied and virtual orbitals of A.
    t_ra
        The induction amplitudes of A being differentiated: the CPSCF ones if `coupled`,
        else the uncoupled ones, ``cache["tB_ra"]``.
    coupled
        True for the CPSCF amplitudes: returns the (ra) right-hand side of the CPSCF
        equations of the derivative, which the caller solves. False for the uncoupled
        amplitudes, which see the electron interaction only through the Fock operator:
        returns the derivative itself, without the derivative of the two-electron
        orbital Hessian and with the orbital energy denominators in place of the solve.

    Returns
    -------
    numpy.ndarray
        (ra) array, see `coupled`.
    """

    Qaa, Qar, Qbs, Qra, Qrr, Qsb = (cache[key] for key in ("Qaa", "Qar", "Qbs", "Qra", "Qrr", "Qsb"))
    eps_ra, omegaB_aa, omegaB_rr = (cache[key] for key in ("eps_ra", "omegaB_aa", "omegaB_rr"))

    xt_B_bs = xt_B_sb.T

    # X_A -> X_A + [W_A, C_A - C_A^\dagger]_1
    gx_A_aa = (prop_A_aa + 4 * contract("ra,Qar,Qcd->cd", xt_A_ra, Qar, Qaa) -
               contract("ra,Qcr,Qad->cd", xt_A_ra, Qar, Qaa) - contract("ra,Qdr,Qac->cd", xt_A_ra, Qar, Qaa))
    gx_A_rr = (prop_A_rr + 4 * contract("ra,Qar,Qcd->cd", xt_A_ra, Qar, Qrr) -
               contract("ra,Qcr,Qad->cd", xt_A_ra, Qrr, Qar) - contract("ra,Qdr,Qac->cd", xt_A_ra, Qrr, Qar))

    u_ra = (contract("rR,Ra->ra", gx_A_rr, t_ra) - contract("Aa,rA->ra", gx_A_aa, t_ra) +
            2 * contract("bs,Qra,Qsb->ra", xt_B_bs, Qra, Qsb) + contract("Ra,rR->ra", xt_A_ra, omegaB_rr) -
            contract("rA,Aa->ra", xt_A_ra, omegaB_aa) + 2 * contract("sb,Qar,Qbs->ra", xt_B_sb, Qar, Qbs))

    if not coupled:
        # uncoupled amplitudes: no two-electron Hessian, denominators instead of CPSCF
        return u_ra * eps_ra

    # two-electron Hessian derivative acting on C_A^(1):
    #   sum_RA t_RA d/dF [4(ra|RA) - (rR|aA) - (rA|aR)]
    u_ra += (-4 * contract("rd,Qda,Qsc,sc->ra", xt_A_ra, Qaa, Qra, t_ra) +
             4 * contract("ta,Qrt,Qsc,sc->ra", xt_A_ra, Qrr, Qra, t_ra) -
             4 * contract("sd,Qra,Qdc,sc->ra", xt_A_ra, Qra, Qaa, t_ra) +
             4 * contract("tc,Qra,Qst,sc->ra", xt_A_ra, Qra, Qrr, t_ra) +
             contract("rd,Qds,Qac,sc->ra", xt_A_ra, Qar, Qaa, t_ra) +
             contract("sd,Qrd,Qac,sc->ra", xt_A_ra, Qra, Qaa, t_ra) -
             contract("ta,Qrs,Qtc,sc->ra", xt_A_ra, Qrr, Qra, t_ra) -
             contract("tc,Qrs,Qat,sc->ra", xt_A_ra, Qrr, Qar, t_ra) +
             contract("rd,Qdc,Qas,sc->ra", xt_A_ra, Qaa, Qar, t_ra) -
             contract("tc,Qrt,Qas,sc->ra", xt_A_ra, Qrr, Qar, t_ra) -
             contract("ta,Qrc,Qts,sc->ra", xt_A_ra, Qra, Qrr, t_ra) +
             contract("sd,Qrc,Qad,sc->ra", xt_A_ra, Qra, Qaa, t_ra))

    return u_ra


def u_sb_amplitudes(cache: dict, xt_A_ra: np.ndarray, xt_B_sb: np.ndarray, prop_B_bb: np.ndarray,
                    prop_B_ss: np.ndarray, t_sb: np.ndarray, coupled: bool) -> np.ndarray:
    """
    Field derivative of the induction amplitudes of monomer B, or the CPSCF right-hand
    side it solves: the monomer B counterpart of :py:func:`u_ra_amplitudes`, with the
    uncoupled amplitudes ``cache["tA_sb"]``.

    Returns
    -------
    numpy.ndarray
        (sb) array.
    """

    Qar, Qbb, Qbs, Qra, Qsb, Qss = (cache[key] for key in ("Qar", "Qbb", "Qbs", "Qra", "Qsb", "Qss"))
    eps_sb, omegaA_bb, omegaA_ss = (cache[key] for key in ("eps_sb", "omegaA_bb", "omegaA_ss"))

    xt_A_ar = xt_A_ra.T

    # X_B -> X_B + [W_B, C_B - C_B^\dagger]_1
    gx_B_bb = (prop_B_bb + 4 * contract("sb,Qbs,Qcd->cd", xt_B_sb, Qbs, Qbb) -
               contract("sb,Qcs,Qbd->cd", xt_B_sb, Qbs, Qbb) - contract("sb,Qds,Qbc->cd", xt_B_sb, Qbs, Qbb))
    gx_B_ss = (prop_B_ss + 4 * contract("sb,Qbs,Qcd->cd", xt_B_sb, Qbs, Qss) -
               contract("sb,Qcs,Qbd->cd", xt_B_sb, Qss, Qbs) - contract("sb,Qds,Qbc->cd", xt_B_sb, Qss, Qbs))

    u_sb = (contract("sS,Sb->sb", gx_B_ss, t_sb) - contract("Bb,sB->sb", gx_B_bb, t_sb) +
            2 * contract("ar,Qra,Qsb->sb", xt_A_ar, Qra, Qsb) + contract("Sb,sS->sb", xt_B_sb, omegaA_ss) -
            contract("sB,Bb->sb", xt_B_sb, omegaA_bb) + 2 * contract("ra,Qar,Qbs->sb", xt_A_ra, Qar, Qbs))

    if not coupled:
        # uncoupled amplitudes: no two-electron Hessian, denominators instead of CPSCF
        return u_sb * eps_sb

    # two-electron Hessian derivative acting on C_B^(1):
    #   sum_SB t_SB d/dF [4(sb|SB) - (sS|bB) - (sB|bS)]
    u_sb += (-4 * contract("sd,Qdb,Quc,uc->sb", xt_B_sb, Qbb, Qsb, t_sb) +
             4 * contract("tb,Qst,Quc,uc->sb", xt_B_sb, Qss, Qsb, t_sb) -
             4 * contract("ud,Qsb,Qdc,uc->sb", xt_B_sb, Qsb, Qbb, t_sb) +
             4 * contract("tc,Qsb,Qut,uc->sb", xt_B_sb, Qsb, Qss, t_sb) +
             contract("sd,Qdu,Qbc,uc->sb", xt_B_sb, Qbs, Qbb, t_sb) +
             contract("ud,Qsd,Qbc,uc->sb", xt_B_sb, Qsb, Qbb, t_sb) -
             contract("tb,Qsu,Qtc,uc->sb", xt_B_sb, Qss, Qsb, t_sb) -
             contract("tc,Qsu,Qbt,uc->sb", xt_B_sb, Qss, Qbs, t_sb) +
             contract("sd,Qdc,Qbu,uc->sb", xt_B_sb, Qbb, Qbs, t_sb) -
             contract("tc,Qst,Qbu,uc->sb", xt_B_sb, Qss, Qbs, t_sb) -
             contract("tb,Qsc,Qtu,uc->sb", xt_B_sb, Qsb, Qss, t_sb) +
             contract("ud,Qsc,Qbd,uc->sb", xt_B_sb, Qsb, Qbb, t_sb))

    return u_sb


def exch_ind2_s2_property(cache: dict, xt_A_ra: np.ndarray, xt_B_sb: np.ndarray, u_ra: np.ndarray, u_sb: np.ndarray,
                          t_ra: np.ndarray, t_sb: np.ndarray) -> float:
    """
    Field derivative of the S^2 exchange-induction energy.

    The energy is linear in the induction amplitudes, so the same expression serves both
    kinds: the coupled ones with their derivatives, or the uncoupled ones with theirs
    (:py:func:`u_ra_amplitudes`, :py:func:`u_sb_amplitudes`).

    Parameters
    ----------
    cache
        Intermediates with the fitted 3-index integrals.
    xt_A_ra, xt_B_sb
        Responses of monomers A and B to the property operator.
    u_ra, u_sb
        Field derivatives of `t_ra` and `t_sb`.
    t_ra, t_sb
        Induction amplitudes of monomers A and B.

    Returns
    -------
    float
        The property term.
    """

    Qaa, Qar, Qbb, Qbs, Qrr, Qss = (cache[key] for key in ("Qaa", "Qar", "Qbb", "Qbs", "Qrr", "Qss"))
    omegaA_bb, omegaA_bs, omegaA_ss, omegaB_aa, omegaB_ar, omegaB_rr = (cache[key]
                                                                        for key in ("omegaA_bb", "omegaA_bs",
                                                                                    "omegaA_ss", "omegaB_aa",
                                                                                    "omegaB_ar", "omegaB_rr"))

    xt_A_ar = xt_A_ra.T
    xt_B_bs = xt_B_sb.T

    s_ab = cache["s_ab"]
    s_ba = s_ab.T
    s_as = cache["s_as"]
    s_sa = s_as.T
    s_rs = cache["s_rs"]
    s_sr = s_rs.T
    s_rb = cache["s_rb"]
    s_br = s_rb.T

    x2_exch_ind_resp_s2 = (
        # < R(X) | V P2 R(V) >
        +2 * contract("bs,sr,ac,ra,cb", xt_B_bs, s_sr, s_ab, t_ra, omegaA_bb) -
        2 * contract("bs,ca,ad,sc,db", xt_B_bs, s_ba, s_ab, t_sb, omegaA_bb) +
        2 * contract("bs,ac,sa,ce,eb", xt_B_bs, s_as, s_sa, t_sb, omegaA_bb) +
        2 * contract("ar,sa,cb,rc,bs", xt_A_ar, s_sa, s_ab, t_ra, omegaA_bs) -
        2 * contract("ar,sc,rb,ca,bs", xt_A_ar, s_sr, s_rb, t_ra, omegaA_bs) +
        2 * contract("bs,sa,ac,db,cd", xt_B_bs, s_sa, s_ab, t_sb, omegaA_bs) +
        2 * contract("bs,ca,ab,sd,dc", xt_B_bs, s_sa, s_ab, t_sb, omegaA_bs) +
        2 * contract("ar,ba,rc,sb,cs", xt_A_ar, s_ba, s_rb, t_sb, omegaA_bs) -
        2 * contract("ar,rs,ca,sb,bc", xt_A_ar, s_rs, s_sa, t_sb, omegaA_bs) -
        2 * contract("bs,cr,ab,ra,sc", xt_B_bs, s_sr, s_ab, t_ra, omegaA_ss) +
        2 * contract("bs,ca,ab,dc,sd", xt_B_bs, s_ba, s_ab, t_sb, omegaA_ss) -
        2 * contract("bs,ac,da,cb,sd", xt_B_bs, s_as, s_sa, t_sb, omegaA_ss) -
        2 * contract("ar,bc,db,rd,ca", xt_A_ar, s_ba, s_ab, t_ra, omegaB_aa) +
        2 * contract("ar,bc,rb,ce,ea", xt_A_ar, s_br, s_rb, t_ra, omegaB_aa) +
        2 * contract("ar,rs,bc,sb,ca", xt_A_ar, s_rs, s_ba, t_sb, omegaB_aa) +
        2 * contract("bs,sa,cb,rc,ar", xt_B_bs, s_sa, s_ab, t_ra, omegaB_ar) +
        2 * contract("ar,ba,cb,rd,dc", xt_A_ar, s_ba, s_rb, t_ra, omegaB_ar) +
        2 * contract("ar,bc,rb,da,cd", xt_A_ar, s_ba, s_rb, t_ra, omegaB_ar) -
        2 * contract("bs,sr,cb,ra,ac", xt_B_bs, s_sr, s_rb, t_ra, omegaB_ar) +
        2 * contract("bs,ca,rb,sc,ar", xt_B_bs, s_ba, s_rb, t_sb, omegaB_ar) -
        2 * contract("bs,rc,sa,cb,ar", xt_B_bs, s_rs, s_sa, t_sb, omegaB_ar) +
        2 * contract("ar,ba,cb,dc,rd", xt_A_ar, s_ba, s_ab, t_ra, omegaB_rr) -
        2 * contract("ar,bc,db,ca,rd", xt_A_ar, s_br, s_rb, t_ra, omegaB_rr) -
        2 * contract("ar,cs,ba,sb,rc", xt_A_ar, s_rs, s_ba, t_sb, omegaB_rr) -
        2 * contract("ar,sc,db,rd,Qca,Qbs", xt_A_ar, s_sa, s_ab, t_ra, Qaa, Qbs) +
        4 * contract("ar,sc,cb,re,Qea,Qbs", xt_A_ar, s_sa, s_ab, t_ra, Qaa, Qbs) +
        2 * contract("ar,sc,rb,ce,Qea,Qbs", xt_A_ar, s_sr, s_rb, t_ra, Qaa, Qbs) -
        2 * contract("ar,bc,rd,sb,Qca,Qds", xt_A_ar, s_ba, s_rb, t_sb, Qaa, Qbs) +
        4 * contract("ar,bc,rb,sd,Qca,Qds", xt_A_ar, s_ba, s_rb, t_sb, Qaa, Qbs) +
        2 * contract("ar,rs,cd,sb,Qda,Qbc", xt_A_ar, s_rs, s_sa, t_sb, Qaa, Qbs) -
        2 * contract("bs,sa,cd,rc,Qar,Qdb", xt_B_bs, s_sa, s_ab, t_ra, Qar, Qbb) +
        4 * contract("bs,sa,ac,rd,Qdr,Qcb", xt_B_bs, s_sa, s_ab, t_ra, Qar, Qbb) +
        2 * contract("bs,sr,cd,ra,Qac,Qdb", xt_B_bs, s_sr, s_rb, t_ra, Qar, Qbb) -
        2 * contract("bs,ca,rd,sc,Qar,Qdb", xt_B_bs, s_ba, s_rb, t_sb, Qar, Qbb) +
        4 * contract("bs,ca,rc,se,Qar,Qeb", xt_B_bs, s_ba, s_rb, t_sb, Qar, Qbb) +
        2 * contract("bs,rc,sa,ce,Qar,Qeb", xt_B_bs, s_rs, s_sa, t_sb, Qar, Qbb) -
        4 * contract("ar,sa,rb,cd,Qdc,Qbs", xt_A_ar, s_sa, s_rb, t_ra, Qar, Qbs) +
        2 * contract("ar,sa,cb,rd,Qdc,Qbs", xt_A_ar, s_sa, s_rb, t_ra, Qar, Qbs) +
        2 * contract("ar,sc,rb,da,Qcd,Qbs", xt_A_ar, s_sa, s_rb, t_ra, Qar, Qbs) -
        4 * contract("bs,sa,rb,cd,Qar,Qdc", xt_B_bs, s_sa, s_rb, t_sb, Qar, Qbs) +
        2 * contract("bs,sa,rc,db,Qar,Qcd", xt_B_bs, s_sa, s_rb, t_sb, Qar, Qbs) +
        2 * contract("bs,ca,rb,sd,Qar,Qdc", xt_B_bs, s_sa, s_rb, t_sb, Qar, Qbs) +
        4 * contract("bs,ca,dc,rd,Qar,Qbs", xt_B_bs, s_ba, s_ab, t_ra, Qar, Qbs) -
        4 * contract("bs,cr,dc,ra,Qad,Qbs", xt_B_bs, s_br, s_rb, t_ra, Qar, Qbs) -
        4 * contract("bs,rc,da,cd,Qar,Qbs", xt_B_bs, s_rs, s_ba, t_sb, Qar, Qbs) -
        4 * contract("bs,ca,ab,rd,Qdr,Qsc", xt_B_bs, s_sa, s_ab, t_ra, Qar, Qss) +
        2 * contract("bs,ca,db,rd,Qar,Qsc", xt_B_bs, s_sa, s_ab, t_ra, Qar, Qss) -
        2 * contract("bs,cr,db,ra,Qad,Qsc", xt_B_bs, s_sr, s_rb, t_ra, Qar, Qss) -
        4 * contract("bs,ca,rc,eb,Qar,Qse", xt_B_bs, s_ba, s_rb, t_sb, Qar, Qss) +
        2 * contract("bs,ca,rb,dc,Qar,Qsd", xt_B_bs, s_ba, s_rb, t_sb, Qar, Qss) -
        2 * contract("bs,rc,da,cb,Qar,Qsd", xt_B_bs, s_rs, s_sa, t_sb, Qar, Qss) -
        4 * contract("ar,sc,db,cd,Qar,Qbs", xt_A_ar, s_sr, s_ab, t_ra, Qar, Qbs) +
        4 * contract("ar,bc,ce,sb,Qar,Qes", xt_A_ar, s_ba, s_ab, t_sb, Qar, Qbs) -
        4 * contract("ar,cs,dc,sb,Qar,Qbd", xt_A_ar, s_as, s_sa, t_sb, Qar, Qbs) -
        4 * contract("ar,sc,cb,ea,Qre,Qbs", xt_A_ar, s_sa, s_ab, t_ra, Qrr, Qbs) +
        2 * contract("ar,sa,cb,dc,Qrd,Qbs", xt_A_ar, s_sa, s_ab, t_ra, Qrr, Qbs) -
        2 * contract("ar,sc,db,ca,Qrd,Qbs", xt_A_ar, s_sr, s_rb, t_ra, Qrr, Qbs) -
        4 * contract("ar,ba,cb,sd,Qrc,Qds", xt_A_ar, s_ba, s_rb, t_sb, Qrr, Qbs) +
        2 * contract("ar,ba,cd,sb,Qrc,Qds", xt_A_ar, s_ba, s_rb, t_sb, Qrr, Qbs) -
        2 * contract("ar,cs,da,sb,Qrc,Qbd", xt_A_ar, s_rs, s_sa, t_sb, Qrr, Qbs)
        # < V P2 R([V, R(X)]) > + < V P2 R([X, R(V)]) >
        - 2 * contract("sr,ab,ra,bs", s_sr, s_ab, u_ra, omegaA_bs) + 2 * contract(
            "ba,cb,rc,ar", s_ba, s_ab, u_ra, omegaB_ar) - 2 * contract("br,cb,ra,ac", s_br, s_rb, u_ra, omegaB_ar) -
        4 * contract("sa,ab,rc,Qcr,Qbs", s_sa, s_ab, u_ra, Qar, Qbs) +
        2 * contract("sa,cb,rc,Qar,Qbs", s_sa, s_ab, u_ra, Qar, Qbs) -
        2 * contract("sr,cb,ra,Qac,Qbs", s_sr, s_rb, u_ra, Qar, Qbs) + 2 * contract(
            "ba,ac,sb,cs", s_ba, s_ab, u_sb, omegaA_bs) - 2 * contract("as,ca,sb,bc", s_as, s_sa, u_sb, omegaA_bs) -
        2 * contract("rs,ba,sb,ar", s_rs, s_ba, u_sb, omegaB_ar) -
        4 * contract("ba,rb,sc,Qar,Qcs", s_ba, s_rb, u_sb, Qar, Qbs) +
        2 * contract("ba,rc,sb,Qar,Qcs", s_ba, s_rb, u_sb, Qar, Qbs) -
        2 * contract("rs,ca,sb,Qar,Qbc", s_rs, s_sa, u_sb, Qar, Qbs)
        # < V P2 R(X) R(V) >
        + 2 * contract("sb,br,ac,ra,cs", xt_B_sb, s_br, s_ab, t_ra, omegaA_bs) -
        2 * contract("sb,cr,as,ra,bc", xt_B_sb, s_sr, s_as, t_ra, omegaA_bs) +
        2 * contract("ra,br,ac,sb,cs", xt_A_ra, s_br, s_ab, t_sb, omegaA_bs) -
        2 * contract("ra,sr,ac,cb,bs", xt_A_ra, s_sr, s_as, t_sb, omegaA_bs) +
        2 * contract("sb,as,ca,dc,bd", xt_B_sb, s_as, s_ba, t_sb, omegaA_bs) +
        2 * contract("sb,ac,ba,ce,es", xt_B_sb, s_as, s_ba, t_sb, omegaA_bs) +
        2 * contract("ra,br,cb,dc,ad", xt_A_ra, s_br, s_ab, t_ra, omegaB_ar) +
        2 * contract("ra,bc,ab,ce,er", xt_A_ra, s_br, s_ab, t_ra, omegaB_ar) +
        2 * contract("sb,as,bc,ra,cr", xt_B_sb, s_as, s_ba, t_ra, omegaB_ar) -
        2 * contract("sb,br,cs,ra,ac", xt_B_sb, s_br, s_rs, t_ra, omegaB_ar) +
        2 * contract("ra,as,bc,sb,cr", xt_A_ra, s_as, s_ba, t_sb, omegaB_ar) -
        2 * contract("ra,br,cs,sb,ac", xt_A_ra, s_br, s_rs, t_sb, omegaB_ar) -
        2 * contract("sb,ba,cd,rc,Qar,Qds", xt_B_sb, s_ba, s_ab, t_ra, Qar, Qbs) +
        4 * contract("sb,ba,ac,rd,Qdr,Qcs", xt_B_sb, s_ba, s_ab, t_ra, Qar, Qbs) +
        4 * contract("sb,ca,dc,rd,Qar,Qbs", xt_B_sb, s_ba, s_ab, t_ra, Qar, Qbs) -
        4 * contract("ra,sr,ab,cd,Qdc,Qbs", xt_A_ra, s_sr, s_ab, t_ra, Qar, Qbs) -
        4 * contract("ra,sc,db,cd,Qar,Qbs", xt_A_ra, s_sr, s_ab, t_ra, Qar, Qbs) +
        2 * contract("ra,sr,cb,dc,Qad,Qbs", xt_A_ra, s_sr, s_ab, t_ra, Qar, Qbs) +
        2 * contract("ra,sc,ab,ce,Qer,Qbs", xt_A_ra, s_sr, s_ab, t_ra, Qar, Qbs) -
        4 * contract("sb,cr,dc,ra,Qad,Qbs", xt_B_sb, s_br, s_rb, t_ra, Qar, Qbs) +
        2 * contract("sb,br,cd,ra,Qac,Qds", xt_B_sb, s_br, s_rb, t_ra, Qar, Qbs) -
        2 * contract("sb,cr,ds,ra,Qad,Qbc", xt_B_sb, s_sr, s_rs, t_ra, Qar, Qbs) -
        4 * contract("sb,as,ca,rd,Qdr,Qbc", xt_B_sb, s_as, s_sa, t_ra, Qar, Qbs) +
        2 * contract("sb,as,cd,ra,Qdr,Qbc", xt_B_sb, s_as, s_sa, t_ra, Qar, Qbs) -
        2 * contract("ra,bc,ad,sb,Qcr,Qds", xt_A_ra, s_ba, s_ab, t_sb, Qar, Qbs) +
        4 * contract("ra,bc,ab,sd,Qcr,Qds", xt_A_ra, s_ba, s_ab, t_sb, Qar, Qbs) +
        4 * contract("ra,bc,ce,sb,Qar,Qes", xt_A_ra, s_ba, s_ab, t_sb, Qar, Qbs) -
        4 * contract("sb,rs,ba,cd,Qar,Qdc", xt_B_sb, s_rs, s_ba, t_sb, Qar, Qbs) -
        4 * contract("sb,rc,da,cd,Qar,Qbs", xt_B_sb, s_rs, s_ba, t_sb, Qar, Qbs) +
        2 * contract("sb,rs,ca,dc,Qar,Qbd", xt_B_sb, s_rs, s_ba, t_sb, Qar, Qbs) +
        2 * contract("sb,rc,ba,ce,Qar,Qes", xt_B_sb, s_rs, s_ba, t_sb, Qar, Qbs) -
        4 * contract("ra,br,cb,sd,Qac,Qds", xt_A_ra, s_br, s_rb, t_sb, Qar, Qbs) +
        2 * contract("ra,br,cd,sb,Qac,Qds", xt_A_ra, s_br, s_rb, t_sb, Qar, Qbs) -
        2 * contract("ra,sr,cd,db,Qac,Qbs", xt_A_ra, s_sr, s_rs, t_sb, Qar, Qbs) -
        4 * contract("ra,cs,dc,sb,Qar,Qbd", xt_A_ra, s_as, s_sa, t_sb, Qar, Qbs) +
        2 * contract("ra,as,cd,sb,Qdr,Qbc", xt_A_ra, s_as, s_sa, t_sb, Qar, Qbs))

    return float(x2_exch_ind_resp_s2)
