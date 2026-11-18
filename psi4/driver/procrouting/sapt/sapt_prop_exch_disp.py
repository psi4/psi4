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
Second-order exchange-dispersion term of the interaction-induced SAPT properties, S^2.

Ported from propSAPT v0.5.0 (https://github.com/btyrcha/propSAPT),
``prop_sapt/properties_components/exch_disp2.py``. A term is the derivative of the
exchange-dispersion energy with respect to the strength of a field coupled to a
one-electron operator; the contractions keep propSAPT's index strings. The dispersion is
uncoupled. Orbital indices: a, b occupied and r, s virtual orbitals of monomers A and B.
The S^inf term is not ported: propSAPT's is not yet the derivative of its energy
(propSAPT issue #1).
"""

from functools import partial

import numpy as np

__all__ = ["u_rsab_amplitudes", "exch_disp2_s2_property"]

# the contractions keep propSAPT's index strings
contract = partial(np.einsum, optimize=True)


def u_rsab_amplitudes(cache: dict, xt_A_ra: np.ndarray, xt_B_sb: np.ndarray, prop_A_aa: np.ndarray,
                      prop_A_rr: np.ndarray, prop_B_bb: np.ndarray, prop_B_ss: np.ndarray) -> np.ndarray:
    """
    Field derivative of the uncoupled dispersion amplitudes ``cache["t_rsab"]``.

    The monomers relax in the field: `xt_A_ra` and `xt_B_sb` are their CPSCF responses
    to the property operator, and the Fock operators change with them.

    Parameters
    ----------
    cache
        Intermediates from :py:func:`~psi4.driver.procrouting.sapt.sapt_prop_cache.build_prop_cache`,
        with the fitted 3-index integrals and the dispersion amplitudes.
    xt_A_ra, xt_B_sb
        Responses of monomers A and B to the property operator.
    prop_A_aa, prop_A_rr, prop_B_bb, prop_B_ss
        The property operator in the occupied and virtual orbitals of A and B.

    Returns
    -------
    numpy.ndarray
        (rsab) array.
    """

    Qaa, Qar, Qbb, Qbs, Qrr, Qss = (cache[key] for key in ("Qaa", "Qar", "Qbb", "Qbs", "Qrr", "Qss"))
    eps_rsab, t_rsab = (cache[key] for key in ("eps_rsab", "t_rsab"))

    # X_A -> X_A + [W_A, C_A - C_A^\dagger]_1
    gx_A_aa = (prop_A_aa + 4 * contract("ra,Qar,Qcd->cd", xt_A_ra, Qar, Qaa) -
               contract("ra,Qcr,Qad->cd", xt_A_ra, Qar, Qaa) - contract("ra,Qdr,Qac->cd", xt_A_ra, Qar, Qaa))
    gx_A_rr = (prop_A_rr + 4 * contract("ra,Qar,Qcd->cd", xt_A_ra, Qar, Qrr) -
               contract("ra,Qcr,Qad->cd", xt_A_ra, Qrr, Qar) - contract("ra,Qdr,Qac->cd", xt_A_ra, Qrr, Qar))

    # X_B -> X_B + [W_B, C_B - C_B^\dagger]_1
    gx_B_bb = (prop_B_bb + 4 * contract("sb,Qbs,Qcd->cd", xt_B_sb, Qbs, Qbb) -
               contract("sb,Qcs,Qbd->cd", xt_B_sb, Qbs, Qbb) - contract("sb,Qds,Qbc->cd", xt_B_sb, Qbs, Qbb))
    gx_B_ss = (prop_B_ss + 4 * contract("sb,Qbs,Qcd->cd", xt_B_sb, Qbs, Qss) -
               contract("sb,Qcs,Qbd->cd", xt_B_sb, Qss, Qbs) - contract("sb,Qds,Qbc->cd", xt_B_sb, Qss, Qbs))

    u_rsab = (
        # R([X, R(V)])
        contract("rR,Rsab->rsab", gx_A_rr, t_rsab) - contract("Aa,rsAb->rsab", gx_A_aa, t_rsab) +
        contract("sS,rSab->rsab", gx_B_ss, t_rsab) - contract("Bb,rsaB->rsab", gx_B_bb, t_rsab)
        # R([V, R(X)])
        + contract("Ra,QrR,Qbs->rsab", xt_A_ra, Qrr, Qbs) - contract("rA,QAa,Qbs->rsab", xt_A_ra, Qaa, Qbs) +
        contract("Sb,Qar,QsS->rsab", xt_B_sb, Qar, Qss) - contract("sB,Qar,QBb->rsab", xt_B_sb, Qar, Qbb))

    u_rsab = contract("rsab,rsab->rsab", u_rsab, eps_rsab)

    return u_rsab


def exch_disp2_s2_property(cache: dict, xt_A_ra: np.ndarray, xt_B_sb: np.ndarray, prop_A_aa: np.ndarray,
                           prop_A_rr: np.ndarray, prop_B_bb: np.ndarray, prop_B_ss: np.ndarray) -> float:
    """
    Field derivative of the S^2 exchange-dispersion energy, with the uncoupled
    dispersion amplitudes.

    Parameters
    ----------
    cache
        Intermediates with the fitted 3-index integrals and the dispersion amplitudes.
    xt_A_ra, xt_B_sb
        Responses of monomers A and B to the property operator.
    prop_A_aa, prop_A_rr, prop_B_bb, prop_B_ss
        The property operator in the occupied and virtual orbitals of A and B.

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
    t_rsab = cache["t_rsab"]

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

    u_rsab = u_rsab_amplitudes(cache=cache,
                               xt_A_ra=xt_A_ra,
                               xt_B_sb=xt_B_sb,
                               prop_A_aa=prop_A_aa,
                               prop_A_rr=prop_A_rr,
                               prop_B_bb=prop_B_bb,
                               prop_B_ss=prop_B_ss)

    # < R(X) | V P2 R(V) >
    x2_exch_disp_s2 = (-2 * contract("bs,cr,ad,rsac,db", xt_B_bs, s_br, s_ab, t_rsab, omegaA_bb) +
                       4 * contract("bs,cr,ac,rsae,eb", xt_B_bs, s_br, s_ab, t_rsab, omegaA_bb) +
                       2 * contract("bs,sr,ac,rcae,eb", xt_B_bs, s_sr, s_as, t_rsab, omegaA_bb) -
                       2 * contract("ar,ba,cd,rscb,ds", xt_A_ar, s_ba, s_ab, t_rsab, omegaA_bs) +
                       4 * contract("ar,ba,cb,rsce,es", xt_A_ar, s_ba, s_ab, t_rsab, omegaA_bs) -
                       4 * contract("bs,sr,ab,rcad,dc", xt_B_bs, s_sr, s_ab, t_rsab, omegaA_bs) +
                       2 * contract("bs,sr,ac,rdab,cd", xt_B_bs, s_sr, s_ab, t_rsab, omegaA_bs) +
                       2 * contract("bs,cr,ab,rsad,dc", xt_B_bs, s_sr, s_ab, t_rsab, omegaA_bs) -
                       4 * contract("ar,bc,rb,csae,es", xt_A_ar, s_br, s_rb, t_rsab, omegaA_bs) +
                       2 * contract("ar,bc,rd,csab,ds", xt_A_ar, s_br, s_rb, t_rsab, omegaA_bs) -
                       2 * contract("ar,sc,rd,cdab,bs", xt_A_ar, s_sr, s_rs, t_rsab, omegaA_bs) +
                       2 * contract("ar,cs,da,rscb,bd", xt_A_ar, s_as, s_sa, t_rsab, omegaA_bs) -
                       4 * contract("bs,cr,ac,reab,se", xt_B_bs, s_br, s_ab, t_rsab, omegaA_ss) +
                       2 * contract("bs,cr,ab,rdac,sd", xt_B_bs, s_br, s_ab, t_rsab, omegaA_ss) -
                       2 * contract("bs,cr,ad,rdab,sc", xt_B_bs, s_sr, s_as, t_rsab, omegaA_ss) -
                       2 * contract("ar,cs,bd,rscb,da", xt_A_ar, s_as, s_ba, t_rsab, omegaB_aa) +
                       4 * contract("ar,cs,bc,rseb,ea", xt_A_ar, s_as, s_ba, t_rsab, omegaB_aa) +
                       2 * contract("ar,bc,rs,cseb,ea", xt_A_ar, s_br, s_rs, t_rsab, omegaB_aa) -
                       2 * contract("bs,ca,db,rsdc,ar", xt_B_bs, s_ba, s_ab, t_rsab, omegaB_ar) +
                       4 * contract("bs,ca,ab,rsdc,dr", xt_B_bs, s_ba, s_ab, t_rsab, omegaB_ar) -
                       4 * contract("ar,rs,ba,csdb,dc", xt_A_ar, s_rs, s_ba, t_rsab, omegaB_ar) +
                       2 * contract("ar,rs,bc,dsab,cd", xt_A_ar, s_rs, s_ba, t_rsab, omegaB_ar) +
                       2 * contract("ar,cs,ba,rsdb,dc", xt_A_ar, s_rs, s_ba, t_rsab, omegaB_ar) +
                       2 * contract("bs,cr,db,rsac,ad", xt_B_bs, s_br, s_rb, t_rsab, omegaB_ar) -
                       2 * contract("bs,sr,cd,rdab,ac", xt_B_bs, s_sr, s_rs, t_rsab, omegaB_ar) -
                       4 * contract("bs,ac,sa,rceb,er", xt_B_bs, s_as, s_sa, t_rsab, omegaB_ar) +
                       2 * contract("bs,ac,sd,rcab,dr", xt_B_bs, s_as, s_sa, t_rsab, omegaB_ar) -
                       4 * contract("ar,cs,bc,esab,re", xt_A_ar, s_as, s_ba, t_rsab, omegaB_rr) +
                       2 * contract("ar,cs,ba,dscb,rd", xt_A_ar, s_as, s_ba, t_rsab, omegaB_rr) -
                       2 * contract("ar,bc,ds,csab,rd", xt_A_ar, s_br, s_rs, t_rsab, omegaB_rr) -
                       4 * contract("ar,bc,db,rsdf,Qca,Qfs", xt_A_ar, s_ba, s_ab, t_rsab, Qaa, Qbs) -
                       4 * contract("ar,bc,ce,rsfb,Qfa,Qes", xt_A_ar, s_ba, s_ab, t_rsab, Qaa, Qbs) +
                       2 * contract("ar,bc,de,rsdb,Qca,Qes", xt_A_ar, s_ba, s_ab, t_rsab, Qaa, Qbs) -
                       2 * contract("ar,bc,rd,csfb,Qfa,Qds", xt_A_ar, s_br, s_rb, t_rsab, Qaa, Qbs) +
                       4 * contract("ar,bc,rb,csef,Qea,Qfs", xt_A_ar, s_br, s_rb, t_rsab, Qaa, Qbs) +
                       2 * contract("ar,sc,rd,cdgb,Qga,Qbs", xt_A_ar, s_sr, s_rs, t_rsab, Qaa, Qbs) -
                       2 * contract("ar,cs,de,rscb,Qea,Qbd", xt_A_ar, s_as, s_sa, t_rsab, Qaa, Qbs) +
                       4 * contract("ar,cs,dc,rsfb,Qfa,Qbd", xt_A_ar, s_as, s_sa, t_rsab, Qaa, Qbs) -
                       4 * contract("bs,ca,dc,rsdg,Qar,Qgb", xt_B_bs, s_ba, s_ab, t_rsab, Qar, Qbb) -
                       4 * contract("bs,ca,ad,rsec,Qer,Qdb", xt_B_bs, s_ba, s_ab, t_rsab, Qar, Qbb) +
                       2 * contract("bs,ca,de,rsdc,Qar,Qeb", xt_B_bs, s_ba, s_ab, t_rsab, Qar, Qbb) -
                       2 * contract("bs,cr,de,rsac,Qad,Qeb", xt_B_bs, s_br, s_rb, t_rsab, Qar, Qbb) +
                       4 * contract("bs,cr,dc,rsaf,Qad,Qfb", xt_B_bs, s_br, s_rb, t_rsab, Qar, Qbb) +
                       2 * contract("bs,sr,cd,rdaf,Qac,Qfb", xt_B_bs, s_sr, s_rs, t_rsab, Qar, Qbb) -
                       2 * contract("bs,ac,sd,rcaf,Qdr,Qfb", xt_B_bs, s_as, s_sa, t_rsab, Qar, Qbb) +
                       4 * contract("bs,ac,sa,rcef,Qer,Qfb", xt_B_bs, s_as, s_sa, t_rsab, Qar, Qbb) -
                       2 * contract("bs,sa,cd,recb,Qar,Qde", xt_B_bs, s_sa, s_ab, t_rsab, Qar, Qbs) -
                       2 * contract("bs,ca,db,rsdf,Qar,Qfc", xt_B_bs, s_sa, s_ab, t_rsab, Qar, Qbs) +
                       4 * contract("bs,sa,cb,rdcf,Qar,Qfd", xt_B_bs, s_sa, s_ab, t_rsab, Qar, Qbs) +
                       4 * contract("bs,sa,ac,rdeb,Qer,Qcd", xt_B_bs, s_sa, s_ab, t_rsab, Qar, Qbs) +
                       4 * contract("bs,ca,ab,rsde,Qdr,Qec", xt_B_bs, s_sa, s_ab, t_rsab, Qar, Qbs) -
                       2 * contract("ar,ba,cd,rseb,Qec,Qds", xt_A_ar, s_ba, s_rb, t_rsab, Qar, Qbs) -
                       2 * contract("ar,bc,rd,esab,Qce,Qds", xt_A_ar, s_ba, s_rb, t_rsab, Qar, Qbs) +
                       4 * contract("ar,ba,cb,rsde,Qdc,Qes", xt_A_ar, s_ba, s_rb, t_rsab, Qar, Qbs) +
                       4 * contract("ar,ba,rc,dseb,Qed,Qcs", xt_A_ar, s_ba, s_rb, t_rsab, Qar, Qbs) +
                       4 * contract("ar,bc,rb,dsae,Qcd,Qes", xt_A_ar, s_ba, s_rb, t_rsab, Qar, Qbs) -
                       4 * contract("bs,sr,cb,rdae,Qac,Qed", xt_B_bs, s_sr, s_rb, t_rsab, Qar, Qbs) +
                       2 * contract("bs,sr,cd,reab,Qac,Qde", xt_B_bs, s_sr, s_rb, t_rsab, Qar, Qbs) +
                       2 * contract("bs,cr,db,rsae,Qad,Qec", xt_B_bs, s_sr, s_rb, t_rsab, Qar, Qbs) -
                       4 * contract("ar,rs,ca,dseb,Qed,Qbc", xt_A_ar, s_rs, s_sa, t_rsab, Qar, Qbs) +
                       2 * contract("ar,rs,cd,esab,Qde,Qbc", xt_A_ar, s_rs, s_sa, t_rsab, Qar, Qbs) +
                       2 * contract("ar,cs,da,rseb,Qec,Qbd", xt_A_ar, s_rs, s_sa, t_rsab, Qar, Qbs) -
                       8 * contract("bs,ac,da,rcfd,Qfr,Qbs", xt_B_bs, s_as, s_ba, t_rsab, Qar, Qbs) +
                       4 * contract("bs,ac,de,rcad,Qer,Qbs", xt_B_bs, s_as, s_ba, t_rsab, Qar, Qbs) -
                       4 * contract("bs,cr,de,reac,Qad,Qbs", xt_B_bs, s_br, s_rs, t_rsab, Qar, Qbs) -
                       2 * contract("bs,ca,db,redc,Qar,Qse", xt_B_bs, s_ba, s_ab, t_rsab, Qar, Qss) +
                       4 * contract("bs,ca,ab,rdec,Qer,Qsd", xt_B_bs, s_ba, s_ab, t_rsab, Qar, Qss) +
                       4 * contract("bs,ca,dc,rfdb,Qar,Qsf", xt_B_bs, s_ba, s_ab, t_rsab, Qar, Qss) -
                       4 * contract("bs,cr,dc,rfab,Qad,Qsf", xt_B_bs, s_br, s_rb, t_rsab, Qar, Qss) +
                       2 * contract("bs,cr,db,reac,Qad,Qse", xt_B_bs, s_br, s_rb, t_rsab, Qar, Qss) -
                       2 * contract("bs,cr,de,reab,Qad,Qsc", xt_B_bs, s_sr, s_rs, t_rsab, Qar, Qss) -
                       4 * contract("bs,ac,da,rcfb,Qfr,Qsd", xt_B_bs, s_as, s_sa, t_rsab, Qar, Qss) +
                       2 * contract("bs,ac,de,rcab,Qer,Qsd", xt_B_bs, s_as, s_sa, t_rsab, Qar, Qss) -
                       8 * contract("ar,bc,db,csdg,Qar,Qgs", xt_A_ar, s_br, s_ab, t_rsab, Qar, Qbs) +
                       4 * contract("ar,bc,de,csdb,Qar,Qes", xt_A_ar, s_br, s_ab, t_rsab, Qar, Qbs) -
                       4 * contract("ar,sc,de,cedb,Qar,Qbs", xt_A_ar, s_sr, s_as, t_rsab, Qar, Qbs) -
                       2 * contract("ar,ba,cd,escb,Qre,Qds", xt_A_ar, s_ba, s_ab, t_rsab, Qrr, Qbs) +
                       4 * contract("ar,ba,cb,dscf,Qrd,Qfs", xt_A_ar, s_ba, s_ab, t_rsab, Qrr, Qbs) +
                       4 * contract("ar,bc,ce,fsab,Qrf,Qes", xt_A_ar, s_ba, s_ab, t_rsab, Qrr, Qbs) -
                       4 * contract("ar,bc,db,csaf,Qrd,Qfs", xt_A_ar, s_br, s_rb, t_rsab, Qrr, Qbs) +
                       2 * contract("ar,bc,de,csab,Qrd,Qes", xt_A_ar, s_br, s_rb, t_rsab, Qrr, Qbs) -
                       2 * contract("ar,sc,de,ceab,Qrd,Qbs", xt_A_ar, s_sr, s_rs, t_rsab, Qrr, Qbs) -
                       4 * contract("ar,cs,dc,fsab,Qrf,Qbd", xt_A_ar, s_as, s_sa, t_rsab, Qrr, Qbs) +
                       2 * contract("ar,cs,da,escb,Qre,Qbd", xt_A_ar, s_as, s_sa, t_rsab, Qrr, Qbs))

    # < V P2 R([V, R(X)]) > + < V P2 R([X, R(V)]) >
    x2_exch_disp_s2 += (-4 * contract("br,ab,rsac,cs", s_br, s_ab, u_rsab, omegaA_bs) +
                        2 * contract("br,ac,rsab,cs", s_br, s_ab, u_rsab, omegaA_bs) -
                        2 * contract("sr,ac,rcab,bs", s_sr, s_as, u_rsab, omegaA_bs) -
                        4 * contract("as,ba,rscb,cr", s_as, s_ba, u_rsab, omegaB_ar) +
                        2 * contract("as,bc,rsab,cr", s_as, s_ba, u_rsab, omegaB_ar) -
                        2 * contract("br,cs,rsab,ac", s_br, s_rs, u_rsab, omegaB_ar) -
                        2 * contract("ba,cd,rscb,Qar,Qds", s_ba, s_ab, u_rsab, Qar, Qbs) +
                        4 * contract("ba,cb,rsce,Qar,Qes", s_ba, s_ab, u_rsab, Qar, Qbs) +
                        4 * contract("ba,ac,rsdb,Qdr,Qcs", s_ba, s_ab, u_rsab, Qar, Qbs) -
                        4 * contract("br,cb,rsad,Qac,Qds", s_br, s_rb, u_rsab, Qar, Qbs) +
                        2 * contract("br,cd,rsab,Qac,Qds", s_br, s_rb, u_rsab, Qar, Qbs) -
                        2 * contract("sr,cd,rdab,Qac,Qbs", s_sr, s_rs, u_rsab, Qar, Qbs) -
                        4 * contract("as,ca,rsdb,Qdr,Qbc", s_as, s_sa, u_rsab, Qar, Qbs) +
                        2 * contract("as,cd,rsab,Qdr,Qbc", s_as, s_sa, u_rsab, Qar, Qbs))

    # < V P2 R(X) R(V) >
    x2_exch_disp_s2 += (-4 * contract("sb,br,as,rcad,dc", xt_B_sb, s_br, s_as, t_rsab, omegaA_bs) +
                        2 * contract("sb,br,ac,rcae,es", xt_B_sb, s_br, s_as, t_rsab, omegaA_bs) +
                        2 * contract("sb,cr,as,rdac,bd", xt_B_sb, s_br, s_as, t_rsab, omegaA_bs) -
                        4 * contract("ra,br,as,csdb,dc", xt_A_ra, s_br, s_as, t_rsab, omegaB_ar) +
                        2 * contract("ra,br,cs,dscb,ad", xt_A_ra, s_br, s_as, t_rsab, omegaB_ar) +
                        2 * contract("ra,bc,as,cseb,er", xt_A_ra, s_br, s_as, t_rsab, omegaB_ar) -
                        8 * contract("ra,bc,db,csdg,Qar,Qgs", xt_A_ra, s_br, s_ab, t_rsab, Qar, Qbs) -
                        2 * contract("ra,br,cd,escb,Qae,Qds", xt_A_ra, s_br, s_ab, t_rsab, Qar, Qbs) -
                        2 * contract("ra,bc,ad,csfb,Qfr,Qds", xt_A_ra, s_br, s_ab, t_rsab, Qar, Qbs) +
                        4 * contract("ra,br,cb,dscf,Qad,Qfs", xt_A_ra, s_br, s_ab, t_rsab, Qar, Qbs) +
                        4 * contract("ra,br,ac,dseb,Qed,Qcs", xt_A_ra, s_br, s_ab, t_rsab, Qar, Qbs) +
                        4 * contract("ra,bc,ab,csef,Qer,Qfs", xt_A_ra, s_br, s_ab, t_rsab, Qar, Qbs) +
                        4 * contract("ra,bc,de,csdb,Qar,Qes", xt_A_ra, s_br, s_ab, t_rsab, Qar, Qbs) -
                        4 * contract("ra,sr,ac,dcfb,Qfd,Qbs", xt_A_ra, s_sr, s_as, t_rsab, Qar, Qbs) -
                        4 * contract("ra,sc,de,cedb,Qar,Qbs", xt_A_ra, s_sr, s_as, t_rsab, Qar, Qbs) +
                        2 * contract("ra,sr,cd,edcb,Qae,Qbs", xt_A_ra, s_sr, s_as, t_rsab, Qar, Qbs) +
                        2 * contract("ra,sc,ad,cdgb,Qgr,Qbs", xt_A_ra, s_sr, s_as, t_rsab, Qar, Qbs) -
                        8 * contract("sb,ac,da,rcfd,Qfr,Qbs", xt_B_sb, s_as, s_ba, t_rsab, Qar, Qbs) -
                        2 * contract("sb,as,cd,reac,Qdr,Qbe", xt_B_sb, s_as, s_ba, t_rsab, Qar, Qbs) -
                        2 * contract("sb,ac,bd,rcaf,Qdr,Qfs", xt_B_sb, s_as, s_ba, t_rsab, Qar, Qbs) +
                        4 * contract("sb,as,ca,rdec,Qer,Qbd", xt_B_sb, s_as, s_ba, t_rsab, Qar, Qbs) +
                        4 * contract("sb,as,bc,rdae,Qcr,Qed", xt_B_sb, s_as, s_ba, t_rsab, Qar, Qbs) +
                        4 * contract("sb,ac,ba,rcef,Qer,Qfs", xt_B_sb, s_as, s_ba, t_rsab, Qar, Qbs) +
                        4 * contract("sb,ac,de,rcad,Qer,Qbs", xt_B_sb, s_as, s_ba, t_rsab, Qar, Qbs) -
                        4 * contract("sb,br,cs,rdae,Qac,Qed", xt_B_sb, s_br, s_rs, t_rsab, Qar, Qbs) -
                        4 * contract("sb,cr,de,reac,Qad,Qbs", xt_B_sb, s_br, s_rs, t_rsab, Qar, Qbs) +
                        2 * contract("sb,br,cd,rdaf,Qac,Qfs", xt_B_sb, s_br, s_rs, t_rsab, Qar, Qbs) +
                        2 * contract("sb,cr,ds,reac,Qad,Qbe", xt_B_sb, s_br, s_rs, t_rsab, Qar, Qbs))

    return float(x2_exch_disp_s2)
