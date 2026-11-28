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
Second-order exchange-induction terms of the interaction-induced SAPT properties, S^inf.

Ported from propSAPT v0.5.0 (https://github.com/btyrcha/propSAPT),
``prop_sapt/properties_components/exch_ind2_sinf.py``, the S^inf second-quantized
formulation of B. Tyrcha, T. Gupta, K. Patkowski and P. S. Zuchowski, J. Chem. Theory
Comput. 21, 4562 (2025). The contractions keep propSAPT's index strings. Orbital indices:
a, b occupied and r, s virtual orbitals of monomers A and B.
"""

from functools import partial

import numpy as np

__all__ = ["exch_ind2_sinf_property"]

# the contractions keep propSAPT's index strings
contract = partial(np.einsum, optimize=True)


def exch_ind2_sinf_property(cache: dict, xt_A_ra: np.ndarray, xt_B_sb: np.ndarray, u_ra: np.ndarray, u_sb: np.ndarray,
                            t_ra: np.ndarray, t_sb: np.ndarray) -> float:
    """
    Field derivative of the S^inf exchange-induction energy plus the induction energy.

    Subtract the induction term made with the same amplitudes to get the
    exchange-induction term. The energy is linear in the induction amplitudes, so the
    same expression serves the coupled and the uncoupled ones, with their derivatives
    (see :py:func:`~psi4.driver.procrouting.sapt.sapt_prop_exch_ind.u_ra_amplitudes`).

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
        The induction plus S^inf exchange-induction property term.
    """

    A_bb, B_aa, C_rr, D_ss, E_ar, E_ra = (cache[key] for key in ("A_bb", "B_aa", "C_rr", "D_ss", "E_ar", "E_ra"))
    F_bs, F_sb, G_rs, G_sr, H_ab, H_ba = (cache[key] for key in ("F_bs", "F_sb", "G_rs", "G_sr", "H_ab", "H_ba"))
    I_br, I_rb, J_as, J_sa, Qaa, Qar = (cache[key] for key in ("I_br", "I_rb", "J_as", "J_sa", "Qaa", "Qar"))
    Qbb, Qbs, Qrr, Qss = (cache[key] for key in ("Qbb", "Qbs", "Qrr", "Qss"))
    omegaA_bb, omegaA_bs, omegaA_ss = (cache[key] for key in ("omegaA_bb", "omegaA_bs", "omegaA_ss"))
    omegaB_aa, omegaB_ar, omegaB_rr = (cache[key] for key in ("omegaB_aa", "omegaB_ar", "omegaB_rr"))

    xt_A_ar = xt_A_ra.T
    xt_B_bs = xt_B_sb.T

    # yapf: disable
    vrx_ar = (
        - contract("Ar,aA->ar", xt_A_ar, omegaB_aa)
        + contract("aR,Rr->ar", xt_A_ar, omegaB_rr)
        + 2 * contract("bs,Qar,Qbs->ar", xt_B_bs, Qar, Qbs)
    )

    vrx_bs = (
        - contract("Bs,bB->bs", xt_B_bs, omegaA_bb)
        + contract("bS,Ss->bs", xt_B_bs, omegaA_ss)
        + 2 * contract("ar,Qar,Qbs->bs", xt_A_ar, Qar, Qbs)
    )

    vrx_abrs = (
        contract("aR,QRr,Qbs->abrs", xt_A_ar, Qrr, Qbs)
        - contract("Ar,QaA,Qbs->abrs", xt_A_ar, Qaa, Qbs)
        + contract("bS,Qar,QSs->abrs", xt_B_bs, Qar, Qss)
        - contract("Bs,Qar,QbB->abrs", xt_B_bs, Qar, Qbb)
    )

    x2_exch_ind_resp_sinf = (
        # < R(X) | V P R(V) >
        - 2 * contract("cd,ef,de,fc", G_sr, H_ab, t_ra, vrx_bs)
        - 2 * contract("cd,ef,de,fc", G_rs, H_ba, t_sb, vrx_ar)
        + 2 * contract("cd,ef,fc,de", A_bb, D_ss, t_sb, vrx_bs)
        + 2 * contract("cd,ef,fc,de", B_aa, C_rr, t_ra, vrx_ar)
        - 4 * contract("cd,ef,gh,he,dfcg", E_ra, A_bb, D_ss, t_sb, vrx_abrs)
        - 4 * contract("cd,ef,gh,fc,dheg", B_aa, C_rr, F_sb, t_ra, vrx_abrs)
        - 2 * contract("cd,ef,gh,he,fdcg", I_rb, B_aa, G_sr, t_ra, vrx_abrs)
        - 2 * contract("cd,ef,gh,fc,hdeg", A_bb, G_rs, J_sa, t_sb, vrx_abrs)
        + 2 * contract("cd,ef,gh,fg,hdce", I_rb, D_ss, H_ba, t_sb, vrx_abrs)
        + 2 * contract("cd,ef,gh,de,hfcg", C_rr, H_ab, J_sa, t_ra, vrx_abrs)
        + 4 * contract("cd,ef,gh,fg,dhce", E_ra, G_sr, H_ab, t_ra, vrx_abrs)
        + 4 * contract("cd,ef,gh,fg,hdec", F_sb, G_rs, H_ba, t_sb, vrx_abrs)
        - 2 * contract("ca,rd,ef,ar,de,fc", E_ra, G_rs, H_ba, omegaB_ar, t_sb, xt_A_ar)
        - 2 * contract("rc,de,fa,ar,ef,cd", E_ra, G_rs, H_ba, omegaB_ar, t_sb, xt_A_ar)
        - 2 * contract("cb,de,sf,bs,fd,ec", I_rb, B_aa, G_sr, omegaA_bs, t_ra, xt_A_ar)
        - 2 * contract("rc,da,ef,ar,fd,ce", I_rb, B_aa, G_sr, omegaB_ar, t_ra, xt_B_bs)
        - 2 * contract("cb,de,sf,bs,ec,fd", A_bb, G_rs, J_sa, omegaA_bs, t_sb, xt_A_ar)
        - 2 * contract("cd,re,fa,ar,ec,df", A_bb, G_rs, J_sa, omegaB_ar, t_sb, xt_B_bs)
        - 2 * contract("cb,sd,ef,bs,de,fc", F_sb, G_sr, H_ab, omegaA_bs, t_ra, xt_B_bs)
        - 2 * contract("sc,de,fb,bs,ef,cd", F_sb, G_sr, H_ab, omegaA_bs, t_ra, xt_B_bs)
        + 2 * contract("ca,de,rf,ar,fd,ec", E_ra, B_aa, C_rr, omegaB_ar, t_ra, xt_A_ar)
        + 2 * contract("rc,da,ef,ar,fd,ce", E_ra, B_aa, C_rr, omegaB_ar, t_ra, xt_A_ar)
        + 2 * contract("cb,sd,ef,bs,de,fc", I_rb, D_ss, H_ba, omegaA_bs, t_sb, xt_A_ar)
        + 2 * contract("rc,de,fa,ar,ef,cd", I_rb, D_ss, H_ba, omegaB_ar, t_sb, xt_B_bs)
        + 2 * contract("cb,de,sf,bs,ec,fd", A_bb, D_ss, F_sb, omegaA_bs, t_sb, xt_B_bs)
        + 2 * contract("cd,se,fb,bs,ec,df", A_bb, D_ss, F_sb, omegaA_bs, t_sb, xt_B_bs)
        + 2 * contract("rc,de,fa,ar,cd,ef", C_rr, H_ab, J_sa, omegaB_ar, t_ra, xt_B_bs)
        + 2 * contract("cd,eb,sf,bs,de,fc", C_rr, H_ab, J_sa, omegaA_bs, t_ra, xt_A_ar)
        - 4 * contract("ra,cb,sd,ef,de,Qar,Qbs,fc", E_ra, I_rb, D_ss, H_ba, t_sb, Qar, Qbs, xt_A_ar)
        - 4 * contract("ra,cb,de,sf,ec,Qar,Qbs,fd", E_ra, A_bb, D_ss, F_sb, t_sb, Qar, Qbs, xt_B_bs)
        - 4 * contract("ra,cd,se,fb,ec,Qar,Qbs,df", E_ra, A_bb, D_ss, F_sb, t_sb, Qar, Qbs, xt_B_bs)
        - 4 * contract("ra,cd,eb,sf,de,Qar,Qbs,fc", E_ra, C_rr, H_ab, J_sa, t_ra, Qar, Qbs, xt_A_ar)
        - 4 * contract("ca,rd,eb,sf,fe,Qar,Qbs,dc", E_ra, E_ra, A_bb, D_ss, t_sb, Qar, Qbs, xt_A_ar)
        - 4 * contract("ca,de,rf,sb,fd,Qar,Qbs,ec", E_ra, B_aa, C_rr, F_sb, t_ra, Qar, Qbs, xt_A_ar)
        - 4 * contract("rc,da,ef,sb,fd,Qar,Qbs,ce", E_ra, B_aa, C_rr, F_sb, t_ra, Qar, Qbs, xt_A_ar)
        - 4 * contract("cb,da,re,sf,ed,Qar,Qbs,fc", I_rb, B_aa, C_rr, J_sa, t_ra, Qar, Qbs, xt_A_ar)
        - 4 * contract("rc,db,se,fa,ed,Qar,Qbs,cf", I_rb, A_bb, D_ss, J_sa, t_sb, Qar, Qbs, xt_B_bs)
        - 4 * contract("rc,de,sb,fa,ef,Qar,Qbs,cd", I_rb, D_ss, F_sb, H_ba, t_sb, Qar, Qbs, xt_B_bs)
        - 4 * contract("ca,rd,eb,sf,dc,Qar,Qbs,fe", B_aa, C_rr, F_sb, F_sb, t_ra, Qar, Qbs, xt_B_bs)
        - 4 * contract("rc,sb,de,fa,cd,Qar,Qbs,ef", C_rr, F_sb, H_ab, J_sa, t_ra, Qar, Qbs, xt_B_bs)
        - 2 * contract("ca,rb,de,sf,fd,Qar,Qbs,ec", E_ra, I_rb, B_aa, G_sr, t_ra, Qar, Qbs, xt_A_ar)
        - 2 * contract("ca,db,re,sf,ed,Qar,Qbs,fc", E_ra, A_bb, G_rs, J_sa, t_sb, Qar, Qbs, xt_A_ar)
        - 2 * contract("rc,db,ea,sf,fe,Qar,Qbs,cd", E_ra, I_rb, B_aa, G_sr, t_ra, Qar, Qbs, xt_A_ar)
        - 2 * contract("rc,db,ef,sa,fd,Qar,Qbs,ce", E_ra, A_bb, G_rs, J_sa, t_sb, Qar, Qbs, xt_A_ar)
        - 2 * contract("rb,ca,sd,ef,fc,Qar,Qbs,de", I_rb, B_aa, F_sb, G_sr, t_ra, Qar, Qbs, xt_B_bs)
        - 2 * contract("rb,sc,de,fa,cd,Qar,Qbs,ef", I_rb, G_sr, H_ab, J_sa, t_ra, Qar, Qbs, xt_B_bs)
        - 2 * contract("rb,cd,ea,sf,de,Qar,Qbs,fc", I_rb, G_rs, H_ba, J_sa, t_sb, Qar, Qbs, xt_A_ar)
        - 2 * contract("cb,rd,ef,sa,de,Qar,Qbs,fc", I_rb, G_rs, H_ba, J_sa, t_sb, Qar, Qbs, xt_A_ar)
        - 2 * contract("rc,da,eb,sf,fd,Qar,Qbs,ce", I_rb, B_aa, F_sb, G_sr, t_ra, Qar, Qbs, xt_B_bs)
        - 2 * contract("rc,de,fb,sa,ef,Qar,Qbs,cd", I_rb, G_sr, H_ab, J_sa, t_ra, Qar, Qbs, xt_B_bs)
        - 2 * contract("cb,sd,re,fa,ec,Qar,Qbs,df", A_bb, F_sb, G_rs, J_sa, t_sb, Qar, Qbs, xt_B_bs)
        - 2 * contract("cd,eb,rf,sa,fc,Qar,Qbs,de", A_bb, F_sb, G_rs, J_sa, t_sb, Qar, Qbs, xt_B_bs)
        + 2 * contract("ca,rb,sd,ef,de,Qar,Qbs,fc", E_ra, I_rb, D_ss, H_ba, t_sb, Qar, Qbs, xt_A_ar)
        + 2 * contract("ca,rd,eb,sf,de,Qar,Qbs,fc", E_ra, C_rr, H_ab, J_sa, t_ra, Qar, Qbs, xt_A_ar)
        + 2 * contract("rc,db,se,fa,ef,Qar,Qbs,cd", E_ra, I_rb, D_ss, H_ba, t_sb, Qar, Qbs, xt_A_ar)
        + 2 * contract("rc,de,fb,sa,ef,Qar,Qbs,cd", E_ra, C_rr, H_ab, J_sa, t_ra, Qar, Qbs, xt_A_ar)
        + 2 * contract("rb,cd,se,fa,ec,Qar,Qbs,df", I_rb, A_bb, D_ss, J_sa, t_sb, Qar, Qbs, xt_B_bs)
        + 2 * contract("rb,ca,de,sf,ec,Qar,Qbs,fd", I_rb, B_aa, C_rr, J_sa, t_ra, Qar, Qbs, xt_A_ar)
        + 2 * contract("rb,cd,se,fa,df,Qar,Qbs,ec", I_rb, D_ss, F_sb, H_ba, t_sb, Qar, Qbs, xt_B_bs)
        + 2 * contract("cb,de,rf,sa,fd,Qar,Qbs,ec", I_rb, B_aa, C_rr, J_sa, t_ra, Qar, Qbs, xt_A_ar)
        + 2 * contract("rc,db,ef,sa,fd,Qar,Qbs,ce", I_rb, A_bb, D_ss, J_sa, t_sb, Qar, Qbs, xt_B_bs)
        + 2 * contract("rc,sd,eb,fa,df,Qar,Qbs,ce", I_rb, D_ss, F_sb, H_ba, t_sb, Qar, Qbs, xt_B_bs)
        + 2 * contract("rc,db,ef,sa,ce,Qar,Qbs,fd", C_rr, F_sb, H_ab, J_sa, t_ra, Qar, Qbs, xt_B_bs)
        + 2 * contract("rc,sd,eb,fa,ce,Qar,Qbs,df", C_rr, F_sb, H_ab, J_sa, t_ra, Qar, Qbs, xt_B_bs)
        + 4 * contract("ra,cb,de,sf,fd,Qar,Qbs,ec", E_ra, I_rb, B_aa, G_sr, t_ra, Qar, Qbs, xt_A_ar)
        + 4 * contract("ra,cb,de,sf,ec,Qar,Qbs,fd", E_ra, A_bb, G_rs, J_sa, t_sb, Qar, Qbs, xt_A_ar)
        + 4 * contract("ra,cb,sd,ef,de,Qar,Qbs,fc", E_ra, F_sb, G_sr, H_ab, t_ra, Qar, Qbs, xt_B_bs)
        + 4 * contract("ra,sc,de,fb,ef,Qar,Qbs,cd", E_ra, F_sb, G_sr, H_ab, t_ra, Qar, Qbs, xt_B_bs)
        + 4 * contract("ca,rd,se,fb,ef,Qar,Qbs,dc", E_ra, E_ra, G_sr, H_ab, t_ra, Qar, Qbs, xt_A_ar)
        + 4 * contract("ca,sb,rd,ef,de,Qar,Qbs,fc", E_ra, F_sb, G_rs, H_ba, t_sb, Qar, Qbs, xt_A_ar)
        + 4 * contract("rc,sb,de,fa,ef,Qar,Qbs,cd", E_ra, F_sb, G_rs, H_ba, t_sb, Qar, Qbs, xt_A_ar)
        + 4 * contract("cb,rd,ea,sf,de,Qar,Qbs,fc", I_rb, G_rs, H_ba, J_sa, t_sb, Qar, Qbs, xt_A_ar)
        + 4 * contract("rc,da,sb,ef,fd,Qar,Qbs,ce", I_rb, B_aa, F_sb, G_sr, t_ra, Qar, Qbs, xt_B_bs)
        + 4 * contract("rc,sd,eb,fa,de,Qar,Qbs,cf", I_rb, G_sr, H_ab, J_sa, t_ra, Qar, Qbs, xt_B_bs)
        + 4 * contract("cd,sb,re,fa,ec,Qar,Qbs,df", A_bb, F_sb, G_rs, J_sa, t_sb, Qar, Qbs, xt_B_bs)
        + 4 * contract("cb,sd,re,fa,ef,Qar,Qbs,dc", F_sb, F_sb, G_rs, H_ba, t_sb, Qar, Qbs, xt_B_bs)
        # < V P R([V, R(X)]) > + < V P R([X, R(V)]) >
        - 2 * contract("sc,db,bs,cd", G_sr, H_ab, omegaA_bs, u_ra)
        - 2 * contract("rc,da,ar,cd", G_rs, H_ba, omegaB_ar, u_sb)
        + 2 * contract("cb,sd,bs,dc", A_bb, D_ss, omegaA_bs, u_sb)
        + 2 * contract("ca,rd,ar,dc", B_aa, C_rr, omegaB_ar, u_ra)
        - 4 * contract("ra,cb,sd,dc,Qar,Qbs", E_ra, A_bb, D_ss, u_sb, Qar, Qbs)
        - 4 * contract("ca,rd,sb,dc,Qar,Qbs", B_aa, C_rr, F_sb, u_ra, Qar, Qbs)
        - 2 * contract("rb,ca,sd,dc,Qar,Qbs", I_rb, B_aa, G_sr, u_ra, Qar, Qbs)
        - 2 * contract("cb,rd,sa,dc,Qar,Qbs", A_bb, G_rs, J_sa, u_sb, Qar, Qbs)
        + 2 * contract("rb,sc,da,cd,Qar,Qbs", I_rb, D_ss, H_ba, u_sb, Qar, Qbs)
        + 2 * contract("rc,db,sa,cd,Qar,Qbs", C_rr, H_ab, J_sa, u_ra, Qar, Qbs)
        + 4 * contract("ra,sc,db,cd,Qar,Qbs", E_ra, G_sr, H_ab, u_ra, Qar, Qbs)
        + 4 * contract("sb,rc,da,cd,Qar,Qbs", F_sb, G_rs, H_ba, u_sb, Qar, Qbs)
        # < V P R(X) R(V) >
        - 2 * contract("cd,se,fb,bs,ec,df", E_ar, G_sr, H_ab, omegaA_bs, t_ra, xt_A_ra)
        - 2 * contract("cd,se,fb,bs,df,ec", E_ar, G_sr, H_ab, omegaA_bs, t_ra, xt_A_ra)
        - 2 * contract("cd,ea,rf,ar,fc,de", I_br, B_aa, G_rs, omegaB_ar, t_sb, xt_A_ra)
        - 2 * contract("cd,ea,rf,ar,de,fc", I_br, B_aa, G_rs, omegaB_ar, t_ra, xt_B_sb)
        - 2 * contract("cb,sd,ef,bs,de,fc", A_bb, G_sr, J_as, omegaA_bs, t_ra, xt_B_sb)
        - 2 * contract("cb,sd,ef,bs,fc,de", A_bb, G_sr, J_as, omegaA_bs, t_sb, xt_A_ra)
        - 2 * contract("cd,re,fa,ar,ec,df", F_bs, G_rs, H_ba, omegaB_ar, t_sb, xt_B_sb)
        - 2 * contract("cd,re,fa,ar,df,ec", F_bs, G_rs, H_ba, omegaB_ar, t_sb, xt_B_sb)
        + 2 * contract("cd,ea,rf,ar,fc,de", E_ar, B_aa, C_rr, omegaB_ar, t_ra, xt_A_ra)
        + 2 * contract("cd,ea,rf,ar,de,fc", E_ar, B_aa, C_rr, omegaB_ar, t_ra, xt_A_ra)
        + 2 * contract("cd,se,fb,bs,ec,df", I_br, D_ss, H_ab, omegaA_bs, t_sb, xt_A_ra)
        + 2 * contract("cd,se,fb,bs,df,ec", I_br, D_ss, H_ab, omegaA_bs, t_ra, xt_B_sb)
        + 2 * contract("cb,sd,ef,bs,de,fc", A_bb, D_ss, F_bs, omegaA_bs, t_sb, xt_B_sb)
        + 2 * contract("cb,sd,ef,bs,fc,de", A_bb, D_ss, F_bs, omegaA_bs, t_sb, xt_B_sb)
        + 2 * contract("rc,da,ef,ar,fd,ce", C_rr, H_ba, J_as, omegaB_ar, t_sb, xt_A_ra)
        + 2 * contract("rc,da,ef,ar,ce,fd", C_rr, H_ba, J_as, omegaB_ar, t_ra, xt_B_sb)
        - 4 * contract("cd,ea,rf,sb,fc,Qar,Qbs,de", E_ar, B_aa, C_rr, F_sb, t_ra, Qar, Qbs, xt_A_ra)
        - 4 * contract("cd,ea,rf,sb,de,Qar,Qbs,fc", E_ar, B_aa, C_rr, F_sb, t_ra, Qar, Qbs, xt_A_ra)
        - 4 * contract("ra,cd,se,fb,ec,Qar,Qbs,df", E_ra, I_br, D_ss, H_ab, t_sb, Qar, Qbs, xt_A_ra)
        - 4 * contract("ra,cd,se,fb,df,Qar,Qbs,ec", E_ra, I_br, D_ss, H_ab, t_ra, Qar, Qbs, xt_B_sb)
        - 4 * contract("ra,cb,sd,ef,de,Qar,Qbs,fc", E_ra, A_bb, D_ss, F_bs, t_sb, Qar, Qbs, xt_B_sb)
        - 4 * contract("ra,cb,sd,ef,fc,Qar,Qbs,de", E_ra, A_bb, D_ss, F_bs, t_sb, Qar, Qbs, xt_B_sb)
        - 4 * contract("cb,sd,re,fa,ef,Qar,Qbs,dc", A_bb, D_ss, G_rs, H_ba, t_sb, Qar, Qbs, xt_B_sb)
        - 4 * contract("cb,sd,re,fa,dc,Qar,Qbs,ef", A_bb, D_ss, G_rs, H_ba, t_sb, Qar, Qbs, xt_B_sb)
        - 4 * contract("ca,rd,se,fb,ef,Qar,Qbs,dc", B_aa, C_rr, G_sr, H_ab, t_ra, Qar, Qbs, xt_A_ra)
        - 4 * contract("ca,rd,se,fb,dc,Qar,Qbs,ef", B_aa, C_rr, G_sr, H_ab, t_ra, Qar, Qbs, xt_A_ra)
        - 4 * contract("rc,sb,da,ef,fd,Qar,Qbs,ce", C_rr, F_sb, H_ba, J_as, t_sb, Qar, Qbs, xt_A_ra)
        - 4 * contract("rc,sb,da,ef,ce,Qar,Qbs,fd", C_rr, F_sb, H_ba, J_as, t_ra, Qar, Qbs, xt_B_sb)
        - 2 * contract("cd,rb,ea,sf,fc,Qar,Qbs,de", E_ar, I_rb, B_aa, G_sr, t_ra, Qar, Qbs, xt_A_ra)
        - 2 * contract("cd,rb,ea,sf,de,Qar,Qbs,fc", E_ar, I_rb, B_aa, G_sr, t_ra, Qar, Qbs, xt_A_ra)
        - 2 * contract("cd,re,fb,sa,ec,Qar,Qbs,df", I_br, G_rs, H_ab, J_sa, t_sb, Qar, Qbs, xt_A_ra)
        - 2 * contract("cd,re,fb,sa,df,Qar,Qbs,ec", I_br, G_rs, H_ab, J_sa, t_ra, Qar, Qbs, xt_B_sb)
        - 2 * contract("rb,sc,da,ef,fd,Qar,Qbs,ce", I_rb, G_sr, H_ba, J_as, t_sb, Qar, Qbs, xt_A_ra)
        - 2 * contract("rb,sc,da,ef,ce,Qar,Qbs,fd", I_rb, G_sr, H_ba, J_as, t_ra, Qar, Qbs, xt_B_sb)
        - 2 * contract("cb,da,se,rf,ed,Qar,Qbs,fc", A_bb, B_aa, G_sr, G_rs, t_ra, Qar, Qbs, xt_B_sb)
        - 2 * contract("cb,de,rf,sa,fd,Qar,Qbs,ec", A_bb, F_bs, G_rs, J_sa, t_sb, Qar, Qbs, xt_B_sb)
        - 2 * contract("cb,da,se,rf,fc,Qar,Qbs,ed", A_bb, B_aa, G_sr, G_rs, t_sb, Qar, Qbs, xt_A_ra)
        - 2 * contract("cb,de,rf,sa,ec,Qar,Qbs,fd", A_bb, F_bs, G_rs, J_sa, t_sb, Qar, Qbs, xt_B_sb)
        - 2 * contract("rc,sd,ea,fb,de,Qar,Qbs,cf", C_rr, D_ss, H_ba, H_ab, t_sb, Qar, Qbs, xt_A_ra)
        - 2 * contract("rc,sd,ea,fb,cf,Qar,Qbs,de", C_rr, D_ss, H_ba, H_ab, t_ra, Qar, Qbs, xt_B_sb)
        + 2 * contract("cd,re,fb,sa,ec,Qar,Qbs,df", E_ar, C_rr, H_ab, J_sa, t_ra, Qar, Qbs, xt_A_ra)
        + 2 * contract("cd,re,fb,sa,df,Qar,Qbs,ec", E_ar, C_rr, H_ab, J_sa, t_ra, Qar, Qbs, xt_A_ra)
        + 2 * contract("cd,rb,ea,sf,fc,Qar,Qbs,de", I_br, I_rb, B_aa, D_ss, t_sb, Qar, Qbs, xt_A_ra)
        + 2 * contract("cd,rb,ea,sf,de,Qar,Qbs,fc", I_br, I_rb, B_aa, D_ss, t_ra, Qar, Qbs, xt_B_sb)
        + 2 * contract("rb,sc,de,fa,ef,Qar,Qbs,cd", I_rb, D_ss, F_bs, H_ba, t_sb, Qar, Qbs, xt_B_sb)
        + 2 * contract("rb,sc,de,fa,cd,Qar,Qbs,ef", I_rb, D_ss, F_bs, H_ba, t_sb, Qar, Qbs, xt_B_sb)
        + 2 * contract("cb,rd,ef,sa,de,Qar,Qbs,fc", A_bb, C_rr, J_as, J_sa, t_ra, Qar, Qbs, xt_B_sb)
        + 2 * contract("cb,sd,re,fa,df,Qar,Qbs,ec", A_bb, D_ss, G_rs, H_ba, t_sb, Qar, Qbs, xt_B_sb)
        + 2 * contract("cb,rd,ef,sa,fc,Qar,Qbs,de", A_bb, C_rr, J_as, J_sa, t_sb, Qar, Qbs, xt_A_ra)
        + 2 * contract("cb,sd,re,fa,ec,Qar,Qbs,df", A_bb, D_ss, G_rs, H_ba, t_sb, Qar, Qbs, xt_B_sb)
        + 2 * contract("ca,rd,se,fb,df,Qar,Qbs,ec", B_aa, C_rr, G_sr, H_ab, t_ra, Qar, Qbs, xt_A_ra)
        + 2 * contract("ca,rd,se,fb,ec,Qar,Qbs,df", B_aa, C_rr, G_sr, H_ab, t_ra, Qar, Qbs, xt_A_ra)
        + 4 * contract("cd,ra,se,fb,ec,Qar,Qbs,df", E_ar, E_ra, G_sr, H_ab, t_ra, Qar, Qbs, xt_A_ra)
        + 4 * contract("cd,ra,se,fb,df,Qar,Qbs,ec", E_ar, E_ra, G_sr, H_ab, t_ra, Qar, Qbs, xt_A_ra)
        + 4 * contract("ra,cb,sd,ef,de,Qar,Qbs,fc", E_ra, A_bb, G_sr, J_as, t_ra, Qar, Qbs, xt_B_sb)
        + 4 * contract("ra,cb,sd,ef,fc,Qar,Qbs,de", E_ra, A_bb, G_sr, J_as, t_sb, Qar, Qbs, xt_A_ra)
        + 4 * contract("cd,ea,sb,rf,fc,Qar,Qbs,de", I_br, B_aa, F_sb, G_rs, t_sb, Qar, Qbs, xt_A_ra)
        + 4 * contract("cd,ea,sb,rf,de,Qar,Qbs,fc", I_br, B_aa, F_sb, G_rs, t_ra, Qar, Qbs, xt_B_sb)
        + 4 * contract("cb,da,re,sf,ed,Qar,Qbs,fc", A_bb, B_aa, C_rr, D_ss, t_ra, Qar, Qbs, xt_B_sb)
        + 4 * contract("cb,da,re,sf,fc,Qar,Qbs,ed", A_bb, B_aa, C_rr, D_ss, t_sb, Qar, Qbs, xt_A_ra)
        + 4 * contract("cd,sb,re,fa,ec,Qar,Qbs,df", F_bs, F_sb, G_rs, H_ba, t_sb, Qar, Qbs, xt_B_sb)
        + 4 * contract("cd,sb,re,fa,df,Qar,Qbs,ec", F_bs, F_sb, G_rs, H_ba, t_sb, Qar, Qbs, xt_B_sb)
        + 4 * contract("sc,rd,ea,fb,de,Qar,Qbs,cf", G_sr, G_rs, H_ba, H_ab, t_sb, Qar, Qbs, xt_A_ra)
        + 4 * contract("sc,rd,ea,fb,cf,Qar,Qbs,de", G_sr, G_rs, H_ba, H_ab, t_ra, Qar, Qbs, xt_B_sb)
    )
    # yapf: enable

    return float(x2_exch_ind_resp_sinf)
