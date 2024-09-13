/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include "blas.h"
#include "mrcc.h"

namespace psi {
namespace psimrcc {

void CCMRCC::add_matrices() {
    // O^4
    wfn_->blas()->add_Matrix("<[oo]:[oo]>");
    wfn_->blas()->add_Matrix("<[oo]|[oo]>");

    // O^3V
    wfn_->blas()->add_Matrix("([oo]:[ov])");
    wfn_->blas()->add_Matrix("([oo]|[ov])");
    wfn_->blas()->add_Matrix("<[ooo]:[v]>");
    wfn_->blas()->add_Matrix("<[ooo]|[v]>");
    wfn_->blas()->add_Matrix("<[o]:[oov]>");
    wfn_->blas()->add_Matrix("<[o]|[oov]>");
    wfn_->blas()->add_Matrix("<[o]:[voo]>");
    wfn_->blas()->add_Matrix("<[o]|[voo]>");
    wfn_->blas()->add_Matrix("<[o]|[ovo]>");
    wfn_->blas()->add_Matrix("<[o]:[ovo]>");
    wfn_->blas()->add_Matrix("<[oo]:[ov]>");
    wfn_->blas()->add_Matrix("<[oo]|[ov]>");

    // O^2V^2
    wfn_->blas()->add_Matrix("<[o>o]:[v>v]>");
    wfn_->blas()->add_Matrix("<[oo]:[vv]>");
    wfn_->blas()->add_Matrix("<[oo]|[vv]>");
    wfn_->blas()->add_Matrix("<[v]:[voo]>");
    wfn_->blas()->add_Matrix("<[v]|[voo]>");
    wfn_->blas()->add_Matrix("<[o]:[ovv]>");
    wfn_->blas()->add_Matrix("<[o]|[ovv]>");
    wfn_->blas()->add_Matrix("([ov]|[ov])");
    wfn_->blas()->add_Matrix("([ov]:[ov])");
    wfn_->blas()->add_Matrix("([ov]|[vo])");
    wfn_->blas()->add_Matrix("([ov]:[vo])");
    wfn_->blas()->add_Matrix("<[ov]:[vo]>");
    wfn_->blas()->add_Matrix("<[ov]|[vo]>");
    wfn_->blas()->add_Matrix("<[vo]|[ov]>");
    wfn_->blas()->add_Matrix("<[vo]|[vo]>");
    wfn_->blas()->add_Matrix("<[ov]:[ov]>");
    wfn_->blas()->add_Matrix("<[ov]|[ov]>");

    // OV^3
    wfn_->blas()->add_Matrix("([ov]|[vv])");
    /*
      mrcc_f_int.cc:  wfn_->blas()->append("F_ae[v][v]{c} += #12# ([ov]|[vv]) 1@1 t1[ov]{c} ");
      mrcc_f_int.cc:  wfn_->blas()->append("F_ae[v][v]{o} += #12# ([ov]|[vv]) 1@1 t1[OV]{o} ");
      mrcc_f_int.cc:  wfn_->blas()->append("F_AE[V][V]{o} += #12# ([ov]|[vv]) 1@1 t1[ov]{o} ");
    */

    wfn_->blas()->add_Matrix("([ov]:[vv])");  // Used only three times in ccmrcc_f_int.cpp
                                              /*
                                                mrcc_f_int.cc:  wfn_->blas()->append("F_ae[v][v]{c} += #12# ([ov]:[vv]) 1@1 t1[ov]{c}");
                                                mrcc_f_int.cc:  wfn_->blas()->append("F_ae[v][v]{o} += #12# ([ov]:[vv]) 1@1 t1[ov]{o}");
                                                mrcc_f_int.cc:  wfn_->blas()->append("F_AE[V][V]{o} += #12# ([ov]:[vv]) 1@1 t1[OV]{o}");
                                              */

    wfn_->blas()->add_Matrix("<[vo]|[vv]>");
    /*
      mrcc_z_int.cc:  wfn_->blas()->append("Z_ijam[oov][o]{u} = #1234#   tau[oo][vv]{u} 2@2 <[vo]|[vv]>");
      mrcc_z_int.cc:  wfn_->blas()->append("Z_iJaM[oOv][O]{u} = #1234#   tau[oO][vV]{u} 2@2 <[vo]|[vv]>");
      mrcc_z_int.cc:  wfn_->blas()->append("Z_iJAm[oOV][o]{u} = #1234# - tau[oO][Vv]{u} 2@2 <[vo]|[vv]>");
      mrcc_z_int.cc:  wfn_->blas()->append("Z_IJAM[OOV][O]{u} = #1234#   tau[OO][VV]{u} 2@2 <[vo]|[vv]>");
    */

    wfn_->blas()->add_Matrix("<[v]|[ovv]>");
    /*
      mrcc_t1_amps.cc:  wfn_->blas()->append("t1_eqns[o][v]{c} +=     t2[o][OvV]{c} 2@2 <[v]|[ovv]>");
      mrcc_t1_amps.cc:  wfn_->blas()->append("t1_eqns[o][v]{o} +=     t2[o][OvV]{o} 2@2 <[v]|[ovv]>");
      mrcc_t1_amps.cc:  wfn_->blas()->append("t1_eqns[O][V]{o} +=     t2[O][oVv]{o} 2@2 <[v]|[ovv]>");
      mrcc_t2_amps.cc:  wfn_->blas()->append("t2_eqns[oO][vV]{c} += #1234#   t1[o][v]{c} 2@1 <[v]|[ovv]>");
      mrcc_t2_amps.cc:  wfn_->blas()->append("t2_eqns[oO][vV]{c} += #2143#   t1[O][V]{c} 2@1 <[v]|[ovv]>");
      mrcc_t2_amps.cc:  wfn_->blas()->append("t2_eqns[oO][vV]{o} += #1234#   t1[o][v]{o} 2@1 <[v]|[ovv]>");
      mrcc_t2_amps.cc:  wfn_->blas()->append("t2_eqns[oO][vV]{o} += #2143#   t1[O][V]{o} 2@1 <[v]|[ovv]>");
      mrcc_w_int.cc:  wfn_->blas()->append("W_jbme[ov][ov]{u} += #3241#   <[v]|[ovv]> 1@2 t1[o][v]{u}");
      mrcc_w_int.cc:  wfn_->blas()->append("W_JBme[OV][ov]{o} += #3241# <[v]|[ovv]> 1@2 t1[O][V]{o}");
      mrcc_w_int.cc:  wfn_->blas()->append("W_jbME[ov][OV]{u} += #3241# <[v]|[ovv]> 1@2 t1[o][v]{u}");
      mrcc_w_int.cc:  wfn_->blas()->append("W_JBME[OV][OV]{o} += #3241#   <[v]|[ovv]> 1@2 t1[O][V]{o}");
      mrcc_w_t3_int.cc:   wfn_->blas()->solve("W'_aBIc[vVO][v]{u}  = #4312# <[v]|[ovv]>");
      mrcc_w_t3_int.cc:   wfn_->blas()->solve("W'_AbiC[Vvo][V]{u}  = #4312# <[v]|[ovv]>");
      mrcc_w_t3_int.cc:  wfn_->blas()->solve("W_aIbC[v][OvV]{u}  = <[v]|[ovv]>");
      mrcc_w_t3_int.cc:  wfn_->blas()->solve("W_AiBc[V][oVv]{u}  = <[v]|[ovv]>");
    */

    wfn_->blas()->add_Matrix("<[v]:[ovv]>");  // Used several times
                                              /*
                                                mrcc_t1_amps.cc:  wfn_->blas()->append("t1_eqns[o][v]{c} += 1/2 t2[o][ovv]{c} 2@2 <[v]:[ovv]>");
                                                mrcc_t1_amps.cc:  wfn_->blas()->append("t1_eqns[o][v]{o} += 1/2 t2[o][ovv]{o} 2@2 <[v]:[ovv]>");
                                                mrcc_t1_amps.cc:  wfn_->blas()->append("t1_eqns[O][V]{o} += 1/2 t2[O][OVV]{o} 2@2 <[v]:[ovv]>");
                                                mrcc_t2_amps.cc:    wfn_->blas()->append("t2_eqns[oo][vv]{c} += #1234#   t1[o][v]{c} 2@1 <[v]:[ovv]>");
                                                mrcc_t2_amps.cc:    wfn_->blas()->append("t2_eqns[oo][vv]{c} += #2134# - t1[o][v]{c} 2@1 <[v]:[ovv]>");
                                                mrcc_t2_amps.cc:  wfn_->blas()->append("t2_eqns[oo][vv]{o} += #1234#   t1[o][v]{o} 2@1 <[v]:[ovv]>");
                                                mrcc_t2_amps.cc:  wfn_->blas()->append("t2_eqns[oo][vv]{o} += #2134# - t1[o][v]{o} 2@1 <[v]:[ovv]>");
                                                mrcc_t2_amps.cc:  wfn_->blas()->append("t2_eqns[OO][VV]{o} += #1234#   t1[O][V]{o} 2@1 <[v]:[ovv]>");
                                                mrcc_t2_amps.cc:  wfn_->blas()->append("t2_eqns[OO][VV]{o} += #2134# - t1[O][V]{o} 2@1 <[v]:[ovv]>");
                                                mrcc_w_t3_int.cc:   wfn_->blas()->solve("W'_abic[vvo][v]{u}  = #4312# <[v]:[ovv]>");
                                                mrcc_w_t3_int.cc:   wfn_->blas()->solve("W'_ABIC[VVO][V]{u}  = #4312# <[v]:[ovv]>");
                                                mrcc_w_t3_int.cc:  wfn_->blas()->solve("W_aibc[v][ovv]{u}  = <[v]:[ovv]>");
                                                mrcc_w_t3_int.cc:  wfn_->blas()->solve("W_AIBC[V][OVV]{u}  = <[v]:[ovv]>");
                                              */

    wfn_->blas()->add_Matrix("([vvo]|[v])");  // Used four times in ccmrcc_w_int.cpp
                                              /*
                                                mrcc_w_int.cc:  wfn_->blas()->append("W_jbme[ov][ov]{u} += #2431# - ([vvo]|[v]) 2@2 t1[o][v]{u}");
                                                mrcc_w_int.cc:  wfn_->blas()->append("W_jBmE[oV][oV]{u} += #2431# - ([vvo]|[v]) 2@2 t1[o][v]{u}");
                                                mrcc_w_int.cc:  wfn_->blas()->append("W_JbMe[Ov][Ov]{o} += #2431# - ([vvo]|[v]) 2@2 t1[O][V]{o}");
                                                mrcc_w_int.cc:  wfn_->blas()->append("W_JBME[OV][OV]{o} += #2431# - ([vvo]|[v]) 2@2 t1[O][V]{o}");
                                              */

    // V^4
    wfn_->blas()->add_Matrix("<[v>v]:[v>v]>");
    wfn_->blas()->add_Matrix("<[vv]|[v>=v]>");

    // Fock Matrix
    wfn_->blas()->add_Matrix("fock[o][o]{u}");
    wfn_->blas()->add_Matrix("fock[o][v]{u}");
    wfn_->blas()->add_Matrix("fock[v][v]{u}");
    wfn_->blas()->add_Matrix("fock[O][O]{u}");
    wfn_->blas()->add_Matrix("fock[O][V]{u}");
    wfn_->blas()->add_Matrix("fock[V][V]{u}");

    wfn_->blas()->add_Matrix("fock[v][o]{u}");
    wfn_->blas()->add_Matrix("fock[V][O]{u}");

    wfn_->blas()->add_Matrix("fock[oo]{u}");
    wfn_->blas()->add_Matrix("fock[ov]{u}");
    wfn_->blas()->add_Matrix("fock[vv]{u}");
    wfn_->blas()->add_Matrix("fock[OO]{u}");
    wfn_->blas()->add_Matrix("fock[OV]{u}");
    wfn_->blas()->add_Matrix("fock[VV]{u}");

    wfn_->blas()->add_Matrix("fock[vo]{u}");
    wfn_->blas()->add_Matrix("fock[VO]{u}");

    // Denominators
    wfn_->blas()->add_Matrix("d1[ov]{u}");
    wfn_->blas()->add_Matrix("d1[OV]{u}");
    wfn_->blas()->add_Matrix("d1[o][v]{u}");
    wfn_->blas()->add_Matrix("d1[O][V]{u}");
    wfn_->blas()->add_Matrix("d2[oo][vv]{u}");
    wfn_->blas()->add_Matrix("d2[oO][vV]{u}");
    wfn_->blas()->add_Matrix("d2[OO][VV]{u}");
    wfn_->blas()->add_Matrix("d2[o>o][v>v]{u}");
    wfn_->blas()->add_Matrix("d2[oO][vV]{u}");
    wfn_->blas()->add_Matrix("d2[O>O][V>V]{u}");
    // Shifted denominators
    wfn_->blas()->add_Matrix("d'1[o][v]{u}");
    wfn_->blas()->add_Matrix("d'1[O][V]{u}");
    wfn_->blas()->add_Matrix("d'2[oo][vv]{u}");
    wfn_->blas()->add_Matrix("d'2[oO][vV]{u}");
    wfn_->blas()->add_Matrix("d'2[OO][VV]{u}");

    // Amplitudes
    wfn_->blas()->add_Matrix("t1[ov]{u}");
    wfn_->blas()->add_Matrix("t1[OV]{u}");
    wfn_->blas()->add_Matrix("t1[o][v]{u}");
    wfn_->blas()->add_Matrix("t1[O][V]{u}");

    wfn_->blas()->add_Matrix("t2[oO][vV]{u}");
    wfn_->blas()->add_Matrix("t2[oo][vv]{u}");
    wfn_->blas()->add_Matrix("t2[OO][VV]{u}");

    wfn_->blas()->add_Matrix("t2[ov][OV]{u}");
    wfn_->blas()->add_Matrix("t2[ov][ov]{u}");
    wfn_->blas()->add_Matrix("t2[OV][OV]{u}");
    wfn_->blas()->add_Matrix("t2[oV][Ov]{u}");

    wfn_->blas()->add_Matrix("t2[o][ovv]{u}");
    wfn_->blas()->add_Matrix("t2[o][OvV]{u}");
    wfn_->blas()->add_Matrix("t2[O][oVv]{u}");
    wfn_->blas()->add_Matrix("t2[O][OVV]{u}");

    wfn_->blas()->add_Matrix("t2[v][voo]{u}");
    wfn_->blas()->add_Matrix("t2[v][VoO]{u}");
    wfn_->blas()->add_Matrix("t2[V][vOo]{u}");
    wfn_->blas()->add_Matrix("t2[V][VOO]{u}");

    wfn_->blas()->add_Matrix("t2[oo][v>v]{u}");
    wfn_->blas()->add_Matrix("t2[OO][V>V]{u}");
    wfn_->blas()->add_Matrix("t2[oO][v>=V]{u}");
    wfn_->blas()->add_Matrix("t2[oO][V>=v]{u}");

    wfn_->blas()->add_Matrix("t1_delta[o][v]{u}");
    wfn_->blas()->add_Matrix("t1_delta[O][V]{u}");
    wfn_->blas()->add_Matrix("t2_delta[oo][vv]{u}");
    wfn_->blas()->add_Matrix("t2_delta[oO][vV]{u}");
    wfn_->blas()->add_Matrix("t2_delta[OO][VV]{u}");

    // Similarity transformed Hamiltonian
    wfn_->blas()->add_Matrix("t1_eqns[o][v]{u}");
    wfn_->blas()->add_Matrix("t1_eqns[O][V]{u}");
    wfn_->blas()->add_Matrix("t2_eqns[oo][vv]{u}");
    wfn_->blas()->add_Matrix("t2_eqns[oO][vV]{u}");
    wfn_->blas()->add_Matrix("t2_eqns[OO][VV]{u}");
    wfn_->blas()->add_Matrix("t2_eqns[oo][v>v]{u}");
    wfn_->blas()->add_Matrix("t2_eqns[OO][V>V]{u}");

    // F intermediates
    wfn_->blas()->add_Matrix("F_ae[v][v]{u}");
    wfn_->blas()->add_Matrix("F_AE[V][V]{u}");
    wfn_->blas()->add_Matrix("F_mi[o][o]{u}");
    wfn_->blas()->add_Matrix("F_MI[O][O]{u}");

    // W intermediates
    wfn_->blas()->add_Matrix("W_mnij[oo][oo]{u}");
    wfn_->blas()->add_Matrix("W_mNiJ[oO][oO]{u}");
    wfn_->blas()->add_Matrix("W_MNIJ[OO][OO]{u}");

    wfn_->blas()->add_Matrix("W_jbme[ov][ov]{u}");
    wfn_->blas()->add_Matrix("W_jbME[ov][OV]{u}");
    wfn_->blas()->add_Matrix("W_JBme[OV][ov]{o}");
    wfn_->blas()->add_Matrix("W_JBME[OV][OV]{o}");
    wfn_->blas()->add_Matrix("W_jBmE[oV][oV]{u}");
    wfn_->blas()->add_Matrix("W_JbMe[Ov][Ov]{o}");

    wfn_->blas()->add_Matrix("t1_old[o][v]{u}");
    wfn_->blas()->add_Matrix("t1_old[O][V]{u}");
    wfn_->blas()->add_Matrix("t2_old[oO][vV]{u}");
    wfn_->blas()->add_Matrix("t2_old[oo][vv]{u}");
    wfn_->blas()->add_Matrix("t2_old[OO][VV]{u}");

    // MRPT2 Intermediates
    wfn_->blas()->add_Matrix("ERef{u}");
    wfn_->blas()->add_Matrix("EPT2{u}");
    wfn_->blas()->add_Matrix("Eaa{u}");
    wfn_->blas()->add_Matrix("Ebb{u}");
    wfn_->blas()->add_Matrix("Eaaaa{u}");
    wfn_->blas()->add_Matrix("Eabab{u}");
    wfn_->blas()->add_Matrix("Ebbbb{u}");

    wfn_->blas()->add_Matrix("t1[a][v]{u}");
    wfn_->blas()->add_Matrix("t1_eqns[a][a]{u}");
    wfn_->blas()->add_Matrix("t1_eqns[A][A]{o}");
    wfn_->blas()->add_Matrix("fock[a][a]{u}");
    wfn_->blas()->add_Matrix("fock[a][v]{u}");
    wfn_->blas()->add_Matrix("fock[A][A]{o}");

    wfn_->blas()->add_Matrix("factor_mk{u}");
    wfn_->blas()->add_Matrix("neg_factor_mk{u}");

    //   wfn_->blas()->add_Matrix("F_ae[v][v]{u}");
    //   wfn_->blas()->add_Matrix("F_AE[V][V]{o}");
    //
    //   wfn_->blas()->add_Matrix("F_mi[o][o]{u}");
    //   wfn_->blas()->add_Matrix("F_MI[O][O]{o}");

    wfn_->blas()->add_Matrix("F_me[o][v]{u}");
    wfn_->blas()->add_Matrix("F_me[ov]{u}");

    wfn_->blas()->add_Matrix("F_ME[O][V]{o}");
    wfn_->blas()->add_Matrix("F_ME[OV]{o}");

    wfn_->blas()->add_Matrix("F'_ae[v][v]{u}");
    wfn_->blas()->add_Matrix("F'_AE[V][V]{o}");

    wfn_->blas()->add_Matrix("F'_mi[o][o]{u}");
    wfn_->blas()->add_Matrix("F'_MI[O][O]{o}");

    wfn_->blas()->add_Matrix("tau[oo][vv]{u}");
    wfn_->blas()->add_Matrix("tau[oO][vV]{u}");
    wfn_->blas()->add_Matrix("tau[oO][Vv]{u}");
    wfn_->blas()->add_Matrix("tau[OO][VV]{u}");

    // Tau tilde intermediates (8 * O^2V^2)
    wfn_->blas()->add_Matrix("tau2[v][voo]{u}");
    wfn_->blas()->add_Matrix("tau2[o][ovv]{u}");

    wfn_->blas()->add_Matrix("tau2[V][VOO]{u}");
    wfn_->blas()->add_Matrix("tau2[O][OVV]{u}");

    wfn_->blas()->add_Matrix("tau2[v][VoO]{u}");
    wfn_->blas()->add_Matrix("tau2[o][OvV]{u}");

    wfn_->blas()->add_Matrix("tau2[V][vOo]{u}");
    wfn_->blas()->add_Matrix("tau2[O][oVv]{u}");

    wfn_->blas()->add_Matrix("tau[oo][v>v]{u}");
    wfn_->blas()->add_Matrix("tau[OO][V>V]{u}");
    wfn_->blas()->add_Matrix("tau[oO][v>=V]{u}");
    wfn_->blas()->add_Matrix("tau[oO][V>=v]{u}");

    // Amplitudes
    wfn_->blas()->add_Matrix("ECCSD{u}");
    wfn_->blas()->add_Matrix("t1_norm{u}");
    wfn_->blas()->add_Matrix("||Delta_t1||{u}");
    wfn_->blas()->add_Matrix("||Delta_t2||{u}");

    wfn_->blas()->add_Matrix("t1t1_iame[ov][ov]{u}");
    wfn_->blas()->add_Matrix("t1t1_iAMe[oV][Ov]{u}");
    wfn_->blas()->add_Matrix("t1t1_IAME[OV][OV]{u}");

    wfn_->blas()->add_Matrix("tau3[ov][ov]{u}");
    wfn_->blas()->add_Matrix("tau3[OV][OV]{u}");
    wfn_->blas()->add_Matrix("tau3[oV][vO]{u}");
    wfn_->blas()->add_Matrix("tau3[Ov][Vo]{u}");

    wfn_->blas()->add_Matrix("Z_ijam[oov][o]{u}");
    wfn_->blas()->add_Matrix("Z_iJaM[oOv][O]{u}");
    wfn_->blas()->add_Matrix("Z_iJAm[oOV][o]{u}");
    wfn_->blas()->add_Matrix("Z_IJAM[OOV][O]{u}");

    // Mukherjee Terms
    wfn_->blas()->add_Matrix("Mk1[o][v]{u}");
    wfn_->blas()->add_Matrix("Mk1[O][V]{u}");
    wfn_->blas()->add_Matrix("Mk2[oo][vv]{u}");
    wfn_->blas()->add_Matrix("Mk2[oO][vV]{u}");
    wfn_->blas()->add_Matrix("Mk2[OO][VV]{u}");

    // Triples
    if (triples_type >= ccsd_t) {
        wfn_->blas()->add_Matrix("epsilon[o][o]{u}");
        wfn_->blas()->add_Matrix("epsilon[O][O]{u}");
        wfn_->blas()->add_Matrix("epsilon[v][v]{u}");
        wfn_->blas()->add_Matrix("epsilon[V][V]{u}");

        wfn_->blas()->add_Matrix("t2[Oo][Vv]{u}");
        wfn_->blas()->add_Matrix("F_ME[O][V]{u}");

        wfn_->blas()->add_Matrix("W_ijka[oo][ov]{u}");
        wfn_->blas()->add_Matrix("W_iJkA[oO][oV]{u}");
        wfn_->blas()->add_Matrix("W_IjKa[Oo][Ov]{u}");
        wfn_->blas()->add_Matrix("W_IJKA[OO][OV]{u}");

        wfn_->blas()->add_Matrix("W_aibc[v][ovv]{u}");
        wfn_->blas()->add_Matrix("W_aIbC[v][OvV]{u}");
        wfn_->blas()->add_Matrix("W_AiBc[V][oVv]{u}");
        wfn_->blas()->add_Matrix("W_AIBC[V][OVV]{u}");
    }

    if (triples_type > ccsd_t) {  // TODO: ccsd_t should not require storage

        wfn_->blas()->add_Matrix("ERROR{u}");
        wfn_->blas()->add_Matrix("<[oo]:[ov]>");
        wfn_->blas()->add_Matrix("<[oo]|[ov]>");

        // Required by T2 * W_bcek
        wfn_->blas()->add_Matrix("t2[oov][v]{u}");
        wfn_->blas()->add_Matrix("t2[oOv][V]{u}");
        wfn_->blas()->add_Matrix("t2[OoV][v]{u}");
        wfn_->blas()->add_Matrix("t2[OOV][V]{u}");

        /*
            wfn_->blas()->add_Matrix("t2_test[oo][vv]{u}");
            wfn_->blas()->add_Matrix("t2_test[oO][vV]{u}");
            wfn_->blas()->add_Matrix("t2_test[OO][VV]{u}");*/

        wfn_->blas()->add_Matrix("W'_abic[vvo][v]{u}");
        wfn_->blas()->add_Matrix("W'_aBIc[vVO][v]{u}");
        wfn_->blas()->add_Matrix("W'_AbiC[Vvo][V]{u}");
        wfn_->blas()->add_Matrix("W'_ABIC[VVO][V]{u}");

        wfn_->blas()->add_Matrix("W'_ajki[voo][o]{u}");
        wfn_->blas()->add_Matrix("W'_AjKi[VoO][o]{u}");
        wfn_->blas()->add_Matrix("W'_aJkI[vOo][O]{u}");
        wfn_->blas()->add_Matrix("W'_AJKI[VOO][O]{u}");

        wfn_->blas()->add_Matrix("W_kija[o][oov]{u}");
        wfn_->blas()->add_Matrix("W_kiJA[o][oOV]{u}");
        wfn_->blas()->add_Matrix("W_KIja[O][Oov]{u}");
        wfn_->blas()->add_Matrix("W_KIJA[O][OOV]{u}");

        wfn_->blas()->add_Matrix("W_aibc[v][ovv]{u}");
        wfn_->blas()->add_Matrix("W_aIbC[v][OvV]{u}");
        wfn_->blas()->add_Matrix("W_AiBc[V][oVv]{u}");
        wfn_->blas()->add_Matrix("W_AIBC[V][OVV]{u}");

        //
        wfn_->blas()->add_Matrix("DELTA_t1[o][v]");
        wfn_->blas()->add_Matrix("DELTA_t1[O][V]");
        wfn_->blas()->add_Matrix("DELTA_t2[oo][vv]");
        wfn_->blas()->add_Matrix("DELTA_t2[oO][vV]");
        wfn_->blas()->add_Matrix("DELTA_t2[OO][VV]");

        // Required by T2*W_mcjk
        wfn_->blas()->add_Matrix("t2[ovv][o]{u}");
        wfn_->blas()->add_Matrix("t2[OvV][o]{u}");
        wfn_->blas()->add_Matrix("t2[oVv][O]{u}");
        wfn_->blas()->add_Matrix("t2[OVV][O]{u}");

        wfn_->blas()->add_Matrix("t3[ooo][vvv]{u}");
        wfn_->blas()->add_Matrix("t3[ooO][vvV]{u}");
        wfn_->blas()->add_Matrix("t3[oOO][vVV]{u}");
        wfn_->blas()->add_Matrix("t3[OOO][VVV]{u}");

        wfn_->blas()->add_Matrix("t3_eqns[ooo][vvv]{u}");
        wfn_->blas()->add_Matrix("t3_eqns[ooO][vvV]{u}");
        wfn_->blas()->add_Matrix("t3_eqns[oOO][vVV]{u}");
        wfn_->blas()->add_Matrix("t3_eqns[OOO][VVV]{u}");
    }

    wfn_->blas()->add_Matrix("fock[ff]{u}");
    wfn_->blas()->add_Matrix("fock[FF]{u}");
}

}  // namespace psimrcc
}  // namespace psi
