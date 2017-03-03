/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
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

extern FILE* outfile;

namespace psi{ namespace psimrcc{

void CCMRCC::add_matrices()
{
  // O^4
  blas->add_Matrix("<[oo]:[oo]>");
  blas->add_Matrix("<[oo]|[oo]>");

  // O^3V
  blas->add_Matrix("([oo]:[ov])");
  blas->add_Matrix("([oo]|[ov])");
  blas->add_Matrix("<[ooo]:[v]>");
  blas->add_Matrix("<[ooo]|[v]>");
  blas->add_Matrix("<[o]:[oov]>");
  blas->add_Matrix("<[o]|[oov]>");
  blas->add_Matrix("<[o]:[voo]>");
  blas->add_Matrix("<[o]|[voo]>");
  blas->add_Matrix("<[o]|[ovo]>");
  blas->add_Matrix("<[o]:[ovo]>");
  blas->add_Matrix("<[oo]:[ov]>");
  blas->add_Matrix("<[oo]|[ov]>");

  // O^2V^2
  blas->add_Matrix("<[o>o]:[v>v]>");
  blas->add_Matrix("<[oo]:[vv]>");
  blas->add_Matrix("<[oo]|[vv]>");
  blas->add_Matrix("<[v]:[voo]>");
  blas->add_Matrix("<[v]|[voo]>");
  blas->add_Matrix("<[o]:[ovv]>");
  blas->add_Matrix("<[o]|[ovv]>");
  blas->add_Matrix("([ov]|[ov])");
  blas->add_Matrix("([ov]:[ov])");
  blas->add_Matrix("([ov]|[vo])");
  blas->add_Matrix("([ov]:[vo])");
  blas->add_Matrix("<[ov]:[vo]>");
  blas->add_Matrix("<[ov]|[vo]>");
  blas->add_Matrix("<[vo]|[ov]>");
  blas->add_Matrix("<[vo]|[vo]>");
  blas->add_Matrix("<[ov]:[ov]>");
  blas->add_Matrix("<[ov]|[ov]>");

  // OV^3
  blas->add_Matrix("([ov]|[vv])");
/*
  mrcc_f_int.cc:  blas->append("F_ae[v][v]{c} += #12# ([ov]|[vv]) 1@1 t1[ov]{c} ");
  mrcc_f_int.cc:  blas->append("F_ae[v][v]{o} += #12# ([ov]|[vv]) 1@1 t1[OV]{o} ");
  mrcc_f_int.cc:  blas->append("F_AE[V][V]{o} += #12# ([ov]|[vv]) 1@1 t1[ov]{o} ");
*/

  blas->add_Matrix("([ov]:[vv])");  // Used only three times in ccmrcc_f_int.cpp
/*
  mrcc_f_int.cc:  blas->append("F_ae[v][v]{c} += #12# ([ov]:[vv]) 1@1 t1[ov]{c}");
  mrcc_f_int.cc:  blas->append("F_ae[v][v]{o} += #12# ([ov]:[vv]) 1@1 t1[ov]{o}");
  mrcc_f_int.cc:  blas->append("F_AE[V][V]{o} += #12# ([ov]:[vv]) 1@1 t1[OV]{o}");
*/

  blas->add_Matrix("<[vo]|[vv]>");
/*
  mrcc_z_int.cc:  blas->append("Z_ijam[oov][o]{u} = #1234#   tau[oo][vv]{u} 2@2 <[vo]|[vv]>");
  mrcc_z_int.cc:  blas->append("Z_iJaM[oOv][O]{u} = #1234#   tau[oO][vV]{u} 2@2 <[vo]|[vv]>");
  mrcc_z_int.cc:  blas->append("Z_iJAm[oOV][o]{u} = #1234# - tau[oO][Vv]{u} 2@2 <[vo]|[vv]>");
  mrcc_z_int.cc:  blas->append("Z_IJAM[OOV][O]{u} = #1234#   tau[OO][VV]{u} 2@2 <[vo]|[vv]>");
*/

  blas->add_Matrix("<[v]|[ovv]>");
/*
  mrcc_t1_amps.cc:  blas->append("t1_eqns[o][v]{c} +=     t2[o][OvV]{c} 2@2 <[v]|[ovv]>");
  mrcc_t1_amps.cc:  blas->append("t1_eqns[o][v]{o} +=     t2[o][OvV]{o} 2@2 <[v]|[ovv]>");
  mrcc_t1_amps.cc:  blas->append("t1_eqns[O][V]{o} +=     t2[O][oVv]{o} 2@2 <[v]|[ovv]>");
  mrcc_t2_amps.cc:  blas->append("t2_eqns[oO][vV]{c} += #1234#   t1[o][v]{c} 2@1 <[v]|[ovv]>");
  mrcc_t2_amps.cc:  blas->append("t2_eqns[oO][vV]{c} += #2143#   t1[O][V]{c} 2@1 <[v]|[ovv]>");
  mrcc_t2_amps.cc:  blas->append("t2_eqns[oO][vV]{o} += #1234#   t1[o][v]{o} 2@1 <[v]|[ovv]>");
  mrcc_t2_amps.cc:  blas->append("t2_eqns[oO][vV]{o} += #2143#   t1[O][V]{o} 2@1 <[v]|[ovv]>");
  mrcc_w_int.cc:  blas->append("W_jbme[ov][ov]{u} += #3241#   <[v]|[ovv]> 1@2 t1[o][v]{u}");
  mrcc_w_int.cc:  blas->append("W_JBme[OV][ov]{o} += #3241# <[v]|[ovv]> 1@2 t1[O][V]{o}");
  mrcc_w_int.cc:  blas->append("W_jbME[ov][OV]{u} += #3241# <[v]|[ovv]> 1@2 t1[o][v]{u}");
  mrcc_w_int.cc:  blas->append("W_JBME[OV][OV]{o} += #3241#   <[v]|[ovv]> 1@2 t1[O][V]{o}");
  mrcc_w_t3_int.cc:   blas->solve("W'_aBIc[vVO][v]{u}  = #4312# <[v]|[ovv]>");
  mrcc_w_t3_int.cc:   blas->solve("W'_AbiC[Vvo][V]{u}  = #4312# <[v]|[ovv]>");
  mrcc_w_t3_int.cc:  blas->solve("W_aIbC[v][OvV]{u}  = <[v]|[ovv]>");
  mrcc_w_t3_int.cc:  blas->solve("W_AiBc[V][oVv]{u}  = <[v]|[ovv]>");
*/

  blas->add_Matrix("<[v]:[ovv]>");  // Used several times
/*
  mrcc_t1_amps.cc:  blas->append("t1_eqns[o][v]{c} += 1/2 t2[o][ovv]{c} 2@2 <[v]:[ovv]>");
  mrcc_t1_amps.cc:  blas->append("t1_eqns[o][v]{o} += 1/2 t2[o][ovv]{o} 2@2 <[v]:[ovv]>");
  mrcc_t1_amps.cc:  blas->append("t1_eqns[O][V]{o} += 1/2 t2[O][OVV]{o} 2@2 <[v]:[ovv]>");
  mrcc_t2_amps.cc:    blas->append("t2_eqns[oo][vv]{c} += #1234#   t1[o][v]{c} 2@1 <[v]:[ovv]>");
  mrcc_t2_amps.cc:    blas->append("t2_eqns[oo][vv]{c} += #2134# - t1[o][v]{c} 2@1 <[v]:[ovv]>");
  mrcc_t2_amps.cc:  blas->append("t2_eqns[oo][vv]{o} += #1234#   t1[o][v]{o} 2@1 <[v]:[ovv]>");
  mrcc_t2_amps.cc:  blas->append("t2_eqns[oo][vv]{o} += #2134# - t1[o][v]{o} 2@1 <[v]:[ovv]>");
  mrcc_t2_amps.cc:  blas->append("t2_eqns[OO][VV]{o} += #1234#   t1[O][V]{o} 2@1 <[v]:[ovv]>");
  mrcc_t2_amps.cc:  blas->append("t2_eqns[OO][VV]{o} += #2134# - t1[O][V]{o} 2@1 <[v]:[ovv]>");
  mrcc_w_t3_int.cc:   blas->solve("W'_abic[vvo][v]{u}  = #4312# <[v]:[ovv]>");
  mrcc_w_t3_int.cc:   blas->solve("W'_ABIC[VVO][V]{u}  = #4312# <[v]:[ovv]>");
  mrcc_w_t3_int.cc:  blas->solve("W_aibc[v][ovv]{u}  = <[v]:[ovv]>");
  mrcc_w_t3_int.cc:  blas->solve("W_AIBC[V][OVV]{u}  = <[v]:[ovv]>");
*/

  blas->add_Matrix("([vvo]|[v])");  // Used four times in ccmrcc_w_int.cpp
/*
  mrcc_w_int.cc:  blas->append("W_jbme[ov][ov]{u} += #2431# - ([vvo]|[v]) 2@2 t1[o][v]{u}");
  mrcc_w_int.cc:  blas->append("W_jBmE[oV][oV]{u} += #2431# - ([vvo]|[v]) 2@2 t1[o][v]{u}");
  mrcc_w_int.cc:  blas->append("W_JbMe[Ov][Ov]{o} += #2431# - ([vvo]|[v]) 2@2 t1[O][V]{o}");
  mrcc_w_int.cc:  blas->append("W_JBME[OV][OV]{o} += #2431# - ([vvo]|[v]) 2@2 t1[O][V]{o}");
*/

  // V^4
  blas->add_Matrix("<[v>v]:[v>v]>");
  blas->add_Matrix("<[vv]|[v>=v]>");

  // Fock Matrix
  blas->add_Matrix("fock[o][o]{u}");
  blas->add_Matrix("fock[o][v]{u}");
  blas->add_Matrix("fock[v][v]{u}");
  blas->add_Matrix("fock[O][O]{u}");
  blas->add_Matrix("fock[O][V]{u}");
  blas->add_Matrix("fock[V][V]{u}");

  blas->add_Matrix("fock[v][o]{u}");
  blas->add_Matrix("fock[V][O]{u}");

  blas->add_Matrix("fock[oo]{u}");
  blas->add_Matrix("fock[ov]{u}");
  blas->add_Matrix("fock[vv]{u}");
  blas->add_Matrix("fock[OO]{u}");
  blas->add_Matrix("fock[OV]{u}");
  blas->add_Matrix("fock[VV]{u}");

  blas->add_Matrix("fock[vo]{u}");
  blas->add_Matrix("fock[VO]{u}");

  // Denominators
  blas->add_Matrix("d1[ov]{u}");
  blas->add_Matrix("d1[OV]{u}");
  blas->add_Matrix("d1[o][v]{u}");
  blas->add_Matrix("d1[O][V]{u}");
  blas->add_Matrix("d2[oo][vv]{u}");
  blas->add_Matrix("d2[oO][vV]{u}");
  blas->add_Matrix("d2[OO][VV]{u}");
  blas->add_Matrix("d2[o>o][v>v]{u}");
  blas->add_Matrix("d2[oO][vV]{u}");
  blas->add_Matrix("d2[O>O][V>V]{u}");
  // Shifted denominators
  blas->add_Matrix("d'1[o][v]{u}");
  blas->add_Matrix("d'1[O][V]{u}");
  blas->add_Matrix("d'2[oo][vv]{u}");
  blas->add_Matrix("d'2[oO][vV]{u}");
  blas->add_Matrix("d'2[OO][VV]{u}");

  // Amplitudes
  blas->add_Matrix("t1[ov]{u}");
  blas->add_Matrix("t1[OV]{u}");
  blas->add_Matrix("t1[o][v]{u}");
  blas->add_Matrix("t1[O][V]{u}");

  blas->add_Matrix("t2[oO][vV]{u}");
  blas->add_Matrix("t2[oo][vv]{u}");
  blas->add_Matrix("t2[OO][VV]{u}");

  blas->add_Matrix("t2[ov][OV]{u}");
  blas->add_Matrix("t2[ov][ov]{u}");
  blas->add_Matrix("t2[OV][OV]{u}");
  blas->add_Matrix("t2[oV][Ov]{u}");

  blas->add_Matrix("t2[o][ovv]{u}");
  blas->add_Matrix("t2[o][OvV]{u}");
  blas->add_Matrix("t2[O][oVv]{u}");
  blas->add_Matrix("t2[O][OVV]{u}");

  blas->add_Matrix("t2[v][voo]{u}");
  blas->add_Matrix("t2[v][VoO]{u}");
  blas->add_Matrix("t2[V][vOo]{u}");
  blas->add_Matrix("t2[V][VOO]{u}");

  blas->add_Matrix("t2[oo][v>v]{u}");
  blas->add_Matrix("t2[OO][V>V]{u}");
  blas->add_Matrix("t2[oO][v>=V]{u}");
  blas->add_Matrix("t2[oO][V>=v]{u}");

  blas->add_Matrix("t1_delta[o][v]{u}");
  blas->add_Matrix("t1_delta[O][V]{u}");
  blas->add_Matrix("t2_delta[oo][vv]{u}");
  blas->add_Matrix("t2_delta[oO][vV]{u}");
  blas->add_Matrix("t2_delta[OO][VV]{u}");

  // Similarity transformed Hamiltonian
  blas->add_Matrix("t1_eqns[o][v]{u}");
  blas->add_Matrix("t1_eqns[O][V]{u}");
  blas->add_Matrix("t2_eqns[oo][vv]{u}");
  blas->add_Matrix("t2_eqns[oO][vV]{u}");
  blas->add_Matrix("t2_eqns[OO][VV]{u}");
  blas->add_Matrix("t2_eqns[oo][v>v]{u}");
  blas->add_Matrix("t2_eqns[OO][V>V]{u}");

  // F intermediates
  blas->add_Matrix("F_ae[v][v]{u}");
  blas->add_Matrix("F_AE[V][V]{u}");
  blas->add_Matrix("F_mi[o][o]{u}");
  blas->add_Matrix("F_MI[O][O]{u}");

  // W intermediates
  blas->add_Matrix("W_mnij[oo][oo]{u}");
  blas->add_Matrix("W_mNiJ[oO][oO]{u}");
  blas->add_Matrix("W_MNIJ[OO][OO]{u}");

  blas->add_Matrix("W_jbme[ov][ov]{u}");
  blas->add_Matrix("W_jbME[ov][OV]{u}");
  blas->add_Matrix("W_JBme[OV][ov]{o}");
  blas->add_Matrix("W_JBME[OV][OV]{o}");
  blas->add_Matrix("W_jBmE[oV][oV]{u}");
  blas->add_Matrix("W_JbMe[Ov][Ov]{o}");

  blas->add_Matrix("t1_old[o][v]{u}");
  blas->add_Matrix("t1_old[O][V]{u}");
  blas->add_Matrix("t2_old[oO][vV]{u}");
  blas->add_Matrix("t2_old[oo][vv]{u}");
  blas->add_Matrix("t2_old[OO][VV]{u}");

  // MRPT2 Intermediates
  blas->add_Matrix("ERef{u}");
  blas->add_Matrix("EPT2{u}");
  blas->add_Matrix("Eaa{u}");
  blas->add_Matrix("Ebb{u}");
  blas->add_Matrix("Eaaaa{u}");
  blas->add_Matrix("Eabab{u}");
  blas->add_Matrix("Ebbbb{u}");

  blas->add_Matrix("t1[a][v]{u}");
  blas->add_Matrix("t1_eqns[a][a]{u}");
  blas->add_Matrix("t1_eqns[A][A]{o}");
  blas->add_Matrix("fock[a][a]{u}");
  blas->add_Matrix("fock[a][v]{u}");
  blas->add_Matrix("fock[A][A]{o}");

  blas->add_Matrix("factor_mk{u}");
  blas->add_Matrix("neg_factor_mk{u}");

//   blas->add_Matrix("F_ae[v][v]{u}");
//   blas->add_Matrix("F_AE[V][V]{o}");
//
//   blas->add_Matrix("F_mi[o][o]{u}");
//   blas->add_Matrix("F_MI[O][O]{o}");

  blas->add_Matrix("F_me[o][v]{u}");
  blas->add_Matrix("F_me[ov]{u}");

  blas->add_Matrix("F_ME[O][V]{o}");
  blas->add_Matrix("F_ME[OV]{o}");

  blas->add_Matrix("F'_ae[v][v]{u}");
  blas->add_Matrix("F'_AE[V][V]{o}");

  blas->add_Matrix("F'_mi[o][o]{u}");
  blas->add_Matrix("F'_MI[O][O]{o}");

  blas->add_Matrix("tau[oo][vv]{u}");
  blas->add_Matrix("tau[oO][vV]{u}");
  blas->add_Matrix("tau[oO][Vv]{u}");
  blas->add_Matrix("tau[OO][VV]{u}");

  // Tau tilde intermediates (8 * O^2V^2)
  blas->add_Matrix("tau2[v][voo]{u}");
  blas->add_Matrix("tau2[o][ovv]{u}");

  blas->add_Matrix("tau2[V][VOO]{u}");
  blas->add_Matrix("tau2[O][OVV]{u}");

  blas->add_Matrix("tau2[v][VoO]{u}");
  blas->add_Matrix("tau2[o][OvV]{u}");

  blas->add_Matrix("tau2[V][vOo]{u}");
  blas->add_Matrix("tau2[O][oVv]{u}");

  blas->add_Matrix("tau[oo][v>v]{u}");
  blas->add_Matrix("tau[OO][V>V]{u}");
  blas->add_Matrix("tau[oO][v>=V]{u}");
  blas->add_Matrix("tau[oO][V>=v]{u}");

  // Amplitudes
  blas->add_Matrix("ECCSD{u}");
  blas->add_Matrix("t1_norm{u}");
  blas->add_Matrix("||Delta_t1||{u}");
  blas->add_Matrix("||Delta_t2||{u}");


  blas->add_Matrix("t1t1_iame[ov][ov]{u}");
  blas->add_Matrix("t1t1_iAMe[oV][Ov]{u}");
  blas->add_Matrix("t1t1_IAME[OV][OV]{u}");

  blas->add_Matrix("tau3[ov][ov]{u}");
  blas->add_Matrix("tau3[OV][OV]{u}");
  blas->add_Matrix("tau3[oV][vO]{u}");
  blas->add_Matrix("tau3[Ov][Vo]{u}");

  blas->add_Matrix("Z_ijam[oov][o]{u}");
  blas->add_Matrix("Z_iJaM[oOv][O]{u}");
  blas->add_Matrix("Z_iJAm[oOV][o]{u}");
  blas->add_Matrix("Z_IJAM[OOV][O]{u}");

  // Mukherjee Terms
  blas->add_Matrix("Mk1[o][v]{u}");
  blas->add_Matrix("Mk1[O][V]{u}");
  blas->add_Matrix("Mk2[oo][vv]{u}");
  blas->add_Matrix("Mk2[oO][vV]{u}");
  blas->add_Matrix("Mk2[OO][VV]{u}");

  // Triples
  if(triples_type >= ccsd_t){
    blas->add_Matrix("epsilon[o][o]{u}");
    blas->add_Matrix("epsilon[O][O]{u}");
    blas->add_Matrix("epsilon[v][v]{u}");
    blas->add_Matrix("epsilon[V][V]{u}");

    blas->add_Matrix("t2[Oo][Vv]{u}");
    blas->add_Matrix("F_ME[O][V]{u}");

    blas->add_Matrix("W_ijka[oo][ov]{u}");
    blas->add_Matrix("W_iJkA[oO][oV]{u}");
    blas->add_Matrix("W_IjKa[Oo][Ov]{u}");
    blas->add_Matrix("W_IJKA[OO][OV]{u}");

    blas->add_Matrix("W_aibc[v][ovv]{u}");
    blas->add_Matrix("W_aIbC[v][OvV]{u}");
    blas->add_Matrix("W_AiBc[V][oVv]{u}");
    blas->add_Matrix("W_AIBC[V][OVV]{u}");
  }

  if(triples_type > ccsd_t){  // TODO: ccsd_t should not require storage

    blas->add_Matrix("ERROR{u}");
    blas->add_Matrix("<[oo]:[ov]>");
    blas->add_Matrix("<[oo]|[ov]>");

    // Required by T2 * W_bcek
    blas->add_Matrix("t2[oov][v]{u}");
    blas->add_Matrix("t2[oOv][V]{u}");
    blas->add_Matrix("t2[OoV][v]{u}");
    blas->add_Matrix("t2[OOV][V]{u}");

/*
    blas->add_Matrix("t2_test[oo][vv]{u}");
    blas->add_Matrix("t2_test[oO][vV]{u}");
    blas->add_Matrix("t2_test[OO][VV]{u}");*/



    blas->add_Matrix("W'_abic[vvo][v]{u}");
    blas->add_Matrix("W'_aBIc[vVO][v]{u}");
    blas->add_Matrix("W'_AbiC[Vvo][V]{u}");
    blas->add_Matrix("W'_ABIC[VVO][V]{u}");

    blas->add_Matrix("W'_ajki[voo][o]{u}");
    blas->add_Matrix("W'_AjKi[VoO][o]{u}");
    blas->add_Matrix("W'_aJkI[vOo][O]{u}");
    blas->add_Matrix("W'_AJKI[VOO][O]{u}");

    blas->add_Matrix("W_kija[o][oov]{u}");
    blas->add_Matrix("W_kiJA[o][oOV]{u}");
    blas->add_Matrix("W_KIja[O][Oov]{u}");
    blas->add_Matrix("W_KIJA[O][OOV]{u}");

    blas->add_Matrix("W_aibc[v][ovv]{u}");
    blas->add_Matrix("W_aIbC[v][OvV]{u}");
    blas->add_Matrix("W_AiBc[V][oVv]{u}");
    blas->add_Matrix("W_AIBC[V][OVV]{u}");

    //
    blas->add_Matrix("DELTA_t1[o][v]");
    blas->add_Matrix("DELTA_t1[O][V]");
    blas->add_Matrix("DELTA_t2[oo][vv]");
    blas->add_Matrix("DELTA_t2[oO][vV]");
    blas->add_Matrix("DELTA_t2[OO][VV]");

    // Required by T2*W_mcjk
    blas->add_Matrix("t2[ovv][o]{u}");
    blas->add_Matrix("t2[OvV][o]{u}");
    blas->add_Matrix("t2[oVv][O]{u}");
    blas->add_Matrix("t2[OVV][O]{u}");


    blas->add_Matrix("t3[ooo][vvv]{u}");
    blas->add_Matrix("t3[ooO][vvV]{u}");
    blas->add_Matrix("t3[oOO][vVV]{u}");
    blas->add_Matrix("t3[OOO][VVV]{u}");

    blas->add_Matrix("t3_eqns[ooo][vvv]{u}");
    blas->add_Matrix("t3_eqns[ooO][vvV]{u}");
    blas->add_Matrix("t3_eqns[oOO][vVV]{u}");
    blas->add_Matrix("t3_eqns[OOO][VVV]{u}");
  }

  blas->add_Matrix("fock[ff]{u}");
  blas->add_Matrix("fock[FF]{u}");

  if(pert_cbs){
    blas->add_Matrix("t2_eqns[oo][vf]{u}");
    blas->add_Matrix("t2_eqns[oO][vF]{u}");
    blas->add_Matrix("t2_eqns[OO][VF]{u}");

    blas->add_Matrix("t2_eqns[oo][fv]{u}");
    blas->add_Matrix("t2_eqns[oO][fV]{u}");
    blas->add_Matrix("t2_eqns[OO][FV]{u}");

    blas->add_Matrix("t2_eqns[oo][ff]{u}");
    blas->add_Matrix("t2_eqns[oO][fF]{u}");
    blas->add_Matrix("t2_eqns[OO][FF]{u}");

    blas->add_Matrix("t2_1[oo][vv]{u}");
    blas->add_Matrix("t2_1[oO][vV]{u}");
    blas->add_Matrix("t2_1[OO][VV]{u}");

    blas->add_Matrix("t2_1[oo][vf]{u}");
    blas->add_Matrix("t2_1[oO][vF]{u}");
    blas->add_Matrix("t2_1[OO][VF]{u}");

    blas->add_Matrix("t2_1[oo][fv]{u}");
    blas->add_Matrix("t2_1[oO][fV]{u}");
    blas->add_Matrix("t2_1[OO][FV]{u}");

    blas->add_Matrix("t2_1[oo][ff]{u}");
    blas->add_Matrix("t2_1[oO][fF]{u}");
    blas->add_Matrix("t2_1[OO][FF]{u}");

    blas->add_Matrix("<[oo]:[fv]>");
    blas->add_Matrix("<[oo]|[fv]>");

    blas->add_Matrix("<[oo]:[vf]>");
    blas->add_Matrix("<[oo]|[vf]>");

    blas->add_Matrix("<[oo]:[ff]>");
    blas->add_Matrix("<[oo]|[ff]>");

    blas->add_Matrix("d2[oo][vf]{u}");
    blas->add_Matrix("d2[oO][vF]{u}");
    blas->add_Matrix("d2[OO][VF]{u}");

    blas->add_Matrix("d2[oo][fv]{u}");
    blas->add_Matrix("d2[oO][fV]{u}");
    blas->add_Matrix("d2[OO][FV]{u}");

    blas->add_Matrix("d2[oo][ff]{u}");
    blas->add_Matrix("d2[oO][fF]{u}");
    blas->add_Matrix("d2[OO][FF]{u}");


    // THIRD-ORDER
    // Ladder terms
    blas->add_Matrix("<[vf]:[vv]>");
    blas->add_Matrix("<[vf]|[vv]>");

    blas->add_Matrix("<[ff]:[vv]>");
    blas->add_Matrix("<[ff]|[vv]>");

    // Ring terms
    blas->add_Matrix("([fo]:[ov])");
    blas->add_Matrix("([fo]|[ov])");

    blas->add_Matrix("<[of]|[ov]>");

    // T1 Contributions
    blas->add_Matrix("<[v]:[ovf]>");
    blas->add_Matrix("<[v]|[ovf]>");
    blas->add_Matrix("<[v]|[ofv]>");

    blas->add_Matrix("<[o]:[foo]>");
    blas->add_Matrix("<[o]|[foo]>");

    blas->add_Matrix("<[v]:[off]>");
    blas->add_Matrix("<[v]|[off]>");

    blas->add_Matrix("t2_2[oo][vf]{u}");
    blas->add_Matrix("t2_2[oO][vF]{u}");
    blas->add_Matrix("t2_2[OO][VF]{u}");

    blas->add_Matrix("t2_2[oo][ff]{u}");
    blas->add_Matrix("t2_2[oO][fF]{u}");
    blas->add_Matrix("t2_2[OO][FF]{u}");

    if(pert_cbs_coupling){
      blas->add_Matrix("<[fv]|[vv]>");

      blas->add_Matrix("([vo]:[of])");
      blas->add_Matrix("<[ov]|[of]>");
      blas->add_Matrix("<[ov]|[of]>");
      blas->add_Matrix("([vo]|[of])");
      blas->add_Matrix("([vo]:[of])");

      blas->add_Matrix("t2_1[ov][of]{u}");
      blas->add_Matrix("t2_1[ov][OF]{u}");
      blas->add_Matrix("t2_1[oV][Of]{u}");
      blas->add_Matrix("t2_1[oF][Ov]{u}");
      blas->add_Matrix("t2_1[OV][OF]{u}");
      blas->add_Matrix("t2_1[of][OV]{u}");

      blas->add_Matrix("t2_1[o][ovf]{u}");
      blas->add_Matrix("t2_1[o][OvF]{u}");
      blas->add_Matrix("t2_1[o][OfV]{u}");
      blas->add_Matrix("t2_1[o][off]{u}");
      blas->add_Matrix("t2_1[o][OfF]{u}");
      blas->add_Matrix("t2_1[v][foo]{u}");
      blas->add_Matrix("t2_1[v][FoO]{u}");
    }
  }
}

}} /* End Namespaces */
