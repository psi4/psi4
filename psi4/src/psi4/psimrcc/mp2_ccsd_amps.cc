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

#include "blas.h"
#include "mp2_ccsd.h"

namespace psi{ namespace psimrcc{

void MP2_CCSD::synchronize_amps()
{


  blas->solve("t1[ov]{u}     = #12# t1[o][v]{u}");
  blas->solve("t1[OV]{u}     = #12# t1[O][V]{u}");

  blas->reduce_spaces("t1_ov[a][v]{u}","t1[o][v]{u}");
  blas->reduce_spaces("t1_OV[A][V]{u}","t1[O][V]{u}");

  blas->reduce_spaces("t1_ov[o][a]{u}","t1[o][v]{u}");
  blas->reduce_spaces("t1_OV[O][A]{u}","t1[O][V]{u}");

  blas->solve("t2[o][ovv]{u} = #1234# t2[oo][vv]{u}");
  blas->solve("t2[o][OvV]{u} = #1234# t2[oO][vV]{u}");
  blas->solve("t2[O][oVv]{u} = #2143# t2[oO][vV]{u}");
  blas->solve("t2[O][OVV]{u} = #1234# t2[OO][VV]{u}");

  blas->solve("t2[v][voo]{u} = #3412# t2[oo][vv]{u}");
  blas->solve("t2[v][VoO]{u} = #3412# t2[oO][vV]{u}");
  blas->solve("t2[V][vOo]{u} = #4321# t2[oO][vV]{u}");
  blas->solve("t2[V][VOO]{u} = #3412# t2[OO][VV]{u}");

  blas->solve("t2[ov][OV]{u} = #1324# t2[oO][vV]{u}");
  blas->solve("t2[ov][ov]{u} = #1324# t2[oo][vv]{u}");

  blas->reduce_spaces("t2_oovv[o][aaa]{u}","t2[o][ovv]{u}");
  blas->reduce_spaces("t2_OoVv[O][aAa]{u}","t2[O][oVv]{u}");
  blas->reduce_spaces("t2_oOvV[o][AaA]{u}","t2[o][OvV]{u}");
  blas->reduce_spaces("t2_OOVV[O][AAA]{u}","t2[O][OVV]{u}");

  blas->reduce_spaces("t2_oovv[oo][aa]{u}","t2[oo][vv]{u}");
  blas->reduce_spaces("t2_oOvV[oO][aA]{u}","t2[oO][vV]{u}");
  blas->reduce_spaces("t2_OOVV[OO][AA]{u}","t2[OO][VV]{u}");

  blas->reduce_spaces("t2_oovv[a][ovv]{u}","t2[o][ovv]{u}");
  blas->reduce_spaces("t2_oOvV[a][OvV]{u}","t2[o][OvV]{u}");
  blas->reduce_spaces("t2_vvoo[a][voo]{u}","t2[v][voo]{u}");
  blas->reduce_spaces("t2_vVoO[a][VoO]{u}","t2[v][VoO]{u}");
  blas->reduce_spaces("t2_OOVV[A][OVV]{u}","t2[O][OVV]{u}");
  blas->reduce_spaces("t2_OoVv[A][oVv]{u}","t2[O][oVv]{u}");
  blas->reduce_spaces("t2_OoVv[O][oAa]{u}","t2[O][oVv]{u}");
  blas->reduce_spaces("t2_VVOO[A][VOO]{u}","t2[V][VOO]{u}");
  blas->reduce_spaces("t2_VvOo[A][vOo]{u}","t2[V][vOo]{u}");

  blas->reduce_spaces("t2_oovv[aa][vv]{u}","t2[oo][vv]{u}");
  blas->reduce_spaces("t2_oOvV[aA][vV]{u}","t2[oO][vV]{u}");
  blas->reduce_spaces("t2_OOVV[AA][VV]{u}","t2[OO][VV]{u}");

  blas->reduce_spaces("t2_vvoo[v][aaa]{u}","t2[v][voo]{u}");

  blas->reduce_spaces("t2_VvOo[V][aAa]{u}","t2[V][vOo]{u}");
  blas->reduce_spaces("t2_VvOo[V][aAo]{u}","t2[V][vOo]{u}");

  blas->reduce_spaces("t2_vVoO[v][AaA]{u}","t2[v][VoO]{u}");
  blas->reduce_spaces("t2_vVoO[v][AoA]{u}","t2[v][VoO]{u}");
  blas->reduce_spaces("t2_VVOO[V][AAA]{u}","t2[V][VOO]{u}");


  blas->solve("t2_ovOV[oa][OV]{u} = #2413# t2_vVoO[a][VoO]{u}");

  blas->solve("t2_oVOv[oA][Ov]{u} = #2431# t2_VvOo[A][vOo]{u}");


  blas->reduce_spaces("t2_oovv[ao][av]{u}","t2[oo][vv]{u}");
  blas->solve("t2_ovov[aa][ov]{u} = #1324# t2_oovv[ao][av]{u}");

  blas->solve("t2_ovov[oa][ov]{u} = #2413# t2_vvoo[a][voo]{u}");

  blas->solve("t2_ovov[av][ov]{u} = #1324# t2_oovv[a][ovv]{u}");

  blas->reduce_spaces("t2_oOvV[oA][vA]{u}","t2[oO][vV]{u}");
  blas->solve("t2_ovOV[ov][AA]{u} = #1324# t2_oOvV[oA][vA]{u}");

  blas->reduce_spaces("t2_oOvV[aO][aV]{u}","t2[oO][vV]{u}");
  blas->solve("t2_ovOV[aa][OV]{u} = #1324# t2_oOvV[aO][aV]{u}");

  blas->solve("t2_ovOV[av][OV]{u} = #1324# t2_oOvV[a][OvV]{u}");

//   blas->reduce_spaces("t2_oOvV[aO][vV]{u}","t2[oO][vV]{u}");
//   blas->solve("t2_ovOV[av][OV]{u} = #1324# t2_oOvV[aO][vV]{u}");

  blas->reduce_spaces("t2_OOVV[AO][AV]{u}","t2[OO][VV]{u}");
  blas->solve("t2_OVOV[AA][OV]{u} = #1324# t2_OOVV[AO][AV]{u}");

  blas->reduce_spaces("t2_oOvV[aO][vA]{u}","t2[oO][vV]{u}");
  blas->solve("t2_oVOv[aA][Ov]{u} = #1342# t2_oOvV[aO][vA]{u}");

  blas->reduce_spaces("t2_oOvV[oA][aV]{u}","t2[oO][vV]{u}");
  blas->solve("t2_oVOv[oV][Aa]{u} = #1342# t2_oOvV[oA][aV]{u}");

  blas->solve("t2_oVOv[oV][Av]{u} = #3124# t2_OoVv[A][oVv]{u}");
  

  blas->solve("t2_VvOo[V][vAa]{u} = #4321# t2_oOvV[aA][vV]{u}");

  blas->solve("t2_OoVv[O][aAv]{u} = #2143# t2_oOvV[aO][vA]{u}");

  blas->solve("t2_oOvV[o][AvA]{u} = #1234# t2_oOvV[oA][vA]{u}");

}

void MP2_CCSD::build_tau()
{
//   // t1t1[ov][ov]{u}, Ok
//   blas->solve("t1t1_iame[ov][ov]{u} = #1432#   t1[o][v]{u} X t1[o][v]{u}");  
//   blas->solve("t1t1_IAME[OV][OV]{u} = #1432#   t1[O][V]{u} X t1[O][V]{u}");  
//   blas->solve("t1t1_iAMe[oV][Ov]{u} = #1432#   t1[o][v]{u} X t1[O][V]{u}");  

  // tau[oo][vv]{u}, Ok
  blas->solve("tau[oo][vv]{u}  = t2[oo][vv]{u}");
  blas->solve("tau[oo][vv]{u} += #1324#   t1[o][v]{u} X t1[o][v]{u}");  
  blas->solve("tau[oo][vv]{u} += #2314# - t1[o][v]{u} X t1[o][v]{u}"); 
  // tau[oO][vV]{u}, Ok
  blas->solve("tau[oO][vV]{u}  = t2[oO][vV]{u}");
  blas->solve("tau[oO][vV]{u} += #1324#   t1[o][v]{u} X t1[O][V]{u}");  
  // tau[OO][VV]{u}, Ok
  blas->solve("tau[OO][VV]{u}  = t2[OO][VV]{u}");
  blas->solve("tau[OO][VV]{u} += #1324#   t1[O][V]{u} X t1[O][V]{u}");  
  blas->solve("tau[OO][VV]{u} += #2314# - t1[O][V]{u} X t1[O][V]{u}"); 

//   // tau[oO][Vv]{u}, Ok
//   blas->solve("tau[oO][Vv]{u}  = #1243#   tau[oO][vV]{u}");  

  // tau2[v][voo]{u}, Ok
  blas->solve("tau2[v][voo]{u}  = #3412# t2[oo][vv]{u}");
  blas->solve("tau2[v][voo]{u} += #3142# 1/2 t1[o][v]{u} X t1[o][v]{u}");  
  blas->solve("tau2[v][voo]{u} += #4132# -1/2 t1[o][v]{u} X t1[o][v]{u}"); 

  // tau2[V][VOO]{u}, Ok
  blas->solve("tau2[V][VOO]{u}  = #3412# t2[OO][VV]{u}");
  blas->solve("tau2[V][VOO]{u} += #3142# 1/2 t1[O][V]{u} X t1[O][V]{u}");  
  blas->solve("tau2[V][VOO]{u} += #4132# -1/2 t1[O][V]{u} X t1[O][V]{u}");

  // tau2[v][VoO]{u}, Ok
  blas->solve("tau2[v][VoO]{u}  = #3412# t2[oO][vV]{u}");
  blas->solve("tau2[v][VoO]{u} += #3142# 1/2 t1[o][v]{u} X t1[O][V]{u}");

  // tau2[V][vOo]{u}, Ok
  blas->solve("tau2[V][vOo]{u}  = #4321# t2[oO][vV]{u}");
  blas->solve("tau2[V][vOo]{u} += #4231# 1/2 t1[o][v]{u} X t1[O][V]{u}");  

  // tau2[o][ovv]{u}, Ok
  blas->solve("tau2[o][ovv]{u}  = #1234# t2[oo][vv]{u}");
  blas->solve("tau2[o][ovv]{u} += #1324# 1/2 t1[o][v]{u} X t1[o][v]{u}");  
  blas->solve("tau2[o][ovv]{u} += #2314# -1/2 t1[o][v]{u} X t1[o][v]{u}"); 

  // tau2[O][OVV]{u}, Ok
  blas->solve("tau2[O][OVV]{u}  = #1234# t2[OO][VV]{u}");
  blas->solve("tau2[O][OVV]{u} += #1324# 1/2 t1[O][V]{u} X t1[O][V]{u}");  
  blas->solve("tau2[O][OVV]{u} += #2314# -1/2 t1[O][V]{u} X t1[O][V]{u}"); 

  // tau2[o][OvV]{u}, Ok
  blas->solve("tau2[o][OvV]{u}  = #1234# t2[oO][vV]{u}");
  blas->solve("tau2[o][OvV]{u} += #1324# 1/2 t1[o][v]{u} X t1[O][V]{u}");  

  // tau2[O][oVv]{u}, Ok
  blas->solve("tau2[O][oVv]{u}  = #2143# t2[oO][vV]{u}");
  blas->solve("tau2[O][oVv]{u} += #2413# 1/2 t1[o][v]{u} X t1[O][V]{u}");  


  blas->reduce_spaces("tau_oOvV[oO][aA]{u}","tau[oO][vV]{u}");
  blas->reduce_spaces("tau_oOvV[oA][vV]{u}","tau[oO][vV]{u}");
  blas->reduce_spaces("tau_oOvV[oO][vA]{u}","tau[oO][vV]{u}");


  blas->reduce_spaces("tau_oOvV[aA][vV]{u}","tau[oO][vV]{u}");

  blas->solve("tau_oOVv[aA][Vv]{u} = #1243# tau_oOvV[aA][vV]{u}");  
  blas->solve("tau_oOVv[oA][Vv]{u} = #1243# tau_oOvV[oA][vV]{u}");



  blas->solve("t1t1_iame[aa][ov]{u} = #1432#   t1_ov[a][v]{u} X t1_ov[o][a]{u}");
  blas->solve("t1t1_iame[oa][ov]{u} = #1432#   t1[o][v]{u} X t1_ov[o][a]{u}");

  blas->solve("t1t1_iame[av][ov]{u} = #1432#   t1_ov[a][v]{u} X t1[o][v]{u}");

  blas->solve("t1t1_IAME[AA][OV]{u} = #1432#   t1_OV[A][V]{u} X t1_OV[O][A]{u}");
  blas->solve("t1t1_iAMe[aA][Ov]{u} = #1432#   t1_ov[a][v]{u} X t1_OV[O][A]{u}");
  blas->solve("t1t1_iAMe[oA][Ov]{u} = #1432#   t1[o][v]{u} X t1_OV[O][A]{u}");
  blas->solve("t1t1_iAMe[oV][Aa]{u} = #1432#   t1_ov[o][a]{u} X t1_OV[A][V]{u}");

  blas->solve("t1t1_iAMe[oV][Av]{u} = #1432#   t1[o][v]{u} X t1_OV[A][V]{u}");

  // tau3[ov][ov]{u}
  blas->reduce_spaces("t2_oovv[ao][va]{u}","t2[oo][vv]{u}");
  blas->solve("tau3_ovov[aa][ov]{u}  = #1342# 1/2 t2_oovv[ao][va]{u}");
  blas->solve("tau3_ovov[aa][ov]{u} += #1432# t1_ov[a][v]{u} X t1_ov[o][a]{u}");

  blas->reduce_spaces("t2_oovv[oo][va]{u}","t2[oo][vv]{u}");
  blas->solve("tau3_ovov[oa][ov]{u}  = #1342# 1/2 t2_oovv[oo][va]{u}");
  blas->solve("tau3_ovov[oa][ov]{u} += #1432# t1[o][v]{u} X t1_ov[o][a]{u}");


  blas->reduce_spaces("t2_oovv[ao][vv]{u}","t2[oo][vv]{u}");
  blas->solve("tau3_ovov[av][ov]{u}  = #1342# 1/2 t2_oovv[ao][vv]{u}");
  blas->solve("tau3_ovov[av][ov]{u} += #1432# t1_ov[a][v]{u} X t1[o][v]{u}");

  // tau3[oV][vO]{u}
  blas->reduce_spaces("t2_oOvV[aO][vA]{u}","t2[oO][vV]{u}");
  blas->solve("tau3_oVvO[aA][vO]{u}  = #1432# 1/2 t2_oOvV[aO][vA]{u}");
  blas->solve("tau3_oVvO[aA][vO]{u} += #1342# t1_ov[a][v]{u} X t1_OV[O][A]{u}");

  blas->reduce_spaces("t2_oOvV[oO][vA]{u}","t2[oO][vV]{u}");
  blas->solve("tau3_oVvO[oA][vO]{u}  = #1432# 1/2 t2_oOvV[oO][vA]{u}");
  blas->solve("tau3_oVvO[oA][vO]{u} += #1342# t1[o][v]{u} X t1_OV[O][A]{u}");

  blas->reduce_spaces("t2_oOvV[aO][vV]{u}","t2[oO][vV]{u}");
  blas->solve("tau3_oVvO[aV][vO]{u}  = #1432# 1/2 t2_oOvV[aO][vV]{u}");
  blas->solve("tau3_oVvO[aV][vO]{u} += #1342# t1_ov[a][v]{u} X t1[O][V]{u}");

/*


  // tau3[Ov][Vo]{u}
  blas->solve("tau3[Ov][Vo]{u}  = #4123# 1/2 t2[oO][vV]{u}");
  blas->solve("tau3[Ov][Vo]{u} += #4213# t1[o][v]{u} X t1[O][V]{u}");*/
}

}}
