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
#include "index_iterator.h"
#include "mrccsd_t.h"
#include "special_matrices.h"

namespace psi{ namespace psimrcc{

double MRCCSD_T::compute_AB_ooO_contribution_to_Heff(int u_abs,int V_abs,int x_abs,int Y_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3)
{
  double value = 0.0;
  int    i_sym  = o->get_tuple_irrep(i_abs);
  int    j_sym  = o->get_tuple_irrep(j_abs);
  int    k_sym  = o->get_tuple_irrep(k_abs);

  int  ijk_sym  = i_sym xor j_sym xor k_sym;

  size_t i_rel  = o->get_tuple_rel_index(i_abs);

  int  x_sym    = v->get_tuple_irrep(x_abs);
  int  y_sym    = v->get_tuple_irrep(Y_abs);
  int ij_sym    = oo->get_tuple_irrep(i_abs,j_abs);
  int jk_sym    = oo->get_tuple_irrep(j_abs,k_abs);
  int uv_sym    = oo->get_tuple_irrep(u_abs,V_abs);
  int xy_sym    = vv->get_tuple_irrep(x_abs,Y_abs);

  size_t  x_rel = v->get_tuple_rel_index(x_abs);
  size_t  y_rel = v->get_tuple_rel_index(Y_abs);

  size_t ij_rel = oo->get_tuple_rel_index(i_abs,j_abs);
  size_t kj_rel = oo->get_tuple_rel_index(k_abs,j_abs);
  size_t xy_rel = vv->get_tuple_rel_index(x_abs,Y_abs);

  if((j_abs == u_abs) and (k_abs == V_abs)){
    CCIndexIterator  e("[v]",i_sym);
    for(e.first(); !e.end(); e.next()){
      int    e_sym  = v->get_tuple_irrep(e.ind_abs<0>());
      size_t e_abs  = e.ind_abs<0>();
      size_t e_rel  = v->get_tuple_rel_index(e_abs);
      if(uv_sym == xy_sym){
        value += T3->get(e_sym,e_rel,xy_rel) * F2_ov[mu][i_sym][i_rel][e_rel];
      }
    }
  }
  if(i_abs == u_abs){
    CCIndexIterator  e("[v]",ijk_sym xor xy_sym);
    for(e.first(); !e.end(); e.next()){
      int    e_sym  = v->get_tuple_irrep(e.ind_abs<0>());
      size_t e_rel  = v->get_tuple_rel_index(e.ind_abs<0>());
      int    ve_sym = ov->get_tuple_irrep(V_abs,e.ind_abs<0>());
      size_t ve_rel = ov->get_tuple_rel_index(V_abs,e.ind_abs<0>());
      if(jk_sym == ve_sym){
        value += T3->get(e_sym,e_rel,xy_rel) * W_OoOv[mu][jk_sym][kj_rel][ve_rel];
      }
    }
  }
  if(k_abs == V_abs){
    CCIndexIterator  e("[v]",ijk_sym xor xy_sym);
    for(e.first(); !e.end(); e.next()){
      int    e_sym  = v->get_tuple_irrep(e.ind_abs<0>());
      size_t e_rel  = v->get_tuple_rel_index(e.ind_abs<0>());
      int    ue_sym = ov->get_tuple_irrep(u_abs,e.ind_abs<0>());
      size_t ue_rel = ov->get_tuple_rel_index(u_abs,e.ind_abs<0>());
      if(ij_sym == ue_sym){
        value += 0.5 * T3->get(e_sym,e_rel,xy_rel) * W_ooov[mu][ij_sym][ij_rel][ue_rel];
      }
    }
  }
  if((j_abs == u_abs) and (k_abs == V_abs)){
    CCIndexIterator  ef("[vv]",ijk_sym xor x_sym);
    for(ef.first(); !ef.end(); ef.next()){
      int   ief_sym  = ovv->get_tuple_irrep(i_abs,ef.ind_abs<0>(),ef.ind_abs<1>());
      size_t fe_rel  = vv->get_tuple_rel_index(ef.ind_abs<1>(),ef.ind_abs<0>());
      size_t ief_rel = ovv->get_tuple_rel_index(i_abs,ef.ind_abs<0>(),ef.ind_abs<1>());

      if(y_sym == ief_sym){
        value -= T3->get(x_sym,x_rel,fe_rel) * W_VoVv[mu][y_sym][y_rel][ief_rel];
      }
    }
  }
  if((j_abs == u_abs) and (k_abs == V_abs)){
    CCIndexIterator  ef("[vv]",ijk_sym xor y_sym);
    for(ef.first(); !ef.end(); ef.next()){
      int      e_sym =   v->get_tuple_irrep(ef.ind_abs<0>());
      int    ief_sym = ovv->get_tuple_irrep(i_abs,ef.ind_abs<0>(),ef.ind_abs<1>());
      size_t   e_rel =   v->get_tuple_rel_index(ef.ind_abs<0>());
      size_t  fy_rel =  vv->get_tuple_rel_index(ef.ind_abs<1>(),Y_abs);
      size_t ief_rel = ovv->get_tuple_rel_index(i_abs,ef.ind_abs<0>(),ef.ind_abs<1>());

      if(x_sym == ief_sym){
        value -= 0.5 * T3->get(e_sym,e_rel,fy_rel) * W_vovv[mu][x_sym][x_rel][ief_rel];
      }
    }
  }
  return value;
}

double MRCCSD_T::compute_AB_oOO_contribution_to_Heff(int u_abs,int V_abs,int x_abs,int Y_abs,int i_abs,int j_abs,int k_abs,int mu,BlockMatrix* T3)
{
  double value = 0.0;
  int    i_sym  = o->get_tuple_irrep(i_abs);
  int    j_sym  = o->get_tuple_irrep(j_abs);
  int    k_sym  = o->get_tuple_irrep(k_abs);

  int  ijk_sym  = i_sym xor j_sym xor k_sym;

  size_t k_rel  = o->get_tuple_rel_index(k_abs);


  int  x_sym    = v->get_tuple_irrep(x_abs);
  int  y_sym    = v->get_tuple_irrep(Y_abs);
  int ij_sym    = oo->get_tuple_irrep(i_abs,j_abs);
  int jk_sym    = oo->get_tuple_irrep(j_abs,k_abs);
  int uv_sym    = oo->get_tuple_irrep(u_abs,V_abs);
  int xy_sym    = vv->get_tuple_irrep(x_abs,Y_abs);

  size_t  x_rel = v->get_tuple_rel_index(x_abs);
  size_t  y_rel = v->get_tuple_rel_index(Y_abs);

  size_t ij_rel = oo->get_tuple_rel_index(i_abs,j_abs);
  size_t jk_rel = oo->get_tuple_rel_index(j_abs,k_abs);

  if((i_abs == u_abs) and (j_abs == V_abs)){
    CCIndexIterator  e("[v]",k_sym);
    for(e.first(); !e.end(); e.next()){
      size_t  e_rel  = v->get_tuple_rel_index(e.ind_abs<0>());
      size_t ye_rel  = vv->get_tuple_rel_index(Y_abs,e.ind_abs<0>());
      if(uv_sym == xy_sym){
        value += T3->get(x_sym,x_rel,ye_rel) * F2_OV[mu][k_sym][k_rel][e_rel];
      }
    }
  }
  if(i_abs == u_abs){
    CCIndexIterator  e("[v]",ijk_sym xor xy_sym);
    for(e.first(); !e.end(); e.next()){
      int    ve_sym = ov->get_tuple_irrep(V_abs,e.ind_abs<0>());
      size_t ve_rel = ov->get_tuple_rel_index(V_abs,e.ind_abs<0>());
      size_t ye_rel  = vv->get_tuple_rel_index(Y_abs,e.ind_abs<0>());
      if(jk_sym == ve_sym){
        value -= 0.5 * T3->get(x_sym,x_rel,ye_rel) * W_OOOV[mu][jk_sym][jk_rel][ve_rel];
      }
    }
  }
  if(k_abs == V_abs){
    CCIndexIterator  e("[v]",ijk_sym xor xy_sym);
    for(e.first(); !e.end(); e.next()){
      int    ue_sym = ov->get_tuple_irrep(u_abs,e.ind_abs<0>());
      size_t ue_rel = ov->get_tuple_rel_index(u_abs,e.ind_abs<0>());
      size_t ye_rel = vv->get_tuple_rel_index(Y_abs,e.ind_abs<0>());
      if(ij_sym == ue_sym){
        value += T3->get(x_sym,x_rel,ye_rel) * W_oOoV[mu][ij_sym][ij_rel][ue_rel];
      }
    }
  }
  if((i_abs == u_abs) and (j_abs == V_abs)){
    CCIndexIterator  ef("[vv]",ijk_sym xor x_sym);
    for(ef.first(); !ef.end(); ef.next()){
      int   kef_sym  = ovv->get_tuple_irrep(k_abs,ef.ind_abs<0>(),ef.ind_abs<1>());
      size_t ef_rel  = vv->get_tuple_rel_index(ef.ind_abs<0>(),ef.ind_abs<1>());
      size_t kef_rel = ovv->get_tuple_rel_index(k_abs,ef.ind_abs<0>(),ef.ind_abs<1>());
      if(y_sym == kef_sym){
        value += 0.5 * T3->get(x_sym,x_rel,ef_rel) * W_VOVV[mu][y_sym][y_rel][kef_rel];
      }
    }
  }
  if((i_abs == u_abs) and (j_abs == V_abs)){
    CCIndexIterator  ef("[vv]",ijk_sym xor y_sym);
    for(ef.first(); !ef.end(); ef.next()){
      int      e_sym =   v->get_tuple_irrep(ef.ind_abs<0>());
      int   kef_sym  = ovv->get_tuple_irrep(k_abs,ef.ind_abs<0>(),ef.ind_abs<1>());
      size_t   e_rel =   v->get_tuple_rel_index(ef.ind_abs<0>());
      size_t  yf_rel =  vv->get_tuple_rel_index(Y_abs,ef.ind_abs<1>());
      size_t kef_rel = ovv->get_tuple_rel_index(k_abs,ef.ind_abs<0>(),ef.ind_abs<1>());

      if(x_sym == kef_sym){
        value += T3->get(e_sym,e_rel,yf_rel) * W_vOvV[mu][x_sym][x_rel][kef_rel];
      }
    }
  }
  return value;
}


}} /* End Namespaces */
