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

/**
 *  @file ccmrcc_pert_triples_form_matrices.cpp
 *  @ingroup (PSIMRCC)
 *  @brief Computes the (T) correction
*/

#include <cstdlib>

#include "psi4/liboptions/liboptions.h"
#include "psi4/libmoinfo/libmoinfo.h"

#include "blas.h"
#include "index_iterator.h"
#include "matrix.h"
#include "mrccsd_t.h"
#include "special_matrices.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;

// Creates the matrix T2[ij][a][b] with irreps ordered according to b
void MRCCSD_T::form_T2_ij_a_b(IndexMatrix* T2_ij_a_b,bool spin1,bool spin2,bool transpose)
{
  CCIndexIterator  ij("[oo]");

  // Copy the matrix elements
  for(int ref = 0; ref < nrefs; ++ref){
    int unique_ref = moinfo->get_ref_number(ref,AllRefs);

    double*** Tijab;

    if(unique_ref == ref){
      if(spin1 && spin2){
        Tijab = blas->get_MatTmp("t2[oo][vv]",ref,none)->get_matrix();
      }else if (!spin1 && !spin2){
        Tijab = blas->get_MatTmp("t2[OO][VV]",ref,none)->get_matrix();
      }else if (spin1 && !spin2){
        Tijab = blas->get_MatTmp("t2[oO][vV]",ref,none)->get_matrix();
      }
    }else{
      if(spin1 && spin2){
        Tijab = blas->get_MatTmp("t2[OO][VV]",unique_ref,none)->get_matrix();
      }else if (!spin1 && !spin2){
        Tijab = blas->get_MatTmp("t2[oo][vv]",unique_ref,none)->get_matrix();
      }else if (spin1 && !spin2){
        Tijab = blas->get_MatTmp("t2[Oo][Vv]",unique_ref,none)->get_matrix();
      }
    }

    for(ij.first(); !ij.end(); ij.next()){
      BlockMatrix* block_matrix = new BlockMatrix(nirreps,v->get_tuplespi(),v->get_tuplespi(),ij.sym());
      CCIndexIterator  ab("[vv]",ij.sym());
//      ab.reset();
//      ab.set_irrep(ij.sym());
      for(ab.first(); !ab.end(); ab.next()){
        int    a_sym = v->get_tuple_irrep(ab.ind_abs<0>());
        int    b_sym = v->get_tuple_irrep(ab.ind_abs<1>());
        size_t a_rel = v->get_tuple_rel_index(ab.ind_abs<0>());
        size_t b_rel = v->get_tuple_rel_index(ab.ind_abs<1>());
        if(!transpose)
          block_matrix->set(a_sym,a_rel,b_rel,Tijab[ij.sym()][ij.rel()][ab.rel()]);
        else
          block_matrix->set(b_sym,b_rel,a_rel,Tijab[ij.sym()][ij.rel()][ab.rel()]);
      }
      T2_ij_a_b->add_block_matrix(ij.abs(),ref,block_matrix);
    }
  }
}

void MRCCSD_T::form_T2_i_ab_j(IndexMatrix* T2_i_ab_j,bool spin1,bool spin2,bool transpose)
{
  CCIndexIterator  i("[o]");

  // Copy the matrix elements
  for(int ref = 0; ref < nrefs; ++ref){
    int unique_ref = moinfo->get_ref_number(ref,AllRefs);

    double*** Tijab;

    if(unique_ref == ref){
      if(spin1 && spin2){
        Tijab = blas->get_MatTmp("t2[oo][vv]",ref,none)->get_matrix();
      }else if (!spin1 && !spin2){
        Tijab = blas->get_MatTmp("t2[OO][VV]",ref,none)->get_matrix();
      }else if (spin1 && !spin2){
        Tijab = blas->get_MatTmp("t2[oO][vV]",ref,none)->get_matrix();
      }
    }else{
      if(spin1 && spin2){
        Tijab = blas->get_MatTmp("t2[OO][VV]",unique_ref,none)->get_matrix();
      }else if (!spin1 && !spin2){
        Tijab = blas->get_MatTmp("t2[oo][vv]",unique_ref,none)->get_matrix();
      }else if (spin1 && !spin2){
        Tijab = blas->get_MatTmp("t2[Oo][Vv]",unique_ref,none)->get_matrix();
      }
    }

    for(i.first(); !i.end(); i.next()){
      BlockMatrix* block_matrix = new BlockMatrix(nirreps,vv->get_tuplespi(),o->get_tuplespi(),i.sym());
//      abj.reset();
//      abj.set_irrep(i.sym());
      CCIndexIterator  abj("[vvo]",i.sym());
      for(abj.first(); !abj.end(); abj.next()){
        int ij_sym = oo->get_tuple_irrep(i.ind_abs<0>(),abj.ind_abs<2>());
        size_t ij_rel = oo->get_tuple_rel_index(i.ind_abs<0>(),abj.ind_abs<2>());
        size_t ji_rel = oo->get_tuple_rel_index(abj.ind_abs<2>(),i.ind_abs<0>());
        int ab_sym = vv->get_tuple_irrep(abj.ind_abs<0>(),abj.ind_abs<1>());
        size_t ab_rel = vv->get_tuple_rel_index(abj.ind_abs<0>(),abj.ind_abs<1>());
        size_t j_rel = o->get_tuple_rel_index(abj.ind_abs<2>());
        if(!transpose)
          block_matrix->set(ab_sym,ab_rel,j_rel,Tijab[ij_sym][ij_rel][ab_rel]);
        else
          block_matrix->set(ab_sym,ab_rel,j_rel,Tijab[ij_sym][ji_rel][ab_rel]);
      }
      T2_i_ab_j->add_block_matrix(i.abs(),ref,block_matrix);
    }
  }
}

void MRCCSD_T::form_V_k_bc_e(IndexMatrix* V_k_bc_e,double direct,double exchange)
{
  // Build the matrices
  // (v_k)_{bc,e} = <bc||ek> and (v_k)_{bc,e} = <bc|ek>
  // from the integrals <ek||bc>
  // <bc||ek> = <ek||bc>
  CCIndexIterator  k("[o]");

  double*** V = blas->get_MatTmp("<[vo]|[vv]>",none)->get_matrix();
  for(k.first(); !k.end(); k.next()){
    BlockMatrix* block_matrix = new BlockMatrix(nirreps,vv->get_tuplespi(),v->get_tuplespi(),k.sym());
//    ebc.reset();
//    ebc.set_irrep(k.sym());
    CCIndexIterator  ebc("[vvv]",k.sym());
    for(ebc.first(); !ebc.end(); ebc.next()){
      size_t e_rel   = v->get_tuple_rel_index(ebc.ind_abs<0>());

      int    ek_sym  = vo->get_tuple_irrep(ebc.ind_abs<0>(),k.ind_abs<0>());
      size_t ek_rel  = vo->get_tuple_rel_index(ebc.ind_abs<0>(),k.ind_abs<0>());

      int    bc_sym  = vv->get_tuple_irrep(ebc.ind_abs<1>(),ebc.ind_abs<2>());
      size_t bc_rel  = vv->get_tuple_rel_index(ebc.ind_abs<1>(),ebc.ind_abs<2>());
      size_t cb_rel  = vv->get_tuple_rel_index(ebc.ind_abs<2>(),ebc.ind_abs<1>());

      double value = direct * V[ek_sym][ek_rel][bc_rel] + exchange * V[ek_sym][ek_rel][cb_rel];
      block_matrix->set(bc_sym,bc_rel,e_rel,value);
    }
    V_k_bc_e->add_block_matrix(k.abs(),int(0),block_matrix);
  }
}

void MRCCSD_T::form_V_jk_c_m(IndexMatrix* V_jk_c_m,double direct,double exchange)
{
  CCIndexIterator  jk("[oo]");

  double*** V = blas->get_MatTmp("<[oo]|[ov]>",none)->get_matrix();

  for(jk.first(); !jk.end(); jk.next()){
    BlockMatrix* block_matrix = new BlockMatrix(nirreps,v->get_tuplespi(),o->get_tuplespi(),jk.sym());
//    mc.reset();
//    mc.set_irrep(jk.sym());
    CCIndexIterator  mc("[ov]",jk.sym());
    for(mc.first(); !mc.end(); mc.next()){
      size_t m_rel   = o->get_tuple_rel_index(mc.ind_abs<0>());
      int    c_sym   = v->get_tuple_irrep(mc.ind_abs<1>());
      size_t c_rel   = v->get_tuple_rel_index(mc.ind_abs<1>());
      size_t kj_rel  = oo->get_tuple_rel_index(jk.ind_abs<1>(),jk.ind_abs<0>());

      double value = direct * V[jk.sym()][jk.rel()][mc.rel()] + exchange * V[jk.sym()][kj_rel][mc.rel()];
      block_matrix->set(c_sym,c_rel,m_rel,value);
    }
    V_jk_c_m->add_block_matrix(jk.abs(),int(0),block_matrix);
  }
}

}}  /* End Namespaces */
