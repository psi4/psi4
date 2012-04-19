/**
 *  @file ccmrcc_pert_triples.cpp
 *  @ingroup (PSIMRCC)
 *  @brief Computes the (T) correction
*/

#include <cstdlib>

#include <liboptions/liboptions.h>
#include <libmoinfo/libmoinfo.h>
#include <libchkpt/chkpt.hpp>

#include "blas.h"
#include "debugging.h"
#include "index_iterator.h"
#include "mrcc.h"
#include "mrccsd_t.h"
#include "special_matrices.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

void MRCCSD_T::compute_spin_adapted()
{
  fprintf(outfile,"\n\n  Computing (T) correction using the spin-adapted algorithm.\n");
  fflush(outfile);

  compute_ooO_triples_spin_adapted();

  fprintf(outfile,"\n\n  Mk-MRCCSD(T) diagonal contributions to the effective Hamiltonian:\n");
  fprintf(outfile,"\n   Ref         E[4]              E_T[4]            E_ST[4]           E_DT[4]");
  fprintf(outfile,"\n  ------------------------------------------------------------------------------");
  if(nrefs < 100){
    for(int mu = 0; mu < nrefs; ++mu){
      fprintf(outfile,"\n   %2d  ",mu);
      fprintf(outfile," %17.12lf",E4_ooo[mu] + E4_ooO[mu] + E4_oOO[mu] + E4_OOO[mu]);
      fprintf(outfile," %17.12lf",E4T_ooo[mu] + E4T_ooO[mu] + E4T_oOO[mu] + E4T_OOO[mu]);
      fprintf(outfile," %17.12lf",E4ST_ooo[mu] + E4ST_ooO[mu] + E4ST_oOO[mu] + E4ST_OOO[mu]);
      fprintf(outfile," %17.12lf",E4DT_ooo[mu] + E4DT_ooO[mu] + E4DT_oOO[mu] + E4DT_OOO[mu]);
    }
  }
  fprintf(outfile,"\n   Tot ");
  double E4 = 0.0;
  for(int mu = 0; mu < nrefs; ++mu)
    E4 += (E4_ooo[mu] + E4_ooO[mu] + E4_oOO[mu] + E4_OOO[mu]) * h_eff->get_left_eigenvector(mu) * h_eff->get_right_eigenvector(mu);
  fprintf(outfile," %17.12lf",E4);
  double E4T = 0.0;
  for(int mu = 0; mu < nrefs; ++mu)
    E4T += (E4T_ooo[mu] + E4T_ooO[mu] + E4T_oOO[mu] + E4T_OOO[mu])  * h_eff->get_left_eigenvector(mu) * h_eff->get_right_eigenvector(mu);
  fprintf(outfile," %17.12lf",E4T);
  double E4ST = 0.0;
  for(int mu = 0; mu < nrefs; ++mu)
    E4ST += (E4ST_ooo[mu] + E4ST_ooO[mu] + E4ST_oOO[mu] + E4ST_OOO[mu])  * h_eff->get_left_eigenvector(mu) * h_eff->get_right_eigenvector(mu);
  fprintf(outfile," %17.12lf",E4ST);
  double E4DT = 0.0;
  for(int mu = 0; mu < nrefs; ++mu)
    E4DT += (E4DT_ooo[mu] + E4DT_ooO[mu] + E4DT_oOO[mu] + E4DT_OOO[mu]) * h_eff->get_left_eigenvector(mu) * h_eff->get_right_eigenvector(mu);
  fprintf(outfile," %17.12lf",E4DT);
  fprintf(outfile,"\n  ------------------------------------------------------------------------------");

  fprintf(outfile,"\n\n  Mk-MRCCSD(T) off-diagonal contributions to the effective Hamiltonian:\n");
  for(int mu = 0; mu < nrefs; ++mu){
    fprintf(outfile,"\n");
    for(int nu = 0; nu < nrefs; ++nu){
      fprintf(outfile," %17.12lf",d_h_eff[mu][nu]);
    }
  }

  for(int mu = 0; mu < nrefs; ++mu){
    for(int nu = 0; nu < nrefs; ++nu){
      if(mu != nu){
        if(options_.get_bool("OFFDIAGONAL_CCSD_T")){  // Option to add the diagonal correction
          h_eff->add_matrix(mu,nu,2.0 * d_h_eff[mu][nu]);
        }
      }else{
        if(options_.get_bool("DIAGONAL_CCSD_T")){  // Option to add the off-diagonal correction
          h_eff->add_matrix(mu,nu,2.0 * E4_ooo[mu] + 2.0 * E4_ooO[mu]);
        }
      }
    }
  }
  h_eff->print_matrix();
}

void MRCCSD_T::compute_ooO_triples_spin_adapted()
{
  CCIndexIterator  ijk(ooo);

  for(ijk.first(); !ijk.end(); ijk.next()){

    size_t i_abs = o->get_tuple_abs_index(ijk.ind_abs<0>());
    size_t j_abs = o->get_tuple_abs_index(ijk.ind_abs<1>());
    size_t k_abs = o->get_tuple_abs_index(ijk.ind_abs<2>());

    if((i_abs <= j_abs) and (j_abs <= k_abs)){
      int i_sym     = o->get_tuple_irrep(ijk.ind_abs<0>());
      int j_sym     = o->get_tuple_irrep(ijk.ind_abs<1>());
      int k_sym     = o->get_tuple_irrep(ijk.ind_abs<2>());

      size_t i_rel = o->get_tuple_rel_index(ijk.ind_abs<0>());
      size_t j_rel = o->get_tuple_rel_index(ijk.ind_abs<1>());
      size_t k_rel = o->get_tuple_rel_index(ijk.ind_abs<2>());

      size_t ij_abs = oo->get_tuple_abs_index(ijk.ind_abs<0>(),ijk.ind_abs<1>());
      size_t ji_abs = oo->get_tuple_abs_index(ijk.ind_abs<1>(),ijk.ind_abs<0>());
      size_t ik_abs = oo->get_tuple_abs_index(ijk.ind_abs<0>(),ijk.ind_abs<2>());
      size_t ki_abs = oo->get_tuple_abs_index(ijk.ind_abs<2>(),ijk.ind_abs<0>());
      size_t jk_abs = oo->get_tuple_abs_index(ijk.ind_abs<1>(),ijk.ind_abs<2>());
      size_t kj_abs = oo->get_tuple_abs_index(ijk.ind_abs<2>(),ijk.ind_abs<1>());

      int    ij_sym = oo->get_tuple_irrep(ijk.ind_abs<0>(),ijk.ind_abs<1>());
      int    ik_sym = oo->get_tuple_irrep(ijk.ind_abs<0>(),ijk.ind_abs<2>());
      int    jk_sym = oo->get_tuple_irrep(ijk.ind_abs<1>(),ijk.ind_abs<2>());
      size_t ij_rel = oo->get_tuple_rel_index(ijk.ind_abs<0>(),ijk.ind_abs<1>());
      size_t ji_rel = oo->get_tuple_rel_index(ijk.ind_abs<1>(),ijk.ind_abs<0>());
      size_t ik_rel = oo->get_tuple_rel_index(ijk.ind_abs<0>(),ijk.ind_abs<2>());
      size_t ki_rel = oo->get_tuple_rel_index(ijk.ind_abs<2>(),ijk.ind_abs<0>());
      size_t jk_rel = oo->get_tuple_rel_index(ijk.ind_abs<1>(),ijk.ind_abs<2>());
      size_t kj_rel = oo->get_tuple_rel_index(ijk.ind_abs<2>(),ijk.ind_abs<1>());

      int ijk_sym = ijk.sym();

      // Compute W for all unique references (d N^7)
      for(int mu = 0; mu < nrefs; ++mu){
        // Check if ijk belong to the occupied space of mu
        if(is_aocc[mu][i_abs] && is_aocc[mu][j_abs] && is_aocc[mu][k_abs]){
            // 1. x_ijk^abc
            Z[mu][ijk_sym]->contract(T2_iJ_a_B->get_block_matrix(ij_abs,mu),V_K_bC_e->get_block_matrix(k_abs),1.0,0.0);
            // + x_ijk^abc - x_ijk^bac                                           abc  acb  bac  bca  cab  cba
            W_ijk[mu][ijk_sym]->add(Z[mu][ijk_sym],0.0,1.0);
            // + x_ijk^acb - x_ijk^bca                                           abc  acb  bac  bca  cab  cba
            W_ikj[mu][ijk_sym]->add_acb(0.0,Z[mu][ijk_sym],vvv,v,vv,1.0);
            // + x_ijk^cab - x_ijk^cba                                           abc  acb  bac  bca  cab  cba
            W_jki[mu][ijk_sym]->add_cab(0.0,Z[mu][ijk_sym],vvv,v,vv,1.0);//->add_permutation_1_2(0.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

            // 2. x_jik^abc
            Z[mu][ijk_sym]->contract(T2_iJ_a_B->get_block_matrix(ji_abs,mu),V_K_bC_e->get_block_matrix(k_abs),1.0,0.0);
            // - x_jik^abc + x_jik^bac                                           abc  acb  bac  bca  cab  cba
            W_ijk[mu][ijk_sym]->add(Z[mu][ijk_sym],1.0,-1.0);
            // + x_jik^cab - x_jik^cba                                           abc  acb  bac  bca  cab  cba
            W_ikj[mu][ijk_sym]->add_cab(1.0,Z[mu][ijk_sym],vvv,v,vv,1.0);
            // + x_jik^acb - x_jik^bca                                           abc  acb  bac  bca  cab  cba
            W_jki[mu][ijk_sym]->add_acb(1.0,Z[mu][ijk_sym],vvv,v,vv,1.0);//->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0);

            // 3. x_jki^abc
            Z[mu][ijk_sym]->contract(T2_iJ_a_B->get_block_matrix(jk_abs,mu),V_K_bC_e->get_block_matrix(i_abs),1.0,0.0);
            // - x_jki^acb + x_jki^bca                                           abc  acb  bac  bca  cab  cba
            W_ijk[mu][ijk_sym]->add_acb(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0);//->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0);
            // - x_jki^cab + x_jki^cba                                           abc  acb  bac  bca  cab  cba
            W_ikj[mu][ijk_sym]->add_cab(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0);
            // + x_jki^abc - x_jki^bac                                           abc  acb  bac  bca  cab  cba
            W_jki[mu][ijk_sym]->add(Z[mu][ijk_sym],1.0,1.0);

            // 4. x_kji^abc
            Z[mu][ijk_sym]->contract(T2_iJ_a_B->get_block_matrix(kj_abs,mu),V_K_bC_e->get_block_matrix(i_abs),1.0,0.0);
            // - x_kji^cab + x_kji^cba                                           abc  acb  bac  bca  cab  cba
            W_ijk[mu][ijk_sym]->add_cab(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0);
            // - x_kji^acb + x_kji^bac                                           abc  acb  bac  bca  cab  cba
            W_ikj[mu][ijk_sym]->add_acb(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0);//->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0);
            // - x_kji^abc + x_kji^bac                                           abc  acb  bac  bca  cab  cba
            W_jki[mu][ijk_sym]->add(Z[mu][ijk_sym],1.0,-1.0);

            // 5. x_kij^abc
            Z[mu][ijk_sym]->contract(T2_iJ_a_B->get_block_matrix(ki_abs,mu),V_K_bC_e->get_block_matrix(j_abs),1.0,0.0);
            // + x_kij^cab - x_kij^cba                                           abc  acb  bac  bca  cab  cba
            W_ijk[mu][ijk_sym]->add_cab(1.0,Z[mu][ijk_sym],vvv,v,vv,1.0);
            // - x_kij^abc + x_kij^bac                                           abc  acb  bac  bca  cab  cba
            W_ikj[mu][ijk_sym]->add(Z[mu][ijk_sym],1.0,-1.0);
            // - x_kij^acb + x_kij^bca                                           abc  acb  bac  bca  cab  cba
            W_jki[mu][ijk_sym]->add_acb(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0);//->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0);

            // 6. x_ikj^abc
            Z[mu][ijk_sym]->contract(T2_iJ_a_B->get_block_matrix(ik_abs,mu),V_K_bC_e->get_block_matrix(j_abs),1.0,0.0);
            // + x_ikj^acb - x_ikj^bca                                           abc  acb  bac  bca  cab  cba
            W_ijk[mu][ijk_sym]->add_acb(1.0,Z[mu][ijk_sym],vvv,v,vv,1.0);//->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
            // + x_ikj^abc - x_ikj^bac                                           abc  acb  bac  bca  cab  cba
            W_ikj[mu][ijk_sym]->add(Z[mu][ijk_sym],1.0,1.0);
            // - x_ikj^cab + x_ikj^cba                                           abc  acb  bac  bca  cab  cba
            W_jki[mu][ijk_sym]->add_cab(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0);

            // 7. y_ijk^abc
            Z[mu][ijk_sym]->contract(V_jK_C_m->get_block_matrix(ij_abs),T2_i_aB_J->get_block_matrix(k_abs,mu),1.0,0.0);
            // + y_ijk^acb - y_ijk^bca                                           abc  acb  bac  bca  cab  cba
            W_ijk[mu][ijk_sym]->add_acb(1.0,Z[mu][ijk_sym],vvv,v,vv,1.0);//->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
            // + y_ijk^cab - y_ijk^cba                                           abc  acb  bac  bca  cab  cba
            W_ikj[mu][ijk_sym]->add_cab(1.0,Z[mu][ijk_sym],vvv,v,vv,1.0);
            // - y_ijk^abc + y_ijk^bac                                           abc  acb  bac  bca  cab  cba
            W_jki[mu][ijk_sym]->add(Z[mu][ijk_sym],1.0,-1.0);

            // 8. y_jik^acb
            Z[mu][ijk_sym]->contract(V_jK_C_m->get_block_matrix(ji_abs),T2_i_aB_J->get_block_matrix(k_abs,mu),1.0,0.0);
            // - y_jik^acb + y_jik^bca                                           abc  acb  bac  bca  cab  cba
            W_ijk[mu][ijk_sym]->add_acb(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0);//->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0);
            // - y_jik^abc + y_jik^bac                                           abc  acb  bac  bca  cab  cba
            W_ikj[mu][ijk_sym]->add(Z[mu][ijk_sym],1.0,-1.0);
            // + y_jik^cab - y_jik^cba                                           abc  acb  bac  bca  cab  cba
            W_jki[mu][ijk_sym]->add_cab(1.0,Z[mu][ijk_sym],vvv,v,vv,1.0);

            // 9. y_jki^abc
            Z[mu][ijk_sym]->contract(V_jK_C_m->get_block_matrix(jk_abs),T2_i_aB_J->get_block_matrix(i_abs,mu),1.0,0.0);
            // - y_jki^cab + y_jki^cba                                           abc  acb  bac  bca  cab  cba
            W_ijk[mu][ijk_sym]->add_cab(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0);
            // + y_jki^abc - y_jki^bac                                           abc  acb  bac  bca  cab  cba
            W_ikj[mu][ijk_sym]->add(Z[mu][ijk_sym],1.0,1.0);
            // + y_jki^acb - y_jki^bca                                           abc  acb  bac  bca  cab  cba
            W_jki[mu][ijk_sym]->add_acb(1.0,Z[mu][ijk_sym],vvv,v,vv,1.0);//->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0);

            // 10. y_kji^abc
            Z[mu][ijk_sym]->contract(V_jK_C_m->get_block_matrix(kj_abs),T2_i_aB_J->get_block_matrix(i_abs,mu),1.0,0.0);
            // + y_kji^abc - y_kji^bac                                           abc  acb  bac  bca  cab  cba
            W_ijk[mu][ijk_sym]->add(Z[mu][ijk_sym],1.0,1.0);
            // - y_kji^cab + y_kji^cba                                           abc  acb  bac  bca  cab  cba
            W_ikj[mu][ijk_sym]->add_cab(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0);
            // - y_kji^acb + y_kji^bca                                           abc  acb  bac  bca  cab  cba
            W_jki[mu][ijk_sym]->add_acb(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0);//->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0);

            // 11. y_kij^abc
            Z[mu][ijk_sym]->contract(V_jK_C_m->get_block_matrix(ki_abs),T2_i_aB_J->get_block_matrix(j_abs,mu),1.0,0.0);
            // - y_kij^abc + y_kij^bac                                           abc  acb  bac  bca  cab  cba
            W_ijk[mu][ijk_sym]->add(Z[mu][ijk_sym],1.0,-1.0);
            // - y_kij^acb + y_kij^bca                                           abc  acb  bac  bca  cab  cba
            W_ikj[mu][ijk_sym]->add_acb(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0);//->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0);
            // - y_kij^cab + y_kij^cba                                           abc  acb  bac  bca  cab  cba
            W_jki[mu][ijk_sym]->add_cab(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0);

            // 12. y_ikj^abc
            Z[mu][ijk_sym]->contract(V_jK_C_m->get_block_matrix(ik_abs),T2_i_aB_J->get_block_matrix(j_abs,mu),1.0,0.0);
            // + y_ikj^cab - y_ikj^cba                                           abc  acb  bac  bca  cab  cba
            W_ijk[mu][ijk_sym]->add_cab(1.0,Z[mu][ijk_sym],vvv,v,vv,1.0);
            // + y_ikj^acb - y_ikj^bca                                           abc  acb  bac  bca  cab  cba
            W_ikj[mu][ijk_sym]->add_acb(1.0,Z[mu][ijk_sym],vvv,v,vv,1.0);//->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
            // + y_ikj^abc - y_ikj^bac                                           abc  acb  bac  bca  cab  cba
            W_jki[mu][ijk_sym]->add(Z[mu][ijk_sym],1.0,1.0);

            W_ijk[mu][ijk_sym]->a_b_permutation(vvv,v,vv);
            W_ikj[mu][ijk_sym]->a_b_permutation(vvv,v,vv);
            W_jki[mu][ijk_sym]->a_b_permutation(vvv,v,vv);
        }
      }

      //////////////
      // IJK CASE //
      //////////////
      if((i_abs < j_abs) and (j_abs <= k_abs)){
        for(int mu = 0; mu < nrefs; ++mu){
          T[mu][ijk_sym]->zero();
        }

        // Compute T (d^2 N^6)
        int    cycle = 0;
        double oldE  = 1.0;
        double newE  = 0.0;
        while(fabs(oldE-newE) > threshold){
          cycle++;
          oldE = newE;
          newE = 0.0;
          // Iterate the Mk-MRCCSD(T) Equations
          for(int mu = 0; mu < nrefs; ++mu){
            e4T[mu] = e4ST[mu] = e4DT[mu] = 0.0;
            // Check if ijk belong to the occupied space of mu
            if(is_aocc[mu][i_abs] && is_aocc[mu][j_abs] && is_bocc[mu][k_abs]){

              double***  F_ov_mu    = F_ov[mu];
              double***  F_OV_mu    = F_OV[mu];
              double***  T1_ov_mu   = T1_ov[mu];
              double***  T1_OV_mu   = T1_OV[mu];
              double***  T2_oovv_mu = T2_oovv[mu];
              double***  T2_oOvV_mu = T2_oOvV[mu];

              double D_ijK = e_oo[mu][i_abs] + e_oo[mu][j_abs] + e_OO[mu][k_abs];

              // Add W
              Z[mu][ijk_sym]->add(W_ijk[mu][ijk_sym],0.0,1.0);

              // Add the coupling terms
              for(int nu = 0; nu < nrefs; ++nu){
                if(nu != mu){
                  Z[mu][ijk_sym]->add(T[nu][ijk_sym],1.0,Mk_factor[mu][nu]);
                }
              }

              // Divide by the denominator
              std::vector<double>& e_vv_mu  = e_vv[mu];
              std::vector<double>& e_VV_mu  = e_VV[mu];
              std::vector<bool>& is_avir_mu = is_avir[mu];
              std::vector<bool>& is_bvir_mu = is_bvir[mu];

              CCIndexIterator  abc(vvv,ijk_sym);
              // Loop over abc
              for(abc.first(); !abc.end(); abc.next()){
                size_t a_abs = v->get_tuple_abs_index(abc.ind_abs<0>());
                size_t b_abs = v->get_tuple_abs_index(abc.ind_abs<1>());
                size_t c_abs = v->get_tuple_abs_index(abc.ind_abs<2>());
                if(is_avir_mu[a_abs] && is_avir_mu[b_abs] && is_bvir_mu[c_abs]){
                  int     a_sym = v->get_tuple_irrep(abc.ind_abs<0>());
                  int     c_sym = v->get_tuple_irrep(abc.ind_abs<2>());
                  int    ab_sym = vv->get_tuple_irrep(abc.ind_abs<0>(),abc.ind_abs<1>());
                  int    bc_sym = vv->get_tuple_irrep(abc.ind_abs<1>(),abc.ind_abs<2>());
                  size_t  a_rel = v->get_tuple_rel_index(abc.ind_abs<0>());
                  size_t  c_rel = v->get_tuple_rel_index(abc.ind_abs<2>());
                  size_t ab_rel = vv->get_tuple_rel_index(abc.ind_abs<0>(),abc.ind_abs<1>());
                  size_t bc_rel = vv->get_tuple_rel_index(abc.ind_abs<1>(),abc.ind_abs<2>());

                  double D_abC = e_vv_mu[a_abs] + e_vv_mu[b_abs] + e_VV_mu[c_abs];

                  // Update T
                  T[mu][ijk_sym]->set(a_sym,a_rel,bc_rel,Z[mu][ijk_sym]->get(a_sym,a_rel,bc_rel)/(Mk_shift[mu] + D_ijK - D_abC));

                  // Compute the energy
                  e4T[mu] += W_ijk[mu][ijk_sym]->get(a_sym,a_rel,bc_rel) * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel) / 2.0;
                  if((i_sym == a_sym) & (jk_sym == bc_sym)){
                    e4ST[mu] += T1_ov_mu[i_sym][i_rel][a_rel] *  V_oOvV[jk_sym][jk_rel][bc_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                    e4DT[mu] +=  F_ov_mu[i_sym][i_rel][a_rel] * T2_oOvV_mu[jk_sym][jk_rel][bc_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                  }
                  if((j_sym == a_sym) & (ik_sym == bc_sym)){
                    e4ST[mu] -= T1_ov_mu[j_sym][j_rel][a_rel] *  V_oOvV[ik_sym][ik_rel][bc_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                    e4DT[mu] -=  F_ov_mu[j_sym][j_rel][a_rel] * T2_oOvV_mu[ik_sym][ik_rel][bc_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                  }
                  if((k_sym == c_sym) & (ij_sym == ab_sym)){
                    e4ST[mu] += 0.5 * T1_OV_mu[k_sym][k_rel][c_rel] *  V_oovv[ij_sym][ij_rel][ab_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                    e4DT[mu] += 0.5 *  F_OV_mu[k_sym][k_rel][c_rel] * T2_oovv_mu[ij_sym][ij_rel][ab_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                  }
                }
              }  // End loop over abc
              newE += std::fabs(e4T[mu]) + std::fabs(e4ST[mu]) + std::fabs(e4DT[mu]);
            }  // End loop over allowed ijk
          }  // End of iterations
        }

        //      // Compute the contributions to the off-diagonal elements of Heff
        //      for(int mu = 0; mu < nrefs; ++mu){
        //        compute_ooO_contribution_to_Heff_spin_adapted(i_abs,j_abs,k_abs,mu,T[mu][ijk_sym]);
        //      }

        // Add the energy contributions from ijk
        for(int mu = 0; mu < nrefs; ++mu){
          E4T_ooO[mu]  += e4T[mu];
          E4ST_ooO[mu] += e4ST[mu];
          E4DT_ooO[mu] += e4DT[mu];
          E4T_oOO[mu]  += e4T[mu];
          E4ST_oOO[mu] += e4ST[mu];
          E4DT_oOO[mu] += e4DT[mu];
        }
      }


      //////////////
      // IKJ CASE //
      //////////////
      if((i_abs <= j_abs) and (j_abs < k_abs)){

        for(int mu = 0; mu < nrefs; ++mu){
          T[mu][ijk_sym]->zero();
        }

        // Compute T (d^2 N^6)
        int    cycle = 0;
        double oldE  = 1.0;
        double newE  = 0.0;
        while(fabs(oldE-newE) > threshold){
          cycle++;
          oldE = newE;
          newE = 0.0;
          // Iterate the Mk-MRCCSD(T) Equations
          for(int mu = 0; mu < nrefs; ++mu){
            e4T[mu] = e4ST[mu] = e4DT[mu] = 0.0;
            // Check if ijk belong to the occupied space of mu
            if(is_aocc[mu][i_abs] && is_aocc[mu][j_abs] && is_bocc[mu][k_abs]){

              double***  F_ov_mu    = F_ov[mu];
              double***  F_OV_mu    = F_OV[mu];
              double***  T1_ov_mu   = T1_ov[mu];
              double***  T1_OV_mu   = T1_OV[mu];
              double***  T2_oovv_mu = T2_oovv[mu];
              double***  T2_oOvV_mu = T2_oOvV[mu];

              double D_ijK = e_oo[mu][i_abs] + e_oo[mu][j_abs] + e_OO[mu][k_abs];

              // Add W
              Z[mu][ijk_sym]->add(W_ikj[mu][ijk_sym],0.0,1.0);

              // Add the coupling terms
              for(int nu = 0; nu < nrefs; ++nu){
                if(nu != mu){
                  Z[mu][ijk_sym]->add(T[nu][ijk_sym],1.0,Mk_factor[mu][nu]);
                }
              }

              // Divide by the denominator
              std::vector<double>& e_vv_mu  = e_vv[mu];
              std::vector<double>& e_VV_mu  = e_VV[mu];
              std::vector<bool>& is_avir_mu = is_avir[mu];
              std::vector<bool>& is_bvir_mu = is_bvir[mu];

              CCIndexIterator  abc(vvv,ijk_sym);
              // Loop over abc
              for(abc.first(); !abc.end(); abc.next()){
                size_t a_abs = v->get_tuple_abs_index(abc.ind_abs<0>());
                size_t b_abs = v->get_tuple_abs_index(abc.ind_abs<1>());
                size_t c_abs = v->get_tuple_abs_index(abc.ind_abs<2>());
                if(is_avir_mu[a_abs] && is_avir_mu[b_abs] && is_bvir_mu[c_abs]){
                  int     a_sym = v->get_tuple_irrep(abc.ind_abs<0>());
                  int     c_sym = v->get_tuple_irrep(abc.ind_abs<2>());
                  int    ab_sym = vv->get_tuple_irrep(abc.ind_abs<0>(),abc.ind_abs<1>());
                  int    bc_sym = vv->get_tuple_irrep(abc.ind_abs<1>(),abc.ind_abs<2>());
                  size_t  a_rel = v->get_tuple_rel_index(abc.ind_abs<0>());
                  size_t  c_rel = v->get_tuple_rel_index(abc.ind_abs<2>());
                  size_t ab_rel = vv->get_tuple_rel_index(abc.ind_abs<0>(),abc.ind_abs<1>());
                  size_t bc_rel = vv->get_tuple_rel_index(abc.ind_abs<1>(),abc.ind_abs<2>());

                  double D_abC = e_vv_mu[a_abs] + e_vv_mu[b_abs] + e_VV_mu[c_abs];

                  // Update T
                  T[mu][ijk_sym]->set(a_sym,a_rel,bc_rel,Z[mu][ijk_sym]->get(a_sym,a_rel,bc_rel)/(Mk_shift[mu] + D_ijK - D_abC));

                  // Compute the energy
                  e4T[mu] += W_ikj[mu][ijk_sym]->get(a_sym,a_rel,bc_rel) * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel) / 2.0;
                  if((i_sym == a_sym) & (jk_sym == bc_sym)){
                    e4ST[mu] += T1_ov_mu[i_sym][i_rel][a_rel] *  V_oOvV[jk_sym][kj_rel][bc_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel); //
                    e4DT[mu] +=  F_ov_mu[i_sym][i_rel][a_rel] * T2_oOvV_mu[jk_sym][kj_rel][bc_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                  }
                  if((k_sym == a_sym) & (ij_sym == bc_sym)){
                    e4ST[mu] -= T1_ov_mu[k_sym][k_rel][a_rel] *  V_oOvV[ij_sym][ij_rel][bc_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel); //
                    e4DT[mu] -=  F_ov_mu[k_sym][k_rel][a_rel] * T2_oOvV_mu[ij_sym][ij_rel][bc_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                  }
                  if((j_sym == c_sym) & (ik_sym == ab_sym)){
                    e4ST[mu] += 0.5 * T1_OV_mu[j_sym][j_rel][c_rel] *  V_oovv[ik_sym][ik_rel][ab_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                    e4DT[mu] += 0.5 *  F_OV_mu[j_sym][j_rel][c_rel] * T2_oovv_mu[ik_sym][ik_rel][ab_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                  }
                }
              }  // End loop over abc
              newE += std::fabs(e4T[mu]) + std::fabs(e4ST[mu]) + std::fabs(e4DT[mu]);
            }  // End loop over allowed ijk
          }  // End of iterations
        }

        //      // Compute the contributions to the off-diagonal elements of Heff
        //      for(int mu = 0; mu < nrefs; ++mu){
        //        compute_ooO_contribution_to_Heff_spin_adapted(i_abs,j_abs,k_abs,mu,T[mu][ijk_sym]);
        //      }

        // Add the energy contributions from ijk
        for(int mu = 0; mu < nrefs; ++mu){
          E4T_ooO[mu]  += e4T[mu];
          E4ST_ooO[mu] += e4ST[mu];
          E4DT_ooO[mu] += e4DT[mu];
          E4T_oOO[mu]  += e4T[mu];
          E4ST_oOO[mu] += e4ST[mu];
          E4DT_oOO[mu] += e4DT[mu];
        }
      }


      //////////////
      // JKI CASE //
      //////////////
      if((i_abs < j_abs) and (j_abs < k_abs)){

        for(int mu = 0; mu < nrefs; ++mu){
          T[mu][ijk_sym]->zero();
        }

        // Compute T (d^2 N^6)
        int    cycle = 0;
        double oldE  = 1.0;
        double newE  = 0.0;
        while(fabs(oldE-newE) > threshold){
          cycle++;
          oldE = newE;
          newE = 0.0;
          // Iterate the Mk-MRCCSD(T) Equations
          for(int mu = 0; mu < nrefs; ++mu){
            e4T[mu] = e4ST[mu] = e4DT[mu] = 0.0;
            // Check if ijk belong to the occupied space of mu
            if(is_aocc[mu][i_abs] && is_aocc[mu][j_abs] && is_bocc[mu][k_abs]){

              double***  F_ov_mu    = F_ov[mu];
              double***  F_OV_mu    = F_OV[mu];
              double***  T1_ov_mu   = T1_ov[mu];
              double***  T1_OV_mu   = T1_OV[mu];
              double***  T2_oovv_mu = T2_oovv[mu];
              double***  T2_oOvV_mu = T2_oOvV[mu];

              double D_ijK = e_oo[mu][i_abs] + e_oo[mu][j_abs] + e_OO[mu][k_abs];

              // Add W
              Z[mu][ijk_sym]->add(W_jki[mu][ijk_sym],0.0,1.0);

              // Add the coupling terms
              for(int nu = 0; nu < nrefs; ++nu){
                if(nu != mu){
                  Z[mu][ijk_sym]->add(T[nu][ijk_sym],1.0,Mk_factor[mu][nu]);
                }
              }

              // Divide by the denominator
              std::vector<double>& e_vv_mu  = e_vv[mu];
              std::vector<double>& e_VV_mu  = e_VV[mu];
              std::vector<bool>& is_avir_mu = is_avir[mu];
              std::vector<bool>& is_bvir_mu = is_bvir[mu];

              CCIndexIterator  abc(vvv,ijk_sym);
              // Loop over abc
              for(abc.first(); !abc.end(); abc.next()){
                size_t a_abs = v->get_tuple_abs_index(abc.ind_abs<0>());
                size_t b_abs = v->get_tuple_abs_index(abc.ind_abs<1>());
                size_t c_abs = v->get_tuple_abs_index(abc.ind_abs<2>());
                if(is_avir_mu[a_abs] && is_avir_mu[b_abs] && is_bvir_mu[c_abs]){
                  int     a_sym = v->get_tuple_irrep(abc.ind_abs<0>());
                  int     c_sym = v->get_tuple_irrep(abc.ind_abs<2>());
                  int    ab_sym = vv->get_tuple_irrep(abc.ind_abs<0>(),abc.ind_abs<1>());
                  int    bc_sym = vv->get_tuple_irrep(abc.ind_abs<1>(),abc.ind_abs<2>());
                  size_t  a_rel = v->get_tuple_rel_index(abc.ind_abs<0>());
                  size_t  c_rel = v->get_tuple_rel_index(abc.ind_abs<2>());
                  size_t ab_rel = vv->get_tuple_rel_index(abc.ind_abs<0>(),abc.ind_abs<1>());
                  size_t bc_rel = vv->get_tuple_rel_index(abc.ind_abs<1>(),abc.ind_abs<2>());

                  double D_abC = e_vv_mu[a_abs] + e_vv_mu[b_abs] + e_VV_mu[c_abs];

                  // Update T
                  T[mu][ijk_sym]->set(a_sym,a_rel,bc_rel,Z[mu][ijk_sym]->get(a_sym,a_rel,bc_rel)/(Mk_shift[mu] + D_ijK - D_abC));

                  // Compute the energy
                  e4T[mu] += W_jki[mu][ijk_sym]->get(a_sym,a_rel,bc_rel) * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel) / 2.0;
                  if((j_sym == a_sym) & (ik_sym == bc_sym)){
                    e4ST[mu] += T1_ov_mu[j_sym][j_rel][a_rel] *  V_oOvV[ik_sym][ki_rel][bc_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                    e4DT[mu] +=  F_ov_mu[j_sym][j_rel][a_rel] * T2_oOvV_mu[ik_sym][ki_rel][bc_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                  }
                  if((k_sym == a_sym) & (ij_sym == bc_sym)){
                    e4ST[mu] -= T1_ov_mu[k_sym][k_rel][a_rel] *  V_oOvV[ij_sym][ji_rel][bc_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                    e4DT[mu] -=  F_ov_mu[k_sym][k_rel][a_rel] * T2_oOvV_mu[ij_sym][ji_rel][bc_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                  }
                  if((i_sym == c_sym) & (jk_sym == ab_sym)){
                    e4ST[mu] += 0.5 * T1_OV_mu[i_sym][i_rel][c_rel] *  V_oovv[jk_sym][jk_rel][ab_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                    e4DT[mu] += 0.5 *  F_OV_mu[i_sym][i_rel][c_rel] * T2_oovv_mu[jk_sym][jk_rel][ab_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                  }
                }
              }  // End loop over abc
              newE += std::fabs(e4T[mu]) + std::fabs(e4ST[mu]) + std::fabs(e4DT[mu]);
            }  // End loop over allowed ijk
          }  // End of iterations
        }

        //      // Compute the contributions to the off-diagonal elements of Heff
        //      for(int mu = 0; mu < nrefs; ++mu){
        //        compute_ooO_contribution_to_Heff_spin_adapted(i_abs,j_abs,k_abs,mu,T[mu][ijk_sym]);
        //      }

        // Add the energy contributions from ijk
        for(int mu = 0; mu < nrefs; ++mu){
          E4T_ooO[mu]  += e4T[mu];
          E4ST_ooO[mu] += e4ST[mu];
          E4DT_ooO[mu] += e4DT[mu];
          E4T_oOO[mu]  += e4T[mu];
          E4ST_oOO[mu] += e4ST[mu];
          E4DT_oOO[mu] += e4DT[mu];
        }
      }


      //////////////
      // AAA CASE //
      //////////////
      if((i_abs < j_abs) and (j_abs < k_abs)){

        for(int mu = 0; mu < nrefs; ++mu){
          T[mu][ijk_sym]->zero();
        }

        // Compute T (d^2 N^6)
        int    cycle = 0;
        double oldE  = 1.0;
        double newE  = 0.0;
        while(fabs(oldE-newE) > threshold){
          cycle++;
          oldE = newE;
          newE = 0.0;
          // Iterate the Mk-MRCCSD(T) Equations
          for(int mu = 0; mu < nrefs; ++mu){
            e4T[mu] = e4ST[mu] = e4DT[mu] = 0.0;
            // Check if ijk belong to the occupied space of mu
            if(is_aocc[mu][i_abs] && is_aocc[mu][j_abs] && is_bocc[mu][k_abs]){

              double***  F_ov_mu    = F_ov[mu];
              double***  F_OV_mu    = F_OV[mu];
              double***  T1_ov_mu   = T1_ov[mu];
              double***  T1_OV_mu   = T1_OV[mu];
              double***  T2_oovv_mu = T2_oovv[mu];
              double***  T2_oOvV_mu = T2_oOvV[mu];

              double D_ijK = e_oo[mu][i_abs] + e_oo[mu][j_abs] + e_OO[mu][k_abs];

              // Add W
              Z[mu][ijk_sym]->add(W_ijk[mu][ijk_sym],0.0, 1.0);
              Z[mu][ijk_sym]->add(W_jki[mu][ijk_sym],1.0, 1.0);
              Z[mu][ijk_sym]->add(W_ikj[mu][ijk_sym],1.0,-1.0);

              // Add the coupling terms
              for(int nu = 0; nu < nrefs; ++nu){
                if(nu != mu){
                  Z[mu][ijk_sym]->add(T[nu][ijk_sym],1.0,Mk_factor[mu][nu]);
                }
              }

              // Divide by the denominator
              std::vector<double>& e_vv_mu  = e_vv[mu];
              std::vector<double>& e_VV_mu  = e_VV[mu];
              std::vector<bool>& is_avir_mu = is_avir[mu];
              std::vector<bool>& is_bvir_mu = is_bvir[mu];

              CCIndexIterator  abc(vvv,ijk_sym);
              // Loop over abc
              for(abc.first(); !abc.end(); abc.next()){
                size_t a_abs = v->get_tuple_abs_index(abc.ind_abs<0>());
                size_t b_abs = v->get_tuple_abs_index(abc.ind_abs<1>());
                size_t c_abs = v->get_tuple_abs_index(abc.ind_abs<2>());
                if(is_avir_mu[a_abs] && is_avir_mu[b_abs] && is_bvir_mu[c_abs]){
                  int     a_sym = v->get_tuple_irrep(abc.ind_abs<0>());
                  int     c_sym = v->get_tuple_irrep(abc.ind_abs<2>());
                  int    ab_sym = vv->get_tuple_irrep(abc.ind_abs<0>(),abc.ind_abs<1>());
                  int    bc_sym = vv->get_tuple_irrep(abc.ind_abs<1>(),abc.ind_abs<2>());
                  size_t  a_rel = v->get_tuple_rel_index(abc.ind_abs<0>());
                  size_t  c_rel = v->get_tuple_rel_index(abc.ind_abs<2>());
                  size_t ab_rel = vv->get_tuple_rel_index(abc.ind_abs<0>(),abc.ind_abs<1>());
                  size_t bc_rel = vv->get_tuple_rel_index(abc.ind_abs<1>(),abc.ind_abs<2>());

                  double D_abC = e_vv_mu[a_abs] + e_vv_mu[b_abs] + e_VV_mu[c_abs];

                  // Update T
                  T[mu][ijk_sym]->set(a_sym,a_rel,bc_rel,Z[mu][ijk_sym]->get(a_sym,a_rel,bc_rel)/(Mk_shift[mu] + D_ijK - D_abC));

                  // Compute the energy
                  e4T[mu] += W_ijk[mu][ijk_sym]->get(a_sym,a_rel,bc_rel) * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel) / 6.0;
                  e4T[mu] += W_jki[mu][ijk_sym]->get(a_sym,a_rel,bc_rel) * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel) / 6.0;
                  e4T[mu] -= W_ikj[mu][ijk_sym]->get(a_sym,a_rel,bc_rel) * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel) / 6.0;

                  if((i_sym == a_sym) & (jk_sym == bc_sym)){
                    e4ST[mu] += T1_ov_mu[i_sym][i_rel][a_rel] *  V_oOvV[jk_sym][jk_rel][bc_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                    e4DT[mu] +=  F_ov_mu[i_sym][i_rel][a_rel] * T2_oOvV_mu[jk_sym][jk_rel][bc_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                  }
                  if((j_sym == a_sym) & (ik_sym == bc_sym)){
                    e4ST[mu] -= T1_ov_mu[j_sym][j_rel][a_rel] *  V_oOvV[ik_sym][ik_rel][bc_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                    e4DT[mu] -=  F_ov_mu[j_sym][j_rel][a_rel] * T2_oOvV_mu[ik_sym][ik_rel][bc_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                  }
                  if((k_sym == c_sym) & (ij_sym == ab_sym)){
                    e4ST[mu] += 0.5 * T1_OV_mu[k_sym][k_rel][c_rel] *  V_oovv[ij_sym][ij_rel][ab_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                    e4DT[mu] += 0.5 *  F_OV_mu[k_sym][k_rel][c_rel] * T2_oovv_mu[ij_sym][ij_rel][ab_rel] * T[mu][ijk_sym]->get(a_sym,a_rel,bc_rel);
                  }
                }
              }  // End loop over abc
              newE += std::fabs(e4T[mu]) + std::fabs(e4ST[mu]) + std::fabs(e4DT[mu]);
            }  // End loop over allowed ijk
          }  // End of iterations
        }

        //      // Compute the contributions to the off-diagonal elements of Heff
        //      for(int mu = 0; mu < nrefs; ++mu){
        //        compute_ooO_contribution_to_Heff_spin_adapted(i_abs,j_abs,k_abs,mu,T[mu][ijk_sym]);
        //      }

        // Add the energy contributions from ijk
        for(int mu = 0; mu < nrefs; ++mu){
          E4T_ooo[mu]  += e4T[mu];
          E4ST_ooo[mu] += e4ST[mu];
          E4DT_ooo[mu] += e4DT[mu];
          E4T_OOO[mu]  += e4T[mu];
          E4ST_OOO[mu] += e4ST[mu];
          E4DT_OOO[mu] += e4DT[mu];
        }
      }

    }
  } // End loop over ijk

  for(int mu = 0; mu < nrefs; ++mu){
    fprintf(outfile,"\n  E_T[4]  (aaa) = %20.15lf (%d)",E4T_ooo[mu],mu);
    fprintf(outfile,"\n  E_ST[4] (aaa) = %20.15lf (%d)",E4ST_ooo[mu],mu);
    fprintf(outfile,"\n  E_DT[4] (aaa) = %20.15lf (%d)",E4DT_ooo[mu],mu);

    fprintf(outfile,"\n  E_T[4]  (aab) = %20.15lf (%d)",E4T_ooO[mu],mu);
    fprintf(outfile,"\n  E_ST[4] (aab) = %20.15lf (%d)",E4ST_ooO[mu],mu);
    fprintf(outfile,"\n  E_DT[4] (aab) = %20.15lf (%d)",E4DT_ooO[mu],mu);
    E4_ooo[mu] = E4T_ooo[mu] + E4ST_ooo[mu] + E4DT_ooo[mu];
    E4_ooO[mu] = E4T_ooO[mu] + E4ST_ooO[mu] + E4DT_ooO[mu];
    E4_oOO[mu] = E4T_oOO[mu] + E4ST_oOO[mu] + E4DT_oOO[mu];
    E4_OOO[mu] = E4T_OOO[mu] + E4ST_OOO[mu] + E4DT_OOO[mu];
  }
}

}} /* End Namespaces */












//// 1. x_ijk^abc
//Z[mu][ijk_sym]->contract(T2_iJ_a_B->get_block_matrix(ij_abs,mu),V_K_bC_e->get_block_matrix(k_abs),1.0,0.0);
//// + x_ijk^abc - x_ijk^bac                                           abc  acb  bac  bca  cab  cba
//W_ijk[mu][ijk_sym]->add_permutation_1_2(0.0,Z[mu][ijk_sym],vvv,v,vv, 1.0, 0.0,-1.0, 0.0, 0.0, 0.0);
//// + x_ijk^acb - x_ijk^bca                                           abc  acb  bac  bca  cab  cba
//W_ikj[mu][ijk_sym]->add_permutation_1_2(0.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 1.0, 0.0,-1.0, 0.0, 0.0);
//// + x_ijk^cab - x_ijk^cba                                           abc  acb  bac  bca  cab  cba
//W_jki[mu][ijk_sym]->add_permutation_1_2(0.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0);
//
//// 2. x_jik^abc
//Z[mu][ijk_sym]->contract(T2_iJ_a_B->get_block_matrix(ji_abs,mu),V_K_bC_e->get_block_matrix(k_abs),1.0,0.0);
//// - x_jik^abc + x_jik^bac                                           abc  acb  bac  bca  cab  cba
//W_ijk[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0, 0.0, 1.0, 0.0, 0.0, 0.0);
//// + x_jik^cab - x_jik^cba                                           abc  acb  bac  bca  cab  cba
//W_ikj[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0);
//// + x_jik^acb - x_jik^bca                                           abc  acb  bac  bca  cab  cba
//W_jki[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 1.0, 0.0,-1.0, 0.0, 0.0);
//
//// 3. x_jki^abc
//Z[mu][ijk_sym]->contract(T2_iJ_a_B->get_block_matrix(jk_abs,mu),V_K_bC_e->get_block_matrix(i_abs),1.0,0.0);
//// - x_jki^acb + x_jki^bca                                           abc  acb  bac  bca  cab  cba
//W_ijk[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0,-1.0, 0.0, 1.0, 0.0, 0.0);
//// - x_jki^cab + x_jki^cba                                           abc  acb  bac  bca  cab  cba
//W_ikj[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 0.0, 0.0, 0.0,-1.0, 1.0);
//// + x_jki^abc - x_jki^bac                                           abc  acb  bac  bca  cab  cba
//W_jki[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 1.0, 0.0,-1.0, 0.0, 0.0, 0.0);
//
//// 4. x_kji^abc
//Z[mu][ijk_sym]->contract(T2_iJ_a_B->get_block_matrix(kj_abs,mu),V_K_bC_e->get_block_matrix(i_abs),1.0,0.0);
//// - x_kji^cab + x_kji^cba                                           abc  acb  bac  bca  cab  cba
//W_ijk[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 0.0, 0.0, 0.0,-1.0, 1.0);
//// - x_kji^acb + x_kji^bac                                           abc  acb  bac  bca  cab  cba
//W_ikj[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0,-1.0, 0.0, 1.0, 0.0, 0.0);
//// - x_kji^abc + x_kji^bac                                           abc  acb  bac  bca  cab  cba
//W_jki[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0, 0.0, 1.0, 0.0, 0.0, 0.0);
//
//// 5. x_kij^abc
//Z[mu][ijk_sym]->contract(T2_iJ_a_B->get_block_matrix(ki_abs,mu),V_K_bC_e->get_block_matrix(j_abs),1.0,0.0);
//// + x_kij^cab - x_kij^cba                                           abc  acb  bac  bca  cab  cba
//W_ijk[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0);
//// - x_kij^abc + x_kij^bac                                           abc  acb  bac  bca  cab  cba
//W_ikj[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0, 0.0, 1.0, 0.0, 0.0, 0.0);
//// - x_kij^acb + x_kij^bca                                           abc  acb  bac  bca  cab  cba
//W_jki[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0,-1.0, 0.0, 1.0, 0.0, 0.0);
//
//// 6. x_ikj^abc
//Z[mu][ijk_sym]->contract(T2_iJ_a_B->get_block_matrix(ik_abs,mu),V_K_bC_e->get_block_matrix(j_abs),1.0,0.0);
//// + x_ikj^acb - x_ikj^bca                                           abc  acb  bac  bca  cab  cba
//W_ijk[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 1.0, 0.0,-1.0, 0.0, 0.0);
//// + x_ikj^abc - x_ikj^bac                                           abc  acb  bac  bca  cab  cba
//W_ikj[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 1.0, 0.0,-1.0, 0.0, 0.0, 0.0);
//// - x_ikj^cab + x_ikj^cba                                           abc  acb  bac  bca  cab  cba
//W_jki[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 0.0, 0.0, 0.0,-1.0, 1.0);
//
//// 7. y_ijk^abc
//Z[mu][ijk_sym]->contract(V_jK_C_m->get_block_matrix(ij_abs),T2_i_aB_J->get_block_matrix(k_abs,mu),1.0,0.0);
//// + y_ijk^acb - y_ijk^bca                                           abc  acb  bac  bca  cab  cba
//W_ijk[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 1.0, 0.0,-1.0, 0.0, 0.0);
//// + y_ijk^cab - y_ijk^cba                                           abc  acb  bac  bca  cab  cba
//W_ikj[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0);
//// - y_ijk^abc + y_ijk^bac                                           abc  acb  bac  bca  cab  cba
//W_jki[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0, 0.0, 1.0, 0.0, 0.0, 0.0);
//
//// 8. y_jik^acb
//Z[mu][ijk_sym]->contract(V_jK_C_m->get_block_matrix(ji_abs),T2_i_aB_J->get_block_matrix(k_abs,mu),1.0,0.0);
//// - y_jik^acb + y_jik^bca                                           abc  acb  bac  bca  cab  cba
//W_ijk[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0,-1.0, 0.0, 1.0, 0.0, 0.0);
//// - y_jik^abc + y_jik^bac                                           abc  acb  bac  bca  cab  cba
//W_ikj[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0, 0.0, 1.0, 0.0, 0.0, 0.0);
//// + y_jik^cab - y_jik^cba                                           abc  acb  bac  bca  cab  cba
//W_jki[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0);
//
//// 9. y_jki^abc
//Z[mu][ijk_sym]->contract(V_jK_C_m->get_block_matrix(jk_abs),T2_i_aB_J->get_block_matrix(i_abs,mu),1.0,0.0);
//// - y_jki^cab + y_jki^cba                                           abc  acb  bac  bca  cab  cba
//W_ijk[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 0.0, 0.0, 0.0,-1.0, 1.0);
//// + y_jki^abc - y_jki^bac                                           abc  acb  bac  bca  cab  cba
//W_ikj[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 1.0, 0.0,-1.0, 0.0, 0.0, 0.0);
//// + y_jki^acb - y_jki^bca                                           abc  acb  bac  bca  cab  cba
//W_jki[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 1.0, 0.0,-1.0, 0.0, 0.0);
//
//// 10. y_kji^abc
//Z[mu][ijk_sym]->contract(V_jK_C_m->get_block_matrix(kj_abs),T2_i_aB_J->get_block_matrix(i_abs,mu),1.0,0.0);
//// + y_kji^abc - y_kji^bac                                           abc  acb  bac  bca  cab  cba
//W_ijk[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 1.0, 0.0,-1.0, 0.0, 0.0, 0.0);
//// - y_kji^cab + y_kji^cba                                           abc  acb  bac  bca  cab  cba
//W_ikj[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 0.0, 0.0, 0.0,-1.0, 1.0);
//// - y_kji^acb + y_kji^bca                                           abc  acb  bac  bca  cab  cba
//W_jki[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0,-1.0, 0.0, 1.0, 0.0, 0.0);
//
//// 11. y_kij^abc
//Z[mu][ijk_sym]->contract(V_jK_C_m->get_block_matrix(ki_abs),T2_i_aB_J->get_block_matrix(j_abs,mu),1.0,0.0);
//// - y_kij^abc + y_kij^bac                                           abc  acb  bac  bca  cab  cba
//W_ijk[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv,-1.0, 0.0, 1.0, 0.0, 0.0, 0.0);
//// - y_kij^acb + y_kij^bca                                           abc  acb  bac  bca  cab  cba
//W_ikj[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0,-1.0, 0.0, 1.0, 0.0, 0.0);
//// - y_kij^cab + y_kij^cba                                           abc  acb  bac  bca  cab  cba
//W_jki[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 0.0, 0.0, 0.0,-1.0, 1.0);
//
//// 12. y_ikj^abc
//Z[mu][ijk_sym]->contract(V_jK_C_m->get_block_matrix(ik_abs),T2_i_aB_J->get_block_matrix(j_abs,mu),1.0,0.0);
//// + y_ikj^cab - y_ikj^cba                                           abc  acb  bac  bca  cab  cba
//W_ijk[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0);
//// + y_ikj^acb - y_ikj^bca                                           abc  acb  bac  bca  cab  cba
//W_ikj[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 0.0, 1.0, 0.0,-1.0, 0.0, 0.0);
//// + y_ikj^abc - y_ikj^bac                                           abc  acb  bac  bca  cab  cba
//W_jki[mu][ijk_sym]->add_permutation_1_2(1.0,Z[mu][ijk_sym],vvv,v,vv, 1.0, 0.0,-1.0, 0.0, 0.0, 0.0);

