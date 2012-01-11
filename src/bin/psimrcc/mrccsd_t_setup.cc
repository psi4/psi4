#include <algorithm>
#include <cmath>
#include <utility>

#include <libmoinfo/libmoinfo.h>
#include <liboptions/liboptions.h>

#include "blas.h"
#include "heff.h"
#include "index_iterator.h"
#include "matrix.h"
#include "mrccsd_t.h"
#include "special_matrices.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;
    extern MemoryManager* memory_manager;


void MRCCSD_T::startup()
{
  if(options_.get_str("TRIPLES_ALGORITHM") == "SPIN_ADAPTED"){
    triples_algorithm = SpinAdaptedTriples;
  }else if(options_.get_str("TRIPLES_ALGORITHM") == "RESTRICTED"){
    triples_algorithm = RestrictedTriples;
  }else{
    triples_algorithm = UnrestrictedTriples;
  }

  nirreps   = moinfo->get_nirreps();
  nrefs     = moinfo->get_ref_size(AllRefs);
  threshold = 0.1 * options_.get_double("E_CONVERGENCE");

  build_W_intermediates();

  o   = blas->get_index("[o]");
  oo  = blas->get_index("[oo]");
  ov  = blas->get_index("[ov]");
  v   = blas->get_index("[v]");
  vo  = blas->get_index("[vo]");
  vv  = blas->get_index("[vv]");
  vvv = blas->get_index("[vvv]");
  ovv = blas->get_index("[ovv]");
  ooo = blas->get_index("[ooo]");

  T2_ij_a_b = new IndexMatrix();
  T2_iJ_a_B = new IndexMatrix();
  T2_iJ_B_a = new IndexMatrix();
  T2_IJ_A_B = new IndexMatrix();

  T2_i_ab_j = new IndexMatrix();
  T2_i_aB_J = new IndexMatrix();
  T2_J_aB_i = new IndexMatrix();
  T2_I_AB_J = new IndexMatrix();

  V_k_bc_e = new IndexMatrix();
  V_k_bC_E = new IndexMatrix();
  V_K_bC_e = new IndexMatrix();

  V_jk_c_m = new IndexMatrix();
  V_jK_c_M = new IndexMatrix();
  V_jK_C_m = new IndexMatrix();

  form_T2_ij_a_b(T2_ij_a_b,true,true,false);
  form_T2_ij_a_b(T2_iJ_a_B,true,false,false);
  form_T2_ij_a_b(T2_iJ_B_a,true,false,true);   // T2_iJ_B_a = t2_{iJ}^{aB}
  form_T2_ij_a_b(T2_IJ_A_B,false,false,false);

  form_T2_i_ab_j(T2_i_ab_j,true,true,false);   // T2_i_ab_j = t2_{ij}^{ab}
  form_T2_i_ab_j(T2_i_aB_J,true,false,false);  // T2_i_aB_J = t2_{iJ}^{aB}
  form_T2_i_ab_j(T2_J_aB_i,true,false,true);   // T2_J_aB_i = t2_{iJ}^{aB}
  form_T2_i_ab_j(T2_I_AB_J,false,false,false); // T2_I_AB_J = t2_{IJ}^{AB}

  form_V_k_bc_e(V_k_bc_e,1.0,-1.0); // = <bc:ek>
  form_V_k_bc_e(V_k_bC_E,0.0,1.0);  // = <bC|kE> = <Ek|Cb>
  form_V_k_bc_e(V_K_bC_e,1.0,0.0);  // = <bC|eK> = <eK|bC>

  form_V_jk_c_m(V_jk_c_m,1.0,-1.0); // = <jk:mc>
  form_V_jk_c_m(V_jK_c_M,0.0,1.0);  // = <jk|cm>
  form_V_jk_c_m(V_jK_C_m,1.0,0.0);  // = <jk|mc>

  if(options_.get_bool("FAVG_CCSD_T")){
    fprintf(outfile,"\n\n  Using the average Fock matrix for the all references\n");
    for(int mu = 0; mu < nrefs; ++mu){
      int unique_mu = moinfo->get_ref_number(mu,AllRefs);
      double c_mu_2 = h_eff->get_zeroth_order_eigenvector(unique_mu)
                    * h_eff->get_zeroth_order_eigenvector(unique_mu);
      std::string factor = to_string(c_mu_2);
      blas->solve("epsilon[o][o]{u} += " + factor + " fock[o][o]{" +  to_string(unique_mu) + "}");
      blas->solve("epsilon[O][O]{u} += " + factor + " fock[O][O]{" +  to_string(unique_mu) + "}");
      blas->solve("epsilon[v][v]{u} += " + factor + " fock[v][v]{" +  to_string(unique_mu) + "}");
      blas->solve("epsilon[V][V]{u} += " + factor + " fock[V][V]{" +  to_string(unique_mu) + "}");
    }
  }else{
    blas->solve("epsilon[o][o]{u} = fock[o][o]{u}");
    blas->solve("epsilon[O][O]{u} = fock[O][O]{u}");
    blas->solve("epsilon[v][v]{u} = fock[v][v]{u}");
    blas->solve("epsilon[V][V]{u} = fock[V][V]{u}");
  }

  for(int mu = 0; mu < nrefs; ++mu){
    int unique_mu = moinfo->get_ref_number(mu,AllRefs);

    //  Unique references
    if(mu == unique_mu){
      // Setup the denominators
      double*** F_oo = blas->get_MatTmp("epsilon[o][o]",mu,none)->get_matrix();
      std::vector<double>  e_oo_mu;
      {
        CCIndexIterator i("[o]");
        for(i.first(); !i.end(); i.next()){
          int    i_sym = o->get_tuple_irrep(i.ind_abs<0>());
          size_t i_rel = o->get_tuple_rel_index(i.ind_abs<0>());
          e_oo_mu.push_back(F_oo[i_sym][i_rel][i_rel]);
        }
      }
      e_oo.push_back(e_oo_mu);

      double*** F_OO = blas->get_MatTmp("epsilon[O][O]",mu,none)->get_matrix();
      std::vector<double>  e_OO_mu;
      {
        CCIndexIterator i("[o]");
        for(i.first(); !i.end(); i.next()){
          int    i_sym = o->get_tuple_irrep(i.ind_abs<0>());
          size_t i_rel = o->get_tuple_rel_index(i.ind_abs<0>());
          e_OO_mu.push_back(F_OO[i_sym][i_rel][i_rel]);
        }
      }
      e_OO.push_back(e_OO_mu);


      double*** F_vv = blas->get_MatTmp("epsilon[v][v]",mu,none)->get_matrix();
      std::vector<double>  e_vv_mu;
      {
        CCIndexIterator a("[v]");
        for(a.first(); !a.end(); a.next()){
          int    a_sym = v->get_tuple_irrep(a.ind_abs<0>());
          size_t a_rel = v->get_tuple_rel_index(a.ind_abs<0>());
          e_vv_mu.push_back(F_vv[a_sym][a_rel][a_rel]);
        }
      }
      e_vv.push_back(e_vv_mu);

      double*** F_VV = blas->get_MatTmp("epsilon[V][V]",mu,none)->get_matrix();
      std::vector<double>  e_VV_mu;
      {
        CCIndexIterator a("[v]");
        for(a.first(); !a.end(); a.next()){
          int    a_sym = v->get_tuple_irrep(a.ind_abs<0>());
          size_t a_rel = v->get_tuple_rel_index(a.ind_abs<0>());
          e_VV_mu.push_back(F_VV[a_sym][a_rel][a_rel]);
        }
      }
      e_VV.push_back(e_VV_mu);

      F_ov.push_back(blas->get_MatTmp("fock[o][v]",unique_mu,none)->get_matrix());
      F_OV.push_back(blas->get_MatTmp("fock[O][V]",unique_mu,none)->get_matrix());

      if(options_.get_bool("HEFF4")){
        F2_ov.push_back(blas->get_MatTmp("F_me[o][v]",unique_mu,none)->get_matrix());
        F2_OV.push_back(blas->get_MatTmp("F_ME[O][V]",unique_mu,none)->get_matrix());
      }else{
        F2_ov.push_back(blas->get_MatTmp("fock[o][v]",unique_mu,none)->get_matrix());
        F2_OV.push_back(blas->get_MatTmp("fock[O][V]",unique_mu,none)->get_matrix());
      }

      W_ooov.push_back(blas->get_MatTmp("W_ijka[oo][ov]",unique_mu,none)->get_matrix());
      W_oOoV.push_back(blas->get_MatTmp("W_iJkA[oO][oV]",unique_mu,none)->get_matrix());
      W_OoOv.push_back(blas->get_MatTmp("W_IjKa[Oo][Ov]",unique_mu,none)->get_matrix());
      W_OOOV.push_back(blas->get_MatTmp("W_IJKA[OO][OV]",unique_mu,none)->get_matrix());

      W_vovv.push_back(blas->get_MatTmp("W_aibc[v][ovv]",unique_mu,none)->get_matrix());
      W_vOvV.push_back(blas->get_MatTmp("W_aIbC[v][OvV]",unique_mu,none)->get_matrix());
      W_VoVv.push_back(blas->get_MatTmp("W_AiBc[V][oVv]",unique_mu,none)->get_matrix());
      W_VOVV.push_back(blas->get_MatTmp("W_AIBC[V][OVV]",unique_mu,none)->get_matrix());

      T1_ov.push_back(blas->get_MatTmp("t1[o][v]",mu,none)->get_matrix());
      T1_OV.push_back(blas->get_MatTmp("t1[O][V]",mu,none)->get_matrix());

      T2_oovv.push_back(blas->get_MatTmp("t2[oo][vv]",mu,none)->get_matrix());
      T2_oOvV.push_back(blas->get_MatTmp("t2[oO][vV]",mu,none)->get_matrix());
      T2_OOVV.push_back(blas->get_MatTmp("t2[OO][VV]",mu,none)->get_matrix());

      is_aocc.push_back(moinfo->get_is_aocc(mu,AllRefs));
      is_bocc.push_back(moinfo->get_is_bocc(mu,AllRefs));
      is_avir.push_back(moinfo->get_is_avir(mu,AllRefs));
      is_bvir.push_back(moinfo->get_is_bvir(mu,AllRefs));
    }else{
      // Setup the denominators
      double*** F_oo = blas->get_MatTmp("epsilon[O][O]",unique_mu,none)->get_matrix();
      std::vector<double>  e_oo_mu;
      {
        CCIndexIterator i("[o]");
        for(i.first(); !i.end(); i.next()){
          int    i_sym = o->get_tuple_irrep(i.ind_abs<0>());
          size_t i_rel = o->get_tuple_rel_index(i.ind_abs<0>());
          e_oo_mu.push_back(F_oo[i_sym][i_rel][i_rel]);
        }
      }
      e_oo.push_back(e_oo_mu);

      double*** F_OO = blas->get_MatTmp("epsilon[o][o]",unique_mu,none)->get_matrix();
      std::vector<double>  e_OO_mu;
      {
        CCIndexIterator i("[o]");
        for(i.first(); !i.end(); i.next()){
          int    i_sym = o->get_tuple_irrep(i.ind_abs<0>());
          size_t i_rel = o->get_tuple_rel_index(i.ind_abs<0>());
          e_OO_mu.push_back(F_OO[i_sym][i_rel][i_rel]);
        }
      }
      e_OO.push_back(e_OO_mu);

      double*** F_vv = blas->get_MatTmp("epsilon[V][V]",unique_mu,none)->get_matrix();
      std::vector<double>  e_vv_mu;
      {
        CCIndexIterator a("[v]");
        for(a.first(); !a.end(); a.next()){
          int    a_sym = v->get_tuple_irrep(a.ind_abs<0>());
          size_t a_rel = v->get_tuple_rel_index(a.ind_abs<0>());
          e_vv_mu.push_back(F_vv[a_sym][a_rel][a_rel]);
        }
      }
      e_vv.push_back(e_vv_mu);

      double*** F_VV = blas->get_MatTmp("epsilon[v][v]",unique_mu,none)->get_matrix();
      std::vector<double>  e_VV_mu;
      {
        CCIndexIterator a("[v]");
        for(a.first(); !a.end(); a.next()){
          int    a_sym = v->get_tuple_irrep(a.ind_abs<0>());
          size_t a_rel = v->get_tuple_rel_index(a.ind_abs<0>());
          e_VV_mu.push_back(F_VV[a_sym][a_rel][a_rel]);
        }
      }
      e_VV.push_back(e_VV_mu);

      F_ov.push_back(blas->get_MatTmp("fock[O][V]",unique_mu,none)->get_matrix());
      F_OV.push_back(blas->get_MatTmp("fock[o][v]",unique_mu,none)->get_matrix());

      if(options_.get_bool("HEFF4")){
        F2_ov.push_back(blas->get_MatTmp("F_ME[O][V]",unique_mu,none)->get_matrix());
        F2_OV.push_back(blas->get_MatTmp("F_me[o][v]",unique_mu,none)->get_matrix());
      }else{
        F2_ov.push_back(blas->get_MatTmp("fock[O][V]",unique_mu,none)->get_matrix());
        F2_OV.push_back(blas->get_MatTmp("fock[o][v]",unique_mu,none)->get_matrix());
      }

      W_ooov.push_back(blas->get_MatTmp("W_IJKA[OO][OV]",unique_mu,none)->get_matrix());
      W_oOoV.push_back(blas->get_MatTmp("W_IjKa[Oo][Ov]",unique_mu,none)->get_matrix());
      W_OoOv.push_back(blas->get_MatTmp("W_iJkA[oO][oV]",unique_mu,none)->get_matrix());
      W_OOOV.push_back(blas->get_MatTmp("W_ijka[oo][ov]",unique_mu,none)->get_matrix());

      W_vovv.push_back(blas->get_MatTmp("W_AIBC[V][OVV]",unique_mu,none)->get_matrix());
      W_vOvV.push_back(blas->get_MatTmp("W_AiBc[V][oVv]",unique_mu,none)->get_matrix());
      W_VoVv.push_back(blas->get_MatTmp("W_aIbC[v][OvV]",unique_mu,none)->get_matrix());
      W_VOVV.push_back(blas->get_MatTmp("W_aibc[v][ovv]",unique_mu,none)->get_matrix());

      T1_ov.push_back(blas->get_MatTmp("t1[O][V]",unique_mu,none)->get_matrix());
      T1_OV.push_back(blas->get_MatTmp("t1[o][v]",unique_mu,none)->get_matrix());

      T2_oovv.push_back(blas->get_MatTmp("t2[OO][VV]",unique_mu,none)->get_matrix());
      T2_oOvV.push_back(blas->get_MatTmp("t2[Oo][Vv]",unique_mu,none)->get_matrix());
      T2_OOVV.push_back(blas->get_MatTmp("t2[oo][vv]",unique_mu,none)->get_matrix());

      is_bocc.push_back(moinfo->get_is_aocc(unique_mu,AllRefs));
      is_aocc.push_back(moinfo->get_is_bocc(unique_mu,AllRefs));
      is_bvir.push_back(moinfo->get_is_avir(unique_mu,AllRefs));
      is_avir.push_back(moinfo->get_is_bvir(unique_mu,AllRefs));
    }

    if(options_.get_str("CORR_CCSD_T") == "STANDARD"){
        std::vector<double> factor_row;
      for(int nu = 0; nu < nrefs; ++nu){
        double c_mu   = h_eff->get_right_eigenvector(mu);
        double c_nu   = h_eff->get_right_eigenvector(nu);
        double factor = h_eff->get_matrix(mu,nu) * c_nu / c_mu;
        if(options_.get_bool("TIKHONOW_TRIPLES")){
          double omega  = static_cast<double>(options_.get_int("TIKHONOW_OMEGA")) / 1000.0;
          factor = h_eff->get_matrix(mu,nu) * c_nu * c_mu / (pow(c_mu,2.0) + pow(omega,2.0));
        }
        factor_row.push_back(factor);
      }
      Mk_factor.push_back(factor_row);

      Mk_shift.push_back(h_eff->get_eigenvalue() - h_eff->get_matrix(mu,mu));
    }else if(options_.get_str("CORR_CCSD_T") == "PITTNER"){
        std::vector<double> factor_row;
      for(int nu = 0; nu < nrefs; ++nu){
        factor_row.push_back(0.0);
      }
      Mk_factor.push_back(factor_row);

      Mk_shift.push_back(0.0);
    }

    std::vector<double> d_h_eff_row;
    for(int nu = 0; nu < nrefs; ++nu){
      d_h_eff_row.push_back(0.0);
    }
    d_h_eff.push_back(d_h_eff_row);
  }

  V_oovv = blas->get_MatTmp("<[oo]:[vv]>",none)->get_matrix();
  V_oOvV = blas->get_MatTmp("<[oo]|[vv]>",none)->get_matrix();

  // Allocate Z, this will hold the results
  allocate2(BlockMatrix**,Z,nrefs,nirreps);
  for(int mu = 0; mu < nrefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      Z[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
    }
  }

  // Allocate W
  if((triples_algorithm == UnrestrictedTriples) or (triples_algorithm == RestrictedTriples)){
    allocate2(BlockMatrix**,W,nrefs,nirreps);
    for(int mu = 0; mu < nrefs; ++mu){
      for(int h = 0; h < nirreps; ++h){
        W[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
      }
    }
  }else if(triples_algorithm == SpinAdaptedTriples){
    allocate2(BlockMatrix**,W_ijk,nrefs,nirreps);
    allocate2(BlockMatrix**,W_ikj,nrefs,nirreps);
    allocate2(BlockMatrix**,W_jki,nrefs,nirreps);
    for(int mu = 0; mu < nrefs; ++mu){
      for(int h = 0; h < nirreps; ++h){
        W_ijk[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
        W_ikj[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
        W_jki[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
      }
    }
  }

  // Allocate T
  allocate2(BlockMatrix**,T,nrefs,nirreps);
  for(int mu = 0; mu < nrefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      T[mu][h] = new BlockMatrix(nirreps,v->get_tuplespi(),vv->get_tuplespi(),h);
    }
  }

  e4T.assign(nrefs,0.0);
  e4ST.assign(nrefs,0.0);
  e4DT.assign(nrefs,0.0);

  E4T_ooo.assign(nrefs,0.0);
  E4T_ooO.assign(nrefs,0.0);
  E4T_oOO.assign(nrefs,0.0);
  E4T_OOO.assign(nrefs,0.0);
  E4ST_ooo.assign(nrefs,0.0);
  E4ST_ooO.assign(nrefs,0.0);
  E4ST_oOO.assign(nrefs,0.0);
  E4ST_OOO.assign(nrefs,0.0);
  E4DT_ooo.assign(nrefs,0.0);
  E4DT_ooO.assign(nrefs,0.0);
  E4DT_oOO.assign(nrefs,0.0);
  E4DT_OOO.assign(nrefs,0.0);
  E4_ooo.assign(nrefs,0.0);
  E4_ooO.assign(nrefs,0.0);
  E4_oOO.assign(nrefs,0.0);
  E4_OOO.assign(nrefs,0.0);
}

void MRCCSD_T::check_intruders()
{
    std::vector<int> occ_to_mo =  moinfo->get_occ_to_mo();
    std::vector<int> vir_to_mo =  moinfo->get_vir_to_mo();
  // Identify intruders
  for(int mu = 0; mu < nrefs; ++mu){
    std::vector<std::pair<double,std::vector<short> > > aaa_sample;
    std::vector<std::pair<double,std::vector<short> > > aab_sample;
    // Loop over ijk
    CCIndexIterator  ijk("[ooo]");
    for(ijk.first(); !ijk.end(); ijk.next()){
      size_t i_abs = o->get_tuple_abs_index(ijk.ind_abs<0>());
      size_t j_abs = o->get_tuple_abs_index(ijk.ind_abs<1>());
      size_t k_abs = o->get_tuple_abs_index(ijk.ind_abs<2>());
      int ijk_sym = ijk.sym();
      // Loop over abc
      CCIndexIterator  abc("[vvv]",ijk_sym);
      for(abc.first(); !abc.end(); abc.next()){
        size_t a_abs = v->get_tuple_abs_index(abc.ind_abs<0>());
        size_t b_abs = v->get_tuple_abs_index(abc.ind_abs<1>());
        size_t c_abs = v->get_tuple_abs_index(abc.ind_abs<2>());

        // AAA Case
        if(is_aocc[mu][i_abs] and is_aocc[mu][j_abs] and is_aocc[mu][k_abs]){
          if(is_avir[mu][a_abs] and is_avir[mu][b_abs] and is_avir[mu][c_abs]){
            if((i_abs < j_abs) and (j_abs < k_abs) and (a_abs < b_abs) and (b_abs < c_abs)){
              double D_ijk = e_oo[mu][i_abs] + e_oo[mu][j_abs] + e_oo[mu][k_abs];
              double D_abc = e_vv[mu][a_abs] + e_vv[mu][b_abs] + e_vv[mu][c_abs];
              double denominator = D_ijk - D_abc;
              if(abs(denominator) < 0.1){
                std::vector<short> T3;
                T3.push_back(i_abs);
                T3.push_back(j_abs);
                T3.push_back(k_abs);
                T3.push_back(a_abs);
                T3.push_back(b_abs);
                T3.push_back(c_abs);
                aaa_sample.push_back(make_pair(denominator,T3));
              }
            }
          }
        }

        // AAB Case
        if(is_aocc[mu][i_abs] and is_aocc[mu][j_abs] and is_bocc[mu][k_abs]){
          if(is_avir[mu][a_abs] and is_avir[mu][b_abs] and is_bvir[mu][c_abs]){
            if((i_abs < j_abs) and (a_abs < b_abs)){
              double D_ijk = e_oo[mu][i_abs] + e_oo[mu][j_abs] + e_OO[mu][k_abs];
              double D_abc = e_vv[mu][a_abs] + e_vv[mu][b_abs] + e_VV[mu][c_abs];
              double denominator = D_ijk - D_abc;
              if(abs(denominator) < 0.1){
                std::vector<short> T3;
                T3.push_back(i_abs);
                T3.push_back(j_abs);
                T3.push_back(k_abs);
                T3.push_back(a_abs);
                T3.push_back(b_abs);
                T3.push_back(c_abs);
                aab_sample.push_back(make_pair(denominator,T3));
              }
            }
          }
        }
      }  // End loop over abc
    }  // End loop over allowed ijk

    int max_aaa = std::min(10,static_cast<int>(aaa_sample.size()));
    if(max_aaa > 0){
      fprintf(outfile,"\n\n  Intruders diagnostics for reference %d, AAA triple excitations",mu);
      fprintf(outfile,"\n  has found the following denominators with absolute value < 0.1\n");
      std::sort(aaa_sample.begin(),aaa_sample.end());
      for(int n = 0; n < max_aaa; ++n){
        fprintf(outfile,"\n  [%3d][%3d][%3d] -> [%3d][%3d][%3d] = %12.8f",
            aaa_sample[n].second[0],aaa_sample[n].second[1],aaa_sample[n].second[2],
            aaa_sample[n].second[3],aaa_sample[n].second[4],aaa_sample[n].second[5],
            aaa_sample[n].first);
      }
      fprintf(outfile,"\n  please check your results.");
    }

    int max_aab = std::min(10,static_cast<int>(aab_sample.size()));
    if(max_aab > 0){
      fprintf(outfile,"\n\n  Intruders diagnostics for reference %d, AAB triple excitations",mu);
      fprintf(outfile,"\n  has found the following denominators with absolute value < 0.1\n");

      std::sort(aab_sample.begin(),aab_sample.end());
      for(int n = 0; n < max_aab; ++n){
        fprintf(outfile,"\n  [%3d][%3d][%3d] -> [%3d][%3d][%3d] = %12.8f",
            aab_sample[n].second[0],aab_sample[n].second[1],aab_sample[n].second[2],
            aab_sample[n].second[3],aab_sample[n].second[4],aab_sample[n].second[5],
            aab_sample[n].first);
      }
      fprintf(outfile,"\n\n  please check your results.");
    }
    if((max_aaa + max_aab) > 0 ){
      fprintf(outfile,"\n  Orbital labels used to describe the intruder triple excitations refer to"
                      "\n  the occupied (docc + actv) and the virtual (actv + extr) spaces.");

      fprintf(outfile,"\n\n  Printing occupied orbital energies for reference %d",mu);
      CCIndexIterator  i("[o]");
      fprintf(outfile,"\n   OCC   MO      e(alpha)         e(beta)");
      for(i.first(); !i.end();i.next()){
        fprintf(outfile,"\n  %4d %4d",i.ind_abs<0>(),occ_to_mo[i.ind_abs<0>()]);
        if(is_aocc[mu][i.ind_abs<0>()]){
          fprintf(outfile,"%15.9f  ",e_oo[mu][i.ind_abs<0>()]);
        }else{
          fprintf(outfile,"         ---   ");
        }
        if(is_bocc[mu][i.ind_abs<0>()]){
          fprintf(outfile,"%15.9f",e_OO[mu][i.ind_abs<0>()]);
        }else{
          fprintf(outfile,"           ---");
        }
      }
      fprintf(outfile,"\n\n  Printing virtual orbital energies for reference %d",mu);
      CCIndexIterator  a("[v]");
      fprintf(outfile,"\n   VIR   MO      e(alpha)         e(beta)");
      for(a.first(); !a.end();a.next()){
        fprintf(outfile,"\n  %4d %4d",a.ind_abs<0>(),vir_to_mo[a.ind_abs<0>()]);
        if(is_avir[mu][a.ind_abs<0>()]){
          fprintf(outfile,"%15.9f  ",e_vv[mu][a.ind_abs<0>()]);
        }else{
          fprintf(outfile,"         ---   ");
        }
        if(is_bvir[mu][a.ind_abs<0>()]){
          fprintf(outfile,"%15.9f",e_VV[mu][a.ind_abs<0>()]);
        }else{
          fprintf(outfile,"           ---");
        }
      }


    }
  }

  /*
  for(int mu = 0; mu < nrefs; ++mu){
    fprintf(outfile,"\n  @=@%d %.9f",mu,Mk_shift[mu]);
  }

  for(int mu = 0; mu < nrefs; ++mu){
    fprintf(outfile,"\n  @@@%d ",mu);
    // Loop over ijk
    CCIndexIterator  ijk("[ooo]");
    for(ijk.first(); !ijk.end(); ijk.next()){
      size_t i_abs = o->get_tuple_abs_index(ijk.ind_abs<0>());
      size_t j_abs = o->get_tuple_abs_index(ijk.ind_abs<1>());
      size_t k_abs = o->get_tuple_abs_index(ijk.ind_abs<2>());
      int ijk_sym = ijk.sym();
      // Loop over abc
      CCIndexIterator  abc("[vvv]",ijk_sym);
      for(abc.first(); !abc.end(); abc.next()){
        size_t a_abs = v->get_tuple_abs_index(abc.ind_abs<0>());
        size_t b_abs = v->get_tuple_abs_index(abc.ind_abs<1>());
        size_t c_abs = v->get_tuple_abs_index(abc.ind_abs<2>());

        // AAA Case
        if(is_aocc[mu][i_abs] and is_aocc[mu][j_abs] and is_aocc[mu][k_abs]){
          if(is_avir[mu][a_abs] and is_avir[mu][b_abs] and is_avir[mu][c_abs]){
            if((i_abs < j_abs) and (j_abs < k_abs) and (a_abs < b_abs) and (b_abs < c_abs)){
              double D_ijk = e_oo[mu][i_abs] + e_oo[mu][j_abs] + e_oo[mu][k_abs];
              double D_abc = e_vv[mu][a_abs] + e_vv[mu][b_abs] + e_vv[mu][c_abs];
              double denominator = D_ijk - D_abc;
              fprintf(outfile," %.6f",denominator);
            }
          }
        }

      }  // End loop over abc
    }  // End loop over allowed ijk
  }


  for(int mu = 0; mu < nrefs; ++mu){
    fprintf(outfile,"\n  @@#%d ",mu);
    // Loop over ijk
    CCIndexIterator  ijk("[ooo]");
    for(ijk.first(); !ijk.end(); ijk.next()){
      size_t i_abs = o->get_tuple_abs_index(ijk.ind_abs<0>());
      size_t j_abs = o->get_tuple_abs_index(ijk.ind_abs<1>());
      size_t k_abs = o->get_tuple_abs_index(ijk.ind_abs<2>());
      int ijk_sym = ijk.sym();
      // Loop over abc
      CCIndexIterator  abc("[vvv]",ijk_sym);
      for(abc.first(); !abc.end(); abc.next()){
        size_t a_abs = v->get_tuple_abs_index(abc.ind_abs<0>());
        size_t b_abs = v->get_tuple_abs_index(abc.ind_abs<1>());
        size_t c_abs = v->get_tuple_abs_index(abc.ind_abs<2>());


        // AAB Case
        if(is_aocc[mu][i_abs] and is_aocc[mu][j_abs] and is_bocc[mu][k_abs]){
          if(is_avir[mu][a_abs] and is_avir[mu][b_abs] and is_bvir[mu][c_abs]){
            if((i_abs < j_abs) and (a_abs < b_abs)){
              double D_ijk = e_oo[mu][i_abs] + e_oo[mu][j_abs] + e_OO[mu][k_abs];
              double D_abc = e_vv[mu][a_abs] + e_vv[mu][b_abs] + e_VV[mu][c_abs];
              double denominator = D_ijk - D_abc;
              fprintf(outfile," %.6f",denominator);
            }
          }
        }
      }  // End loop over abc
    }  // End loop over allowed ijk
  }
  */
}

void MRCCSD_T::build_W_intermediates()
{
  blas->solve("W_ijka[oo][ov]{u}  = <[oo]:[ov]>");
  if(options_.get_bool("HEFF4"))
    blas->solve("W_ijka[oo][ov]{u} += #4123# <[v]:[voo]> 1@2 t1[o][v]{u}");

  blas->solve("W_iJkA[oO][oV]{u}  = <[oo]|[ov]>");
  if(options_.get_bool("HEFF4"))
    blas->solve("W_iJkA[oO][oV]{u} += #4123# <[v]|[voo]> 1@2 t1[o][v]{u}");

  blas->solve("W_IjKa[Oo][Ov]{u}  = <[oo]|[ov]>");
  if(options_.get_bool("HEFF4"))
    blas->solve("W_IjKa[Oo][Ov]{u} += #4123# <[v]|[voo]> 1@2 t1[O][V]{u}");

  blas->solve("W_IJKA[OO][OV]{u}  = <[oo]:[ov]>");
  if(options_.get_bool("HEFF4"))
    blas->solve("W_IJKA[OO][OV]{u} += #4123# <[v]:[voo]> 1@2 t1[O][V]{u}");

  blas->solve("W_aibc[v][ovv]{u}  = <[v]:[ovv]>");
  if(options_.get_bool("HEFF4"))
    blas->solve("W_aibc[v][ovv]{u} += - t1[o][v]{u} 1@1 <[o]:[ovv]>");

  blas->solve("W_aIbC[v][OvV]{u}  = <[v]|[ovv]>");
  if(options_.get_bool("HEFF4"))
    blas->solve("W_aIbC[v][OvV]{u} += - t1[o][v]{u} 1@1 <[o]|[ovv]>");

  blas->solve("W_AiBc[V][oVv]{u}  = <[v]|[ovv]>");
  if(options_.get_bool("HEFF4"))
    blas->solve("W_AiBc[V][oVv]{u} += - t1[O][V]{u} 1@1 <[o]|[ovv]>");

  blas->solve("W_AIBC[V][OVV]{u}  = <[v]:[ovv]>");
  if(options_.get_bool("HEFF4"))
    blas->solve("W_AIBC[V][OVV]{u} += - t1[O][V]{u} 1@1 <[o]:[ovv]>");
}

void MRCCSD_T::cleanup()
{
  delete T2_ij_a_b;
  delete T2_iJ_a_B;
  delete T2_iJ_B_a;
  delete T2_IJ_A_B;

  delete T2_i_ab_j;
  delete T2_i_aB_J;
  delete T2_J_aB_i;
  delete T2_I_AB_J;

  delete V_k_bc_e;
  delete V_k_bC_E;
  delete V_K_bC_e;

  delete V_jk_c_m;
  delete V_jK_c_M;
  delete V_jK_C_m;

  // Deallocate Z
  for(int mu = 0; mu < nrefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      delete Z[mu][h];
    }
  }
  release2(Z);

  // Deallocate W
  if((triples_algorithm == UnrestrictedTriples) or (triples_algorithm == RestrictedTriples)){
    for(int mu = 0; mu < nrefs; ++mu){
      for(int h = 0; h < nirreps; ++h){
        delete W[mu][h];
      }
    }
    release2(W);
  }else if(triples_algorithm == SpinAdaptedTriples){
    for(int mu = 0; mu < nrefs; ++mu){
      for(int h = 0; h < nirreps; ++h){
        delete W_ijk[mu][h];  delete W_ikj[mu][h];  delete W_jki[mu][h];
      }
    }
    release2(W_ijk);  release2(W_ikj);  release2(W_jki);
  }

  // Deallocate T
  for(int mu = 0; mu < nrefs; ++mu){
    for(int h = 0; h < nirreps; ++h){
      delete T[mu][h];
    }
  }
  release2(T);
}

}} /* End Namespaces */
