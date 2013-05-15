/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include <libmoinfo/libmoinfo.h>
#include "mrcc.h"
#include "blas.h"
#include <libutil/libutil.h>
#include "algebra_interface.h"
#include "memory_manager.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

/**
 * @brief Computes the contraction
 * \f[ -P(ij/k) \sum_m t_{ijm}^{abc} f_{mk} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000)
 */
void CCMRCC::build_t3_ijkabc_amplitudes_diagram3()
{

//    Algorithm
//    Loop over references
//      Loop over irreps of t3_eqns[ijk][abc]
//        Loop over [abc]
//          t3_eqns[ijk][abc] <- ddot t[ij][m](abc) 2@2 Fock[k][m]
//        End loop over [abc]
//      End loop irrep
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);
    string ref_str  = to_string(unique_ref);

    // Grab the temporary matrices
    CCMatTmp  TijkabcMatTmp = blas->get_MatTmp("t3[ooo][vvv]",unique_ref,none);
    CCMatTmp  HijkabcMatTmp = blas->get_MatTmp("t3_eqns[ooo][vvv]",unique_ref,none);
    CCMatTmp  FockmeMatTmp  = blas->get_MatTmp("fock[o][o]",unique_ref,none);

    // Grab the indexing for t3[ijk][abc]
    short**   ooo_tuples  = TijkabcMatTmp->get_left()->get_tuples();

    double*** Tijkabc_matrix = TijkabcMatTmp->get_matrix();
    double*** Fockme_matrix  = FockmeMatTmp->get_matrix();
    CCIndex* oo_indexing = blas->get_index("[oo]");
    CCIndex* o_indexing  = blas->get_index("[o]");

    double ***H_ijk;
    double ***t_ijm;
//  Loop over irrep
    allocate1(double**,H_ijk,moinfo->get_nirreps());
    allocate1(double**,t_ijm,moinfo->get_nirreps());
    for(int h =0; h < moinfo->get_nirreps();h++){
      for(int ij_sym = 0; ij_sym < moinfo->get_nirreps(); ij_sym++){
        int m_sym = h^ij_sym;
        allocate2(double,H_ijk[ij_sym],oo_indexing->get_pairpi(ij_sym),o_indexing ->get_pairpi(m_sym));
        allocate2(double,t_ijm[ij_sym],oo_indexing->get_pairpi(ij_sym),o_indexing ->get_pairpi(m_sym));
      }
      size_t ooo_offset  = TijkabcMatTmp->get_left()->get_first(h);
      for(int abc = 0;abc<TijkabcMatTmp->get_right_pairpi(h);abc++){
        // Copy T_ijk(abc) into t_[ij][k]
        for(int ijk = 0;ijk<TijkabcMatTmp->get_left_pairpi(h);ijk++){
          int i = ooo_tuples[ooo_offset + ijk][0];
          int j = ooo_tuples[ooo_offset + ijk][1];
          int k = ooo_tuples[ooo_offset + ijk][2];
          int ij = oo_indexing->get_tuple_rel_index(i,j);
          int m  = o_indexing ->get_tuple_rel_index(k);
          int ij_sym = oo_indexing->get_tuple_irrep(i,j);
          t_ijm[ij_sym][ij][m] = Tijkabc_matrix[h][ijk][abc];
        }

        // Loop over irreps
        for(int ij_sym = 0; ij_sym < moinfo->get_nirreps(); ij_sym++){
          int m = oo_indexing->get_pairpi(ij_sym);
          int n = o_indexing ->get_pairpi(ij_sym^h);
          int k = n;
          if(m*n*k){
            double alpha = 1.0;
            double beta  = 0.0;
            C_DGEMM_22(m, n, k, alpha,&(t_ijm[ij_sym][0][0]), n, 
                    &(Fockme_matrix[ij_sym^h][0][0]), n, beta, &(H_ijk[ij_sym][0][0]),n);
          }
        }

        for(int ijk = 0;ijk<HijkabcMatTmp->get_left_pairpi(h);ijk++){
          int i = ooo_tuples[ooo_offset + ijk][0];
          int j = ooo_tuples[ooo_offset + ijk][1];
          int k = ooo_tuples[ooo_offset + ijk][2];
          int ij = oo_indexing->get_tuple_rel_index(i,j);
          int m  = o_indexing ->get_tuple_rel_index(k);
          int ij_sym = oo_indexing->get_tuple_irrep(i,j);
          HijkabcMatTmp->add_six_address_element_Pij_k(i,j,k,abc,-H_ijk[ij_sym][ij][m]);
        }
      } // End loop over abc
      for(int ij_sym = 0; ij_sym < moinfo->get_nirreps(); ij_sym++){
        int m_sym = h^ij_sym;
        release2(H_ijk[ij_sym]);
        release2(t_ijm[ij_sym]);
      }
    } // End loop over irrep
    release1(H_ijk);
    release1(t_ijm);

  } // End loop over references
//   blas->print("t3_test[ooo][vvv]{u}");
//   blas->solve("ERROR{u} = 1000.0 t3_test[ooo][vvv]{u} . t3_test[ooo][vvv]{u}");
//   blas->print("ERROR{u}");
}

/**
 * @brief Computes the contraction
 * \f[ -P(ij/K) \sum_m t_{ijm}^{abC} f_{mK} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000)
 * This algorithm hence performs
 * \f[ \sum_m t_{ijM}^{abC} f_{MK} \rightarrow \left \{ -\bar{H}_{ijK}^{abC} \right \} \f]
 * \f[ P(ij)\sum_m t_{imK}^{abC} f_{mj} \rightarrow \left \{ -\bar{H}_{ijK}^{abC},\bar{H}_{jiK}^{abC} \right \} \f]
 */
void CCMRCC::build_t3_ijKabC_amplitudes_diagram3()
{

//    Algorithm
//    Loop over references
//      Loop over irreps of t3_eqns[ijk][abc]
//        Loop over [abc]
//          t3_eqns[ijk][abc] <- ddot t[ij][m](abc) 2@2 Fock[k][m]
//        End loop over [abc]
//      End loop irrep
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);
    string ref_str  = to_string(unique_ref);

    // Grab the temporary matrices
    CCMatTmp  TijKabCMatTmp = blas->get_MatTmp("t3[ooO][vvV]",unique_ref,none);
    CCMatTmp  HijKabCMatTmp = blas->get_MatTmp("t3_eqns[ooO][vvV]",unique_ref,none);
    CCMatTmp  FockmeMatTmp  = blas->get_MatTmp("fock[o][o]",unique_ref,none);
    CCMatTmp  FockMEMatTmp  = blas->get_MatTmp("fock[O][O]",unique_ref,none);

    // Grab the indexing for t3[ijk][abc]
    short**   ooo_tuples  = TijKabCMatTmp->get_left()->get_tuples();

    double*** TijKabC_matrix = TijKabCMatTmp->get_matrix();
    double*** Fockme_matrix  = FockmeMatTmp->get_matrix();
    double*** FockME_matrix  = FockMEMatTmp->get_matrix();
    CCIndex* oo_indexing = blas->get_index("[oo]");
    CCIndex* o_indexing  = blas->get_index("[o]");

    double ***H_ijk;
    double ***t_ijm;
//  Loop over irrep
    for(int h =0; h < moinfo->get_nirreps();h++){
      allocate1(double**,H_ijk,moinfo->get_nirreps());
      allocate1(double**,t_ijm,moinfo->get_nirreps());
      for(int ij_sym = 0; ij_sym < moinfo->get_nirreps(); ij_sym++){
        int m_sym = h^ij_sym;
        allocate2(double,H_ijk[ij_sym],oo_indexing->get_pairpi(ij_sym),o_indexing ->get_pairpi(m_sym));
        allocate2(double,t_ijm[ij_sym],oo_indexing->get_pairpi(ij_sym),o_indexing ->get_pairpi(m_sym));
      }
      size_t ooo_offset  = TijKabCMatTmp->get_left()->get_first(h);

      for(int abc = 0;abc<TijKabCMatTmp->get_right_pairpi(h);abc++){
        // Copy T_ijm(abc) into t_[ij][m]
        for(int ijk = 0;ijk<TijKabCMatTmp->get_left_pairpi(h);ijk++){
          int i = ooo_tuples[ooo_offset + ijk][0];
          int j = ooo_tuples[ooo_offset + ijk][1];
          int k = ooo_tuples[ooo_offset + ijk][2];
          int ij = oo_indexing->get_tuple_rel_index(i,j);
          int m  = o_indexing ->get_tuple_rel_index(k);
          int ij_sym = oo_indexing->get_tuple_irrep(i,j);
          t_ijm[ij_sym][ij][m] = TijKabC_matrix[h][ijk][abc];
        }

        // Loop over irreps
        for(int ij_sym = 0; ij_sym < moinfo->get_nirreps(); ij_sym++){
          int m = oo_indexing->get_pairpi(ij_sym);
          int n = o_indexing ->get_pairpi(ij_sym^h);
          int k = n;
          if(m*n*k){
            double alpha = 1.0;
            double beta  = 0.0;
            C_DGEMM_22(m, n, k, alpha,&(t_ijm[ij_sym][0][0]), n, 
                    &(FockME_matrix[ij_sym^h][0][0]), n, beta, &(H_ijk[ij_sym][0][0]),n);
          }
        }

        for(int ijk = 0;ijk<HijKabCMatTmp->get_left_pairpi(h);ijk++){
          int i = ooo_tuples[ooo_offset + ijk][0];
          int j = ooo_tuples[ooo_offset + ijk][1];
          int k = ooo_tuples[ooo_offset + ijk][2];
          int ij = oo_indexing->get_tuple_rel_index(i,j);
          int m  = o_indexing ->get_tuple_rel_index(k);
          int ij_sym = oo_indexing->get_tuple_irrep(i,j);
          HijKabCMatTmp->add_six_address_element_abc(i,j,k,abc,-H_ijk[ij_sym][ij][m]);
        }

        // Copy T_ijk(abc) into t_[ik][j]
        for(int ijk = 0;ijk<TijKabCMatTmp->get_left_pairpi(h);ijk++){
          int i = ooo_tuples[ooo_offset + ijk][0];
          int j = ooo_tuples[ooo_offset + ijk][1];
          int k = ooo_tuples[ooo_offset + ijk][2];
          int ik = oo_indexing->get_tuple_rel_index(i,k);
          int m  = o_indexing ->get_tuple_rel_index(j);
          int ik_sym = oo_indexing->get_tuple_irrep(i,k);
          t_ijm[ik_sym][ik][m] = TijKabC_matrix[h][ijk][abc];
        }

        // Loop over irreps
        for(int ij_sym = 0; ij_sym < moinfo->get_nirreps(); ij_sym++){
          int m = oo_indexing->get_pairpi(ij_sym);
          int n = o_indexing ->get_pairpi(ij_sym^h);
          int k = n;
          if(m*n*k){
            double alpha = 1.0;
            double beta  = 0.0;
            C_DGEMM_22(m, n, k, alpha,&(t_ijm[ij_sym][0][0]), n, 
                    &(Fockme_matrix[ij_sym^h][0][0]), n, beta, &(H_ijk[ij_sym][0][0]),n);
          }
        }

        for(int ijk = 0;ijk<HijKabCMatTmp->get_left_pairpi(h);ijk++){
          int i = ooo_tuples[ooo_offset + ijk][0];
          int j = ooo_tuples[ooo_offset + ijk][1];
          int k = ooo_tuples[ooo_offset + ijk][2];
          int ik = oo_indexing->get_tuple_rel_index(i,k);
          int m  = o_indexing ->get_tuple_rel_index(j);
          int ik_sym = oo_indexing->get_tuple_irrep(i,k);
          HijKabCMatTmp->add_six_address_element_Pij_abc(i,j,k,abc,-H_ijk[ik_sym][ik][m]);
        }


      } // End loop over abc
      for(int ij_sym = 0; ij_sym < moinfo->get_nirreps(); ij_sym++){
        int m_sym = h^ij_sym;
        release2(H_ijk[ij_sym]);
        release2(t_ijm[ij_sym]);
      }
      release1(H_ijk);
      release1(t_ijm);
    } // End loop over irrep

  } // End loop over references
//   blas->print("t3_test[ooO][vvV]{u}");
//   blas->solve("ERROR{u} = 1000.0 t3_test[ooO][vvV]{u} . t3_test[ooO][vvV]{u}");
//   blas->print("ERROR{u}");
}


/**
 * @brief Computes the contraction
 * \f[ -P(iJ/K) \sum_m t_{iJm}^{aBC} f_{mK} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000)
 * This algorithm hence performs
 * \f[ \sum_m t_{iJM}^{abC} f_{MK} \rightarrow \left \{ -\bar{H}_{ijK}^{abC} \right \} \f]
 * \f[ P(ij)\sum_m t_{imK}^{abC} f_{mj} \rightarrow \left \{ -\bar{H}_{ijK}^{abC},\bar{H}_{jiK}^{abC} \right \} \f]
 */
void CCMRCC::build_t3_iJKaBC_amplitudes_diagram3()
{

//    Algorithm
//    Loop over references
//      Loop over irreps of t3_eqns[ijk][abc]
//        Loop over [abc]
//          t3_eqns[ijk][abc] <- ddot t[ij][m](abc) 2@2 Fock[k][m]
//        End loop over [abc]
//      End loop irrep
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);
    string ref_str  = to_string(unique_ref);

    // Grab the temporary matrices
    CCMatTmp  TiJKaBCMatTmp = blas->get_MatTmp("t3[oOO][vVV]",unique_ref,none);
    CCMatTmp  HiJKaBCMatTmp = blas->get_MatTmp("t3_eqns[oOO][vVV]",unique_ref,none);
    CCMatTmp  FockmeMatTmp  = blas->get_MatTmp("fock[o][o]",unique_ref,none);
    CCMatTmp  FockMEMatTmp  = blas->get_MatTmp("fock[O][O]",unique_ref,none);

    // Grab the indexing for t3[ijk][abc]
    short**   ooo_tuples  = TiJKaBCMatTmp->get_left()->get_tuples();

    double*** TiJKaBC_matrix = TiJKaBCMatTmp->get_matrix();
    double*** Fockme_matrix  = FockmeMatTmp->get_matrix();
    double*** FockME_matrix  = FockMEMatTmp->get_matrix();
    CCIndex* oo_indexing = blas->get_index("[oo]");
    CCIndex* o_indexing  = blas->get_index("[o]");

    double ***H_ijk;
    double ***t_ijm;
//  Loop over irrep
    for(int h =0; h < moinfo->get_nirreps();h++){
      allocate1(double**,H_ijk,moinfo->get_nirreps());
      allocate1(double**,t_ijm,moinfo->get_nirreps());
      for(int ij_sym = 0; ij_sym < moinfo->get_nirreps(); ij_sym++){
        int m_sym = h^ij_sym;
        allocate2(double,H_ijk[ij_sym],oo_indexing->get_pairpi(ij_sym),o_indexing ->get_pairpi(m_sym));
        allocate2(double,t_ijm[ij_sym],oo_indexing->get_pairpi(ij_sym),o_indexing ->get_pairpi(m_sym));
      }
      size_t ooo_offset  = TiJKaBCMatTmp->get_left()->get_first(h);

      for(int abc = 0;abc<TiJKaBCMatTmp->get_right_pairpi(h);abc++){
        // Copy T_ijm(abc) into t_[ij][m]
        for(int ijk = 0;ijk<TiJKaBCMatTmp->get_left_pairpi(h);ijk++){
          int i = ooo_tuples[ooo_offset + ijk][0];
          int j = ooo_tuples[ooo_offset + ijk][1];
          int k = ooo_tuples[ooo_offset + ijk][2];
          int ij = oo_indexing->get_tuple_rel_index(i,j);
          int m  = o_indexing ->get_tuple_rel_index(k);
          int ij_sym = oo_indexing->get_tuple_irrep(i,j);
          t_ijm[ij_sym][ij][m] = TiJKaBC_matrix[h][ijk][abc];
        }

        // Loop over irreps
        for(int ij_sym = 0; ij_sym < moinfo->get_nirreps(); ij_sym++){
          int m = oo_indexing->get_pairpi(ij_sym);
          int n = o_indexing ->get_pairpi(ij_sym^h);
          int k = n;
          if(m*n*k){
            double alpha = 1.0;
            double beta  = 0.0;
            C_DGEMM_22(m, n, k, alpha,&(t_ijm[ij_sym][0][0]), n, 
                    &(FockME_matrix[ij_sym^h][0][0]), n, beta, &(H_ijk[ij_sym][0][0]),n);
          }
        }

        for(int ijk = 0;ijk<HiJKaBCMatTmp->get_left_pairpi(h);ijk++){
          int i = ooo_tuples[ooo_offset + ijk][0];
          int j = ooo_tuples[ooo_offset + ijk][1];
          int k = ooo_tuples[ooo_offset + ijk][2];
          int ij = oo_indexing->get_tuple_rel_index(i,j);
          int m  = o_indexing ->get_tuple_rel_index(k);
          int ij_sym = oo_indexing->get_tuple_irrep(i,j);
          HiJKaBCMatTmp->add_six_address_element_Pjk_abc(i,j,k,abc,-H_ijk[ij_sym][ij][m]);
        }

        // Copy T_ijk(abc) into t_[jk][i]
        for(int ijk = 0;ijk<TiJKaBCMatTmp->get_left_pairpi(h);ijk++){
          int i = ooo_tuples[ooo_offset + ijk][0];
          int j = ooo_tuples[ooo_offset + ijk][1];
          int k = ooo_tuples[ooo_offset + ijk][2];
          int jk = oo_indexing->get_tuple_rel_index(j,k);
          int m  = o_indexing ->get_tuple_rel_index(i);
          int jk_sym = oo_indexing->get_tuple_irrep(j,k);
          t_ijm[jk_sym][jk][m] = TiJKaBC_matrix[h][ijk][abc];
        }

        // Loop over irreps
        for(int ij_sym = 0; ij_sym < moinfo->get_nirreps(); ij_sym++){
          int m = oo_indexing->get_pairpi(ij_sym);
          int n = o_indexing ->get_pairpi(ij_sym^h);
          int k = n;
          if(m*n*k){
            double alpha = 1.0;
            double beta  = 0.0;
            C_DGEMM_22(m, n, k, alpha,&(t_ijm[ij_sym][0][0]), n, 
                    &(Fockme_matrix[ij_sym^h][0][0]), n, beta, &(H_ijk[ij_sym][0][0]),n);
          }
        }

        for(int ijk = 0;ijk<HiJKaBCMatTmp->get_left_pairpi(h);ijk++){
          int i = ooo_tuples[ooo_offset + ijk][0];
          int j = ooo_tuples[ooo_offset + ijk][1];
          int k = ooo_tuples[ooo_offset + ijk][2];
          int jk = oo_indexing->get_tuple_rel_index(j,k);
          int m  = o_indexing ->get_tuple_rel_index(i);
          int jk_sym = oo_indexing->get_tuple_irrep(j,k);
          HiJKaBCMatTmp->add_six_address_element_abc(i,j,k,abc,-H_ijk[jk_sym][jk][m]);
        }


      } // End loop over abc
      for(int ij_sym = 0; ij_sym < moinfo->get_nirreps(); ij_sym++){
        int m_sym = h^ij_sym;
        release2(H_ijk[ij_sym]);
        release2(t_ijm[ij_sym]);
      }
      release1(H_ijk);
      release1(t_ijm);
    } // End loop over irrep

  } // End loop over references
//   blas->print("t3_test[oOO][vVV]{u}");
//   blas->solve("ERROR{u} = 1000.0 t3_test[oOO][vVV]{u} . t3_test[oOO][vVV]{u}");
//   blas->print("ERROR{u}");
}

/**
 * @brief Computes the contraction
 * \f[ -P(IJ/K) \sum_M t_{IJM}^{ABC} f_{MK} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000)
 */
void CCMRCC::build_t3_IJKABC_amplitudes_diagram3()
{

//    Algorithm
//    Loop over references
//      Loop over irreps of t3_eqns[ijk][abc]
//        Loop over [abc]
//          t3_eqns[ijk][abc] <- ddot t[ij][m](abc) 2@2 Fock[k][m]
//        End loop over [abc]
//      End loop irrep
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);
    string ref_str  = to_string(unique_ref);

    // Grab the temporary matrices
    CCMatTmp  TIJKABCMatTmp = blas->get_MatTmp("t3[OOO][VVV]",unique_ref,none);
    CCMatTmp  HIJKABCMatTmp = blas->get_MatTmp("t3_eqns[OOO][VVV]",unique_ref,none);
    CCMatTmp  FockMEMatTmp  = blas->get_MatTmp("fock[O][O]",unique_ref,none);

    // Grab the indexing for t3[ijk][abc]
    short**   ooo_tuples  = TIJKABCMatTmp->get_left()->get_tuples();

    double*** TIJKABC_matrix = TIJKABCMatTmp->get_matrix();
    double*** FockME_matrix  = FockMEMatTmp->get_matrix();
    CCIndex* oo_indexing = blas->get_index("[oo]");
    CCIndex* o_indexing  = blas->get_index("[o]");

    double ***H_ijk;
    double ***t_ijm;
//  Loop over irrep
    for(int h =0; h < moinfo->get_nirreps();h++){
      allocate1(double**,H_ijk,moinfo->get_nirreps());
      allocate1(double**,t_ijm,moinfo->get_nirreps());
      for(int ij_sym = 0; ij_sym < moinfo->get_nirreps(); ij_sym++){
        int m_sym = h^ij_sym;
        allocate2(double,H_ijk[ij_sym],oo_indexing->get_pairpi(ij_sym),o_indexing ->get_pairpi(m_sym));
        allocate2(double,t_ijm[ij_sym],oo_indexing->get_pairpi(ij_sym),o_indexing ->get_pairpi(m_sym));
      }
      size_t ooo_offset  = TIJKABCMatTmp->get_left()->get_first(h);
      for(int abc = 0;abc<TIJKABCMatTmp->get_right_pairpi(h);abc++){
        // Copy T_ijm(abc) into t_[ij][m]
        for(int ijk = 0;ijk<TIJKABCMatTmp->get_left_pairpi(h);ijk++){
          int i = ooo_tuples[ooo_offset + ijk][0];
          int j = ooo_tuples[ooo_offset + ijk][1];
          int k = ooo_tuples[ooo_offset + ijk][2];
          int ij = oo_indexing->get_tuple_rel_index(i,j);
          int m  = o_indexing ->get_tuple_rel_index(k);
          int ij_sym = oo_indexing->get_tuple_irrep(i,j);
          t_ijm[ij_sym][ij][m] = TIJKABC_matrix[h][ijk][abc];
        }

        // Loop over irreps
        for(int ij_sym = 0; ij_sym < moinfo->get_nirreps(); ij_sym++){
          int m = oo_indexing->get_pairpi(ij_sym);
          int n = o_indexing ->get_pairpi(ij_sym^h);
          int k = n;
          if(m*n*k){
            double alpha = 1.0;
            double beta  = 0.0;
            C_DGEMM_22(m, n, k, alpha,&(t_ijm[ij_sym][0][0]), n, 
                    &(FockME_matrix[ij_sym^h][0][0]), n, beta, &(H_ijk[ij_sym][0][0]),n);
          }
        }

        for(int ijk = 0;ijk<HIJKABCMatTmp->get_left_pairpi(h);ijk++){
          int i = ooo_tuples[ooo_offset + ijk][0];
          int j = ooo_tuples[ooo_offset + ijk][1];
          int k = ooo_tuples[ooo_offset + ijk][2];
          int ij = oo_indexing->get_tuple_rel_index(i,j);
          int m  = o_indexing ->get_tuple_rel_index(k);
          int ij_sym = oo_indexing->get_tuple_irrep(i,j);
          HIJKABCMatTmp->add_six_address_element_Pij_k(i,j,k,abc,-H_ijk[ij_sym][ij][m]);
        }
      } // End loop over abc
      for(int ij_sym = 0; ij_sym < moinfo->get_nirreps(); ij_sym++){
        int m_sym = h^ij_sym;
        release2(H_ijk[ij_sym]);
        release2(t_ijm[ij_sym]);
      }
      release1(H_ijk);
      release1(t_ijm);
    } // End loop over irrep

  } // End loop over references
//   blas->print("t3_test[OOO][VVV]{u}");
//   blas->solve("ERROR{u} = 1000.0 t3_test[OOO][VVV]{u} . t3_test[OOO][VVV]{u}");
//   blas->print("ERROR{u}");
}

/**
 * @brief Computes the contraction
 * \f[ P(ab/c) \sum_e t_{ijk}^{abe} f_{ce} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000)
 */
void CCMRCC::build_t3_ijkabc_amplitudes_diagram4()
{

//    Algorithm
//    Loop over references
//      Loop over irreps of t3_eqns[ijk][abc]
//        Loop over [abc]
//          t3_eqns[ijk][abc] <- ddot t[ij][m](abc) 2@2 Fock[k][m]
//        End loop over [abc]
//      End loop irrep
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);
    string ref_str  = to_string(unique_ref);

    // Grab the temporary matrices
    CCMatTmp  TijkabcMatTmp = blas->get_MatTmp("t3[ooo][vvv]",unique_ref,none);
    CCMatTmp  HijkabcMatTmp = blas->get_MatTmp("t3_eqns[ooo][vvv]",unique_ref,none);
    CCMatTmp  FockaeMatTmp  = blas->get_MatTmp("fock[v][v]",unique_ref,none);

    // Grab the indexing for t3[ijk][abc]
    short**   vvv_tuples  = TijkabcMatTmp->get_right()->get_tuples();

    double*** Tijkabc_matrix = TijkabcMatTmp->get_matrix();
    double*** Fockae_matrix  = FockaeMatTmp->get_matrix();
    CCIndex* vv_indexing = blas->get_index("[vv]");
    CCIndex* v_indexing  = blas->get_index("[v]");

    double ***H_abc;
    double ***t_abc;
//  Loop over irrep
    for(int h =0; h < moinfo->get_nirreps();h++){
      allocate1(double**,H_abc,moinfo->get_nirreps());
      allocate1(double**,t_abc,moinfo->get_nirreps());
      for(int ab_sym = 0; ab_sym < moinfo->get_nirreps(); ab_sym++){
        int e_sym = h^ab_sym;
        allocate2(double,H_abc[ab_sym],vv_indexing->get_pairpi(ab_sym),v_indexing ->get_pairpi(e_sym));
        allocate2(double,t_abc[ab_sym],vv_indexing->get_pairpi(ab_sym),v_indexing ->get_pairpi(e_sym));
      }
      size_t vvv_offset  = TijkabcMatTmp->get_right()->get_first(h);
      for(int ijk = 0;ijk<TijkabcMatTmp->get_left_pairpi(h);ijk++){
        // Copy T_ijk(abc) into t_[ij][k]
        for(int abc = 0;abc<TijkabcMatTmp->get_right_pairpi(h);abc++){
          int a = vvv_tuples[vvv_offset + abc][0];
          int b = vvv_tuples[vvv_offset + abc][1];
          int c = vvv_tuples[vvv_offset + abc][2];
          int ab = vv_indexing->get_tuple_rel_index(a,b);
          int e  = v_indexing ->get_tuple_rel_index(c);
          int ab_sym = vv_indexing->get_tuple_irrep(a,b);
          t_abc[ab_sym][ab][e] = Tijkabc_matrix[h][ijk][abc];
        }

        // Loop over irreps
        for(int ab_sym = 0; ab_sym < moinfo->get_nirreps(); ab_sym++){
          int m = vv_indexing->get_pairpi(ab_sym);
          int n = v_indexing ->get_pairpi(ab_sym^h);
          int k = n;
          if(m*n*k){
            double alpha = 1.0;
            double beta  = 0.0;
            C_DGEMM_22(m, n, k, alpha,&(t_abc[ab_sym][0][0]), n, 
                    &(Fockae_matrix[ab_sym^h][0][0]), n, beta, &(H_abc[ab_sym][0][0]),n);
          }
        }

        for(int abc = 0;abc<TijkabcMatTmp->get_right_pairpi(h);abc++){
          int a = vvv_tuples[vvv_offset + abc][0];
          int b = vvv_tuples[vvv_offset + abc][1];
          int c = vvv_tuples[vvv_offset + abc][2];
          int ab = vv_indexing->get_tuple_rel_index(a,b);
          int e  = v_indexing ->get_tuple_rel_index(c);
          int ab_sym = vv_indexing->get_tuple_irrep(a,b);       
          HijkabcMatTmp->add_six_address_element_Pab_c(ijk,a,b,c,H_abc[ab_sym][ab][e]);
        }
      } // End loop over abc
      for(int ab_sym = 0; ab_sym < moinfo->get_nirreps(); ab_sym++){
        int e_sym = h^ab_sym;
        release2(H_abc[ab_sym]);
        release2(t_abc[ab_sym]);
      }
      release1(H_abc);
      release1(t_abc);
    } // End loop over irrep

  } // End loop over references
//   blas->print("t3_test[ooo][vvv]{u}");
//   blas->solve("ERROR{u} = 1000.0 t3_test[ooo][vvv]{u} . t3_test[ooo][vvv]{u}");
//   blas->print("ERROR{u}");
}

/**
 * @brief Computes the contraction
 * \f[ P(ab/C) \sum_e t_{ijK}^{abe} f_{Ce} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000)
 */
void CCMRCC::build_t3_ijKabC_amplitudes_diagram4()
{
//    Algorithm
//    Loop over references
//      Loop over irreps of t3_eqns[ijk][abc]
//        Loop over [abc]
//          t3_eqns[ijk][abc] <- ddot t[ij][m](abc) 2@2 Fock[k][m]
//        End loop over [abc]
//      End loop irrep
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);
    string ref_str  = to_string(unique_ref);

    // Grab the temporary matrices
    CCMatTmp  TijKabCMatTmp = blas->get_MatTmp("t3[ooO][vvV]",unique_ref,none);
    CCMatTmp  HijKabCMatTmp = blas->get_MatTmp("t3_eqns[ooO][vvV]",unique_ref,none);
    CCMatTmp  FockaeMatTmp  = blas->get_MatTmp("fock[v][v]",unique_ref,none);
    CCMatTmp  FockAEMatTmp  = blas->get_MatTmp("fock[V][V]",unique_ref,none);

    // Grab the indexing for t3[ijk][abc]
    short**   vvv_tuples  = TijKabCMatTmp->get_right()->get_tuples();

    double*** TijKabC_matrix = TijKabCMatTmp->get_matrix();
    double*** Fockae_matrix  = FockaeMatTmp->get_matrix();
    double*** FockAE_matrix  = FockAEMatTmp->get_matrix();
    CCIndex* vv_indexing = blas->get_index("[vv]");
    CCIndex* v_indexing  = blas->get_index("[v]");

    double ***H_abc;
    double ***t_abc;
//  Loop over irrep
    for(int h =0; h < moinfo->get_nirreps();h++){
      allocate1(double**,H_abc,moinfo->get_nirreps());
      allocate1(double**,t_abc,moinfo->get_nirreps());
      for(int ab_sym = 0; ab_sym < moinfo->get_nirreps(); ab_sym++){
        int e_sym = h^ab_sym;
        allocate2(double,H_abc[ab_sym],vv_indexing->get_pairpi(ab_sym),v_indexing ->get_pairpi(e_sym));
        allocate2(double,t_abc[ab_sym],vv_indexing->get_pairpi(ab_sym),v_indexing ->get_pairpi(e_sym));
      }
      size_t vvv_offset  = TijKabCMatTmp->get_right()->get_first(h);
      for(int ijk = 0;ijk<TijKabCMatTmp->get_left_pairpi(h);ijk++){
        // Copy T_ijk(abc) into t_[ij][k]
        for(int abc = 0;abc<TijKabCMatTmp->get_right_pairpi(h);abc++){
          int a = vvv_tuples[vvv_offset + abc][0];
          int b = vvv_tuples[vvv_offset + abc][1];
          int c = vvv_tuples[vvv_offset + abc][2];
          int ab = vv_indexing->get_tuple_rel_index(a,b);
          int e  = v_indexing ->get_tuple_rel_index(c);
          int ab_sym = vv_indexing->get_tuple_irrep(a,b);
          t_abc[ab_sym][ab][e] = TijKabC_matrix[h][ijk][abc];
        }

        // Loop over irreps
        for(int ab_sym = 0; ab_sym < moinfo->get_nirreps(); ab_sym++){
          int m = vv_indexing->get_pairpi(ab_sym);
          int n = v_indexing ->get_pairpi(ab_sym^h);
          int k = n;
          if(m*n*k){
            double alpha = 1.0;
            double beta  = 0.0;
            C_DGEMM_22(m, n, k, alpha,&(t_abc[ab_sym][0][0]), n, 
                    &(FockAE_matrix[ab_sym^h][0][0]), n, beta, &(H_abc[ab_sym][0][0]),n);
          }
        }

        for(int abc = 0;abc<TijKabCMatTmp->get_right_pairpi(h);abc++){
          int a = vvv_tuples[vvv_offset + abc][0];
          int b = vvv_tuples[vvv_offset + abc][1];
          int c = vvv_tuples[vvv_offset + abc][2];
          int ab = vv_indexing->get_tuple_rel_index(a,b);
          int e  = v_indexing ->get_tuple_rel_index(c);
          int ab_sym = vv_indexing->get_tuple_irrep(a,b);       
          HijKabCMatTmp->add_six_address_element_ijk(ijk,a,b,c,H_abc[ab_sym][ab][e]);
        }

        // Copy T_ijk(abc) into t_[ij][k]
        for(int abc = 0;abc<TijKabCMatTmp->get_right_pairpi(h);abc++){
          int a = vvv_tuples[vvv_offset + abc][0];
          int b = vvv_tuples[vvv_offset + abc][1];
          int c = vvv_tuples[vvv_offset + abc][2];
          int ac = vv_indexing->get_tuple_rel_index(a,c);
          int e  = v_indexing ->get_tuple_rel_index(b);
          int ac_sym = vv_indexing->get_tuple_irrep(a,c);
          t_abc[ac_sym][ac][e] = TijKabC_matrix[h][ijk][abc];
        }

        // Loop over irreps
        for(int ab_sym = 0; ab_sym < moinfo->get_nirreps(); ab_sym++){
          int m = vv_indexing->get_pairpi(ab_sym);
          int n = v_indexing ->get_pairpi(ab_sym^h);
          int k = n;
          if(m*n*k){
            double alpha = 1.0;
            double beta  = 0.0;
            C_DGEMM_22(m, n, k, alpha,&(t_abc[ab_sym][0][0]), n, 
                    &(Fockae_matrix[ab_sym^h][0][0]), n, beta, &(H_abc[ab_sym][0][0]),n);
          }
        }

        for(int abc = 0;abc<TijKabCMatTmp->get_right_pairpi(h);abc++){
          int a = vvv_tuples[vvv_offset + abc][0];
          int b = vvv_tuples[vvv_offset + abc][1];
          int c = vvv_tuples[vvv_offset + abc][2];
          int ac = vv_indexing->get_tuple_rel_index(a,c);
          int e  = v_indexing ->get_tuple_rel_index(b);
          int ac_sym = vv_indexing->get_tuple_irrep(a,c);
          HijKabCMatTmp->add_six_address_element_Pab_ijk(ijk,a,b,c,H_abc[ac_sym][ac][e]);
        }

      } // End loop over abc

      for(int ab_sym = 0; ab_sym < moinfo->get_nirreps(); ab_sym++){
        int e_sym = h^ab_sym;
        release2(H_abc[ab_sym]);
        release2(t_abc[ab_sym]);
      }
      release1(H_abc);
      release1(t_abc);
    } // End loop over irrep

  } // End loop over references
//   blas->print("t3_test[ooO][vvV]{u}");
//   blas->solve("ERROR{u} = 1000.0 t3_test[ooO][vvV]{u} . t3_test[ooO][vvV]{u}");
//   blas->print("ERROR{u}");
}

/**
 * @brief Computes the contraction
 * \f[ P(aB/C) \sum_e t_{iJK}^{aBe} f_{Ce} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000)
 */
void CCMRCC::build_t3_iJKaBC_amplitudes_diagram4()
{
//    Algorithm
//    Loop over references
//      Loop over irreps of t3_eqns[ijk][abc]
//        Loop over [abc]
//          t3_eqns[ijk][abc] <- ddot t[ij][m](abc) 2@2 Fock[k][m]
//        End loop over [abc]
//      End loop irrep
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);
    string ref_str  = to_string(unique_ref);

    // Grab the temporary matrices
    CCMatTmp  TiJKaBCMatTmp = blas->get_MatTmp("t3[oOO][vVV]",unique_ref,none);
    CCMatTmp  HiJKaBCMatTmp = blas->get_MatTmp("t3_eqns[oOO][vVV]",unique_ref,none);
    CCMatTmp  FockaeMatTmp  = blas->get_MatTmp("fock[v][v]",unique_ref,none);
    CCMatTmp  FockAEMatTmp  = blas->get_MatTmp("fock[V][V]",unique_ref,none);

    // Grab the indexing for t3[ijk][abc]
    short**   vvv_tuples  = TiJKaBCMatTmp->get_right()->get_tuples();

    double*** TiJKaBC_matrix = TiJKaBCMatTmp->get_matrix();
    double*** Fockae_matrix  = FockaeMatTmp->get_matrix();
    double*** FockAE_matrix  = FockAEMatTmp->get_matrix();
    CCIndex* vv_indexing = blas->get_index("[vv]");
    CCIndex* v_indexing  = blas->get_index("[v]");

    double ***H_abc;
    double ***t_abc;
//  Loop over irrep
    for(int h =0; h < moinfo->get_nirreps();h++){
      allocate1(double**,H_abc,moinfo->get_nirreps());
      allocate1(double**,t_abc,moinfo->get_nirreps());
      for(int ab_sym = 0; ab_sym < moinfo->get_nirreps(); ab_sym++){
        int e_sym = h^ab_sym;
        allocate2(double,H_abc[ab_sym],vv_indexing->get_pairpi(ab_sym),v_indexing ->get_pairpi(e_sym));
        allocate2(double,t_abc[ab_sym],vv_indexing->get_pairpi(ab_sym),v_indexing ->get_pairpi(e_sym));
      }
      size_t vvv_offset  = TiJKaBCMatTmp->get_right()->get_first(h);
      for(int ijk = 0;ijk<TiJKaBCMatTmp->get_left_pairpi(h);ijk++){
        // Copy T_ijk(abc) into t_[ij][k]
        for(int abc = 0;abc<TiJKaBCMatTmp->get_right_pairpi(h);abc++){
          int a = vvv_tuples[vvv_offset + abc][0];
          int b = vvv_tuples[vvv_offset + abc][1];
          int c = vvv_tuples[vvv_offset + abc][2];
          int ab = vv_indexing->get_tuple_rel_index(a,b);
          int e  = v_indexing ->get_tuple_rel_index(c);
          int ab_sym = vv_indexing->get_tuple_irrep(a,b);
          t_abc[ab_sym][ab][e] = TiJKaBC_matrix[h][ijk][abc];
        }

        // Loop over irreps
        for(int ab_sym = 0; ab_sym < moinfo->get_nirreps(); ab_sym++){
          int m = vv_indexing->get_pairpi(ab_sym);
          int n = v_indexing ->get_pairpi(ab_sym^h);
          int k = n;
          if(m*n*k){
            double alpha = 1.0;
            double beta  = 0.0;
            C_DGEMM_22(m, n, k, alpha,&(t_abc[ab_sym][0][0]), n, 
                    &(FockAE_matrix[ab_sym^h][0][0]), n, beta, &(H_abc[ab_sym][0][0]),n);
          }
        }

        for(int abc = 0;abc<TiJKaBCMatTmp->get_right_pairpi(h);abc++){
          int a = vvv_tuples[vvv_offset + abc][0];
          int b = vvv_tuples[vvv_offset + abc][1];
          int c = vvv_tuples[vvv_offset + abc][2];
          int ab = vv_indexing->get_tuple_rel_index(a,b);
          int e  = v_indexing ->get_tuple_rel_index(c);
          int ab_sym = vv_indexing->get_tuple_irrep(a,b);       
          HiJKaBCMatTmp->add_six_address_element_Pbc_ijk(ijk,a,b,c,H_abc[ab_sym][ab][e]);
        }

        // Copy T_ijk(abc) into t_[ij][k]
        for(int abc = 0;abc<TiJKaBCMatTmp->get_right_pairpi(h);abc++){
          int a = vvv_tuples[vvv_offset + abc][0];
          int b = vvv_tuples[vvv_offset + abc][1];
          int c = vvv_tuples[vvv_offset + abc][2];
          int bc = vv_indexing->get_tuple_rel_index(b,c);
          int e  = v_indexing ->get_tuple_rel_index(a);
          int bc_sym = vv_indexing->get_tuple_irrep(b,c);
          t_abc[bc_sym][bc][e] = TiJKaBC_matrix[h][ijk][abc];
        }

        // Loop over irreps
        for(int ab_sym = 0; ab_sym < moinfo->get_nirreps(); ab_sym++){
          int m = vv_indexing->get_pairpi(ab_sym);
          int n = v_indexing ->get_pairpi(ab_sym^h);
          int k = n;
          if(m*n*k){
            double alpha = 1.0;
            double beta  = 0.0;
            C_DGEMM_22(m, n, k, alpha,&(t_abc[ab_sym][0][0]), n, 
                    &(Fockae_matrix[ab_sym^h][0][0]), n, beta, &(H_abc[ab_sym][0][0]),n);
          }
        }

        for(int abc = 0;abc<TiJKaBCMatTmp->get_right_pairpi(h);abc++){
          int a = vvv_tuples[vvv_offset + abc][0];
          int b = vvv_tuples[vvv_offset + abc][1];
          int c = vvv_tuples[vvv_offset + abc][2];
          int bc = vv_indexing->get_tuple_rel_index(b,c);
          int e  = v_indexing ->get_tuple_rel_index(a);
          int bc_sym = vv_indexing->get_tuple_irrep(b,c);
          HiJKaBCMatTmp->add_six_address_element_ijk(ijk,a,b,c,H_abc[bc_sym][bc][e]);
        }

      } // End loop over abc

      for(int ab_sym = 0; ab_sym < moinfo->get_nirreps(); ab_sym++){
        int e_sym = h^ab_sym;
        release2(H_abc[ab_sym]);
        release2(t_abc[ab_sym]);
      }
      release1(H_abc);
      release1(t_abc);
    } // End loop over irrep

  } // End loop over references
//   blas->print("t3_test[oOO][vVV]{u}");
//   blas->solve("ERROR{u} = 1000.0 t3_test[oOO][vVV]{u} . t3_test[oOO][vVV]{u}");
//   blas->print("ERROR{u}");
}


/**
 * @brief Computes the contraction
 * \f[ P(AB/C) \sum_E t_{IJK}^{ABE} f_{CE} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000)
 */
void CCMRCC::build_t3_IJKABC_amplitudes_diagram4()
{

//    Algorithm
//    Loop over references
//      Loop over irreps of t3_eqns[ijk][abc]
//        Loop over [abc]
//          t3_eqns[ijk][abc] <- ddot t[ij][m](abc) 2@2 Fock[k][m]
//        End loop over [abc]
//      End loop irrep
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);
    string ref_str  = to_string(unique_ref);

    // Grab the temporary matrices
    CCMatTmp  TIJKABCMatTmp = blas->get_MatTmp("t3[OOO][VVV]",unique_ref,none);
    CCMatTmp  HIJKABCMatTmp = blas->get_MatTmp("t3_eqns[OOO][VVV]",unique_ref,none);
    CCMatTmp  FockAEMatTmp  = blas->get_MatTmp("fock[V][V]",unique_ref,none);

    // Grab the indexing for t3[ijk][abc]
    short**   vvv_tuples  = TIJKABCMatTmp->get_right()->get_tuples();

    double*** TIJKABC_matrix = TIJKABCMatTmp->get_matrix();
    double*** FockAE_matrix  = FockAEMatTmp->get_matrix();
    CCIndex* vv_indexing = blas->get_index("[vv]");
    CCIndex* v_indexing  = blas->get_index("[v]");

    double ***H_abc;
    double ***t_abc;
//  Loop over irrep
    for(int h =0; h < moinfo->get_nirreps();h++){
      allocate1(double**,H_abc,moinfo->get_nirreps());
      allocate1(double**,t_abc,moinfo->get_nirreps());
      for(int ab_sym = 0; ab_sym < moinfo->get_nirreps(); ab_sym++){
        int e_sym = h^ab_sym;
        allocate2(double,H_abc[ab_sym],vv_indexing->get_pairpi(ab_sym),v_indexing ->get_pairpi(e_sym));
        allocate2(double,t_abc[ab_sym],vv_indexing->get_pairpi(ab_sym),v_indexing ->get_pairpi(e_sym));
      }
      size_t vvv_offset  = TIJKABCMatTmp->get_right()->get_first(h);
      for(int ijk = 0;ijk<TIJKABCMatTmp->get_left_pairpi(h);ijk++){
        // Copy T_ijk(abc) into t_[ij][k]
        for(int abc = 0;abc<TIJKABCMatTmp->get_right_pairpi(h);abc++){
          int a = vvv_tuples[vvv_offset + abc][0];
          int b = vvv_tuples[vvv_offset + abc][1];
          int c = vvv_tuples[vvv_offset + abc][2];
          int ab = vv_indexing->get_tuple_rel_index(a,b);
          int e  = v_indexing ->get_tuple_rel_index(c);
          int ab_sym = vv_indexing->get_tuple_irrep(a,b);
          t_abc[ab_sym][ab][e] = TIJKABC_matrix[h][ijk][abc];
        }

        // Loop over irreps
        for(int ab_sym = 0; ab_sym < moinfo->get_nirreps(); ab_sym++){
          int m = vv_indexing->get_pairpi(ab_sym);
          int n = v_indexing ->get_pairpi(ab_sym^h);
          int k = n;
          if(m*n*k){
            double alpha = 1.0;
            double beta  = 0.0;
            C_DGEMM_22(m, n, k, alpha,&(t_abc[ab_sym][0][0]), n, 
                    &(FockAE_matrix[ab_sym^h][0][0]), n, beta, &(H_abc[ab_sym][0][0]),n);
          }
        }

        for(int abc = 0;abc<TIJKABCMatTmp->get_right_pairpi(h);abc++){
          int a = vvv_tuples[vvv_offset + abc][0];
          int b = vvv_tuples[vvv_offset + abc][1];
          int c = vvv_tuples[vvv_offset + abc][2];
          int ab = vv_indexing->get_tuple_rel_index(a,b);
          int e  = v_indexing ->get_tuple_rel_index(c);
          int ab_sym = vv_indexing->get_tuple_irrep(a,b);       
          HIJKABCMatTmp->add_six_address_element_Pab_c(ijk,a,b,c,H_abc[ab_sym][ab][e]);
        }
      } // End loop over abc
      for(int ab_sym = 0; ab_sym < moinfo->get_nirreps(); ab_sym++){
        int e_sym = h^ab_sym;
        release2(H_abc[ab_sym]);
        release2(t_abc[ab_sym]);
      }
      release1(H_abc);
      release1(t_abc);
    } // End loop over irrep

  } // End loop over references
//   blas->print("t3_test[OOO][VVV]{u}");
//   blas->solve("ERROR{u} = 1000.0 t3_test[OOO][VVV]{u} . t3_test[OOO][VVV]{u}");
//   blas->print("ERROR{u}");
}

}} /* End Namespaces */
