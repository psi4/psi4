/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include <libmoinfo/libmoinfo.h>
#include "mrcc.h"
#include "blas.h"
#include "debugging.h"
#include <libutil/libutil.h>
#include "algebra_interface.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void CCMRCC::build_t3_amplitudes()
{
  if(triples_type>ccsd_t){
    Timer timer;
    DEBUGGING(1,
      fprintf(outfile,"\n\tBuilding the T3 Amplitudes        ...");
      fflush(outfile);
    )
    DEBUGGING(1,
      fprintf(outfile,"\n\tDiagram 1 and 2 ...");
      fflush(outfile);
    );
    Timer timer1;
    build_t3_ijkabc_amplitudes_diagrams12();
    build_t3_ijKabC_amplitudes_diagrams12();
    build_t3_iJKaBC_amplitudes_diagrams12();
    build_t3_IJKABC_amplitudes_diagrams12();
    DEBUGGING(1,
      fprintf(outfile," done. Timing %20.6f s",timer1.get());
      fflush(outfile);
    )

    DEBUGGING(1,
      fprintf(outfile,"\n\tDiagram 3     ...");
      fflush(outfile);
    );
    Timer timer2;
    build_t3_ijkabc_amplitudes_diagram3();
    build_t3_ijKabC_amplitudes_diagram3();
    build_t3_iJKaBC_amplitudes_diagram3();
    build_t3_IJKABC_amplitudes_diagram3();
    DEBUGGING(1,
      fprintf(outfile," done. Timing %20.6f s",timer2.get());
      fflush(outfile);
    )

    DEBUGGING(1,
      fprintf(outfile,"\n\tDiagram 4     ...");
      fflush(outfile);
    )
    Timer timer3;
    build_t3_ijkabc_amplitudes_diagram4();
    build_t3_ijKabC_amplitudes_diagram4();
    build_t3_iJKaBC_amplitudes_diagram4();
    build_t3_IJKABC_amplitudes_diagram4();
    DEBUGGING(1,
      fprintf(outfile," done. Timing %20.6f s",timer3.get());
      fflush(outfile);
    )

    DEBUGGING(1,
      fprintf(outfile," done. Timing %20.6f s",timer.get());
      fflush(outfile);
    )
  }
}

/**
 * @brief Computes the contraction
 * \f[ P(ij/k) P(bc/a) \sum_e t_{ij}^{ae} W'_{bcek} - P(jk/i) P(ab/c) \sum_m t_{im}^{ab} W'_{mcjk} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000). N.B. There is a typo
 * in the original article, in the first equation P(bc/a) was mistyped as P(ab/c).
 */
void CCMRCC::build_t3_ijkabc_amplitudes_diagrams12()
{

  //  Algorithm
  //  Loop over references
  //    Loop over irreps of t3_eqns[ijk][abc]
  //     Loop over [ijk]
  //       unpack [ijk] -> i j k
  //       Loop over [abc]
  //        unpack [abc] -> a b c
  //        form [aij], [bck]
  //        t3_eqns[aij][bck] <- ddot t[aij] and W[bck]
  //      End loop over [abc]
  //    End loop over [ijk]
  //  End loop irrep

  int sym_left,rel_left,sym_right,rel_right,vec_length;
  int incx = 1;
  int incy = 1;
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);
    string ref_str  = to_string(unique_ref);

    // Grab the temporary matrices
    CCMatTmp  MatTmp      = blas->get_MatTmp("t3_eqns[ooo][vvv]",unique_ref,none);
    CCMatTmp  TijabMatTmp = blas->get_MatTmp("t2[oov][v]",unique_ref,none);
    CCMatTmp  WabicMatTmp = blas->get_MatTmp("W'_abic[vvo][v]",unique_ref,none);
    CCMatTmp  TiabjMatTmp = blas->get_MatTmp("t2[ovv][o]",unique_ref,none);
    CCMatTmp  WajkiMatTmp = blas->get_MatTmp("W'_ajki[voo][o]",unique_ref,none);

    // Zero the T3_eqns matrix
    MatTmp->zero_matrix();

    // Grab the indexing for t3[ijk][abc]
    short**   left_tuples  = MatTmp->get_left()->get_tuples();
    short**   right_tuples = MatTmp->get_right()->get_tuples();

    double*** h_matrix     = MatTmp->get_matrix();
    double*** Tijab_matrix = TijabMatTmp->get_matrix();
    double*** Wabic_matrix = WabicMatTmp->get_matrix();
    double*** Tiabj_matrix = TiabjMatTmp->get_matrix();
    double*** Wajki_matrix = WajkiMatTmp->get_matrix();

    // Loop over irrep
    for(int h =0; h < moinfo->get_nirreps();h++){
      size_t left_offset  = MatTmp->get_left()->get_first(h);
      size_t right_offset = MatTmp->get_right()->get_first(h);
      
      for(int ijk = 0;ijk<MatTmp->get_left_pairpi(h);ijk++){
        int i = left_tuples[left_offset + ijk][0];
        int j = left_tuples[left_offset + ijk][1];
        int k = left_tuples[left_offset + ijk][2];
        if((i!=j) && (j!=k) && (i!=k)){
          for(int abc = 0;abc<MatTmp->get_right_pairpi(h);abc++){
            int a = right_tuples[right_offset + abc][0];
            int b = right_tuples[right_offset + abc][1];
            int c = right_tuples[right_offset + abc][2];
            if((a!=b) && (b!=c) && (a!=c)){
              double value=0.0;

              // 1. T[ija][e] W[bck][e]
              sym_left   = TijabMatTmp->get_left()->get_tuple_irrep(i,j,a);
              rel_left   = TijabMatTmp->get_left()->get_tuple_rel_index(i,j,a);
              vec_length = TijabMatTmp->get_right_pairpi(sym_left);
              sym_right  = WabicMatTmp->get_left()->get_tuple_irrep(b,c,k);
              rel_right  = WabicMatTmp->get_left()->get_tuple_rel_index(b,c,k);
              if(vec_length>0){
                value = F_DDOT(&vec_length,&(Tijab_matrix[sym_left][rel_left][0]),&incx,&(Wabic_matrix[sym_right][rel_right][0]),&incy);
                MatTmp->add_six_address_element_Pij_k_Pa_bc(i,j,k,a,b,c,value);
              }
    
              // 2. T[iab][m] W[cjk][m]
              sym_left  = TiabjMatTmp->get_left()->get_tuple_irrep(i,a,b);
              rel_left  = TiabjMatTmp->get_left()->get_tuple_rel_index(i,a,b);
              vec_length= TiabjMatTmp->get_right_pairpi(sym_left);
              sym_right = WajkiMatTmp->get_left()->get_tuple_irrep(c,j,k);
              rel_right = WajkiMatTmp->get_left()->get_tuple_rel_index(c,j,k);
              if(vec_length>0){
                value = F_DDOT(&vec_length,&(Tiabj_matrix[sym_left][rel_left][0]),&incx,&(Wajki_matrix[sym_right][rel_right][0]),&incy);
                MatTmp->add_six_address_element_Pi_jk_Pab_c(i,j,k,a,b,c,-value);
              }

            } // End a!=b!=c test 
          } // End loop over abc
        } // End i!=j!=k test 
      } // End loop over ijk
    } // End loop over irrep

  } // End loop over references
//   blas->print("t3_eqns[ooo][vvv]{u}");
//   blas->solve("ERROR{u} = 1000.0 t3_eqns[ooo][vvv]{u} . t3_eqns[ooo][vvv]{u}");
//   blas->print("ERROR{u}");
}


/**
 * @brief Computes the contraction
 * \f[ P(ij/K) P(bC/a) \sum_e t_{ij}^{ae} W'_{bCeK} - P(jK/i)P(ab/C) \sum_m t_{im}^{ab} W'_{mCjK} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000). N.B. There is a typo
 * in the original article, in the first equation P(bc/a) was mistyped as P(ab/c).
 * This algorithm hence performs
 * \f[ t_{jK}^{eC} W'_{baei} \rightarrow \left \{ \bar{H}_{ijK}^{abC}, -\bar{H}_{jiK}^{abC} \right \} \f]
 * \f[ t_{ij}^{ae} W'_{bCeK} \rightarrow \left \{ \bar{H}_{ijK}^{abC}, -\bar{H}_{ijK}^{baC} \right \} \f]
 * \f[ t_{jK}^{aE} W'_{CbEi} \rightarrow \left \{ - \bar{H}_{ijK}^{abC}, \bar{H}_{ijK}^{baC}, \bar{H}_{jiK}^{abC}, - \bar{H}_{jiK}^{baC} \right \} \f]
 * \f[ t_{im}^{ab} W'_{mCjK} \rightarrow \left \{ -\bar{H}_{ijK}^{abC}, \bar{H}_{jiK}^{abC} \right \} \f]
 * \f[ t_{mK}^{bC} W'_{maji} \rightarrow \left \{ - \bar{H}_{ijK}^{abC}, \bar{H}_{ijK}^{baC} \right \} \f]
 * \f[ t_{jM}^{aC} W'_{MbKi} \rightarrow \left \{ \bar{H}_{ijK}^{abC}, - \bar{H}_{ijK}^{baC}, - \bar{H}_{jiK}^{abC},  \bar{H}_{jiK}^{baC} \right \} \f]
 */
void CCMRCC::build_t3_ijKabC_amplitudes_diagrams12()
{

  //  Algorithm
  //  Loop over references
  //    Loop over irreps of t3_eqns[ijK][abC]
  //     Loop over [ijK]
  //       unpack [ijK] -> i j K
  //       Loop over [abC]
  //        unpack [abC] -> a b C
  //        form [KjC], [bai]
  //        t3_eqns[ijK][abC] <- ddot t[KjC] and W[bai]
  //      End loop over [abC]
  //    End loop over [ijK]
  //  End loop irrep

  int sym_left,rel_left,sym_right,rel_right,vec_length;
  int incx = 1;
  int incy = 1;
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);
    string ref_str  = to_string(unique_ref);

    // Grab the temporary matrices
    CCMatTmp  MatTmp      = blas->get_MatTmp("t3_eqns[ooO][vvV]",unique_ref,none);
    // For diagram 1
    CCMatTmp  TijabMatTmp = blas->get_MatTmp("t2[oov][v]",unique_ref,none);
    CCMatTmp  TJiBaMatTmp = blas->get_MatTmp("t2[OoV][v]",unique_ref,none);
    CCMatTmp  TiJaBMatTmp = blas->get_MatTmp("t2[oOv][V]",unique_ref,none);
    CCMatTmp  WabicMatTmp = blas->get_MatTmp("W'_abic[vvo][v]",unique_ref,none);
    CCMatTmp  WaBIcMatTmp = blas->get_MatTmp("W'_aBIc[vVO][v]",unique_ref,none);
    CCMatTmp  WAbiCMatTmp = blas->get_MatTmp("W'_AbiC[Vvo][V]",unique_ref,none);
    // For diagram 2
    CCMatTmp  TiabjMatTmp = blas->get_MatTmp("t2[ovv][o]",unique_ref,none);
    CCMatTmp  TJaBiMatTmp = blas->get_MatTmp("t2[OvV][o]",unique_ref,none);
    CCMatTmp  TiBaJMatTmp = blas->get_MatTmp("t2[oVv][O]",unique_ref,none);
    CCMatTmp  WajkiMatTmp = blas->get_MatTmp("W'_ajki[voo][o]",unique_ref,none);
    CCMatTmp  WaJkIMatTmp = blas->get_MatTmp("W'_aJkI[vOo][O]",unique_ref,none);
    CCMatTmp  WAjKiMatTmp = blas->get_MatTmp("W'_AjKi[VoO][o]",unique_ref,none);


    // Zero the T3_eqns matrix
    MatTmp->zero_matrix();


    // Grab the indexing for t3[ijK][abC]
    short**   left_tuples  = MatTmp->get_left()->get_tuples();
    short**   right_tuples = MatTmp->get_right()->get_tuples();

    double*** h_matrix     = MatTmp->get_matrix();
    // For diagram 1
    double*** Tijab_matrix = TijabMatTmp->get_matrix();
    double*** TiJaB_matrix = TiJaBMatTmp->get_matrix();
    double*** TJiBa_matrix = TJiBaMatTmp->get_matrix();
    double*** Wabic_matrix = WabicMatTmp->get_matrix();
    double*** WaBIc_matrix = WaBIcMatTmp->get_matrix();
    double*** WAbiC_matrix = WAbiCMatTmp->get_matrix();
    // For diagram 2
    double*** Tiabj_matrix = TiabjMatTmp->get_matrix();
    double*** TJaBi_matrix = TJaBiMatTmp->get_matrix();
    double*** TiBaJ_matrix = TiBaJMatTmp->get_matrix();
    double*** Wajki_matrix = WajkiMatTmp->get_matrix();
    double*** WaJkI_matrix = WaJkIMatTmp->get_matrix();
    double*** WAjKi_matrix = WAjKiMatTmp->get_matrix();

    // Loop over irrep
    for(int h =0; h < moinfo->get_nirreps();h++){
      size_t left_offset  = MatTmp->get_left()->get_first(h);
      size_t right_offset = MatTmp->get_right()->get_first(h);
      
      for(int ijk = 0;ijk<MatTmp->get_left_pairpi(h);ijk++){
        int i = left_tuples[left_offset + ijk][0];
        int j = left_tuples[left_offset + ijk][1];
        int K = left_tuples[left_offset + ijk][2];
        if(i!=j){
          for(int abc = 0;abc<MatTmp->get_right_pairpi(h);abc++){
            int a = right_tuples[right_offset + abc][0];
            int b = right_tuples[right_offset + abc][1];
            int C = right_tuples[right_offset + abc][2];
            if(a!=b){
              double value=0.0;

              // 1. T[KjC][e] W[bai][e]
              sym_left   = TJiBaMatTmp->get_left()->get_tuple_irrep(K,j,C);
              rel_left   = TJiBaMatTmp->get_left()->get_tuple_rel_index(K,j,C);
              vec_length = TJiBaMatTmp->get_right_pairpi(sym_left);
              sym_right  = WabicMatTmp->get_left()->get_tuple_irrep(b,a,i);
              rel_right  = WabicMatTmp->get_left()->get_tuple_rel_index(b,a,i);
              if(vec_length>0){
                value = F_DDOT(&vec_length,&(TJiBa_matrix[sym_left][rel_left][0]),&incx,&(Wabic_matrix[sym_right][rel_right][0]),&incy);
                MatTmp->add_six_address_element_Pij(i,j,K,a,b,C,value);
              }

              // 2. T[ija][e] W[bCK][e]
              sym_left   = TijabMatTmp->get_left()->get_tuple_irrep(i,j,a);
              rel_left   = TijabMatTmp->get_left()->get_tuple_rel_index(i,j,a);
              vec_length = TijabMatTmp->get_right_pairpi(sym_left);
              sym_right  = WaBIcMatTmp->get_left()->get_tuple_irrep(b,C,K);
              rel_right  = WaBIcMatTmp->get_left()->get_tuple_rel_index(b,C,K);
              if(vec_length>0){
                value = F_DDOT(&vec_length,&(Tijab_matrix[sym_left][rel_left][0]),&incx,&(WaBIc_matrix[sym_right][rel_right][0]),&incy);
                MatTmp->add_six_address_element_Pab(i,j,K,a,b,C,value);
              }

              // 3. T[jKa][E] W[Cbi][E]
              sym_left   = TiJaBMatTmp->get_left()->get_tuple_irrep(j,K,a);
              rel_left   = TiJaBMatTmp->get_left()->get_tuple_rel_index(j,K,a);
              vec_length = TiJaBMatTmp->get_right_pairpi(sym_left);
              sym_right  = WAbiCMatTmp->get_left()->get_tuple_irrep(C,b,i);
              rel_right  = WAbiCMatTmp->get_left()->get_tuple_rel_index(C,b,i);
              if(vec_length>0){
                value = F_DDOT(&vec_length,&(TiJaB_matrix[sym_left][rel_left][0]),&incx,&(WAbiC_matrix[sym_right][rel_right][0]),&incy);
                MatTmp->add_six_address_element_Pij_Pab(i,j,K,a,b,C,-value);
              }

              // 4. T[iab][m] W[CjK][m]
              sym_left  = TiabjMatTmp->get_left()->get_tuple_irrep(i,a,b);
              rel_left  = TiabjMatTmp->get_left()->get_tuple_rel_index(i,a,b);
              vec_length= TiabjMatTmp->get_right_pairpi(sym_left);
              sym_right = WAjKiMatTmp->get_left()->get_tuple_irrep(C,j,K);
              rel_right = WAjKiMatTmp->get_left()->get_tuple_rel_index(C,j,K);
              if(vec_length>0){
                value = F_DDOT(&vec_length,&(Tiabj_matrix[sym_left][rel_left][0]),&incx,&(WAjKi_matrix[sym_right][rel_right][0]),&incy);
                MatTmp->add_six_address_element_Pij(i,j,K,a,b,C,-value);
              }

              // 5. T[KbC][m] W[aji][m]
              sym_left  = TJaBiMatTmp->get_left()->get_tuple_irrep(K,b,C);
              rel_left  = TJaBiMatTmp->get_left()->get_tuple_rel_index(K,b,C);
              vec_length= TJaBiMatTmp->get_right_pairpi(sym_left);
              sym_right = WajkiMatTmp->get_left()->get_tuple_irrep(a,j,i);
              rel_right = WajkiMatTmp->get_left()->get_tuple_rel_index(a,j,i);
              if(vec_length>0){
                value = F_DDOT(&vec_length,&(TJaBi_matrix[sym_left][rel_left][0]),&incx,&(Wajki_matrix[sym_right][rel_right][0]),&incy);
                MatTmp->add_six_address_element_Pab(i,j,K,a,b,C,-value);
              }

              // 6. T[jCa][M] W[bKi][M]
              sym_left  = TiBaJMatTmp->get_left()->get_tuple_irrep(j,C,a);
              rel_left  = TiBaJMatTmp->get_left()->get_tuple_rel_index(j,C,a);
              vec_length= TiBaJMatTmp->get_right_pairpi(sym_left);
              sym_right = WaJkIMatTmp->get_left()->get_tuple_irrep(b,K,i);
              rel_right = WaJkIMatTmp->get_left()->get_tuple_rel_index(b,K,i);
              if(vec_length>0){
                value = F_DDOT(&vec_length,&(TiBaJ_matrix[sym_left][rel_left][0]),&incx,&(WaJkI_matrix[sym_right][rel_right][0]),&incy);
                MatTmp->add_six_address_element_Pij_Pab(i,j,K,a,b,C,value);
              }

            } // End a!=b!=c test 
          } // End loop over abc
        } // End i!=j!=k test 
      } // End loop over ijk
    } // End loop over irrep

  } // End loop over references
//   blas->print("t3_eqns[ooO][vvV]{u}");
//   blas->solve("ERROR{u} = 1000.0 t3_eqns[ooO][vvV]{u} . t3_eqns[ooO][vvV]{u}");
//   blas->print("ERROR{u}");
}


/**
 * @brief Computes the contraction
 * \f[ P(iJ/K) P(BC/a) \sum_e t_{iJ}^{ae} W'_{BCeK} - P(JK/i)P(aB/C) \sum_m t_{im}^{aB} W'_{mCJK} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000). N.B. There is a typo
 * in the original article, in the first equation P(bc/a) was mistyped as P(ab/c).
 * This algorithm hence performs
 * \f[ t_{iJ}^{aE} W'_{BCEK} \rightarrow \left \{ \bar{H}_{iJK}^{aBC}, -\bar{H}_{iKJ}^{aBC} \right \} \f]
 * \f[ t_{KJ}^{CE} W'_{BaEi} \rightarrow \left \{ \bar{H}_{iJK}^{aBC}, -\bar{H}_{iJK}^{aCB} \right \} \f]
 * \f[ t_{iJ}^{eB} W'_{aCeK} \rightarrow \left \{ \bar{H}_{iJK}^{aBC}, -\bar{H}_{iJK}^{aCB}, -\bar{H}_{iKJ}^{aBC},  \bar{H}_{iKJ}^{aCB} \right \} \f]
 * \f[ t_{KM}^{CB} W'_{MaJi} \rightarrow \left \{ -\bar{H}_{iJK}^{aBC}, \bar{H}_{iKJ}^{aBC} \right \} \f]
 * \f[ t_{iM}^{aB} W'_{MCJK} \rightarrow \left \{ -\bar{H}_{iJK}^{aBC}, \bar{H}_{iJK}^{aCB} \right \} \f]
 * \f[ t_{mK}^{aB} W'_{mCiJ} \rightarrow \left \{ \bar{H}_{iJK}^{aBC}, -\bar{H}_{iJK}^{aCB}, -\bar{H}_{iKJ}^{aBC}, \bar{H}_{iKJ}^{aCB} \right \} \f]
 */
void CCMRCC::build_t3_iJKaBC_amplitudes_diagrams12()
{

  //  Algorithm
  //  Loop over references
  //    Loop over irreps of t3_eqns[iJK][aBC]
  //     Loop over [iJK]
  //       unpack [iJK] -> i J K
  //       Loop over [aBC]
  //        unpack [aBC] -> a B C
  //        form [KjC], [bai]
  //        t3_eqns[iJK][aBC] <- ddot t[KJC] and W[Bai]
  //      End loop over [aBC]
  //    End loop over [iJK]
  //  End loop irrep

  int sym_left,rel_left,sym_right,rel_right,vec_length;
  int incx = 1;
  int incy = 1;
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);
    string ref_str  = to_string(unique_ref);

    // Grab the temporary matrices
    CCMatTmp  MatTmp      = blas->get_MatTmp("t3_eqns[oOO][vVV]",unique_ref,none);
    // For diagram 1
    CCMatTmp  TiJaBMatTmp = blas->get_MatTmp("t2[oOv][V]",unique_ref,none);
    CCMatTmp  TIJABMatTmp = blas->get_MatTmp("t2[OOV][V]",unique_ref,none);
    CCMatTmp  TJiBaMatTmp = blas->get_MatTmp("t2[OoV][v]",unique_ref,none);
    CCMatTmp  WaBIcMatTmp = blas->get_MatTmp("W'_aBIc[vVO][v]",unique_ref,none);
    CCMatTmp  WAbiCMatTmp = blas->get_MatTmp("W'_AbiC[Vvo][V]",unique_ref,none);
    CCMatTmp  WABICMatTmp = blas->get_MatTmp("W'_ABIC[VVO][V]",unique_ref,none); 
    // For diagram 2
    CCMatTmp  TiBaJMatTmp = blas->get_MatTmp("t2[oVv][O]",unique_ref,none);
    CCMatTmp  TJaBiMatTmp = blas->get_MatTmp("t2[OvV][o]",unique_ref,none);
    CCMatTmp  TIABJMatTmp = blas->get_MatTmp("t2[OVV][O]",unique_ref,none);
    CCMatTmp  WaJkIMatTmp = blas->get_MatTmp("W'_aJkI[vOo][O]",unique_ref,none);
    CCMatTmp  WAjKiMatTmp = blas->get_MatTmp("W'_AjKi[VoO][o]",unique_ref,none);
    CCMatTmp  WAJKIMatTmp = blas->get_MatTmp("W'_AJKI[VOO][O]",unique_ref,none);


    // Zero the T3_eqns matrix
    MatTmp->zero_matrix();


    // Grab the indexing for t3[ijk][abc]
    short**   left_tuples  = MatTmp->get_left()->get_tuples();
    short**   right_tuples = MatTmp->get_right()->get_tuples();

    double*** h_matrix     = MatTmp->get_matrix();
    // For diagram 1
    double*** TJiBa_matrix = TJiBaMatTmp->get_matrix();
    double*** TIJAB_matrix = TIJABMatTmp->get_matrix();
    double*** TiJaB_matrix = TiJaBMatTmp->get_matrix();
    double*** WaBIc_matrix = WaBIcMatTmp->get_matrix();
    double*** WAbiC_matrix = WAbiCMatTmp->get_matrix();
    double*** WABIC_matrix = WABICMatTmp->get_matrix();
    // For diagram 2
    double*** TiBaJ_matrix = TiBaJMatTmp->get_matrix();
    double*** TJaBi_matrix = TJaBiMatTmp->get_matrix();
    double*** TIABJ_matrix = TIABJMatTmp->get_matrix();
    double*** WaJkI_matrix = WaJkIMatTmp->get_matrix();
    double*** WAjKi_matrix = WAjKiMatTmp->get_matrix();
    double*** WAJKI_matrix = WAJKIMatTmp->get_matrix();

    // Loop over irrep
    for(int h =0; h < moinfo->get_nirreps();h++){
      size_t left_offset  = MatTmp->get_left()->get_first(h);
      size_t right_offset = MatTmp->get_right()->get_first(h);
      
      for(int ijk = 0;ijk<MatTmp->get_left_pairpi(h);ijk++){
        int i = left_tuples[left_offset + ijk][0];
        int J = left_tuples[left_offset + ijk][1];
        int K = left_tuples[left_offset + ijk][2];
        if(J!=K){
          for(int abc = 0;abc<MatTmp->get_right_pairpi(h);abc++){
            int a = right_tuples[right_offset + abc][0];
            int B = right_tuples[right_offset + abc][1];
            int C = right_tuples[right_offset + abc][2];
            if(B!=C){
              double value=0.0;

              // 1. T[iJa][E] W[BCK][E]
              sym_left   = TiJaBMatTmp->get_left()->get_tuple_irrep(i,J,a);
              rel_left   = TiJaBMatTmp->get_left()->get_tuple_rel_index(i,J,a);
              vec_length = TiJaBMatTmp->get_right_pairpi(sym_left);
              sym_right  = WABICMatTmp->get_left()->get_tuple_irrep(B,C,K);
              rel_right  = WABICMatTmp->get_left()->get_tuple_rel_index(B,C,K);
              if(vec_length>0){
                value = F_DDOT(&vec_length,&(TiJaB_matrix[sym_left][rel_left][0]),&incx,&(WABIC_matrix[sym_right][rel_right][0]),&incy);
                MatTmp->add_six_address_element_Pjk(i,J,K,a,B,C,value);
              }

               // 2. T[KJC][E] W[Bai][E]
              sym_left   = TIJABMatTmp->get_left()->get_tuple_irrep(K,J,C);
              rel_left   = TIJABMatTmp->get_left()->get_tuple_rel_index(K,J,C);
              vec_length = TIJABMatTmp->get_right_pairpi(sym_left);
              sym_right  = WAbiCMatTmp->get_left()->get_tuple_irrep(B,a,i);
              rel_right  = WAbiCMatTmp->get_left()->get_tuple_rel_index(B,a,i);
              if(vec_length>0){
                value = F_DDOT(&vec_length,&(TIJAB_matrix[sym_left][rel_left][0]),&incx,&(WAbiC_matrix[sym_right][rel_right][0]),&incy);
                MatTmp->add_six_address_element_Pbc(i,J,K,a,B,C,value);
              }

              // 3. T[JiB][e] W[aCK][e]
              sym_left   = TJiBaMatTmp->get_left()->get_tuple_irrep(J,i,B);
              rel_left   = TJiBaMatTmp->get_left()->get_tuple_rel_index(J,i,B);
              vec_length = TJiBaMatTmp->get_right_pairpi(sym_left);
              sym_right  = WaBIcMatTmp->get_left()->get_tuple_irrep(a,C,K);
              rel_right  = WaBIcMatTmp->get_left()->get_tuple_rel_index(a,C,K);
              if(vec_length>0){
                value = F_DDOT(&vec_length,&(TJiBa_matrix[sym_left][rel_left][0]),&incx,&(WaBIc_matrix[sym_right][rel_right][0]),&incy);
                MatTmp->add_six_address_element_Pjk_Pbc(i,J,K,a,B,C,value);
              }

              // 4. T[KCB][M] W[aJi][M]
              sym_left   = TIABJMatTmp->get_left()->get_tuple_irrep(K,C,B);
              rel_left   = TIABJMatTmp->get_left()->get_tuple_rel_index(K,C,B);
              vec_length = TIABJMatTmp->get_right_pairpi(sym_left);
              sym_right  = WaJkIMatTmp->get_left()->get_tuple_irrep(a,J,i);
              rel_right  = WaJkIMatTmp->get_left()->get_tuple_rel_index(a,J,i);
              if(vec_length>0){
                value = F_DDOT(&vec_length,&(TIABJ_matrix[sym_left][rel_left][0]),&incx,&(WaJkI_matrix[sym_right][rel_right][0]),&incy);
                MatTmp->add_six_address_element_Pjk(i,J,K,a,B,C,-value);
              }

              // 5. T[iBa][M] W[CJK][M]
              sym_left   = TiBaJMatTmp->get_left()->get_tuple_irrep(i,B,a);
              rel_left   = TiBaJMatTmp->get_left()->get_tuple_rel_index(i,B,a);
              vec_length = TiBaJMatTmp->get_right_pairpi(sym_left);
              sym_right  = WAJKIMatTmp->get_left()->get_tuple_irrep(C,J,K);
              rel_right  = WAJKIMatTmp->get_left()->get_tuple_rel_index(C,J,K);
              if(vec_length>0){
                value = F_DDOT(&vec_length,&(TiBaJ_matrix[sym_left][rel_left][0]),&incx,&(WAJKI_matrix[sym_right][rel_right][0]),&incy);
                MatTmp->add_six_address_element_Pbc(i,J,K,a,B,C,-value);
              }

              // 6. T[KaB][M] W[CJK][M]
              sym_left   = TJaBiMatTmp->get_left()->get_tuple_irrep(K,a,B);
              rel_left   = TJaBiMatTmp->get_left()->get_tuple_rel_index(K,a,B);
              vec_length = TJaBiMatTmp->get_right_pairpi(sym_left);
              sym_right  = WAjKiMatTmp->get_left()->get_tuple_irrep(C,i,J);
              rel_right  = WAjKiMatTmp->get_left()->get_tuple_rel_index(C,i,J);
              if(vec_length>0){
                value = F_DDOT(&vec_length,&(TJaBi_matrix[sym_left][rel_left][0]),&incx,&(WAjKi_matrix[sym_right][rel_right][0]),&incy);
                MatTmp->add_six_address_element_Pjk_Pbc(i,J,K,a,B,C,value);
              }
            } // End a!=b!=c test 
          } // End loop over abc
        } // End i!=j!=k test 
      } // End loop over ijk
    } // End loop over irrep

  } // End loop over references
//   blas->print("t3_eqns[oOO][vVV]{u}");
//   blas->solve("ERROR{u} = 1000.0 t3_eqns[oOO][vVV]{u} . t3_eqns[oOO][vVV]{u}");
//   blas->print("ERROR{u}");
}


/**
 * @brief Computes the contraction
 * \f[ P(IJ/K) P(BC/A) \sum_E t_{IJ}^{AE} W'_{BCEK} - P(JK/I) P(AB/C) \sum_m t_{IM}^{AB} W'_{MCJK} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000). N.B. There is a typo
 * in the original article, in the first equation P(BC/A) was mistyped as P(AB/C).
 */
void CCMRCC::build_t3_IJKABC_amplitudes_diagrams12()
{

  //  Algorithm
  //  Loop over references
  //    Loop over irreps of t3_eqns[IJK][ABC]
  //     Loop over [IJK]
  //       unpack [IJK] -> I J K
  //       Loop over [ABC]
  //        unpack [ABC] -> A B C
  //        form [AIJ], [BCK]
  //        t3_eqns[AIJ][BCK] <- ddot t[AIJ] and W[BCK]
  //      End loop over [ABC]
  //    End loop over [IJK]
  //  End loop irrep

  int sym_left,rel_left,sym_right,rel_right,vec_length;
  int incx = 1;
  int incy = 1;
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);
    string ref_str  = to_string(unique_ref);

    // Grab the temporary matrices
    CCMatTmp  MatTmp      = blas->get_MatTmp("t3_eqns[OOO][VVV]",unique_ref,none);
    CCMatTmp  TIJABMatTmp = blas->get_MatTmp("t2[OOV][V]",unique_ref,none);
    CCMatTmp  WABICMatTmp = blas->get_MatTmp("W'_ABIC[VVO][V]",unique_ref,none);
    CCMatTmp  TIABJMatTmp = blas->get_MatTmp("t2[OVV][O]",unique_ref,none);
    CCMatTmp  WAJKIMatTmp = blas->get_MatTmp("W'_AJKI[VOO][O]",unique_ref,none);


    // Zero the T3_eqns matrix
    MatTmp->zero_matrix();


    // Grab the indexing for t3[ijk][abc]
    short**   left_tuples  = MatTmp->get_left()->get_tuples();
    short**   right_tuples = MatTmp->get_right()->get_tuples();

    double*** h_matrix     = MatTmp->get_matrix();
    double*** TIJAB_matrix = TIJABMatTmp->get_matrix();
    double*** WABIC_matrix = WABICMatTmp->get_matrix();
    double*** TIABJ_matrix = TIABJMatTmp->get_matrix();
    double*** WAJKI_matrix = WAJKIMatTmp->get_matrix();

    // Loop over irrep
    for(int h =0; h < moinfo->get_nirreps();h++){
      size_t left_offset  = MatTmp->get_left()->get_first(h);
      size_t right_offset = MatTmp->get_right()->get_first(h);
      
      for(int ijk = 0;ijk<MatTmp->get_left_pairpi(h);ijk++){
        int I = left_tuples[left_offset + ijk][0];
        int J = left_tuples[left_offset + ijk][1];
        int K = left_tuples[left_offset + ijk][2];
        if((I!=J) && (J!=K) && (I!=K)){
          for(int abc = 0;abc<MatTmp->get_right_pairpi(h);abc++){
            int A = right_tuples[right_offset + abc][0];
            int B = right_tuples[right_offset + abc][1];
            int C = right_tuples[right_offset + abc][2];
            if((A!=B) && (B!=C) && (A!=C)){
              double value=0.0;

              // 1. T[IJA][E] W[BCK][E]
              sym_left   = TIJABMatTmp->get_left()->get_tuple_irrep(I,J,A);
              rel_left   = TIJABMatTmp->get_left()->get_tuple_rel_index(I,J,A);
              vec_length = TIJABMatTmp->get_right_pairpi(sym_left);
              sym_right  = WABICMatTmp->get_left()->get_tuple_irrep(B,C,K);
              rel_right  = WABICMatTmp->get_left()->get_tuple_rel_index(B,C,K);
              if(vec_length>0){
                value = F_DDOT(&vec_length,&(TIJAB_matrix[sym_left][rel_left][0]),&incx,&(WABIC_matrix[sym_right][rel_right][0]),&incy);
                MatTmp->add_six_address_element_Pij_k_Pa_bc(I,J,K,A,B,C,value);
              }

              // 2. T[IAB][M] W[CJK][M]
              sym_left  = TIABJMatTmp->get_left()->get_tuple_irrep(I,A,B);
              rel_left  = TIABJMatTmp->get_left()->get_tuple_rel_index(I,A,B);
              vec_length= TIABJMatTmp->get_right_pairpi(sym_left);
              sym_right = WAJKIMatTmp->get_left()->get_tuple_irrep(C,J,K);
              rel_right = WAJKIMatTmp->get_left()->get_tuple_rel_index(C,J,K);
              if(vec_length>0){
                value = F_DDOT(&vec_length,&(TIABJ_matrix[sym_left][rel_left][0]),&incx,&(WAJKI_matrix[sym_right][rel_right][0]),&incy);
                MatTmp->add_six_address_element_Pi_jk_Pab_c(I,J,K,A,B,C,-value);
              }

            } // End a!=b!=c test 
          } // End loop over abc
        } // End i!=j!=k test 
      } // End loop over ijk
    } // End loop over irrep

  } // End loop over references
//   blas->print("t3_eqns[OOO][VVV]{u}");
//   blas->solve("ERROR{u} = 1000.0 t3_eqns[OOO][VVV]{u} . t3_eqns[OOO][VVV]{u}");
//   blas->print("ERROR{u}");
}




}} /* End Namespaces */
