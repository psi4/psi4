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
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include "blas.h"
#include "index.h"
#include "matrix.h"
#include "mrcc.h"
#include "debugging.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;

using namespace std;

void CCMRCC::build_t1_amplitudes()
{
  build_t1_ia_amplitudes();
  build_t1_IA_amplitudes();
}

void CCMRCC::build_t1_ia_amplitudes()
{
  Timer timer;
  DEBUGGING(1,
    outfile->Printf("\n\tBuilding the t1_ia Amplitudes     ...");

  )
  // Closed-shell
  blas->append("t1_eqns[o][v]{c} = fock[o][v]{c}");
  blas->append("t1_eqns[o][v]{c} += t1[o][v]{c} 2@2 F_ae[v][v]{c}");
  blas->append("t1_eqns[o][v]{c} += - F_mi[o][o]{c} 1@1 t1[o][v]{c}");
  blas->append("t1_eqns[o][v]{c} += #12# t2[ov][ov]{c} 2@1 F_me[ov]{c}");
  blas->append("t1_eqns[o][v]{c} += #12# t2[ov][OV]{c} 2@1 F_me[ov]{c}");

  blas->append("t1_eqns[o][v]{c} += #12# - <[ov]|[ov]> 2@1 t1[ov]{c}");
  blas->append("t1_eqns[o][v]{c} += #21# 2 ([ov]|[vo]) 1@1 t1[ov]{c}");

  blas->append("t1_eqns[o][v]{c} += 1/2 t2[o][ovv]{c} 2@2 <[v]:[ovv]>");
  blas->append("t1_eqns[o][v]{c} +=     t2[o][OvV]{c} 2@2 <[v]|[ovv]>");

  blas->append("t1_eqns[o][v]{c} += -1/2 <[o]:[voo]> 2@2 t2[v][voo]{c}");
  blas->append("t1_eqns[o][v]{c} += - <[o]|[voo]> 2@2 t2[v][VoO]{c}");

  // Open-shell
  blas->append("t1_eqns[o][v]{o} = fock[o][v]{o}");
  blas->append("t1_eqns[o][v]{o} += t1[o][v]{o} 2@2 F_ae[v][v]{o}");
  blas->append("t1_eqns[o][v]{o} += - F_mi[o][o]{o} 1@1 t1[o][v]{o}");
  blas->append("t1_eqns[o][v]{o} += #12# t2[ov][ov]{o} 2@1 F_me[ov]{o}");
  blas->append("t1_eqns[o][v]{o} += #12# t2[ov][OV]{o} 2@1 F_ME[OV]{o}");

  blas->append("t1_eqns[o][v]{o} += #12# - <[ov]|[ov]> 2@1 t1[ov]{o}");
  blas->append("t1_eqns[o][v]{o} += #21#  ([ov]|[vo]) 1@1 t1[ov]{o}");
  blas->append("t1_eqns[o][v]{o} += #21#  ([ov]|[vo]) 1@1 t1[OV]{o}");

  blas->append("t1_eqns[o][v]{o} += 1/2 t2[o][ovv]{o} 2@2 <[v]:[ovv]>");
  blas->append("t1_eqns[o][v]{o} +=     t2[o][OvV]{o} 2@2 <[v]|[ovv]>");

  blas->append("t1_eqns[o][v]{o} += -1/2 <[o]:[voo]> 2@2 t2[v][voo]{o}");
  blas->append("t1_eqns[o][v]{o} += - <[o]|[voo]> 2@2 t2[v][VoO]{o}");


  if(pert_cbs && pert_cbs_coupling){
    outfile->Printf("\n Computing frozen-virtual contribution to H(ia)");
    blas->append("t1_eqns[o][v]{u} +=     t2_1[o][ovf]{u} 2@2 <[v]:[ovf]>");
    blas->append("t1_eqns[o][v]{u} +=     t2_1[o][OvF]{u} 2@2 <[v]|[ovf]>");
    blas->append("t1_eqns[o][v]{u} +=     t2_1[o][OfV]{u} 2@2 <[v]|[ofv]>");

    blas->append("t1_eqns[o][v]{u} += 1/2 t2_1[o][off]{u} 2@2 <[v]:[off]>");
    blas->append("t1_eqns[o][v]{u} +=     t2_1[o][OfF]{u} 2@2 <[v]|[off]>");

    blas->append("t1_eqns[o][v]{u} += -1/2 <[o]:[foo]> 2@2 t2_1[v][foo]{u}");
    blas->append("t1_eqns[o][v]{u} += -    <[o]|[foo]> 2@2 t2_1[v][FoO]{u}");
  }

  DEBUGGING(3,blas->print("t1_eqns[o][v]{u}"););


  DEBUGGING(1,
    outfile->Printf(" done. Timing %20.6f s",timer.get());

  );
}

void CCMRCC::build_t1_IA_amplitudes()
{
  Timer timer;
  DEBUGGING(1,
    outfile->Printf("\n\tBuilding the t1_IA Amplitudes     ...");

  );
  // Closed-shell
  blas->append("t1_eqns[O][V]{c} = t1_eqns[o][v]{c}");

  // Open-shell
  blas->append("t1_eqns[O][V]{o} = fock[O][V]{o}");
  blas->append("t1_eqns[O][V]{o} += t1[O][V]{o} 2@2 F_AE[V][V]{o}");
  blas->append("t1_eqns[O][V]{o} += - F_MI[O][O]{o} 1@1 t1[O][V]{o}");
  blas->append("t1_eqns[O][V]{o} += #12# t2[OV][OV]{o} 2@1 F_ME[OV]{o}");
  blas->append("t1_eqns[O][V]{o} += #12# t2[ov][OV]{o} 1@1 F_me[ov]{o}");

  blas->append("t1_eqns[O][V]{o} += #12# - <[ov]|[ov]> 2@1 t1[OV]{o}");
  blas->append("t1_eqns[O][V]{o} += #21#  ([ov]|[vo]) 1@1 t1[OV]{o}");
  blas->append("t1_eqns[O][V]{o} += #21#  ([ov]|[vo]) 1@1 t1[ov]{o}");

  blas->append("t1_eqns[O][V]{o} += 1/2 t2[O][OVV]{o} 2@2 <[v]:[ovv]>");
  blas->append("t1_eqns[O][V]{o} +=     t2[O][oVv]{o} 2@2 <[v]|[ovv]>");

  blas->append("t1_eqns[O][V]{o} += -1/2 <[o]:[voo]> 2@2 t2[V][VOO]{o}");
  blas->append("t1_eqns[O][V]{o} += - <[o]|[voo]> 2@2 t2[V][vOo]{o}");

  DEBUGGING(3,blas->print("t1_eqns[O][V]{u}"););

  DEBUGGING(1,
    outfile->Printf(" done. Timing %20.6f s",timer.get());

  );
}


void CCMRCC::build_t1_amplitudes_triples()
{
  build_t1_ia_amplitudes_triples();
  build_t1_IA_amplitudes_triples();
}

/**
 * @brief Computes the contraction
 * \f[ \frac{1}{4} \sum_{mnef} t_{imn}^{aef} <mn||ef> \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000).
 * This algorithm hence performs
 * \f[ \frac{1}{4} \sum_{mnef} t_{imn}^{aef} <mn||ef> + \sum_{mNeF} t_{imN}^{aeF} <mN|eF> + \frac{1}{4} \sum_{MNEF} t_{iMN}^{aEF} <MN||EF>\rightarrow  \{ \bar{H}_{i}^{a} \} \f]
 */
void CCMRCC::build_t1_ia_amplitudes_triples()
{
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);

    // Grab the temporary matrices
    CCMatTmp  HiaMatTmp     = blas->get_MatTmp("t1_eqns[o][v]",unique_ref,none);
    CCMatTmp  TijkabcMatTmp = blas->get_MatTmp("t3[ooo][vvv]",unique_ref,none);
    CCMatTmp  TijKabCMatTmp = blas->get_MatTmp("t3[ooO][vvV]",unique_ref,none);
    CCMatTmp  TiJKaBCMatTmp = blas->get_MatTmp("t3[oOO][vVV]",unique_ref,none);
    CCMatTmp  ImnefMatTmp   = blas->get_MatTmp("<[oo]:[vv]>",none);
    CCMatTmp  ImNeFMatTmp   = blas->get_MatTmp("<[oo]|[vv]>",none);

    // Grab the indexing for t3[ijk][abc]
    short**   i_tuples = HiaMatTmp->get_left()->get_tuples();
    short**   a_tuples = HiaMatTmp->get_right()->get_tuples();
    short**   mn_tuples = ImnefMatTmp->get_left()->get_tuples();
    short**   ef_tuples = ImnefMatTmp->get_right()->get_tuples();

    double*** Tijkabc_matrix = TijkabcMatTmp->get_matrix();
    double*** TijKabC_matrix = TijKabCMatTmp->get_matrix();
    double*** TiJKaBC_matrix = TiJKaBCMatTmp->get_matrix();
    double*** Hia_matrix     = HiaMatTmp->get_matrix();
    double*** Imnef_matrix   = ImnefMatTmp->get_matrix();
    double*** ImNeF_matrix   = ImNeFMatTmp->get_matrix();
    CCIndex*  ijkIndex       = blas->get_index("[ooo]");
    CCIndex*  abcIndex       = blas->get_index("[vvv]");

    for(int h =0; h < moinfo->get_nirreps();h++){
      size_t i_offset  = HiaMatTmp->get_left()->get_first(h);
      size_t a_offset = HiaMatTmp->get_right()->get_first(h);
      for(int a = 0;a <HiaMatTmp->get_right_pairpi(h);a++){
        int a_abs = a + a_offset;
        for(int i = 0;i<HiaMatTmp->get_left_pairpi(h);i++){
          int i_abs = i + i_offset;
          // <mn||ef> contribution
          for(int mn_sym =0; mn_sym < moinfo->get_nirreps();mn_sym++){
            size_t mn_offset = ImnefMatTmp->get_left()->get_first(mn_sym);
            size_t ef_offset = ImnefMatTmp->get_right()->get_first(mn_sym);
            for(int ef = 0;ef <ImnefMatTmp->get_right_pairpi(mn_sym);ef++){
              int e = ef_tuples[ef_offset + ef][0];
              int f = ef_tuples[ef_offset + ef][1];
              size_t aef  = abcIndex->get_tuple_rel_index(a_abs,e,f);
              int aef_sym = abcIndex->get_tuple_irrep(a_abs,e,f);
              for(int mn = 0;mn <ImnefMatTmp->get_left_pairpi(mn_sym);mn++){
                int m = mn_tuples[mn_offset + mn][0];
                int n = mn_tuples[mn_offset + mn][1];
                size_t imn  = ijkIndex->get_tuple_rel_index(i_abs,m,n);
                Hia_matrix[h][i][a] += 0.25 * Tijkabc_matrix[aef_sym][imn][aef] * Imnef_matrix[mn_sym][mn][ef];
                Hia_matrix[h][i][a] += 0.25 * TiJKaBC_matrix[aef_sym][imn][aef] * Imnef_matrix[mn_sym][mn][ef];
                Hia_matrix[h][i][a] +=        TijKabC_matrix[aef_sym][imn][aef] * ImNeF_matrix[mn_sym][mn][ef];
              }
            }
          }
        }
      }
    }
  }
//   blas->print("t1_eqns[o][v]{u}");
//   blas->solve("ERROR{u} = 100000000000.0 t1_eqns[o][v]{u} . t1_eqns[o][v]{u}");
//   blas->print("ERROR{u}");
}

/**
 * @brief Computes the contraction
 * \f[ \frac{1}{4} \sum_{mnef} t_{imn}^{aef} <mn||ef> \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000).
 * This algorithm hence performs
 * \f[ \frac{1}{4} \sum_{mnef} t_{mnI}^{efA} <mn||ef> + \sum_{mNeF} t_{mNI}^{eFA} <mN|eF> + \frac{1}{4} \sum_{MNEF} t_{MNI}^{EFA} <MN||EF>\rightarrow  \{ \bar{H}_{I}^{A} \} \f]
 */
void CCMRCC::build_t1_IA_amplitudes_triples()
{
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);

    // Grab the temporary matrices
    CCMatTmp  HIAMatTmp     = blas->get_MatTmp("t1_eqns[O][V]",unique_ref,none);
    CCMatTmp  TijKabCMatTmp = blas->get_MatTmp("t3[ooO][vvV]",unique_ref,none);
    CCMatTmp  TiJKaBCMatTmp = blas->get_MatTmp("t3[oOO][vVV]",unique_ref,none);
    CCMatTmp  TIJKABCMatTmp = blas->get_MatTmp("t3[OOO][VVV]",unique_ref,none);
    CCMatTmp  ImnefMatTmp   = blas->get_MatTmp("<[oo]:[vv]>",none);
    CCMatTmp  ImNeFMatTmp   = blas->get_MatTmp("<[oo]|[vv]>",none);

    // Grab the indexing for t3[ijk][abc]
    short**   I_tuples = HIAMatTmp->get_left()->get_tuples();
    short**   A_tuples = HIAMatTmp->get_right()->get_tuples();
    short**   mn_tuples = ImnefMatTmp->get_left()->get_tuples();
    short**   ef_tuples = ImnefMatTmp->get_right()->get_tuples();

    double*** TijKabC_matrix = TijKabCMatTmp->get_matrix();
    double*** TiJKaBC_matrix = TiJKaBCMatTmp->get_matrix();
    double*** TIJKABC_matrix = TIJKABCMatTmp->get_matrix();
    double*** HIA_matrix     = HIAMatTmp->get_matrix();
    double*** Imnef_matrix   = ImnefMatTmp->get_matrix();
    double*** ImNeF_matrix   = ImNeFMatTmp->get_matrix();
    CCIndex*  ijkIndex       = blas->get_index("[ooo]");
    CCIndex*  abcIndex       = blas->get_index("[vvv]");

    for(int h =0; h < moinfo->get_nirreps();h++){
      size_t i_offset  = HIAMatTmp->get_left()->get_first(h);
      size_t a_offset = HIAMatTmp->get_right()->get_first(h);
      for(int a = 0;a <HIAMatTmp->get_right_pairpi(h);a++){
        int a_abs = a + a_offset;
        for(int i = 0;i<HIAMatTmp->get_left_pairpi(h);i++){
          int i_abs = i + i_offset;
          // <mn||ef> contribution
          for(int mn_sym =0; mn_sym < moinfo->get_nirreps();mn_sym++){
            size_t mn_offset = ImnefMatTmp->get_left()->get_first(mn_sym);
            size_t ef_offset = ImnefMatTmp->get_right()->get_first(mn_sym);
            for(int ef = 0;ef <ImnefMatTmp->get_right_pairpi(mn_sym);ef++){
              int e = ef_tuples[ef_offset + ef][0];
              int f = ef_tuples[ef_offset + ef][1];
              size_t efa  = abcIndex->get_tuple_rel_index(e,f,a_abs);
              int efa_sym = abcIndex->get_tuple_irrep(e,f,a_abs);
              for(int mn = 0;mn <ImnefMatTmp->get_left_pairpi(mn_sym);mn++){
                int m = mn_tuples[mn_offset + mn][0];
                int n = mn_tuples[mn_offset + mn][1];
                size_t mni  = ijkIndex->get_tuple_rel_index(m,n,i_abs);
                HIA_matrix[h][i][a] += 0.25 * TijKabC_matrix[efa_sym][mni][efa] * Imnef_matrix[mn_sym][mn][ef];
                HIA_matrix[h][i][a] += 0.25 * TIJKABC_matrix[efa_sym][mni][efa] * Imnef_matrix[mn_sym][mn][ef];
                HIA_matrix[h][i][a] +=        TiJKaBC_matrix[efa_sym][mni][efa] * ImNeF_matrix[mn_sym][mn][ef];
              }
            }
          }
        }
      }
    }
  }
//   blas->print("t1_eqns[O][V]{u}");
//   blas->solve("ERROR{u} = 100000000000.0 t1_eqns[O][V]{u} . t1_eqns[O][V]{u}");
//   blas->print("ERROR{u}");
}

}} /* End Namespaces */
