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
#include "debugging.h"
#include "index.h"
#include "mrcc.h"
#include "matrix.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;

using namespace std;

void CCMRCC::synchronize_amps()
{
  blas->solve("t1[ov]{u}     = #12# t1[o][v]{u}");
  blas->solve("t1[OV]{u}     = #12# t1[O][V]{u}");
  blas->solve("t2[ov][OV]{u} = #1324# t2[oO][vV]{u}");
  blas->solve("t2[ov][ov]{u} = #1324# t2[oo][vv]{u}");
  blas->solve("t2[OV][OV]{u} = #1324# t2[OO][VV]{u}");
  blas->solve("t2[oV][Ov]{u} = #1342# t2[oO][vV]{u}");

  blas->solve("t2[o][ovv]{u} = #1234# t2[oo][vv]{u}");
  blas->solve("t2[o][OvV]{u} = #1234# t2[oO][vV]{u}");
  blas->solve("t2[O][oVv]{u} = #2143# t2[oO][vV]{u}");
  blas->solve("t2[O][OVV]{u} = #1234# t2[OO][VV]{u}");

  blas->solve("t2[v][voo]{u} = #3412# t2[oo][vv]{u}");
  blas->solve("t2[v][VoO]{u} = #3412# t2[oO][vV]{u}");
  blas->solve("t2[V][vOo]{u} = #4321# t2[oO][vV]{u}");
  blas->solve("t2[V][VOO]{u} = #3412# t2[OO][VV]{u}");

  if(triples_type == ccsd_t){
    blas->solve("t2[Oo][Vv]{u} = #2143# t2[oO][vV]{u}");
  }

  if(triples_type > ccsd_t){
    blas->solve("t2[ovv][o]{u} = #1423# t2[oo][vv]{u}"); // T2_iab_j
    blas->solve("t2[oVv][O]{u} = #1432# t2[oO][vV]{u}"); // T2_iBa_J
    blas->solve("t2[OvV][o]{u} = #4123# t2[oO][vV]{u}"); // T2_JaB_i
    blas->solve("t2[OVV][O]{u} = #1423# t2[OO][VV]{u}"); // T2_IAB_J

    blas->solve("t2[oov][v]{u} = #1234# t2[oo][vv]{u}"); // T2_ija_b
    blas->solve("t2[oOv][V]{u} = #1234# t2[oO][vV]{u}"); // T2_iJa_B
    blas->solve("t2[OoV][v]{u} = #2143# t2[oO][vV]{u}"); // T2_JiB_a
    blas->solve("t2[OOV][V]{u} = #1234# t2[OO][VV]{u}"); // T2_IJA_B
  }
}

void CCMRCC::compute_delta_amps()
{
  blas->solve("||Delta_t1||{u}  = t1_delta[o][v]{u} . t1_delta[o][v]{u}");
  blas->solve("||Delta_t1||{u} += t1_delta[O][V]{u} . t1_delta[O][V]{u}");

  blas->solve("||Delta_t2||{u}  = t2_delta[oo][vv]{u} . t2_delta[oo][vv]{u}");
  blas->solve("||Delta_t2||{u} += t2_delta[oO][vV]{u} . t2_delta[oO][vV]{u}");
  blas->solve("||Delta_t2||{u} += t2_delta[OO][VV]{u} . t2_delta[OO][VV]{u}");

  // Compute the T-AMPS difference
  delta_t1_amps=0.0;
  delta_t2_amps=0.0;
  for(int n = 0; n < moinfo->get_ref_size(AllRefs); n++){
    double c_n2 = std::pow(h_eff.get_right_eigenvector(n),2.0);
    delta_t1_amps += c_n2 * blas->get_scalar("||Delta_t1||",moinfo->get_ref_number(n));
    delta_t2_amps += c_n2 * blas->get_scalar("||Delta_t2||",moinfo->get_ref_number(n));
  }
  delta_t1_amps = std::sqrt(delta_t1_amps);
  delta_t2_amps = std::sqrt(delta_t2_amps);
}

void CCMRCC::update_t3_amps()
{
  update_t3_ijkabc_amps();
  update_t3_ijKabC_amps();
  update_t3_iJKaBC_amps();
  update_t3_IJKABC_amps();
}

void CCMRCC::update_t3_ijkabc_amps()
{
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);

    // Grab the temporary matrices
    CCMatTmp  TijkabcMatTmp = blas->get_MatTmp("t3[ooo][vvv]",unique_ref,none);
    CCMatTmp  HijkabcMatTmp = blas->get_MatTmp("t3_eqns[ooo][vvv]",unique_ref,none);

    double*** Tijkabc_matrix = TijkabcMatTmp->get_matrix();
    double*** Hijkabc_matrix = HijkabcMatTmp->get_matrix();

    for(int h =0; h < moinfo->get_nirreps();h++){
      for(size_t abc = 0;abc<TijkabcMatTmp->get_right_pairpi(h);abc++){
        double delta_abc = d3_vvv[ref][h][abc];
        for(size_t ijk = 0;ijk<TijkabcMatTmp->get_left_pairpi(h);ijk++){
          double delta_ijk = d3_ooo[ref][h][ijk];
          Tijkabc_matrix[h][ijk][abc]+=Hijkabc_matrix[h][ijk][abc]/(delta_ijk-delta_abc);
        }
      }
    }
  }
//   blas->print("t3[ooo][vvv]{u}");
//   blas->solve("ERROR{u} = 1000000.0 t3[ooo][vvv]{u} . t3[ooo][vvv]{u}");
//   blas->print("ERROR{u}");
}

void CCMRCC::update_t3_ijKabC_amps()
{
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);

    // Grab the temporary matrices
    CCMatTmp  TijKabCMatTmp = blas->get_MatTmp("t3[ooO][vvV]",unique_ref,none);
    CCMatTmp  HijKabCMatTmp = blas->get_MatTmp("t3_eqns[ooO][vvV]",unique_ref,none);

    double*** TijKabC_matrix = TijKabCMatTmp->get_matrix();
    double*** HijKabC_matrix = HijKabCMatTmp->get_matrix();

    for(int h =0; h < moinfo->get_nirreps();h++){
      for(size_t abc = 0;abc<TijKabCMatTmp->get_right_pairpi(h);abc++){
        double delta_abc = d3_vvV[ref][h][abc];
        for(size_t ijk = 0;ijk<TijKabCMatTmp->get_left_pairpi(h);ijk++){
          double delta_ijk = d3_ooO[ref][h][ijk];
          TijKabC_matrix[h][ijk][abc]+=HijKabC_matrix[h][ijk][abc]/(delta_ijk-delta_abc);
        }
      }
    }
  }
//   blas->print("t3[ooO][vvV]{u}");
//   blas->solve("ERROR{u} = 1000000.0 t3[ooO][vvV]{u} . t3[ooO][vvV]{u}");
//   blas->print("ERROR{u}");
}


void CCMRCC::update_t3_iJKaBC_amps()
{
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);

    // Grab the temporary matrices
    CCMatTmp  TiJKaBCMatTmp = blas->get_MatTmp("t3[oOO][vVV]",unique_ref,none);
    CCMatTmp  HiJKaBCMatTmp = blas->get_MatTmp("t3_eqns[oOO][vVV]",unique_ref,none);

    // Grab the indexing for t3[ijk][abc]
    short**   left_tuples  = TiJKaBCMatTmp->get_left()->get_tuples();
    short**   right_tuples = TiJKaBCMatTmp->get_right()->get_tuples();

    double*** TiJKaBC_matrix = TiJKaBCMatTmp->get_matrix();
    double*** HiJKaBC_matrix = HiJKaBCMatTmp->get_matrix();

    for(int h =0; h < moinfo->get_nirreps();h++){
      for(int abc = 0;abc<TiJKaBCMatTmp->get_right_pairpi(h);abc++){
        double delta_abc = d3_vVV[ref][h][abc];
        for(int ijk = 0;ijk<TiJKaBCMatTmp->get_left_pairpi(h);ijk++){
          double delta_ijk = d3_oOO[ref][h][ijk];
          TiJKaBC_matrix[h][ijk][abc]+=HiJKaBC_matrix[h][ijk][abc]/(delta_ijk-delta_abc);
        }
      }
    }
  }
//   blas->print("t3[oOO][vVV]{u}");
//   blas->solve("ERROR{u} = 1000000.0 t3[oOO][vVV]{u} . t3[oOO][vVV]{u}");
//   blas->print("ERROR{u}");
}

void CCMRCC::update_t3_IJKABC_amps()
{
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int unique_ref  = moinfo->get_ref_number(ref,UniqueRefs);

    // Grab the temporary matrices
    CCMatTmp  TIJKABCMatTmp = blas->get_MatTmp("t3[OOO][VVV]",unique_ref,none);
    CCMatTmp  HIJKABCMatTmp = blas->get_MatTmp("t3_eqns[OOO][VVV]",unique_ref,none);

    // Grab the indexing for t3[ijk][abc]
    short**   left_tuples  = TIJKABCMatTmp->get_left()->get_tuples();
    short**   right_tuples = TIJKABCMatTmp->get_right()->get_tuples();

    double*** TIJKABC_matrix = TIJKABCMatTmp->get_matrix();
    double*** HIJKABC_matrix = HIJKABCMatTmp->get_matrix();

    for(int h =0; h < moinfo->get_nirreps();h++){
      for(int abc = 0;abc<TIJKABCMatTmp->get_right_pairpi(h);abc++){
        double delta_abc = d3_VVV[ref][h][abc];
        for(int ijk = 0;ijk<TIJKABCMatTmp->get_left_pairpi(h);ijk++){
          double delta_ijk = d3_OOO[ref][h][ijk];
          TIJKABC_matrix[h][ijk][abc]+=HIJKABC_matrix[h][ijk][abc]/(delta_ijk-delta_abc);
        }
      }
    }
  }
//   blas->print("t3[OOO][VVV]{u}");
//   blas->solve("ERROR{u} = 1000000.0 t3[OOO][VVV]{u} . t3[OOO][VVV]{u}");
//   blas->print("ERROR{u}");
}


// void CCMRCC::update_t3_ijkabc_amps()
// {
//   // Loop over references
//   for(int ref=0;ref<moinfo->get_nunique();ref++){
//     int unique_ref  = moinfo->get_ref_number("u",ref);
//
//     // Grab the temporary matrices
//     CCMatTmp  TijkabcMatTmp = blas->get_MatTmp("t3[ooo][vvv]",unique_ref,none);
//     CCMatTmp  HijkabcMatTmp = blas->get_MatTmp("t3_eqns[ooo][vvv]",unique_ref,none);
//     CCMatTmp  FockmnMatTmp  = blas->get_MatTmp("fock[oo]",unique_ref,none);
//     CCMatTmp  FockefMatTmp  = blas->get_MatTmp("fock[vv]",unique_ref,none);
//
//     // Grab the indexing for t3[ijk][abc]
//     short**   left_tuples  = TijkabcMatTmp->get_left()->get_tuples();
//     short**   right_tuples = TijkabcMatTmp->get_right()->get_tuples();
//
//     double*** Tijkabc_matrix = TijkabcMatTmp->get_matrix();
//     double*** Hijkabc_matrix = HijkabcMatTmp->get_matrix();
//     double*** Fockmn_matrix  = FockmnMatTmp->get_matrix();
//     double*** Fockef_matrix  = FockefMatTmp->get_matrix();
//
//     for(int h =0; h < moinfo->get_nirreps();h++){
//       size_t left_offset  = TijkabcMatTmp->get_left()->get_first(h);
//       size_t right_offset = TijkabcMatTmp->get_right()->get_first(h);
//       for(int abc = 0;abc<TijkabcMatTmp->get_right_pairpi(h);abc++){
//         int a = right_tuples[right_offset + abc][0];
//         int b = right_tuples[right_offset + abc][1];
//         int c = right_tuples[right_offset + abc][2];
//         double delta_abc = FockefMatTmp->get_two_address_element(a,a) +
//                            FockefMatTmp->get_two_address_element(b,b) +
//                            FockefMatTmp->get_two_address_element(c,c);
//         double delta_abc_new = d3_vvv[ref][h][abc];
//         if(delta_abc_new!=delta_abc)
//           outfile->Printf("\nDenominators disagree %d",abc);
//         for(int ijk = 0;ijk<TijkabcMatTmp->get_left_pairpi(h);ijk++){
//           int i = left_tuples[left_offset + ijk][0];
//           int j = left_tuples[left_offset + ijk][1];
//           int k = left_tuples[left_offset + ijk][2];
//           double delta_ijk = FockmnMatTmp->get_two_address_element(i,i) +
//                              FockmnMatTmp->get_two_address_element(j,j) +
//                              FockmnMatTmp->get_two_address_element(k,k);
//           Tijkabc_matrix[h][ijk][abc]+=Hijkabc_matrix[h][ijk][abc]/(delta_ijk-delta_abc);
//         }
//       }
//     }
//   }
// //   blas->print("t3[ooo][vvv]{u}");
// //   blas->solve("ERROR{u} = 1000000.0 t3[ooo][vvv]{u} . t3[ooo][vvv]{u}");
// //   blas->print("ERROR{u}");
// }

// void CCMRCC::update_t3_ijKabC_amps()
// {
//   // Loop over references
//   for(int ref=0;ref<moinfo->get_nunique();ref++){
//     int unique_ref  = moinfo->get_ref_number("u",ref);
//
//     // Grab the temporary matrices
//     CCMatTmp  TijKabCMatTmp = blas->get_MatTmp("t3[ooO][vvV]",unique_ref,none);
//     CCMatTmp  HijKabCMatTmp = blas->get_MatTmp("t3_eqns[ooO][vvV]",unique_ref,none);
//     CCMatTmp  FockmnMatTmp  = blas->get_MatTmp("fock[oo]",unique_ref,none);
//     CCMatTmp  FockMNMatTmp  = blas->get_MatTmp("fock[OO]",unique_ref,none);
//     CCMatTmp  FockefMatTmp  = blas->get_MatTmp("fock[vv]",unique_ref,none);
//     CCMatTmp  FockEFMatTmp  = blas->get_MatTmp("fock[VV]",unique_ref,none);
//
//     // Grab the indexing for t3[ijk][abc]
//     short**   left_tuples  = TijKabCMatTmp->get_left()->get_tuples();
//     short**   right_tuples = TijKabCMatTmp->get_right()->get_tuples();
//
//     double*** TijKabC_matrix = TijKabCMatTmp->get_matrix();
//     double*** HijKabC_matrix = HijKabCMatTmp->get_matrix();
//     double*** Fockmn_matrix  = FockmnMatTmp->get_matrix();
//     double*** FockMN_matrix  = FockMNMatTmp->get_matrix();
//     double*** Fockef_matrix  = FockefMatTmp->get_matrix();
//     double*** FockEF_matrix  = FockEFMatTmp->get_matrix();
//
//     for(int h =0; h < moinfo->get_nirreps();h++){
//       size_t left_offset  = TijKabCMatTmp->get_left()->get_first(h);
//       size_t right_offset = TijKabCMatTmp->get_right()->get_first(h);
//       for(int abc = 0;abc<TijKabCMatTmp->get_right_pairpi(h);abc++){
//         int a = right_tuples[right_offset + abc][0];
//         int b = right_tuples[right_offset + abc][1];
//         int c = right_tuples[right_offset + abc][2];
//         double delta_abc = FockefMatTmp->get_two_address_element(a,a) +
//                            FockefMatTmp->get_two_address_element(b,b) +
//                            FockEFMatTmp->get_two_address_element(c,c);
//         for(int ijk = 0;ijk<TijKabCMatTmp->get_left_pairpi(h);ijk++){
//           int i = left_tuples[left_offset + ijk][0];
//           int j = left_tuples[left_offset + ijk][1];
//           int k = left_tuples[left_offset + ijk][2];
//           double delta_ijk = FockmnMatTmp->get_two_address_element(i,i) +
//                              FockmnMatTmp->get_two_address_element(j,j) +
//                              FockMNMatTmp->get_two_address_element(k,k);
//           TijKabC_matrix[h][ijk][abc]+=HijKabC_matrix[h][ijk][abc]/(delta_ijk-delta_abc);
//         }
//       }
//     }
//   }
// //   blas->print("t3[ooO][vvV]{u}");
// //   blas->solve("ERROR{u} = 1000000.0 t3[ooO][vvV]{u} . t3[ooO][vvV]{u}");
// //   blas->print("ERROR{u}");
// }
//
//
// void CCMRCC::update_t3_iJKaBC_amps()
// {
//   // Loop over references
//   for(int ref=0;ref<moinfo->get_nunique();ref++){
//     int unique_ref  = moinfo->get_ref_number("u",ref);
//
//     // Grab the temporary matrices
//     CCMatTmp  TiJKaBCMatTmp = blas->get_MatTmp("t3[oOO][vVV]",unique_ref,none);
//     CCMatTmp  HiJKaBCMatTmp = blas->get_MatTmp("t3_eqns[oOO][vVV]",unique_ref,none);
//     CCMatTmp  FockmnMatTmp  = blas->get_MatTmp("fock[oo]",unique_ref,none);
//     CCMatTmp  FockMNMatTmp  = blas->get_MatTmp("fock[OO]",unique_ref,none);
//     CCMatTmp  FockefMatTmp  = blas->get_MatTmp("fock[vv]",unique_ref,none);
//     CCMatTmp  FockEFMatTmp  = blas->get_MatTmp("fock[VV]",unique_ref,none);
//
//     // Grab the indexing for t3[ijk][abc]
//     short**   left_tuples  = TiJKaBCMatTmp->get_left()->get_tuples();
//     short**   right_tuples = TiJKaBCMatTmp->get_right()->get_tuples();
//
//     double*** TiJKaBC_matrix = TiJKaBCMatTmp->get_matrix();
//     double*** HiJKaBC_matrix = HiJKaBCMatTmp->get_matrix();
//     double*** Fockmn_matrix  = FockmnMatTmp->get_matrix();
//     double*** FockMN_matrix  = FockMNMatTmp->get_matrix();
//     double*** Fockef_matrix  = FockefMatTmp->get_matrix();
//     double*** FockEF_matrix  = FockEFMatTmp->get_matrix();
//
//     for(int h =0; h < moinfo->get_nirreps();h++){
//       size_t left_offset  = TiJKaBCMatTmp->get_left()->get_first(h);
//       size_t right_offset = TiJKaBCMatTmp->get_right()->get_first(h);
//       for(int abc = 0;abc<TiJKaBCMatTmp->get_right_pairpi(h);abc++){
//         int a = right_tuples[right_offset + abc][0];
//         int b = right_tuples[right_offset + abc][1];
//         int c = right_tuples[right_offset + abc][2];
//         double delta_abc = FockefMatTmp->get_two_address_element(a,a) +
//                            FockEFMatTmp->get_two_address_element(b,b) +
//                            FockEFMatTmp->get_two_address_element(c,c);
//         for(int ijk = 0;ijk<TiJKaBCMatTmp->get_left_pairpi(h);ijk++){
//           int i = left_tuples[left_offset + ijk][0];
//           int j = left_tuples[left_offset + ijk][1];
//           int k = left_tuples[left_offset + ijk][2];
//           double delta_ijk = FockmnMatTmp->get_two_address_element(i,i) +
//                              FockMNMatTmp->get_two_address_element(j,j) +
//                              FockMNMatTmp->get_two_address_element(k,k);
//           TiJKaBC_matrix[h][ijk][abc]+=HiJKaBC_matrix[h][ijk][abc]/(delta_ijk-delta_abc); // TODO, check this when doing MR-CCSD
//         }
//       }
//     }
//   }
// //   blas->print("t3[oOO][vVV]{u}");
// //   blas->solve("ERROR{u} = 1000000.0 t3[oOO][vVV]{u} . t3[oOO][vVV]{u}");
// //   blas->print("ERROR{u}");
// }
//
// void CCMRCC::update_t3_IJKABC_amps()
// {
//   // Loop over references
//   for(int ref=0;ref<moinfo->get_nunique();ref++){
//     int unique_ref  = moinfo->get_ref_number("u",ref);
//
//     // Grab the temporary matrices
//     CCMatTmp  TIJKABCMatTmp = blas->get_MatTmp("t3[OOO][VVV]",unique_ref,none);
//     CCMatTmp  HIJKABCMatTmp = blas->get_MatTmp("t3_eqns[OOO][VVV]",unique_ref,none);
//     CCMatTmp  FockMNMatTmp  = blas->get_MatTmp("fock[OO]",unique_ref,none);
//     CCMatTmp  FockEFMatTmp  = blas->get_MatTmp("fock[VV]",unique_ref,none);
//
//     // Grab the indexing for t3[ijk][abc]
//     short**   left_tuples  = TIJKABCMatTmp->get_left()->get_tuples();
//     short**   right_tuples = TIJKABCMatTmp->get_right()->get_tuples();
//
//     double*** TIJKABC_matrix = TIJKABCMatTmp->get_matrix();
//     double*** HIJKABC_matrix = HIJKABCMatTmp->get_matrix();
//     double*** FockMN_matrix  = FockMNMatTmp->get_matrix();
//     double*** FockEF_matrix  = FockEFMatTmp->get_matrix();
//
//     for(int h =0; h < moinfo->get_nirreps();h++){
//       size_t left_offset  = TIJKABCMatTmp->get_left()->get_first(h);
//       size_t right_offset = TIJKABCMatTmp->get_right()->get_first(h);
//       for(int abc = 0;abc<TIJKABCMatTmp->get_right_pairpi(h);abc++){
//         int a = right_tuples[right_offset + abc][0];
//         int b = right_tuples[right_offset + abc][1];
//         int c = right_tuples[right_offset + abc][2];
//         double delta_abc = FockEFMatTmp->get_two_address_element(a,a) +
//                            FockEFMatTmp->get_two_address_element(b,b) +
//                            FockEFMatTmp->get_two_address_element(c,c);
//         for(int ijk = 0;ijk<TIJKABCMatTmp->get_left_pairpi(h);ijk++){
//           int i = left_tuples[left_offset + ijk][0];
//           int j = left_tuples[left_offset + ijk][1];
//           int k = left_tuples[left_offset + ijk][2];
//           double delta_ijk = FockMNMatTmp->get_two_address_element(i,i) +
//                              FockMNMatTmp->get_two_address_element(j,j) +
//                              FockMNMatTmp->get_two_address_element(k,k);
//           TIJKABC_matrix[h][ijk][abc]+=HIJKABC_matrix[h][ijk][abc]/(delta_ijk-delta_abc);
//         }
//       }
//     }
//   }
// //   blas->print("t3[OOO][VVV]{u}");
// //   blas->solve("ERROR{u} = 1000000.0 t3[OOO][VVV]{u} . t3[OOO][VVV]{u}");
// //   blas->print("ERROR{u}");
// }



}} /* End Namespaces */
