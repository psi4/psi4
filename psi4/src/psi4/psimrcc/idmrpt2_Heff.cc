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
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include "blas.h"
#include "debugging.h"
#include "idmrpt2.h"
#include "matrix.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;

using namespace std;

void IDMRPT2::build_Heff_mrpt2_offdiagonal()
{
  build_Heff_uv();
  build_Heff_UV();
  build_Heff_uVxY();
  build_Heff_uvxy();
  build_Heff_UVXY();

  intvec occ_to_act = moinfo->get_occ_to_actv();
  intvec vir_to_act = moinfo->get_vir_to_actv();

  for(int i = 0; i < moinfo->get_ref_size(AllRefs); ++i){
    int i_unique = moinfo->get_ref_number(i);
    // Find the off_diagonal elements for reference i
    // Loop over reference j (in a safe way)
    for(int j=0;j<moinfo->get_ref_size(AllRefs);j++){
      if(i!=j){
        vector<pair<int,int> >  alpha_internal_excitation = moinfo->get_alpha_internal_excitation(i,j);
        vector<pair<int,int> >   beta_internal_excitation = moinfo->get_beta_internal_excitation(i,j);
        double                   sign_internal_excitation = moinfo->get_sign_internal_excitation(i,j);

        double element = 0.0;
        if(i==i_unique){
          // Set alpha-alpha single excitations
          if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==0))
            element=sign_internal_excitation * blas->get_MatTmp("Hia[a][a]",i_unique,none)->get_two_address_element(
                                               occ_to_act[alpha_internal_excitation[0].first],
                                               vir_to_act[alpha_internal_excitation[0].second]);

          // Set beta-beta single excitations
          if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==1))
            element=sign_internal_excitation * blas->get_MatTmp("HIA[A][A]",i_unique,none)->get_two_address_element(
                                               occ_to_act[beta_internal_excitation[0].first],
                                               vir_to_act[beta_internal_excitation[0].second]);

          // Set (alpha,alpha)->(alpha,alpha) double excitations
          if((alpha_internal_excitation.size()==2)&&(beta_internal_excitation.size()==0))
            element=sign_internal_excitation * blas->get_MatTmp("Hijab[aa][aa]",i_unique,none)->get_four_address_element(
                                               occ_to_act[alpha_internal_excitation[0].first],
                                               occ_to_act[alpha_internal_excitation[1].first],
                                               vir_to_act[alpha_internal_excitation[0].second],
                                               vir_to_act[alpha_internal_excitation[1].second]);

          // Set (alpha,beta)->(alpha,beta) double excitations
          if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==1))
            element=sign_internal_excitation * blas->get_MatTmp("HiJaB[aA][aA]",i_unique,none)->get_four_address_element(
                                               occ_to_act[alpha_internal_excitation[0].first],
                                               occ_to_act[beta_internal_excitation[0].first],
                                               vir_to_act[alpha_internal_excitation[0].second],
                                               vir_to_act[beta_internal_excitation[0].second]);

          // Set (beta,beta)->(beta,beta) double excitations
          if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==2))
            element=sign_internal_excitation * blas->get_MatTmp("HIJAB[AA][AA]",i_unique,none)->get_four_address_element(
                                               occ_to_act[beta_internal_excitation[0].first],
                                               occ_to_act[beta_internal_excitation[1].first],
                                               vir_to_act[beta_internal_excitation[0].second],
                                               vir_to_act[beta_internal_excitation[1].second]);
        }else{
          // Set alpha-alpha single excitations
          if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==0))
            element=sign_internal_excitation * blas->get_MatTmp("HIA[A][A]",i_unique,none)->get_two_address_element(
                                               occ_to_act[alpha_internal_excitation[0].first],
                                               vir_to_act[alpha_internal_excitation[0].second]);

          // Set beta-beta single excitations
          if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==1))
            element=sign_internal_excitation * blas->get_MatTmp("Hia[a][a]",i_unique,none)->get_two_address_element(
                                               occ_to_act[beta_internal_excitation[0].first],
                                               vir_to_act[beta_internal_excitation[0].second]);

          // Set (alpha,alpha)->(alpha,alpha) double excitations
          if((alpha_internal_excitation.size()==2)&&(beta_internal_excitation.size()==0))
            element=sign_internal_excitation * blas->get_MatTmp("HIJAB[AA][AA]",i_unique,none)->get_four_address_element(
                                               occ_to_act[alpha_internal_excitation[0].first],
                                               occ_to_act[alpha_internal_excitation[1].first],
                                               vir_to_act[alpha_internal_excitation[0].second],
                                               vir_to_act[alpha_internal_excitation[1].second]);

          // Set (alpha,beta)->(alpha,beta) double excitations
          if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==1))
            element=sign_internal_excitation * blas->get_MatTmp("HiJaB[aA][aA]",i_unique,none)->get_four_address_element(
                                               occ_to_act[beta_internal_excitation[0].first],
                                               occ_to_act[alpha_internal_excitation[0].first],
                                               vir_to_act[beta_internal_excitation[0].second],
                                               vir_to_act[alpha_internal_excitation[0].second]);

          // Set (beta,beta)->(beta,beta) double excitations
          if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==2))
            element=sign_internal_excitation * blas->get_MatTmp("Hijab[aa][aa]",i_unique,none)->get_four_address_element(
                                               occ_to_act[beta_internal_excitation[0].first],
                                               occ_to_act[beta_internal_excitation[1].first],
                                               vir_to_act[beta_internal_excitation[0].second],
                                               vir_to_act[beta_internal_excitation[1].second]);
        }
        Heff_mrpt2[j][i]=element;
      }
    }
  }
}




void IDMRPT2::build_Heff_ijkabc()
{
  Timer timer;
  DEBUGGING(1,
    outfile->Printf("\n\tBuilding the H_ijkabc Matrix Elements   ...");

  );

  blas->reduce_spaces("t2_oovv[aaa][v]{u}","t2[oov][v]{u}");
  blas->reduce_spaces("t2_ovvo[aaa][o]{u}","t2[ovv][o]{u}");

  blas->solve("Hijkabc[aaa][aaa]{u}  = #124653#   t2_oovv[aaa][v]{u} 2@2 <[aaa]:[v]>");
  blas->solve("Hijkabc[aaa][aaa]{u}  = #324651# - t2_oovv[aaa][v]{u} 2@2 <[aaa]:[v]>");
  blas->solve("Hijkabc[aaa][aaa]{u}  = #134652# - t2_oovv[aaa][v]{u} 2@2 <[aaa]:[v]>");

  blas->solve("Hijkabc[aaa][aaa]{u}  = #126453# - t2_oovv[aaa][v]{u} 2@2 <[aaa]:[v]>");
  blas->solve("Hijkabc[aaa][aaa]{u}  = #326451#   t2_oovv[aaa][v]{u} 2@2 <[aaa]:[v]>");
  blas->solve("Hijkabc[aaa][aaa]{u}  = #136452#   t2_oovv[aaa][v]{u} 2@2 <[aaa]:[v]>");

  blas->solve("Hijkabc[aaa][aaa]{u}  = #125643# - t2_oovv[aaa][v]{u} 2@2 <[aaa]:[v]>");
  blas->solve("Hijkabc[aaa][aaa]{u}  = #325641#   t2_oovv[aaa][v]{u} 2@2 <[aaa]:[v]>");
  blas->solve("Hijkabc[aaa][aaa]{u}  = #135642#   t2_oovv[aaa][v]{u} 2@2 <[aaa]:[v]>");


  blas->solve("Hijkabc[aaa][aaa]{u}  = #145623#   t2_ovvo[aaa][o]{u} 2@1 <[o]:[aaa]>");
  blas->solve("Hijkabc[aaa][aaa]{u}  = #245613# - t2_ovvo[aaa][o]{u} 2@1 <[o]:[aaa]>");
  blas->solve("Hijkabc[aaa][aaa]{u}  = #345621# - t2_ovvo[aaa][o]{u} 2@1 <[o]:[aaa]>");

  blas->solve("Hijkabc[aaa][aaa]{u}  = #165423# - t2_ovvo[aaa][o]{u} 2@1 <[o]:[aaa]>");
  blas->solve("Hijkabc[aaa][aaa]{u}  = #265413#   t2_ovvo[aaa][o]{u} 2@1 <[o]:[aaa]>");
  blas->solve("Hijkabc[aaa][aaa]{u}  = #365421#   t2_ovvo[aaa][o]{u} 2@1 <[o]:[aaa]>");

  blas->solve("Hijkabc[aaa][aaa]{u}  = #146523# - t2_ovvo[aaa][o]{u} 2@1 <[o]:[aaa]>");
  blas->solve("Hijkabc[aaa][aaa]{u}  = #246513#   t2_ovvo[aaa][o]{u} 2@1 <[o]:[aaa]>");
  blas->solve("Hijkabc[aaa][aaa]{u}  = #346521#   t2_ovvo[aaa][o]{u} 2@1 <[o]:[aaa]>");


  DEBUGGING(3,
    blas->print("Hijkabc[aaa][aaa]{u}");
  );

  DEBUGGING(1,
    outfile->Printf(" done. Timing %20.6f s",timer.get());

  );
}

void IDMRPT2::build_Heff_IJKABC()
{
  Timer timer;
  DEBUGGING(1,
    outfile->Printf("\n\tBuilding the H_IJKABC Matrix Elements   ...");

  );

  blas->reduce_spaces("t2_OOVV[AAA][V]{u}","t2[oov][v]{u}");

  blas->solve("HIJKABC[AAA][AAA]{u}  = #124653#   t2_OOVV[AAA][V]{u} 2@2 <[aaa]:[v]>");
  blas->solve("HIJKABC[AAA][AAA]{u}  = #324651# - t2_OOVV[AAA][V]{u} 2@2 <[aaa]:[v]>");
  blas->solve("HIJKABC[AAA][AAA]{u}  = #134652# - t2_OOVV[AAA][V]{u} 2@2 <[aaa]:[v]>");

  blas->solve("HIJKABC[AAA][AAA]{u}  = #126453# - t2_OOVV[AAA][V]{u} 2@2 <[aaa]:[v]>");
  blas->solve("HIJKABC[AAA][AAA]{u}  = #326451#   t2_OOVV[AAA][V]{u} 2@2 <[aaa]:[v]>");
  blas->solve("HIJKABC[AAA][AAA]{u}  = #136452#   t2_OOVV[AAA][V]{u} 2@2 <[aaa]:[v]>");

  blas->solve("HIJKABC[AAA][AAA]{u}  = #125643# - t2_OOVV[AAA][V]{u} 2@2 <[aaa]:[v]>");
  blas->solve("HIJKABC[AAA][AAA]{u}  = #325641#   t2_OOVV[AAA][V]{u} 2@2 <[aaa]:[v]>");
  blas->solve("HIJKABC[AAA][AAA]{u}  = #135642#   t2_OOVV[AAA][V]{u} 2@2 <[aaa]:[v]>");


  blas->solve("HIJKABC[AAA][AAA]{u}  = #145623#   t2_OVVO[AAA][O]{u} 2@1 <[o]:[aaa]>");
  blas->solve("HIJKABC[AAA][AAA]{u}  = #245613# - t2_OVVO[AAA][O]{u} 2@1 <[o]:[aaa]>");
  blas->solve("HIJKABC[AAA][AAA]{u}  = #345621# - t2_OVVO[AAA][O]{u} 2@1 <[o]:[aaa]>");

  blas->solve("HIJKABC[AAA][AAA]{u}  = #165423# - t2_OVVO[AAA][O]{u} 2@1 <[o]:[aaa]>");
  blas->solve("HIJKABC[AAA][AAA]{u}  = #265413#   t2_OVVO[AAA][O]{u} 2@1 <[o]:[aaa]>");
  blas->solve("HIJKABC[AAA][AAA]{u}  = #365421#   t2_OVVO[AAA][O]{u} 2@1 <[o]:[aaa]>");

  blas->solve("HIJKABC[AAA][AAA]{u}  = #146523# - t2_OVVO[AAA][O]{u} 2@1 <[o]:[aaa]>");
  blas->solve("HIJKABC[AAA][AAA]{u}  = #246513#   t2_OVVO[AAA][O]{u} 2@1 <[o]:[aaa]>");
  blas->solve("HIJKABC[AAA][AAA]{u}  = #346521#   t2_OVVO[AAA][O]{u} 2@1 <[o]:[aaa]>");

  DEBUGGING(3,
    blas->print("HIJKABC[AAA][AAA]{u}");
  );

  DEBUGGING(1,
    outfile->Printf(" done. Timing %20.6f s",timer.get());

  );
}

}} /* End Namespaces */
