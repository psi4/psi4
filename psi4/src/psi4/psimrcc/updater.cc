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
 *  @file updater.cc
 *  @ingroup (PSIMRCC)
 *  @brief Contains methods for updating the CC equations
*/

#include <vector>
//#include <string>
//
#include <cstdio>
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/liboptions/liboptions.h"

#include "blas.h"
#include "matrix.h"
#include "matrixtmp.h"
#include "updater.h"

namespace psi{

    namespace psimrcc{
    extern MOInfo *moinfo;


Updater::Updater(Options &options):
        options_(options)
{
}

Updater::~Updater()
{
}

void Updater::zero_internal_amps()
{
  if(options_.get_bool("ZERO_INTERNAL_AMPS")){
    // Zero internal amplitudes for unique reference i
    for(int i=0;i<moinfo->get_nunique();i++){
      int unique_i = moinfo->get_ref_number(i,UniqueRefs);
      // Loop over reference j
      for(int j=0;j<moinfo->get_ref_size(AllRefs);j++){
        std::vector<std::pair<int,int> >  alpha_internal_excitation = moinfo->get_alpha_internal_excitation(unique_i,j);
        std::vector<std::pair<int,int> >   beta_internal_excitation = moinfo->get_beta_internal_excitation(unique_i,j);

        // Zero alpha-alpha single excitations
        if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==0)){
          blas->get_MatTmp("t1[o][v]",unique_i,none)->set_two_address_element(
                                            alpha_internal_excitation[0].first,
                                            alpha_internal_excitation[0].second,
                                            0.0);
        }

        // Zero beta-beta single excitations
        if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==1))
          blas->get_MatTmp("t1[O][V]",unique_i,none)->set_two_address_element(
                                            beta_internal_excitation[0].first,
                                            beta_internal_excitation[0].second,
                                            0.0);

        // Zero (alpha,alpha)->(alpha,alpha) double excitations (all permutations)
        if((alpha_internal_excitation.size()==2)&&(beta_internal_excitation.size()==0)){
          blas->get_MatTmp("t2[oo][vv]",unique_i,none)->set_four_address_element(
                                            alpha_internal_excitation[0].first,
                                            alpha_internal_excitation[1].first,
                                            alpha_internal_excitation[0].second,
                                            alpha_internal_excitation[1].second,
                                            0.0);
          blas->get_MatTmp("t2[oo][vv]",unique_i,none)->set_four_address_element(
                                            alpha_internal_excitation[0].first,
                                            alpha_internal_excitation[1].first,
                                            alpha_internal_excitation[1].second,
                                            alpha_internal_excitation[0].second,
                                            0.0);
          blas->get_MatTmp("t2[oo][vv]",unique_i,none)->set_four_address_element(
                                            alpha_internal_excitation[1].first,
                                            alpha_internal_excitation[0].first,
                                            alpha_internal_excitation[0].second,
                                            alpha_internal_excitation[1].second,
                                            0.0);
          blas->get_MatTmp("t2[oo][vv]",unique_i,none)->set_four_address_element(
                                            alpha_internal_excitation[1].first,
                                            alpha_internal_excitation[0].first,
                                            alpha_internal_excitation[1].second,
                                            alpha_internal_excitation[0].second,
                                            0.0);
        }

        // Zero (alpha,beta)->(alpha,beta) double excitations
        if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==1)){
          blas->get_MatTmp("t2[oO][vV]",unique_i,none)->set_four_address_element(
                                            alpha_internal_excitation[0].first,
                                            beta_internal_excitation[0].first,
                                            alpha_internal_excitation[0].second,
                                            beta_internal_excitation[0].second,
                                            0.0);
        }

        // Zero (beta,beta)->(beta,beta) double excitations (all permutations)
        if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==2)){
          blas->get_MatTmp("t2[OO][VV]",unique_i,none)->set_four_address_element(
                                            beta_internal_excitation[0].first,
                                            beta_internal_excitation[1].first,
                                            beta_internal_excitation[0].second,
                                            beta_internal_excitation[1].second,
                                            0.0);
          blas->get_MatTmp("t2[OO][VV]",unique_i,none)->set_four_address_element(
                                            beta_internal_excitation[0].first,
                                            beta_internal_excitation[1].first,
                                            beta_internal_excitation[1].second,
                                            beta_internal_excitation[0].second,
                                            0.0);
          blas->get_MatTmp("t2[OO][VV]",unique_i,none)->set_four_address_element(
                                            beta_internal_excitation[1].first,
                                            beta_internal_excitation[0].first,
                                            beta_internal_excitation[0].second,
                                            beta_internal_excitation[1].second,
                                            0.0);
          blas->get_MatTmp("t2[OO][VV]",unique_i,none)->set_four_address_element(
                                            beta_internal_excitation[1].first,
                                            beta_internal_excitation[0].first,
                                            beta_internal_excitation[1].second,
                                            beta_internal_excitation[0].second,
                                            0.0);
        }
      }
    }

    // Print the t-amplitudes
//    DEBUGGING(3,
//      blas->print("t1[o][v]{u}");
//      blas->print("t1[O][V]{u}");
//      blas->print("t2[oo][vv]{u}");
//      blas->print("t2[oO][vV]{u}");
//      blas->print("t2[OO][VV]{u}");
//    )
  }else{
    outfile->Printf("\n  Warning: the internal amplitudes are not zeroed.\n  This is not proper Mk-MRCC. Size-extensivity might be lost\n");
  }
}


void Updater::zero_t1_internal_amps()
{
  if(options_.get_bool("ZERO_INTERNAL_AMPS")){
    // Zero internal amplitudes for unique reference i
    for(int i=0;i<moinfo->get_nunique();i++){
      int unique_i = moinfo->get_ref_number(i,UniqueRefs);
      // Loop over reference j
      for(int j=0;j<moinfo->get_ref_size(AllRefs);j++){
        std::vector<std::pair<int,int> >  alpha_internal_excitation = moinfo->get_alpha_internal_excitation(unique_i,j);
        std::vector<std::pair<int,int> >   beta_internal_excitation = moinfo->get_beta_internal_excitation(unique_i,j);

        // Zero alpha-alpha single excitations
        if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==0))
          blas->get_MatTmp("t1[o][v]",unique_i,none)->set_two_address_element(
                                            alpha_internal_excitation[0].first,
                                            alpha_internal_excitation[0].second,
                                            0.0);

        // Zero beta-beta single excitations
        if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==1))
          blas->get_MatTmp("t1[O][V]",unique_i,none)->set_two_address_element(
                                            beta_internal_excitation[0].first,
                                            beta_internal_excitation[0].second,
                                            0.0);
      }
    }

    // Print the t-amplitudes
//    DEBUGGING(3,
//      blas->print("t1[o][v]{u}");
//      blas->print("t1[O][V]{u}");
//    )
  }else{
    outfile->Printf("\n  Warning: the internal amplitudes are not zeroed.\n  This is not proper Mk-MRCC. Size-extensivity might be lost\n");
  }
}

void Updater::zero_internal_delta_amps()
{
  if(options_.get_bool("ZERO_INTERNAL_AMPS")){
    // Zero internal amplitudes for unique reference i
    for(int i=0;i<moinfo->get_nunique();i++){
      int unique_i = moinfo->get_ref_number(i,UniqueRefs);
      // Loop over reference j
      for(int j=0;j<moinfo->get_ref_size(AllRefs);j++){
        std::vector<std::pair<int,int> >  alpha_internal_excitation = moinfo->get_alpha_internal_excitation(unique_i,j);
        std::vector<std::pair<int,int> >   beta_internal_excitation = moinfo->get_beta_internal_excitation(unique_i,j);

        // Zero alpha-alpha single excitations
        if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==0))
          blas->get_MatTmp("t1_delta[o][v]",unique_i,none)->set_two_address_element(
                                            alpha_internal_excitation[0].first,
                                            alpha_internal_excitation[0].second,
                                            0.0);

        // Zero beta-beta single excitations
        if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==1))
          blas->get_MatTmp("t1_delta[O][V]",unique_i,none)->set_two_address_element(
                                            beta_internal_excitation[0].first,
                                            beta_internal_excitation[0].second,
                                            0.0);

        // Zero (alpha,alpha)->(alpha,alpha) double excitations (all permutations)
        if((alpha_internal_excitation.size()==2)&&(beta_internal_excitation.size()==0)){
          blas->get_MatTmp("t2_delta[oo][vv]",unique_i,none)->set_four_address_element(
                                            alpha_internal_excitation[0].first,
                                            alpha_internal_excitation[1].first,
                                            alpha_internal_excitation[0].second,
                                            alpha_internal_excitation[1].second,
                                            0.0);
          blas->get_MatTmp("t2_delta[oo][vv]",unique_i,none)->set_four_address_element(
                                            alpha_internal_excitation[0].first,
                                            alpha_internal_excitation[1].first,
                                            alpha_internal_excitation[1].second,
                                            alpha_internal_excitation[0].second,
                                            0.0);
          blas->get_MatTmp("t2_delta[oo][vv]",unique_i,none)->set_four_address_element(
                                            alpha_internal_excitation[1].first,
                                            alpha_internal_excitation[0].first,
                                            alpha_internal_excitation[0].second,
                                            alpha_internal_excitation[1].second,
                                            0.0);
          blas->get_MatTmp("t2_delta[oo][vv]",unique_i,none)->set_four_address_element(
                                            alpha_internal_excitation[1].first,
                                            alpha_internal_excitation[0].first,
                                            alpha_internal_excitation[1].second,
                                            alpha_internal_excitation[0].second,
                                            0.0);
        }

        // Zero (alpha,beta)->(alpha,beta) double excitations
        if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==1)){
          blas->get_MatTmp("t2_delta[oO][vV]",unique_i,none)->set_four_address_element(
                                            alpha_internal_excitation[0].first,
                                            beta_internal_excitation[0].first,
                                            alpha_internal_excitation[0].second,
                                            beta_internal_excitation[0].second,
                                            0.0);
        }

        // Zero (beta,beta)->(beta,beta) double excitations (all permutations)
        if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==2)){
          blas->get_MatTmp("t2_delta[OO][VV]",unique_i,none)->set_four_address_element(
                                            beta_internal_excitation[0].first,
                                            beta_internal_excitation[1].first,
                                            beta_internal_excitation[0].second,
                                            beta_internal_excitation[1].second,
                                            0.0);
          blas->get_MatTmp("t2_delta[OO][VV]",unique_i,none)->set_four_address_element(
                                            beta_internal_excitation[0].first,
                                            beta_internal_excitation[1].first,
                                            beta_internal_excitation[1].second,
                                            beta_internal_excitation[0].second,
                                            0.0);
          blas->get_MatTmp("t2_delta[OO][VV]",unique_i,none)->set_four_address_element(
                                            beta_internal_excitation[1].first,
                                            beta_internal_excitation[0].first,
                                            beta_internal_excitation[0].second,
                                            beta_internal_excitation[1].second,
                                            0.0);
          blas->get_MatTmp("t2_delta[OO][VV]",unique_i,none)->set_four_address_element(
                                            beta_internal_excitation[1].first,
                                            beta_internal_excitation[0].first,
                                            beta_internal_excitation[1].second,
                                            beta_internal_excitation[0].second,
                                            0.0);
        }
      }
    }

    // Print the t-amplitudes
//    DEBUGGING(3,
//      blas->print("t1_delta[o][v]{u}");
//      blas->print("t1_delta[O][V]{u}");
//      blas->print("t2_delta[oo][vv]{u}");
//      blas->print("t2_delta[oO][vV]{u}");
//      blas->print("t2_delta[OO][VV]{u}");
//    )
  }else{
    outfile->Printf("\n  Warning: the internal amplitudes are not zeroed.\n  This is not proper Mk-MRCC. Size-extensivity might be lost\n");
  }
}

}}
