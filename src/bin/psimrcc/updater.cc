/**
 *  @file updater.cc
 *  @ingroup (PSIMRCC)
 *  @brief Contains methods for updating the CC equations
*/

#include <vector>
//#include <string>
//
#include <cstdio>
#include <libmoinfo/libmoinfo.h>
#include <liboptions/liboptions.h>

#include "blas.h"
#include "matrix.h"
#include "matrixtmp.h"
#include "updater.h"

extern FILE *outfile;

namespace psi{ namespace psimrcc{


Updater::Updater()
{
}

Updater::~Updater()
{
}

void Updater::zero_internal_amps()
{
  if(options_get_bool("ZERO_INTERNAL_AMPS")){
    // Zero internal amplitudes for unique reference i
    for(int i=0;i<moinfo->get_nunique();i++){
      int unique_i = moinfo->get_ref_number(i,UniqueRefs);
      // Loop over reference j
      for(int j=0;j<moinfo->get_ref_size(AllRefs);j++){
        vector<pair<int,int> >  alpha_internal_excitation = moinfo->get_alpha_internal_excitation(unique_i,j);
        vector<pair<int,int> >   beta_internal_excitation = moinfo->get_beta_internal_excitation(unique_i,j);

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
    fprintf(outfile,"\n  Warning: the internal amplitudes are not zeroed.\n  This is not proper Mk-MRCC. Size-extensivity might be lost\n");
  }
}


void Updater::zero_t1_internal_amps()
{
  if(options_get_bool("ZERO_INTERNAL_AMPS")){
    // Zero internal amplitudes for unique reference i
    for(int i=0;i<moinfo->get_nunique();i++){
      int unique_i = moinfo->get_ref_number(i,UniqueRefs);
      // Loop over reference j
      for(int j=0;j<moinfo->get_ref_size(AllRefs);j++){
        vector<pair<int,int> >  alpha_internal_excitation = moinfo->get_alpha_internal_excitation(unique_i,j);
        vector<pair<int,int> >   beta_internal_excitation = moinfo->get_beta_internal_excitation(unique_i,j);

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
    fprintf(outfile,"\n  Warning: the internal amplitudes are not zeroed.\n  This is not proper Mk-MRCC. Size-extensivity might be lost\n");
  }
}

void Updater::zero_internal_delta_amps()
{
  if(options_get_bool("ZERO_INTERNAL_AMPS")){
    // Zero internal amplitudes for unique reference i
    for(int i=0;i<moinfo->get_nunique();i++){
      int unique_i = moinfo->get_ref_number(i,UniqueRefs);
      // Loop over reference j
      for(int j=0;j<moinfo->get_ref_size(AllRefs);j++){
        vector<pair<int,int> >  alpha_internal_excitation = moinfo->get_alpha_internal_excitation(unique_i,j);
        vector<pair<int,int> >   beta_internal_excitation = moinfo->get_beta_internal_excitation(unique_i,j);

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
    fprintf(outfile,"\n  Warning: the internal amplitudes are not zeroed.\n  This is not proper Mk-MRCC. Size-extensivity might be lost\n");
  }
}

}}
