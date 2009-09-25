#include <libmoinfo/libmoinfo.h>

#include "mrccsd_t.h"

namespace psi{ namespace psimrcc{

void MRCCSD_T::compute_ooo_contribution_to_Heff_restricted(int i,int j,int k,int mu,BlockMatrix* T3){
  // Find the off_diagonal elements for reference mu
  // Loop over reference nu (in a safe way)
  for(int nu = 0; nu < nrefs; nu++){
    if(nu != mu){
      vector<pair<int,int> >  alpha_internal_excitation = moinfo->get_alpha_internal_excitation(mu,nu);
      vector<pair<int,int> >   beta_internal_excitation = moinfo->get_beta_internal_excitation(mu,nu);
      double                   sign_internal_excitation = moinfo->get_sign_internal_excitation(mu,nu);

      // Set (alpha)->(alpha) single excitations
      if((alpha_internal_excitation.size() == 1) and (beta_internal_excitation.size() == 0)){
        d_h_eff[nu][mu] += sign_internal_excitation * compute_A_ooo_contribution_to_Heff_restricted(alpha_internal_excitation[0].first,alpha_internal_excitation[0].second,i,j,k,mu,T3);
      }
    }
  }
}

void MRCCSD_T::compute_ooO_contribution_to_Heff_restricted(int i,int j,int k,int mu,BlockMatrix* T3){
  // Find the off_diagonal elements for reference mu
  // Loop over reference nu (in a safe way)
  for(int nu = 0; nu < nrefs; nu++){
    if(nu != mu){
      vector<pair<int,int> >  alpha_internal_excitation = moinfo->get_alpha_internal_excitation(mu,nu);
      vector<pair<int,int> >   beta_internal_excitation = moinfo->get_beta_internal_excitation(mu,nu);
      double                   sign_internal_excitation = moinfo->get_sign_internal_excitation(mu,nu);

      // Set (alpha)->(alpha) single excitations
      if((alpha_internal_excitation.size() == 1) and (beta_internal_excitation.size() == 0)){
        d_h_eff[nu][mu] += sign_internal_excitation * compute_A_ooO_contribution_to_Heff_restricted(alpha_internal_excitation[0].first,alpha_internal_excitation[0].second,i,j,k,mu,T3);
      }
      // Set (beta)->(beta) single excitations
      if((alpha_internal_excitation.size() == 0) and (beta_internal_excitation.size() == 1)){
        d_h_eff[nu][mu] += sign_internal_excitation * compute_B_ooO_contribution_to_Heff_restricted(beta_internal_excitation[0].first,beta_internal_excitation[0].second,i,j,k,mu,T3);
      }
      // Set (alpha,beta)->(alpha,beta) double excitations
      if((alpha_internal_excitation.size() == 1) and (beta_internal_excitation.size() == 1)){
        d_h_eff[nu][mu] += sign_internal_excitation * compute_AB_ooO_contribution_to_Heff_restricted(alpha_internal_excitation[0].first, beta_internal_excitation[0].first,
                                                               alpha_internal_excitation[0].second,beta_internal_excitation[0].second,i,j,k,mu,T3);
      }
    }
  }
}

void MRCCSD_T::compute_oOO_contribution_to_Heff_restricted(int i,int j,int k,int mu,BlockMatrix* T3){
  // Find the off_diagonal elements for reference mu
  // Loop over reference nu (in a safe way)
  for(int nu = 0; nu < nrefs; nu++){
    if(nu != mu){
      vector<pair<int,int> >  alpha_internal_excitation = moinfo->get_alpha_internal_excitation(mu,nu);
      vector<pair<int,int> >   beta_internal_excitation = moinfo->get_beta_internal_excitation(mu,nu);
      double                   sign_internal_excitation = moinfo->get_sign_internal_excitation(mu,nu);

      // Set (alpha)->(alpha) single excitations
      if((alpha_internal_excitation.size() == 1) and (beta_internal_excitation.size() == 0)){
        d_h_eff[nu][mu] += sign_internal_excitation * compute_A_oOO_contribution_to_Heff_restricted(alpha_internal_excitation[0].first,alpha_internal_excitation[0].second,i,j,k,mu,T3);
      }
      // Set (beta)->(beta) single excitations
      if((alpha_internal_excitation.size() == 0) and (beta_internal_excitation.size() == 1)){
        d_h_eff[nu][mu] += sign_internal_excitation * compute_B_oOO_contribution_to_Heff_restricted(beta_internal_excitation[0].first,beta_internal_excitation[0].second,i,j,k,mu,T3);
      }
      // Set (alpha,beta)->(alpha,beta) double excitations
      if((alpha_internal_excitation.size() == 1) and (beta_internal_excitation.size() == 1)){
        d_h_eff[nu][mu] += sign_internal_excitation * compute_AB_oOO_contribution_to_Heff_restricted(alpha_internal_excitation[0].first, beta_internal_excitation[0].first,
                                                               alpha_internal_excitation[0].second,beta_internal_excitation[0].second,i,j,k,mu,T3);
      }
    }
  }
}

void MRCCSD_T::compute_OOO_contribution_to_Heff_restricted(int i,int j,int k,int mu,BlockMatrix* T3){
  // Find the off_diagonal elements for reference mu
  // Loop over reference nu (in a safe way)
  for(int nu = 0; nu < nrefs; nu++){
    if(nu != mu){
      vector<pair<int,int> >  alpha_internal_excitation = moinfo->get_alpha_internal_excitation(mu,nu);
      vector<pair<int,int> >   beta_internal_excitation = moinfo->get_beta_internal_excitation(mu,nu);
      double                   sign_internal_excitation = moinfo->get_sign_internal_excitation(mu,nu);

      // Set (beta)->(beta) single excitations
      if((alpha_internal_excitation.size() == 0) and (beta_internal_excitation.size() == 1)){
        d_h_eff[nu][mu] += sign_internal_excitation * compute_B_OOO_contribution_to_Heff_restricted(beta_internal_excitation[0].first,beta_internal_excitation[0].second,i,j,k,mu,T3);
      }
    }
  }
}

}} /* End Namespaces */
