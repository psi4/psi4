/**
 *  @file updater_mk.cc
 *  @ingroup (PSIMRCC)
 *  @brief Contains methods for updating the CC equations
*/

#include <cstdio>
#include <cmath>
#include <string>

#include <libmoinfo/libmoinfo.h>
#include <liboptions/liboptions.h>

#include "blas.h"
#include "heff.h"
#include "updater.h"

extern FILE *outfile;

namespace psi{ namespace psimrcc{

MkUpdater::MkUpdater() : Updater()
{
}

MkUpdater::~MkUpdater()
{
}

void MkUpdater::update(int cycle,Hamiltonian* heff)
{
  // Setup the Tikhonow omega parameter
  double omega = 0;
  int    tikhonow_max   = options_get_int("TIKHONOW_MAX");
  double tikhonow_omega = static_cast<double>(options_get_int("TIKHONOW_OMEGA")) / 1000.0;

  if(tikhonow_max == 0){  // Tikhonow always turned on
    omega = tikhonow_omega;
  }else{
    if(cycle < tikhonow_max){
      omega = tikhonow_omega;
      fprintf(outfile,"\n  Tikhonow regularization turned on.  Omega = %6.3e",omega);
    }
  }

  blas->solve("d'1[o][v]{u}    = d1[o][v]{u}");
  blas->solve("d'1[O][V]{u}    = d1[O][V]{u}");
  blas->solve("d'2[oo][vv]{u}  = d2[oo][vv]{u}");
  blas->solve("d'2[oO][vV]{u}  = d2[oO][vV]{u}");
  blas->solve("d'2[OO][VV]{u}  = d2[OO][VV]{u}");

  // Shift the denominators
  if(options_get_bool("COUPLING_TERMS")){
    for(int mu = 0; mu < moinfo->get_nunique(); ++mu){
      int mu_unique = moinfo->get_ref_number(mu,UniqueRefs);
      double denominator_shift = heff->get_eigenvalue() - heff->get_matrix(mu_unique,mu_unique);
      std::string shift  = to_string(denominator_shift);
      std::string mu_str = to_string(mu_unique);
      blas->solve("d'1[o][v]{" + mu_str + "} += " + shift);
      blas->solve("d'1[O][V]{" + mu_str + "} += " + shift);
      blas->solve("d'2[oo][vv]{" + mu_str + "} += " + shift);
      blas->solve("d'2[oO][vV]{" + mu_str + "} += " + shift);
      blas->solve("d'2[OO][VV]{" + mu_str + "} += " + shift);
    }
  }

  for(int i=0;i<moinfo->get_nunique();i++){
    int unique_i = moinfo->get_ref_number(i,UniqueRefs);
    string i_str = to_string(unique_i);
    // Form the coupling terms
    if(options_get_bool("COUPLING_TERMS")){
      for(int j=0;j<moinfo->get_nrefs();j++){
        int unique_j = moinfo->get_ref_number(j);
        string j_str = to_string(unique_j);

        double term = heff->get_right_eigenvector(j) * heff->get_right_eigenvector(unique_i) /
                      (std::pow(heff->get_right_eigenvector(unique_i),2.0) + std::pow(omega,2.0));

        if(fabs(term) > 100.0) {
          fprintf(outfile,"\n  Warning: c_nu/c_mu = %e ."
              "\n  1) turn on Tikhonow regularization or increase omega (TIKHONOW_OMEGA > 0)",term);
        }

        blas->set_scalar("factor_mk",unique_j,heff->get_matrix(unique_i,j) * term);
        if(unique_i!=j){
          if(j==unique_j){
            blas->solve("t1_eqns[o][v]{" + i_str + "} += factor_mk{" + j_str + "} t1[o][v]{" + j_str + "}");
            blas->solve("t1_eqns[O][V]{" + i_str + "} += factor_mk{" + j_str + "} t1[O][V]{" + j_str + "}");
          }else{
            blas->solve("t1_eqns[o][v]{" + i_str + "} += factor_mk{" + j_str + "} t1[O][V]{" + j_str + "}");
            blas->solve("t1_eqns[O][V]{" + i_str + "} += factor_mk{" + j_str + "} t1[o][v]{" + j_str + "}");
          }
        }
      }
    }

    // Update t1 for reference i
    blas->solve("t1_delta[o][v]{" + i_str + "}  =   t1_eqns[o][v]{" + i_str + "} / d'1[o][v]{" + i_str + "} - t1[o][v]{" + i_str + "}");
    blas->solve("t1_delta[O][V]{" + i_str + "}  =   t1_eqns[O][V]{" + i_str + "} / d'1[O][V]{" + i_str + "} - t1[O][V]{" + i_str + "}");

    blas->solve("t1[o][v]{" + i_str + "} = t1_eqns[o][v]{" + i_str + "} / d'1[o][v]{" + i_str + "}");
    blas->solve("t1[O][V]{" + i_str + "} = t1_eqns[O][V]{" + i_str + "} / d'1[O][V]{" + i_str + "}");

    zero_internal_amps();

    if(options_get_bool("COUPLING_TERMS")){
      // Add the contribution from the other references
      for(int j=0;j<moinfo->get_nrefs();j++){
        int unique_j = moinfo->get_ref_number(j);
        string j_str = to_string(unique_j);

        double term = heff->get_right_eigenvector(j) * heff->get_right_eigenvector(unique_i) /
                      (std::pow(heff->get_right_eigenvector(unique_i),2.0) + std::pow(omega,2.0));

        if(fabs(term) > 100.0) {
          fprintf(outfile,"\n  Warning: c_nu/c_mu = %e ."
              "\n  1) turn on Tikhonow regularization or increase omega (TIKHONOW_OMEGA > 0)",term);
        }

        blas->set_scalar("factor_mk",unique_j,heff->get_matrix(unique_i,j) * term);
        if(unique_i!=j){
          if(j==unique_j){
            // aaaa case
            // + t_ij^ab(nu/mu)
            blas->solve("Mk2[oo][vv]{" + i_str + "}  = t2[oo][vv]{" + j_str + "}");

            // P(ij)t_i^a(nu/mu)t_j^b(nu/mu)
            blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324#   t1[o][v]{" + j_str + "} X t1[o][v]{" + j_str + "}");
            blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314# - t1[o][v]{" + j_str + "} X t1[o][v]{" + j_str + "}");

            // -P(ij)P(ab)t_i^a(mu)t_j^b(nu/mu)
            blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324# - t1[o][v]{" + i_str + "} X t1[o][v]{" + j_str + "}");
            blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314#   t1[o][v]{" + i_str + "} X t1[o][v]{" + j_str + "}");
            blas->solve("Mk2[oo][vv]{" + i_str + "} += #1423#   t1[o][v]{" + i_str + "} X t1[o][v]{" + j_str + "}");
            blas->solve("Mk2[oo][vv]{" + i_str + "} += #2413# - t1[o][v]{" + i_str + "} X t1[o][v]{" + j_str + "}");

            // P(ij)t_i^a(mu)t_j^b(mu)
            blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324#   t1[o][v]{" + i_str + "} X t1[o][v]{" + i_str + "}");
            blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314# - t1[o][v]{" + i_str + "} X t1[o][v]{" + i_str + "}");

            blas->solve("t2_eqns[oo][vv]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[oo][vv]{" + i_str + "}");

            // abab case
            // + t_ij^ab(nu/mu)
            blas->solve("Mk2[oO][vV]{" + i_str + "}  = t2[oO][vV]{" + j_str + "}");

            // P(ij)t_i^a(nu/mu)t_J^B(nu/mu)
            blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324#   t1[o][v]{" + j_str + "} X t1[O][V]{" + j_str + "}");

            // -P(iJ)P(aB)t_i^a(mu)t_J^B(nu/mu)
            blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324# - t1[o][v]{" + i_str + "} X t1[O][V]{" + j_str + "}");
            blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324# - t1[o][v]{" + j_str + "} X t1[O][V]{" + i_str + "}");

            // P(iJ)t_i^a(mu)t_J^B(mu)
            blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324#   t1[o][v]{" + i_str + "} X t1[O][V]{" + i_str + "}");

            blas->solve("t2_eqns[oO][vV]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[oO][vV]{" + i_str + "}");

            // bbbb case
            // + t_ij^ab(nu/mu)
            blas->solve("Mk2[OO][VV]{" + i_str + "}  = t2[OO][VV]{" + j_str + "}");

            // P(ij)t_i^a(nu/mu)t_j^b(nu/mu)
            blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324#   t1[O][V]{" + j_str + "} X t1[O][V]{" + j_str + "}");
            blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314# - t1[O][V]{" + j_str + "} X t1[O][V]{" + j_str + "}");

            // -P(ij)P(ab)t_i^a(mu)t_j^b(nu/mu)
            blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324# - t1[O][V]{" + i_str + "} X t1[O][V]{" + j_str + "}");
            blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314#   t1[O][V]{" + i_str + "} X t1[O][V]{" + j_str + "}");
            blas->solve("Mk2[OO][VV]{" + i_str + "} += #1423#   t1[O][V]{" + i_str + "} X t1[O][V]{" + j_str + "}");
            blas->solve("Mk2[OO][VV]{" + i_str + "} += #2413# - t1[O][V]{" + i_str + "} X t1[O][V]{" + j_str + "}");

            // P(ij)t_i^a(mu)t_j^b(mu)
            blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324#   t1[O][V]{" + i_str + "} X t1[O][V]{" + i_str + "}");
            blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314# - t1[O][V]{" + i_str + "} X t1[O][V]{" + i_str + "}");

            blas->solve("t2_eqns[OO][VV]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[OO][VV]{" + i_str + "}");
          }else{
            // aaaa case
            // + t_ij^ab(nu/mu)
            blas->solve("Mk2[oo][vv]{" + i_str + "}  = t2[OO][VV]{" + j_str + "}");

            // P(ij)t_i^a(nu/mu)t_j^b(nu/mu)
            blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324#   t1[O][V]{" + j_str + "} X t1[O][V]{" + j_str + "}");
            blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314# - t1[O][V]{" + j_str + "} X t1[O][V]{" + j_str + "}");

            // -P(ij)P(ab)t_i^a(mu)t_j^b(nu/mu)
            blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324# - t1[o][v]{" + i_str + "} X t1[O][V]{" + j_str + "}");
            blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314#   t1[o][v]{" + i_str + "} X t1[O][V]{" + j_str + "}");
            blas->solve("Mk2[oo][vv]{" + i_str + "} += #1423#   t1[o][v]{" + i_str + "} X t1[O][V]{" + j_str + "}");
            blas->solve("Mk2[oo][vv]{" + i_str + "} += #2413# - t1[o][v]{" + i_str + "} X t1[O][V]{" + j_str + "}");

            // P(ij)t_i^a(mu)t_j^b(mu)
            blas->solve("Mk2[oo][vv]{" + i_str + "} += #1324#   t1[o][v]{" + i_str + "} X t1[o][v]{" + i_str + "}");
            blas->solve("Mk2[oo][vv]{" + i_str + "} += #2314# - t1[o][v]{" + i_str + "} X t1[o][v]{" + i_str + "}");

            blas->solve("t2_eqns[oo][vv]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[oo][vv]{" + i_str + "}");

            // abab case
            // + t_ij^ab(nu/mu)
            blas->solve("Mk2[oO][vV]{" + i_str + "}  = #2143# t2[oO][vV]{" + j_str + "}");

            // P(ij)t_i^a(nu/mu)t_J^B(nu/mu)
            blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324#   t1[O][V]{" + j_str + "} X t1[o][v]{" + j_str + "}");

            // -P(iJ)P(aB)t_i^a(mu)t_J^B(nu/mu)
            blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324# - t1[o][v]{" + i_str + "} X t1[o][v]{" + j_str + "}");
            blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324# - t1[O][V]{" + j_str + "} X t1[O][V]{" + i_str + "}");

            // P(iJ)t_i^a(mu)t_J^B(mu)
            blas->solve("Mk2[oO][vV]{" + i_str + "} += #1324#   t1[o][v]{" + i_str + "} X t1[O][V]{" + i_str + "}");

            blas->solve("t2_eqns[oO][vV]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[oO][vV]{" + i_str + "}");

            // bbbb case
            // + t_ij^ab(nu/mu)
            blas->solve("Mk2[OO][VV]{" + i_str + "}  = t2[oo][vv]{" + j_str + "}");

            // P(ij)t_i^a(nu/mu)t_j^b(nu/mu)
            blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324#   t1[o][v]{" + j_str + "} X t1[o][v]{" + j_str + "}");
            blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314# - t1[o][v]{" + j_str + "} X t1[o][v]{" + j_str + "}");

            // -P(ij)P(ab)t_i^a(mu)t_j^b(nu/mu)
            blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324# - t1[O][V]{" + i_str + "} X t1[o][v]{" + j_str + "}");
            blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314#   t1[O][V]{" + i_str + "} X t1[o][v]{" + j_str + "}");
            blas->solve("Mk2[OO][VV]{" + i_str + "} += #1423#   t1[O][V]{" + i_str + "} X t1[o][v]{" + j_str + "}");
            blas->solve("Mk2[OO][VV]{" + i_str + "} += #2413# - t1[O][V]{" + i_str + "} X t1[o][v]{" + j_str + "}");

            // P(ij)t_i^a(mu)t_j^b(mu)
            blas->solve("Mk2[OO][VV]{" + i_str + "} += #1324#   t1[O][V]{" + i_str + "} X t1[O][V]{" + i_str + "}");
            blas->solve("Mk2[OO][VV]{" + i_str + "} += #2314# - t1[O][V]{" + i_str + "} X t1[O][V]{" + i_str + "}");

            blas->solve("t2_eqns[OO][VV]{" + i_str + "} += factor_mk{" + j_str + "} Mk2[OO][VV]{" + i_str + "}");
          }

        }
      }
    }
    blas->solve("t2_delta[oo][vv]{" + i_str + "} = t2_eqns[oo][vv]{" + i_str + "} / d'2[oo][vv]{" + i_str + "} - t2[oo][vv]{" + i_str + "}");
    blas->solve("t2_delta[oO][vV]{" + i_str + "} = t2_eqns[oO][vV]{" + i_str + "} / d'2[oO][vV]{" + i_str + "} - t2[oO][vV]{" + i_str + "}");
    blas->solve("t2_delta[OO][VV]{" + i_str + "} = t2_eqns[OO][VV]{" + i_str + "} / d'2[OO][VV]{" + i_str + "} - t2[OO][VV]{" + i_str + "}");

    string damp = to_string(double(options_get_int("DAMPING_FACTOR"))/1000.0);
    string one_minus_damp = to_string(1.0-double(options_get_int("DAMPING_FACTOR"))/1000.0);
    blas->solve("t2[oo][vv]{" + i_str + "} = " + one_minus_damp + " t2_eqns[oo][vv]{" + i_str + "} / d'2[oo][vv]{" + i_str + "}");
    blas->solve("t2[oO][vV]{" + i_str + "} = " + one_minus_damp + " t2_eqns[oO][vV]{" + i_str + "} / d'2[oO][vV]{" + i_str + "}");
    blas->solve("t2[OO][VV]{" + i_str + "} = " + one_minus_damp + " t2_eqns[OO][VV]{" + i_str + "} / d'2[OO][VV]{" + i_str + "}");
    blas->solve("t2[oo][vv]{" + i_str + "} += " + damp + " t2_old[oo][vv]{" + i_str + "}");
    blas->solve("t2[oO][vV]{" + i_str + "} += " + damp + " t2_old[oO][vV]{" + i_str + "}");
    blas->solve("t2[OO][VV]{" + i_str + "} += " + damp + " t2_old[OO][VV]{" + i_str + "}");
    zero_internal_amps();
    blas->solve("t2_old[oo][vv]{" + i_str + "} = t2[oo][vv]{" + i_str + "}");
    blas->solve("t2_old[oO][vV]{" + i_str + "} = t2[oO][vV]{" + i_str + "}");
    blas->solve("t2_old[OO][VV]{" + i_str + "} = t2[OO][VV]{" + i_str + "}");
  }
}

}} /* End Namespaces */


//#include <psifiles.h>

  //#include <libpsio/psio.hpp>
  //#include <libutil/libutil.h>
  //

  //#include "debugging.h"
  //#include "index.h"
  //#include "mrcc.h"
  //#include "matrix.h"



