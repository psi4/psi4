#include <libmoinfo/libmoinfo.h>
#include <liboptions/liboptions.h>
#include "mrcc.h"
#include "matrix.h"
#include "blas.h"
#include "debugging.h"
#include <libutil/libutil.h>
#include <psi4-dec.h>

#include <boost/shared_ptr.hpp>
#include <libchkpt/chkpt.hpp>

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;

using namespace std;

bool CCMRCC::build_diagonalize_Heff(int cycle, double time)
{
  total_time     = time;
  bool converged = false;

  build_Heff_diagonal();
  build_Heff_offdiagonal();

  if(cycle==0){
    current_energy=diagonalize_Heff(moinfo->get_root(),moinfo->get_nrefs(),Heff,right_eigenvector,left_eigenvector,true);
    old_energy=current_energy;
    print_eigensystem(moinfo->get_nrefs(),Heff,right_eigenvector);
  }
  if((cycle>0) || (cycle==-1)){
    // Compute the energy difference
    old_energy=current_energy;
    current_energy=diagonalize_Heff(moinfo->get_root(),moinfo->get_nrefs(),Heff,right_eigenvector,left_eigenvector,false);

    if(options_.get_bool("HEFF_PRINT"))
      print_eigensystem(moinfo->get_nrefs(),Heff,right_eigenvector);
    DEBUGGING(3,
      print_eigensystem(moinfo->get_nrefs(),Heff,right_eigenvector);
    )
    double delta_energy = current_energy-old_energy;
    if(fabs(delta_energy) < options_.get_double("CONVERGENCE"))
      converged = true;


///    TODO fix this code which is temporarly not working
//     if(options_get_int("DAMPING_FACTOR")>0){
//       if(fabs(delta_energy) < moinfo->get_no_damp_convergence()){
//         double damping_factor = moinfo->get_damping_factor();
//         damping_factor *= 0.95;
//         fprintf(outfile,"\n\t# Scaling damp factor to zero, damping_factor = %lf",moinfo->get_damping_factor());
//         moinfo->set_damping_factor(damping_factor);
//       }
//     }
  }
  print_mrccsd_energy(cycle);
  if(converged){
    print_eigensystem(moinfo->get_nrefs(),Heff,right_eigenvector);
    Process::environment.globals["CURRENT ENERGY"]    = current_energy;
    Process::environment.globals["MRCC TOTAL ENERGY"] = current_energy;
  }
  return(converged);
}

void CCMRCC::build_Heff_diagonal()
{
  // Compute the diagonal elements of the effective Hamiltonian
  // using a simple UCCSD energy expression
  blas->solve("Eaa{u}   = t1[o][v]{u} . fock[o][v]{u}");
  blas->solve("Ebb{u}   = t1[O][V]{u} . fock[O][V]{u}");

  blas->solve("Eaaaa{u} = 1/4 tau[oo][vv]{u} . <[oo]:[vv]>");
  blas->solve("Eabab{u} =     tau[oO][vV]{u} . <[oo]|[vv]>");
  blas->solve("Ebbbb{u} = 1/4 tau[OO][VV]{u} . <[oo]:[vv]>");

  blas->solve("ECCSD{u}  = Eaa{u} + Ebb{u} + Eaaaa{u} + Eabab{u} + Ebbbb{u} + ERef{u}");

  for(int n=0;n<moinfo->get_nrefs();n++)
    Heff[n][n] = blas->get_scalar("ECCSD",moinfo->get_ref_number(n));
}

void CCMRCC::build_Heff_offdiagonal()
{
  for(int i=0;i<moinfo->get_ref_size(AllRefs);i++){
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
            element=sign_internal_excitation * blas->get_MatTmp("t1_eqns[o][v]",i_unique,none)->get_two_address_element(
                                               alpha_internal_excitation[0].first,
                                               alpha_internal_excitation[0].second);

          // Set beta-beta single excitations
          if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==1))
            element=sign_internal_excitation * blas->get_MatTmp("t1_eqns[O][V]",i_unique,none)->get_two_address_element(
                                               beta_internal_excitation[0].first,
                                               beta_internal_excitation[0].second);

          // Set (alpha,alpha)->(alpha,alpha) double excitations
          if((alpha_internal_excitation.size()==2)&&(beta_internal_excitation.size()==0))
            element=sign_internal_excitation * blas->get_MatTmp("t2_eqns[oo][vv]",i_unique,none)->get_four_address_element(
                                               alpha_internal_excitation[0].first,
                                               alpha_internal_excitation[1].first,
                                               alpha_internal_excitation[0].second,
                                               alpha_internal_excitation[1].second);

          // Set (alpha,beta)->(alpha,beta) double excitations
          if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==1))
            element=sign_internal_excitation * blas->get_MatTmp("t2_eqns[oO][vV]",i_unique,none)->get_four_address_element(
                                               alpha_internal_excitation[0].first,
                                                beta_internal_excitation[0].first,
                                               alpha_internal_excitation[0].second,
                                                beta_internal_excitation[0].second);

          // Set (beta,beta)->(beta,beta) double excitations
          if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==2))
            element=sign_internal_excitation * blas->get_MatTmp("t2_eqns[OO][VV]",i_unique,none)->get_four_address_element(
                                               beta_internal_excitation[0].first,
                                               beta_internal_excitation[1].first,
                                               beta_internal_excitation[0].second,
                                               beta_internal_excitation[1].second);
        }else{
          // Set alpha-alpha single excitations
          if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==0))
            element=sign_internal_excitation * blas->get_MatTmp("t1_eqns[O][V]",i_unique,none)->get_two_address_element(
                                               alpha_internal_excitation[0].first,
                                               alpha_internal_excitation[0].second);

          // Set beta-beta single excitations
          if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==1))
            element=sign_internal_excitation * blas->get_MatTmp("t1_eqns[o][v]",i_unique,none)->get_two_address_element(
                                               beta_internal_excitation[0].first,
                                               beta_internal_excitation[0].second);

          // Set (alpha,alpha)->(alpha,alpha) double excitations
          if((alpha_internal_excitation.size()==2)&&(beta_internal_excitation.size()==0))
            element=sign_internal_excitation * blas->get_MatTmp("t2_eqns[OO][VV]",i_unique,none)->get_four_address_element(
                                               alpha_internal_excitation[0].first,
                                               alpha_internal_excitation[1].first,
                                               alpha_internal_excitation[0].second,
                                               alpha_internal_excitation[1].second);

          // Set (alpha,beta)->(alpha,beta) double excitations
          if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==1))
            element=sign_internal_excitation * blas->get_MatTmp("t2_eqns[oO][vV]",i_unique,none)->get_four_address_element(
                                                beta_internal_excitation[0].first,
                                               alpha_internal_excitation[0].first,
                                                beta_internal_excitation[0].second,
                                               alpha_internal_excitation[0].second);

          // Set (beta,beta)->(beta,beta) double excitations
          if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==2))
            element=sign_internal_excitation * blas->get_MatTmp("t2_eqns[oo][vv]",i_unique,none)->get_four_address_element(
                                               beta_internal_excitation[0].first,
                                               beta_internal_excitation[1].first,
                                               beta_internal_excitation[0].second,
                                               beta_internal_excitation[1].second);
        }
        Heff[j][i]=element;
      }
    }
  }
}

}} /* End Namespaces */
