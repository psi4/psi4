/* 
 *  SAPT0.CC 
 *
 */
#ifdef _MKL
#include <mkl.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>
#include <time.h>

#include <psifiles.h>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include <libmints/mints.h>

#include "sapt3bn7.h"

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace sapt {

SAPT3BN7::SAPT3BN7(Options& options, shared_ptr<PSIO> psio, 
    shared_ptr<Chkpt> chkpt) : SAPT3BN6(options, psio, chkpt)
{
}

SAPT3BN7::~SAPT3BN7()
{
}

double SAPT3BN7::compute_energy()
{  
    print_header(); 
    compute_integrals();
    compute_amplitudes();
    cphf_induction();

    exch100_s2();
    exch100_s3();
    exch100_s4();
    ind110();
    ind210();
    ind111();
    exch_ind200_s2();
    exch_ind110_s2();
    exch_ind200_s3();
    exch_ind110_s3();
    exch_disp200_s2();
    exch_disp110_s2();
    exch_disp200_s3();
    exch_disp110_s3();
    ind_disp210();
    disp111();
    disp3100();
    disp211_D();
    disp211_T();
    disp220_S();
    disp220_D();
    disp220_T();
    disp220_Q();

    return print_results();
}

void SAPT3BN7::print_header()
{
     fprintf(outfile,"       3B-SAPT  \n");
     fprintf(outfile,"    Ed Hohenstein\n") ;
     fprintf(outfile,"   20 October 2010\n") ;
     fprintf(outfile,"\n");
     fprintf(outfile,"    Orbital Information\n");
     fprintf(outfile,"  -----------------------\n");
     fprintf(outfile,"    NSO     = %9d\n",calc_info_.nso);
     fprintf(outfile,"    NMO     = %9d\n",calc_info_.nmo);
     fprintf(outfile,"    NRI     = %9d\n",calc_info_.nrio);
     fprintf(outfile,"    NOCC_A  = %9d\n",calc_info_.noccA);
     fprintf(outfile,"    NOCC_B  = %9d\n",calc_info_.noccB);
     fprintf(outfile,"    NOCC_C  = %9d\n",calc_info_.noccC);
     fprintf(outfile,"    NVIR_A  = %9d\n",calc_info_.nvirA);
     fprintf(outfile,"    NVIR_B  = %9d\n",calc_info_.nvirB);
     fprintf(outfile,"    NVIR_C  = %9d\n\n",calc_info_.nvirC);
     #ifdef _OPENMP
     fprintf(outfile,"Running SAPT with %d OMP threads\n\n",
       omp_get_max_threads());
     #endif 
     fprintf(outfile,"Notation: E_( V_AB V_AC V_BC ; W_A W_B W_C )\n\n");
     fflush(outfile);
}   

double SAPT3BN7::print_results()
{
  double exch_1_s2 = results_.exch100_s2 + results_.exch010_s2 + results_.exch001_s2;
  double exch_1_s3 = results_.exch100_s3 + results_.exch010_s3 + results_.exch001_s3;
  double exch_1_s4 = results_.exch100_s4 + results_.exch010_s4 + results_.exch001_s4;
  double ind2r = results_.ind110 + results_.ind101 + results_.ind011;
  double ind3 = results_.ind210 + results_.ind201 + results_.ind120 
    + results_.ind021 + results_.ind102 + results_.ind012 + results_.ind111;
  double exchind_s2 = results_.exch_ind200_s2 + results_.exch_ind020_s2 + results_.exch_ind002_s2 
    + results_.exch_ind110_s2 + results_.exch_ind101_s2 + results_.exch_ind011_s2;
  double exchind_s3 = results_.exch_ind200_s3 + results_.exch_ind020_s3 + results_.exch_ind002_s3 
    + results_.exch_ind110_s3 + results_.exch_ind101_s3 + results_.exch_ind011_s3;
  double dHF = results_.e_HF - exch_1_s2 - exch_1_s3 - exch_1_s4 - ind2r - ind3 - exchind_s2 - exchind_s3;
  double exchdisp_s2 = results_.exch_disp200_s2 + results_.exch_disp020_s2 + results_.exch_disp002_s2 
    + results_.exch_disp110_s2 + results_.exch_disp101_s2 + results_.exch_disp011_s2;
  double exchdisp_s3 = results_.exch_disp200_s3 + results_.exch_disp020_s3 + results_.exch_disp002_s3 
    + results_.exch_disp110_s3 + results_.exch_disp101_s3 + results_.exch_disp011_s3;
  double ind_disp = results_.ind_disp210 + results_.ind_disp201 + results_.ind_disp120 
    + results_.ind_disp021 + results_.ind_disp102 + results_.ind_disp012;

  double disp31 = results_.disp3100 + results_.disp3010 + results_.disp3001;
  double disp4sdq = results_.disp211d + results_.disp121d + results_.disp112d
    + results_.disp220s + results_.disp202s + results_.disp022s 
    + results_.disp220d + results_.disp202d + results_.disp022d 
    + results_.disp220q + results_.disp202q + results_.disp022q;

  double disp4t = results_.disp211t + results_.disp121t + results_.disp112t
    + results_.disp220t + results_.disp202t + results_.disp022t;

  double n5sapt = results_.e_HF + exchdisp_s2 + exchdisp_s3 + ind_disp + results_.disp111; 
  double n6sapt = n5sapt + disp31 + disp4sdq; 
  double n7sapt = n6sapt + disp4t; 

  fprintf(outfile,"    SAPT Results  \n");
  fprintf(outfile,"  ------------------------------------------------------------------------\n");
  fprintf(outfile,"    E_HF               %16.8lf mH %16.8lf kcal mol^-1\n",results_.e_HF*1000.0,results_.e_HF*627.5095);
  fprintf(outfile,"    Exch10 (S^2)       %16.8lf mH %16.8lf kcal mol^-1\n",exch_1_s2*1000.0,exch_1_s2*627.5095);
  fprintf(outfile,"    Exch10 (S^3)       %16.8lf mH %16.8lf kcal mol^-1\n",exch_1_s3*1000.0,exch_1_s3*627.5095);
  fprintf(outfile,"    Exch10 (S^4)       %16.8lf mH %16.8lf kcal mol^-1\n",exch_1_s4*1000.0,exch_1_s4*627.5095);
  fprintf(outfile,"    Ind20,r            %16.8lf mH %16.8lf kcal mol^-1\n",ind2r*1000.0,ind2r*627.5095);
  fprintf(outfile,"    Ind30              %16.8lf mH %16.8lf kcal mol^-1\n",ind3*1000.0,ind3*627.5095);
  fprintf(outfile,"    Exch-Ind20,r (S^2) %16.8lf mH %16.8lf kcal mol^-1\n",exchind_s2*1000.0,exchind_s2*627.5095);
  fprintf(outfile,"    Exch-Ind20,r (S^3) %16.8lf mH %16.8lf kcal mol^-1\n",exchind_s3*1000.0,exchind_s3*627.5095);
  fprintf(outfile,"    delta HF,r         %16.8lf mH %16.8lf kcal mol^-1\n",dHF*1000.0,dHF*627.5095);
  fprintf(outfile,"    Exch-Disp20 (S^2)  %16.8lf mH %16.8lf kcal mol^-1\n",exchdisp_s2*1000.0,exchdisp_s2*627.5095);
  fprintf(outfile,"    Exch-Disp20 (S^3)  %16.8lf mH %16.8lf kcal mol^-1\n",exchdisp_s3*1000.0,exchdisp_s3*627.5095);
  fprintf(outfile,"    Ind-Disp30         %16.8lf mH %16.8lf kcal mol^-1\n",ind_disp*1000.0,ind_disp*627.5095);
  fprintf(outfile,"    Disp30             %16.8lf mH %16.8lf kcal mol^-1\n",results_.disp111*1000.0,results_.disp111*627.5095);
  fprintf(outfile,"    Disp31             %16.8lf mH %16.8lf kcal mol^-1\n",disp31*1000.0,disp31*627.5095);
  fprintf(outfile,"    Disp40 (SDQ)       %16.8lf mH %16.8lf kcal mol^-1\n",disp4sdq*1000.0,disp4sdq*627.5095);
  fprintf(outfile,"    Disp40 (T)         %16.8lf mH %16.8lf kcal mol^-1\n\n",disp4t*1000.0,disp4t*627.5095);
  fprintf(outfile,"    3B-SAPT (N^5)      %16.8lf mH %16.8lf kcal mol^-1\n",n5sapt*1000.0,n5sapt*627.5095);
  fprintf(outfile,"    3B-SAPT (N^6)      %16.8lf mH %16.8lf kcal mol^-1\n",n6sapt*1000.0,n6sapt*627.5095);
  fprintf(outfile,"    3B-SAPT (N^7)      %16.8lf mH %16.8lf kcal mol^-1\n",n7sapt*1000.0,n7sapt*627.5095);

  double tot_exch = exch_1_s2 + exch_1_s3 + exch_1_s4;
  double tot_ind = ind2r + ind3 + exchind_s2 + exchind_s3 + dHF;
  double tot_disp = exchdisp_s2 + exchdisp_s3 + results_.disp111 + disp31
    + disp4sdq + disp4t;
  double tot_ind_disp = ind_disp;

  Process::environment.globals["SAPT EXCH ENERGY"] = tot_exch;
  Process::environment.globals["SAPT IND ENERGY"] = tot_ind;
  Process::environment.globals["SAPT DISP ENERGY"] = tot_disp;
  Process::environment.globals["SAPT IND-DISP ENERGY"] = tot_ind_disp;
  Process::environment.globals["SAPT 3B-SAPTN5 ENERGY"] = n5sapt;
  Process::environment.globals["SAPT 3B-SAPTN6 ENERGY"] = n6sapt;
  Process::environment.globals["SAPT 3B-SAPTN7 ENERGY"] = n7sapt;
  Process::environment.globals["SAPT ENERGY"] = n7sapt;

  return(n7sapt);
}

}}
