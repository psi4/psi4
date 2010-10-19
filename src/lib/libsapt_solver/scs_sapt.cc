/* 
 *  SAPT0.CC 
 *
 */

#include "scs_sapt.h"

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

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace sapt {

SCS_SAPT::SCS_SAPT(Options& options, shared_ptr<PSIO> psio, 
  shared_ptr<Chkpt> chkpt) : SAPT0(options, psio, chkpt)
{
}

SCS_SAPT::~SCS_SAPT()
{
}

double SCS_SAPT::compute_energy()
{  
    print_header(); 
    compute_integrals();
    compute_amplitudes();
    elst10();
    exch10();
    disp20();
    exch_disp20();
    cphf_induction();
    ind20();
    exch_ind20();
    scale_spin_components();

    return print_results();
}

void SCS_SAPT::print_header()
{
     fprintf(outfile,"      SCS-SAPT\n");
     fprintf(outfile,"    Ed Hohenstein\n");
     fprintf(outfile,"   16 October 2010\n");
     fprintf(outfile,"\n");
     fprintf(outfile,"    Orbital Information\n");
     fprintf(outfile,"  -----------------------\n");
     fprintf(outfile,"    NSO     = %9d\n",calc_info_.nso);
     fprintf(outfile,"    NMO     = %9d\n",calc_info_.nmo);
     fprintf(outfile,"    NRI     = %9d\n",calc_info_.nri);
     fprintf(outfile,"    NOCC_A  = %9d\n",calc_info_.noccA);
     fprintf(outfile,"    NOCC_B  = %9d\n",calc_info_.noccB);
     fprintf(outfile,"    NVIR_A  = %9d\n",calc_info_.nvirA);
     fprintf(outfile,"    NVIR_B  = %9d\n\n",calc_info_.nvirB);
    
     #ifdef _OPENMP
     fprintf(outfile,"Running SAPT with %d OMP threads\n\n",
       omp_get_max_threads());
     #endif 
    
     fflush(outfile);
}   

double SCS_SAPT::print_results()
{
  double eHF = calc_info_.eHF_D - calc_info_.eHF_A - calc_info_.eHF_B;
  double sapt0 = eHF + results_.disp20 + results_.exch_disp20;
  double dHF = eHF - (results_.elst10 + results_.exch10 + results_.ind20 + 
    results_.exch_ind20);
  double scs_disp = results_.disp20_os + results_.disp20_ss;
  double scs_exch_disp = results_.exch_disp20_os + results_.exch_disp20_ss;
  double scs_sapt = eHF + scs_disp + scs_exch_disp;

  fprintf(outfile,"    SAPT Results  \n");
  fprintf(outfile,"  ------------------------------------------------------------------\n");
  fprintf(outfile,"    E_HF             %16.8lf mH %16.8lf kcal mol^-1\n",
          eHF*1000.0,eHF*627.5095);
  fprintf(outfile,"    Elst10           %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.elst10*1000.0,results_.elst10*627.5095);
  fprintf(outfile,"    Exch10(S^2)      %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.exch10*1000.0,results_.exch10*627.5095);
  fprintf(outfile,"    Ind20,r          %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.ind20*1000.0,results_.ind20*627.5095);
  fprintf(outfile,"    Exch-Ind20,r     %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.exch_ind20*1000.0,results_.exch_ind20*627.5095);
  fprintf(outfile,"    delta HF,r       %16.8lf mH %16.8lf kcal mol^-1\n",
          dHF*1000.0,dHF*627.5095);
  fprintf(outfile,"    Disp20           %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.disp20*1000.0,results_.disp20*627.5095);
  fprintf(outfile,"    Exch-Disp20      %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.exch_disp20*1000.0,results_.exch_disp20*627.5095);
  fprintf(outfile,"    SCS-Disp20       %16.8lf mH %16.8lf kcal mol^-1\n",
          scs_disp*1000.0,scs_disp*627.5095);
  fprintf(outfile,"    SCS-Exch-Disp20  %16.8lf mH %16.8lf kcal mol^-1\n\n",
          scs_exch_disp*1000.0,scs_exch_disp*627.5095);
  fprintf(outfile,"    Total SAPT0      %16.8lf mH %16.8lf kcal mol^-1\n",
          sapt0*1000.0,sapt0*627.5095);
  fprintf(outfile,"    Total SCS-SAPT   %16.8lf mH %16.8lf kcal mol^-1\n",
          scs_sapt*1000.0,scs_sapt*627.5095);

  Process::environment.globals["SAPT ELST10 ENERGY"] = results_.elst10;
  Process::environment.globals["SAPT EXCH10 ENERGY"] = results_.exch10;
  Process::environment.globals["SAPT IND20 ENERGY"] = results_.ind20;
  Process::environment.globals["SAPT EXCH-IND20 ENERGY"] = results_.exch_ind20;
  Process::environment.globals["SAPT DELTA-HF ENERGY"] = dHF;
  Process::environment.globals["SAPT DISP20 ENERGY"] = results_.disp20;
  Process::environment.globals["SAPT EXCH-DISP20 ENERGY"] = 
    results_.exch_disp20;
  Process::environment.globals["SAPT SCS-DISP20 ENERGY"] = scs_disp;
  Process::environment.globals["SAPT SCS-EXCH-DISP20 ENERGY"] = scs_exch_disp;
  Process::environment.globals["SAPT SAPT0 ENERGY"] = sapt0;
  Process::environment.globals["SAPT SCS-SAPT ENERGY"] = scs_sapt;
  Process::environment.globals["SAPT ENERGY"] = scs_sapt;

  return(scs_sapt);
}

void SCS_SAPT::scale_spin_components()
{
  results_.disp20_os = 0.5*results_.disp20;
  results_.disp20_ss = 0.5*results_.disp20;

  double OS_SAPT = results_.disp20_os+results_.exch_disp20_os;
  double SS_SAPT = results_.disp20_ss+results_.exch_disp20_ss;

  fprintf(outfile,"    SCS-SAPT Decomposition  \n");
  fprintf(outfile,"  ------------------------------------------------------------------\n");
  fprintf(outfile,"    Disp20 (OS)          %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.disp20_os*1000.0,results_.disp20_os*627.5095);
  fprintf(outfile,"    Disp20 (SS)          %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.disp20_ss*1000.0,results_.disp20_ss*627.5095);
  fprintf(outfile,"    Exch-Disp20 (OS)     %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.exch_disp20_os*1000.0,results_.exch_disp20_os*627.5095);
  fprintf(outfile,"    Exch-Disp20 (SS)     %16.8lf mH %16.8lf kcal mol^-1\n\n",
          results_.exch_disp20_ss*1000.0,results_.exch_disp20_ss*627.5095);
  fprintf(outfile,"    Total OS SAPT        %16.8lf mH %16.8lf kcal mol^-1\n",
          OS_SAPT*1000.0,OS_SAPT*627.5095);
  fprintf(outfile,"    Total SS SAPT        %16.8lf mH %16.8lf kcal mol^-1\n\n",
          SS_SAPT*1000.0,SS_SAPT*627.5095);

  double sOS = 1.0;
  double sSS = 1.0;

  results_.disp20_os *= sOS;
  results_.exch_disp20_os *= sOS;
  results_.disp20_ss *= sSS;
  results_.exch_disp20_ss *= sSS;

  OS_SAPT = results_.disp20_os+results_.exch_disp20_os;
  SS_SAPT = results_.disp20_ss+results_.exch_disp20_ss;

  fprintf(outfile,"    SCS-Disp20 (OS)      %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.disp20_os*1000.0,results_.disp20_os*627.5095);
  fprintf(outfile,"    SCS-Disp20 (SS)      %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.disp20_ss*1000.0,results_.disp20_ss*627.5095);
  fprintf(outfile,"    SCS-Exch-Disp20 (OS) %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.exch_disp20_os*1000.0,results_.exch_disp20_os*627.5095);
  fprintf(outfile,"    SCS-Exch-Disp20 (SS) %16.8lf mH %16.8lf kcal mol^-1\n\n",
          results_.exch_disp20_ss*1000.0,results_.exch_disp20_ss*627.5095);
  fprintf(outfile,"    Total OS SCS-SAPT    %16.8lf mH %16.8lf kcal mol^-1\n",
          OS_SAPT*1000.0,OS_SAPT*627.5095);
  fprintf(outfile,"    Total SS SCS-SAPT    %16.8lf mH %16.8lf kcal mol^-1\n\n",
          SS_SAPT*1000.0,SS_SAPT*627.5095);
}

}}
