/* 
 *  SAPT0.CC 
 *
 */

#ifdef HAVEMKL
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

#include "sapt2.h"

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace sapt {

SAPT2::SAPT2(Options& options, shared_ptr<PSIO> psio, 
  shared_ptr<Chkpt> chkpt) : SAPT0(options, psio, chkpt)
{
}

SAPT2::~SAPT2()
{
}

double SAPT2::compute_energy()
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
    elst12();
    exch11();
    exch12();
    ind22();

    return print_results();
}

void SAPT2::print_header()
{
     fprintf(outfile,"        SAPT2  \n");
     fprintf(outfile,"    Ed Hohenstein\n");
     fprintf(outfile,"     6 June 2009\n");
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

double SAPT2::print_results()
{
  double eHF = calc_info_.eHF_D - calc_info_.eHF_A - calc_info_.eHF_B;
  double exch_ind22 = results_.ind22*results_.exch_ind20/results_.ind20;
  double sapt0 = eHF + results_.disp20 + results_.exch_disp20;
  double sapt2 = eHF + results_.disp20 + results_.exch_disp20 + 
    results_.elst12 + results_.exch11 + results_.exch12 + results_.ind22 +
    exch_ind22;
  double dHF = eHF - (results_.elst10 + results_.exch10 + results_.ind20 + 
    results_.exch_ind20);

  fprintf(outfile,"    SAPT Results  \n");
  fprintf(outfile,"  ------------------------------------------------------------------\n");
  fprintf(outfile,"    E_HF             %16.8lf mH %16.8lf kcal mol^-1\n",
          eHF*1000.0,eHF*627.5095);
  fprintf(outfile,"    Elst10           %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.elst10*1000.0,results_.elst10*627.5095);
  fprintf(outfile,"    Elst12,r         %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.elst12*1000.0,results_.elst12*627.5095);
  fprintf(outfile,"    Exch10(S^2)      %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.exch10*1000.0,results_.exch10*627.5095);
  fprintf(outfile,"    Exch11(S^2)      %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.exch11*1000.0,results_.exch11*627.5095);
  fprintf(outfile,"    Exch12(S^2)      %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.exch12*1000.0,results_.exch12*627.5095);
  fprintf(outfile,"    Ind20,r          %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.ind20*1000.0,results_.ind20*627.5095);
  fprintf(outfile,"    Ind22            %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.ind22*1000.0,results_.ind22*627.5095);
  fprintf(outfile,"    Exch-Ind20,r     %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.exch_ind20*1000.0,results_.exch_ind20*627.5095);
  fprintf(outfile,"    Exch-Ind22       %16.8lf mH %16.8lf kcal mol^-1\n",
          exch_ind22*1000.0,exch_ind22*627.5095);
  fprintf(outfile,"    delta HF,r       %16.8lf mH %16.8lf kcal mol^-1\n",
          dHF*1000.0,dHF*627.5095);
  fprintf(outfile,"    Disp20           %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.disp20*1000.0,results_.disp20*627.5095);
  fprintf(outfile,"    Exch-Disp20      %16.8lf mH %16.8lf kcal mol^-1\n\n",
          results_.exch_disp20*1000.0,results_.exch_disp20*627.5095);
  fprintf(outfile,"    Total SAPT0      %16.8lf mH %16.8lf kcal mol^-1\n",
          sapt0*1000.0,sapt0*627.5095);
  fprintf(outfile,"    Total SAPT2      %16.8lf mH %16.8lf kcal mol^-1\n",
          sapt2*1000.0,sapt2*627.5095);

  Process::environment.globals["SAPT ELST10 ENERGY"] = results_.elst10;
  Process::environment.globals["SAPT ELST12 ENERGY"] = results_.elst12;
  Process::environment.globals["SAPT EXCH10 ENERGY"] = results_.exch10;
  Process::environment.globals["SAPT EXCH11 ENERGY"] = results_.exch11;
  Process::environment.globals["SAPT EXCH12 ENERGY"] = results_.exch12;
  Process::environment.globals["SAPT IND20 ENERGY"] = results_.ind20;
  Process::environment.globals["SAPT EXCH-IND20 ENERGY"] = results_.exch_ind20;
  Process::environment.globals["SAPT IND22 ENERGY"] = results_.ind22;
  Process::environment.globals["SAPT EXCH-IND22 ENERGY"] = exch_ind22;
  Process::environment.globals["SAPT DELTA-HF ENERGY"] = dHF;
  Process::environment.globals["SAPT DISP20 ENERGY"] = results_.disp20;
  Process::environment.globals["SAPT EXCH-DISP20 ENERGY"] = 
    results_.exch_disp20;
  Process::environment.globals["SAPT SAPT0 ENERGY"] = sapt0;
  Process::environment.globals["SAPT SAPT2 ENERGY"] = sapt2;
  Process::environment.globals["SAPT ENERGY"] = sapt2;

  return(sapt2);
}

}}
