/*! \file
    \ingroup LMP2
    \brief Enter brief description of file here
*/
//#include "mpi.h"
//#include <iostream>
//#include <fstream>              // file I/O support
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
#include <libparallel/parallel.h>
#include <cmath>
//#include <libipv1/ip_lib.h>
//#include <libciomr/libciomr.h>
//#include <libchkpt/chkpt.hpp>
//#include <libpsio/psio.hpp>
//#include <libqt/qt.h>
//#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{

namespace lmp2{

void LMP2::get_params(Options &options) {

  int errcod, iconv, rconv, fs;
  long int max_bytes;
  char *cachetype = NULL;
  char *junk;


  /* Default reference is RHF */
  if(options.get_str("REFERENCE") == "RHF") ref = 0;

  print = options.get_int("PRINT");

  maxiter = options.get_int("MAXITER");

  iconv = options.get_int("ENERGY_CONV");
  econv = 1.0*pow(10.0,(double) -iconv);

  iconv = options.get_int("SCREENING");
  escreen = 1.0*pow(10, (double) -iconv);

  screen_int = options.get_bool("SCREEN_INTS");

  rconv = options.get_int("RMS_CONV");
  rmsconv = 1.0*pow(10.0,(double) -rconv);

  fs = options.get_int("FSKIP");
  fskip = 1.0*pow(10.0,(double) -fs);

  diis = options.get_bool("USE_DIIS");

  it_diis = options.get_int("DIISSTART");
  if(it_diis < 3) {
    if(Communicator::world->me() == 0) {
      fprintf(outfile, "\n\t*** WARNING ***\n");
      fprintf(outfile, "\tDIISSTART can not be less than 3\n");
      fprintf(outfile, "\tReseting DIISSTART to 3\n");
      fprintf(outfile, "\t***************\n");
    }
    it_diis = 3;
  }
  ndiis = options.get_int("NDIIS");

  memory = options.get_int("MEMORY");
  wfn = const_cast<char*>(options.get_cstr("WFN"));

  cutoff = options.get_double("LOCAL_CUTOFF");
  neglectdp = options.get_bool("NEGLECT_DP");
  dpcutoff = options.get_double("DISTANT_PAIR");

  schwarz_tol = options.get_int("SCHWARTZ_TOL");
  tol = 1.0*pow(10.0,(double) -schwarz_tol);

  RILMP2 = 0;
  if(options.get_bool("RI_LMP2"))
    RILMP2 = 1;

  scs_scale_os = options.get_double("SCALE_OS");
  scs_scale_ss = options.get_double("SCALE_SS");

  if(Communicator::world->me() == 0)
    print_params();
}

void LMP2::print_params(){

    fprintf(outfile, "\n");
    fprintf(outfile, "\tInput parameters:\n");
    fprintf(outfile, "\t-----------------\n");
    fprintf(outfile, "\tWave function \t\t= %s\n", wfn);
    fprintf(outfile, "\tReference WFN \t\t= %s\n", (ref==0)?"RHF":((ref==1)?"ROHF":"UHF"));
    fprintf(outfile, "\tMemory (MB)   \t\t= %.1f\n",memory/1e6);
    fprintf(outfile, "\tNum Procs     \t\t= %d\n", Communicator::world->nproc());
    fprintf(outfile, "\tMaxiter       \t\t= %d\n", maxiter);
    fprintf(outfile, "\tEnergy Convergence   \t= %3.1e\n", econv);
    fprintf(outfile, "\tRMS Convergence   \t= %3.1e\n", rmsconv);
    fprintf(outfile, "\tIntegral Prescreen \t= %3.1e\n", escreen);
    fprintf(outfile, "\tF-Skip   \t\t= %3.1e\n", fskip);
    fprintf(outfile, "\tPrint Level   \t\t= %d\n", print);
    fprintf(outfile, "\tLocal Cutoff  \t\t= %3.1e\n", cutoff);
    fprintf(outfile, "\tUse RI Approximation \t= %s\n", RILMP2 ? "Yes" : "No");
    fprintf(outfile, "\tNeglect Distant Pairs \t= %s\n", neglectdp ? "Yes" : "No");
    if(neglectdp)
        fprintf(outfile, "\tDistace Pair Cutoff \t= %2.1f\n", dpcutoff);
    fprintf(outfile,"\tOpposite-spin scaled by %10.4lf\n",scs_scale_os);
    fprintf(outfile,"\tSame-spin scaled by     %10.4lf\n",scs_scale_ss);
    fprintf(outfile, "\n\tUse DIIS extrapolation = %s\n", diis ? "Yes" : "No");
    if(diis == 1) {
      fprintf(outfile, "\titerations before DIIS extrapolation = %d\n", it_diis);
      fprintf(outfile, "\t%d error matrices will be kept\n", ndiis);
    }

}

}} /* End namespace */

