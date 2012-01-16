/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void get_params(Options &options)
{
  int errcod;
  std::string cachetype = "NULL";
  std::string ref;
  
  params.wfn = options.get_str("WFN");
  ref = options.get_str("REFERENCE");

  /* Default reference is RHF */
  params.ref = 0;
  params.semicanonical = 0;
  if(ref == "RHF") params.ref = 0;
  else if(ref == "ROHF" && params.wfn == "MP2") {
    params.ref = 2;
    params.semicanonical = 1;
  }
  else if(ref == "ROHF") params.ref = 1;
  else if(ref == "UHF") params.ref = 2;
  else {
    throw PsiException("Invalid Reference", __FILE__, __LINE__);
  }
  
  /* Default Jobtype */
  params.jobtype = options.get_str("JOBTYPE");

  /* Default Dertype */
  params.dertype = options.get_str("DERTYPE");

  params.gradient = options.get_str("DERTYPE") == "FIRST";
  params.relax_opdm = options.get_bool("OPDM_RELAX");
  params.opdm = options.get_bool("OPDM");
//  if(params.jobtype == "SP") {
//    params.opdm = 0;
//    params.relax_opdm = 0;
//    params.gradient = 0;
//  }
//  else if(params.jobtype == "OEPROP" && params.dertype == "NONE") {
//    params.opdm = 1;
//    params.relax_opdm = 0;
//    params.gradient = 0;
//  }
//  else if(params.jobtype == "OEPROP" && params.dertype == "FIRST") {
//    params.opdm = 1;
//    params.relax_opdm = 1;
//    params.gradient = 0;
//  }
//  else if(params.jobtype == "OPT" && params.dertype == "NONE") {
//    params.opdm = 0;
//    params.relax_opdm = 0;
//    params.gradient = 0;
//  }
//  else if(params.jobtype == "OPT" && params.dertype == "FIRST") {
//    params.opdm = 0;
//    params.relax_opdm = 0;
//    params.gradient = 1;
//  }
//  else if(params.jobtype == "OPT_FC" && params.dertype == "FIRST") {
//    params.opdm = 0;
//    params.relax_opdm = 0;
//    params.gradient = 1;
//  }
//  else if(params.jobtype == "SYMM_FC" && params.dertype == "FIRST") {
//    params.opdm = 0;
//    params.relax_opdm = 0;
//    params.gradient = 1;
//  }
//  else if(params.jobtype == "FREQ" && params.dertype == "NONE") {
//    params.opdm = 0;
//    params.relax_opdm = 0;
//    params.gradient = 0;
//  }
//  else if(params.jobtype == "FREQ" && params.dertype == "FIRST") {
//    params.opdm = 0;
//    params.relax_opdm = 0;
//    params.gradient = 1;
//  }
//  else {
//    throw PsiException("Invalid combination of JOBTYPE and DERTYPE", __FILE__, __LINE__);
//  }

  if((params.relax_opdm || params.gradient) && 
     (mo.nfzdocc != 0 || mo.nfzvirt != 0)) {
    throw PsiException("The Z-vector equations DO NOT work with frozen orbitals ... yet", __FILE__, __LINE__);
  }

  params.print = options.get_int("PRINT");
  
  params.cachelev = options.get_int("CACHELEVEL");
  
  cachetype = options.get_str("CACHETYPE");
  if(cachetype != "NULL") {
    if(cachetype == "LOW")
      params.cachetype = 1;
    else if(cachetype == "LRU")
      params.cachetype = 0;
    else {
      fprintf(outfile, "Invalide CACHETYPE = %s\n",cachetype.c_str());
      abort();
    }
  }
  
  /* get parameters related to SCS-MP2 or SCS-N-MP2 */
  /* see papers by S. Grimme or J. Platz */
  params.scs = options.get_int("SCS_N");
  if (params.scs == 1) {
    params.scs_scale_os = 0.0;
    params.scs_scale_ss = 1.76;
  }
  if(options.get_int("SCS") == 1) 
    params.scs = 1;
  if (params.scs == 1) { 
    params.scs_scale_os = options.get_double("SCALE_OS");
    params.scs_scale_ss = options.get_double("SCALE_SS");
  }

  params.memory = Process::environment.get_memory();
 
  fprintf(outfile, "\n");
  fprintf(outfile, "\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tWave function \t=\t%s\n", params.wfn.c_str());
  if(params.semicanonical) {
  fprintf(outfile, "\tReference WFN \t=\tROHF changed to UHF for Semicanonical Orbitals\n");
  }
  else {
  fprintf(outfile, "\tReference WFN \t=\t%s\n", (params.ref==0)?"RHF":((params.ref==1)?"ROHF":"UHF"));
  } 
  fprintf(outfile, "\tDerivative    \t=\t%s\n", params.dertype.c_str());
  fprintf(outfile, "\tCache Level   \t=\t%d\n", params.cachelev);
  fprintf(outfile, "\tCache Type    \t=\t%s\n", params.cachetype ? "LOW":"LRU");
  fprintf(outfile, "\tMemory (MB)   \t=\t%.1f\n",params.memory/1e6);
  fprintf(outfile, "\tPrint Level   \t=\t%d\n", params.print);
  fprintf(outfile, "\tOPDM          \t=\t%s\n", params.opdm ? "YES":"NO");
  fprintf(outfile, "\tSCS           \t=\t%s\n", (params.scs == 1) ? "True" : "False");
  fprintf(outfile, "\tSCALE_OS      \t=\t%.6f\n",params.scs_scale_os);
  fprintf(outfile, "\tSCALE_SS      \t=\t%.6f\n",params.scs_scale_ss);

  if (params.scs && params.dertype != "NONE") {
    fprintf(outfile,"\nWarning: SCS-MP2 computation requested but\n");
    fprintf(outfile,"derivative will be evaluated for standard MP2 energy.\n");
  }

}

}} /* End namespaces */
