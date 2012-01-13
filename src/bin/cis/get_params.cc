/*! \file
    \ingroup CIS
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cmath>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cis {

void get_params(Options & options)
{
  int h, i, j, errcod, ref, iconv, max_dim;
  char *cachetype = NULL;
  std::string junk;

  params.wfn = options.get_str("WFN"); //default 0? options CCSD CCSD_T EOM_CCSD CIS

  junk = options.get_str("REFERENCE"); //default RHF options RHF ==> 0, ROHF ==> 1, UHF ==> 2
  if(junk == "RHF") ref = 0;
  else if(junk == "ROHF") ref = 1;
  else if(junk == "UHF") ref = 2;

  /* Make sure the value of ref matches that from CC_INFO */
  if(params.ref != ref) {
    printf("Value of REFERENCE from input.dat (%1d) and CC_INFO (%1d) do not match!\n", 
	   ref, params.ref);
    throw PsiException("cis input comparison error REFERENCE and params.ref", __FILE__, __LINE__);
  }

  local.amp_print_cutoff = options.get_double("LOCAL_AMPS_PRINT_CUTOFF"); //default 0.60
  params.print = options.get_int("PRINT"); //default 0
  params.maxiter = options.get_int("MAXITER"); //default 500

  params.convergence = options.get_double("R_CONVERGENCE"); //default 7 ==> 1e-7

  if (options["STATES_PER_IRREP"].size() > 0) {
    i = options["STATES_PER_IRREP"].size();
    if (i != moinfo.nirreps) {
      fprintf(outfile,"Dim. of states_per_irrep vector must be %d\n", moinfo.nirreps);
      throw PsiException("cis input comparison error STATES_PER_IRREP and moinfo.nirreps", __FILE__, __LINE__);
    }
    // possibly replace the following
    //for (i=0; i<moinfo.irreps; ++i)
    //  params.rpi[i] = options["STATES_PER_IRREP"][i].to_integer();
    // with the following line
    params.rpi = options.get_int_array("STATES_PER_IRREP");
  }
  else {
    params.rpi = init_int_array(moinfo.nirreps);
    for (i=0; i<moinfo.nirreps; i++) {
      params.rpi[i] = 2;
    }
  } 

/* start pre- porting
  params.rpi = (int *) malloc(moinfo.nirreps * sizeof(int));
  if (ip_exist("STATES_PER_IRREP",0)) {  
    ip_count("STATES_PER_IRREP", &i, 0);
    if (i != moinfo.nirreps) {
      fprintf(outfile,"Dim. of states_per_irrep vector must be %d\n", 
        moinfo.nirreps) ;
      throw PsiException("cis input comparison error STATES_PER_IRREP and moinfo.nirreps", __FILE__, __LINE__);
    }
    for (i=0;i<moinfo.nirreps;++i)
      errcod = ip_data("STATES_PER_IRREP","%d",&(params.rpi[i]),1,i);
  }
  else { 
    params.rpi = init_int_array(moinfo.nirreps);
    for (i=0;i<moinfo.nirreps;i++) {
      params.rpi[i] = 2;
    }
  } 
end pre-porting */

  /* Test STATES_PER_IRREP vector for correct dimensions */
  for(h=0; h < moinfo.nirreps; h++) {
    max_dim = 0;
    for(i=0; i < moinfo.nirreps; i++) {
      j = i^h;
      if(params.ref == 0) max_dim += moinfo.occpi[i] * moinfo.virtpi[j];
      else max_dim += moinfo.aoccpi[i] * moinfo.avirtpi[j] + moinfo.boccpi[i] * moinfo.bvirtpi[j];
    }
    if(params.rpi[h] > max_dim) {
      fprintf(outfile, "\n\tRequested no. of roots for irrep %d exceeds basis limit.\n", h);
      fprintf(outfile, "\tSetting rpi[%d] = %d.\n", h, max_dim);
      params.rpi[h] = max_dim;
    }
  }

  params.diag_method = options.get_str("DIAG_METHOD");  //default DAVIDSON options DAVIDSON FULL
  params.local = options.get_bool("LOCAL"); //default 0
  local.cutoff = options.get_double("LOCAL_CUTOFF"); //default 0.02
  local.method = options.get_str("LOCAL_METHOD"); //default WERNER options AOBASIS WERNER
  local.weakp = options.get_str("LOCAL_WEAKP"); //default MP2 options MP2 NEGLECT NONE
  local.ghost = options.get_int("LOCAL_GHOST"); //default -1

  params.memory = module.get_memory();
//  fndcor(&(params.memory),infile,outfile);

  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tWave function          = %6s\n", params.wfn.c_str());
  fprintf(outfile, "\tReference wfn          = %5s\n",
	  (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
  fprintf(outfile, "\tMemory (Mbytes)        = %5.1f\n",params.memory/1e6);
  fprintf(outfile, "\tRoots sought per irrep =  ");
  for (i=0;i<moinfo.nirreps;++i) fprintf(outfile, "%3d", params.rpi[i]);
  fprintf(outfile,"\n");
  fprintf(outfile, "\tDiagonalization        = %s\n", params.diag_method.c_str());
  if(params.diag_method == "DAVIDSON") {
    fprintf(outfile, "\tMaxiter                = %4d\n", params.maxiter);
    fprintf(outfile, "\tConvergence            = %3.1e\n", params.convergence);
  }
  fprintf(outfile, "\tPrint Level            =    %1d\n",  params.print);
  fprintf(outfile, "\tLocal CC        =     %s\n", params.local ? "Yes" : "No");
  if(params.local) {
    fprintf(outfile, "\tLocal Cutoff    = %3.1e\n", local.cutoff);
    fprintf(outfile, "\tLocal Method    =    %s\n", local.method.c_str());
    fprintf(outfile, "\tWeak pairs      =    %s\n", local.weakp.c_str());
    fprintf(outfile, "\tGhost atom      =    %d\n", local.ghost);
    fprintf(outfile, "\tAmp Print Cutoff =   %f\n", local.amp_print_cutoff);
  }
  fprintf(outfile, "\n");

}


}} // namespace psi::cis
