/*! \file
    \ingroup CIS
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cis {

void get_params()
{
  int h, i, j, errcod, ref, iconv, max_dim;
  char *cachetype = NULL;
  char *junk;

  errcod = ip_string("WFN", &(params.wfn), 0);
  if(strcmp(params.wfn, "CCSD") && strcmp(params.wfn, "CCSD_T") &&
     strcmp(params.wfn, "EOM_CCSD") && strcmp(params.wfn, "CIS")) {
    fprintf(outfile, "Invalid value of input keyword WFN: %s\n", params.wfn);
    exit(PSI_RETURN_FAILURE);
  }

  errcod = ip_string("REFERENCE", &(junk),0);
  /* if no reference is given, assume rhf */
  if (errcod != IPE_OK) {
    ref = 0;
  }
  else {
    if(!strcmp(junk, "RHF")) ref = 0;
    else if(!strcmp(junk, "ROHF")) ref = 1;
    else if(!strcmp(junk, "UHF")) ref = 2;
    else { 
      printf("Invalid value of input keyword REFERENCE: %s\n", junk);
      exit(PSI_RETURN_FAILURE); 
    }
    free(junk);
  }

  /* Make sure the value of ref matches that from CC_INFO */
  if(params.ref != ref) {
    printf("Value of REFERENCE from input.dat (%1d) and CC_INFO (%1d) do not match!\n", 
	   ref, params.ref);
    exit(PSI_RETURN_FAILURE);
  }

  local.amp_print_cutoff = 0.60;
  errcod = ip_data("LOCAL_AMP_PRINT_CUTOFF","%lf",&(local.amp_print_cutoff),0);

  params.print = 0;
  errcod = ip_data("PRINT", "%d", &(params.print),0);

  params.maxiter = 500;
  errcod = ip_data("MAXITER","%d",&(params.maxiter),0);
  params.convergence = 1e-7;
  errcod = ip_data("CONVERGENCE","%d",&(iconv),0);
  if(errcod == IPE_OK) params.convergence = 1.0*pow(10.0,(double) -iconv);

  params.rpi = (int *) malloc(moinfo.nirreps * sizeof(int));
  if (ip_exist("STATES_PER_IRREP",0)) {
    ip_count("STATES_PER_IRREP", &i, 0);
    if (i != moinfo.nirreps) {
      fprintf(outfile,"Dim. of states_per_irrep vector must be %d\n", 
        moinfo.nirreps) ;
      exit(PSI_RETURN_FAILURE);
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

  if(ip_exist("DIAG_METHOD",0)) {
    errcod = ip_string("DIAG_METHOD", &(params.diag_method), 0);
    if(strcmp(params.diag_method,"DAVIDSON") && strcmp(params.diag_method,"FULL")) {
      fprintf(outfile, "Invalid diagonalization method requested: %s\n", params.diag_method);
      exit(PSI_RETURN_FAILURE);
    }
  }
  else {
    params.diag_method = (char *) malloc(9 * sizeof(char));
    sprintf(params.diag_method, "%s", "DAVIDSON");
  }

  params.local = 0;
  errcod = ip_boolean("LOCAL", &(params.local),0);
  local.cutoff = 0.02;
  errcod = ip_data("LOCAL_CUTOFF", "%lf", &(local.cutoff), 0);

  if(ip_exist("LOCAL_METHOD",0)) {
    errcod = ip_string("LOCAL_METHOD", &(local.method), 0);
    if(strcmp(local.method,"AOBASIS") && strcmp(local.method,"WERNER")) {
      fprintf(outfile, "Invalid local correlation method: %s\n", local.method);
      exit(PSI_RETURN_FAILURE);
    }
  }
  else if(params.local) {
    local.method = (char *) malloc(7 * sizeof(char));
    sprintf(local.method, "%s", "WERNER");
  }

  if(ip_exist("LOCAL_WEAKP",0)) {
    errcod = ip_string("LOCAL_WEAKP", &(local.weakp), 0);
    if(strcmp(local.weakp,"MP2") && strcmp(local.weakp,"NEGLECT") && strcmp(local.weakp,"NONE")) {
      fprintf(outfile, "Invalid method for treating local pairs: %s\n", local.weakp);
      exit(PSI_RETURN_FAILURE);
    }
  }
  else if(params.local) {
    local.weakp = (char *) malloc(4 * sizeof(char));
    sprintf(local.weakp, "%s", "MP2");
  }

  local.ghost = -1;
  if(ip_exist("LOCAL_GHOST",0))
    errcod = ip_data("LOCAL_GHOST", "%d", &(local.ghost), 0);

  fndcor(&(params.memory),infile,outfile);

  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tWave function          = %6s\n", params.wfn);
  fprintf(outfile, "\tReference wfn          = %5s\n",
	  (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
  fprintf(outfile, "\tMemory (Mbytes)        = %5.1f\n",params.memory/1e6);
  fprintf(outfile, "\tRoots sought per irrep =  ");
  for (i=0;i<moinfo.nirreps;++i) fprintf(outfile, "%3d", params.rpi[i]);
  fprintf(outfile,"\n");
  fprintf(outfile, "\tDiagonalization        = %s\n", params.diag_method);
  if(!strcmp(params.diag_method,"DAVIDSON")) {
    fprintf(outfile, "\tMaxiter                = %4d\n", params.maxiter);
    fprintf(outfile, "\tConvergence            = %3.1e\n", params.convergence);
  }
  fprintf(outfile, "\tPrint Level            =    %1d\n",  params.print);
  fprintf(outfile, "\tLocal CC        =     %s\n", params.local ? "Yes" : "No");
  if(params.local) {
    fprintf(outfile, "\tLocal Cutoff    = %3.1e\n", local.cutoff);
    fprintf(outfile, "\tLocal Method    =    %s\n", local.method);
    fprintf(outfile, "\tWeak pairs      =    %s\n", local.weakp);
    fprintf(outfile, "\tGhost atom      =    %d\n", local.ghost);
    fprintf(outfile, "\tAmp Print Cutoff =   %f\n", local.amp_print_cutoff);
  }
  fprintf(outfile, "\n");

}


}} // namespace psi::cis
