/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void get_params()
{
  int errcod;
  char *cachetype = NULL;
  char *junk;
  
  errcod = ip_string("WFN", &(params.wfn), 0);

  errcod = ip_string("REFERENCE", &(junk),0);

  /* Default reference is RHF */
  params.ref = 0;
  params.semicanonical = 0;
  if(!strcmp(junk,"RHF")) params.ref = 0;
  else if(!strcmp(junk,"ROHF") && !strcmp(params.wfn,"MP2")) {
    params.ref = 2;
    params.semicanonical = 1;
  }
  else if(!strcmp(junk,"ROHF")) params.ref = 1;
  else if(!strcmp(junk,"UHF")) params.ref = 2;
  else {
    fprintf(outfile,"\nInvalid Reference: %s\n",junk);
    exit(PSI_RETURN_FAILURE);
  }
  free(junk);
  
  /* Default Jobtype */
  if(ip_exist("JOBTYPE",0)) {
    ip_string("JOBTYPE", &(params.jobtype),0);
  }
  else {
    params.jobtype = strdup("SP");
  }
    
  /* Default Dertype */
  if(ip_exist("DERTYPE",0)) {
    ip_string("DERTYPE", &(params.dertype),0);
  }
  else {
    params.dertype = strdup("NONE");
  }

  if(!strcmp(params.jobtype,"SP")) {
    params.opdm = 0;
    params.relax_opdm = 0;
    params.gradient = 0;
  }
  else if(!strcmp(params.jobtype,"OEPROP") && !strcmp(params.dertype,"NONE")) {
    params.opdm = 1;
    params.relax_opdm = 0;
    params.gradient = 0;
  }
  else if(!strcmp(params.jobtype,"OEPROP") && !strcmp(params.dertype,"FIRST")) {
    params.opdm = 1;
    params.relax_opdm = 1;
    params.gradient = 0;
  }
  else if(!strcmp(params.jobtype,"OPT") && !strcmp(params.dertype,"NONE")) {
    params.opdm = 0;
    params.relax_opdm = 0;
    params.gradient = 0;
  }
  else if(!strcmp(params.jobtype,"OPT") && !strcmp(params.dertype,"FIRST")) {
    params.opdm = 0;
    params.relax_opdm = 0;
    params.gradient = 1;
  }
  else if(!strcmp(params.jobtype,"OPT_FC") && !strcmp(params.dertype,"FIRST")) {
    params.opdm = 0;
    params.relax_opdm = 0;
    params.gradient = 1;
  }
  else if(!strcmp(params.jobtype,"SYMM_FC") && !strcmp(params.dertype,"FIRST")) {
    params.opdm = 0;
    params.relax_opdm = 0;
    params.gradient = 1;
  }
  else if(!strcmp(params.jobtype,"FREQ") && !strcmp(params.dertype,"NONE")) {
    params.opdm = 0;
    params.relax_opdm = 0;
    params.gradient = 0;
  }
  else if(!strcmp(params.jobtype,"FREQ") && !strcmp(params.dertype,"FIRST")) {
    params.opdm = 0;
    params.relax_opdm = 0;
    params.gradient = 1;
  }
  else {
    printf("Invalid combination of JOBTYPE and DERTYPE\n");
    exit(PSI_RETURN_FAILURE);
  }

  if((params.relax_opdm || params.gradient) && 
     (mo.nfzdocc != 0 || mo.nfzvirt != 0)) {
    fprintf(outfile,"\n\tThe Z-vector equations DO NOT work with frozen orbitals ... yet\n");
    exit(PSI_RETURN_FAILURE);
  }

  params.print = 0;
  ip_data("PRINT", "%d", &(params.print),0);
  
  params.cachelev = 2;
  ip_data("CACHELEV", "%d", &(params.cachelev),0);
  
  params.cachetype = 1;
  errcod = ip_string("CACHETYPE", &(cachetype),0);
  if (cachetype != NULL && strlen(cachetype)) {
    if (!strcmp(cachetype,"LOW")) 
      params.cachetype = 1;
    else if (!strcmp(cachetype,"LRU")) 
      params.cachetype = 0;
    else {
      fprintf(outfile, "Invalide CACHETYPE = %s\n",cachetype);
      abort();
    }
    free(cachetype);
  }
  
  /* get parameters related to SCS-MP2 or SCS-N-MP2 */
  /* see papers by S. Grimme or J. Platz */
  params.scs = 0;
  params.scs_scale_os = 6.0/5.0;
  params.scs_scale_ss = 1.0/3.0;
  errcod = ip_boolean("SCS_N",&(params.scs),0);
  if (params.scs == 1) {
    params.scs_scale_os = 0.0;
    params.scs_scale_ss = 1.76;
  }
  errcod = ip_boolean("SCS",&(params.scs),0);
  if (params.scs == 1) { 
    errcod = ip_data("SCALE_OS","%lf",&(params.scs_scale_os),0); 
    errcod = ip_data("SCALE_SS","%lf",&(params.scs_scale_ss),0); 
    }

  fndcor(&(params.memory),infile,outfile);
 
  fprintf(outfile, "\n");
  fprintf(outfile, "\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tWave function \t=\t%s\n", params.wfn);
  if(params.semicanonical) {
  fprintf(outfile, "\tReference WFN \t=\tROHF changed to UHF for Semicanonical Orbitals\n");
  }
  else {
  fprintf(outfile, "\tReference WFN \t=\t%s\n", (params.ref==0)?"RHF":((params.ref==1)?"ROHF":"UHF"));
  } 
  fprintf(outfile, "\tDerivative    \t=\t%s\n", params.dertype);
  fprintf(outfile, "\tCache Level   \t=\t%d\n", params.cachelev);
  fprintf(outfile, "\tCache Type    \t=\t%s\n", params.cachetype ? "LOW":"LRU");
  fprintf(outfile, "\tMemory (MB)   \t=\t%.1f\n",params.memory/1e6);
  fprintf(outfile, "\tPrint Level   \t=\t%d\n", params.print);
  fprintf(outfile, "\tOPDM          \t=\t%s\n", params.opdm ? "YES":"NO");
  fprintf(outfile, "\tSCS           \t=\t%s\n", (params.scs == 1) ? "True" : "False");
  fprintf(outfile, "\tSCALE_OS      \t=\t%.6f\n",params.scs_scale_os);
  fprintf(outfile, "\tSCALE_SS      \t=\t%.6f\n",params.scs_scale_ss);

  if (params.scs && (strcmp(params.dertype,"NONE")!=0)) {
    fprintf(outfile,"\nWarning: SCS-MP2 computation requested but\n");
    fprintf(outfile,"derivative will be evaluated for standard MP2 energy.\n");
  }

}

}} /* End namespaces */
