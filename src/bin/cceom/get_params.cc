/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void get_params(void)
{
  int errcod, iconv;
  char *cachetype = NULL;
  char *read_ref, *read_eom_ref;

  errcod = ip_string("WFN", &(params.wfn), 0);

  if(!strcmp(params.wfn,"EOM_CC2")) {
    psio_read_entry(CC_INFO, "CC2 Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,"\tCC2 energy          (file100) = %20.15f\n",moinfo.ecc);
  }
  else if(!strcmp(params.wfn,"EOM_CCSD")) {
    psio_read_entry(CC_INFO, "CCSD Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,"\tCCSD energy         (file100) = %20.15f\n",moinfo.ecc);
  }
  else if(!strcmp(params.wfn,"EOM_CC3")) {
    psio_read_entry(CC_INFO, "CC3 Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,"\tCC3 energy          (file100) = %20.15f\n",moinfo.ecc);
  }

  params.semicanonical = 0;
  errcod = ip_string("REFERENCE", &(read_ref),0);
  if(!strcmp(read_ref, "RHF")) params.ref = 0;
  else if(!strcmp(read_ref,"ROHF") && (!strcmp(params.wfn,"EOM_CC3"))) {
    params.ref = 2;
    params.semicanonical = 1;
  }
  else if(!strcmp(read_ref, "ROHF")) params.ref = 1;
  else if(!strcmp(read_ref, "UHF")) params.ref = 2;
  else { 
    fprintf(outfile,
	    "\nInvalid value of input keyword REFERENCE: %s\n", read_ref);
    exit(2); 
  }

  if (params.ref == 0) { /* for RHF refs, allow CCEOM to do RHF, ROHF, UHF modes */
    errcod = ip_string("EOM_REFERENCE", &(read_eom_ref),0);
    if (errcod == IPE_OK) {
      if(!strcmp(read_eom_ref, "RHF")) params.eom_ref = 0;
      else if(!strcmp(read_eom_ref, "ROHF")) params.eom_ref = 1;
      else if(!strcmp(read_eom_ref, "UHF")) params.eom_ref = 2;
      else { 
        fprintf(outfile,
		"\nInvalid value of input keyword EOM_REFERENCE: %s\n", read_eom_ref);
        exit(2); 
      }
    }
    else {
      params.eom_ref = 0;
      read_eom_ref = (char *) malloc(10*sizeof(char));
      sprintf(read_eom_ref,"%s","RHF"); /* just for printing below */
    }
  }
  else if (params.ref == 1) { /* for ROHF refs, allow CCEOM to do ROHF & UHF modes */
    errcod = ip_string("EOM_REFERENCE", &(read_eom_ref),0);
    if (errcod == IPE_OK) {
      if(!strcmp(read_eom_ref, "ROHF")) params.eom_ref = 1;
      else if(!strcmp(read_eom_ref, "UHF")) params.eom_ref = 2;
      else { 
        fprintf(outfile,
		"\nInvalid value of input keyword EOM_REFERENCE: %s\n", read_eom_ref);
        exit(2); 
      }
    }
    else {
      params.eom_ref = 1;
      read_eom_ref = (char *) malloc(10*sizeof(char));
      sprintf(read_eom_ref,"%s","ROHF"); /* just for printing below */
    }
  }
  else { /* run in UHF mode - ignore EOM_REFERENCE */
    params.eom_ref = 2;
    read_eom_ref = (char *) malloc(10*sizeof(char));
    sprintf(read_eom_ref,"%s","UHF"); /* just for printing below */
  }

  /*  fprintf(outfile, "\nCCEOM not yet UHF capable\n"); */

  fndcor(&(params.memory),infile,outfile);

  params.aobasis = 0;
  errcod = ip_boolean("AO_BASIS", &(params.aobasis),0);

  params.full_matrix = 0;
  errcod = ip_boolean("FULL_MATRIX", &(params.full_matrix),0);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);

  params.cachetype = 0;
  errcod = ip_string("CACHETYPE", &(cachetype),0);
  if(cachetype != NULL && strlen(cachetype)) {
    if(!strcmp(cachetype,"LOW")) params.cachetype = 1;
    else if(!strcmp(cachetype,"LRU")) params.cachetype = 0;
    else {
      fprintf(outfile, "Error in input: invalid CACHETYPE = %s\n",
	      cachetype);
      exit(1);
    }
    free(cachetype);
  }
  if(params.ref == 2) /* No LRU cacheing yet for UHF references */
    params.cachetype = 0;

  params.nthreads = 1;
  errcod = ip_data("NTHREADS", "%d", &(params.nthreads),0);

  if(ip_exist("ABCD",0)) {
    errcod = ip_string("ABCD", &(params.abcd), 0);
    if(strcmp(params.abcd,"NEW") && strcmp(params.abcd,"OLD")) {
      fprintf(outfile, "Invalid ABCD algorithm: %s\n", params.abcd);
      exit(PSI_RETURN_FAILURE);
    }
  }
  else params.abcd = strdup("NEW");

  params.t3_Ws_incore = 0;
  errcod = ip_boolean("T3_WS_INCORE", &(params.t3_Ws_incore),0);

  params.local = 0;
  errcod = ip_boolean("LOCAL", &(params.local),0);

  local.cutoff = 0.02;
  errcod = ip_data("LOCAL_CUTOFF", "%lf", &(local.cutoff), 0);

  if(ip_exist("LOCAL_METHOD",0)) {
    errcod = ip_string("LOCAL_METHOD", &(local.method), 0);
    if(strcmp(local.method,"AOBASIS") && strcmp(local.method,"WERNER")) {
      fprintf(outfile, "\nInvalid local correlation method: %s\n", local.method);
      exit(2);
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
      exit(2);
    }
  }
  else if(params.local) {
    local.weakp = (char *) malloc(5 * sizeof(char));
    sprintf(local.weakp, "%s", "NONE");
  }

  if(ip_exist("LOCAL_PRECONDITIONER",0) && params.local) {
    errcod = ip_string("LOCAL_PRECONDITIONER", &(local.precon), 0);
    if(strcmp(local.precon,"FOCK") && strcmp(local.precon,"HBAR")) {
      fprintf(outfile, "Invalid choice of local-pair preconditioner: %s\n", local.precon);
      exit(2);
    }
  }
  else if(params.local){
    local.precon = (char *) malloc(5 * sizeof(char));
    sprintf(local.precon, "%s", "HBAR");
  }

  local.ghost = -1;
  if(ip_exist("LOCAL_GHOST",0))
    errcod = ip_data("LOCAL_GHOST", "%d", &(local.ghost), 0);

  local.do_singles = 1;
  errcod = ip_boolean("LOCAL_DO_SINGLES", &(local.do_singles), 0);

  local.filter_singles = 1;
  ip_boolean("LOCAL_FILTER_SINGLES", &(local.filter_singles), 0);

  params.newtrips = 1;
  ip_boolean("NEWTRIPS", &params.newtrips, 0);

  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  if(params.semicanonical)
    fprintf(outfile, "\tReference wfn   = ROHF changed to UHF for Semicanonical Orbitals\n");
  else 
    fprintf(outfile, "\tReference wfn   =    %4s\n", read_ref);
  fprintf(outfile, "\tReference EOM wfn=    %4s\n", read_eom_ref);
  fprintf(outfile, "\tMemory (Mbytes) =  %5.1f\n",params.memory/1e6);
  fprintf(outfile, "\tAO Basis        =     %s\n", 
	  params.aobasis ? "Yes" : "No");
  fprintf(outfile, "\tABCD            =     %s\n", params.abcd);
  fprintf(outfile, "\tCache Level     =    %1d\n", 
	  params.cachelev);
  fprintf(outfile, "\tCache Type      =    %4s\n", 
	  params.cachetype ? "LOW" : "LRU");
  if ( !strcmp(params.wfn,"EOM_CC3") )
    fprintf(outfile, "\tT3 Ws incore  =    %4s\n", params.t3_Ws_incore ? "Yes" : "No");
  fprintf(outfile, "\tNum. of threads =     %d\n",params.nthreads);
  fprintf(outfile, "\tLocal CC        =     %s\n", params.local ? "Yes" : "No");
  if(params.local) {
    fprintf(outfile, "\tLocal Cutoff    = %3.1e\n", local.cutoff);
    fprintf(outfile, "\tLocal Method    =    %s\n", local.method);
    fprintf(outfile, "\tWeak pairs      =    %s\n", local.weakp);
    fprintf(outfile, "\tLocal precon.   =    %s\n", local.precon);
    fprintf(outfile, "\tGhost atom      =    %d\n", local.ghost);
    fprintf(outfile, "\tLocal guess     =    %s\n", 
	    local.do_singles ? "HBAR_SS" : "UNIT VECTORS" );
    fprintf(outfile, "\tFilter singles  =    %s\n", local.filter_singles ? "Yes" : "No");
  }
  fprintf(outfile, "\n");

  free(read_ref);
  free(read_eom_ref);
}


}} // namespace psi::cceom
