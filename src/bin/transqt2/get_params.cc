/*! \file
    \ingroup TRANSQT2
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi {
  namespace transqt2 {

void get_params()
{
  int errcod, tol;
  char *junk;
  int *mu_irreps, tmp;

  errcod = ip_string("WFN", &(params.wfn), 0);
  /*
  if(strcmp(params.wfn, "MP2") && strcmp(params.wfn, "CCSD") &&
     strcmp(params.wfn, "CCSD_T") && strcmp(params.wfn, "EOM_CCSD") &&
     strcmp(params.wfn, "LEOM_CCSD") && strcmp(params.wfn, "BCCD") &&
     strcmp(params.wfn,"BCCD_T") && strcmp(params.wfn, "SCF") &&
     strcmp(params.wfn,"CIS") && strcmp(params.wfn,"RPA") &&
     strcmp(params.wfn,"CC2") && strcmp(params.wfn,"CC3") &&
     strcmp(params.wfn,"EOM_CC3") && strcmp(params.wfn,"EOM_CC2") &&
     strcmp(params.wfn,"CCSD_MVD") ) {
    printf("Invalid value of input keyword WFN: %s\n", params.wfn);
    exit(PSI_RETURN_FAILURE);
  }
  */

  /* NB: SCF wfns are allowed because, at present, ccsort is needed for
     RPA-type calculations */

  params.semicanonical = 0;
  errcod = ip_string("REFERENCE", &(junk),0);
  if (errcod != IPE_OK)
    params.ref = 0; /* if no reference is given, assume rhf */
  else {
    if(!strcmp(junk, "RHF")) params.ref = 0;
    else if(!strcmp(junk,"ROHF") &&
	    (!strcmp(params.wfn,"MP2") || !strcmp(params.wfn,"CCSD_T") ||
	     !strcmp(params.wfn,"CC3") || !strcmp(params.wfn, "EOM_CC3") ||
	     !strcmp(params.wfn,"CC2") || !strcmp(params.wfn, "EOM_CC2"))) {
      params.ref = 2;
      params.semicanonical = 1;
    }
    else if(!strcmp(junk, "ROHF")) params.ref = 1;
    else if(!strcmp(junk, "UHF")) params.ref = 2;
    else if(!strcmp(junk, "TWOCON")) params.ref = 1; /* Treat twocon as rhf */
    else {
      printf("Invalid value of input keyword REFERENCE: %s\n", junk);
      exit(PSI_RETURN_FAILURE);
    }
    free(junk);
  }

  params.dertype = 0;
  if(ip_exist("DERTYPE",0)) {
    errcod = ip_string("DERTYPE", &(junk),0);
    if(errcod != IPE_OK) params.dertype = 0; /* assume no derivative */
    else if(!strcmp(junk,"NONE")) params.dertype = 0;
    else if(!strcmp(junk,"FIRST")) params.dertype = 1;
    else if(!strcmp(junk,"SECOND")) params.dertype = 2;
    else if(!strcmp(junk,"RESPONSE")) params.dertype = 3; /* linear response */
    else {
      printf("Invalid value of input keyword DERTYPE: %s\n", junk);
      exit(PSI_RETURN_FAILURE);
    }
    free(junk);
  }
  else {  /* if jobtype=opt and dertype is absent, assume dertype = 1 */
    if(ip_exist("JOBTYPE",0)) {
      errcod = ip_string("JOBTYPE", &(junk),0);
      if(!strcmp(junk,"OPT")) params.dertype = 1;
      free(junk);
    }
  }

  /*   params.print_lvl = 1; */ /* default set in init_io() */
  errcod = ip_data("PRINT","%d",&(params.print_lvl),0);

  params.print_tei = 0;
  errcod = ip_boolean("PRINT_TEI", &params.print_tei, 0);

  params.tolerance = 1e-14;
  errcod = ip_data("TOLERANCE", "%d", &(tol),0);
  if(errcod == IPE_OK) params.tolerance = 1.0*pow(10.0,(double) -tol);

  fndcor(&(params.memory), infile, outfile);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);

  params.delete_tei = 1;
  /* If AO-basis chosen, keep the SO_TEI file */
  {
    char* aobasis;
    if(ip_exist("AO_BASIS",0)) {
      errcod = ip_string("AO_BASIS", &aobasis,0);
    }
    else aobasis = strdup("NONE");
    if(!strcmp(aobasis,"DISK")) {
      params.delete_tei = 0;
    }
  }
  // any MCSCF-type wavefunction needs multiple transforms so don't delete
  // the AO two-electron ints
  if ((strcmp(params.wfn,"OOCCD")==0 || strcmp(params.wfn,"DETCAS")==0 ||
       strcmp(params.wfn,"CASSCF")==0|| strcmp(params.wfn,"RASSCF")==0))
    params.delete_tei = 0;

  errcod = ip_boolean("DELETE_TEI", &params.delete_tei, 0);

  if(params.print_lvl) {
    fprintf(outfile, "\n\tInput parameters:\n");
    fprintf(outfile, "\t-----------------\n");
    fprintf(outfile, "\tWave function   =\t%s\n", params.wfn);
    fprintf(outfile, "\tBacktransform   =\t%s\n", params.backtr ? "Yes" : "No");
    fprintf(outfile, "\tPrint Level     =\t%d\n", params.print_lvl);
    fprintf(outfile, "\tPrint TEIs      =\t%s\n", params.print_tei ? "Yes" : "No");
    if(params.semicanonical) {
      fprintf(outfile, "\tReference wfn   =\tROHF (using UHF for semicanonical orbitals)\n");
    }
    else {
      fprintf(outfile, "\tReference wfn   =\t%s\n",
	      (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
    }
    if(params.dertype == 0) fprintf(outfile, "\tDerivative      =\tNone\n");
    else if(params.dertype == 1) fprintf(outfile, "\tDerivative      =\tFirst\n");
    else if(params.dertype == 2) fprintf(outfile, "\tDerivative      =\tSecond\n");
    else if(params.dertype == 3) fprintf(outfile, "\tDerivative      =\tResponse\n");
    fprintf(outfile, "\tDelete TEI File =\t%s\n", params.delete_tei ? "Yes" : "No");
    fprintf(outfile, "\tMemory (Mbytes) =\t%.1f\n", params.memory/1e6);
    fprintf(outfile, "\tCache Level     =\t%d\n", params.cachelev);
    fprintf(outfile, "\tCache Type      =\t%s\n", "LRU");
    fflush(outfile);
  }
}

  } // namespace transqt2
} // namespace psi
