/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include <physconst.h>
#include "Params.h"
#include "MOInfo.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

void get_params()
{
  int errcod, tol, count, i;
  char *junk, units[20];
  int *mu_irreps;
  
  errcod = ip_string("WFN", &(params.wfn), 0);
  if(strcmp(params.wfn, "MP2") && strcmp(params.wfn, "CCSD") && 
     strcmp(params.wfn, "CCSD_T") && strcmp(params.wfn, "EOM_CCSD") && 
     strcmp(params.wfn, "LEOM_CCSD") && strcmp(params.wfn, "BCCD") && 
     strcmp(params.wfn,"BCCD_T") && strcmp(params.wfn, "SCF") &&
     strcmp(params.wfn,"CIS") && strcmp(params.wfn,"RPA") &&
     strcmp(params.wfn,"CC2") && strcmp(params.wfn,"CC3") &&
     strcmp(params.wfn,"EOM_CC3") && strcmp(params.wfn,"EOM_CC2") &&
     strcmp(params.wfn,"CCSD_MVD")) {
    fprintf(outfile, "Invalid value of input keyword WFN: %s\n", params.wfn);
    exit(PSI_RETURN_FAILURE);
  }

  /* NB: SCF wfns are allowed because, at present, ccsort is needed for 
     RPA calculations */
  
  params.semicanonical = 0;
  errcod = ip_string("REFERENCE", &(junk),0);
  if (errcod != IPE_OK) {
    /* if no reference is given, assume rhf */
    params.ref = 0;
  }
  else {
    if(!strcmp(junk, "RHF")) params.ref = 0;
    else if(!strcmp(junk,"ROHF") && (!strcmp(params.wfn,"MP2") || !strcmp(params.wfn,"CCSD_T") || 
				     !strcmp(params.wfn,"CC3") || !strcmp(params.wfn, "EOM_CC3") ||
				     !strcmp(params.wfn,"CC2") || !strcmp(params.wfn, "EOM_CC2"))) {
      params.ref = 2;
      params.semicanonical = 1;
    }
    else if(!strcmp(junk, "ROHF")) params.ref = 1;
    else if(!strcmp(junk, "UHF")) params.ref = 2;
    else { 
      printf("Invalid value of input keyword REFERENCE: %s\n", junk);
      exit(PSI_RETURN_FAILURE); 
    }
    free(junk);
  }


  params.dertype = 0;
  if(ip_exist("DERTYPE",0)) {
    errcod = ip_string("DERTYPE", &(junk),0);
    if(errcod != IPE_OK) params.dertype = 0;
    else if(!strcmp(junk,"NONE")) params.dertype = 0;
    else if(!strcmp(junk,"FIRST")) params.dertype = 1;
    else if(!strcmp(junk,"RESPONSE")) params.dertype = 3; /* linear response */
    else {
      printf("Invalid value of input keyword DERTYPE: %s\n", junk);
      exit(PSI_RETURN_FAILURE); 
    }
    free(junk);
  }
	else {  /* if jobtype=opt and dertype is absent, infer dertype = 1 */
    if(ip_exist("JOBTYPE",0)) {
      errcod = ip_string("JOBTYPE", &(junk),0);
      if(!strcmp(junk,"OPT")) params.dertype = 1;
			free(junk);
		}
	}

  if(ip_exist("PROPERTY",0)) {
    errcod = ip_string("PROPERTY", &(params.prop), 0);
    if(strcmp(params.prop,"POLARIZABILITY") && strcmp(params.prop,"ROTATION") &&
       strcmp(params.prop,"ALL") && strcmp(params.prop,"MAGNETIZABILITY")
       && strcmp(params.prop,"ROA")) {
      fprintf(outfile, "Invalid choice of response property: %s\n", params.prop);
      exit(PSI_RETURN_FAILURE);
    }
  }
  else params.prop = strdup("POLARIZABILITY");
                                                                                                              
  params.local = 0;
  errcod = ip_boolean("LOCAL", &(params.local),0);
  local.cutoff = 0.02;
  errcod = ip_data("LOCAL_CUTOFF", "%lf", &(local.cutoff), 0);

  local.cphf_cutoff = 0.10;
  ip_data("LOCAL_CPHF_CUTOFF", "%lf", &(local.cphf_cutoff), 0);

  local.core_cutoff = 0.05;
  ip_data("LOCAL_CORE_CUTOFF", "%lf", &(local.core_cutoff), 0);

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
    sprintf(local.weakp, "%s", "NONE");
  }

  local.freeze_core = NULL;
  ip_string("FREEZE_CORE", &local.freeze_core, 0);
  if(local.freeze_core == NULL) local.freeze_core = strdup("FALSE");

  if(ip_exist("LOCAL_PAIRDEF",0)){
    errcod = ip_string("LOCAL_PAIRDEF", &(local.pairdef), 0);
    if(strcmp(local.pairdef,"BP") && strcmp(local.pairdef,"RESPONSE")) {
      fprintf(outfile, "Invalid keyword for strong/weak pair definition: %s\n", local.pairdef);
      exit(PSI_RETURN_FAILURE);
    }
  }
  else if(params.local && params.dertype == 3)
    local.pairdef = strdup("RESPONSE");
  else if(params.local)
    local.pairdef = strdup("BP");

  if(params.local && (strcmp(local.pairdef,"BP"))) {
    if(params.dertype == 3) {
      if(!(strcmp(params.prop,"POLARIZABILITY"))) {
	local.domain_polar = 1;
	local.domain_mag = 0;
      }
      else if(!(strcmp(params.prop,"MAGNETIZABILITY"))) {
	local.domain_polar = 0;
	local.domain_mag = 1;
      }
      else if(!(strcmp(params.prop,"ROTATION"))) {
	local.domain_polar = 1;
	local.domain_mag = 1;
      }
    }
    else {
      local.domain_polar = 0;
      local.domain_mag = 0;
    }
  }
  ip_boolean("LOCAL_DOMAIN_POLAR", &(local.domain_polar), 0);
  ip_boolean("LOCAL_DOMAIN_MAG", &(local.domain_mag), 0);

  local.domain_sep = 0;
  errcod = ip_boolean("LOCAL_DOMAIN_SEP", &(local.domain_sep),0);

  if(params.dertype == 3)
    local.filter_singles = 0;
  else
    local.filter_singles = 1;
  ip_boolean("LOCAL_FILTER_SINGLES", &(local.filter_singles), 0);

  if(ip_exist("AO_BASIS",0)) {
    errcod = ip_string("AO_BASIS", &(params.aobasis),0);
  }
  else params.aobasis = strdup("NONE");

  /* Do we need MO-basis <ab|cd> integrals? */
  if(!strcmp(params.wfn,"MP2") || !strcmp(params.aobasis,"DISK") ||
     !strcmp(params.aobasis,"DIRECT")) {
    params.make_abcd = 0;
  }
  else {
    params.make_abcd = 1;
  }

  /* If we need the MO-basis <ab|cd> integrals, do we need the fully unpacked list? */
  params.make_unpacked_abcd = 0;
  if(params.make_abcd) {
    if(params.ref != 0 || params.dertype == 1 || !strcmp(params.wfn,"EOM_CC2") ||
       !strcmp(params.wfn,"CC3") || !strcmp(params.wfn,"EOM_CC3") || !strcmp(params.wfn,"CCSD_MVD")) {
      params.make_unpacked_abcd = 1;
    }
    errcod = ip_string("EOM_REFERENCE", &(junk), 0);
    if (errcod == IPE_OK) {
      if(!strcmp(junk,"ROHF")) params.make_unpacked_abcd = 1;
      free(junk);
    }
  }

  /* for now, generate <ai|bc> ordering if CC gradient, ROHF-CC, CC2, or CC3 */
  if(params.dertype == 1 || params.ref == 1 || !strcmp(params.wfn,"CC2") ||
     !strcmp(params.wfn,"CC3") || !strcmp(params.wfn,"EOM_CC3") ||
     !strcmp(params.wfn,"EOM_CC2"))
    params.make_aibc = 1;
  else params.make_aibc = 0;
  errcod = ip_string("EOM_REFERENCE", &(junk), 0);
  if (errcod == IPE_OK) {
    if(!strcmp(junk,"ROHF")) params.make_aibc = 1;
    free(junk);
  }
  
  params.print_lvl = 1;
  errcod = ip_data("PRINT","%d",&(params.print_lvl),0);

  params.keep_TEIFile = 0;
  errcod = ip_boolean("KEEP_TEIFILE",&(params.keep_TEIFile),0);

  params.keep_OEIFile = 0;
  errcod = ip_boolean("KEEP_OEIFILE",&(params.keep_OEIFile),0);

  params.tolerance = 1e-14;
  errcod = ip_data("TOLERANCE", "%d", &(tol),0);
  if(errcod == IPE_OK) params.tolerance = 1.0*pow(10.0,(double) -tol);

  fndcor(&(params.memory), infile, outfile);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);

  params.local = 0;
  errcod = ip_boolean("LOCAL", &(params.local),0);

  /* grab the field frequencies from input -- a few different units are converted to E_h */
  if(ip_exist("OMEGA",0)) {
    errcod = ip_count("OMEGA", &count, 0);

    if(errcod == IPE_NOT_AN_ARRAY || count == 1) { /* assume Hartrees and only one frequency */
      params.nomega = 1;
      params.omega = init_array(1);
      errcod = ip_data("OMEGA", "%lf", &(params.omega[0]), 0);
    }
    else if(count >= 2) {
      params.nomega = count-1;
      params.omega = init_array(params.nomega);

      errcod = ip_data("OMEGA", "%s", units, 1, count-1);
      for(junk = units; *junk != '\0'; junk++)
	if(*junk>='a' && *junk <= 'z') *junk += 'A' - 'a';

      for(i=0; i < count-1; i++) {
	errcod = ip_data("OMEGA", "%lf", &(params.omega[i]), 1, i);

	if(!strcmp(units, "HZ")) params.omega[i] *= _h / _hartree2J;
	else if(!strcmp(units, "AU")) 1; /* do nothing */
	else if(!strcmp(units, "NM")) params.omega[i] = (_c*_h*1e9)/(params.omega[i]*_hartree2J);
	else if(!strcmp(units, "EV")) params.omega[i] /= _hartree2ev;
	else {
	  fprintf(outfile, "\n\tError in unit for input field frequencies.  Must use one of:\n");
	  fprintf(outfile,   "\tau, hz, nm, or ev.\n");
	  exit(PSI_RETURN_FAILURE);
	}
      }
    }
    else {
      fprintf(outfile, "\n\tError reading input field frequencies.  Please use the format:\n");
      fprintf(outfile,   "\t  omega = (value1 value2 ... units)\n");
      fprintf(outfile,   "\twhere units = hartrees, hz, nm, or ev.\n");
      exit(PSI_RETURN_FAILURE);
    }
  }
  else { /* assume static field by default */
    params.omega = init_array(1);
    params.omega[0] = 0.0;
    params.nomega = 1;
  }

  mu_irreps = init_int_array(3);
  errcod = ip_int_array("MU_IRREPS", mu_irreps, 3);
  if(errcod == IPE_OK) {
    moinfo.irrep_x = mu_irreps[0];
    moinfo.irrep_y = mu_irreps[1];
    moinfo.irrep_z = mu_irreps[2];
  }
  else {
    if(params.dertype == 3) {
      fprintf(outfile, "\nYou must supply the irreps of x, y, and z with the MU_IRREPS keyword.\n");
      exit(PSI_RETURN_FAILURE);
    }
  }
  free(mu_irreps);

  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tWave function   =\t%s\n", params.wfn);
  if(params.semicanonical) {
    fprintf(outfile, "\tReference wfn   =\tROHF changed to UHF for Semicanonical Orbitals\n"); 
  }
  else {
    fprintf(outfile, "\tReference wfn   =\t%s\n", 
	    (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
  }	  
  if(params.dertype == 0) fprintf(outfile, "\tDerivative      =\tNone\n");
  else if(params.dertype == 1) fprintf(outfile, "\tDerivative      =\tFirst\n");
  else if(params.dertype == 3) fprintf(outfile, "\tDerivative      =\tResponse\n");
  fprintf(outfile, "\tMemory (Mbytes) =\t%.1f\n", params.memory/1e6);
  fprintf(outfile, "\tAO Basis        =\t%s\n", params.aobasis);
  fprintf(outfile, "\tMake (ab|cd)    =\t%s\n", 
	  (params.make_abcd == 1) ? "True" : "False");
  fprintf(outfile, "\tMake unpacked (ab|cd) =\t%s\n", 
	  (params.make_unpacked_abcd == 1) ? "True" : "False");
  fprintf(outfile, "\tCache Level     =\t%d\n", params.cachelev);
  fprintf(outfile, "\tCache Type      =\t%s\n", "LRU");
  fprintf(outfile, "\tLocal CC        =     %s\n", params.local ? "Yes" : "No");
  if(params.local) {
    fprintf(outfile, "\tLocal Cutoff       = %3.1e\n", local.cutoff);
    fprintf(outfile, "\tLocal Core Cutoff  = %3.1e\n", local.core_cutoff);
    fprintf(outfile, "\tLocal Method      =    %s\n", local.method);
    fprintf(outfile, "\tWeak pairs        =    %s\n", local.weakp);
    fprintf(outfile, "\tFilter singles    =    %s\n", local.filter_singles ? "Yes" : "No");
    fprintf(outfile, "\tLocal pairs       =    %s\n", local.pairdef);
    fprintf(outfile, "\tLocal CPHF cutoff =  %3.1e\n", local.cphf_cutoff);
  }
  fprintf(outfile, "\n");
  fflush(outfile);
}

}} // namespace psi::ccsort
