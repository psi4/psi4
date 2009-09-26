/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
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
  
  params.wfn = options.get_cstr("WFN");
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
  junk = options.get_cstr("REFERENCE");
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

  junk = options.get_cstr("DERTYPE");
  else if(!strcmp(junk,"NONE")) params.dertype = 0;
  else if(!strcmp(junk,"FIRST")) params.dertype = 1;
  else if(!strcmp(junk,"RESPONSE")) params.dertype = 3; /* linear response */
  else {
    printf("Invalid value of input keyword DERTYPE: %s\n", junk);
    exit(PSI_RETURN_FAILURE); 
  }

  params.prop = options.get_cstr("PROPERTY");
  if(strcmp(params.prop,"POLARIZABILITY") && strcmp(params.prop,"ROTATION") &&
     strcmp(params.prop,"ALL") && strcmp(params.prop,"MAGNETIZABILITY")
     && strcmp(params.prop,"ROA")) {
     fprintf(outfile, "Invalid choice of response property: %s\n", params.prop);
     throw PsiException("ccsort error", __FILE__, __LINE__);
  }

  params.local = options.get_bool("LOCAL");
  local.cutoff = options.get_double("LOCAL_CUTOFF");
  local.cphf_cutoff = options.get_double("LOCAL_CPHF_CUTOFF");
  local.core_cutoff = options.get_double("LOCAL_CORE_CUTOFF");

  sprintf(local.method, "%s", "WERNER");
  local.method = options.get_cstr("LOCAL_METHOD");
  if(strcmp(local.method,"AOBASIS") && strcmp(local.method,"WERNER")) {
    fprintf(outfile, "Invalid local correlation method: %s\n", local.method);
    throw PsiException("ccsort error", __FILE__, __LINE__);
  }

  sprintf(local.weakp, "%s", "NONE");
  local.weakp = options.get_cstr("LOCAL_WEAKP");
  if(strcmp(local.weakp,"MP2") && strcmp(local.weakp,"NEGLECT") && strcmp(local.weakp,"NONE")) {
    fprintf(outfile, "Invalid method for treating local pairs: %s\n", local.weakp);
    exit(PSI_RETURN_FAILURE);
  }

  local.freeze_core = NULL;
  local.freeze_core = options.get_cstr("FREEZE_CORE");

  local.pairdef = options.get_cstr("LOCAL_PAIRDEF");
  if(strcmp(local.pairdef,"BP") && strcmp(local.pairdef,"RESPONSE")) {
    fprintf(outfile, "Invalid keyword for strong/weak pair definition: %s\n", local.pairdef);
    exit(PSI_RETURN_FAILURE);
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
  local.domain_polar = options.get_bool("LOCAL_DOMAIN_POLAR");
  local.domain_mag = options.get_bool("LOCAL_DOMAIN_MAG");

  local.domain_sep = 0;
  local.domain_sep = options.get_bool("LOCAL_DOMAIN_SEP");

  if(params.dertype == 3)
    local.filter_singles = 0;
  else
    local.filter_singles = 1;
  local.filter_singles = options.get_bool("LOCAL_FILTER_SINGLES");

  params.aobasis = strdup("NONE");
  params.aobasis = options.get_cstr("AO_BASIS");

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
   junk = options.get_cstr("EOM_REFERENCE");
   if(!strcmp(junk,"ROHF")) params.make_unpacked_abcd = 1;
  }

  /* for now, generate <ai|bc> ordering if CC gradient, ROHF-CC, CC2, or CC3 */
  if(params.dertype == 1 || params.ref == 1 || !strcmp(params.wfn,"CC2") ||
     !strcmp(params.wfn,"CC3") || !strcmp(params.wfn,"EOM_CC3") ||
     !strcmp(params.wfn,"EOM_CC2"))
    params.make_aibc = 1;
  else params.make_aibc = 0;
  junk = options.get_cstr("EOM_REFERENCE");
  if(!strcmp(junk,"ROHF")) params.make_aibc = 1;
  
  params.print_lvl = 1;
  params.print_lvl = options.get_int("PRINT");

  params.keep_TEIFile = 0;
  params.ket_TEIFile = options.get_bool("KEEP_TEIFILE");

  params.keep_OEIFile = 0;
  params.keep_OEIFile = options.get_bool("KEEP_OEIFILE");

  params.tolerance = 1e-14;
  tol = options.get_int("TOLERANCE");
  params.tolerance = 1.0*pow(10.0,(double) -tol);

  params.memory = module.get_memory();

  params.cachelev = 2;
  params.cachelev = options.get_int("CACHELEV");

  params.local = 0;
  params.local = options.get_bool("LOCAL");

  /* grab the field frequencies from input -- a few different units are converted to E_h */
  count = options["OMEGA"].size();
  if(count == 1) { /* assume Hartrees and only one frequency */
    params.nomega = 1;
    params.omega = init_array(1);
    params.omega[0] = options["OMEGA"][0].to_double();
  }
  else if(count >= 2) {
    params.nomega = count-1;
    params.omega = init_array(params.nomega);

    units = options["OMEGA"][count-1].to_cstr();
    for(junk = units; *junk != '\0'; junk++)
      if(*junk>='a' && *junk <= 'z') *junk += 'A' - 'a';

    for(i=0; i < count-1; i++) {
      params.omega[i] = options["OMEGA"][i].to_double();

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
    throw PsiException("Failure in ccsort.", __FILE__, __LINE__);
  }

  else { /* assume static field by default */
    params.omega = init_array(1);
    params.omega[0] = 0.0;
    params.nomega = 1;
  }

  mu_irreps = init_int_array(3);
  moinfo.irrep_x = options["MU_IRREPS"][0].to_integer();
  moinfo.irrep_y = options["MU_IRREPS"][1].to_integer();
  moinfo.irrep_z = options["MU_IRREPS"][2].to_integer();

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
