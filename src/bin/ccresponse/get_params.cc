/*! \file
    \ingroup CCRESPONSE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include <physconst.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

void get_params()
{
  int i, errcod, ref, count, iconv, *tmpi;
  char *junk, units[20];

  errcod = ip_string("WFN", &(params.wfn), 0);
  if(strcmp(params.wfn, "CCSD") && strcmp(params.wfn, "CC2")) {
    fprintf(outfile, "Invalid value of input keyword WFN: %s\n", params.wfn);
    exit(PSI_RETURN_FAILURE);
  }

  params.print = 1;
  errcod = ip_data("PRINT","%d",&(params.print),0);

  fndcor(&(params.memory), infile, outfile);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);
  params.cachelev = 0;

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
    fprintf(outfile, "Value of REFERENCE from input.dat (%1d) and CC_INFO (%1d) do not match!\n", 
	    ref, params.ref);
    fprintf(outfile, "Is this what you want to do?\n");
    params.ref = ref;
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

  if(ip_exist("GAUGE",0)) {
    errcod = ip_string("GAUGE", &(params.gauge), 0);
    if(strcmp(params.gauge,"LENGTH") && strcmp(params.gauge,"VELOCITY") &&
       strcmp(params.gauge,"BOTH")) {
      printf("Invalid choice of gauge: %s\n", params.gauge);
      exit(PSI_RETURN_FAILURE);
    }
  }
  else
    params.gauge = strdup("LENGTH");

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

  moinfo.mu_irreps = init_int_array(3);
  errcod = ip_int_array("MU_IRREPS", moinfo.mu_irreps, 3);
  if(errcod != IPE_OK) {
    fprintf(outfile, "\nYou must supply the irreps of x, y, and z with the MU_IRREPS keyword.\n");
    exit(PSI_RETURN_FAILURE);
  }

  /* compute the irreps of the angular momentum operator while we're here */
  moinfo.l_irreps = init_int_array(3);
  for(i=0; i < 3; i++)
    moinfo.l_irreps[i] = moinfo.mu_irreps[(int) (i+1)%3] ^ moinfo.mu_irreps[(int) (i+2)%3];

  params.maxiter = 50;
  errcod = ip_data("MAXITER","%d",&(params.maxiter),0);
  params.convergence = 1e-7;
  errcod = ip_data("CONVERGENCE","%d",&(iconv),0);
  if(errcod == IPE_OK) params.convergence = 1.0*pow(10.0,(double) -iconv);
  params.diis = 1;
  errcod = ip_boolean("DIIS", &(params.diis), 0);

  if(ip_exist("PROPERTY",0)) {
    errcod = ip_string("PROPERTY", &(params.prop), 0);
    if(strcmp(params.prop,"POLARIZABILITY") && strcmp(params.prop,"ROTATION") 
       && strcmp(params.prop,"ROA") && strcmp(params.prop,"ALL")) {
      fprintf(outfile, "Invalid choice of resp. property: %s\n", params.prop);
      exit(PSI_RETURN_FAILURE);
    }
  }
  else params.prop = strdup("POLARIZABILITY");

  if(ip_exist("ABCD",0)) {
    errcod = ip_string("ABCD", &(params.abcd), 0);
    if(strcmp(params.abcd,"NEW") && strcmp(params.abcd,"OLD")) {
      fprintf(outfile, "Invalid ABCD algorithm: %s\n", params.abcd);
      exit(PSI_RETURN_FAILURE);
    }
  }
  else params.abcd = strdup("NEW");

  params.restart = 1;
  errcod = ip_boolean("RESTART", &params.restart, 0);

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
    sprintf(local.weakp, "%s", "NONE");
  }

  if(params.dertype == 3)
    local.filter_singles = 0;
  else
    local.filter_singles = 1;
  ip_boolean("LOCAL_FILTER_SINGLES", &(local.filter_singles), 0);

  local.cphf_cutoff = 0.10;
  ip_data("LOCAL_CPHF_CUTOFF", "%lf", &(local.cphf_cutoff), 0);

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

  params.analyze = 0;
  ip_boolean("ANALYZE", &(params.analyze), 0);

  params.num_amps = 5;
  if(ip_exist("NUM_AMPS", 0)) {
    errcod = ip_data("NUM_AMPS", "%d", &(params.num_amps),0);
  }

  params.sekino = 0;
  if(ip_exist("SEKINO",0)) {
    errcod = ip_boolean("SEKINO", &params.sekino, 0);
    if(errcod != IPE_OK) params.sekino = 0;
  }

  params.linear = 0;
  if(ip_exist("LINEAR",0)) {
    errcod = ip_boolean("LINEAR", &params.linear, 0);
    if(errcod != IPE_OK) params.linear = 0;
  }

  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  if(!strcmp(params.prop,"ALL"))
    fprintf(outfile, "\tProperty         =    POLARIZABILITY + ROTATION\n");
  else
    fprintf(outfile, "\tProperty         =    %s\n", params.prop);
  fprintf(outfile, "\tReference wfn    =    %5s\n",
	  (params.ref == 0) ? "RHF" : ((params.ref == 1) ? "ROHF" : "UHF"));
  fprintf(outfile, "\tMemory (Mbytes)  =  %5.1f\n",params.memory/1e6);
  fprintf(outfile, "\tCache Level      =    %1d\n", params.cachelev);
  fprintf(outfile, "\tPrint Level      =    %1d\n",  params.print);
  fprintf(outfile, "\tMaxiter          =    %3d\n",  params.maxiter);
  fprintf(outfile, "\tConvergence      = %3.1e\n", params.convergence);
  fprintf(outfile, "\tRestart          =     %s\n", params.restart ? "Allowed" : "Not Allowed");
  fprintf(outfile, "\tDIIS             =     %s\n", params.diis ? "Yes" : "No");
  fprintf(outfile, "\tModel III        =     %s\n", params.sekino ? "Yes" : "No");
  fprintf(outfile, "\tLinear Model     =     %s\n", params.linear ? "Yes" : "No");
  fprintf(outfile, "\tABCD             =     %s\n", params.abcd);
  fprintf(outfile, "\tIrrep X          =    %3s\n", moinfo.labels[moinfo.mu_irreps[0]]);
  fprintf(outfile, "\tIrrep Y          =    %3s\n", moinfo.labels[moinfo.mu_irreps[1]]);
  fprintf(outfile, "\tIrrep Z          =    %3s\n", moinfo.labels[moinfo.mu_irreps[2]]);
  fprintf(outfile, "\tIrrep RX         =    %3s\n", moinfo.labels[moinfo.l_irreps[0]]);
  fprintf(outfile, "\tIrrep RY         =    %3s\n", moinfo.labels[moinfo.l_irreps[1]]);
  fprintf(outfile, "\tIrrep RZ         =    %3s\n", moinfo.labels[moinfo.l_irreps[2]]);
  fprintf(outfile, "\tGauge            =    %s\n", params.gauge);
  for(i=0; i < params.nomega; i++) {
    if(params.omega[i] == 0.0) 
      fprintf(outfile, "\tApplied field %2d =  0.000\n", i);
    else 
      fprintf(outfile, "\tApplied field %2d =    %5.3f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", i, params.omega[i],
	      (_c*_h*1e9)/(_hartree2J*params.omega[i]), _hartree2ev*params.omega[i],
	      _hartree2wavenumbers*params.omega[i]);
  }
  fprintf(outfile, "\tAnalyze X2 Amps  =    %s\n", params.analyze ? "Yes" : "No");
  fprintf(outfile, "\tLocal CC         =    %s\n", params.local ? "Yes" : "No");
  if(params.local) {
    fprintf(outfile, "\tLocal Cutoff      = %3.1e\n", local.cutoff);
    fprintf(outfile, "\tLocal Method      =    %s\n", local.method);
    fprintf(outfile, "\tWeak pairs        =    %s\n", local.weakp);
    fprintf(outfile, "\tFilter singles    =    %s\n", local.filter_singles ? "Yes" : "No");
    fprintf(outfile, "\tLocal pairs       =    %s\n", local.pairdef);
    fprintf(outfile, "\tLocal CPHF cutoff =  %3.1e\n", local.cphf_cutoff);
  }
  fprintf(outfile, "\n");
}


}} // namespace psi::ccresponse
