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
#include <libchkpt/chkpt.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void get_eom_params(Options &options)
{
  int errcod, i, j, sym, iconv,exist, state_irrep, c_irrep;

  chkpt_init(PSIO_OPEN_OLD);
  if (chkpt_rd_override_occ()) {
    eom_params.states_per_irrep = chkpt_rd_statespi();
  }
  else if (ip_exist("STATES_PER_IRREP",0)) {
    eom_params.states_per_irrep = (int *) malloc(moinfo.nirreps * sizeof(int));
    ip_count("STATES_PER_IRREP", &i, 0);
    if (i != moinfo.nirreps) {
      fprintf(outfile,"Dim. of states_per_irrep vector must be %d\n", moinfo.nirreps) ;
      exit(0);
    }
    for (i=0;i<moinfo.nirreps;++i)
      errcod = ip_data("STATES_PER_IRREP","%d",&(eom_params.states_per_irrep[i]),1,i);
  }
  else { fprintf(outfile,"Must have states_per_irrep vector in input.\n"); exit(PSI_RETURN_FAILURE); } 
  chkpt_close();


  eom_params.cs_per_irrep = (int *) malloc(moinfo.nirreps * sizeof(int));
  eom_params.number_of_states = 0;
  for (state_irrep=0; state_irrep<moinfo.nirreps; ++state_irrep) {
    eom_params.cs_per_irrep[state_irrep^moinfo.sym] = eom_params.states_per_irrep[state_irrep];
	eom_params.number_of_states += eom_params.states_per_irrep[state_irrep];
  }
  eom_params.state_energies = (double*) malloc(eom_params.number_of_states * sizeof(double));

  eom_params.max_iter = 80 * moinfo.nirreps;
  errcod = ip_data("MAX_ITER","%d",&(eom_params.max_iter),0);

  /* Use prop_sym and prop_root only to determine what energy to write the file32 */
  if (ip_exist("PROP_SYM",0)) {
    ip_data("PROP_SYM","%d",&(eom_params.prop_sym),0);
    eom_params.prop_sym = (eom_params.prop_sym - 1)^moinfo.sym;
  }
  else {
    for (i=0;i<moinfo.nirreps;++i)
      if (eom_params.states_per_irrep[i])
        eom_params.prop_sym = i^moinfo.sym;
  }
  if (ip_exist("PROP_ROOT",0)) {
    ip_data("PROP_ROOT","%d",&(eom_params.prop_root),0);
    if (eom_params.prop_root > eom_params.states_per_irrep[eom_params.prop_sym^moinfo.sym]) {
      fprintf(outfile,"prop_root is too big\n");
      exit(1);
    }
  }
  else {
    eom_params.prop_root = eom_params.states_per_irrep[eom_params.prop_sym^moinfo.sym];
  }
  --eom_params.prop_root;

  // by default for CC3 use follow_root to root-following based on CCSD soln
  if ( (!strcmp(params.wfn,"EOM_CC3")) && (eom_params.prop_root != 0) ) {
    eom_params.follow_root = true;
  }
  // allow user to explicitly turn off following which may help in bizarre cases
  // in this case prop_root is used to determine which residuals are used
  if (ip_exist("CC3_FOLLOW_ROOT",0)) {
    i = 0;
    errcod = ip_boolean("CC3_FOLLOW_ROOT",&(i),0);
    if (errcod == IPE_OK) eom_params.follow_root = i;
  }

  /* so far, all R's are always kept so this is not used */
  eom_params.save_all = 0;

  eom_params.mult = 1;
  eom_params.rhf_triplets = 0;
  errcod = ip_data("RHF_TRIPLETS","%d",&(eom_params.rhf_triplets),0);
  if (eom_params.rhf_triplets != 0) eom_params.mult = 3;

  eom_params.excitation_range = 2;
  errcod = ip_data("EXCITATION_RANGE","%d",&(eom_params.excitation_range),0);

  eom_params.print_singles = 0;
  errcod = ip_data("PRINT_SINGLES","%d",&(eom_params.print_singles),0);

  eom_params.vectors_per_root_SS = 5;
  errcod = ip_data("VECTORS_PER_ROOT_SS","%d",&(eom_params.vectors_per_root_SS),0);

  eom_params.vectors_per_root = 12;
  errcod = ip_data("VECTORS_PER_ROOT","%d",&(eom_params.vectors_per_root),0);

  eom_params.vectors_cc3 = eom_params.prop_root + 9;
  errcod = ip_data("VECTORS_CC3","%d",&(eom_params.vectors_cc3),0);
  if (eom_params.vectors_cc3 > eom_params.vectors_per_root)
    eom_params.vectors_cc3 = eom_params.vectors_per_root;

  eom_params.collapse_with_last = 1;
  errcod = ip_boolean("COLLAPSE_WITH_LAST",&(eom_params.collapse_with_last),0);

  eom_params.complex_tol = 1E-12;
  errcod = ip_data("COMPLEX_TOL","%d",&(iconv),0);
  if(errcod == IPE_OK) eom_params.complex_tol = 1.0*pow(10.0,(double) -iconv);

  eom_params.residual_tol = 1E-6;
  errcod = ip_data("RESIDUAL_TOL","%d",&(iconv),0);
  if(errcod == IPE_OK) eom_params.residual_tol = 1.0*pow(10.0,(double) -iconv);

  eom_params.residual_tol_SS = 1E-6;
  errcod = ip_data("RESIDUAL_TOL_SS","%d",&(iconv),0);
  if(errcod == IPE_OK) eom_params.residual_tol_SS = 1.0*pow(10.0,(double) -iconv);

  eom_params.eval_tol = 1E-8;
  errcod = ip_data("EVAL_TOL","%d",&(iconv),0);
  if(errcod == IPE_OK) eom_params.eval_tol = 1.0*pow(10.0,(double) -iconv);

  eom_params.eval_tol_SS = 1E-6;
  errcod = ip_data("EVAL_TOL_SS","%d",&(iconv),0);
  if(errcod == IPE_OK) eom_params.eval_tol_SS = 1.0*pow(10.0,(double) -iconv);

  eom_params.amps_to_print = 5;
  errcod = ip_data("AMPS_TO_PRINT","%d",&i,0);
  if(errcod == IPE_OK) eom_params.amps_to_print = i;

  eom_params.schmidt_add_residual_tol = 1E-3;
  errcod = ip_data("SCHMIDT_ADD_RESIDUAL_TOL","%d",&(iconv),0);
  if(errcod == IPE_OK) eom_params.schmidt_add_residual_tol = 1.0*pow(10.0,(double) -iconv);

  eom_params.skip_diagSS = 0;
  i = 0;
  errcod = ip_boolean("SKIP_DIAGSS",&(i),0);
  if (i) eom_params.skip_diagSS = 1;

  eom_params.restart_eom_cc3 = 0;
  i = 0;
  errcod = ip_boolean("RESTART_EOM_CC3",&(i),0);
  if (i) eom_params.restart_eom_cc3 = 1;

  eom_params.max_iter_SS = 500;

  if(eom_params.guess == NULL) { /* we didn't use the cmdline arg --reuse */
    if(ip_exist("EOM_GUESS",0))
      errcod = ip_string("EOM_GUESS", &(eom_params.guess), 0);
    else
      eom_params.guess = strdup("SINGLES");
  }

  fprintf(outfile, "\n\tCCEOM parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tStates sought per irrep     =");
  if(strcmp(eom_params.guess,"INPUT")) {
    for (i=0;i<moinfo.nirreps;++i) fprintf(outfile, " %s %d,", moinfo.irr_labs[i],
					   eom_params.states_per_irrep[i]);
  }
  else fprintf(outfile, " USER INPUT");

  fprintf(outfile,"\n");
  fprintf(outfile, "\tMax. number of iterations   = %5d\n", eom_params.max_iter);
  fprintf(outfile, "\tVectors stored per root     = %5d\n", eom_params.vectors_per_root);
  fprintf(outfile, "\tPrint HbarSS iterations?    = %5d\n", eom_params.print_singles);
  fprintf(outfile, "\tExcitation range for HBarSS = %5d\n", eom_params.excitation_range);
  fprintf(outfile, "\tEigenvalue tolerance        = %5.1e\n", eom_params.eval_tol);
  fprintf(outfile, "\tEigenvalue toleranceSS      = %5.1e\n", eom_params.eval_tol_SS);
  fprintf(outfile, "\tResidual vector tolerance   = %5.1e\n", eom_params.residual_tol);
  fprintf(outfile, "\tResidual vector toleranceSS = %5.1e\n", eom_params.residual_tol_SS);
  fprintf(outfile, "\tComplex tolerance           = %5.1e\n", eom_params.complex_tol);
  fprintf(outfile, "\tRoot for properties         = %5d\n", eom_params.prop_root + 1);
  fprintf(outfile, "\tSym of state for properties = %6s\n", moinfo.irr_labs[eom_params.prop_sym]);
  fprintf(outfile, "\tGuess vectors taken from    = %s\n", eom_params.guess);
  fprintf(outfile, "\tRestart EOM CC3             = %s\n", eom_params.restart_eom_cc3?"YES":"NO");
  fprintf(outfile, "\tCollapse with last vector   = %s\n", eom_params.collapse_with_last ? "YES":"NO");
  if (eom_params.follow_root)
    fprintf(outfile, "\tRoot following for CC3 turned on.\n");
  fprintf(outfile, "\n\n");
}


}} // namespace psi::cceom
