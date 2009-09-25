/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include <libqt/qt.h>
#include <libchkpt/chkpt.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void get_params(void)
{
  int errcod, iconv,i,j,k,l,prop_sym,prop_root, excited_method=0;
	int *states_per_irrep, prop_all, lambda_and_Ls = 0;
  char lbl[32];
  char *junk;

  /* check WFN keyword in input */
  errcod = ip_string("WFN", &(params.wfn), 0);
	excited_method = cc_excited(params.wfn);

  if(!strcmp(params.wfn,"CC2") || !strcmp(params.wfn,"EOM_CC2")) {
    psio_read_entry(CC_INFO, "CC2 Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,  "\tCC2 energy          (CC_INFO) = %20.15f\n",moinfo.ecc);
    fprintf(outfile,  "\tTotal CC2 energy    (CC_INFO) = %20.15f\n",
            moinfo.eref+moinfo.ecc);
  }
  else if(!strcmp(params.wfn,"CCSD") || !strcmp(params.wfn,"EOM_CCSD")) {
    psio_read_entry(CC_INFO, "CCSD Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,  "\tCCSD energy         (CC_INFO) = %20.15f\n",moinfo.ecc);
    fprintf(outfile,  "\tTotal CCSD energy   (CC_INFO) = %20.15f\n",
            moinfo.eref+moinfo.ecc);
  }
  else if(!strcmp(params.wfn,"CC3") || !strcmp(params.wfn,"EOM_CC3")) {
    psio_read_entry(CC_INFO, "CC3 Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,  "\tCC3 energy          (CC_INFO) = %20.15f\n",moinfo.ecc);
    fprintf(outfile,  "\tTotal CC3 energy    (CC_INFO) = %20.15f\n",
            moinfo.eref+moinfo.ecc);
  }

  /* read in the easy-to-understand parameters */

  params.convergence = 1e-7;
  errcod = ip_data("CONVERGENCE","%d",&(iconv),0);
  if(errcod == IPE_OK) params.convergence = 1.0*pow(10.0,(double) -iconv);

  params.restart = 1;
  errcod = ip_boolean("RESTART", &(params.restart),0);
  if(!moinfo.phase) params.restart = 0;

  fndcor(&(params.memory),infile,outfile);

  params.print = 0;
  errcod = ip_data("PRINT", "%d", &(params.print),0);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);

  params.sekino = 0;
  if(ip_exist("SEKINO",0)) {
    errcod = ip_boolean("SEKINO", &params.sekino, 0);
    if(errcod != IPE_OK) params.sekino = 0;
  }

  params.diis = 1;
  errcod = ip_boolean("DIIS", &params.diis, 0);

  params.aobasis = 0;
  errcod = ip_boolean("AO_BASIS", &(params.aobasis),0);
  params.aobasis = 0;  /* AO basis code not yet working for lambda */

  if(ip_exist("ABCD",0)) {
    errcod = ip_string("ABCD", &(params.abcd), 0);
    if(strcmp(params.abcd,"NEW") && strcmp(params.abcd,"OLD")) {
      fprintf(outfile, "Invalid ABCD algorithm: %s\n", params.abcd);
      exit(PSI_RETURN_FAILURE);
    }
  }
  else params.abcd = strdup("NEW");

  params.num_amps = 10;
  if(ip_exist("NUM_AMPS",0)) {
    errcod = ip_data("NUM_AMPS", "%d", &(params.num_amps), 0);
  }

  /* Determine DERTYPE */
  params.dertype = 0;
  if(ip_exist("DERTYPE",0)) {
    errcod = ip_string("DERTYPE", &(junk),0);
    if(errcod != IPE_OK) {
      printf("Invalid value of input keyword DERTYPE: %s\n", junk);
      exit(PSI_RETURN_FAILURE); 
		}
    else if(!strcmp(junk,"NONE")) params.dertype = 0;
    else if(!strcmp(junk,"FIRST")) params.dertype = 1;
    else if(!strcmp(junk,"RESPONSE")) params.dertype = 3; /* linear response */
    else {
      printf("Invalid value of input keyword DERTYPE: %s\n", junk);
      exit(PSI_RETURN_FAILURE); 
    }
    free(junk);
  }
	else { /* DERTYPE is absent, assume 1 if jobtype=opt; 0 if jobtype=oeprop */
    ip_string("JOBTYPE", &(junk),0);
    if(!strcmp(junk,"OEPROP")) params.dertype = 0;
    else if(!strcmp(junk,"OPT")) params.dertype = 1;
    else {
      printf("Don't know what to do with DERTYPE missing and jobtype: %s\n", junk);
      exit(PSI_RETURN_FAILURE); 
    }
    free(junk);
	}

  /* begin local parameters */
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

	/* Now setup the structure which determines what will be solved */ 
	/* if --zeta, use Xi and solve for Zeta */
	/* if (DERTYPE == FIRST) determine ground vs. excited from wfn.
	    if ground, do only lambda.
			if excited, compute only one L chosen as described below.
			*/
	/* if (DERTYPE == RESPONSE), determine ground vs. excited from wfn.
	    Compute lambda.
			if excited, also do L(s) chosen as described below */
	/* if (DERTYPE == NONE) determine ground vs. excited from wfn.
	    Compute lambda.
			if excited, also do L(s) chosen as described below */
/* To determine which L(s) to compute for multiple L(s):
	  Check PROP_ALL in input
		 - If (PROP_ALL == true), compute L for all excited states.
		 - If false, check PROP_SYM for irrep desired, and PROP_ROOT
			     for root desired, as in cceom. */
/* To determine which L(s) to compute for single L(s) 
		 - Check PROP_SYM for irrep desired, and PROP_ROOT
			     for root desired, as in cceom. */

  /* setup property variables for excited states */
  if (cc_excited(params.wfn)) {
    chkpt_init(PSIO_OPEN_OLD);
    if (chkpt_rd_override_occ()) {
      states_per_irrep = chkpt_rd_statespi();
    }
    else {
      ip_count("STATES_PER_IRREP", &i, 0);
	    states_per_irrep = (int *) malloc(moinfo.nirreps * sizeof(int));
      for (i=0;i<moinfo.nirreps;++i)
        errcod = ip_data("STATES_PER_IRREP","%d",&(states_per_irrep[i]),1,i);
    }
    chkpt_close();

	  prop_all = 0;
	  if (ip_exist("PROP_ALL",0)) ip_boolean("PROP_ALL",&prop_all,0);
		/* command-line overrides this keyword (at least for now) */
		if (params.all) prop_all = 1;

	  if (ip_exist("PROP_SYM",0)) {  /* read symmetry of state for properties */
      ip_data("PROP_SYM","%d",&(prop_sym),0);
      prop_sym -= 1;
      prop_sym = moinfo.sym^prop_sym;
	  }
    else { /* just use last irrep of states requested for symmetry of states */
      for (i=0;i<moinfo.nirreps;++i) {
		    if (states_per_irrep[i] > 0)
          prop_sym = i^moinfo.sym;
		  }
    }

    if (ip_exist("PROP_ROOT",0)) { /* read prop_root */
      ip_data("PROP_ROOT","%d",&(prop_root),0);
      prop_root -= 1;
	  }
	  else { /* just use highest root, if you need only one of them */
	    prop_root = states_per_irrep[prop_sym^moinfo.sym];
			prop_root -= 1;
	  }
	}
	

  if (params.zeta) { /* only use Xi to solve for Zeta */
	  params.nstates = 1;
    pL_params = (struct L_Params *) malloc(params.nstates * sizeof(struct L_Params));
    psio_read_entry(CC_INFO, "XI Irrep", (char *) &i,sizeof(int));
    fprintf(outfile,"\tIrrep of Zeta       (CC_INFO) = %d\n", i);
    pL_params[0].irrep = prop_sym = i; /* is this always A1? I forget */
    pL_params[0].root = prop_root = 0;
    pL_params[0].ground = 0;
    pL_params[0].cceom_energy = 0.0;
    pL_params[0].R0 = 0.0; /* <Zeta0|R0> = 0, since zeta_0 = 0 */
    sprintf(pL_params[0].L1A_lbl,"ZIA");
    sprintf(pL_params[0].L1B_lbl,"Zia");
    sprintf(pL_params[0].L2AA_lbl,"ZIJAB");
    sprintf(pL_params[0].L2BB_lbl,"Zijab");
    sprintf(pL_params[0].L2AB_lbl,"ZIjAb");
    sprintf(pL_params[0].L2RHF_lbl,"2ZIjAb - ZIjbA");
  }
	else if (params.dertype == 1) { /* analytic gradient, ignore prop_all */
	  if (!cc_excited(params.wfn)) { /* do only lambda for ground state */
	    params.nstates = 1;
      pL_params = (struct L_Params *) malloc(params.nstates * sizeof(struct L_Params));
      pL_params[0].irrep = 0;
      pL_params[0].root = -1;
      pL_params[0].ground = 1;
      pL_params[0].R0 = 1.0;
      pL_params[0].cceom_energy = 0.0;
      sprintf(pL_params[0].L1A_lbl,"LIA %d %d",0, -1);
      sprintf(pL_params[0].L1B_lbl,"Lia %d %d",0, -1);
      sprintf(pL_params[0].L2AA_lbl,"LIJAB %d %d",0, -1);
      sprintf(pL_params[0].L2BB_lbl,"Lijab %d %d",0, -1);
      sprintf(pL_params[0].L2AB_lbl,"LIjAb %d %d",0, -1);
      sprintf(pL_params[0].L2RHF_lbl,"2LIjAb - LIjbA %d %d",0, -1);
		}
		else { /* do only one L for excited state */
		  params.nstates = 1;
      pL_params = (struct L_Params *) malloc(params.nstates * sizeof(struct L_Params));
      pL_params[0].irrep = prop_sym;
      pL_params[0].root = prop_root;
      pL_params[0].ground = 0;
      if(!strcmp(params.wfn,"EOM_CC2")) {
        sprintf(lbl,"EOM CC2 Energy for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].cceom_energy),sizeof(double));
        sprintf(lbl,"EOM CC2 R0 for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].R0),sizeof(double));
      }
      else if(!strcmp(params.wfn,"EOM_CCSD")) {
        sprintf(lbl,"EOM CCSD Energy for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].cceom_energy),sizeof(double));
        sprintf(lbl,"EOM CCSD R0 for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].R0),sizeof(double));
      }
      else if(!strcmp(params.wfn,"EOM_CC3")) {
        sprintf(lbl,"EOM CC3 Energy for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].cceom_energy),sizeof(double));
        sprintf(lbl,"EOM CC3 R0 for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].R0),sizeof(double));
      }
      sprintf(pL_params[0].L1A_lbl,"LIA %d %d",prop_sym, prop_root);
      sprintf(pL_params[0].L1B_lbl,"Lia %d %d",prop_sym, prop_root);
      sprintf(pL_params[0].L2AA_lbl,"LIJAB %d %d",prop_sym, prop_root);
      sprintf(pL_params[0].L2BB_lbl,"Lijab %d %d",prop_sym, prop_root);
      sprintf(pL_params[0].L2AB_lbl,"LIjAb %d %d",prop_sym, prop_root);
      sprintf(pL_params[0].L2RHF_lbl,"2LIjAb - LIjbA %d %d",prop_sym, prop_root);
		}
  }
	else if (params.dertype == 3) { /* response calculation */
	  if (!cc_excited(params.wfn)) { /* ground state */
	    params.nstates = 1;
      pL_params = (struct L_Params *) malloc(params.nstates * sizeof(struct L_Params));
      pL_params[0].irrep = 0;
      pL_params[0].root = -1;
      pL_params[0].ground = 1;
      pL_params[0].R0 = 1.0;
      pL_params[0].cceom_energy = 0.0;
      sprintf(pL_params[0].L1A_lbl,"LIA %d %d",0, -1);
      sprintf(pL_params[0].L1B_lbl,"Lia %d %d",0, -1);
      sprintf(pL_params[0].L2AA_lbl,"LIJAB %d %d",0, -1);
      sprintf(pL_params[0].L2BB_lbl,"Lijab %d %d",0, -1);
      sprintf(pL_params[0].L2AB_lbl,"LIjAb %d %d",0, -1);
      sprintf(pL_params[0].L2RHF_lbl,"2LIjAb - LIjbA %d %d",0, -1);
		}
		else { /* excited state */
		  lambda_and_Ls = 1; /* code is below */
		}
	}
	else if (params.dertype == 0) {
	  if (!cc_excited(params.wfn)) { /* ground state */
	    params.nstates = 1;
      pL_params = (struct L_Params *) malloc(params.nstates * sizeof(struct L_Params));
      pL_params[0].irrep = 0;
      pL_params[0].root = -1;
      pL_params[0].ground = 1;
      pL_params[0].R0 = 1.0;
      pL_params[0].cceom_energy = 0.0;
      sprintf(pL_params[0].L1A_lbl,"LIA %d %d",0, -1);
      sprintf(pL_params[0].L1B_lbl,"Lia %d %d",0, -1);
      sprintf(pL_params[0].L2AA_lbl,"LIJAB %d %d",0, -1);
      sprintf(pL_params[0].L2BB_lbl,"Lijab %d %d",0, -1);
      sprintf(pL_params[0].L2AB_lbl,"LIjAb %d %d",0, -1);
      sprintf(pL_params[0].L2RHF_lbl,"2LIjAb - LIjbA %d %d",0, -1);
		}
		else { /* excited state */
		  lambda_and_Ls = 1; /* code is below */
		}
	}


  /* do lambda for ground state AND do L(s) for excited states */
  if (lambda_and_Ls) {
    /* determine number of states to converge */
	  params.nstates = 1; /* for ground state */
		if (prop_all) {
		  for (i=0; i<moinfo.nirreps; ++i) 
			  params.nstates += states_per_irrep[i]; /* do all L(s) */
		}
		else {
		  params.nstates += 1; /* do only one L */
		}

    pL_params = (struct L_Params *) malloc(params.nstates * sizeof(struct L_Params));

		/* ground state */
    pL_params[0].irrep = 0;
    pL_params[0].root = -1;
    pL_params[0].ground = 1;
    pL_params[0].R0 = 1.0;
    pL_params[0].cceom_energy = 0.0;
    sprintf(pL_params[0].L1A_lbl,"LIA %d %d",0, -1);
    sprintf(pL_params[0].L1B_lbl,"Lia %d %d",0, -1);
    sprintf(pL_params[0].L2AA_lbl,"LIJAB %d %d",0, -1);
    sprintf(pL_params[0].L2BB_lbl,"Lijab %d %d",0, -1);
    sprintf(pL_params[0].L2AB_lbl,"LIjAb %d %d",0, -1);
    sprintf(pL_params[0].L2RHF_lbl,"2LIjAb - LIjbA %d %d",0, -1);

		if (prop_all) { /* do all L(s) */
		  k=1;
		    for (i=0; i<moinfo.nirreps; ++i)  { /* look over irrep of L(s) */
				  for (j=0; j < states_per_irrep[i^moinfo.sym]; ++j) {
            pL_params[k].irrep = i;
            pL_params[k].root = j;
            pL_params[k].ground = 0;

            if(!strcmp(params.wfn,"EOM_CC2")) {
              sprintf(lbl,"EOM CC2 Energy for root %d %d", i, j);
              psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[k].cceom_energy),sizeof(double));
              sprintf(lbl,"EOM CC2 R0 for root %d %d", i, j);
              psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[k].R0),sizeof(double));
            }
            else if(!strcmp(params.wfn,"EOM_CCSD")) {
              sprintf(lbl,"EOM CCSD Energy for root %d %d", i, j);
              psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[k].cceom_energy),sizeof(double));
              sprintf(lbl,"EOM CCSD R0 for root %d %d", i, j);
              psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[k].R0),sizeof(double));
            }
            else if(!strcmp(params.wfn,"EOM_CC3")) {
              sprintf(lbl,"EOM CC3 Energy for root %d %d", i, j);
              psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[k].cceom_energy),sizeof(double));
              sprintf(lbl,"EOM CC3 R0 for root %d %d", i, j);
              psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[k].R0),sizeof(double));
            }

            sprintf(pL_params[k].L1A_lbl,"LIA %d %d",i, j);
            sprintf(pL_params[k].L1B_lbl,"Lia %d %d",i, j);
            sprintf(pL_params[k].L2AA_lbl,"LIJAB %d %d",i, j);
            sprintf(pL_params[k].L2BB_lbl,"Lijab %d %d",i, j);
            sprintf(pL_params[k].L2AB_lbl,"LIjAb %d %d",i, j);
            sprintf(pL_params[k].L2RHF_lbl,"2LIjAb - LIjbA %d %d",i, j);
						k++;
					}
			  }
		}
		else { /* use prop_sym and prop_root determined above from input or inferrence */
      pL_params[1].irrep = prop_sym;
      pL_params[1].root = prop_root;
      pL_params[1].ground = 0;

      if(!strcmp(params.wfn,"EOM_CC2")) {
        sprintf(lbl,"EOM CC2 Energy for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[1].cceom_energy),sizeof(double));
        sprintf(lbl,"EOM CC2 R0 for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[1].R0),sizeof(double));
      }
      else if(!strcmp(params.wfn,"EOM_CCSD")) {
        sprintf(lbl,"EOM CCSD Energy for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[1].cceom_energy),sizeof(double));
        sprintf(lbl,"EOM CCSD R0 for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[1].R0),sizeof(double));
      }
      else if(!strcmp(params.wfn,"EOM_CC3")) {
        sprintf(lbl,"EOM CC3 Energy for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[1].cceom_energy),sizeof(double));
        sprintf(lbl,"EOM CC3 R0 for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[1].R0),sizeof(double));
      }

      sprintf(pL_params[1].L1A_lbl,"LIA %d %d", prop_sym, prop_root);
      sprintf(pL_params[1].L1B_lbl,"Lia %d %d", prop_sym, prop_root);
      sprintf(pL_params[1].L2AA_lbl,"LIJAB %d %d", prop_sym, prop_root);
      sprintf(pL_params[1].L2BB_lbl,"Lijab %d %d", prop_sym, prop_root);
      sprintf(pL_params[1].L2AB_lbl,"LIjAb %d %d", prop_sym, prop_root);
      sprintf(pL_params[1].L2RHF_lbl,"2LIjAb - LIjbA %d %d", prop_sym, prop_root);
		}
  }


  params.maxiter = 50 * params.nstates;
  errcod = ip_data("MAXITER","%d",&(params.maxiter),0);

  fprintf(outfile, "\n\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tMaxiter       =    %4d\n", params.maxiter);
  fprintf(outfile, "\tConvergence   = %3.1e\n", params.convergence);
  fprintf(outfile, "\tRestart       =     %s\n", params.restart ? "Yes" : "No");
  fprintf(outfile, "\tCache Level   =     %1d\n", params.cachelev);
  fprintf(outfile, "\tModel III     =     %s\n", params.sekino ? "Yes" : "No");
  fprintf(outfile, "\tDIIS          =     %s\n", params.diis ? "Yes" : "No");
  fprintf(outfile, "\tAO Basis      =     %s\n", 
          params.aobasis ? "Yes" : "No");
  fprintf(outfile, "\tABCD            =     %s\n", params.abcd);
  fprintf(outfile, "\tLocal CC        =     %s\n", params.local ? "Yes" : "No");
  if(params.local) {
    fprintf(outfile, "\tLocal Cutoff    = %3.1e\n", local.cutoff);
    fprintf(outfile, "\tLocal Method    =    %s\n", local.method);
    fprintf(outfile, "\tWeak pairs      =    %s\n", local.weakp);
    fprintf(outfile, "\tFilter singles  =    %s\n", local.filter_singles ? "Yes" : "No");
    fprintf(outfile, "\tLocal pairs       =    %s\n", local.pairdef);
    fprintf(outfile, "\tLocal CPHF cutoff =  %3.1e\n", local.cphf_cutoff);
  }

  fprintf(outfile,"\tParamaters for left-handed eigenvectors:\n");
  fprintf(outfile,"\t    Irr   Root  Ground-State?    EOM energy        R0\n");
  for (i=0; i<params.nstates; ++i) {
    fprintf(outfile,"\t%3d %3d %5d %10s %18.10lf %14.10lf\n", i+1, pL_params[i].irrep, pL_params[i].root+1,
	    (pL_params[i].ground ? "Yes":"No"), pL_params[i].cceom_energy, pL_params[i].R0);
  }

  for (i=0; i<params.nstates; ++i) {
    fprintf(outfile,"\tLabels for eigenvector %d:\n\t%s, %s, %s, %s, %s, %s\n",
	    i+1,pL_params[i].L1A_lbl,pL_params[i].L1B_lbl,pL_params[i].L2AA_lbl,pL_params[i].L2BB_lbl,
	    pL_params[i].L2AB_lbl, pL_params[i].L2RHF_lbl);
  }

  fflush(outfile);
  return;
}

}} // namespace psi::cclambda
