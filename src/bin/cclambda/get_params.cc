/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cmath>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include <libqt/qt.h>
#include <liboptions/liboptions.h>
#include <libchkpt/chkpt.h>
#include <psi4-dec.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void get_params(Options& options)
{
  int errcod, iconv,i,j,k,l,prop_sym,prop_root, excited_method=0;
        int *states_per_irrep, prop_all, lambda_and_Ls = 0;
  char lbl[32];
  std::string junk;

  /* check WFN keyword in input */
  params.wfn = options.get_str("WFN");
  excited_method = cc_excited(params.wfn);

  if(params.wfn == "CC2" || params.wfn == "EOM_CC2") {
    psio_read_entry(CC_INFO, "CC2 Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,  "\tCC2 energy          (CC_INFO) = %20.15f\n",moinfo.ecc);
    fprintf(outfile,  "\tTotal CC2 energy    (CC_INFO) = %20.15f\n",
            moinfo.eref+moinfo.ecc);
  }
  else if(params.wfn == "CCSD" || params.wfn == "EOM_CCSD") {
    psio_read_entry(CC_INFO, "CCSD Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,  "\tCCSD energy         (CC_INFO) = %20.15f\n",moinfo.ecc);
    fprintf(outfile,  "\tTotal CCSD energy   (CC_INFO) = %20.15f\n",
            moinfo.eref+moinfo.ecc);
  }
  else if(params.wfn == "CC3" || params.wfn == "EOM_CC3") {
    psio_read_entry(CC_INFO, "CC3 Energy", (char *) &(moinfo.ecc),
                    sizeof(double));
    fprintf(outfile,  "\tCC3 energy          (CC_INFO) = %20.15f\n",moinfo.ecc);
    fprintf(outfile,  "\tTotal CC3 energy    (CC_INFO) = %20.15f\n",
            moinfo.eref+moinfo.ecc);
  }

  /* read in the easy-to-understand parameters */

  params.convergence = 1e-7;
  params.convergence = options.get_double("R_CONVERGENCE");

  params.restart = 1;
  params.restart = options.get_bool("RESTART");
  if(!moinfo.phase) params.restart = 0;

  params.memory = Process::environment.get_memory();

  params.print = 0;
  params.print = options.get_int("PRINT");

  params.cachelev = 2;
  params.cachelev  = options.get_int("CACHELEVEL");

  params.sekino = 0;
  params.sekino = options.get_bool("SEKINO");

  params.diis = 1;
  params.diis = options.get_bool("DIIS");

  params.aobasis = 0;
  params.aobasis = options.get_bool("AO_BASIS");
  params.aobasis = 0;  /* AO basis code not yet working for lambda */

  params.abcd = options.get_str("ABCD");
  if(params.abcd == "NEW" && params.abcd == "OLD") {
    fprintf(outfile, "Invalid ABCD algorithm: %s\n", params.abcd.c_str());
    throw PsiException("cclambda: error", __FILE__, __LINE__);
  }

  params.num_amps = 10;
   params.num_amps = options.get_int("NUM_AMPS_PRINT");

  /* Determine DERTYPE */
  params.dertype = 0;
  if(options["DERTYPE"].has_changed()) {
    junk = options.get_str("DERTYPE");
    if(junk == "NONE") params.dertype = 0;
    else if(junk == "FIRST") params.dertype = 1;
    else if(junk == "RESPONSE") params.dertype = 3; /* linear response */
    else {
      printf("Invalid value of input keyword DERTYPE: %s\n", junk.c_str());
      throw PsiException("cclambda: error", __FILE__, __LINE__);
    }
  }

  /* begin local parameters */
  params.local = 0;
  params.local = options.get_bool("LOCAL");
  local.cutoff = 0.02;
  local.cutoff = options.get_double("LOCAL_CUTOFF");
  if(options["LOCAL_METHOD"].has_changed()) {
    local.method = options.get_str("LOCAL_METHOD");
    if(local.method == "AOBASIS" && local.method == "WERNER") {
      fprintf(outfile, "Invalid local correlation method: %s\n", local.method.c_str());
      throw PsiException("cclambda: error", __FILE__, __LINE__);
    }
  }
  else if(params.local) {
    local.method = "WERNER";
  }

  if(options["LOCAL_WEAKP"].has_changed()) {
    local.weakp = options.get_str("LOCAL_WEAKP");
    if(local.weakp != "MP2" && local.weakp != "NEGLECT" && local.weakp != "NONE") {
      fprintf(outfile, "Invalid method for treating local pairs: %s\n", local.weakp.c_str());
      throw PsiException("cclambda: error", __FILE__, __LINE__);
    }
  }
  else if(params.local) {
    local.weakp = "NONE";
  }

  if(params.dertype == 3)
    local.filter_singles = 0;
  else
    local.filter_singles = 1;

  local.filter_singles = options.get_bool("LOCAL_FILTER_SINGLES");

  local.cphf_cutoff = 0.10;
  local.cphf_cutoff = options.get_double("LOCAL_CPHF_CUTOFF");

  local.freeze_core = "FALSE";
  local.freeze_core = options.get_str("FREEZE_CORE");

  if(options["LOCAL_PAIRDEF"].has_changed()){
    local.pairdef = options.get_str("LOCAL_PAIRDEF");
    if(local.pairdef != "BP" && local.pairdef != "RESPONSE") {
      fprintf(outfile, "Invalid keyword for strong/weak pair definition: %s\n", local.pairdef.c_str());
      throw PsiException("cclambda: error", __FILE__, __LINE__);
    }
  }
  else if(params.local && params.dertype == 3)
    local.pairdef = "RESPONSE";
  else if(params.local)
    local.pairdef = "BP";

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
      states_per_irrep = options.get_int_array("ROOTS_PER_IRREP");
    }
    chkpt_close();

    prop_all = 1;
    prop_all = options.get_bool("PROP_ALL");
    /* command-line overrides this keyword (at least for now) */
    if (params.all) prop_all = 1;

    if (options["PROP_SYM"].has_changed()) {  /* read symmetry of state for properties */
      prop_sym = options.get_int("PROP_SYM");
      prop_sym -= 1;
      prop_sym = moinfo.sym^prop_sym;
          }
    else { /* just use last irrep of states requested for symmetry of states */
      for (i=0;i<moinfo.nirreps;++i) {
          if (states_per_irrep[i] > 0)
          prop_sym = i^moinfo.sym;
      }
    }

    if (options["PROP_ROOT"].has_changed()) { /* read prop_root */
      prop_root = options.get_int("PROP_ROOT");
      prop_root -= 1;
    }
    else { /* just use highest root, if you need only one of them */
            prop_root = states_per_irrep[prop_sym^moinfo.sym];
                        prop_root -= 1;
    }
  }

  params.zeta = options.get_bool("ZETA");

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
      if(params.wfn == "EOM_CC2") {
        sprintf(lbl,"EOM CC2 Energy for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].cceom_energy),sizeof(double));
        sprintf(lbl,"EOM CC2 R0 for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].R0),sizeof(double));
      }
      else if(params.wfn == "EOM_CCSD") {
        sprintf(lbl,"EOM CCSD Energy for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].cceom_energy),sizeof(double));
        sprintf(lbl,"EOM CCSD R0 for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[0].R0),sizeof(double));
      }
      else if(params.wfn == "EOM_CC3") {
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

            if(params.wfn == "EOM_CC2") {
              sprintf(lbl,"EOM CC2 Energy for root %d %d", i, j);
              psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[k].cceom_energy),sizeof(double));
              sprintf(lbl,"EOM CC2 R0 for root %d %d", i, j);
              psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[k].R0),sizeof(double));
            }
            else if(params.wfn == "EOM_CCSD") {
              sprintf(lbl,"EOM CCSD Energy for root %d %d", i, j);
              psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[k].cceom_energy),sizeof(double));
              sprintf(lbl,"EOM CCSD R0 for root %d %d", i, j);
              psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[k].R0),sizeof(double));
            }
            else if(params.wfn == "EOM_CC3") {
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

      if(params.wfn == "EOM_CC2") {
        sprintf(lbl,"EOM CC2 Energy for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[1].cceom_energy),sizeof(double));
        sprintf(lbl,"EOM CC2 R0 for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[1].R0),sizeof(double));
      }
      else if(params.wfn == "EOM_CCSD") {
        sprintf(lbl,"EOM CCSD Energy for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[1].cceom_energy),sizeof(double));
        sprintf(lbl,"EOM CCSD R0 for root %d %d", prop_sym, prop_root);
        psio_read_entry(CC_INFO, lbl, (char *) &(pL_params[1].R0),sizeof(double));
      }
      else if(params.wfn == "EOM_CC3") {
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
  params.maxiter = options.get_int("MAXITER");

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
  fprintf(outfile, "\tABCD            =     %s\n", params.abcd.c_str());
  fprintf(outfile, "\tLocal CC        =     %s\n", params.local ? "Yes" : "No");
  if(params.local) {
    fprintf(outfile, "\tLocal Cutoff    = %3.1e\n", local.cutoff);
    fprintf(outfile, "\tLocal Method    =    %s\n", local.method.c_str());
    fprintf(outfile, "\tWeak pairs      =    %s\n", local.weakp.c_str());
    fprintf(outfile, "\tFilter singles  =    %s\n", local.filter_singles ? "Yes" : "No");
    fprintf(outfile, "\tLocal pairs       =    %s\n", local.pairdef.c_str());
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
