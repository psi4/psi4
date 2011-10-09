/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <libchkpt/chkpt.h>
#include <liboptions/liboptions.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void get_rho_params(Options& options)
{
  int i,j,k,l,prop_sym,prop_root, lambda_and_Ls=0, errcod, prop_all,cnt;
  char lbl[32];
  int *states_per_irrep;

  /* check WFN keyword in input */
  params.wfn = options.get_str("WFN");
  if ( !cc_excited(params.wfn) ) params.ground = 1;
  else params.ground = 0;

  /* setup propery variables for excited states */
  if (!params.ground) {
    chkpt_init(PSIO_OPEN_OLD);
    if (chkpt_rd_override_occ()) {
      states_per_irrep = chkpt_rd_statespi();
    }
    else {
      states_per_irrep = options.get_int_array("STATES_PER_IRREP");
    }
    chkpt_close();

    prop_all = 0;
    prop_all = options.get_bool("PROP_ALL");
    /* can also be turned on by command-line (at least for now) */
    if (params.prop_all) prop_all = 1;
    if (params.dertype == 1) prop_all = 0; /* don't do relaxed, multiple excited states */
    params.prop_all = prop_all;

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

  /* setup density parameters */
  if (params.ground) { /* just compute ground state density */
    params.nstates = 1;
    rho_params = (struct RHO_Params *) malloc(params.nstates * sizeof(struct RHO_Params));
    rho_params[0].L_irr = 0;
    rho_params[0].R_irr = 0;
    rho_params[0].L_root = -1;
    rho_params[0].R_root = -1;
    rho_params[0].L_ground = 1;
    rho_params[0].R_ground = 1;
    rho_params[0].R0 = 1;
    rho_params[0].L0 = 0;
    rho_params[0].cceom_energy = 0;
  }
  else if (params.calc_xi) { /* just compute Xi for excited-state density */
    params.nstates = 1;
    rho_params = (struct RHO_Params *) malloc(params.nstates * sizeof(struct RHO_Params));
    rho_params[0].L_irr = prop_sym;
    rho_params[0].R_irr = prop_sym;
    rho_params[0].L_root = prop_root;
    rho_params[0].R_root = prop_root;
    rho_params[0].L_ground = 0;
    rho_params[0].R_ground = 0;
  }
  else { /* excited state density(ies) are involved */
    if (params.prop_all)  { /* do all roots */
      params.nstates = 1;
      for(i=0; i<moinfo.nirreps; i++) {
//        cnt = 0;
//        ip_data("STATES_PER_IRREP","%d",&cnt, 1, i);
        params.nstates += states_per_irrep[i];
      }
      rho_params = (struct RHO_Params *) malloc(params.nstates * sizeof(struct RHO_Params));
      rho_params[0].L_irr = 0;
      rho_params[0].R_irr = 0;
      rho_params[0].L_root = -1;
      rho_params[0].R_root = -1;
      rho_params[0].L_ground = 1;
      rho_params[0].R_ground = 1;

      cnt = 0;
      for(i=0; i<moinfo.nirreps; i++) { /* loop over irrep of R */
//        ip_data("STATES_PER_IRREP","%d",&j, 1, i);
//        for (k=0; k<j; ++k) {
        for (k=0; k<states_per_irrep[i]; ++k) {
          ++cnt;
          rho_params[cnt].L_irr = i^moinfo.sym;
          rho_params[cnt].R_irr = i^moinfo.sym;
          rho_params[cnt].L_root = k;
          rho_params[cnt].R_root = k;
          rho_params[cnt].L_ground = 0;
          rho_params[cnt].R_ground = 0;
        }
      }
    }
    else { /* do only one root */
      params.nstates = 1;
      rho_params = (struct RHO_Params *) malloc(params.nstates * sizeof(struct RHO_Params));
      rho_params[0].L_irr = prop_sym;
      rho_params[0].R_irr = prop_sym;
      rho_params[0].L_root = prop_root;
      rho_params[0].R_root = prop_root;
      rho_params[0].L_ground = 0;
      rho_params[0].R_ground = 0;
    }
    if(params.onepdm) params.transition = 1;
  }

  /* for each state, determine G_irr, cceom_energy, R0, and labels for files */
  for(i=0; i<params.nstates; i++) {
    rho_params[i].G_irr = (rho_params[i].L_irr) ^ (rho_params[i].R_irr);

    if (rho_params[i].R_ground) {
      rho_params[i].cceom_energy = 0.0;
      rho_params[i].R0 = 1.0;
    }
    else {
      if(params.wfn == "EOM_CC2") {
        sprintf(lbl,"EOM CC2 Energy for root %d %d", rho_params[i].R_irr, rho_params[i].R_root);
        psio_read_entry(CC_INFO,lbl,(char*)&(rho_params[i].cceom_energy), sizeof(double));
        sprintf(lbl,"EOM CC2 R0 for root %d %d",rho_params[i].R_irr, rho_params[i].R_root);
        psio_read_entry(CC_INFO,lbl,(char*)&(rho_params[i].R0),sizeof(double));
      }
      else if(params.wfn == "EOM_CCSD") {
        sprintf(lbl,"EOM CCSD Energy for root %d %d", rho_params[i].R_irr, rho_params[i].R_root);
        psio_read_entry(CC_INFO,lbl,(char*)&(rho_params[i].cceom_energy), sizeof(double));
        sprintf(lbl,"EOM CCSD R0 for root %d %d",rho_params[i].R_irr, rho_params[i].R_root);
        psio_read_entry(CC_INFO,lbl,(char*)&(rho_params[i].R0),sizeof(double));
      }
      else if( params.wfn == "EOM_CC3") {
        sprintf(lbl,"EOM CC3 Energy for root %d %d", rho_params[i].R_irr, rho_params[i].R_root);
        psio_read_entry(CC_INFO,lbl,(char*)&(rho_params[i].cceom_energy), sizeof(double));
        sprintf(lbl,"EOM CC3 R0 for root %d %d",rho_params[i].R_irr, rho_params[i].R_root);
        psio_read_entry(CC_INFO,lbl,(char*)&(rho_params[i].R0),sizeof(double));
      }
    }
    if (rho_params[i].L_ground)
      rho_params[i].L0 = 1.0;
    else
      rho_params[i].L0 = 0.0;

    sprintf(rho_params[i].R1A_lbl, "RIA %d %d",  rho_params[i].R_irr, rho_params[i].R_root);
    sprintf(rho_params[i].R1B_lbl, "Ria %d %d",  rho_params[i].R_irr, rho_params[i].R_root);
    sprintf(rho_params[i].R2AA_lbl,"RIJAB %d %d",rho_params[i].R_irr, rho_params[i].R_root);
    sprintf(rho_params[i].R2BB_lbl,"Rijab %d %d",rho_params[i].R_irr, rho_params[i].R_root);
    sprintf(rho_params[i].R2AB_lbl,"RIjAb %d %d",rho_params[i].R_irr, rho_params[i].R_root);
    sprintf(rho_params[i].L1A_lbl, "LIA %d %d",  rho_params[i].L_irr, rho_params[i].L_root);
    sprintf(rho_params[i].L1B_lbl, "Lia %d %d",  rho_params[i].L_irr, rho_params[i].L_root);
    sprintf(rho_params[i].L2AA_lbl,"LIJAB %d %d",rho_params[i].L_irr, rho_params[i].L_root);
    sprintf(rho_params[i].L2BB_lbl,"Lijab %d %d",rho_params[i].L_irr, rho_params[i].L_root);
    sprintf(rho_params[i].L2AB_lbl,"LIjAb %d %d",rho_params[i].L_irr, rho_params[i].L_root);

    sprintf(rho_params[i].DIJ_lbl, "DIJ %d", i-1); /* change to a different label ? */
    sprintf(rho_params[i].Dij_lbl, "Dij %d", i-1);
    sprintf(rho_params[i].DAB_lbl, "DAB %d", i-1);
    sprintf(rho_params[i].Dab_lbl, "Dab %d", i-1);
    sprintf(rho_params[i].DAI_lbl, "DAI %d", i-1);
    sprintf(rho_params[i].Dai_lbl, "Dai %d", i-1);
    sprintf(rho_params[i].DIA_lbl, "DIA %d", i-1);
    sprintf(rho_params[i].Dia_lbl, "Dia %d", i-1);

    if (params.ref == 0) {
      if (i == 0) sprintf(rho_params[i].opdm_lbl, "MO-basis OPDM");
      else sprintf(rho_params[i].opdm_lbl, "MO-basis OPDM Root %d", i);
    }
    else if (params.ref == 1) { /* ROHF */
      if (i == 0) sprintf(rho_params[i].opdm_lbl, "MO-basis OPDM");
      else sprintf(rho_params[i].opdm_lbl, "MO-basis OPDM Root %d", i);
    }
    else if (params.ref == 2) { /* UHF */
      if (i == 0) {
        sprintf(rho_params[i].opdm_a_lbl, "MO-basis Alpha OPDM");
        sprintf(rho_params[i].opdm_b_lbl, "MO-basis Beta OPDM");
      }
      else {
        sprintf(rho_params[i].opdm_a_lbl, "MO-basis Alpha OPDM Root %d", i);
        sprintf(rho_params[i].opdm_b_lbl, "MO-basis Beta OPDM Root %d", i);
      }
    }
  }

  psio_write_entry(CC_INFO, "Num. of CC densities", (char *) &(params.nstates), sizeof(int));

  fprintf(outfile,"\tNumber of States = %-d\n",params.nstates);

  /*
  fprintf(outfile,"\n\tGround? L_root L_irr R_root R_irr G_irr     EOM Energy        R0\n");
  for(i=0; i<params.nstates; i++) {
    fprintf(outfile,"\t%5s %6d %7s %4d %7s %4s %15.10lf %12.8lf\n",
            (rho_params[i].R_ground ? "Yes":"No"),
            rho_params[i].L_root+1, moinfo.labels[rho_params[i].L_irr],
            rho_params[i].R_root+1, moinfo.labels[rho_params[i].R_irr],
            moinfo.labels[rho_params[i].G_irr],
            rho_params[i].cceom_energy,
            rho_params[i].R0);
  }
  */

  fprintf(outfile,"\n\tGround?  State     EOM Energy       R0\n");
  for(i=0; i<params.nstates; i++) {
    fprintf(outfile,"\t%5s     %d%3s %15.10lf %12.8lf\n",
            (rho_params[i].R_ground ? "Yes":"No"),
            rho_params[i].L_root+1, moinfo.labels[rho_params[i].L_irr],
            rho_params[i].cceom_energy,
            rho_params[i].R0);
  }

  /* set variables for xi and twopdm code in the old non-state specific structure */
  params.G_irr = rho_params[0].G_irr;
  params.R_irr = rho_params[0].R_irr;
  params.L_irr = rho_params[0].L_irr;
  params.R0 = rho_params[0].R0;
  params.L0 = rho_params[0].L0;
  params.cceom_energy = rho_params[0].cceom_energy;

  return;
}

}} // namespace psi::ccdensity
