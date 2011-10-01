/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libdpd/dpd.h>
#include <liboptions/liboptions.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void get_td_params(Options& options)
{
  int i,j,k,l;
  char lbl[32];

  params.nstates = 0;

  if(options["PROP_SYM"].has_changed() && options["PROP_ROOT"].has_changed()) {
    params.prop_sym = options.get_int("PROP_SYM");
    params.prop_root = options.get_int("PROP_ROOT");
    /*User input counts from 1*/
    params.prop_sym -= 1;
    params.prop_root -= 1;
    params.nstates = 1;
  }
  else if(options["STATES_PER_IRREP"].has_changed()) {
    i = options["STATES_PER_IRREP"].size();
    if(i != moinfo.nirreps) {
      fprintf(outfile,"Dim. of states_per_irrep vector must be %d\n",
              moinfo.nirreps) ;
      throw PsiException("ccdensity: error", __FILE__, __LINE__);
    }
    for(i=0;i<moinfo.nirreps;++i) {
      params.nstates += options["STATES_PER_IRREP"][i].to_integer();
    }
  }
  else {
    fprintf(outfile,"\nUse STATES_PER_IRREP or PROP_SYM and PROP_ROOT\n");
    throw PsiException("ccdensity: error", __FILE__, __LINE__);
  }

  /*
  fprintf(outfile,"\tNumber of States = %d\n",params.nstates);
  fflush(outfile);
  */

  td_params = (struct TD_Params *)malloc(params.nstates*sizeof(struct TD_Params));

  l=0;
  if(options["PROP_SYM"].has_changed() && options["PROP_ROOT"].has_changed()) {
    td_params[0].irrep = params.prop_sym ^ moinfo.sym;
    k = td_params[0].root = params.prop_root;

    if(params.wfn == "CC2" || params.wfn == "EOM_CC2") {
      sprintf(lbl,"EOM CC2 Energy for root %d %d", td_params[0].irrep, k);
      psio_read_entry(CC_INFO,lbl,(char*)&(td_params[0].cceom_energy),
                      sizeof(double));
      sprintf(lbl,"EOM CC2 R0 for root %d %d",td_params[0].irrep, k);
      psio_read_entry(CC_INFO,lbl,(char*)&(td_params[0].R0),sizeof(double));
    }
    else if(params.wfn == "CCSD" || params.wfn == "EOM_CCSD") {
      sprintf(lbl,"EOM CCSD Energy for root %d %d", td_params[0].irrep, k);
      psio_read_entry(CC_INFO,lbl,(char*)&(td_params[0].cceom_energy),
                      sizeof(double));
      sprintf(lbl,"EOM CCSD R0 for root %d %d",td_params[0].irrep, k);
      psio_read_entry(CC_INFO,lbl,(char*)&(td_params[0].R0),sizeof(double));
    }
    else if(params.wfn == "CC3" || params.wfn == "EOM_CC3") {
      sprintf(lbl,"EOM CC3 Energy for root %d %d", td_params[0].irrep, k);
      psio_read_entry(CC_INFO,lbl,(char*)&(td_params[0].cceom_energy),
                      sizeof(double));
      sprintf(lbl,"EOM CC3 R0 for root %d %d",td_params[0].irrep, k);
      psio_read_entry(CC_INFO,lbl,(char*)&(td_params[0].R0),sizeof(double));
    }

    sprintf(td_params[l].L1A_lbl,"LIA %d %d",td_params[0].irrep, k);
    sprintf(td_params[l].L1B_lbl,"Lia %d %d",td_params[0].irrep, k);
    sprintf(td_params[l].L2AA_lbl,"LIJAB %d %d",td_params[0].irrep, k);
    sprintf(td_params[l].L2BB_lbl,"Lijab %d %d",td_params[0].irrep, k);
    sprintf(td_params[l].L2AB_lbl,"LIjAb %d %d",td_params[0].irrep, k);
    sprintf(td_params[l].R1A_lbl,"RIA %d %d",td_params[0].irrep, k);
    sprintf(td_params[l].R1B_lbl,"Ria %d %d",td_params[0].irrep, k);
    sprintf(td_params[l].R2AA_lbl,"RIJAB %d %d",td_params[0].irrep, k);
    sprintf(td_params[l].R2BB_lbl,"Rijab %d %d",td_params[0].irrep, k);
    sprintf(td_params[l].R2AB_lbl,"RIjAb %d %d",td_params[0].irrep, k);
  }
  else if(options["STATES_PER_IRREP"].has_changed()) {
    for(i=0;i<moinfo.nirreps;++i) {
      j = options["STATES_PER_IRREP"][i].to_integer();
      for (k=0;k<j;++k) {
        td_params[l].irrep = i^moinfo.sym;
        td_params[l].root = k;

        if(params.wfn == "CC2" || params.wfn == "EOM_CC2") {
          sprintf(lbl,"EOM CC2 Energy for root %d %d", td_params[l].irrep, k);
          psio_read_entry(CC_INFO,lbl,(char*)&(td_params[l].cceom_energy),
                        sizeof(double));
          sprintf(lbl,"EOM CC2 R0 for root %d %d",td_params[l].irrep, k);
          psio_read_entry(CC_INFO,lbl,(char*)&(td_params[l].R0),sizeof(double));
        }
        else if(params.wfn == "CCSD" || params.wfn == "EOM_CCSD") {
          sprintf(lbl,"EOM CCSD Energy for root %d %d", td_params[l].irrep, k);
          psio_read_entry(CC_INFO,lbl,(char*)&(td_params[l].cceom_energy),
                        sizeof(double));
          sprintf(lbl,"EOM CCSD R0 for root %d %d",td_params[l].irrep, k);
          psio_read_entry(CC_INFO,lbl,(char*)&(td_params[l].R0),sizeof(double));
        }
        else if(params.wfn == "CC3" || params.wfn == "EOM_CC3") {
          sprintf(lbl,"EOM CC3 Energy for root %d %d", td_params[l].irrep, k);
          psio_read_entry(CC_INFO,lbl,(char*)&(td_params[l].cceom_energy),
                        sizeof(double));
          sprintf(lbl,"EOM CC3 R0 for root %d %d",td_params[l].irrep, k);
          psio_read_entry(CC_INFO,lbl,(char*)&(td_params[l].R0),sizeof(double));
        }

        sprintf(td_params[l].L1A_lbl,"LIA %d %d",td_params[l].irrep, k);
        sprintf(td_params[l].L1B_lbl,"Lia %d %d",td_params[l].irrep, k);
        sprintf(td_params[l].L2AA_lbl,"LIJAB %d %d",td_params[l].irrep, k);
        sprintf(td_params[l].L2BB_lbl,"Lijab %d %d",td_params[l].irrep, k);
        sprintf(td_params[l].L2AB_lbl,"LIjAb %d %d",td_params[l].irrep, k);
        sprintf(td_params[l].R1A_lbl,"RIA %d %d",td_params[l].irrep, k);
        sprintf(td_params[l].R1B_lbl,"Ria %d %d",td_params[l].irrep, k);
        sprintf(td_params[l].R2AA_lbl,"RIJAB %d %d",td_params[l].irrep, k);
        sprintf(td_params[l].R2BB_lbl,"Rijab %d %d",td_params[l].irrep, k);
        sprintf(td_params[l].R2AB_lbl,"RIjAb %d %d",td_params[l].irrep, k);
        l++;
      }
    }
  }
  /*
  fprintf(outfile,"\n\tState\t  EOM Energy\t    R0\n");
  for(i=0; i<params.nstates; i++) {
    fprintf(outfile,"\t %d%3s %15.10lf %12.8lf\n",
            td_params[i].root+1,moinfo.labels[td_params[i].irrep],
            td_params[i].cceom_energy,td_params[i].R0);
  }
  */

  return;
}

}} // namespace psi::ccdensity
