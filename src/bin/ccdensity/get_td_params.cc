/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libdpd/dpd.h>
#include <libipv1/ip_lib.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void get_td_params(void)
{
  int i,j,k,l;
  char lbl[32];

  params.nstates = 0;

  if(ip_exist("PROP_SYM",0) && ip_exist("PROP_ROOT",0)) {
    ip_data("PROP_SYM","%d",&(params.prop_sym),0);
    ip_data("PROP_ROOT","%d",&(params.prop_root),0);
    /*User input counts from 1*/ 
    params.prop_sym -= 1;  
    params.prop_root -= 1; 
    params.nstates = 1;
  }
  else if(ip_exist("STATES_PER_IRREP",0)) {
    ip_count("STATES_PER_IRREP", &i, 0);
    if(i != moinfo.nirreps) {
      fprintf(outfile,"Dim. of states_per_irrep vector must be %d\n",
              moinfo.nirreps) ;
      exit(0);
    }
    for(i=0;i<moinfo.nirreps;++i) {
      ip_data("STATES_PER_IRREP","%d",&j,1,i);
      params.nstates += j;
    }
  }
  else {
    fprintf(outfile,"\nUse STATES_PER_IRREP or PROP_SYM and PROP_ROOT\n");
    exit(0);
  }

  /*
  fprintf(outfile,"\tNumber of States = %d\n",params.nstates);
  fflush(outfile);
  */

  td_params = (struct TD_Params *)malloc(params.nstates*sizeof(struct TD_Params));

  l=0; 
  if(ip_exist("PROP_SYM",0) && ip_exist("PROP_ROOT",0)) {
    td_params[0].irrep = params.prop_sym ^ moinfo.sym;
    k = td_params[0].root = params.prop_root;

    if(!strcmp(params.wfn,"CC2") || !strcmp(params.wfn,"EOM_CC2")) {
      sprintf(lbl,"EOM CC2 Energy for root %d %d", td_params[0].irrep, k);
      psio_read_entry(CC_INFO,lbl,(char*)&(td_params[0].cceom_energy),
                      sizeof(double));
      sprintf(lbl,"EOM CC2 R0 for root %d %d",td_params[0].irrep, k);
      psio_read_entry(CC_INFO,lbl,(char*)&(td_params[0].R0),sizeof(double));
    }
    else if(!strcmp(params.wfn,"CCSD") || !strcmp(params.wfn,"EOM_CCSD")) {
      sprintf(lbl,"EOM CCSD Energy for root %d %d", td_params[0].irrep, k);
      psio_read_entry(CC_INFO,lbl,(char*)&(td_params[0].cceom_energy),
                      sizeof(double));
      sprintf(lbl,"EOM CCSD R0 for root %d %d",td_params[0].irrep, k);
      psio_read_entry(CC_INFO,lbl,(char*)&(td_params[0].R0),sizeof(double));
    }
    else if(!strcmp(params.wfn,"CC3") || !strcmp(params.wfn,"EOM_CC3")) {
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
  else if(ip_exist("STATES_PER_IRREP",0)) {
    for(i=0;i<moinfo.nirreps;++i) {
      ip_data("STATES_PER_IRREP","%d",&j,1,i);
      for (k=0;k<j;++k) {
        td_params[l].irrep = i^moinfo.sym;
        td_params[l].root = k;

        if(!strcmp(params.wfn,"CC2") || !strcmp(params.wfn,"EOM_CC2")) {
          sprintf(lbl,"EOM CC2 Energy for root %d %d", td_params[l].irrep, k);
          psio_read_entry(CC_INFO,lbl,(char*)&(td_params[l].cceom_energy),
                        sizeof(double));
          sprintf(lbl,"EOM CC2 R0 for root %d %d",td_params[l].irrep, k);
          psio_read_entry(CC_INFO,lbl,(char*)&(td_params[l].R0),sizeof(double));
        }
        else if(!strcmp(params.wfn,"CCSD") || !strcmp(params.wfn,"EOM_CCSD")) {
          sprintf(lbl,"EOM CCSD Energy for root %d %d", td_params[l].irrep, k);
          psio_read_entry(CC_INFO,lbl,(char*)&(td_params[l].cceom_energy),
                        sizeof(double));
          sprintf(lbl,"EOM CCSD R0 for root %d %d",td_params[l].irrep, k);
          psio_read_entry(CC_INFO,lbl,(char*)&(td_params[l].R0),sizeof(double));
        }
        else if(!strcmp(params.wfn,"CC3") || !strcmp(params.wfn,"EOM_CC3")) {
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
