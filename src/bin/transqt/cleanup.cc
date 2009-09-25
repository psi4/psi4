/*! \file
    \ingroup TRANSQT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"

namespace psi { namespace transqt {

void destruct_evects(int nirreps, double ***evects);


void cleanup(void)
{
  int i;

  /* Free moinfo Arrays */
  free(moinfo.sopi);
  free(moinfo.orbspi);
  free(moinfo.clsdpi);
  free(moinfo.openpi);
  free(moinfo.virtpi);
  free(moinfo.sosym);
  free(moinfo.orbsym);
  free(moinfo.order);
  free(moinfo.order_alpha);
  free(moinfo.order_beta);
  free(moinfo.corr2pitz);
  free(moinfo.corr2pitz_a);
  free(moinfo.corr2pitz_b);
  free(moinfo.frdocc);
  free(moinfo.fruocc);
  free(moinfo.rstrdocc);
  free(moinfo.rstruocc);
  free(moinfo.sloc);
  free(moinfo.stype);
  free(moinfo.snuc);
  if (params.backtr) {
    free(moinfo.corr2pitz_nofzv);
    free(moinfo.corr2pitz_nofzv_a);
    free(moinfo.corr2pitz_nofzv_b);
  }
  free(moinfo.first_so);
  free(moinfo.last_so);
  free(moinfo.first);
  free(moinfo.last);
  free(moinfo.fstact);
  free(moinfo.lstact);
  for(i=0; i < moinfo.nirreps; i++)
    free(moinfo.labels[i]);
  free(moinfo.labels);
  if(!strcmp(params.ref,"UHF")) {
    destruct_evects(params.backtr ? moinfo.backtr_nirreps : moinfo.nirreps, 
		    moinfo.evects_alpha);
    destruct_evects(params.backtr ? moinfo.backtr_nirreps : moinfo.nirreps, 
		    moinfo.evects_beta);
  }
  else {
    destruct_evects(params.backtr ? moinfo.backtr_nirreps : moinfo.nirreps, 
		    moinfo.evects);
  }
  free(moinfo.active);
  if(!strcmp(params.ref,"UHF")) {
    free_block(moinfo.scf_vector_alpha);
    free_block(moinfo.scf_vector_beta);
  }
  else free_block(moinfo.scf_vector);
  /* free(moinfo.evals); */
  free(moinfo.oe_ints);
  if(!strcmp(params.ref,"UHF")) {
    free(moinfo.fzc_operator_alpha);
    free(moinfo.fzc_operator_beta);
  }
  else free(moinfo.fzc_operator);
  free(moinfo.S);
  if (params.reorder) free(params.moorder);
  free(moinfo.backtr_mo_first);
  free(moinfo.backtr_mo_last);
  free(moinfo.backtr_mo_fstact);
  free(moinfo.backtr_mo_lstact);
  free(moinfo.backtr_mo_orbspi);
  free(moinfo.backtr_mo_active);
  free(moinfo.backtr_ao_first);
  free(moinfo.backtr_ao_last);
  free(moinfo.backtr_ao_orbspi);
  free(moinfo.backtr_ao_orbsym);

  /* Free ioff Array */
  free(ioff);

  /* Free params Arrays */
  free(params.wfn);
}

}} // end namespace psi::transqt
