/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

namespace psi {

int dpd_close(int dpd_num)
{
  int h,i,j,k,cnt;
  int nirreps,num_subspaces,num_pairs;
  dpd_data *this_dpd;

  /*  dpd_file2_cache_print(stdout); */
  dpd_file2_cache_close();
  /*  dpd_file4_cache_print(stdout);*/
  dpd_file4_cache_close();

  this_dpd = &(dpd_list[dpd_num]);

  num_subspaces = this_dpd->num_subspaces;
  nirreps = this_dpd->nirreps;
  num_pairs = this_dpd->num_pairs;

  for(i=0; i < num_pairs; i++)
    for(j=0; j < num_pairs; j++)
      free_int_matrix(this_dpd->params4[i][j].start13);

  for(i=0; i < num_subspaces; i++)  free(this_dpd->orboff[i]);
  free(this_dpd->orboff);

  for(i=0; i < num_subspaces; i++) {
    for(j=0; j < 5; j++) {
      free_int_matrix(this_dpd->pairidx[5*i+j]);
      for(k=0; k < this_dpd->nirreps; k++) 
	if(this_dpd->pairtot[5*i+j][k])
	  free_int_matrix(this_dpd->pairorb[5*i+j][k]);
      free(this_dpd->pairorb[5*i+j]);
    }
  }
  for(i=0,cnt=5*num_subspaces; i < num_subspaces; i++) {
    for(j=i+1; j < num_subspaces; j++,cnt+=2) {
      free_int_matrix(this_dpd->pairidx[cnt]);
      free_int_matrix(this_dpd->pairidx[cnt+1]);
      for(k=0; k < nirreps; k++) {
	if(this_dpd->pairtot[cnt][k])
	  free_int_matrix(this_dpd->pairorb[cnt][k]);
	if(this_dpd->pairtot[cnt+1][k])
	  free_int_matrix(this_dpd->pairorb[cnt+1][k]);
      }
      free(this_dpd->pairorb[cnt]);
      free(this_dpd->pairorb[cnt+1]);
    }
  }
  free(this_dpd->pairidx); free(this_dpd->pairorb);

  for(i=0; i < num_subspaces; i++) {
    free(this_dpd->orbidx2[i]);
    for(j=0; j < nirreps; j++) {
      if(this_dpd->orbspi[i][j]) free(this_dpd->orbs2[i][j]);
    }
    free(this_dpd->orbs2[i]);
  }
  free(this_dpd->orbidx2); free(this_dpd->orbs2);

  for(i=0; i < num_subspaces; i++) {
    free(this_dpd->orbspi[i]);
    free(this_dpd->orbsym[i]);
  }
  free(this_dpd->orbspi); free(this_dpd->orbsym);

  free_int_matrix(this_dpd->pairtot);

  free(this_dpd->numorbs);

  for(i=0; i < num_pairs; i++)
    free(this_dpd->params4[i]);
  free(this_dpd->params4);
  for(i=0; i < num_subspaces; i++)
    free(this_dpd->params2[i]);
  free(this_dpd->params2);

  /*
    printf("memory = %d; memfree = %d\n",
    dpd_main.memory, dpd_main.memfree);
  */

  return 0;
}

}
