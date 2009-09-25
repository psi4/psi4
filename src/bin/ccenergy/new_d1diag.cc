/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

/* Computes a modified D1 diagnostic developed by T.J. Lee, but not yet
 * published.
 * */

double d1diag_t1_rhf(void);

static double new_d1diag_t1_rohf(void)
{
  int h, nirreps, i, j;
  int nclsd, nuocc, nopen;
  double **T1_hp, **T1_hx, **T1_xp, **T1_sq;
  double *E, **C;
  double max_hp=0.0, max_xp=0.0, max_hx=0.0, max;
  dpdfile2 T1_a, T1_b;

  nirreps = moinfo.nirreps;

  dpd_file2_init(&T1_a, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_mat_init(&T1_a);
  dpd_file2_mat_rd(&T1_a);
      
  dpd_file2_init(&T1_b, CC_OEI, 0, 0, 1, "tia");
  dpd_file2_mat_init(&T1_b);
  dpd_file2_mat_rd(&T1_b);

  for(h=0; h < nirreps; h++) {
    nclsd = moinfo.clsdpi[h];
    nuocc = moinfo.uoccpi[h];
    nopen = moinfo.openpi[h];

    if(nclsd && nuocc) {
      T1_hp = block_matrix(nclsd, nuocc);
      for (i=0; i < nclsd; i++)
	for (j=0; j < nuocc; j++)
	  T1_hp[i][j] = (T1_a.matrix[h][i][j] + T1_b.matrix[h][i][j])/2.;

      T1_sq = block_matrix(nclsd, nclsd);
      C_DGEMM('n','t',nclsd,nclsd,nuocc,1.0,&(T1_hp[0][0]),nuocc, 
	     &(T1_hp[0][0]),nuocc,0.0,&(T1_sq[0][0]),nclsd);

      E = init_array(nclsd);
      C = block_matrix(nclsd, nclsd);
      sq_rsp(nclsd, nclsd, T1_sq, E, 0, C, 1e-12);
      for(i=0; i < nclsd; i++) if(E[i] > max_hp) max_hp = E[i];
      free(E);
      free_block(C);
      free_block(T1_sq);
      free_block(T1_hp);
    }

    if(nclsd && nopen) {
      T1_hx = block_matrix(nclsd, nopen);
      for (i=0; i < nclsd; i++)
	for (j=0; j < nopen; j++)
	  T1_hx[i][j] = T1_b.matrix[h][i][nuocc+j]/sqrt(2.);

      T1_sq = block_matrix(nclsd, nclsd);
      C_DGEMM('n','t',nclsd,nclsd,nopen,1.0,&(T1_hx[0][0]),nopen, 
	     &(T1_hx[0][0]),nopen,0.0,&(T1_sq[0][0]),nclsd);

      E = init_array(nclsd);
      C = block_matrix(nclsd, nclsd);
      sq_rsp(nclsd, nclsd, T1_sq, E, 0, C, 1e-12);
      for(i=0; i < nclsd; i++) if(E[i] > max_hx) max_hx = E[i];
      free(E);
      free_block(C);
      free_block(T1_sq);
      free_block(T1_hx);
    }

    if(nopen && nuocc) {
      T1_xp = block_matrix(nopen, nuocc);
      for (i=0; i < nopen; i++)
	for (j=0; j < nuocc; j++)
	  T1_xp[i][j] = T1_a.matrix[h][nclsd+i][j]/sqrt(2.);

      T1_sq = block_matrix(nopen, nopen);
      C_DGEMM('n','t',nopen,nopen,nuocc,1.0,&(T1_xp[0][0]),nuocc, 
	     &(T1_xp[0][0]),nuocc,0.0,&(T1_sq[0][0]),nopen);

      E = init_array(nopen);
      C = block_matrix(nopen, nopen);
      sq_rsp(nopen, nopen, T1_sq, E, 0, C, 1e-12);
      for(i=0; i < nopen; i++) if(E[i] > max_xp) max_xp = E[i];
      free(E);
      free_block(C);
      free_block(T1_sq);
      free_block(T1_xp);
    }
  }

  dpd_file2_mat_close(&T1_a);
  dpd_file2_close(&T1_a);

  dpd_file2_mat_close(&T1_b);
  dpd_file2_close(&T1_b);

  max_hp = sqrt(max_hp);
  max_hx = sqrt(max_hx);
  max_xp = sqrt(max_xp);

  /*
  fprintf(outfile, "ND1: hp=%8.6f hx=%8.6f xp=%8.6f\n", max_hp, max_hx, max_xp);
  */

  max = max_hp;
  if (max_hx > max) max = max_hx;
  if (max_xp > max) max = max_xp;

  return max;
}

double new_d1diag(void)
{
  double norm = 0.0;

  if(params.ref == 0) { /** RHF **/
    norm = d1diag_t1_rhf();
  }
  else if (params.ref == 1) { /** ROHF **/
    norm = new_d1diag_t1_rohf();
  }
  return norm;
}
}} // namespace psi::ccenergy
