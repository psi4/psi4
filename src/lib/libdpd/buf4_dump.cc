/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libiwl/iwl.h>
#include "dpd.h"

namespace psi {
	
int dpd_buf4_dump(dpdbuf4 *DPDBuf, struct iwlbuf *IWLBuf,
 		  int *prel, int *qrel, int *rrel, int *srel, 
		  int bk_pack, int swap23)
{
  int h, row, col, p, q, r, s, P, Q, R, S, my_irrep;
  int soccs;
  double value;

  my_irrep = DPDBuf->file.my_irrep;

  for(h=0; h < DPDBuf->params->nirreps; h++) {
    dpd_buf4_mat_irrep_init(DPDBuf, h);
    dpd_buf4_mat_irrep_rd(DPDBuf, h);
    for(row=0; row < DPDBuf->params->rowtot[h]; row++) {
      p = DPDBuf->params->roworb[h][row][0]; P = prel[p];
      q = DPDBuf->params->roworb[h][row][1]; Q = qrel[q];
      if(bk_pack) {
	for(col=0; col <= row; col++) {
	  r = DPDBuf->params->colorb[h^my_irrep][col][0]; R = rrel[r];
	  s = DPDBuf->params->colorb[h^my_irrep][col][1]; S = srel[s];
		  
	  value = DPDBuf->matrix[h][row][col];

	  if(swap23)
	    iwl_buf_wrt_val(IWLBuf, P, R, Q, S, value, 0, 
			    (FILE *) NULL, 0);
	  else
	    iwl_buf_wrt_val(IWLBuf, P, Q, R, S, value, 0, 
			    (FILE *) NULL, 0);
	}
      }
      else {
	for(col=0; col < DPDBuf->params->coltot[h^my_irrep]; col++) {
	  r = DPDBuf->params->colorb[h^my_irrep][col][0]; R = rrel[r];
	  s = DPDBuf->params->colorb[h^my_irrep][col][1]; S = srel[s];

	  value = DPDBuf->matrix[h][row][col];

	  if(swap23)
	    iwl_buf_wrt_val(IWLBuf, P, R, Q, S, value, 0,
			    (FILE *) NULL, 0);
	  else
	    iwl_buf_wrt_val(IWLBuf, P, Q, R, S, value, 0,
			    (FILE *) NULL, 0);
	}
      }
    }
    dpd_buf4_mat_irrep_close(DPDBuf, h);
  }

  return 0;
}

}

