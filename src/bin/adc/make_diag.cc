
/*
 *  make_diag.cc
 *  
 *
 *  Created by M.Saitow on 10/07/28.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libiwl/iwl.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {namespace adc {
	
    void make_diag(void)
    {
      char lbl[32];
      int a, m, ma, asym, ia, im, irrep;
      dpdfile2 D;
      dpdbuf4 A;
      
      dpd_buf4_init(&A, CC_MISC, 0, 11, 11, 11, 11, 0, "dA(AI, BJ)");
      for(irrep = 0;irrep < moinfo.nirreps;irrep++){
	dpd_buf4_mat_irrep_init(&A, irrep);
	dpd_buf4_mat_irrep_rd(&A, irrep);
	sprintf(lbl, "Dia[%d]", irrep);
	dpd_file2_init(&D, CC_MISC, irrep, 0, 1, lbl);
	dpd_file2_mat_init(&D);
	for(ma = 0;ma < A.params->rowtot[irrep];ma++){
	  m = A.params->roworb[irrep][ma][0];
	  a = A.params->roworb[irrep][ma][1];
	  
	  ia = D.params->rowidx[a];
	  im = D.params->colidx[m];
	  asym = D.params->psym[a];
	  D.matrix[asym][ia][im] = A.matrix[irrep][ma][ma];
	}
	dpd_buf4_mat_irrep_close(&A, irrep);
	dpd_file2_mat_wrt(&D);
	dpd_file2_mat_close(&D);
	dpd_file2_close(&D);
      }
      dpd_buf4_close(&A);
    }
    
    void shift_diag(double lambda, int root, int irrep)
    {
      char lbl[32];
      int asym, a, m;
      double denom;
      dpdfile2 D, L;
      
      sprintf(lbl, "Dia[%d]", irrep);
      dpd_file2_init(&D, CC_MISC, irrep, 0, 1, lbl);
      dpd_file2_mat_init(&D);
      dpd_file2_mat_rd(&D);
      
      sprintf(lbl, "Lia^(%d)[%d]", root, irrep);
      dpd_file2_init(&L, CC_MISC, irrep, 0, 1, lbl);
      dpd_file2_mat_init(&L);
      dpd_file2_mat_rd(&L);
      
      for(asym = 0;asym < moinfo.nirreps;asym++){
	for(a = 0;a < D.params->rowtot[asym];a++){
	  for(m = 0;m < D.params->coltot[asym^irrep];m++){
	    denom = lambda - D.matrix[asym][a][m];
	    if(fabs(denom) > 1e-6)
	      L.matrix[asym][a][m] = 1 / denom;
	    else 
	      L.matrix[asym][a][m] = 0.0;
	  }
	}
      }
      dpd_file2_mat_wrt(&D);
      dpd_file2_mat_close(&D);
      dpd_file2_mat_wrt(&L);
      dpd_file2_mat_close(&L);
    }
	
}}

