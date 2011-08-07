
/*
 *  denom.cc
 *  
 *
 *  Created by M.Saitow on 10/08/01.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace adc {

    void denom(int irrep, double omega)
    {
      char lbl[32];
      int Gij, Gab;
      int ij, ab, i, j, a, b, I, J, A, B, isym, jsym, asym, bsym;
      int nirreps;
      int *occpi, *virpi, *occ_off, *vir_off;
      double fii, fjj, faa, fbb;
      dpdfile2 fij, fab;
      dpdbuf4 D;
      
      nirreps = moinfo.nirreps;
      
      occpi = moinfo.occpi;
      virpi = moinfo.virpi;
      occ_off = moinfo.occ_off;
      vir_off = moinfo.vir_off;
      
      dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fij");
      dpd_file2_mat_init(&fij);
      dpd_file2_mat_rd(&fij);
      
      dpd_file2_init(&fab, CC_OEI, 0, 1, 1, "fab");
      dpd_file2_mat_init(&fab);
      dpd_file2_mat_rd(&fab);
      
      sprintf(lbl, "Dijab[%d]", irrep);
      dpd_buf4_init(&D, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
      for(Gij = 0; Gij < nirreps; Gij++) {
	dpd_buf4_mat_irrep_init(&D, Gij);
	
	for(ij = 0; ij < D.params->rowtot[Gij]; ij++) {
	  i = D.params->roworb[Gij][ij][0];
	  j = D.params->roworb[Gij][ij][1];
	  isym = D.params->psym[i];
	  jsym = D.params->qsym[j];
	  
	  I = i - occ_off[isym];
	  J = j - occ_off[jsym];
	  
	  fii = fij.matrix[isym][I][I];
	  fjj = fij.matrix[jsym][J][J];
	  
	  for(ab = 0; ab < D.params->coltot[Gij^irrep]; ab++) {
	    a = D.params->colorb[Gij^irrep][ab][0];
	    b = D.params->colorb[Gij^irrep][ab][1];
	    asym = D.params->rsym[a];
	    bsym = D.params->ssym[b];
	    
	    A = a - vir_off[asym];
	    B = b - vir_off[bsym];
	    
	    faa = fab.matrix[asym][A][A];
	    fbb = fab.matrix[bsym][B][B];
	    
	    D.matrix[Gij][ij][ab] = 1.0 / (omega+fii+fjj-faa-fbb);
	  }
	}
	dpd_buf4_mat_irrep_wrt(&D, Gij);
	dpd_buf4_mat_irrep_close(&D, Gij);
      }
      dpd_buf4_close(&D);
      dpd_file2_mat_close(&fij);
      dpd_file2_mat_close(&fab);
      dpd_file2_close(&fij);
      dpd_file2_close(&fab);
    }

}}

