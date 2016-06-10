/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

void denom_rhf(void);
void denom_uhf(void);

void denom(void)
{
  if(params.ref == 2) denom_uhf();
  else denom_rhf();
}

void denom_uhf(void)
{
  int h, nirreps;
  int i, j, a, b, ij, ab, I, J, A, B, isym, jsym, asym, bsym;
  int *aoccpi, *avirtpi, *boccpi, *bvirtpi;
  int *aocc_off, *bocc_off, *avir_off, *bvir_off;
  double fii, fjj, faa, fbb;
  dpdfile2 dIA, fIJ, fij, fAB, fab;
  dpdfile4 dIJAB;

  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi;
  avirtpi = moinfo.avirtpi;
  boccpi = moinfo.boccpi;
  bvirtpi = moinfo.bvirtpi;
  aocc_off = moinfo.aocc_off;
  bocc_off = moinfo.bocc_off;
  avir_off = moinfo.avir_off;
  bvir_off = moinfo.bvir_off;

  global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_mat_init(&fIJ);
  global_dpd_->file2_mat_rd(&fIJ);
  
  global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
  global_dpd_->file2_mat_init(&fij);
  global_dpd_->file2_mat_rd(&fij);
  
  global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_mat_init(&fAB);
  global_dpd_->file2_mat_rd(&fAB);
  
  global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");
  global_dpd_->file2_mat_init(&fab);
  global_dpd_->file2_mat_rd(&fab);

  global_dpd_->file2_init(&dIA, PSIF_CC_OEI, 0, 0, 1, "dIA");
  global_dpd_->file2_mat_init(&dIA);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < aoccpi[h]; i++) {
	  fii = fIJ.matrix[h][i][i];
	  for(a=0; a < avirtpi[h]; a++) {
	      faa = fAB.matrix[h][a][a];
	      dIA.matrix[h][i][a] = 1.0/(fii - faa);
	    }
	}
    }
  global_dpd_->file2_mat_wrt(&dIA);
  global_dpd_->file2_mat_close(&dIA);
  global_dpd_->file2_close(&dIA);

  global_dpd_->file2_init(&dIA, PSIF_CC_OEI, 0, 2, 3, "dia");
  global_dpd_->file2_mat_init(&dIA);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < boccpi[h]; i++) {
	  fii = fij.matrix[h][i][i];
	  for(a=0; a < bvirtpi[h]; a++) {
	      faa = fab.matrix[h][a][a];
	      dIA.matrix[h][i][a] = 1.0/(fii - faa);
	    }
	}
    }
  global_dpd_->file2_mat_wrt(&dIA);
  global_dpd_->file2_mat_close(&dIA);
  global_dpd_->file2_close(&dIA);

  global_dpd_->file4_init(&dIJAB, PSIF_CC_DENOM, 0, 1, 6, "dIJAB");
  for(h=0; h < nirreps; h++) {
      global_dpd_->file4_mat_irrep_init(&dIJAB, h);
      for(ij=0; ij < dIJAB.params->rowtot[h]; ij++) {
	  i = dIJAB.params->roworb[h][ij][0];
	  j = dIJAB.params->roworb[h][ij][1];
	  isym = dIJAB.params->psym[i];
	  jsym = dIJAB.params->qsym[j];
	  I = i - aocc_off[isym];
	  J = j - aocc_off[jsym];
	  fii = fIJ.matrix[isym][I][I];
	  fjj = fIJ.matrix[jsym][J][J];

	  for(ab=0; ab < dIJAB.params->coltot[h]; ab++) {
	      a = dIJAB.params->colorb[h][ab][0];
	      b = dIJAB.params->colorb[h][ab][1];
	      asym = dIJAB.params->rsym[a];
	      bsym = dIJAB.params->ssym[b];
	      A = a - avir_off[asym];
	      B = b - avir_off[bsym];
	      faa = fAB.matrix[asym][A][A];
	      fbb = fAB.matrix[bsym][B][B];

	      dIJAB.matrix[h][ij][ab] = 1.0/(fii + fjj - faa - fbb);
	    }
	}
      global_dpd_->file4_mat_irrep_wrt(&dIJAB, h);
      global_dpd_->file4_mat_irrep_close(&dIJAB, h);
    }
  global_dpd_->file4_close(&dIJAB);

  global_dpd_->file4_init(&dIJAB, PSIF_CC_DENOM, 0, 11, 16, "dijab");
  for(h=0; h < nirreps; h++) {
      global_dpd_->file4_mat_irrep_init(&dIJAB, h);
      for(ij=0; ij < dIJAB.params->rowtot[h]; ij++) {
	  i = dIJAB.params->roworb[h][ij][0];
	  j = dIJAB.params->roworb[h][ij][1];
	  isym = dIJAB.params->psym[i];
	  jsym = dIJAB.params->qsym[j];
	  I = i - bocc_off[isym];
	  J = j - bocc_off[jsym];
	  fii = fij.matrix[isym][I][I];
	  fjj = fij.matrix[jsym][J][J];

	  for(ab=0; ab < dIJAB.params->coltot[h]; ab++) {
	      a = dIJAB.params->colorb[h][ab][0];
	      b = dIJAB.params->colorb[h][ab][1];
	      asym = dIJAB.params->rsym[a];
	      bsym = dIJAB.params->ssym[b];
	      A = a - bvir_off[asym];
	      B = b - bvir_off[bsym];
	      faa = fab.matrix[asym][A][A];
	      fbb = fab.matrix[bsym][B][B];

	      dIJAB.matrix[h][ij][ab] = 1.0/(fii + fjj - faa - fbb);
	    }
	}
      global_dpd_->file4_mat_irrep_wrt(&dIJAB, h);
      global_dpd_->file4_mat_irrep_close(&dIJAB, h);
    }
  global_dpd_->file4_close(&dIJAB);

  global_dpd_->file4_init(&dIJAB, PSIF_CC_DENOM, 0, 22, 28, "dIjAb");
  for(h=0; h < nirreps; h++) {
      global_dpd_->file4_mat_irrep_init(&dIJAB, h);
      for(ij=0; ij < dIJAB.params->rowtot[h]; ij++) {
	  i = dIJAB.params->roworb[h][ij][0];
	  j = dIJAB.params->roworb[h][ij][1];
	  isym = dIJAB.params->psym[i];
	  jsym = dIJAB.params->qsym[j];
	  I = i - aocc_off[isym];
	  J = j - bocc_off[jsym];
	  fii = fIJ.matrix[isym][I][I];
	  fjj = fij.matrix[jsym][J][J];

	  for(ab=0; ab < dIJAB.params->coltot[h]; ab++) {
	      a = dIJAB.params->colorb[h][ab][0];
	      b = dIJAB.params->colorb[h][ab][1];
	      asym = dIJAB.params->rsym[a];
	      bsym = dIJAB.params->ssym[b];
	      A = a - avir_off[asym];
	      B = b - bvir_off[bsym];
	      faa = fAB.matrix[asym][A][A];
	      fbb = fab.matrix[bsym][B][B];

	      dIJAB.matrix[h][ij][ab] = 1.0/(fii + fjj - faa - fbb);
	    }
	}
      global_dpd_->file4_mat_irrep_wrt(&dIJAB, h);
      global_dpd_->file4_mat_irrep_close(&dIJAB, h);
    }
  global_dpd_->file4_close(&dIJAB);

  global_dpd_->file2_mat_close(&fIJ);
  global_dpd_->file2_mat_close(&fij);
  global_dpd_->file2_mat_close(&fAB);
  global_dpd_->file2_mat_close(&fab);
  global_dpd_->file2_close(&fIJ);
  global_dpd_->file2_close(&fij);
  global_dpd_->file2_close(&fAB);
  global_dpd_->file2_close(&fab);
}

void denom_rhf(void)
{
  int nirreps;
  int h, i, j, a, b, ij, ab;
  int I, J, A, B;
  int isym, jsym, asym, bsym;
  int *occpi, *virtpi;
  int *occ_off, *vir_off;
  int *openpi;
  double fii, fjj, faa, fbb;
  dpdfile2 fIJ, fij, fAB, fab;
  dpdfile2 dIA, dia;
  dpdfile4 dIJAB, dijab, dIjAb;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  openpi = moinfo.openpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;

  /* Grab Fock matrices from disk */
  global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_mat_init(&fIJ);
  global_dpd_->file2_mat_rd(&fIJ);
  
  global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fij");
  global_dpd_->file2_mat_init(&fij);
  global_dpd_->file2_mat_rd(&fij);
  
  global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_mat_init(&fAB);
  global_dpd_->file2_mat_rd(&fAB);
  
  global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fab");
  global_dpd_->file2_mat_init(&fab);
  global_dpd_->file2_mat_rd(&fab);

  /* Alpha one-electron denominator */
  global_dpd_->file2_init(&dIA, PSIF_CC_OEI, 0, 0, 1, "dIA");
  global_dpd_->file2_mat_init(&dIA);

  for(h=0; h < nirreps; h++) {
      
      for(i=0; i < occpi[h]; i++) {
	  fii = fIJ.matrix[h][i][i];

	  for(a=0; a < (virtpi[h] - openpi[h]); a++) {
	      faa = fAB.matrix[h][a][a];

	      dIA.matrix[h][i][a] = 1.0/(fii - faa);
	    }
	}
      
    }

  global_dpd_->file2_mat_wrt(&dIA);
  global_dpd_->file2_mat_close(&dIA);
  global_dpd_->file2_close(&dIA);

  /* Beta one-electron denominator */
  global_dpd_->file2_init(&dia, PSIF_CC_OEI, 0, 0, 1, "dia");
  global_dpd_->file2_mat_init(&dia);

  for(h=0; h < nirreps; h++) {
      
      for(i=0; i < (occpi[h] - openpi[h]); i++) {
	  fii = fij.matrix[h][i][i];
 
	  for(a=0; a < virtpi[h]; a++) {
	      faa = fab.matrix[h][a][a];

	      dia.matrix[h][i][a] = 1.0/(fii - faa);
	    }
	}
      
    }
  
  global_dpd_->file2_mat_wrt(&dia);
  global_dpd_->file2_mat_close(&dia);
  global_dpd_->file2_close(&dia);

  /* Alpha-alpha two-electron denominator */
  global_dpd_->file4_init(&dIJAB, PSIF_CC_DENOM, 0, 1, 6, "dIJAB");

  for(h=0; h < nirreps; h++) {

      global_dpd_->file4_mat_irrep_init(&dIJAB, h);

      /* Loop over the rows */
      for(ij=0; ij < dIJAB.params->rowtot[h]; ij++) {
	  i = dIJAB.params->roworb[h][ij][0];
	  j = dIJAB.params->roworb[h][ij][1];
	  isym = dIJAB.params->psym[i];
	  jsym = dIJAB.params->qsym[j];

	  /* Convert to relative orbital index */
	  I = i - occ_off[isym];
	  J = j - occ_off[jsym];

	  fii = fIJ.matrix[isym][I][I];
	  fjj = fIJ.matrix[jsym][J][J];
	  
	  /* Loop over the columns */
	  for(ab=0; ab < dIJAB.params->coltot[h]; ab++) {
	      a = dIJAB.params->colorb[h][ab][0];
	      b = dIJAB.params->colorb[h][ab][1];
	      asym = dIJAB.params->rsym[a];
	      bsym = dIJAB.params->ssym[b];

	      /* Convert to relative orbital index */
	      A = a - vir_off[asym];
	      B = b - vir_off[bsym];

	      faa = fAB.matrix[asym][A][A];
	      fbb = fAB.matrix[bsym][B][B];

	      dIJAB.matrix[h][ij][ab] =
                ((A >= (virtpi[asym] - openpi[asym])) ||
		 (B >= (virtpi[bsym] - openpi[bsym])) ?
		 0.0 : 1.0/(fii + fjj - faa - fbb));
	    }
	}

      global_dpd_->file4_mat_irrep_wrt(&dIJAB, h);
      global_dpd_->file4_mat_irrep_close(&dIJAB, h);

    }

  global_dpd_->file4_close(&dIJAB);

  /* Beta-beta two-electron denominator */
  global_dpd_->file4_init(&dijab, PSIF_CC_DENOM, 0, 1, 6, "dijab");

  for(h=0; h < nirreps; h++) {

      global_dpd_->file4_mat_irrep_init(&dijab, h);

      /* Loop over the rows */
      for(ij=0; ij < dijab.params->rowtot[h]; ij++) {
	  i = dijab.params->roworb[h][ij][0];
	  j = dijab.params->roworb[h][ij][1];
	  isym = dijab.params->psym[i];
	  jsym = dijab.params->qsym[j];

	  /* Convert to relative orbital index */
	  I = i - occ_off[isym];
	  J = j - occ_off[jsym];

	  fii = fij.matrix[isym][I][I];
	  fjj = fij.matrix[jsym][J][J];
	  
	  /* Loop over the columns */
	  for(ab=0; ab < dijab.params->coltot[h]; ab++) {
	      a = dijab.params->colorb[h][ab][0];
	      b = dijab.params->colorb[h][ab][1];
	      asym = dijab.params->rsym[a];
	      bsym = dijab.params->ssym[b];

	      /* Convert to relative orbital index */
	      A = a - vir_off[asym];
	      B = b - vir_off[bsym];

	      faa = fab.matrix[asym][A][A];
	      fbb = fab.matrix[bsym][B][B];

	      dijab.matrix[h][ij][ab] =
                ((I >= (occpi[isym] - openpi[isym])) ||
		 (J >= (occpi[jsym] - openpi[jsym])) ?
		 0.0 : 1.0/(fii + fjj - faa - fbb));
	    }
	}

      global_dpd_->file4_mat_irrep_wrt(&dijab, h);
      global_dpd_->file4_mat_irrep_close(&dijab, h);

    }

  global_dpd_->file4_close(&dijab);

  /* Alpha-beta two-electron denominator */
  global_dpd_->file4_init(&dIjAb, PSIF_CC_DENOM, 0, 0, 5, "dIjAb");

  for(h=0; h < nirreps; h++) {

      global_dpd_->file4_mat_irrep_init(&dIjAb, h);

      /* Loop over the rows */
      for(ij=0; ij < dIjAb.params->rowtot[h]; ij++) {
	  i = dIjAb.params->roworb[h][ij][0];
	  j = dIjAb.params->roworb[h][ij][1];
	  isym = dIjAb.params->psym[i];
	  jsym = dIjAb.params->qsym[j];

	  /* Convert to relative orbital index */
	  I = i - occ_off[isym];
	  J = j - occ_off[jsym];

	  fii = fIJ.matrix[isym][I][I];
	  fjj = fij.matrix[jsym][J][J];
	  
	  /* Loop over the columns */
	  for(ab=0; ab < dIjAb.params->coltot[h]; ab++) {
	      a = dIjAb.params->colorb[h][ab][0];
	      b = dIjAb.params->colorb[h][ab][1];
	      asym = dIjAb.params->rsym[a];
	      bsym = dIjAb.params->ssym[b];

	      /* Convert to relative orbital index */
	      A = a - vir_off[asym];
	      B = b - vir_off[bsym];

	      faa = fAB.matrix[asym][A][A];
	      fbb = fab.matrix[bsym][B][B];

	      dIjAb.matrix[h][ij][ab] =
                ((A >= (virtpi[asym] - openpi[asym])) ||
		 (J >= (occpi[jsym] - openpi[jsym])) ?
		 0.0 : 1.0/(fii + fjj - faa - fbb));
	    }
	}

      global_dpd_->file4_mat_irrep_wrt(&dIjAb, h);
      global_dpd_->file4_mat_irrep_close(&dIjAb, h);

    }

  global_dpd_->file4_close(&dIjAb);

  global_dpd_->file2_mat_close(&fIJ);
  global_dpd_->file2_mat_close(&fij);
  global_dpd_->file2_mat_close(&fAB);
  global_dpd_->file2_mat_close(&fab);
  global_dpd_->file2_close(&fIJ);
  global_dpd_->file2_close(&fij);
  global_dpd_->file2_close(&fAB);
  global_dpd_->file2_close(&fab);

}

}} // namespace psi::ccsort