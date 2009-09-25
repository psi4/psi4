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

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_rd(&fIJ);
  
  dpd_file2_init(&fij, CC_OEI, 0, 2, 2, "fij");
  dpd_file2_mat_init(&fij);
  dpd_file2_mat_rd(&fij);
  
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_rd(&fAB);
  
  dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab");
  dpd_file2_mat_init(&fab);
  dpd_file2_mat_rd(&fab);

  dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
  dpd_file2_mat_init(&dIA);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < aoccpi[h]; i++) {
	  fii = fIJ.matrix[h][i][i];
	  for(a=0; a < avirtpi[h]; a++) {
	      faa = fAB.matrix[h][a][a];
	      dIA.matrix[h][i][a] = 1.0/(fii - faa);
	    }
	}
    }
  dpd_file2_mat_wrt(&dIA);
  dpd_file2_mat_close(&dIA);
  dpd_file2_close(&dIA);

  dpd_file2_init(&dIA, CC_OEI, 0, 2, 3, "dia");
  dpd_file2_mat_init(&dIA);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < boccpi[h]; i++) {
	  fii = fij.matrix[h][i][i];
	  for(a=0; a < bvirtpi[h]; a++) {
	      faa = fab.matrix[h][a][a];
	      dIA.matrix[h][i][a] = 1.0/(fii - faa);
	    }
	}
    }
  dpd_file2_mat_wrt(&dIA);
  dpd_file2_mat_close(&dIA);
  dpd_file2_close(&dIA);

  dpd_file4_init(&dIJAB, CC_DENOM, 0, 1, 6, "dIJAB");
  for(h=0; h < nirreps; h++) {
      dpd_file4_mat_irrep_init(&dIJAB, h);
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
      dpd_file4_mat_irrep_wrt(&dIJAB, h);
      dpd_file4_mat_irrep_close(&dIJAB, h);
    }
  dpd_file4_close(&dIJAB);

  dpd_file4_init(&dIJAB, CC_DENOM, 0, 11, 16, "dijab");
  for(h=0; h < nirreps; h++) {
      dpd_file4_mat_irrep_init(&dIJAB, h);
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
      dpd_file4_mat_irrep_wrt(&dIJAB, h);
      dpd_file4_mat_irrep_close(&dIJAB, h);
    }
  dpd_file4_close(&dIJAB);

  dpd_file4_init(&dIJAB, CC_DENOM, 0, 22, 28, "dIjAb");
  for(h=0; h < nirreps; h++) {
      dpd_file4_mat_irrep_init(&dIJAB, h);
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
      dpd_file4_mat_irrep_wrt(&dIJAB, h);
      dpd_file4_mat_irrep_close(&dIJAB, h);
    }
  dpd_file4_close(&dIJAB);

  dpd_file2_mat_close(&fIJ);
  dpd_file2_mat_close(&fij);
  dpd_file2_mat_close(&fAB);
  dpd_file2_mat_close(&fab);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fij);
  dpd_file2_close(&fAB);
  dpd_file2_close(&fab);
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
  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_rd(&fIJ);
  
  dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fij");
  dpd_file2_mat_init(&fij);
  dpd_file2_mat_rd(&fij);
  
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_rd(&fAB);
  
  dpd_file2_init(&fab, CC_OEI, 0, 1, 1, "fab");
  dpd_file2_mat_init(&fab);
  dpd_file2_mat_rd(&fab);

  /* Alpha one-electron denominator */
  dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
  dpd_file2_mat_init(&dIA);

  for(h=0; h < nirreps; h++) {
      
      for(i=0; i < occpi[h]; i++) {
	  fii = fIJ.matrix[h][i][i];

	  for(a=0; a < (virtpi[h] - openpi[h]); a++) {
	      faa = fAB.matrix[h][a][a];

	      dIA.matrix[h][i][a] = 1.0/(fii - faa);
	    }
	}
      
    }

  dpd_file2_mat_wrt(&dIA);
  dpd_file2_mat_close(&dIA);
  dpd_file2_close(&dIA);

  /* Beta one-electron denominator */
  dpd_file2_init(&dia, CC_OEI, 0, 0, 1, "dia");
  dpd_file2_mat_init(&dia);

  for(h=0; h < nirreps; h++) {
      
      for(i=0; i < (occpi[h] - openpi[h]); i++) {
	  fii = fij.matrix[h][i][i];
 
	  for(a=0; a < virtpi[h]; a++) {
	      faa = fab.matrix[h][a][a];

	      dia.matrix[h][i][a] = 1.0/(fii - faa);
	    }
	}
      
    }
  
  dpd_file2_mat_wrt(&dia);
  dpd_file2_mat_close(&dia);
  dpd_file2_close(&dia);

  /* Alpha-alpha two-electron denominator */
  dpd_file4_init(&dIJAB, CC_DENOM, 0, 1, 6, "dIJAB");

  for(h=0; h < nirreps; h++) {

      dpd_file4_mat_irrep_init(&dIJAB, h);

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

      dpd_file4_mat_irrep_wrt(&dIJAB, h);
      dpd_file4_mat_irrep_close(&dIJAB, h);

    }

  dpd_file4_close(&dIJAB);

  /* Beta-beta two-electron denominator */
  dpd_file4_init(&dijab, CC_DENOM, 0, 1, 6, "dijab");

  for(h=0; h < nirreps; h++) {

      dpd_file4_mat_irrep_init(&dijab, h);

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

      dpd_file4_mat_irrep_wrt(&dijab, h);
      dpd_file4_mat_irrep_close(&dijab, h);

    }

  dpd_file4_close(&dijab);

  /* Alpha-beta two-electron denominator */
  dpd_file4_init(&dIjAb, CC_DENOM, 0, 0, 5, "dIjAb");

  for(h=0; h < nirreps; h++) {

      dpd_file4_mat_irrep_init(&dIjAb, h);

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

      dpd_file4_mat_irrep_wrt(&dIjAb, h);
      dpd_file4_mat_irrep_close(&dIjAb, h);

    }

  dpd_file4_close(&dIjAb);

  dpd_file2_mat_close(&fIJ);
  dpd_file2_mat_close(&fij);
  dpd_file2_mat_close(&fAB);
  dpd_file2_mat_close(&fab);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fij);
  dpd_file2_close(&fAB);
  dpd_file2_close(&fab);

}

}} // namespace psi::ccsort
