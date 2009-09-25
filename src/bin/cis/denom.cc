/*! \file
    \ingroup CIS
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cis {

void denom(int irrep, double root)
{
  int Gij, Gab;
  int ij, ab, i, j, a, b, I, J, A, B, isym, jsym, asym, bsym;
  int nirreps;
  int *occpi, *virtpi, *occ_off, *vir_off;
  int *aoccpi, *avirtpi, *aocc_off, *avir_off;
  int *boccpi, *bvirtpi, *bocc_off, *bvir_off;
  double fii, fjj, faa, fbb;
  dpdfile2 fIJ, fij, fAB, fab;
  dpdbuf4 D;
  char lbl[32];

  nirreps = moinfo.nirreps;

  if(params.ref == 0) { /** RHF **/

    occpi = moinfo.occpi;
    virtpi = moinfo.virtpi;
    occ_off = moinfo.occ_off;
    vir_off = moinfo.vir_off;

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

    sprintf(lbl, "dIjAb[%d]", irrep);
    dpd_buf4_init(&D, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
    for(Gij=0; Gij < nirreps; Gij++) {

      dpd_buf4_mat_irrep_init(&D, Gij);

      for(ij=0; ij < D.params->rowtot[Gij]; ij++) {
	i = D.params->roworb[Gij][ij][0];
	j = D.params->roworb[Gij][ij][1];
	isym = D.params->psym[i];
	jsym = D.params->qsym[j];

	I = i - occ_off[isym];
	J = j - occ_off[jsym];

	fii = fIJ.matrix[isym][I][I];
	fjj = fij.matrix[jsym][J][J];
	  
	for(ab=0; ab < D.params->coltot[Gij^irrep]; ab++) {
	  a = D.params->colorb[Gij^irrep][ab][0];
	  b = D.params->colorb[Gij^irrep][ab][1];
	  asym = D.params->rsym[a];
	  bsym = D.params->ssym[b];

	  A = a - vir_off[asym];
	  B = b - vir_off[bsym];

	  faa = fAB.matrix[asym][A][A];
	  fbb = fab.matrix[bsym][B][B];

	  D.matrix[Gij][ij][ab] = 1.0/(fii + fjj - faa - fbb + root);
      }
    }

    dpd_buf4_mat_irrep_wrt(&D, Gij);
    dpd_buf4_mat_irrep_close(&D, Gij);

  }

  dpd_buf4_close(&D);

}
  else if(params.ref == 2) { /** UHF **/

    aoccpi = moinfo.aoccpi;
    boccpi = moinfo.boccpi;
    avirtpi = moinfo.avirtpi;
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

    sprintf(lbl, "dIJAB[%d]", irrep);
    dpd_buf4_init(&D, CC_MISC, irrep, 1, 6, 1, 6, 0, lbl);
    for(Gij=0; Gij < nirreps; Gij++) {
      dpd_buf4_mat_irrep_init(&D, Gij);
      for(ij=0; ij < D.params->rowtot[Gij]; ij++) {
	i = D.params->roworb[Gij][ij][0];
	j = D.params->roworb[Gij][ij][1];
	isym = D.params->psym[i];
	jsym = D.params->qsym[j];
	I = i - aocc_off[isym];
	J = j - aocc_off[jsym];
	fii = fIJ.matrix[isym][I][I];
	fjj = fIJ.matrix[jsym][J][J];

	for(ab=0; ab < D.params->coltot[Gij^irrep]; ab++) {
	  a = D.params->colorb[Gij^irrep][ab][0];
	  b = D.params->colorb[Gij^irrep][ab][1];
	  asym = D.params->rsym[a];
	  bsym = D.params->ssym[b];
	  A = a - avir_off[asym];
	  B = b - avir_off[bsym];
	  faa = fAB.matrix[asym][A][A];
	  fbb = fAB.matrix[bsym][B][B];

	  D.matrix[Gij][ij][ab] = 1.0/(fii + fjj - faa - fbb + root);
	}
      }
      dpd_buf4_mat_irrep_wrt(&D, Gij);
      dpd_buf4_mat_irrep_close(&D, Gij);
    }
    dpd_buf4_close(&D);

    sprintf(lbl, "dijab[%d]", irrep);
    dpd_buf4_init(&D, CC_MISC, irrep, 11, 16, 11, 16, 0, lbl);
    for(Gij=0; Gij < nirreps; Gij++) {
      dpd_buf4_mat_irrep_init(&D, Gij);
      for(ij=0; ij < D.params->rowtot[Gij]; ij++) {
	i = D.params->roworb[Gij][ij][0];
	j = D.params->roworb[Gij][ij][1];
	isym = D.params->psym[i];
	jsym = D.params->qsym[j];
	I = i - bocc_off[isym];
	J = j - bocc_off[jsym];
	fii = fij.matrix[isym][I][I];
	fjj = fij.matrix[jsym][J][J];

	for(ab=0; ab < D.params->coltot[Gij^irrep]; ab++) {
	  a = D.params->colorb[Gij^irrep][ab][0];
	  b = D.params->colorb[Gij^irrep][ab][1];
	  asym = D.params->rsym[a];
	  bsym = D.params->ssym[b];
	  A = a - bvir_off[asym];
	  B = b - bvir_off[bsym];
	  faa = fab.matrix[asym][A][A];
	  fbb = fab.matrix[bsym][B][B];

	  D.matrix[Gij][ij][ab] = 1.0/(fii + fjj - faa - fbb + root);
	}
      }
      dpd_buf4_mat_irrep_wrt(&D, Gij);
      dpd_buf4_mat_irrep_close(&D, Gij);
    }
    dpd_buf4_close(&D);

    sprintf(lbl, "dIjAb[%d]", irrep);
    dpd_buf4_init(&D, CC_MISC, irrep, 22, 28, 22, 28, 0, lbl);
    for(Gij=0; Gij < nirreps; Gij++) {
      dpd_buf4_mat_irrep_init(&D, Gij);
      for(ij=0; ij < D.params->rowtot[Gij]; ij++) {
	i = D.params->roworb[Gij][ij][0];
	j = D.params->roworb[Gij][ij][1];
	isym = D.params->psym[i];
	jsym = D.params->qsym[j];
	I = i - aocc_off[isym];
	J = j - bocc_off[jsym];
	fii = fIJ.matrix[isym][I][I];
	fjj = fij.matrix[jsym][J][J];

	for(ab=0; ab < D.params->coltot[Gij^irrep]; ab++) {
	  a = D.params->colorb[Gij^irrep][ab][0];
	  b = D.params->colorb[Gij^irrep][ab][1];
	  asym = D.params->rsym[a];
	  bsym = D.params->ssym[b];
	  A = a - avir_off[asym];
	  B = b - bvir_off[bsym];
	  faa = fAB.matrix[asym][A][A];
	  fbb = fab.matrix[bsym][B][B];

	  D.matrix[Gij][ij][ab] = 1.0/(fii + fjj - faa - fbb + root);
	}
      }
      dpd_buf4_mat_irrep_wrt(&D, Gij);
      dpd_buf4_mat_irrep_close(&D, Gij);
    }
    dpd_buf4_close(&D);

    dpd_file2_mat_close(&fIJ);
    dpd_file2_mat_close(&fij);
    dpd_file2_mat_close(&fAB);
    dpd_file2_mat_close(&fab);
    dpd_file2_close(&fIJ);
    dpd_file2_close(&fij);
    dpd_file2_close(&fAB);
    dpd_file2_close(&fab);
  }
}

}} // namespace psi::cis
