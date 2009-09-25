/*! \file
    \ingroup CCHBAR
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

void purge(void) {
  dpdfile2 FAE, Fmi, FME, Fme;
  dpdfile4 W;
  int *occpi, *virtpi;
  int h, a, b, e, f, i, j, m, n;
  int    A, B, E, F, I, J, M, N;
  int mn, ei, ma, ef, me, jb, mb, ij, ab;
  int asym, bsym, esym, fsym, isym, jsym, msym, nsym;
  int *occ_off, *vir_off;
  int *occ_sym, *vir_sym;
  int *openpi, nirreps;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  /* Purge FME matrix elements */
  dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_mat_init(&FME);
  dpd_file2_mat_rd(&FME);
  for(h=0; h < nirreps; h++) {
    for(m=0; m<occpi[h]; m++)
      for(e=(virtpi[h]-openpi[h]); e<virtpi[h]; e++)
	FME.matrix[h][m][e] = 0.0;
  }
  dpd_file2_mat_wrt(&FME);
  dpd_file2_mat_close(&FME);
  dpd_file2_close(&FME);

  /* Purge Fme matrix elements */
  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");
  dpd_file2_mat_init(&Fme);
  dpd_file2_mat_rd(&Fme);
  for(h=0; h < nirreps; h++) {
    for(e=0; e<virtpi[h]; e++)
      for(m=(occpi[h]-openpi[h]); m<occpi[h]; m++)
	Fme.matrix[h][m][e] = 0.0;
  }
  dpd_file2_mat_wrt(&Fme);
  dpd_file2_mat_close(&Fme);
  dpd_file2_close(&Fme);

  /* Purge Fmi matrix elements */
  dpd_file2_init(&Fmi, CC_OEI, 0, 0, 0, "Fmi");
  dpd_file2_mat_init(&Fmi);
  dpd_file2_mat_rd(&Fmi);
  for(h=0; h < nirreps; h++) {

    for(i=0; i<occpi[h]; i++)
      for(j=(occpi[h]-openpi[h]); j<occpi[h]; j++)
	Fmi.matrix[h][i][j] = 0.0;

    for(i=(occpi[h]-openpi[h]); i<occpi[h]; i++)
      for(j=0; j<occpi[h]; j++)
	Fmi.matrix[h][i][j] = 0.0;

  }
  dpd_file2_mat_wrt(&Fmi);
  dpd_file2_mat_close(&Fmi);
  dpd_file2_close(&Fmi);

  /* Purge FAE matrix elements */
  dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
  dpd_file2_mat_init(&FAE);
  dpd_file2_mat_rd(&FAE);
  for(h=0; h < nirreps; h++) {

    for(a=0; a<virtpi[h]; a++)
      for(b=(virtpi[h]-openpi[h]); b<virtpi[h]; b++)
	FAE.matrix[h][a][b] = 0.0;

    for(a=(virtpi[h]-openpi[h]); a<virtpi[h]; a++)
      for(b=0; b<virtpi[h]; b++)
	FAE.matrix[h][a][b] = 0.0;

  }
  dpd_file2_mat_wrt(&FAE);
  dpd_file2_mat_close(&FAE);
  dpd_file2_close(&FAE);

  /* Purge Fmit (matrix elements with zero diagonal) */
  dpd_file2_init(&Fmi, CC_OEI, 0, 0, 0, "Fmit");
  dpd_file2_mat_init(&Fmi);
  dpd_file2_mat_rd(&Fmi);
  for(h=0; h < nirreps; h++) {

    for(i=0; i<occpi[h]; i++)
      for(j=(occpi[h]-openpi[h]); j<occpi[h]; j++)
	Fmi.matrix[h][i][j] = 0.0;

    for(i=(occpi[h]-openpi[h]); i<occpi[h]; i++)
      for(j=0; j<occpi[h]; j++)
	Fmi.matrix[h][i][j] = 0.0;

  }
  dpd_file2_mat_wrt(&Fmi);
  dpd_file2_mat_close(&Fmi);
  dpd_file2_close(&Fmi);

  /* Purge FAEt (matrix elements with zero diagonal) */
  dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAEt");
  dpd_file2_mat_init(&FAE);
  dpd_file2_mat_rd(&FAE);
  for(h=0; h < nirreps; h++) {

    for(a=0; a<virtpi[h]; a++)
      for(b=(virtpi[h]-openpi[h]); b<virtpi[h]; b++)
	FAE.matrix[h][a][b] = 0.0;

    for(a=(virtpi[h]-openpi[h]); a<virtpi[h]; a++)
      for(b=0; b<virtpi[h]; b++)
	FAE.matrix[h][a][b] = 0.0;

  }
  dpd_file2_mat_wrt(&FAE);
  dpd_file2_mat_close(&FAE);
  dpd_file2_close(&FAE);

  /* Purge Wmnij matrix elements */
  dpd_file4_init(&W, CC_HBAR, 0, 2, 2,"Wmnij");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mn=0; mn < W.params->rowtot[h]; mn++) {
      m = W.params->roworb[h][mn][0];
      n = W.params->roworb[h][mn][1];
      msym = W.params->psym[m];
      nsym = W.params->qsym[n];
      M = m - occ_off[msym];
      N = n - occ_off[nsym];
      for(ij=0; ij < W.params->coltot[h]; ij++) {
	i = W.params->colorb[h][ij][0];
	j = W.params->colorb[h][ij][1];
	isym = W.params->rsym[i];
	jsym = W.params->ssym[j];
	I = i - occ_off[isym];
	J = j - occ_off[jsym];
	if ((I >= (occpi[isym] - openpi[isym])) ||
	    (J >= (occpi[jsym] - openpi[jsym])) ||
	    (M >= (occpi[msym] - openpi[msym])) ||
	    (N >= (occpi[nsym] - openpi[nsym])) )
	  W.matrix[h][mn][ij] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 0, 0,"WMnIj");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mn=0; mn < W.params->rowtot[h]; mn++) {
      n = W.params->roworb[h][mn][1];
      nsym = W.params->qsym[n];
      N = n - occ_off[nsym];
      for(ij=0; ij < W.params->coltot[h]; ij++) {
	j = W.params->colorb[h][ij][1];
	jsym = W.params->ssym[j];
	J = j - occ_off[jsym];
	if ((J >= (occpi[jsym] - openpi[jsym])) ||
	    (N >= (occpi[nsym] - openpi[nsym])) )
	  W.matrix[h][mn][ij] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  /* Purge Wmbej matrix elements */
  dpd_file4_init(&W, CC_HBAR, 0, 10, 10,"WMBEJ");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(me=0; me < W.params->rowtot[h]; me++) {
      e = W.params->roworb[h][me][1];
      esym = W.params->qsym[e];
      E = e - vir_off[esym];
      for(jb=0; jb < W.params->coltot[h]; jb++) {
	b = W.params->colorb[h][jb][1];
	bsym = W.params->ssym[b];
	B = b - vir_off[bsym];
	if ((E >= (virtpi[esym] - openpi[esym])) ||
	    (B >= (virtpi[bsym] - openpi[bsym])) )
	  W.matrix[h][me][jb] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 10, 10,"Wmbej");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(me=0; me < W.params->rowtot[h]; me++) {
      m = W.params->roworb[h][me][0];
      msym = W.params->psym[m];
      M = m - occ_off[msym];
      for(jb=0; jb < W.params->coltot[h]; jb++) {
	j = W.params->colorb[h][jb][0];
	jsym = W.params->rsym[j];
	J = j - occ_off[jsym];
	if ((M >= (occpi[msym] - openpi[msym])) ||
	    (J >= (occpi[jsym] - openpi[jsym])) )
	  W.matrix[h][me][jb] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 10, 10,"WMbEj");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(me=0; me < W.params->rowtot[h]; me++) {
      e = W.params->roworb[h][me][1];
      esym = W.params->qsym[e];
      E = e - vir_off[esym];
      for(jb=0; jb < W.params->coltot[h]; jb++) {
	j = W.params->colorb[h][jb][0];
	jsym = W.params->rsym[j];
	J = j - occ_off[jsym];
	if ((E >= (virtpi[esym] - openpi[esym])) ||
	    (J >= (occpi[jsym] - openpi[jsym])) )
	  W.matrix[h][me][jb] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 10, 10,"WmBeJ");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(me=0; me < W.params->rowtot[h]; me++) {
      m = W.params->roworb[h][me][0];
      msym = W.params->psym[m];
      M = m - occ_off[msym];
      for(jb=0; jb < W.params->coltot[h]; jb++) {
	b = W.params->colorb[h][jb][1];
	bsym = W.params->ssym[b];
	B = b - vir_off[bsym];
	if ((M >= (occpi[msym] - openpi[msym])) ||
	    (B >= (virtpi[bsym] - openpi[bsym])) )
	  W.matrix[h][me][jb] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 10, 10,"WmBEj");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(me=0; me < W.params->rowtot[h]; me++) {
      m = W.params->roworb[h][me][0];
      e = W.params->roworb[h][me][1];
      msym = W.params->psym[m];
      esym = W.params->qsym[e];
      M = m - occ_off[msym];
      E = e - vir_off[esym];
      for(jb=0; jb < W.params->coltot[h]; jb++) {
	j = W.params->colorb[h][jb][0];
	b = W.params->colorb[h][jb][1];
	jsym = W.params->rsym[j];
	bsym = W.params->ssym[b];
	J = j - occ_off[jsym];
	B = b - vir_off[bsym];
	if ((M >= (occpi[msym] - openpi[msym])) ||
	    (E >= (virtpi[esym] - openpi[esym])) ||
	    (J >= (occpi[jsym] - openpi[jsym])) ||
	    (B >= (virtpi[bsym] - openpi[bsym])) )
	  W.matrix[h][me][jb] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  /* WMbeJ is already OK */

  /* Purge Wamef matrix elements */
  dpd_file4_init(&W, CC_HBAR, 0, 11, 7,"WAMEF");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      a = W.params->roworb[h][ma][0];
      asym = W.params->psym[a];
      A = a - vir_off[asym];
      for(ef=0; ef< W.params->coltot[h]; ef++) {
	e = W.params->colorb[h][ef][0];
	f = W.params->colorb[h][ef][1];
	esym = W.params->rsym[e];
	fsym = W.params->ssym[f];
	E = e - vir_off[esym];
	F = f - vir_off[fsym];
	if ((A >= (virtpi[asym] - openpi[asym])) ||
	    (E >= (virtpi[esym] - openpi[esym])) ||
	    (F >= (virtpi[fsym] - openpi[fsym])) )
	  W.matrix[h][ma][ef] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 11, 7,"Wamef");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      m = W.params->roworb[h][ma][1];
      msym = W.params->qsym[m];
      M = m - occ_off[msym];
      for(ef=0; ef< W.params->coltot[h]; ef++) {
	if (M >=  (occpi[msym] - openpi[msym]))
	  W.matrix[h][ma][ef] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 11, 5,"WAmEf");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      a = W.params->roworb[h][ma][0];
      m = W.params->roworb[h][ma][1];
      asym = W.params->psym[a];
      msym = W.params->qsym[m];
      M = m - occ_off[msym];
      A = a - vir_off[asym];
      for(ef=0; ef< W.params->coltot[h]; ef++) {
	e = W.params->colorb[h][ef][0];
	esym = W.params->rsym[e];
	E = e - vir_off[esym];
	if ((A >= (virtpi[asym] - openpi[asym])) ||
	    (M >=  (occpi[msym] - openpi[msym])) ||
	    (E >= (virtpi[esym] - openpi[esym])) )
	  W.matrix[h][ma][ef] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 11, 5,"WaMeF");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      for(ef=0; ef< W.params->coltot[h]; ef++) {
	f = W.params->colorb[h][ef][1];
	fsym = W.params->ssym[f];
	F = f - vir_off[fsym];
	if (F >= (virtpi[fsym] - openpi[fsym]))
	  W.matrix[h][ma][ef] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);


  /* Purge Wmnie matrix elements */
  /* moving to separate function
  dpd_file4_init(&W, CC_HBAR, 0, 2, 11,"WMNIE");
  for(h=0; h < W.params->nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      for(ei=0; ei<W.params->coltot[h]; ei++) {
	e = W.params->colorb[h][ei][0];
	esym = W.params->rsym[e];
	E = e - vir_off[esym];
	if (E >= (virtpi[esym] - openpi[esym]))
	  W.matrix[h][mn][ei] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 2, 11,"Wmnie");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      m = W.params->roworb[h][mn][0];
      n = W.params->roworb[h][mn][1];
      msym = W.params->psym[m];
      nsym = W.params->qsym[n];
      M = m - occ_off[msym];
      N = n - occ_off[nsym];
      for(ei=0; ei<W.params->coltot[h]; ei++) {
	i = W.params->colorb[h][ei][1];
	isym = W.params->ssym[i];
	I = i - occ_off[isym];
	if ((M >= (occpi[msym] - openpi[msym])) ||
	    (N >= (occpi[nsym] - openpi[nsym])) ||
	    (I >= (occpi[isym] - openpi[isym])) )
	  W.matrix[h][mn][ei] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 0, 11,"WMnIe");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      n = W.params->roworb[h][mn][1];
      nsym = W.params->qsym[n];
      N = n - occ_off[nsym];
      for(ei=0; ei<W.params->coltot[h]; ei++) {
	if (N >= (occpi[nsym] - openpi[nsym]))
	  W.matrix[h][mn][ei] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }

  dpd_file4_init(&W, CC_HBAR, 0, 0, 11,"WmNiE");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      m = W.params->roworb[h][mn][0];
      msym = W.params->psym[m];
      M = m - occ_off[msym];
      for(ei=0; ei<W.params->coltot[h]; ei++) {
	e = W.params->colorb[h][ei][0];
	i = W.params->colorb[h][ei][1];
	esym = W.params->rsym[e];
	isym = W.params->ssym[i];
	E = e - vir_off[esym];
	I = i - occ_off[isym];
	if ((M >= (occpi[msym] - openpi[msym])) ||
	    (E >= (virtpi[esym] - openpi[esym])) ||
	    (I >= (occpi[isym] - openpi[isym])) )
	  W.matrix[h][mn][ei] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);
  */


  /* Purge WMBIJ matrix elements */
  dpd_file4_init(&W, CC_HBAR, 0, 10, 2,"WMBIJ");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      b = W.params->roworb[h][mb][1];
      bsym = W.params->qsym[b];
      B = b - vir_off[bsym];
      for(ij=0; ij<W.params->coltot[h]; ij++) {
	if (B >= (virtpi[bsym] - openpi[bsym]))
	  W.matrix[h][mb][ij] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 10, 2,"Wmbij");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      m = W.params->roworb[h][mb][0];
      msym = W.params->psym[m];
      M = m - occ_off[msym];
      for(ij=0; ij<W.params->coltot[h]; ij++) {
	i = W.params->colorb[h][ij][0];
	j = W.params->colorb[h][ij][1];
	isym = W.params->rsym[i];
	jsym = W.params->ssym[j];
	I = i - occ_off[isym];
	J = j - occ_off[jsym];
	if ((M >= (occpi[msym] - openpi[msym])) ||
	    (I >= (occpi[isym] - openpi[isym])) ||
	    (J >= (occpi[jsym] - openpi[jsym])) )
	  W.matrix[h][mb][ij] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 10, 0,"WMbIj");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      for(ij=0; ij<W.params->coltot[h]; ij++) {
	j = W.params->colorb[h][ij][1];
	jsym = W.params->ssym[j];
	J = j - occ_off[jsym];
	if (J >= (occpi[jsym] - openpi[jsym]))
	  W.matrix[h][mb][ij] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 10, 0,"WmBiJ");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      m = W.params->roworb[h][mb][0];
      b = W.params->roworb[h][mb][1];
      msym = W.params->psym[m];
      bsym = W.params->qsym[b];
      M = m - occ_off[msym];
      B = b - vir_off[bsym];
      for(ij=0; ij<W.params->coltot[h]; ij++) {
	i = W.params->colorb[h][ij][0];
	isym = W.params->rsym[i];
	I = i - occ_off[isym];
	if ((M >= (occpi[msym] - openpi[msym])) ||
	    (B >= (virtpi[bsym] - openpi[bsym])) ||
	    (I >= (occpi[isym] - openpi[isym])) )
	  W.matrix[h][mb][ij] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);





  /* Purge Wabei matrix elements */
  dpd_file4_init(&W, CC_HBAR, 0, 11, 7,"WEIAB");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      e = W.params->roworb[h][ei][0];
      esym = W.params->psym[e];
      E = e - vir_off[esym];
      for(ab=0; ab<W.params->coltot[h]; ab++) {
	a = W.params->colorb[h][ab][0];
	b = W.params->colorb[h][ab][1];
	asym = W.params->rsym[a];
	bsym = W.params->ssym[b];
	A = a - vir_off[asym];
	B = b - vir_off[bsym];
	if ((E >= (virtpi[esym] - openpi[esym])) ||
	    (A >= (virtpi[asym] - openpi[asym])) ||
	    (B >= (virtpi[bsym] - openpi[bsym])) )
	  W.matrix[h][ei][ab] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 11, 7,"Weiab");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      i = W.params->roworb[h][ei][1];
      isym = W.params->qsym[i];
      I = i - occ_off[isym];
      for(ab=0; ab<W.params->coltot[h]; ab++) {
	if (I >= (occpi[isym] - openpi[isym]))
	  W.matrix[h][ei][ab] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 11, 5,"WEiAb");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      e = W.params->roworb[h][ei][0];
      i = W.params->roworb[h][ei][1];
      esym = W.params->psym[e];
      isym = W.params->qsym[i];
      E = e - vir_off[esym];
      I = i - occ_off[isym];
      for(ab=0; ab<W.params->coltot[h]; ab++) {
	a = W.params->colorb[h][ab][0];
	asym = W.params->rsym[a];
	bsym = W.params->ssym[b];
	A = a - vir_off[asym];
	if ((E >= (virtpi[esym] - openpi[esym])) ||
	    (I >= (occpi[isym] - openpi[isym])) ||
	    (A >= (virtpi[asym] - openpi[asym])) )
	  W.matrix[h][ei][ab] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 11, 5,"WeIaB");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      for(ab=0; ab<W.params->coltot[h]; ab++) {
	b = W.params->colorb[h][ab][1];
	bsym = W.params->ssym[b];
	B = b - vir_off[bsym];
	if (B >= (virtpi[bsym] - openpi[bsym]))
	  W.matrix[h][ei][ab] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  return;
}

/* Purge Wmnie matrix elements */
void purge_Wmnie(void) {
  dpdfile4 W;
  int *occpi, *virtpi;
  int h, a, b, e, f, i, j, m, n;
  int    A, B, E, F, I, J, M, N;
  int mn, ei, ma, ef, me, jb, mb, ij, ab;
  int asym, bsym, esym, fsym, isym, jsym, msym, nsym;
  int *occ_off, *vir_off;
  int *occ_sym, *vir_sym;
  int *openpi, nirreps;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  dpd_file4_init(&W, CC_HBAR, 0, 0, 11,"WMnIe (Mn,eI)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      n = W.params->roworb[h][mn][1];
      nsym = W.params->qsym[n];
      N = n - occ_off[nsym];
      for(ei=0; ei<W.params->coltot[h]; ei++) {
	if (N >= (occpi[nsym] - openpi[nsym]))
	  W.matrix[h][mn][ei] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }

  dpd_file4_init(&W, CC_HBAR, 0, 2, 11, "WMNIE (M>N,EI)");
  for(h=0; h < W.params->nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      for(ei=0; ei<W.params->coltot[h]; ei++) {
        e = W.params->colorb[h][ei][0];
        esym = W.params->rsym[e];
        E = e - vir_off[esym];
        if (E >= (virtpi[esym] - openpi[esym]))
          W.matrix[h][mn][ei] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 2, 11,"Wmnie (m>n,ei)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      m = W.params->roworb[h][mn][0];
      n = W.params->roworb[h][mn][1];
      msym = W.params->psym[m];
      nsym = W.params->qsym[n];
      M = m - occ_off[msym];
      N = n - occ_off[nsym];
      for(ei=0; ei<W.params->coltot[h]; ei++) {
        i = W.params->colorb[h][ei][1];
        isym = W.params->ssym[i];
        I = i - occ_off[isym];
        if ((M >= (occpi[msym] - openpi[msym])) ||
          (N >= (occpi[nsym] - openpi[nsym])) ||
          (I >= (occpi[isym] - openpi[isym])) )
          W.matrix[h][mn][ei] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_HBAR, 0, 0, 11,"WmNiE (mN,Ei)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      m = W.params->roworb[h][mn][0];
      msym = W.params->psym[m];
      M = m - occ_off[msym];
      for(ei=0; ei<W.params->coltot[h]; ei++) {
        e = W.params->colorb[h][ei][0];
        i = W.params->colorb[h][ei][1];
        esym = W.params->rsym[e];
        isym = W.params->ssym[i];
        E = e - vir_off[esym];
        I = i - occ_off[isym];
        if ((M >= (occpi[msym] - openpi[msym])) ||
            (E >= (virtpi[esym] - openpi[esym])) ||
            (I >= (occpi[isym] - openpi[isym])) )
          W.matrix[h][mn][ei] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);
  return;
}

}} // namespace psi::cchbar
