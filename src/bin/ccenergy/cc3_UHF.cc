/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void cc3_UHF_AAA_obsolete(void)
{
  int h, nirreps;
  int *occ_off, *occpi;
  int *vir_off, *virtpi;
  int Gi, Gj, Gk, Gijk;
  int Ga, Gb, Gc;
  int i, j, k, I, J, K;
  int a, b, c, A, B, C;
  int Gij, ij, Gab, ab, Gjk, jk;
  double ***W1;
  dpdbuf4 T2, E, F;
  dpdfile2 fIJ, fAB;
  dpdfile2 t1new, t1, d1;
  dpdbuf4 D;
  int nrows, ncols, nlinks;
  dpdbuf4 T2new, D2;
  dpdfile2 Fme;
  dpdbuf4 WMAFE, WMNIE;
  int Gd, d, cd, dc, Gid, id, DD;
  int Gm, m, Gmi, mi, im, mc, M;
  double **Z;

  nirreps = moinfo.nirreps;
  occpi = moinfo.aoccpi;
  occ_off = moinfo.aocc_off;
  virtpi = moinfo.avirtpi;
  vir_off = moinfo.avir_off;

  dpd_file2_init(&t1new, CC_OEI, 0, 0, 1, "CC3 tIA");
  dpd_file2_mat_init(&t1new);

  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&D, h);
    dpd_buf4_mat_irrep_rd(&D, h);
  }

  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_mat_init(&Fme);
  dpd_file2_mat_rd(&Fme);

  dpd_buf4_init(&T2new, CC_MISC, 0, 0, 5, 0, 5, 0, "CC3 tIJAB");
  for(h=0; h < nirreps; h++)
    dpd_buf4_mat_irrep_init(&T2new, h);

  dpd_buf4_init(&WMAFE, CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WAMEF (MA,F>E)");
  dpd_buf4_init(&WMNIE, CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMNIE (M>N,IE)");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&WMNIE, h);
    dpd_buf4_mat_irrep_rd(&WMNIE, h);
  }

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");

  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&F, CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WABEI (IE,B>A)");
  dpd_buf4_init(&E, CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMBIJ (I>J,MB)");

  /* target T3 amplitudes go in here */
  W1 = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      Gij = Gi ^ Gj;
      for(Gk=0; Gk < nirreps; Gk++) {
	Gijk = Gi ^ Gj ^ Gk;
	Gjk = Gj ^ Gk;

	for(Gab=0; Gab < nirreps; Gab++) {
	  /* This will need to change for non-totally-symmetric cases */
	  Gc = Gab ^ Gijk;
	  W1[Gab] = dpd_block_matrix(F.params->coltot[Gab], virtpi[Gc]);
	}

	for(i=0; i < occpi[Gi]; i++) {
	  I = occ_off[Gi] + i;
	  for(j=0; j < occpi[Gj]; j++) {
	    J = occ_off[Gj] + j;
	    for(k=0; k < occpi[Gk]; k++) {
	      K = occ_off[Gk] + k;

	      /* T3_AAA(W1, nirreps, I, Gi, J, Gj, K, Gk, &T2, &F, &E, &fIJ, &fAB, 
		     occpi, occ_off, virtpi, vir_off, 0.0); */

	      /* t_KC <-- 1/4 t_IJKABC <IJ||AB> */

	      Gc = Gk;    /* assumes T1 is totally symmetric */
	      Gab = Gij;  /* assumes <ij||ab> is totally symmetric */

	      ij = D.params->rowidx[I][J];

	      nrows = D.params->coltot[Gij];
	      ncols = virtpi[Gc];

	      if(nrows && ncols)
		C_DGEMV('t', nrows, ncols, 0.25, W1[Gab][0], ncols, D.matrix[Gij][ij], 1,
			1.0, t1new.matrix[Gk][k], 1);

	      /* t_IJAB <-- t_IJKABC F_KC */
	      Gc = Gk;    /* assumes Fme is totally symmetric */
	      Gab = Gij;  /* Assumes tIJAB is totally symmetric */

	      nrows = T2new.params->coltot[Gij];
	      ncols = virtpi[Gc];
	      ij = T2new.params->rowidx[I][J];

	      if(nrows && ncols)
		C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][0], ncols, Fme.matrix[Gk][k], 1,
			1.0, T2new.matrix[Gij][ij], 1);

	      /* t_JKDC <-- +1/2 t_IJKABC W_IDAB */
	      /* t_JKCD <-- -1/2 t_IJKABC W_IDAB */
	      jk = T2new.params->rowidx[J][K];
	      for(Gd=0; Gd < nirreps; Gd++) {
		Gab = Gid = Gi ^ Gd;  /* assumes Wieab is totally symmetric */
		Gc = Gab ^ Gijk;  /* assumes T3 is totally symmetric */

		id = WMAFE.row_offset[Gid][I];

		Z = block_matrix(virtpi[Gc],virtpi[Gd]);
		WMAFE.matrix[Gid] = dpd_block_matrix(virtpi[Gd], WMAFE.params->coltot[Gid]);
		dpd_buf4_mat_irrep_rd_block(&WMAFE, Gid, id, virtpi[Gd]);

		nrows = virtpi[Gc];
		ncols = virtpi[Gd];
		nlinks = WMAFE.params->coltot[Gid];

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 't', nrows, ncols, nlinks, 0.5, W1[Gab][0], nrows, 
			  WMAFE.matrix[Gid][0], nlinks, 0.0, Z[0], ncols);

		for(c=0; c < virtpi[Gc]; c++) {
		  C = vir_off[Gc] + c;
		  for(d=0; d < virtpi[Gd]; d++) {
		    DD = vir_off[Gd] + d;
		    cd = T2new.params->colidx[C][DD];
		    dc = T2new.params->colidx[DD][C];
		    T2new.matrix[Gjk][jk][dc] += Z[c][d];
		    T2new.matrix[Gjk][jk][cd] += -Z[c][d];
		  }
		}
		dpd_free_block(WMAFE.matrix[Gid], virtpi[Gd], WMAFE.params->coltot[Gid]);
		free_block(Z);
	      }

	      /* t_MIAB <-- +1/2 t_IJKABC W_JKMC */
	      /* t_IMAB <-- -1/2 t_IJKABC W_JKMC */
	      jk = WMNIE.params->rowidx[J][K];
	      for(Gm=0; Gm < nirreps; Gm++) {
		Gab = Gmi = Gm ^ Gi;  /* assumes totally symmetric */
		Gc = Gab ^ Gijk;      /* assumes totally symmetric */

		mc = WMNIE.col_offset[Gjk][Gm];

		nrows = F.params->coltot[Gab];
		ncols = occpi[Gm];
		nlinks = virtpi[Gc];

		Z = dpd_block_matrix(nrows, ncols);

		if(nrows && ncols && nlinks) 
		  C_DGEMM('n', 't', nrows, ncols, nlinks, 0.5, W1[Gab][0], nlinks,
			  &(WMNIE.matrix[Gjk][jk][mc]), nlinks, 0.0, Z[0], ncols);

		for(m=0; m < ncols; m++) {
		  M = occ_off[Gm] + m;
		  mi = T2new.params->rowidx[M][I];
		  im = T2new.params->rowidx[I][M];
		  for(ab=0; ab < nrows; ab++) {
		    T2new.matrix[Gmi][mi][ab] += Z[ab][m];
		    T2new.matrix[Gmi][im][ab] -= Z[ab][m];
		  }
		}

		dpd_free_block(Z, nrows, ncols);
	      }

	    } /* k */
	  } /* j */
	} /* i */

	for(Gab=0; Gab < nirreps; Gab++) {
	  /* This will need to change for non-totally-symmetric cases */
	  Gc = Gab ^ Gijk;
	  dpd_free_block(W1[Gab], F.params->coltot[Gab], virtpi[Gc]);
	}

      } /* Gk */
    } /* Gj */
  } /* Gi */

  free(W1);

  dpd_buf4_close(&E);
  dpd_buf4_close(&F);
  dpd_buf4_close(&T2);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);

  dpd_file2_mat_wrt(&t1new);
  dpd_file2_mat_close(&t1new);
  /* divide T1 contributions by denominators and add to T1 */
  dpd_file2_init(&d1, CC_OEI, 0, 0, 1, "dIA");
  dpd_file2_dirprd(&d1, &t1new);
  dpd_file2_close(&d1);
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "New tIA");
  dpd_file2_axpy(&t1new, &t1, 1, 0);
  dpd_file2_close(&t1);
  dpd_file2_close(&t1new);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&D, h);
  }
  dpd_buf4_close(&D);

  dpd_file2_mat_close(&Fme);
  dpd_file2_close(&Fme);

  dpd_buf4_close(&WMAFE);
  for(h=0; h < nirreps; h++)
    dpd_buf4_mat_irrep_close(&WMNIE, h);
  dpd_buf4_close(&WMNIE);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&T2new, h);
    dpd_buf4_mat_irrep_close(&T2new, h);
  }
  dpd_buf4_init(&D2, CC_DENOM, 0, 0, 5, 1, 6, 0, "dIJAB");
  dpd_buf4_dirprd(&D2, &T2new);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "New tIJAB");
  dpd_buf4_axpy(&T2new, &T2, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2new);
}

void cc3_UHF_BBB_obsolete(void)
{
  int h, nirreps;
  int *occ_off, *occpi;
  int *vir_off, *virtpi;
  int Gi, Gj, Gk, Gijk;
  int Ga, Gb, Gc;
  int i, j, k, I, J, K;
  int a, b, c, A, B, C;
  int Gij, ij, Gab, ab, Gjk, jk;
  double ***W1;
  dpdbuf4 T2, E, F;
  dpdfile2 fIJ, fAB;
  dpdfile2 t1new, t1, d1;
  dpdbuf4 D;
  int nrows, ncols, nlinks;
  dpdbuf4 T2new, D2;
  dpdfile2 Fme;
  dpdbuf4 WMAFE, WMNIE;
  int Gd, d, cd, dc, Gid, id, DD;
  int Gm, m, Gmi, mi, im, mc, M;
  double **Z;

  nirreps = moinfo.nirreps;
  occpi = moinfo.boccpi;
  occ_off = moinfo.bocc_off;
  virtpi = moinfo.bvirtpi;
  vir_off = moinfo.bvir_off;

  dpd_file2_init(&t1new, CC_OEI, 0, 2, 3, "CC3 tia");
  dpd_file2_mat_init(&t1new);

  dpd_buf4_init(&D, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&D, h);
    dpd_buf4_mat_irrep_rd(&D, h);
  }

  dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
  dpd_file2_mat_init(&Fme);
  dpd_file2_mat_rd(&Fme);

  dpd_buf4_init(&T2new, CC_MISC, 0, 10, 15, 10, 15, 0, "CC3 tijab");
  for(h=0; h < nirreps; h++) 
    dpd_buf4_mat_irrep_init(&T2new, h);

  dpd_buf4_init(&WMAFE, CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wamef (ma,f>e)");
  dpd_buf4_init(&WMNIE, CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmnie (m>n,ie)");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&WMNIE, h);
    dpd_buf4_mat_irrep_rd(&WMNIE, h);
  }

  dpd_file2_init(&fIJ, CC_OEI, 0, 2, 2, "fij");
  dpd_file2_init(&fAB, CC_OEI, 0, 3, 3, "fab");

  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
  dpd_buf4_init(&F, CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wabei (ie,b>a)");
  dpd_buf4_init(&E, CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmbij (i>j,mb)");

  /* target T3 amplitudes go in here */
  W1 = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      Gij = Gi ^ Gj;
      for(Gk=0; Gk < nirreps; Gk++) {
	Gijk = Gi ^ Gj ^ Gk;
	Gjk = Gj ^ Gk;

	for(Gab=0; Gab < nirreps; Gab++) {
	  /* This will need to change for non-totally-symmetric cases */
	  Gc = Gab ^ Gijk;
	  W1[Gab] = dpd_block_matrix(F.params->coltot[Gab], virtpi[Gc]);
	}

	for(i=0; i < occpi[Gi]; i++) {
	  I = occ_off[Gi] + i;
	  for(j=0; j < occpi[Gj]; j++) {
	    J = occ_off[Gj] + j;
	    for(k=0; k < occpi[Gk]; k++) {
	      K = occ_off[Gk] + k;

	      /* T3_AAA(W1, nirreps, I, Gi, J, Gj, K, Gk, &T2, &F, &E, &fIJ, &fAB, 
		     occpi, occ_off, virtpi, vir_off, 0.0); */

	      /* t_kc <-- 1/4 t_ijkabc <ij||ab> */

	      Gc = Gk;    /* assumes T1 is totally symmetric */
	      Gab = Gij;  /* assumes <ij||ab> is totally symmetric */

	      ij = D.params->rowidx[I][J];

	      nrows = D.params->coltot[Gij];
	      ncols = virtpi[Gc];

	      if(nrows && ncols)
		C_DGEMV('t', nrows, ncols, 0.25, W1[Gab][0], ncols, D.matrix[Gij][ij], 1,
			1.0, t1new.matrix[Gk][k], 1);

	      /* t_ijab <-- t_ijkabc F_kc */
	      Gc = Gk;    /* assumes Fme is totally symmetric */
	      Gab = Gij;  /* Assumes tijab is totally symmetric */

	      nrows = T2new.params->coltot[Gij];
	      ncols = virtpi[Gc];
	      ij = T2new.params->rowidx[I][J];

	      if(nrows && ncols)
		C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][0], ncols, Fme.matrix[Gk][k], 1,
			1.0, T2new.matrix[Gij][ij], 1);

	      /* t_jkdc <-- 1/2 t_ijkabc W_idab */
	      /* t_jkcd <-- -1/2 t_ijkabc W_idab */
	      jk = T2new.params->rowidx[J][K];
	      for(Gd=0; Gd < nirreps; Gd++) {
		Gab = Gid = Gi ^ Gd;  /* assumes Wieab is totally symmetric */
		Gc = Gab ^ Gijk;  /* assumes T3 is totally symmetric */

		id = WMAFE.row_offset[Gid][I];

		Z = block_matrix(virtpi[Gc],virtpi[Gd]);
		WMAFE.matrix[Gid] = dpd_block_matrix(virtpi[Gd], WMAFE.params->coltot[Gid]);
		dpd_buf4_mat_irrep_rd_block(&WMAFE, Gid, id, virtpi[Gd]);

		nrows = virtpi[Gc];
		ncols = virtpi[Gd];
		nlinks = WMAFE.params->coltot[Gid];

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 't', nrows, ncols, nlinks, 0.5, W1[Gab][0], nrows, 
			  WMAFE.matrix[Gid][0], nlinks, 0.0, Z[0], ncols);

		for(c=0; c < virtpi[Gc]; c++) {
		  C = vir_off[Gc] + c;
		  for(d=0; d < virtpi[Gd]; d++) {
		    DD = vir_off[Gd] + d;
		    cd = T2new.params->colidx[C][DD];
		    dc = T2new.params->colidx[DD][C];
		    T2new.matrix[Gjk][jk][dc] += Z[c][d];
		    T2new.matrix[Gjk][jk][cd] += -Z[c][d];
		  }
		}
		dpd_free_block(WMAFE.matrix[Gid], virtpi[Gd], WMAFE.params->coltot[Gid]);
		free_block(Z);
	      }

	      /* t_miab <-- +1/2 t_ijkabc W_jkmc */
	      /* t_imab <-- -1/2 t_ijkabc W_jkmc */
	      jk = WMNIE.params->rowidx[J][K];
	      for(Gm=0; Gm < nirreps; Gm++) {
		Gab = Gmi = Gm ^ Gi;  /* assumes totally symmetric */
		Gc = Gab ^ Gijk;      /* assumes totally symmetric */

		mc = WMNIE.col_offset[Gjk][Gm];

		nrows = F.params->coltot[Gab];
		ncols = occpi[Gm];
		nlinks = virtpi[Gc];

		Z = dpd_block_matrix(nrows, ncols);

		if(nrows && ncols && nlinks) 
		  C_DGEMM('n', 't', nrows, ncols, nlinks, 0.5, W1[Gab][0], nlinks,
			  &(WMNIE.matrix[Gjk][jk][mc]), nlinks, 0.0, Z[0], ncols);

		for(m=0; m < ncols; m++) {
		  M = occ_off[Gm] + m;
		  mi = T2new.params->rowidx[M][I];
		  im = T2new.params->rowidx[I][M];
		  for(ab=0; ab < nrows; ab++) {
		    T2new.matrix[Gmi][mi][ab] += Z[ab][m];
		    T2new.matrix[Gmi][im][ab] -= Z[ab][m];
		  }
		}

		dpd_free_block(Z, nrows, ncols);
	      }

	    } /* k */
	  } /* j */
	} /* i */

	for(Gab=0; Gab < nirreps; Gab++) {
	  /* This will need to change for non-totally-symmetric cases */
	  Gc = Gab ^ Gijk;
	  dpd_free_block(W1[Gab], F.params->coltot[Gab], virtpi[Gc]);
	}

      } /* Gk */
    } /* Gj */
  } /* Gi */

  free(W1);

  dpd_buf4_close(&E);
  dpd_buf4_close(&F);
  dpd_buf4_close(&T2);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);

  dpd_file2_mat_wrt(&t1new);
  dpd_file2_mat_close(&t1new);
  /* divide T1 contributions by denominators and add to T1 */
  dpd_file2_init(&d1, CC_OEI, 0, 2, 3, "dia");
  dpd_file2_dirprd(&d1, &t1new);
  dpd_file2_close(&d1);
  dpd_file2_init(&t1, CC_OEI, 0, 2, 3, "New tia");
  dpd_file2_axpy(&t1new, &t1, 1, 0);
  dpd_file2_close(&t1);
  dpd_file2_close(&t1new);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&D, h);
  }
  dpd_buf4_close(&D);

  dpd_file2_mat_close(&Fme);
  dpd_file2_close(&Fme);

  dpd_buf4_close(&WMAFE);
  for(h=0; h < nirreps; h++)
    dpd_buf4_mat_irrep_close(&WMNIE, h);
  dpd_buf4_close(&WMNIE);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&T2new, h);
    dpd_buf4_mat_irrep_close(&T2new, h);
  }
  dpd_buf4_init(&D2, CC_DENOM, 0, 10, 15, 11, 16, 0, "dijab");
  dpd_buf4_dirprd(&D2, &T2new);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 15, 12, 17, 0, "New tijab");
  dpd_buf4_axpy(&T2new, &T2, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2new);
}

void cc3_UHF_AAB_obsolete(void)
{
  int h, nirreps;
  int *aocc_off, *aoccpi;
  int *avir_off, *avirtpi;
  int *bocc_off, *boccpi;
  int *bvir_off, *bvirtpi;
  int Gi, Gj, Gk, Gijk;
  int Ga, Gb, Gc, Gab;
  int Gij, ij, ji, Gjk, jk, Gbc, bc, Gcb, cb;
  int Gd, Gkd, kd, d, DD, ad, da, dc, Gid, id, Gac, ac, bd;
  int i, j, k, I, J, K;
  int a, b, c, A, B, C;
  int ab;
  double ***W1, ***W2, ***W3;
  dpdbuf4 T2AA, T2AB, T2BA, EAA, EAB, EBA, FAA, FAB, FBA;
  dpdfile2 fIJ, fAB, fij, fab;
  dpdfile2 t1newA, t1newB;
  dpdfile2 t1, d1;
  dpdbuf4 DAA, DAB;
  int nrows, ncols, nlinks;
  int **W_offset, offset;
  dpdfile2 FME, Fme;
  dpdbuf4 T2AAnew, D2, T2;
  dpdbuf4 T2ABnew;
  dpdbuf4 WmAfE, WMAFE, WMaFe;
  double **Z;
  dpdbuf4 WMnIe, WMNIE, WmNiE;
  int Gm, m, Gmi, mi, im, mc, M;
  int Gmk, mk, ma, kj, Gim;

  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi;
  aocc_off = moinfo.aocc_off;
  boccpi = moinfo.boccpi;
  bocc_off = moinfo.bocc_off;
  avirtpi = moinfo.avirtpi;
  avir_off = moinfo.avir_off;
  bvirtpi = moinfo.bvirtpi;
  bvir_off = moinfo.bvir_off;

  W_offset = init_int_matrix(nirreps, nirreps);
  for(Gab=0; Gab < nirreps; Gab++) {
    for(Ga=0,offset=0; Ga < nirreps; Ga++) {
      Gb = Ga ^ Gab;
      W_offset[Gab][Ga] = offset;
      offset += avirtpi[Ga] * avirtpi[Gb];
    }
  }

  dpd_file2_init(&t1newA, CC_OEI, 0, 0, 1, "CC3 tIA");
  dpd_file2_mat_init(&t1newA);
  dpd_file2_init(&t1newB, CC_OEI, 0, 2, 3, "CC3 tia");
  dpd_file2_mat_init(&t1newB);
 
  dpd_buf4_init(&DAA, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
  dpd_buf4_init(&DAB, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&DAA, h);
    dpd_buf4_mat_irrep_rd(&DAA, h);
    dpd_buf4_mat_irrep_init(&DAB, h);
    dpd_buf4_mat_irrep_rd(&DAB, h);
  }

  dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_mat_init(&FME);
  dpd_file2_mat_rd(&FME);

  dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
  dpd_file2_mat_init(&Fme);
  dpd_file2_mat_rd(&Fme);

  dpd_buf4_init(&WmAfE, CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAmEf (mA,fE)");
  dpd_buf4_init(&WMAFE, CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WAMEF (MA,F>E)");
  dpd_buf4_init(&WMaFe, CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaMeF (Ma,Fe)");

  dpd_buf4_init(&WMNIE, CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMNIE (M>N,IE)");
  dpd_buf4_init(&WMnIe, CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMnIe (Mn,Ie)");
  dpd_buf4_init(&WmNiE, CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmNiE (mN,iE)");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&WMnIe, h);
    dpd_buf4_mat_irrep_rd(&WMnIe, h);

    dpd_buf4_mat_irrep_init(&WMNIE, h);
    dpd_buf4_mat_irrep_rd(&WMNIE, h);

    dpd_buf4_mat_irrep_init(&WmNiE, h);
    dpd_buf4_mat_irrep_rd(&WmNiE, h);
  }

  dpd_buf4_init(&T2AAnew, CC_MISC, 0, 0, 5, 0, 5, 0, "CC3 tIJAB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2AAnew, h);
  }

  dpd_buf4_init(&T2ABnew, CC_MISC, 0, 22, 28, 22, 28, 0, "CC3 tIjAb");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2ABnew, h);
  }

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&fij, CC_OEI, 0, 2, 2, "fij");
  dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab");

  dpd_buf4_init(&T2AA, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&T2AB, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_buf4_init(&T2BA, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
  dpd_buf4_init(&FAA, CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WABEI (IE,B>A)");
  dpd_buf4_init(&FAB, CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaBeI (Ie,Ba)");
  dpd_buf4_init(&FBA, CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAbEi (iE,bA)");
  dpd_buf4_init(&EAA, CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMBIJ (I>J,MB)");
  dpd_buf4_init(&EAB, CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMbIj (Ij,Mb)");
  dpd_buf4_init(&EBA, CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmBiJ (iJ,mB)");

  /* target T3 amplitudes go in here */
  W1 = (double ***) malloc(nirreps * sizeof(double **));
  W2 = (double ***) malloc(nirreps * sizeof(double **));
  W3 = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      Gij = Gi ^ Gj;
      for(Gk=0; Gk < nirreps; Gk++) {
	Gijk = Gi ^ Gj ^ Gk;
	Gjk = Gj ^ Gk;

	for(Gab=0; Gab < nirreps; Gab++) {
	  Gc = Gab ^ Gijk; /* assumes totally symmetric */
	  W1[Gab] = dpd_block_matrix(FAA.params->coltot[Gab], bvirtpi[Gc]);
	}
	for(Ga=0; Ga < nirreps; Ga++) {
	  Gcb = Ga ^ Gijk; /* assumes totally symmetric */
	  W2[Ga] = dpd_block_matrix(avirtpi[Ga], WmAfE.params->coltot[Gcb]);  /* alpha-beta-alpha */
	}
	for(Gb=0; Gb < nirreps; Gb++) {
	  Gac = Gb ^ Gijk; /* assumes totally symmetric */
	  W3[Gb] = dpd_block_matrix(avirtpi[Gb], WMaFe.params->coltot[Gac]);  /* alpha-alpha-beta */
	}

	for(i=0; i < aoccpi[Gi]; i++) {
	  I = aocc_off[Gi] + i;
	  for(j=0; j < aoccpi[Gj]; j++) {
	    J = aocc_off[Gj] + j;
	    for(k=0; k < boccpi[Gk]; k++) {
	      K = bocc_off[Gk] + k;

	     /* T3_AAB(W1, nirreps, I, Gi, J, Gj, K, Gk, &T2AA, &T2AB, &T2BA, 
		     &FAA, &FAB, &FBA, &EAA, &EAB, &EBA, &fIJ, &fij, &fAB, &fab,
		     aoccpi, aocc_off, boccpi, bocc_off, avirtpi, avir_off, bvirtpi, bvir_off, 0.0);
*/

	      /* t_kc <-- 1/4 t_IJkABc <IJ||AB> */

	      Gc = Gk;   /* assumes T1 is totally symmetric */
	      Gab = Gij; /* assumes <ij||ab> is totally symmetric */

	      ij = DAA.params->rowidx[I][J];

	      nrows = DAA.params->coltot[Gij];
	      ncols = bvirtpi[Gc];

	      if(nrows && ncols)
		C_DGEMV('t', nrows, ncols, 0.25, W1[Gab][0], ncols, DAA.matrix[Gij][ij], 1,
			1.0, t1newB.matrix[Gk][k], 1);

	      /* t_IA <-- t_IJkABc <Jk|Bc> */

	      Ga = Gi;   /* assumes T1 is totally symmetric */
	      Gbc = Gjk; /* assumes <jk|bc> is totally symmetric */

	      jk = DAB.params->rowidx[J][K];

	      for(Gab=0; Gab < nirreps; Gab++) {
		Gb = Ga ^ Gab;
		Gc = Gb ^ Gbc;

		ab = W_offset[Gab][Ga];
		bc = DAB.col_offset[Gjk][Gb];

		nrows = avirtpi[Ga];
		ncols = avirtpi[Gb] * bvirtpi[Gc];

		if(nrows && ncols)
		  C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][ab], ncols, &(DAB.matrix[Gjk][jk][bc]), 1,
			  1.0, t1newA.matrix[Gi][i], 1);
	      }

	      /* t_IJAB <-- t_IJkABc F_kc */
	      Gc = Gk;    /* assumes Fme is totally symmetric */
	      Gab = Gij;  /* Assumes tIJAB is totally symmetric */

	      nrows = T2AAnew.params->coltot[Gij];
	      ncols = bvirtpi[Gc];
	      ij = T2AAnew.params->rowidx[I][J];
 
	      if(nrows && ncols)
		C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][0], ncols, Fme.matrix[Gk][k], 1,
			1.0, T2AAnew.matrix[Gij][ij], 1);

	      /* t_JkBc <-- t_IJkABc F_IA */
	      Ga = Gi;   /* assumes Fia is totally symmetric */
	      Gbc = Gjk; /* assumes t_jKbC is totally symmetric */

	      jk = T2AB.params->rowidx[J][K];

	      for(Gab=0; Gab < nirreps; Gab++) {
		Gb = Ga ^ Gab;
		Gc = Gb ^ Gbc;

		ab = W_offset[Gab][Ga];
		bc = T2ABnew.col_offset[Gjk][Gb];

		nrows = avirtpi[Ga];
		ncols = avirtpi[Gb] * bvirtpi[Gc];

		if(nrows && ncols)
		  C_DGEMV('t', nrows, ncols, 1.0, W1[Gab][ab], ncols, FME.matrix[Gi][i], 1,
			  1.0, &(T2ABnew.matrix[Gjk][jk][bc]), 1);
	      }

	      /* t_JIDA <-- t_IJkABc W_kDcB */
	      /* sort W1(AB,c) to W2(A,cB) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++) {
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Ga = FAA.params->rsym[A];
		  a = A - avir_off[Ga];
		  for(c=0; c < bvirtpi[Gc]; c++) {
		    C = bvir_off[Gc] + c;
		    cb = WmAfE.params->colidx[C][B];
		    W2[Ga][a][cb] = W1[Gab][ab][c];
		  }
		}
	      }

	      ji = T2AAnew.params->rowidx[J][I];

	      for(Gd=0; Gd < nirreps; Gd++) {
		Gcb = Gkd = Gk ^ Gd; /* assumes totally symmetric */
		Ga = Gd ^ Gij;       /* assumes totally symmetric */

		kd = WmAfE.row_offset[Gkd][K];
		WmAfE.matrix[Gkd] = dpd_block_matrix(avirtpi[Gd], WmAfE.params->coltot[Gkd]);
		dpd_buf4_mat_irrep_rd_block(&WmAfE, Gkd, kd, avirtpi[Gd]);
		Z = block_matrix(avirtpi[Ga], avirtpi[Gd]);

		nrows = avirtpi[Ga];
		ncols = avirtpi[Gd];
		nlinks = WmAfE.params->coltot[Gkd];

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W2[Ga][0], nlinks, 
			  WmAfE.matrix[Gkd][0], nlinks, 0.0, Z[0], ncols);

		for(a=0; a < avirtpi[Ga]; a++) {
		  A = avir_off[Ga] + a;
		  for(d=0; d < avirtpi[Gd]; d++) {
		    DD = avir_off[Gd] + d;
		    ad = T2AAnew.params->colidx[A][DD];
		    da = T2AAnew.params->colidx[DD][A];
		    T2AAnew.matrix[Gij][ji][ad] += -Z[a][d];
		    T2AAnew.matrix[Gij][ji][da] += Z[a][d];
		  }
		}

		dpd_free_block(WmAfE.matrix[Gkd], avirtpi[Gd], WmAfE.params->coltot[Gkd]);
		free_block(Z);
	      }

	      /* t_JkDc <-- 1/2 t_IJkABc W_IDAB */

	      jk = T2ABnew.params->rowidx[J][K];

	      for(Gd=0; Gd < nirreps; Gd++) {
		Gab = Gid = Gi ^ Gd; /* assumes totally symmetric */
		Gc = Gab ^ Gijk;     /* assumes totally symmetric */

		id = WMAFE.row_offset[Gid][I];
		WMAFE.matrix[Gid] = dpd_block_matrix(avirtpi[Gd], WMAFE.params->coltot[Gid]);
		dpd_buf4_mat_irrep_rd_block(&WMAFE, Gid, id, avirtpi[Gd]);
		Z = block_matrix(bvirtpi[Gc], avirtpi[Gd]);

		nrows = bvirtpi[Gc];
		ncols = avirtpi[Gd];
		nlinks = WMAFE.params->coltot[Gid];

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 't', nrows, ncols, nlinks, 0.5, W1[Gab][0], nrows,
			  WMAFE.matrix[Gid][0], nlinks, 0.0, Z[0], ncols);

		for(c=0; c < bvirtpi[Gc]; c++) {
		  C = bvir_off[Gc] + c;
		  for(d=0; d < avirtpi[Gd]; d++) {
		    DD = avir_off[Gd] + d;
		    dc = T2ABnew.params->colidx[DD][C];
		    T2ABnew.matrix[Gjk][jk][dc] += Z[c][d];
		  }
		}

		free_block(Z);
		dpd_free_block(WMAFE.matrix[Gid], avirtpi[Gd], WMAFE.params->coltot[Gid]);
	      }

	      /* t_JkBd <-- t_IJkABc W_IdAc */
	      /* sort W1(AB,c) to W3(B,Ac) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++) {
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Gb = FAA.params->ssym[B];
		  b = B - avir_off[Gb];
		  for(c=0; c < bvirtpi[Gc]; c++) {
		    C = bvir_off[Gc] + c;
		    ac = WMaFe.params->colidx[A][C];
		    W3[Gb][b][ac] = W1[Gab][ab][c];
		  }
		}
	      }

	      jk = T2ABnew.params->rowidx[J][K];

	      for(Gd=0; Gd < nirreps; Gd++) {
		Gac = Gid = Gi ^ Gd; /* assumes totally symmetric */
		Gb = Gac ^ Gijk;     /* assumes totally symmetric */

		id = WMaFe.row_offset[Gid][I]; 
		WMaFe.matrix[Gid] = dpd_block_matrix(bvirtpi[Gd], WMaFe.params->coltot[Gid]);
		dpd_buf4_mat_irrep_rd_block(&WMaFe, Gid, id, bvirtpi[Gd]);
		Z = block_matrix(avirtpi[Gb], bvirtpi[Gd]);

		nrows = avirtpi[Gb];
		ncols = bvirtpi[Gd];
		nlinks = WMaFe.params->coltot[Gid];

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W3[Gb][0], nlinks,
			  WMaFe.matrix[Gid][0], nlinks, 0.0, Z[0], ncols);

		for(b=0; b < avirtpi[Gb]; b++) {
		  B = avir_off[Gb] + b;
		  for(d=0; d < bvirtpi[Gd]; d++) {
		    DD = bvir_off[Gd] + d;
		    bd = T2ABnew.params->colidx[B][DD];
		    T2ABnew.matrix[Gjk][jk][bd] += Z[b][d];
		  }
		}

		dpd_free_block(WMaFe.matrix[Gid], bvirtpi[Gd], WMaFe.params->coltot[Gid]);
		free_block(Z);
	      }

	      /* t_MIAB <--- +t_IJkABc W_JkMc */
	      /* t_IMAB <--- -t_IJkABc W_JkMc */

	      jk = WMnIe.params->rowidx[J][K];

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gab = Gmi = Gm ^ Gi;  /* assumes totally symmetric */
		Gc = Gab ^ Gijk;      /* assumes totally symmetric */

		mc = WMnIe.col_offset[Gjk][Gm];

		nrows = FAA.params->coltot[Gab];
		ncols = aoccpi[Gm];
		nlinks = bvirtpi[Gc];

		Z = dpd_block_matrix(nrows, ncols);

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W1[Gab][0], nlinks,
			  &(WMnIe.matrix[Gjk][jk][mc]), nlinks, 0.0, Z[0], ncols);

		for(m=0; m < ncols; m++) {
		  M = aocc_off[Gm] + m;
		  mi = T2AAnew.params->rowidx[M][I];
		  im = T2AAnew.params->rowidx[I][M];
		  for(ab=0; ab < nrows; ab++) {
		    T2AAnew.matrix[Gmi][mi][ab] += Z[ab][m];
		    T2AAnew.matrix[Gmi][im][ab] -= Z[ab][m];
		  }
		}

		dpd_free_block(Z, nrows, ncols);
	      }

	      /* t_MkBc <-- 1/2 t_IJkABc W_IJMA */
	      /* sort W(AB,c) to W(A,Bc) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;  /* assumes totally symmetric */
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++ ){
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Ga = FAA.params->rsym[A];
		  a = A - avir_off[Ga];
		  for(c=0; c < bvirtpi[Gc]; c++) {
		    C = bvir_off[Gc] + c;
		    bc = T2ABnew.params->colidx[B][C];
		    W3[Ga][a][bc] = W1[Gab][ab][c];
		  }
		}
	      }

	      ij = WMNIE.params->rowidx[I][J];

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gbc = Gmk = Gm ^ Gk;  /* assumes totally symmetric */
		Ga = Gbc ^ Gijk;      /* assumes totally symmetric */

		ma = WMNIE.col_offset[Gij][Gm];

		nrows = T2ABnew.params->coltot[Gmk];
		ncols = aoccpi[Gm];
		nlinks = avirtpi[Ga];

		Z = dpd_block_matrix(nrows, ncols);

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 't', nrows, ncols, nlinks, 0.5, W3[Ga][0], nrows,
			  &(WMNIE.matrix[Gij][ij][ma]), nlinks, 0.0, Z[0], ncols);

		for(m=0; m < aoccpi[Gm]; m++) {
		  M = aocc_off[Gm] + m;
		  mk = T2ABnew.params->rowidx[M][K];
		  for(bc=0; bc < nrows; bc++) {
		    T2ABnew.matrix[Gmk][mk][bc] += Z[bc][m];
		  }
		}

		dpd_free_block(Z, nrows, ncols);
	      }

	      /* t_ImBc <-- t_IJkABc W_kJmA */
	      /* sort W(AB,c) to W(A,Bc) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;  /* assumes totally symmetric */
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++ ){
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Ga = FAA.params->rsym[A];
		  a = A - avir_off[Ga];
		  for(c=0; c < bvirtpi[Gc]; c++) {
		    C = bvir_off[Gc] + c;
		    bc = T2ABnew.params->colidx[B][C];
		    W3[Ga][a][bc] = W1[Gab][ab][c];
		  }
		}
	      }

	      kj = WmNiE.params->rowidx[K][J];

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gbc = Gim = Gi ^ Gm;  /* assumes totally symmetric */
		Ga = Gbc ^ Gijk;      /* assumes totally symmetric */

		ma = WmNiE.col_offset[Gjk][Gm];

		nrows = T2ABnew.params->coltot[Gim];
		ncols = boccpi[Gm];
		nlinks = avirtpi[Ga];

		Z = dpd_block_matrix(nrows, ncols);

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, W3[Ga][0], nrows,
			  &(WmNiE.matrix[Gjk][kj][ma]), nlinks, 0.0, Z[0], ncols);

		for(m=0; m < boccpi[Gm]; m++) {
		  M = bocc_off[Gm] + m;
		  im = T2ABnew.params->rowidx[I][M];
		  for(bc=0; bc < nrows; bc++) {
		    T2ABnew.matrix[Gim][im][bc] += Z[bc][m];
		  }
		}

		dpd_free_block(Z, nrows, ncols);
	      }

	    } /* k */
	  } /* j */
	} /* i */

	for(Gab=0; Gab < nirreps; Gab++) {
	  /* This will need to change for non-totally-symmetric cases */
	  Gc = Gab ^ Gijk;
	  dpd_free_block(W1[Gab], FAA.params->coltot[Gab], bvirtpi[Gc]);
	}
	for(Ga=0; Ga < nirreps; Ga++) {
	  Gcb = Ga ^ Gijk; /* assumes totally symmetric */
	  dpd_free_block(W2[Ga], avirtpi[Ga], WmAfE.params->coltot[Gcb]);
	}
	for(Gb=0; Gb < nirreps; Gb++) {
	  Gac = Gb ^ Gijk; /* assumes totally symmetric */
	  dpd_free_block(W3[Gb], avirtpi[Gb], WMaFe.params->coltot[Gac]);
	}

      } /* Gk */
    } /* Gj */
  } /* Gi */

  free(W1);
  free(W2);
  free(W3);

  dpd_buf4_close(&EAA);
  dpd_buf4_close(&EAB);
  dpd_buf4_close(&EBA);
  dpd_buf4_close(&FAA);
  dpd_buf4_close(&FAB);
  dpd_buf4_close(&FBA);
  dpd_buf4_close(&T2AA);
  dpd_buf4_close(&T2AB);
  dpd_buf4_close(&T2BA);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);
  dpd_file2_close(&fij);
  dpd_file2_close(&fab);

  dpd_file2_mat_wrt(&t1newA);
  dpd_file2_init(&d1, CC_OEI, 0, 0, 1, "dIA");
  dpd_file2_dirprd(&d1, &t1newA);
  dpd_file2_close(&d1);
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "New tIA");
  dpd_file2_axpy(&t1newA, &t1, 1, 0);
  dpd_file2_close(&t1);
  dpd_file2_close(&t1newA);

  dpd_file2_mat_wrt(&t1newB);
  dpd_file2_init(&d1, CC_OEI, 0, 2, 3, "dia");
  dpd_file2_dirprd(&d1, &t1newB);
  dpd_file2_close(&d1);
  dpd_file2_init(&t1, CC_OEI, 0, 2, 3, "New tia");
  dpd_file2_axpy(&t1newB, &t1, 1, 0);
  dpd_file2_close(&t1);
  dpd_file2_close(&t1newB);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&DAA, h);
    dpd_buf4_mat_irrep_close(&DAB, h);
  }
  dpd_buf4_close(&DAA);
  dpd_buf4_close(&DAB);

  free_int_matrix(W_offset);

  dpd_file2_mat_close(&FME);
  dpd_file2_close(&FME);

  dpd_file2_mat_close(&Fme);
  dpd_file2_close(&Fme);

  dpd_buf4_close(&WmAfE);
  dpd_buf4_close(&WMAFE);
  dpd_buf4_close(&WMaFe);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&WMnIe, h);
    dpd_buf4_mat_irrep_close(&WMNIE, h);
    dpd_buf4_mat_irrep_close(&WmNiE, h);
  }
  dpd_buf4_close(&WMnIe);
  dpd_buf4_close(&WMNIE);
  dpd_buf4_close(&WmNiE);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&T2AAnew, h);
    dpd_buf4_mat_irrep_close(&T2AAnew, h);
  }
  dpd_buf4_init(&D2, CC_DENOM, 0, 0, 5, 1, 6, 0, "dIJAB");
  dpd_buf4_dirprd(&D2, &T2AAnew);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "New tIJAB");
  dpd_buf4_axpy(&T2AAnew, &T2, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2AAnew);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&T2ABnew, h);
    dpd_buf4_mat_irrep_close(&T2ABnew, h);
  }
  dpd_buf4_init(&D2, CC_DENOM, 0, 22, 28, 22, 28, 0, "dIjAb");
  dpd_buf4_dirprd(&D2, &T2ABnew);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
  dpd_buf4_axpy(&T2ABnew, &T2, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2ABnew);
}

void cc3_UHF_BBA_obsolete(void)
{
  int h, nirreps;
  int *aocc_off, *aoccpi;
  int *avir_off, *avirtpi;
  int *bocc_off, *boccpi;
  int *bvir_off, *bvirtpi;
  int Gi, Gj, Gk, Gijk;
  int Ga, Gb, Gc, Gab;
  int Gij, ij, ji, Gjk, jk, Gbc, bc, Gcb, cb;
  int Gd, Gkd, kd, d, DD, ad, da, dc, Gid, id, Gac, ac, bd;
  int i, j, k, I, J, K;
  int a, b, c, A, B, C;
  int ab;
  double ***W1, ***W2, ***W3;
  dpdbuf4 T2BB, T2AB, T2BA, EBB, EAB, EBA, FBB, FAB, FBA;
  dpdfile2 fIJ, fAB, fij, fab;
  dpdfile2 t1newA, t1newB;
  dpdfile2 t1, d1;
  dpdbuf4 DBB, DBA;
  int nrows, ncols, nlinks;
  int **W_offset, offset;
  dpdfile2 FME, Fme;
  dpdbuf4 T2BBnew, D2, T2;
  dpdbuf4 T2BAnew;
  dpdbuf4 WMaFe, Wmafe, WmAfE;
  double **Z;
  dpdbuf4 WmNiE, Wmnie, WMnIe;
  int Gm, m, Gmi, mi, im, mc, M;
  int Gmk, mk, ma, kj, Gim;

  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi;
  aocc_off = moinfo.aocc_off;
  boccpi = moinfo.boccpi;
  bocc_off = moinfo.bocc_off;
  avirtpi = moinfo.avirtpi;
  avir_off = moinfo.avir_off;
  bvirtpi = moinfo.bvirtpi;
  bvir_off = moinfo.bvir_off;

  W_offset = init_int_matrix(nirreps, nirreps);
  for(Gab=0; Gab < nirreps; Gab++) {
    for(Ga=0,offset=0; Ga < nirreps; Ga++) {
      Gb = Ga ^ Gab;
      W_offset[Gab][Ga] = offset;
      offset += bvirtpi[Ga] * bvirtpi[Gb];
    }
  }

  dpd_file2_init(&t1newA, CC_OEI, 0, 0, 1, "CC3 tIA");
  dpd_file2_mat_init(&t1newA);
  dpd_file2_init(&t1newB, CC_OEI, 0, 2, 3, "CC3 tia");
  dpd_file2_mat_init(&t1newB);

  dpd_buf4_init(&DBB, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
  dpd_buf4_init(&DBA, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&DBB, h);
    dpd_buf4_mat_irrep_rd(&DBB, h);
    dpd_buf4_mat_irrep_init(&DBA, h);
    dpd_buf4_mat_irrep_rd(&DBA, h);
  }

  dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_mat_init(&FME);
  dpd_file2_mat_rd(&FME);

  dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
  dpd_file2_mat_init(&Fme);
  dpd_file2_mat_rd(&Fme);

  dpd_buf4_init(&WMaFe, CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaMeF (Ma,Fe)");
  dpd_buf4_init(&Wmafe, CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wamef (ma,f>e)");
  dpd_buf4_init(&WmAfE, CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAmEf (mA,fE)");

  dpd_buf4_init(&Wmnie, CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmnie (m>n,ie)");
  dpd_buf4_init(&WmNiE, CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmNiE (mN,iE)");
  dpd_buf4_init(&WMnIe, CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMnIe (Mn,Ie)");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&WmNiE, h);
    dpd_buf4_mat_irrep_rd(&WmNiE, h);

    dpd_buf4_mat_irrep_init(&Wmnie, h);
    dpd_buf4_mat_irrep_rd(&Wmnie, h);

    dpd_buf4_mat_irrep_init(&WMnIe, h);
    dpd_buf4_mat_irrep_rd(&WMnIe, h);
  }

  dpd_buf4_init(&T2BBnew, CC_MISC, 0, 10, 15, 10, 15, 0, "CC3 tijab");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2BBnew, h);
  }

  dpd_buf4_init(&T2BAnew, CC_MISC, 0, 23, 29, 23, 29, 0, "CC3 tiJaB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2BAnew, h);
  }

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&fij, CC_OEI, 0, 2, 2, "fij");
  dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab");

  dpd_buf4_init(&T2BB, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
  dpd_buf4_init(&T2AB, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_buf4_init(&T2BA, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
  dpd_buf4_init(&FBB, CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wabei (ie,b>a)");
  dpd_buf4_init(&FAB, CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaBeI (Ie,Ba)");
  dpd_buf4_init(&FBA, CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAbEi (iE,bA)");
  dpd_buf4_init(&EBB, CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmbij (i>j,mb)");
  dpd_buf4_init(&EAB, CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMbIj (Ij,Mb)");
  dpd_buf4_init(&EBA, CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmBiJ (iJ,mB)");

  /* target T3 amplitudes go in here */
  W1 = (double ***) malloc(nirreps * sizeof(double **));
  W2 = (double ***) malloc(nirreps * sizeof(double **));
  W3 = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      Gij = Gi ^ Gj;
      for(Gk=0; Gk < nirreps; Gk++) {
	Gijk = Gi ^ Gj ^ Gk;
	Gjk = Gj ^ Gk;

	for(Gab=0; Gab < nirreps; Gab++) {
	  /* This will need to change for non-totally-symmetric cases */
	  Gc = Gab ^ Gijk;
	  W1[Gab] = dpd_block_matrix(FBB.params->coltot[Gab], avirtpi[Gc]);
	}
	for(Ga=0; Ga < nirreps; Ga++) {
	  Gcb = Ga ^ Gijk;
	  W2[Ga] = dpd_block_matrix(bvirtpi[Ga], WMaFe.params->coltot[Gcb]);
	}
	for(Gb=0; Gb < nirreps; Gb++) {
	  Gac = Gb ^ Gijk; /* assumes totally symmetric */
	  W3[Gb] = dpd_block_matrix(bvirtpi[Gb], WmAfE.params->coltot[Gac]);
	}

	for(i=0; i < boccpi[Gi]; i++) {
	  I = bocc_off[Gi] + i;
	  for(j=0; j < boccpi[Gj]; j++) {
	    J = bocc_off[Gj] + j;
	    for(k=0; k < aoccpi[Gk]; k++) {
	      K = aocc_off[Gk] + k;

/*
	      T3_AAB(W1, nirreps, I, Gi, J, Gj, K, Gk, &T2BB, &T2BA, &T2AB, 
		     &FBB, &FBA, &FAB, &EBB, &EBA, &EAB, &fij, &fIJ, &fab, &fAB,
		     boccpi, bocc_off, aoccpi, aocc_off, bvirtpi, bvir_off, avirtpi, avir_off, 0.0);
*/

	      /* t_KC <-- 1/4 t_ijKabC <ij||ab> */

	      Gc = Gk;  /* assumes T1 is totally symmetric */
	      Gab = Gij; /* assumes <ij||ab> is totally symmetric */

	      ij = DBB.params->rowidx[I][J];

	      nrows = DBB.params->coltot[Gij];
	      ncols = avirtpi[Gc];

	      if(nrows && ncols)
		C_DGEMV('t', nrows, ncols, 0.25, W1[Gab][0], ncols, DBB.matrix[Gij][ij], 1,
			1.0, t1newA.matrix[Gk][k], 1);

	      /* t_ia <-- t_ijKabC <jK|bC> */

	      Ga = Gi;   /* assumes T1 is totally symmetric */
	      Gbc = Gjk; /* assumes <jk|bc> is totally symmetric */

	      jk = DBA.params->rowidx[J][K];

	      for(Gab=0; Gab < nirreps; Gab++) {
		Gb = Ga ^ Gab;
		Gc = Gb ^ Gbc;

		ab = W_offset[Gab][Ga];
		bc = DBA.col_offset[Gjk][Gb];

		nrows = bvirtpi[Ga];
		ncols = bvirtpi[Gb] * avirtpi[Gc];

		if(nrows && ncols)
		  C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][ab], ncols, &(DBA.matrix[Gjk][jk][bc]), 1,
			  1.0, t1newB.matrix[Gi][i], 1);
	      }

	      /* t_ijab <-- t_ijKabC F_KC */
	      Gc = Gk;    /* assumes Fme is totally symmetric */
	      Gab = Gij;  /* Assumes tIJAB is totally symmetric */

	      nrows = T2BBnew.params->coltot[Gij];
	      ncols = avirtpi[Gc];
	      ij = T2BBnew.params->rowidx[I][J];
 
	      if(nrows && ncols)
		C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][0], ncols, FME.matrix[Gk][k], 1,
			1.0, T2BBnew.matrix[Gij][ij], 1);

	      /* t_jKbC <-- t_ijKabC F_ia */
	      Ga = Gi;   /* assumes Fia is totally symmetric */
	      Gbc = Gjk; /* assumes t_jKbC is totally symmetric */

	      jk = T2BA.params->rowidx[J][K];

	      for(Gab=0; Gab < nirreps; Gab++) {
		Gb = Ga ^ Gab;
		Gc = Gb ^ Gbc;

		ab = W_offset[Gab][Ga];
		bc = T2BAnew.col_offset[Gjk][Gb];

		nrows = bvirtpi[Ga];
		ncols = bvirtpi[Gb] * avirtpi[Gc];

		if(nrows && ncols)
		  C_DGEMV('t', nrows, ncols, 1.0, W1[Gab][ab], ncols, Fme.matrix[Gi][i], 1,
			  1.0, &(T2BAnew.matrix[Gjk][jk][bc]), 1);
	      }

	      /* t_jida <-- t_ijKabC W_KdCb */
	      /* sort W1(ab,C) to W2(a,Cb) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;
		for(ab=0; ab < FBB.params->coltot[Gab]; ab++) {
		  A = FBB.params->colorb[Gab][ab][0];
		  B = FBB.params->colorb[Gab][ab][1];
		  Ga = FBB.params->rsym[A];
		  a = A - bvir_off[Ga];
		  for(c=0; c < avirtpi[Gc]; c++) {
		    C = avir_off[Gc] + c;
		    cb = WMaFe.params->colidx[C][B];
		    W2[Ga][a][cb] = W1[Gab][ab][c];
		  }
		}
	      }

	      ji = T2BBnew.params->rowidx[J][I];

	      for(Gd=0; Gd < nirreps; Gd++) {
		Gcb = Gkd = Gk ^ Gd; /* assumes totally symmetric */
		Ga = Gd ^ Gij;       /* assumes totally symmetric */

		kd = WMaFe.row_offset[Gkd][K];
		WMaFe.matrix[Gkd] = dpd_block_matrix(bvirtpi[Gd], WMaFe.params->coltot[Gkd]);
		dpd_buf4_mat_irrep_rd_block(&WMaFe, Gkd, kd, bvirtpi[Gd]);
		Z = block_matrix(bvirtpi[Ga], bvirtpi[Gd]);

		nrows = bvirtpi[Ga];
		ncols = bvirtpi[Gd];
		nlinks = WMaFe.params->coltot[Gkd];

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W2[Ga][0], nlinks, 
			  WMaFe.matrix[Gkd][0], nlinks, 0.0, Z[0], ncols);

		for(a=0; a < bvirtpi[Ga]; a++) {
		  A = bvir_off[Ga] + a;
		  for(d=0; d < bvirtpi[Gd]; d++) {
		    DD = bvir_off[Gd] + d;
		    ad = T2BBnew.params->colidx[A][DD];
		    da = T2BBnew.params->colidx[DD][A];
		    T2BBnew.matrix[Gij][ji][ad] += -Z[a][d];
		    T2BBnew.matrix[Gij][ji][da] += Z[a][d];
		  }
		}

		dpd_free_block(WMaFe.matrix[Gkd], bvirtpi[Gd], WMaFe.params->coltot[Gkd]);
		free_block(Z);
	      }

	      /* t_jKcD <-- 1/2 t_ijKabC W_idab */

	      jk = T2BAnew.params->rowidx[J][K];

	      for(Gd=0; Gd < nirreps; Gd++) {
		Gab = Gid = Gi ^ Gd; /* assumes totally symmetric */
		Gc = Gab ^ Gijk;     /* assumes totally symmetric */

		id = Wmafe.row_offset[Gid][I];
		Wmafe.matrix[Gid] = dpd_block_matrix(bvirtpi[Gd], Wmafe.params->coltot[Gid]);
		dpd_buf4_mat_irrep_rd_block(&Wmafe, Gid, id, bvirtpi[Gd]);
		Z = block_matrix(avirtpi[Gc], bvirtpi[Gd]);

		nrows = avirtpi[Gc];
		ncols = bvirtpi[Gd];
		nlinks = Wmafe.params->coltot[Gid];

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 't', nrows, ncols, nlinks, 0.5, W1[Gab][0], nrows,
			  Wmafe.matrix[Gid][0], nlinks, 0.0, Z[0], ncols);

		for(c=0; c < avirtpi[Gc]; c++) {
		  C = avir_off[Gc] + c;
		  for(d=0; d < bvirtpi[Gd]; d++) {
		    DD = bvir_off[Gd] + d;
		    dc = T2BAnew.params->colidx[DD][C];
		    T2BAnew.matrix[Gjk][jk][dc] += Z[c][d];
		  }
		}

		free_block(Z);
		dpd_free_block(Wmafe.matrix[Gid], bvirtpi[Gd], Wmafe.params->coltot[Gid]);
	      }

	      /* t_jKbD <-- t_ijKabC W_iDaC */
	      /* sort W1(ab,C) to W3(b,aC) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;
		for(ab=0; ab < FBB.params->coltot[Gab]; ab++) {
		  A = FBB.params->colorb[Gab][ab][0];
		  B = FBB.params->colorb[Gab][ab][1];
		  Gb = FBB.params->ssym[B];
		  b = B - bvir_off[Gb];
		  for(c=0; c < avirtpi[Gc]; c++) {
		    C = avir_off[Gc] + c;
		    ac = WmAfE.params->colidx[A][C];
		    W3[Gb][b][ac] = W1[Gab][ab][c];
		  }
		}
	      }

	      jk = T2BAnew.params->rowidx[J][K];

	      for(Gd=0; Gd < nirreps; Gd++) {
		Gac = Gid = Gi ^ Gd; /* assumes totally symmetric */
		Gb = Gac ^ Gijk;     /* assumes totally symmetric */

		id = WmAfE.row_offset[Gid][I]; 
		WmAfE.matrix[Gid] = dpd_block_matrix(avirtpi[Gd], WmAfE.params->coltot[Gid]);
		dpd_buf4_mat_irrep_rd_block(&WmAfE, Gid, id, avirtpi[Gd]);
		Z = block_matrix(bvirtpi[Gb], avirtpi[Gd]);

		nrows = bvirtpi[Gb];
		ncols = avirtpi[Gd];
		nlinks = WmAfE.params->coltot[Gid];

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W3[Gb][0], nlinks,
			  WmAfE.matrix[Gid][0], nlinks, 0.0, Z[0], ncols);

		for(b=0; b < bvirtpi[Gb]; b++) {
		  B = bvir_off[Gb] + b;
		  for(d=0; d < avirtpi[Gd]; d++) {
		    DD = avir_off[Gd] + d;
		    bd = T2BAnew.params->colidx[B][DD];
		    T2BAnew.matrix[Gjk][jk][bd] += Z[b][d];
		  }
		}

		dpd_free_block(WmAfE.matrix[Gid], avirtpi[Gd], WmAfE.params->coltot[Gid]);
		free_block(Z);
	      }

	      /* t_miab <--- +t_ijKabC W_jKmC */
	      /* t_imab <--- -t_ijKabC W_jKmC */

	      jk = WmNiE.params->rowidx[J][K];

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gab = Gmi = Gm ^ Gi;  /* assumes totally symmetric */
		Gc = Gab ^ Gijk;      /* assumes totally symmetric */

		mc = WmNiE.col_offset[Gjk][Gm];

		nrows = FBB.params->coltot[Gab];
		ncols = boccpi[Gm];
		nlinks = avirtpi[Gc];

		Z = dpd_block_matrix(nrows, ncols);

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W1[Gab][0], nlinks,
			  &(WmNiE.matrix[Gjk][jk][mc]), nlinks, 0.0, Z[0], ncols);

		for(m=0; m < ncols; m++) {
		  M = bocc_off[Gm] + m;
		  mi = T2BBnew.params->rowidx[M][I];
		  im = T2BBnew.params->rowidx[I][M];
		  for(ab=0; ab < nrows; ab++) {
		    T2BBnew.matrix[Gmi][mi][ab] += Z[ab][m];
		    T2BBnew.matrix[Gmi][im][ab] -= Z[ab][m];
		  }
		}

		dpd_free_block(Z, nrows, ncols);
	      }

	      /* t_mKbC <-- 1/2 t_ijKabC W_ijma */
	      /* sort W(ab,C) to W(a,bC) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;  /* assumes totally symmetric */
		for(ab=0; ab < FBB.params->coltot[Gab]; ab++ ){
		  A = FBB.params->colorb[Gab][ab][0];
		  B = FBB.params->colorb[Gab][ab][1];
		  Ga = FBB.params->rsym[A];
		  a = A - bvir_off[Ga];
		  for(c=0; c < avirtpi[Gc]; c++) {
		    C = avir_off[Gc] + c;
		    bc = T2BAnew.params->colidx[B][C];
		    W3[Ga][a][bc] = W1[Gab][ab][c];
		  }
		}
	      }

	      ij = Wmnie.params->rowidx[I][J];

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gbc = Gmk = Gm ^ Gk;  /* assumes totally symmetric */
		Ga = Gbc ^ Gijk;      /* assumes totally symmetric */

		ma = Wmnie.col_offset[Gij][Gm];

		nrows = T2BAnew.params->coltot[Gmk];
		ncols = boccpi[Gm];
		nlinks = bvirtpi[Ga];

		Z = dpd_block_matrix(nrows, ncols);

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 't', nrows, ncols, nlinks, 0.5, W3[Ga][0], nrows,
			  &(Wmnie.matrix[Gij][ij][ma]), nlinks, 0.0, Z[0], ncols);

		for(m=0; m < boccpi[Gm]; m++) {
		  M = bocc_off[Gm] + m;
		  mk = T2BAnew.params->rowidx[M][K];
		  for(bc=0; bc < nrows; bc++) {
		    T2BAnew.matrix[Gmk][mk][bc] += Z[bc][m];
		  }
		}

		dpd_free_block(Z, nrows, ncols);
	      }

	      /* t_iMbC <-- t_ijKabC W_KjMa */
	      /* sort W(ab,C) to W(a,bC) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;  /* assumes totally symmetric */
		for(ab=0; ab < FBB.params->coltot[Gab]; ab++ ){
		  A = FBB.params->colorb[Gab][ab][0];
		  B = FBB.params->colorb[Gab][ab][1];
		  Ga = FBB.params->rsym[A];
		  a = A - bvir_off[Ga];
		  for(c=0; c < avirtpi[Gc]; c++) {
		    C = avir_off[Gc] + c;
		    bc = T2BAnew.params->colidx[B][C];
		    W3[Ga][a][bc] = W1[Gab][ab][c];
		  }
		}
	      }

	      kj = WMnIe.params->rowidx[K][J];

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gbc = Gim = Gi ^ Gm;  /* assumes totally symmetric */
		Ga = Gbc ^ Gijk;      /* assumes totally symmetric */

		ma = WMnIe.col_offset[Gjk][Gm];

		nrows = T2BAnew.params->coltot[Gim];
		ncols = aoccpi[Gm];
		nlinks = bvirtpi[Ga];

		Z = dpd_block_matrix(nrows, ncols);

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, W3[Ga][0], nrows,
			  &(WMnIe.matrix[Gjk][kj][ma]), nlinks, 0.0, Z[0], ncols);

		for(m=0; m < aoccpi[Gm]; m++) {
		  M = aocc_off[Gm] + m;
		  im = T2BAnew.params->rowidx[I][M];
		  for(bc=0; bc < nrows; bc++) {
		    T2BAnew.matrix[Gim][im][bc] += Z[bc][m];
		  }
		}

		dpd_free_block(Z, nrows, ncols);
	      }

	    } /* k */
	  } /* j */
	} /* i */

	for(Gab=0; Gab < nirreps; Gab++) {
	  /* This will need to change for non-totally-symmetric cases */
	  Gc = Gab ^ Gijk;
	  dpd_free_block(W1[Gab], FBB.params->coltot[Gab], avirtpi[Gc]);
	}
	for(Ga=0; Ga < nirreps; Ga++) {
	  Gcb = Ga ^ Gijk; /* assumes totally symmetric */
	  dpd_free_block(W2[Ga], bvirtpi[Ga], WMaFe.params->coltot[Gcb]);
	}
	for(Gb=0; Gb < nirreps; Gb++) {
	  Gac = Gb ^ Gijk; /* assumes totally symmetric */
	  dpd_free_block(W3[Gb], bvirtpi[Gb], WmAfE.params->coltot[Gac]);
	}

      } /* Gk */
    } /* Gj */
  } /* Gi */

  free(W1);
  free(W2);
  free(W3);

  dpd_buf4_close(&EBB);
  dpd_buf4_close(&EAB);
  dpd_buf4_close(&EBA);
  dpd_buf4_close(&FBB);
  dpd_buf4_close(&FAB);
  dpd_buf4_close(&FBA);
  dpd_buf4_close(&T2BB);
  dpd_buf4_close(&T2AB);
  dpd_buf4_close(&T2BA);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);
  dpd_file2_close(&fij);
  dpd_file2_close(&fab);

  dpd_file2_mat_wrt(&t1newA);
  dpd_file2_init(&d1, CC_OEI, 0, 0, 1, "dIA");
  dpd_file2_dirprd(&d1, &t1newA);
  dpd_file2_close(&d1);
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "New tIA");
  dpd_file2_axpy(&t1newA, &t1, 1, 0);
  dpd_file2_close(&t1);
  dpd_file2_close(&t1newA);

  dpd_file2_mat_wrt(&t1newB);
  dpd_file2_init(&d1, CC_OEI, 0, 2, 3, "dia");
  dpd_file2_dirprd(&d1, &t1newB);
  dpd_file2_close(&d1);
  dpd_file2_init(&t1, CC_OEI, 0, 2, 3, "New tia");
  dpd_file2_axpy(&t1newB, &t1, 1, 0);
  dpd_file2_close(&t1);
  dpd_file2_close(&t1newB);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&DBB, h);
    dpd_buf4_mat_irrep_close(&DBA, h);
  }
  dpd_buf4_close(&DBB);
  dpd_buf4_close(&DBA);

  free_int_matrix(W_offset);

  dpd_file2_mat_close(&FME);
  dpd_file2_close(&FME);

  dpd_file2_mat_close(&Fme);
  dpd_file2_close(&Fme);

  dpd_buf4_close(&WMaFe);
  dpd_buf4_close(&Wmafe);
  dpd_buf4_close(&WmAfE);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&WmNiE, h);
    dpd_buf4_mat_irrep_close(&Wmnie, h);
    dpd_buf4_mat_irrep_close(&WMnIe, h);
  }
  dpd_buf4_close(&WmNiE);
  dpd_buf4_close(&Wmnie);
  dpd_buf4_close(&WMnIe);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&T2BBnew, h);
    dpd_buf4_mat_irrep_close(&T2BBnew, h);
  }
  dpd_buf4_init(&D2, CC_DENOM, 0, 10, 15, 11, 16, 0, "dijab");
  dpd_buf4_dirprd(&D2, &T2BBnew);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 15, 12, 17, 0, "New tijab");
  dpd_buf4_axpy(&T2BBnew, &T2, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2BBnew);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&T2BAnew, h);
    dpd_buf4_mat_irrep_close(&T2BAnew, h);
  }
  dpd_buf4_sort(&T2BAnew, CC_MISC, qpsr, 22, 28, "CC3 tIjAb");
  dpd_buf4_close(&T2BAnew);
  dpd_buf4_init(&T2BAnew, CC_MISC, 0, 22, 28, 22, 28, 0, "CC3 tIjAb");
  dpd_buf4_init(&D2, CC_DENOM, 0, 22, 28, 22, 28, 0, "dIjAb");
  dpd_buf4_dirprd(&D2, &T2BAnew);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
  dpd_buf4_axpy(&T2BAnew, &T2, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2BAnew);
}
}} // namespace psi::ccenergy
