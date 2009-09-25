/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void cc3_l3l2_RHF_AAA(void);
void cc3_l3l2_RHF_AAB(void);

void L3_AAA(double ***W1, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk, 
	    dpdbuf4 *T2, dpdbuf4 *F, dpdbuf4 *E, dpdfile2 *fIJ, dpdfile2 *fAB, 
	    dpdbuf4 *D, dpdbuf4 *LIJAB, dpdfile2 *LIA, dpdfile2 *FME,
	    int *occpi, int *occ_off, int *virtpi, int *vir_off);

void L3_AAB(double ***W1, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk, 
	    dpdbuf4 *T2AA, dpdbuf4 *T2AB, dpdbuf4 *T2BA, dpdbuf4 *FAA, dpdbuf4 *FAB, dpdbuf4 *FBA,
	    dpdbuf4 *EAA, dpdbuf4 *EAB, dpdbuf4 *EBA, dpdfile2 *fIJ, dpdfile2 *fij, 
	    dpdfile2 *fAB, dpdfile2 *fab, dpdbuf4 *DAA, dpdbuf4 *DAB, dpdbuf4 *LIJAB, dpdbuf4 *LIjAb,
	    dpdfile2 *LIA, dpdfile2 *Lia, dpdfile2 *FME, dpdfile2 *Fme,
	    int *aoccpi, int *aocc_off, int *boccpi, int *bocc_off,
	    int *avirtpi, int *avir_off, int *bvirtpi, int *bvir_off);

void cc3_l3l2(void)
{
  if(params.ref == 0) {
    cc3_l3l2_RHF_AAA();
    cc3_l3l2_RHF_AAB();
  }
}

void cc3_l3l2_RHF_AAA(void)
{
  int h, nirreps;
  int *occ_off, *occpi;
  int *vir_off, *virtpi;
  int Gijk;
  int Gi, Gj, Gk;
  int Ga, Gb, Gc;
  int Gab, ab;
  int i, j, k, I, J, K;
  int a, b, c, A, B, C;
  double ***W1, ***W2;
  dpdbuf4 L, E, F;
  dpdfile2 fIJ, fAB;
  dpdfile2 FME, LIA;
  dpdbuf4 Dints, LIJAB;
  dpdbuf4 WMAFE, WMNIE;
  dpdbuf4 ZIGDE, T2;
  dpdbuf4 ZDMAE;
  dpdbuf4 ZLMAO;
  dpdbuf4 ZIMLE;
  dpdbuf4 L2new, L2, D2;
  int Gjk, jk, Gid, id, Gik, ik;
  int Gd, d, DD;
  int cd, dc;
  int Gm, m, M;
  int Gmi, mi, im, mc;
  int Gjd, jd;
  int Gij, ij, Gmk, mk, am, Gbc, bc;
  int ac, mb;
  int nrows, ncols, nlinks;
  double **Z;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  occ_off = moinfo.occ_off;
  virtpi = moinfo.virtpi;
  vir_off = moinfo.vir_off;

  dpd_buf4_init(&WMAFE, CC3_HET1, 0, 10, 5, 10, 7, 0, "CC3 WABEI (IE,B>A)");
  dpd_buf4_init(&WMNIE, CC3_HET1, 0, 0, 10, 2, 10, 0, "CC3 WMBIJ (I>J,MB)");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&WMNIE, h);
    dpd_buf4_mat_irrep_rd(&WMNIE, h);
  }

  dpd_buf4_init(&L2new, CC3_MISC, 0, 0, 5, 0, 5, 0, "CC3 LIJAB");
  for(h=0; h < nirreps; h++) dpd_buf4_mat_irrep_init(&L2new, h);

  dpd_buf4_init(&ZIGDE, CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZIGDE");
  dpd_buf4_scm(&ZIGDE, 0.0); /* must be cleared in each iteration */

  dpd_buf4_init(&ZDMAE, CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZDMAE (MD,AE)");
  dpd_buf4_scm(&ZDMAE, 0.0);

  dpd_buf4_init(&ZLMAO, CC3_MISC, 0, 0, 11, 0, 11, 0, "CC3 ZLMAO");
  for(h=0; h < nirreps; h++) dpd_buf4_mat_irrep_init(&ZLMAO, h);

  dpd_buf4_init(&ZIMLE, CC3_MISC, 0, 0, 10, 0, 10, 0, "CC3 ZIMLE");
  for(h=0; h < nirreps; h++) dpd_buf4_mat_irrep_init(&ZIMLE, h);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2, h);
    dpd_buf4_mat_irrep_rd(&T2, h);
  }

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");

  dpd_buf4_init(&L, CC_LAMBDA, 0, 0, 5, 2, 7, 0, "LIJAB");
  dpd_buf4_init(&F, CC3_HET1, 0, 10, 5, 10, 7, 0, "CC3 WAMEF (MA,F>E)");
  dpd_buf4_init(&E, CC3_HET1, 0, 0, 10, 2, 10, 0, "CC3 WMNIE (M>N,IE)");

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
  dpd_buf4_init(&LIJAB, CC_LAMBDA, 0, 0, 5, 2, 7, 0, "LIJAB");
  dpd_file2_init(&LIA, CC_LAMBDA, 0, 0, 1, "LIA");
  dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");

  W1 = (double ***) malloc(nirreps * sizeof(double **));
  W2 = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      Gij = Gi ^ Gj;
      for(Gk=0; Gk < nirreps; Gk++) {
	Gijk = Gi ^ Gj ^ Gk;

	Gjk = Gj ^ Gk;
	Gik = Gi ^ Gk;

	for(Gab=0; Gab < nirreps; Gab++) {
	  Gc = Gab ^ Gijk; /* totally symmetric */
	  W1[Gab] = dpd_block_matrix(F.params->coltot[Gab], virtpi[Gc]);
	}
	for(Ga=0; Ga < nirreps; Ga++) {
	  Gbc = Ga ^ Gijk; /* totally symmetric */
	  W2[Ga] = dpd_block_matrix(virtpi[Ga], F.params->coltot[Gbc]);
	}

	for(i=0; i < occpi[Gi]; i++) {
	  I = occ_off[Gi] + i;
	  for(j=0; j < occpi[Gj]; j++) {
	    J = occ_off[Gj] + j;
	    for(k=0; k < occpi[Gk]; k++) {
	      K = occ_off[Gk] + k;

	      L3_AAA(W1, nirreps, I, Gi, J, Gj, K, Gk, &L, &F, &E, &fIJ, &fAB,
		     &Dints, &LIJAB, &LIA, &FME, occpi, occ_off, virtpi, vir_off);

	      /* L_JKDC <-- +1/2 t_IJKABC W_ABID (IDAB) */
	      /* L_JKCD <-- -1/2 t_IJKABC W_ABID (IDAB) */
	      jk = L2new.params->rowidx[J][K];
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
		    cd = L2new.params->colidx[C][DD];
		    dc = L2new.params->colidx[DD][C];
		    L2new.matrix[Gjk][jk][dc] += Z[c][d];
		    L2new.matrix[Gjk][jk][cd] += -Z[c][d];
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
		  mi = L2new.params->rowidx[M][I];
		  im = L2new.params->rowidx[I][M];
		  for(ab=0; ab < nrows; ab++) {
		    L2new.matrix[Gmi][mi][ab] += Z[ab][m];
		    L2new.matrix[Gmi][im][ab] -= Z[ab][m];
		  }
		}

		dpd_free_block(Z, nrows, ncols);
	      }

	      /* Z_IDAB <-- 1/2 L_IJKABC t_JKDC */

	      jk = T2.params->rowidx[J][K];
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gid = Gab; /* totally symmetric */
		Gc = Gab ^ Gijk; /* totally symmetric */
		Gd = Gi ^ Gid;

		nrows = virtpi[Gd];
		ncols = ZIGDE.params->coltot[Gid];
		nlinks = virtpi[Gc];

		dc = T2.col_offset[Gjk][Gd];
		id = ZIGDE.row_offset[Gid][I];
		ZIGDE.matrix[Gid] = dpd_block_matrix(nrows, ncols);
		dpd_buf4_mat_irrep_rd_block(&ZIGDE, Gid, id, nrows);

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 't', nrows, ncols, nlinks, 0.5, &(T2.matrix[Gjk][jk][dc]), nlinks,
			  W1[Gab][0], nlinks, 1.0, ZIGDE.matrix[Gid][0], ncols);

		dpd_buf4_mat_irrep_wrt_block(&ZIGDE, Gid, id, nrows);
		dpd_free_block(ZIGDE.matrix[Gid], nrows, ncols);
	      }

	      /* Z_JDAB <-- 1/2 L_IJKABC t_IKDC */
	      ik = T2.params->rowidx[I][K];
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gjd = Gab; /* totally symmetric */
		Gc = Gab ^ Gijk; /* totally symmetric */
		Gd = Gj ^ Gjd;

		nrows = virtpi[Gd];
		ncols = ZDMAE.params->coltot[Gjd];
		nlinks = virtpi[Gc];

		dc = T2.col_offset[Gik][Gd];
		jd = ZDMAE.row_offset[Gjd][J];
		ZDMAE.matrix[Gjd] = dpd_block_matrix(nrows, ncols);
		dpd_buf4_mat_irrep_rd_block(&ZDMAE, Gjd, jd, nrows);

		if(nrows && ncols && nlinks) 
		  C_DGEMM('n', 't', nrows, ncols, nlinks, 0.5, &(T2.matrix[Gik][ik][dc]), nlinks,
			  W1[Gab][0], nlinks, 1.0, ZDMAE.matrix[Gjd][0], ncols);

		dpd_buf4_mat_irrep_wrt_block(&ZDMAE, Gjd, jd, nrows);
		dpd_free_block(ZDMAE.matrix[Gjd], nrows, ncols);
	      }

	      /* Z_IJAM <-- -1/2 L_IJKABC t_MKBC */
	      /* sort W(AB,C) to W(A,BC) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;
		for(ab=0; ab < F.params->coltot[Gab]; ab++) {
		  A = F.params->colorb[Gab][ab][0];
		  B = F.params->colorb[Gab][ab][1];
		  Ga = F.params->rsym[A];
		  a = A - vir_off[Ga];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    bc = F.params->colidx[B][C];
		    W2[Ga][a][bc] = W1[Gab][ab][c];
		  }
		}
	      }

	      ij = ZLMAO.params->rowidx[I][J];

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gbc = Gmk = Gm ^ Gk; /* totally symmetric */
		Ga = Gij ^ Gm; /* totally symmetric */

		nrows = virtpi[Ga];
		ncols = T2.params->coltot[Gmk];

		for(m=0; m < occpi[Gm]; m++) {
		  M = occ_off[Gm] + m;
		  mk = T2.params->rowidx[M][K];
		  am = ZLMAO.col_offset[Gij][Ga] + m;

		  if(nrows && ncols)
		    C_DGEMV('n', nrows, ncols, -0.5, W2[Ga][0], ncols, T2.matrix[Gmk][mk], 1, 
			    1.0, &(ZLMAO.matrix[Gij][ij][am]), occpi[Gm]);
		}
	      }

	      /* Z_IJMB <-- -1/2 L_IJKABC t_MKAC */
	      /* sort W(AB,C) to W(B,AC) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;
		for(ab=0; ab < F.params->coltot[Gab]; ab++) {
		  A = F.params->colorb[Gab][ab][0];
		  B = F.params->colorb[Gab][ab][1];
		  Gb = F.params->ssym[B];
		  b = B - vir_off[Gb];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    ac = F.params->colidx[A][C];
		    W2[Gb][b][ac] = W1[Gab][ab][c];
		  }
		}
	      }

	      ij = ZIMLE.params->rowidx[I][J];

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gb = Gm ^ Gij; /* totally symmetric */
		Gmk = Gm ^ Gk;

		nrows = virtpi[Gb];
		ncols = T2.params->coltot[Gmk];

		for(m=0; m < occpi[Gm]; m++) {
		  M = occ_off[Gm] + m;
		  mk = T2.params->rowidx[M][K];
		  mb = ZIMLE.col_offset[Gij][Gm] + m * virtpi[Gb];

		  if(nrows && ncols)
		    C_DGEMV('n', nrows, ncols, -0.5, W2[Gb][0], ncols, T2.matrix[Gmk][mk], 1,
			    1.0, &(ZIMLE.matrix[Gij][ij][mb]), 1);
		}

	      }


	    } /* k */
	  } /* j */
	} /* i */

	for(Gab=0; Gab < nirreps; Gab++) {
	  Gc = Gab ^ Gijk; /* totally symmetric */
	  dpd_free_block(W1[Gab], F.params->coltot[Gab], virtpi[Gc]);
	}
	for(Ga=0; Ga < nirreps; Ga++) {
	  Gbc = Ga ^ Gijk; /* totally symmetric */
	  dpd_free_block(W2[Ga], virtpi[Ga], F.params->coltot[Gbc]);
	}

      } /* Gk */
    } /* Gj */
  } /* Gi */

  free(W1);
  free(W2);

  dpd_buf4_close(&E);
  dpd_buf4_close(&F);
  dpd_buf4_close(&L);

  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);

  dpd_file2_close(&FME);
  dpd_file2_close(&LIA);
  dpd_buf4_close(&Dints);
  dpd_buf4_close(&LIJAB);

  dpd_buf4_close(&WMAFE);
  for(h=0; h < nirreps; h++) dpd_buf4_mat_irrep_close(&WMNIE, h);
  dpd_buf4_close(&WMNIE);

  dpd_buf4_close(&ZIGDE);
  dpd_buf4_close(&ZDMAE);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&ZLMAO, h);
    dpd_buf4_mat_irrep_close(&ZLMAO, h);
  }
  dpd_buf4_close(&ZLMAO);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&ZIMLE, h);
    dpd_buf4_mat_irrep_close(&ZIMLE, h);
  }
  dpd_buf4_close(&ZIMLE);

  for(h=0; h < nirreps; h++) dpd_buf4_mat_irrep_close(&T2, h);
  dpd_buf4_close(&T2);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&L2new, h);
    dpd_buf4_mat_irrep_close(&L2new, h);
  }
  dpd_buf4_init(&D2, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
  dpd_buf4_dirprd(&D2, &L2new);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&L2, CC_LAMBDA, 0, 0, 5, 2, 7, 0, "New LIJAB");
  dpd_buf4_axpy(&L2new, &L2, 1);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&L2new);

}

void cc3_l3l2_RHF_AAB(void)
{
  int h, nirreps;
  int *occ_off, *occpi;
  int *vir_off, *virtpi;
  int Gi, Gj, Gk, Gijk;
  int Ga, Gb, Gc, Gab;
  int i, j, k, I, J, K;
  int a, b, c, A, B, C;
  int ab;
  double ***W1, ***W2;
  dpdbuf4 L2AA, L2AB, L2BA, EAA, EAB, EBA, FAA, FAB, FBA;
  dpdfile2 fIJ, fAB, fij, fab;
  dpdbuf4 DAAints, DABints, LIJAB, LIjAb;
  dpdfile2 LIA, Lia, FME, Fme;
  dpdbuf4 L2AAnew, L2ABnew, L2, D2;
  dpdbuf4 WmAfE, WMnIe, WMAFE, WMaFe, WMNIE, WmNiE;
  dpdbuf4 ZIGDE, T2AB, T2AA, ZIgDe;
  dpdbuf4 ZDMAE, ZDmAe, ZdMAe;
  dpdbuf4 ZLMAO, ZLmAo;
  dpdbuf4 ZIMLE, ZImLe, ZImlE;
  int nrows, ncols, nlinks;
  int Gcb, cb;
  int Gij, ij, Gji, ji, Gjk, jk, kj, Gkj;
  int Gd, d, DD, ad, da, Gkd, kd;
  int Gm, m, M, Gmi, mi, im, mc;
  int Gid, id, dc, cd;
  int Gac, ac, Gca, ca, bd, db;
  int Gbc, bc, Gmk, mk, km, Gim, ma, am;
  int Gik, ik, Gki, ki, Gjd, jd;
  int Gmj, mj, cm, Gjm, jm;
  int mb;
  double **Z;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  occ_off = moinfo.occ_off;
  virtpi = moinfo.virtpi;
  vir_off = moinfo.vir_off;

  dpd_buf4_init(&L2AAnew, CC3_MISC, 0, 0, 5, 0, 5, 0, "CC3 LIJAB");
  for(h=0; h < nirreps; h++) dpd_buf4_mat_irrep_init(&L2AAnew, h);

  dpd_buf4_init(&L2ABnew, CC3_MISC, 0, 0, 5, 0, 5, 0, "CC3 LIjAb");
  for(h=0; h < nirreps; h++) dpd_buf4_mat_irrep_init(&L2ABnew, h);

  dpd_buf4_init(&T2AB, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&T2AA, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2AB, h);
    dpd_buf4_mat_irrep_rd(&T2AB, h);

    dpd_buf4_mat_irrep_init(&T2AA, h);
    dpd_buf4_mat_irrep_rd(&T2AA, h);
  }

  dpd_buf4_init(&ZIGDE, CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZIGDE");
  dpd_buf4_scm(&ZIGDE, 0.0); /* this must be cleared in each iteration */
  dpd_buf4_init(&ZIgDe, CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZIgDe");
  dpd_buf4_scm(&ZIgDe, 0.0); /* this must be cleared in each iteration */

  dpd_buf4_init(&ZDMAE, CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZDMAE (MD,AE)");
  dpd_buf4_scm(&ZDMAE, 0.0); /* must be cleared in each iteration */
  dpd_buf4_init(&ZDmAe, CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZDmAe (mD,Ae)");
  dpd_buf4_scm(&ZDmAe, 0.0); /* must be cleared in each iteration */
  dpd_buf4_init(&ZdMAe, CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZdMAe (Md,Ae)");
  dpd_buf4_scm(&ZdMAe, 0.0); /* must be cleared in each iteration */

  dpd_buf4_init(&ZLMAO, CC3_MISC, 0, 0, 11, 0, 11, 0, "CC3 ZLMAO");
  dpd_buf4_init(&ZLmAo, CC3_MISC, 0, 0, 11, 0, 11, 0, "CC3 ZLmAo");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&ZLMAO, h);
    dpd_buf4_mat_irrep_rd(&ZLMAO, h);

    dpd_buf4_mat_irrep_init(&ZLmAo, h);
  }

  dpd_buf4_init(&ZIMLE, CC3_MISC, 0, 0, 10, 0, 10, 0, "CC3 ZIMLE");
  dpd_buf4_init(&ZImLe, CC3_MISC, 0, 0, 10, 0, 10, 0, "CC3 ZImLe");
  dpd_buf4_init(&ZImlE, CC3_MISC, 0, 0, 10, 0, 10, 0, "CC3 ZImlE");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&ZIMLE, h);
    dpd_buf4_mat_irrep_rd(&ZIMLE, h);

    dpd_buf4_mat_irrep_init(&ZImLe, h);
    dpd_buf4_mat_irrep_init(&ZImlE, h);
  }

  dpd_buf4_init(&WmAfE, CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (iE,bA)");
  dpd_buf4_init(&WMAFE, CC3_HET1, 0, 10, 5, 10, 7, 0, "CC3 WABEI (IE,B>A)");
  dpd_buf4_init(&WMaFe, CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WaBeI (Ie,Ba)");

  dpd_buf4_init(&WMnIe, CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMbIj (Ij,Mb)");
  dpd_buf4_init(&WMNIE, CC3_HET1, 0, 0, 10, 2, 10, 0, "CC3 WMBIJ (I>J,MB)");
  dpd_buf4_init(&WmNiE, CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WmBiJ (iJ,mB)");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&WMnIe, h);
    dpd_buf4_mat_irrep_rd(&WMnIe, h);
    dpd_buf4_mat_irrep_init(&WMNIE, h);
    dpd_buf4_mat_irrep_rd(&WMNIE, h);
    dpd_buf4_mat_irrep_init(&WmNiE, h);
    dpd_buf4_mat_irrep_rd(&WmNiE, h);
  }

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fij");
  dpd_file2_init(&fab, CC_OEI, 0, 1, 1, "fab");

  dpd_buf4_init(&L2AA, CC_LAMBDA, 0, 0, 5, 2, 7, 0, "LIJAB");
  dpd_buf4_init(&L2AB, CC_LAMBDA, 0, 0, 5, 0, 5, 0, "LIjAb");
  dpd_buf4_init(&L2BA, CC_LAMBDA, 0, 0, 5, 0, 5, 0, "LiJaB");
  dpd_buf4_init(&FAA, CC3_HET1, 0, 10, 5, 10, 7, 0, "CC3 WAMEF (MA,F>E)");
  dpd_buf4_init(&FAB, CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WaMeF (Ma,Fe)");
  dpd_buf4_init(&FBA, CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAmEf (mA,fE)");
  dpd_buf4_init(&EAA, CC3_HET1, 0, 0, 10, 2, 10, 0, "CC3 WMNIE (M>N,IE)");
  dpd_buf4_init(&EAB, CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMnIe (Mn,Ie)");
  dpd_buf4_init(&EBA, CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WmNiE (mN,iE)");

  dpd_buf4_init(&DAAints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
  dpd_buf4_init(&DABints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&LIJAB, CC_LAMBDA, 0, 0, 5, 2, 7, 0, "LIJAB");
  dpd_buf4_init(&LIjAb, CC_LAMBDA, 0, 0, 5, 0, 5, 0, "LIjAb");
  dpd_file2_init(&LIA, CC_LAMBDA, 0, 0, 1, "LIA");
  dpd_file2_init(&Lia, CC_LAMBDA, 0, 0, 1, "Lia");
  dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");

  W1 = (double ***) malloc(nirreps * sizeof(double **));
  W2 = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      Gij = Gji = Gi ^ Gj;
      for(Gk=0; Gk < nirreps; Gk++) {
	Gijk = Gi ^ Gj ^ Gk;

	Gjk = Gkj = Gj ^ Gk;
	Gik = Gki = Gi ^ Gk;

	for(Gab=0; Gab < nirreps; Gab++) {
	  Gc = Gab ^ Gijk; /* totally symmetric */
	  W1[Gab] = dpd_block_matrix(FAA.params->coltot[Gab], virtpi[Gc]);
	}
	for(Ga=0; Ga < nirreps; Ga++) {
	  Gcb = Ga ^ Gijk; /* assumes totally symmetric */
	  W2[Ga] = dpd_block_matrix(virtpi[Ga], WmAfE.params->coltot[Gcb]);  /* alpha-beta-alpha */
	}

	for(i=0; i < occpi[Gi]; i++) {
	  I = occ_off[Gi] + i;
	  for(j=0; j < occpi[Gj]; j++) {
	    J = occ_off[Gj] + j;
	    for(k=0; k < occpi[Gk]; k++) {
	      K = occ_off[Gk] + k;

	      L3_AAB(W1, nirreps, I, Gi, J, Gj, K, Gk, &L2AA, &L2AB, &L2BA, 
		     &FAA, &FAB, &FBA, &EAA, &EAB, &EBA, &fIJ, &fij, &fAB, &fab,
		     &DAAints, &DABints, &LIJAB, &LIjAb, &LIA, &Lia, &FME, &Fme,
		     occpi, occ_off, occpi, occ_off, virtpi, vir_off, virtpi, vir_off);

	      /* t_JIDA <-- t_IJkABc W_kDcB */
	      /* sort W1(AB,c) to W2(A,cB) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++) {
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Ga = FAA.params->rsym[A];
		  a = A - vir_off[Ga];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    cb = WmAfE.params->colidx[C][B];
		    W2[Ga][a][cb] = W1[Gab][ab][c];
		  }
		}
	      }

	      ji = L2AAnew.params->rowidx[J][I];

	      for(Gd=0; Gd < nirreps; Gd++) {
		Gcb = Gkd = Gk ^ Gd; /* assumes totally symmetric */
		Ga = Gd ^ Gij;       /* assumes totally symmetric */

		kd = WmAfE.row_offset[Gkd][K];
		WmAfE.matrix[Gkd] = dpd_block_matrix(virtpi[Gd], WmAfE.params->coltot[Gkd]);
		dpd_buf4_mat_irrep_rd_block(&WmAfE, Gkd, kd, virtpi[Gd]);
		Z = block_matrix(virtpi[Ga], virtpi[Gd]);

		nrows = virtpi[Ga];
		ncols = virtpi[Gd];
		nlinks = WmAfE.params->coltot[Gkd];

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W2[Ga][0], nlinks, 
			  WmAfE.matrix[Gkd][0], nlinks, 0.0, Z[0], ncols);

		for(a=0; a < virtpi[Ga]; a++) {
		  A = vir_off[Ga] + a;
		  for(d=0; d < virtpi[Gd]; d++) {
		    DD = vir_off[Gd] + d;
		    ad = L2AAnew.params->colidx[A][DD];
		    da = L2AAnew.params->colidx[DD][A];
		    L2AAnew.matrix[Gij][ji][ad] += -Z[a][d];
		    L2AAnew.matrix[Gij][ji][da] += Z[a][d];
		  }
		}

		dpd_free_block(WmAfE.matrix[Gkd], virtpi[Gd], WmAfE.params->coltot[Gkd]);
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
		ncols = occpi[Gm];
		nlinks = virtpi[Gc];

		Z = dpd_block_matrix(nrows, ncols);

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W1[Gab][0], nlinks,
			  &(WMnIe.matrix[Gjk][jk][mc]), nlinks, 0.0, Z[0], ncols);

		for(m=0; m < ncols; m++) {
		  M = occ_off[Gm] + m;
		  mi = L2AAnew.params->rowidx[M][I];
		  im = L2AAnew.params->rowidx[I][M];
		  for(ab=0; ab < nrows; ab++) {
		    L2AAnew.matrix[Gmi][mi][ab] += Z[ab][m];
		    L2AAnew.matrix[Gmi][im][ab] -= Z[ab][m];
		  }
		}

		dpd_free_block(Z, nrows, ncols);
	      }

	      /* t_JkDc <-- 1/2 t_IJkABc W_IDAB */
	      /* t_KjCd <-- 1/2 t_IJkABc W_IDAB */

	      jk = L2ABnew.params->rowidx[J][K];
	      kj = L2ABnew.params->rowidx[K][J];

	      for(Gd=0; Gd < nirreps; Gd++) {
		Gab = Gid = Gi ^ Gd; /* assumes totally symmetric */
		Gc = Gab ^ Gijk;     /* assumes totally symmetric */

		id = WMAFE.row_offset[Gid][I];
		WMAFE.matrix[Gid] = dpd_block_matrix(virtpi[Gd], WMAFE.params->coltot[Gid]);
		dpd_buf4_mat_irrep_rd_block(&WMAFE, Gid, id, virtpi[Gd]);
		Z = block_matrix(virtpi[Gc], virtpi[Gd]);

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
		    dc = L2ABnew.params->colidx[DD][C];
		    cd = L2ABnew.params->colidx[C][DD];
		    L2ABnew.matrix[Gjk][jk][dc] += Z[c][d];
		    L2ABnew.matrix[Gjk][kj][cd] += Z[c][d];
		  }
		}

		free_block(Z);
		dpd_free_block(WMAFE.matrix[Gid], virtpi[Gd], WMAFE.params->coltot[Gid]);
	      }

	      /* t_JkBd <-- t_IJkABc W_IdAc */
	      /* t_KjBd <-- t_IJkABc W_IdAc */
	      /* sort W1(AB,c) to W2(B,Ac) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++) {
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Gb = FAA.params->ssym[B];
		  b = B - vir_off[Gb];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    ac = WMaFe.params->colidx[A][C];
		    W2[Gb][b][ac] = W1[Gab][ab][c];
		  }
		}
	      }

	      jk = L2ABnew.params->rowidx[J][K];
	      kj = L2ABnew.params->rowidx[K][J];

	      for(Gd=0; Gd < nirreps; Gd++) {
		Gac = Gid = Gi ^ Gd; /* assumes totally symmetric */
		Gb = Gac ^ Gijk;     /* assumes totally symmetric */

		id = WMaFe.row_offset[Gid][I]; 
		WMaFe.matrix[Gid] = dpd_block_matrix(virtpi[Gd], WMaFe.params->coltot[Gid]);
		dpd_buf4_mat_irrep_rd_block(&WMaFe, Gid, id, virtpi[Gd]);
		Z = block_matrix(virtpi[Gb], virtpi[Gd]);

		nrows = virtpi[Gb];
		ncols = virtpi[Gd];
		nlinks = WMaFe.params->coltot[Gid];

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W2[Gb][0], nlinks,
			  WMaFe.matrix[Gid][0], nlinks, 0.0, Z[0], ncols);

		for(b=0; b < virtpi[Gb]; b++) {
		  B = vir_off[Gb] + b;
		  for(d=0; d < virtpi[Gd]; d++) {
		    DD = vir_off[Gd] + d;
		    bd = L2ABnew.params->colidx[B][DD];
		    db = L2ABnew.params->colidx[DD][B];
		    L2ABnew.matrix[Gjk][jk][bd] += Z[b][d];
		    L2ABnew.matrix[Gjk][kj][db] += Z[b][d];
		  }
		}

		dpd_free_block(WMaFe.matrix[Gid], virtpi[Gd], WMaFe.params->coltot[Gid]);
		free_block(Z);
	      }

	      /* t_MkBc <-- 1/2 t_IJkABc W_IJMA */
	      /* sort W(AB,c) to W(A,Bc) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;  /* assumes totally symmetric */
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++ ){
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Ga = FAA.params->rsym[A];
		  a = A - vir_off[Ga];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    bc = L2ABnew.params->colidx[B][C];
		    W2[Ga][a][bc] = W1[Gab][ab][c];
		  }
		}
	      }

	      ij = WMNIE.params->rowidx[I][J];

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gbc = Gmk = Gm ^ Gk;  /* assumes totally symmetric */
		Ga = Gbc ^ Gijk;      /* assumes totally symmetric */

		ma = WMNIE.col_offset[Gij][Gm];

		nrows = L2ABnew.params->coltot[Gmk];
		ncols = occpi[Gm];
		nlinks = virtpi[Ga];

		Z = dpd_block_matrix(nrows, ncols);

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 't', nrows, ncols, nlinks, 0.5, W2[Ga][0], nrows,
			  &(WMNIE.matrix[Gij][ij][ma]), nlinks, 0.0, Z[0], ncols);

		for(m=0; m < occpi[Gm]; m++) {
		  M = occ_off[Gm] + m;
		  mk = L2ABnew.params->rowidx[M][K];
		  km = L2ABnew.params->rowidx[K][M];
		  for(Gb=0; Gb < nirreps; Gb++) {
		    Gc = Gbc ^ Gb;
		    for(b=0; b < virtpi[Gb]; b++) {
		      B = vir_off[Gb] + b;
		      for(c=0; c < virtpi[Gc]; c++) {
			C = vir_off[Gc] + c;
			bc = L2ABnew.params->colidx[B][C];
			cb = L2ABnew.params->colidx[C][B];
			L2ABnew.matrix[Gmk][mk][bc] += Z[bc][m];
			L2ABnew.matrix[Gmk][km][cb] += Z[bc][m];
		      }
		    }
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
		  a = A - vir_off[Ga];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    bc = L2ABnew.params->colidx[B][C];
		    W2[Ga][a][bc] = W1[Gab][ab][c];
		  }
		}
	      }

	      kj = WmNiE.params->rowidx[K][J];

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gbc = Gim = Gi ^ Gm;  /* assumes totally symmetric */
		Ga = Gbc ^ Gijk;      /* assumes totally symmetric */

		ma = WmNiE.col_offset[Gjk][Gm];

		nrows = L2ABnew.params->coltot[Gim];
		ncols = occpi[Gm];
		nlinks = virtpi[Ga];

		Z = dpd_block_matrix(nrows, ncols);

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, W2[Ga][0], nrows,
			  &(WmNiE.matrix[Gjk][kj][ma]), nlinks, 0.0, Z[0], ncols);

		for(m=0; m < occpi[Gm]; m++) {
		  M = occ_off[Gm] + m;
		  im = L2ABnew.params->rowidx[I][M];
		  mi = L2ABnew.params->rowidx[M][I];
		  for(Gb=0; Gb < nirreps; Gb++) {
		    Gc = Gbc ^ Gb;
		    for(b=0; b < virtpi[Gb]; b++) {
		      B = vir_off[Gb] + b;
		      for(c=0; c < virtpi[Gc]; c++) {
			C = vir_off[Gc] + c;
			bc = L2ABnew.params->colidx[B][C];
			cb = L2ABnew.params->colidx[C][B];
			L2ABnew.matrix[Gim][im][bc] += Z[bc][m];
			L2ABnew.matrix[Gim][mi][cb] += Z[bc][m];
		      }
		    }
		  }
		}

		dpd_free_block(Z, nrows, ncols);
	      }

	      /* Z_IDAB <-- L_IJkABc t_JkDc */

	      jk = T2AB.params->rowidx[J][K];
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gid = Gab; /* totally symmetric */
		Gc = Gab ^ Gijk; /* totally symmetric */
		Gd = Gi ^ Gid;

		nrows = virtpi[Gd];
		ncols = ZIGDE.params->coltot[Gid];
		nlinks = virtpi[Gc];

		dc = T2AB.col_offset[Gjk][Gd];
		id = ZIGDE.row_offset[Gid][I];
		ZIGDE.matrix[Gid] = dpd_block_matrix(nrows, ncols);

		if(nrows && ncols && nlinks) {
		  dpd_buf4_mat_irrep_rd_block(&ZIGDE, Gid, id, nrows);

		  C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, &(T2AB.matrix[Gjk][jk][dc]), nlinks,
			  W1[Gab][0], nlinks, 1.0, ZIGDE.matrix[Gid][0], ncols);

		  dpd_buf4_mat_irrep_wrt_block(&ZIGDE, Gid, id, nrows);
		}

		dpd_free_block(ZIGDE.matrix[Gid], nrows, ncols);
	      }

	      /* ZkDCa <-- 1/2 L_ijKabC t_ijdb */

	      ij = T2AA.params->rowidx[I][J];

	      /* sort W(ab,C) to W(b,Ca) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;  /* assumes totally symmetric */
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++ ){
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Gb = FAA.params->ssym[B];
		  b = B - vir_off[Gb];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    ca = ZIgDe.params->colidx[C][A];
		    W2[Gb][b][ca] = W1[Gab][ab][c];
		  }
		}
	      }

	      for(Gb=0; Gb < nirreps; Gb++) {
		Gd = Gb ^ Gij; /* totally symmetric */
		Gca = Gkd = Gk ^ Gd; /* totally symmetric */

		nrows = virtpi[Gd];
		ncols = ZIgDe.params->coltot[Gkd];
		nlinks = virtpi[Gb];

		db = T2AA.col_offset[Gij][Gd];
		kd = ZIgDe.row_offset[Gkd][K];
		ZIgDe.matrix[Gkd] = dpd_block_matrix(nrows, ncols);
		dpd_buf4_mat_irrep_rd_block(&ZIgDe, Gkd, kd, nrows);

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 'n', nrows, ncols, nlinks, 0.5, &(T2AA.matrix[Gij][ij][db]), nlinks,
			  W2[Gb][0], ncols, 1.0, ZIgDe.matrix[Gkd][0], ncols);

		dpd_buf4_mat_irrep_wrt_block(&ZIgDe, Gkd, kd, nrows);
		dpd_free_block(ZIgDe.matrix[Gkd], nrows, ncols);
	      }

	      /* Z_IdAc <-- L_IJkABc t_JkBd */

	      jk = T2AB.params->rowidx[J][K];

	      /* sort W(AB,c) to W(B,Ac) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;  /* assumes totally symmetric */
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++ ){
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Gb = FAA.params->ssym[B];
		  b = B - vir_off[Gb];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    ac = ZIgDe.params->colidx[A][C];
		    W2[Gb][b][ac] = W1[Gab][ab][c];
		  }
		}
	      }

	      for(Gb=0; Gb < nirreps; Gb++) {
		Gd = Gb ^ Gjk; /* totally symmetric */
		Gac = Gid = Gi ^ Gd; /* totally symmetric */

		nrows = virtpi[Gd]; 
		ncols = ZIgDe.params->coltot[Gid];
		nlinks = virtpi[Gb];

		bd = T2AB.col_offset[Gjk][Gb];
		id = ZIgDe.row_offset[Gid][I];
		ZIgDe.matrix[Gid] = dpd_block_matrix(nrows, ncols);
		dpd_buf4_mat_irrep_rd_block(&ZIgDe, Gid, id, nrows);

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, &(T2AB.matrix[Gjk][jk][bd]), nrows,
			  W2[Gb][0], ncols, 1.0, ZIgDe.matrix[Gid][0], ncols);

		dpd_buf4_mat_irrep_wrt_block(&ZIgDe, Gid, id, nrows);
		dpd_free_block(ZIgDe.matrix[Gid], nrows, ncols);
	      }

	      /* Z_JDAB <-- 1/2 L_IJkABc t_IkDc */
	      ik = T2AB.params->rowidx[I][K];
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gjd = Gab; /* totally symmetric */
		Gc = Gab ^ Gijk; /* totally symmetric */
		Gd = Gj ^ Gjd;

		nrows = virtpi[Gd];
		ncols = ZDMAE.params->coltot[Gjd];
		nlinks = virtpi[Gc];

		dc = T2AB.col_offset[Gik][Gd];
		jd = ZDMAE.row_offset[Gjd][J];
		ZDMAE.matrix[Gjd] = dpd_block_matrix(nrows, ncols);
		dpd_buf4_mat_irrep_rd_block(&ZDMAE, Gjd, jd, nrows);

		if(nrows && ncols && nlinks) 
		  C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, &(T2AB.matrix[Gik][ik][dc]), nlinks,
			  W1[Gab][0], nlinks, 1.0, ZDMAE.matrix[Gjd][0], ncols);

		dpd_buf4_mat_irrep_wrt_block(&ZDMAE, Gjd, jd, nrows);
		dpd_free_block(ZDMAE.matrix[Gjd], nrows, ncols);
	      }

	      /* Z_kDAc <-- 1/2 L_IJkABc t_IJDB */
	      ij = T2AA.params->rowidx[I][J];
 	      /* sort W(AB,c) to W(B,Ac) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;  /* assumes totally symmetric */
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++ ){
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Gb = FAA.params->ssym[B];
		  b = B - vir_off[Gb];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    ac = ZDmAe.params->colidx[A][C];
		    W2[Gb][b][ac] = W1[Gab][ab][c];
		  }
		}
	      }

	      for(Gb=0; Gb < nirreps; Gb++) {
		Gd = Gb ^ Gij; /* totally symmetric */
		Gac = Gkd = Gk ^ Gd; /* totally symmetric */

		nrows = virtpi[Gd]; 
		ncols = ZDmAe.params->coltot[Gkd];
		nlinks = virtpi[Gb];

		db = T2AA.col_offset[Gij][Gd];
		kd = ZDmAe.row_offset[Gkd][K];
		ZDmAe.matrix[Gkd] = dpd_block_matrix(nrows, ncols);
		dpd_buf4_mat_irrep_rd_block(&ZDmAe, Gkd, kd, nrows);

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 'n', nrows, ncols, nlinks, 0.5, &(T2AA.matrix[Gij][ij][db]), nlinks,
			  W2[Gb][0], ncols, 1.0, ZDmAe.matrix[Gkd][0], ncols);

		dpd_buf4_mat_irrep_wrt_block(&ZDmAe, Gkd, kd, nrows);
		dpd_free_block(ZDmAe.matrix[Gkd], nrows, ncols);
	      }

	      /* Z_iDCa <-- L_ijKabC t_KjDb */
	      kj = T2AB.params->rowidx[K][J];
 	      /* sort W(AB,c) to W(B,Ca) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;  /* assumes totally symmetric */
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++ ){
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Gb = FAA.params->ssym[B];
		  b = B - vir_off[Gb];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    ca = ZDmAe.params->colidx[C][A];
		    W2[Gb][b][ca] = W1[Gab][ab][c];
		  }
		}
	      }

	      for(Gb=0; Gb < nirreps; Gb++) {
		Gd = Gb ^ Gkj; /* totally symmetric */
		Gca = Gid = Gi ^ Gd; /* totally symmetric */

		nrows = virtpi[Gd]; 
		ncols = ZDmAe.params->coltot[Gid];
		nlinks = virtpi[Gb];

		db = T2AB.col_offset[Gkj][Gd];
		id = ZDmAe.row_offset[Gid][I];
		ZDmAe.matrix[Gid] = dpd_block_matrix(nrows, ncols);
		dpd_buf4_mat_irrep_rd_block(&ZDmAe, Gid, id, nrows);

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0, &(T2AB.matrix[Gkj][kj][db]), nlinks,
			  W2[Gb][0], ncols, 1.0, ZDmAe.matrix[Gid][0], ncols);

		dpd_buf4_mat_irrep_wrt_block(&ZDmAe, Gid, id, nrows);
		dpd_free_block(ZDmAe.matrix[Gid], nrows, ncols);
	      }

	      /* Z_KdCa <-- -1/2 L_ijKabC t_ijdb */
	      ij = T2AA.params->rowidx[I][J];
 	      /* sort W(AB,c) to W(B,Ca) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;  /* assumes totally symmetric */
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++ ){
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Gb = FAA.params->ssym[B];
		  b = B - vir_off[Gb];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    ca = ZdMAe.params->colidx[C][A];
		    W2[Gb][b][ca] = W1[Gab][ab][c];
		  }
		}
	      }

	      for(Gb=0; Gb < nirreps; Gb++) {
		Gd = Gb ^ Gij; /* totally symmetric */
		Gca = Gkd = Gk ^ Gd; /* totally symmetric */

		nrows = virtpi[Gd]; 
		ncols = ZdMAe.params->coltot[Gkd];
		nlinks = virtpi[Gb];

		db = T2AA.col_offset[Gij][Gd];
		kd = ZdMAe.row_offset[Gkd][K];
		ZdMAe.matrix[Gkd] = dpd_block_matrix(nrows, ncols);
		dpd_buf4_mat_irrep_rd_block(&ZdMAe, Gkd, kd, nrows);

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 'n', nrows, ncols, nlinks, -0.5, &(T2AA.matrix[Gij][ij][db]), nlinks,
			  W2[Gb][0], ncols, 1.0, ZdMAe.matrix[Gkd][0], ncols);

		dpd_buf4_mat_irrep_wrt_block(&ZdMAe, Gkd, kd, nrows);
		dpd_free_block(ZdMAe.matrix[Gkd], nrows, ncols);
	      }

	      /* Z_JdAc <-- L_IJkABc t_IkBd */
	      ik = T2AB.params->rowidx[I][K];

 	      /* sort W(AB,c) to W(B,Ca) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;  /* assumes totally symmetric */
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++ ){
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Gb = FAA.params->ssym[B];
		  b = B - vir_off[Gb];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    ac = ZdMAe.params->colidx[A][C];
		    W2[Gb][b][ac] = W1[Gab][ab][c];
		  }
		}
	      }

	      for(Gb=0; Gb < nirreps; Gb++) {
		Gd = Gb ^ Gik; /* totally symmetric */
		Gca = Gjd = Gj ^ Gd; /* totally symmetric */

		nrows = virtpi[Gd]; 
		ncols = ZdMAe.params->coltot[Gjd];
		nlinks = virtpi[Gb];

		bd = T2AB.col_offset[Gik][Gb];
		jd = ZdMAe.row_offset[Gjd][J];
		ZdMAe.matrix[Gjd] = dpd_block_matrix(nrows, ncols);
		dpd_buf4_mat_irrep_rd_block(&ZdMAe, Gjd, jd, nrows);

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, &(T2AB.matrix[Gik][ik][bd]), nrows,
			  W2[Gb][0], ncols, 1.0, ZdMAe.matrix[Gjd][0], ncols);

		dpd_buf4_mat_irrep_wrt_block(&ZdMAe, Gjd, jd, nrows);
		dpd_free_block(ZdMAe.matrix[Gjd], nrows, ncols);
	      }

	      /* Z_IJAM <-- -1/2 L_IJkABc t_MkBc */
	      /* sort W(AB,C) to W(A,BC) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++) {
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Ga = FAA.params->rsym[A];
		  a = A - vir_off[Ga];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    bc = FAA.params->colidx[B][C];
		    W2[Ga][a][bc] = W1[Gab][ab][c];
		  }
		}
	      }

	      ij = ZLMAO.params->rowidx[I][J];

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gbc = Gmk = Gm ^ Gk; /* totally symmetric */
		Ga = Gij ^ Gm; /* totally symmetric */

		nrows = virtpi[Ga];
		ncols = T2AB.params->coltot[Gmk];

		for(m=0; m < occpi[Gm]; m++) {
		  M = occ_off[Gm] + m;
		  mk = T2AB.params->rowidx[M][K];
		  am = ZLMAO.col_offset[Gij][Ga] + m;

		  if(nrows && ncols)
		    C_DGEMV('n', nrows, ncols, -1.0, W2[Ga][0], ncols, T2AB.matrix[Gmk][mk], 1, 
			    1.0, &(ZLMAO.matrix[Gij][ij][am]), occpi[Gm]);
		}
	      }

	      /* Z_KiCm <-- -1/2 L_ijKabC t_mjab */
	      ki = ZLmAo.params->rowidx[K][I];

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gab = Gmj = Gm ^ Gj; /* totally symmetric */
		Gc = Gm ^ Gki; /* totally symmetric */

		nrows = T2AA.params->coltot[Gmj];
		ncols = virtpi[Gc];

		for(m=0; m < occpi[Gm]; m++) {
		  M = occ_off[Gm] + m;
		  mj = T2AA.params->rowidx[M][J];
		  cm = ZLmAo.col_offset[Gki][Gc] + m;

		  if(nrows && ncols)
		    C_DGEMV('t', nrows, ncols, -0.5, W1[Gab][0], ncols, T2AA.matrix[Gmj][mj], 1,
			    1.0, &(ZLmAo.matrix[Gki][ki][cm]), occpi[Gm]);
		}
	      }

	      /* Z_IkAm <-- - L_IJkABc t_mJcB */
	      /* sort W(AB,C) to W(A,CB) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++) {
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Ga = FAA.params->rsym[A];
		  a = A - vir_off[Ga];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    cb = FAA.params->colidx[C][B];
		    W2[Ga][a][cb] = W1[Gab][ab][c];
		  }
		}
	      }

	      ik = ZLmAo.params->rowidx[I][K];

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gbc = Gmj = Gm ^ Gj; /* totally symmetric */
		Ga = Gm ^ Gik; /* totally symmetric */

		nrows = virtpi[Ga];
		ncols = T2AB.params->coltot[Gmj];

		for(m=0; m < occpi[Gm]; m++) {
		  M = occ_off[Gm] + m;
		  mj = T2AB.params->rowidx[M][J];
		  am = ZLmAo.col_offset[Gik][Ga] + m;

		  if(nrows && ncols)
		    C_DGEMV('n', nrows, ncols, -1.0, W2[Ga][0], ncols, T2AB.matrix[Gmj][mj], 1,
			    1.0, &(ZLmAo.matrix[Gik][ik][am]), occpi[Gm]);
		}
	      }

	      /* Z_IJMB <-- - L_IJkABc t_MkAc */
	      /* sort W(AB,C) to W(B,AC) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++) {
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Gb = FAA.params->ssym[B];
		  b = B - vir_off[Gb];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    ac = FAA.params->colidx[A][C];
		    W2[Gb][b][ac] = W1[Gab][ab][c];
		  }
		}
	      }

	      ij = ZIMLE.params->rowidx[I][J];

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gb = Gm ^ Gij; /* totally symmetric */
		Gmk = Gm ^ Gk;

		nrows = virtpi[Gb];
		ncols = T2AB.params->coltot[Gmk];

		for(m=0; m < occpi[Gm]; m++) {
		  M = occ_off[Gm] + m;
		  mk = T2AB.params->rowidx[M][K];
		  mb = ZIMLE.col_offset[Gij][Gm] + m * virtpi[Gb];

		  if(nrows && ncols)
		    C_DGEMV('n', nrows, ncols, -1.0, W2[Gb][0], ncols, T2AB.matrix[Gmk][mk], 1,
			    1.0, &(ZIMLE.matrix[Gij][ij][mb]), 1);
		}

	      }

	      /* Z_IkMc <-- -1/2 L_IJkABc t_MJAB */ 
	      ik = ZImLe.params->rowidx[I][K];
	      for(Gm=0; Gm < nirreps; Gm++) {
		Gc = Gm ^ Gik ; /* totally symmetric */
		Gab = Gmj = Gm ^ Gj; /* totally symmetric */

		nrows = T2AA.params->coltot[Gmj];
		ncols = virtpi[Gc];

		for(m=0; m < occpi[Gm]; m++) {
		  M = occ_off[Gm] + m;
		  mj = T2AA.params->rowidx[M][J];
		  mc = ZImLe.col_offset[Gik][Gm] + m * virtpi[Gc];

		  if(nrows && ncols)
		    C_DGEMV('t', nrows, ncols, -0.5, W1[Gab][0], ncols, T2AA.matrix[Gmj][mj], 1,
			    1.0, &(ZImLe.matrix[Gik][ik][mc]), 1);
		}
	      }

	      /* Z_KiMa <-- - L_ijKabC t_MjCb */
	      /* sort W(AB,C) to W(A,CB) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++) {
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Ga = FAA.params->rsym[A];
		  a = A - vir_off[Ga];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    cb = FAA.params->colidx[C][B];
		    W2[Ga][a][cb] = W1[Gab][ab][c];
		  }
		}
	      }

	      ki = ZImLe.params->rowidx[K][I];
	      for(Gm=0; Gm < nirreps; Gm++) {
		Ga = Gm ^ Gki; /* totally symmetric */
		Gmj = Gm ^ Gj;

		nrows = virtpi[Ga];
		ncols = T2AB.params->coltot[Gmj];

		for(m=0; m < occpi[Gm]; m++) {
		  M = occ_off[Gm] + m;
		  mj = T2AB.params->rowidx[M][J];
		  ma = ZImLe.col_offset[Gki][Gm] + m * virtpi[Ga];

		  if(nrows && ncols)
		    C_DGEMV('n', nrows, ncols, -1.0, W2[Ga][0], ncols, T2AB.matrix[Gmj][mj], 1,
			    1.0, &(ZImLe.matrix[Gki][ki][ma]), 1);
		}
	      }

	      /* Z_KimC <-- 1/2 L_ijKabC t_mjab */
	      ki = ZImlE.params->rowidx[K][I];
	      for(Gm=0; Gm < nirreps; Gm++) {
		Gc = Gm ^ Gki; /* totally symmetric */
		Gab = Gmj = Gm ^ Gj;

		nrows = T2AA.params->coltot[Gmj];
		ncols = virtpi[Gc];

		for(m=0; m < occpi[Gm]; m++) {
		  M = occ_off[Gm] + m;
		  mj = T2AA.params->rowidx[M][J];
		  mc = ZImlE.col_offset[Gki][Gm] + m * virtpi[Gc];

		  if(nrows && ncols)
		    C_DGEMV('t', nrows, ncols, 0.5, W1[Gab][0], ncols, T2AA.matrix[Gmj][mj], 1,
			    1.0, &(ZImlE.matrix[Gki][ki][mc]), 1);
		}
	      }

	      /* Z_IkmB <-- - l_IJkABc t_JmAc */
	      ik = ZImlE.params->rowidx[I][K];
	      /* sort W(AB,C) to W(B,AC) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++) {
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Gb = FAA.params->ssym[B];
		  b = B - vir_off[Gb];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    ac = FAA.params->colidx[A][C];
		    W2[Gb][b][ac] = W1[Gab][ab][c];
		  }
		}
	      }

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gb = Gm ^ Gik; /* totally symmetric */
		Gjm = Gm ^ Gj;

		nrows = virtpi[Gb];
		ncols = T2AB.params->coltot[Gjm];

		for(m=0; m < occpi[Gm]; m++) {
		  M = occ_off[Gm] + m;
		  jm = T2AB.params->rowidx[J][M];
		  mb = ZImlE.col_offset[Gki][Gm] + m * virtpi[Gb];

		  if(nrows && ncols)
		    C_DGEMV('n', nrows, ncols, -1.0, W2[Gb][0], ncols, T2AB.matrix[Gjm][jm], 1,
			    1.0, &(ZImlE.matrix[Gik][ik][mb]), 1);
		}
	      }

	    } /* k */
	  } /* j */
	} /* i */

	for(Gab=0; Gab < nirreps; Gab++) {
	  Gc = Gab ^ Gijk; /* totally symmetric */
	  dpd_free_block(W1[Gab], FAA.params->coltot[Gab], virtpi[Gc]);
	}
	for(Ga=0; Ga < nirreps; Ga++) {
	  Gcb = Ga ^ Gijk; /* assumes totally symmetric */
	  dpd_free_block(W2[Ga], virtpi[Ga], WmAfE.params->coltot[Gcb]);
	}

      } /* Gk */
    } /* Gj */
  } /* Gi */

  free(W1);
  free(W2);

  dpd_buf4_close(&EAA);
  dpd_buf4_close(&EAB);
  dpd_buf4_close(&EBA);
  dpd_buf4_close(&FAA);
  dpd_buf4_close(&FAB);
  dpd_buf4_close(&FBA);
  dpd_buf4_close(&L2AA);
  dpd_buf4_close(&L2AB);
  dpd_buf4_close(&L2BA);

  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);
  dpd_file2_close(&fij);
  dpd_file2_close(&fab);

  dpd_file2_close(&FME);
  dpd_file2_close(&Fme);
  dpd_file2_close(&LIA);
  dpd_file2_close(&Lia);

  dpd_buf4_close(&DAAints);
  dpd_buf4_close(&DABints);
  dpd_buf4_close(&LIJAB);
  dpd_buf4_close(&LIjAb);

  dpd_buf4_close(&WmAfE);
  dpd_buf4_close(&WMAFE);
  dpd_buf4_close(&WMaFe);

  for(h=0; h < nirreps; h++) dpd_buf4_mat_irrep_close(&WMnIe, h);
  for(h=0; h < nirreps; h++) dpd_buf4_mat_irrep_close(&WMNIE, h);
  for(h=0; h < nirreps; h++) dpd_buf4_mat_irrep_close(&WmNiE, h);
  dpd_buf4_close(&WMnIe);
  dpd_buf4_close(&WMNIE);
  dpd_buf4_close(&WmNiE);

  dpd_buf4_close(&ZIgDe);
  dpd_buf4_close(&ZIGDE);

  dpd_buf4_close(&ZDMAE);
  dpd_buf4_close(&ZDmAe);
  dpd_buf4_close(&ZdMAe);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&ZLMAO, h);
    dpd_buf4_mat_irrep_close(&ZLMAO, h);
  }
  dpd_buf4_close(&ZLMAO);
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&ZLmAo, h);
    dpd_buf4_mat_irrep_close(&ZLmAo, h);
  }
  dpd_buf4_close(&ZLmAo);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&ZIMLE, h);
    dpd_buf4_mat_irrep_close(&ZIMLE, h);
  }
  dpd_buf4_close(&ZIMLE);
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&ZImLe, h);
    dpd_buf4_mat_irrep_close(&ZImLe, h);
  }
  dpd_buf4_close(&ZImLe);
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&ZImlE, h);
    dpd_buf4_mat_irrep_close(&ZImlE, h);
  }
  dpd_buf4_close(&ZImlE);

  for(h=0; h < nirreps; h++) dpd_buf4_mat_irrep_close(&T2AB, h);
  for(h=0; h < nirreps; h++) dpd_buf4_mat_irrep_close(&T2AA, h);
  dpd_buf4_close(&T2AB);
  dpd_buf4_close(&T2AA);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&L2AAnew, h);
    dpd_buf4_mat_irrep_close(&L2AAnew, h);
  }
  dpd_buf4_init(&D2, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
  dpd_buf4_dirprd(&D2, &L2AAnew);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&L2, CC_LAMBDA, 0, 0, 5, 2, 7, 0, "New LIJAB");
  dpd_buf4_axpy(&L2AAnew, &L2, 1);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&L2AAnew);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&L2ABnew, h);
    dpd_buf4_mat_irrep_close(&L2ABnew, h);
  }
  dpd_buf4_init(&D2, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
  dpd_buf4_dirprd(&D2, &L2ABnew);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&L2, CC_LAMBDA, 0, 0, 5, 0, 5, 0, "New LIjAb");
  dpd_buf4_axpy(&L2ABnew, &L2, 1);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&L2ABnew);

  /* Spin adaptation will remove this.  And yes, this means that all the above
     calculations for LIJAB were pointless... -TDC */
  dpd_buf4_init(&L2, CC_LAMBDA, 0, 2, 7, 0, 5, 1, "New LIjAb");
  dpd_buf4_copy(&L2, CC_LAMBDA, "New LIJAB");
  dpd_buf4_copy(&L2, CC_LAMBDA, "New Lijab");
  dpd_buf4_close(&L2);

}

}} // namespace psi::cclambda
