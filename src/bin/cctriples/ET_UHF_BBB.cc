/*! \file
    \ingroup CCTRIPLES
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

namespace psi { namespace cctriples {

double ET_UHF_BBB(void)
{
  int cnt;
  int h, nirreps;
  int Gi, Gj, Gk, Ga, Gb, Gc, Gd, Gl;
  int Gji, Gij, Gjk, Gkj, Gik, Gki, Gijk;
  int Gab, Gbc, Gac;
  int Gid, Gjd, Gkd;
  int Gil, Gjl, Gkl;
  int I, J, K, A, B, C;
  int i, j, k, a, b, c;
  int ij, ji, ik, ki, jk, kj;
  int ab, ba, ac, ca, bc, cb;
  int cd, ad, bd;
  int id, jd, kd;
  int il, jl, kl;
  int lc, la, lb;
  int *occpi, *virtpi, *occ_off, *vir_off;
  double value_c, value_d, dijk, denom, ET;
  double t_ia, t_ib, t_ic, t_ja, t_jb, t_jc, t_ka, t_kb, t_kc;
  double f_ia, f_ib, f_ic, f_ja, f_jb, f_jc, f_ka, f_kb, f_kc;
  double D_jkbc, D_jkac, D_jkba, D_ikbc, D_ikac, D_ikba, D_jibc, D_jiac, D_jiba;
  double t_jkbc, t_jkac, t_jkba, t_ikbc, t_ikac, t_ikba, t_jibc, t_jiac, t_jiba;
  int nrows, ncols, nlinks;
  dpdbuf4 T2, Fints, Eints, Dints;
  dpdfile2 fIJ, fAB, fIA, T1;
  double ***WABC, ***WBCA, ***WACB, ***VABC;
  int nijk, mijk;
  FILE *ijkfile;

  nirreps = moinfo.nirreps;
  occpi = moinfo.boccpi; 
  virtpi = moinfo.bvirtpi;
  occ_off = moinfo.bocc_off;
  vir_off = moinfo.bvir_off;

  dpd_file2_init(&fIJ, CC_OEI, 0, 2, 2, "fij");
  dpd_file2_init(&fAB, CC_OEI, 0, 3, 3, "fab");
  dpd_file2_init(&fIA, CC_OEI, 0, 2, 3, "fia");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_init(&fIA);
  dpd_file2_mat_rd(&fIJ);
  dpd_file2_mat_rd(&fAB);
  dpd_file2_mat_rd(&fIA);

  dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
  dpd_file2_mat_init(&T1); 
  dpd_file2_mat_rd(&T1);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 10, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
  dpd_buf4_init(&Dints, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2, h);
    dpd_buf4_mat_irrep_rd(&T2, h);

    dpd_buf4_mat_irrep_init(&Eints, h);
    dpd_buf4_mat_irrep_rd(&Eints, h);

    dpd_buf4_mat_irrep_init(&Dints, h);
    dpd_buf4_mat_irrep_rd(&Dints, h);
  }

  /* Compute the number of IJK combinations in this spin case */
  nijk = 0;
  for(Gi=0; Gi < nirreps; Gi++)
    for(Gj=0; Gj < nirreps; Gj++)
      for(Gk=0; Gk < nirreps; Gk++)
	for(i=0; i < occpi[Gi]; i++) {
	  I = occ_off[Gi] + i;
	  for(j=0; j < occpi[Gj]; j++) {
	    J = occ_off[Gj] + j;
	    for(k=0; k < occpi[Gk]; k++) {
	      K = occ_off[Gk] + k;

	      if(I > J && J > K) nijk++;
	    }
	  }
	}

  ffile(&ijkfile,"ijk.dat",0);
  fprintf(ijkfile, "Spin Case: BBB\n");
  fprintf(ijkfile, "Number of IJK combintions: %d\n", nijk);
  fprintf(ijkfile, "\nCurrent IJK Combination:\n");
  fflush(ijkfile);

  WABC = (double ***) malloc(nirreps * sizeof(double **));
  VABC = (double ***) malloc(nirreps * sizeof(double **));
  WBCA = (double ***) malloc(nirreps * sizeof(double **));
  WACB = (double ***) malloc(nirreps * sizeof(double **));

  mijk = 0;
  ET = 0.0;

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      for(Gk=0; Gk < nirreps; Gk++) {

	Gij = Gji = Gi ^ Gj;
	Gjk = Gkj = Gj ^ Gk;
	Gik = Gki = Gi ^ Gk;

	Gijk = Gi ^ Gj ^ Gk;

	for(i=0; i < occpi[Gi]; i++) {
	  I = occ_off[Gi] + i;
	  for(j=0; j < occpi[Gj]; j++) {
	    J = occ_off[Gj] + j;
	    for(k=0; k < occpi[Gk]; k++) {
	      K = occ_off[Gk] + k;

	      if(I > J && J > K) {

		mijk++;
		fprintf(ijkfile, "%d\n", mijk);
		fflush(ijkfile);

		ij = Eints.params->rowidx[I][J];
		ji = Eints.params->rowidx[J][I];
		jk = Eints.params->rowidx[J][K];
		kj = Eints.params->rowidx[K][J];
		ik = Eints.params->rowidx[I][K];
		ki = Eints.params->rowidx[K][I];

		dijk = 0.0;
		if(fIJ.params->rowtot[Gi])
		  dijk += fIJ.matrix[Gi][i][i];
		if(fIJ.params->rowtot[Gj])
		  dijk += fIJ.matrix[Gj][j][j];
		if(fIJ.params->rowtot[Gk])
		  dijk += fIJ.matrix[Gk][k][k];

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  WABC[Gab] = dpd_block_matrix(Fints.params->coltot[Gab], virtpi[Gc]);
		}

		for(Gd=0; Gd < nirreps; Gd++) {
		  /* -t_jkcd * F_idab */
		  Gab = Gid = Gi ^ Gd;
		  Gc = Gjk ^ Gd;

		  cd = T2.col_offset[Gjk][Gc];
		  id = Fints.row_offset[Gid][I];

		  Fints.matrix[Gid] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gid]);
		  dpd_buf4_mat_irrep_rd_block(&Fints, Gid, id, virtpi[Gd]);

		  nrows = Fints.params->coltot[Gid];
		  ncols = virtpi[Gc];
		  nlinks = virtpi[Gd];

		  if(nrows && ncols && nlinks) 
		    C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
			    &(Fints.matrix[Gid][0][0]), nrows,
			    &(T2.matrix[Gjk][jk][cd]), nlinks, 1.0,
			    &(WABC[Gab][0][0]), ncols);

		  dpd_free_block(Fints.matrix[Gid], virtpi[Gd], Fints.params->coltot[Gid]);

		  /* +t_ikcd * F_jdab */
		  Gab = Gjd = Gj ^ Gd;
		  Gc = Gik ^ Gd;

		  cd = T2.col_offset[Gik][Gc];
		  jd = Fints.row_offset[Gjd][J];

		  Fints.matrix[Gjd] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gjd]);
		  dpd_buf4_mat_irrep_rd_block(&Fints, Gjd, jd, virtpi[Gd]);

		  nrows = Fints.params->coltot[Gjd];
		  ncols = virtpi[Gc];
		  nlinks = virtpi[Gd];

		  if(nrows && ncols && nlinks) 
		    C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			    &(Fints.matrix[Gjd][0][0]), nrows,
			    &(T2.matrix[Gik][ik][cd]), nlinks, 1.0,
			    &(WABC[Gab][0][0]), ncols);

		  dpd_free_block(Fints.matrix[Gjd], virtpi[Gd], Fints.params->coltot[Gjd]);

		  /* +t_jicd * F_kdab */
		  Gab = Gkd = Gk ^ Gd;
		  Gc = Gji ^ Gd;

		  cd = T2.col_offset[Gji][Gc];
		  kd = Fints.row_offset[Gkd][K];

		  Fints.matrix[Gkd] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gkd]);
		  dpd_buf4_mat_irrep_rd_block(&Fints, Gkd, kd, virtpi[Gd]);

		  nrows = Fints.params->coltot[Gkd];
		  ncols = virtpi[Gc];
		  nlinks = virtpi[Gd];

		  if(nrows && ncols && nlinks) 
		    C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			    &(Fints.matrix[Gkd][0][0]), nrows,
			    &(T2.matrix[Gji][ji][cd]), nlinks, 1.0,
			    &(WABC[Gab][0][0]), ncols);

		  dpd_free_block(Fints.matrix[Gkd], virtpi[Gd], Fints.params->coltot[Gkd]);

		}

		for(Gl=0; Gl < nirreps; Gl++) {
		  /* -t_ilab E_jklc */
		  Gab = Gil = Gi ^ Gl;
		  Gc = Gjk ^ Gl;

		  lc = Eints.col_offset[Gjk][Gl];
		  il = T2.row_offset[Gil][I];

		  nrows = T2.params->coltot[Gil];
		  ncols = virtpi[Gc];
		  nlinks = occpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			    &(T2.matrix[Gil][il][0]), nrows,
			    &(Eints.matrix[Gjk][jk][lc]), ncols, 1.0,
			    &(WABC[Gab][0][0]), ncols);

		  /* +t_jlab E_iklc */
		  Gab = Gjl = Gj ^ Gl;
		  Gc = Gik ^ Gl;

		  lc = Eints.col_offset[Gik][Gl];
		  jl = T2.row_offset[Gjl][J];

		  nrows = T2.params->coltot[Gjl];
		  ncols = virtpi[Gc];
		  nlinks = occpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
			    &(T2.matrix[Gjl][jl][0]), nrows,
			    &(Eints.matrix[Gik][ik][lc]), ncols, 1.0,
			    &(WABC[Gab][0][0]), ncols);

		  /* +t_klab E_jilc */
		  Gab = Gkl = Gk ^ Gl;
		  Gc = Gji ^ Gl;

		  lc = Eints.col_offset[Gji][Gl];
		  kl = T2.row_offset[Gkl][K];

		  nrows = T2.params->coltot[Gkl];
		  ncols = virtpi[Gc];
		  nlinks = occpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
			    &(T2.matrix[Gkl][kl][0]), nrows,
			    &(Eints.matrix[Gji][ji][lc]), ncols, 1.0,
			    &(WABC[Gab][0][0]), ncols);
		}

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  WBCA[Gab] = dpd_block_matrix(Fints.params->coltot[Gab], virtpi[Gc]);
		}

		for(Gd=0; Gd < nirreps; Gd++) {
		  /* -t_jkad * F_idbc */
		  Gbc = Gid = Gi ^ Gd;
		  Ga = Gjk ^ Gd;

		  ad = T2.col_offset[Gjk][Ga];
		  id = Fints.row_offset[Gid][I];

		  Fints.matrix[Gid] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gid]);
		  dpd_buf4_mat_irrep_rd_block(&Fints, Gid, id, virtpi[Gd]);

		  nrows = Fints.params->coltot[Gid];
		  ncols = virtpi[Ga];
		  nlinks = virtpi[Gd];

		  if(nrows && ncols && nlinks) 
		    C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
			    &(Fints.matrix[Gid][0][0]), nrows,
			    &(T2.matrix[Gjk][jk][ad]), nlinks, 1.0,
			    &(WBCA[Gbc][0][0]), ncols);

		  dpd_free_block(Fints.matrix[Gid], virtpi[Gd], Fints.params->coltot[Gid]);

		  /* +t_ikad * F_jdbc */
		  Gbc = Gjd = Gj ^ Gd;
		  Ga = Gik ^ Gd;

		  ad = T2.col_offset[Gik][Ga];
		  jd = Fints.row_offset[Gjd][J];

		  Fints.matrix[Gjd] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gjd]);
		  dpd_buf4_mat_irrep_rd_block(&Fints, Gjd, jd, virtpi[Gd]);

		  nrows = Fints.params->coltot[Gjd];
		  ncols = virtpi[Ga];
		  nlinks = virtpi[Gd];

		  if(nrows && ncols && nlinks) 
		    C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			    &(Fints.matrix[Gjd][0][0]), nrows,
			    &(T2.matrix[Gik][ik][ad]), nlinks, 1.0,
			    &(WBCA[Gbc][0][0]), ncols);

		  dpd_free_block(Fints.matrix[Gjd], virtpi[Gd], Fints.params->coltot[Gjd]);

		  /* +t_jiad * F_kdbc */
		  Gbc = Gkd = Gk ^ Gd;
		  Ga = Gji ^ Gd;

		  ad = T2.col_offset[Gji][Ga];
		  kd = Fints.row_offset[Gkd][K];

		  Fints.matrix[Gkd] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gkd]);
		  dpd_buf4_mat_irrep_rd_block(&Fints, Gkd, kd, virtpi[Gd]);

		  nrows = Fints.params->coltot[Gkd];
		  ncols = virtpi[Ga];
		  nlinks = virtpi[Gd];

		  if(nrows && ncols && nlinks) 
		    C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			    &(Fints.matrix[Gkd][0][0]), nrows,
			    &(T2.matrix[Gji][ji][ad]), nlinks, 1.0,
			    &(WBCA[Gbc][0][0]), ncols);

		  dpd_free_block(Fints.matrix[Gkd], virtpi[Gd], Fints.params->coltot[Gkd]);

		}

		for(Gl=0; Gl < nirreps; Gl++) {
		  /* -t_ilbc * E_jkla */
		  Gbc = Gil = Gi ^ Gl;
		  Ga = Gjk ^ Gl;

		  la = Eints.col_offset[Gjk][Gl];
		  il = T2.row_offset[Gil][I];

		  nrows = T2.params->coltot[Gil];
		  ncols = virtpi[Ga];
		  nlinks = occpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			    &(T2.matrix[Gil][il][0]), nrows,
			    &(Eints.matrix[Gjk][jk][la]), ncols, 1.0,
			    &(WBCA[Gbc][0][0]), ncols);

		  /* +t_jlbc E_ikla */
		  Gbc = Gjl = Gj ^ Gl;
		  Ga = Gik ^ Gl;

		  la = Eints.col_offset[Gik][Gl];
		  jl = T2.row_offset[Gjl][J];

		  nrows = T2.params->coltot[Gjl];
		  ncols = virtpi[Ga];
		  nlinks = occpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
			    &(T2.matrix[Gjl][jl][0]), nrows,
			    &(Eints.matrix[Gik][ik][la]), ncols, 1.0,
			    &(WBCA[Gbc][0][0]), ncols);

		  /* +t_klbc E_jila */
		  Gbc = Gkl = Gk ^ Gl;
		  Ga = Gji ^ Gl;

		  la = Eints.col_offset[Gji][Gl];
		  kl = T2.row_offset[Gkl][K];

		  nrows = T2.params->coltot[Gkl];
		  ncols = virtpi[Ga];
		  nlinks = occpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
			    &(T2.matrix[Gkl][kl][0]), nrows,
			    &(Eints.matrix[Gji][ji][la]), ncols, 1.0,
			    &(WBCA[Gbc][0][0]), ncols);
		}

		dpd_3d_sort(WBCA, WABC, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx,
		       Fints.params->colorb, Fints.params->rsym, Fints.params->ssym,
		       vir_off, vir_off, virtpi, vir_off, Fints.params->colidx, cab, 1);

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  dpd_free_block(WBCA[Gab], Fints.params->coltot[Gab], virtpi[Gc]);

		  WACB[Gab] = dpd_block_matrix(Fints.params->coltot[Gab], virtpi[Gc]);
		}

		for(Gd=0; Gd < nirreps; Gd++) {
		  /* +t_jkbd * F_idac */
		  Gac = Gid = Gi ^ Gd;
		  Gb = Gjk ^ Gd;

		  bd = T2.col_offset[Gjk][Gb];
		  id = Fints.row_offset[Gid][I];

		  Fints.matrix[Gid] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gid]);
		  dpd_buf4_mat_irrep_rd_block(&Fints, Gid, id, virtpi[Gd]);

		  nrows = Fints.params->coltot[Gid];
		  ncols = virtpi[Gb];
		  nlinks = virtpi[Gd];

		  if(nrows && ncols && nlinks) 
		    C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			    &(Fints.matrix[Gid][0][0]), nrows,
			    &(T2.matrix[Gjk][jk][bd]), nlinks, 1.0,
			    &(WACB[Gac][0][0]), ncols);

		  dpd_free_block(Fints.matrix[Gid], virtpi[Gd], Fints.params->coltot[Gid]);

		  /* -t_ikbd * F_jdac */
		  Gac = Gjd = Gj ^ Gd;
		  Gb = Gik ^ Gd;

		  bd = T2.col_offset[Gik][Gb];
		  jd = Fints.row_offset[Gjd][J];

		  Fints.matrix[Gjd] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gjd]);
		  dpd_buf4_mat_irrep_rd_block(&Fints, Gjd, jd, virtpi[Gd]);

		  nrows = Fints.params->coltot[Gjd];
		  ncols = virtpi[Gb];
		  nlinks = virtpi[Gd];

		  if(nrows && ncols && nlinks) 
		    C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
			    &(Fints.matrix[Gjd][0][0]), nrows,
			    &(T2.matrix[Gik][ik][bd]), nlinks, 1.0,
			    &(WACB[Gac][0][0]), ncols);

		  dpd_free_block(Fints.matrix[Gjd], virtpi[Gd], Fints.params->coltot[Gjd]);

		  /* -t_jibd * F_kdac */
		  Gac = Gkd = Gk ^ Gd;
		  Gb = Gji ^ Gd;

		  bd = T2.col_offset[Gji][Gb];
		  kd = Fints.row_offset[Gkd][K];

		  Fints.matrix[Gkd] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gkd]);
		  dpd_buf4_mat_irrep_rd_block(&Fints, Gkd, kd, virtpi[Gd]);

		  nrows = Fints.params->coltot[Gkd];
		  ncols = virtpi[Gb];
		  nlinks = virtpi[Gd];

		  if(nrows && ncols && nlinks) 
		    C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
			    &(Fints.matrix[Gkd][0][0]), nrows,
			    &(T2.matrix[Gji][ji][bd]), nlinks, 1.0,
			    &(WACB[Gac][0][0]), ncols);

		  dpd_free_block(Fints.matrix[Gkd], virtpi[Gd], Fints.params->coltot[Gkd]);

		}

		for(Gl=0; Gl < nirreps; Gl++) {
		  /* +t_ilac * E_jklb */
		  Gac = Gil = Gi ^ Gl;
		  Gb = Gjk ^ Gl;

		  lb = Eints.col_offset[Gjk][Gl];
		  il = T2.row_offset[Gil][I];

		  nrows = T2.params->coltot[Gil];
		  ncols = virtpi[Gb];
		  nlinks = occpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
			    &(T2.matrix[Gil][il][0]), nrows,
			    &(Eints.matrix[Gjk][jk][lb]), ncols, 1.0,
			    &(WACB[Gac][0][0]), ncols);

		  /* -t_jlac * E_iklb */
		  Gac = Gjl = Gj ^ Gl;
		  Gb = Gik ^ Gl;

		  lb = Eints.col_offset[Gik][Gl];
		  jl = T2.row_offset[Gjl][J];

		  nrows = T2.params->coltot[Gjl];
		  ncols = virtpi[Gb];
		  nlinks = occpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			    &(T2.matrix[Gjl][jl][0]), nrows,
			    &(Eints.matrix[Gik][ik][lb]), ncols, 1.0,
			    &(WACB[Gac][0][0]), ncols);

		  /* -t_klac * E_jilb */
		  Gac = Gkl = Gk ^ Gl;
		  Gb = Gji ^ Gl;

		  lb = Eints.col_offset[Gji][Gl];
		  kl = T2.row_offset[Gkl][K];

		  nrows = T2.params->coltot[Gkl];
		  ncols = virtpi[Gb];
		  nlinks = occpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			    &(T2.matrix[Gkl][kl][0]), nrows,
			    &(Eints.matrix[Gji][ji][lb]), ncols, 1.0,
			    &(WACB[Gac][0][0]), ncols);
		}

		dpd_3d_sort(WACB, WABC, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx,
		       Fints.params->colorb, Fints.params->rsym, Fints.params->ssym,
		       vir_off, vir_off, virtpi, vir_off, Fints.params->colidx, acb, 1);

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  dpd_free_block(WACB[Gab], Fints.params->coltot[Gab], virtpi[Gc]);
		}

		/* Add disconnected triples and finish W and V */
		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  VABC[Gab] = dpd_block_matrix(Fints.params->coltot[Gab], virtpi[Gc]);

		  for(ab=0; ab < Fints.params->coltot[Gab]; ab++) {
		    A = Fints.params->colorb[Gab][ab][0];
		    Ga = Fints.params->rsym[A];
		    a = A - vir_off[Ga];
		    B = Fints.params->colorb[Gab][ab][1];
		    Gb = Fints.params->ssym[B];
		    b = B - vir_off[Gb];

		    Gbc = Gb ^ Gc;
		    Gac = Ga ^ Gc;

		    for(c=0; c < virtpi[Gc]; c++) {
		      C = vir_off[Gc] + c;

		      bc = Dints.params->colidx[B][C];
		      ac = Dints.params->colidx[A][C];

		      /* +t_ia * D_jkbc + f_ia * t_jkbc */
		      if(Gi == Ga && Gjk == Gbc) {
			t_ia = D_jkbc = f_ia = t_jkbc = 0.0;

			if(T1.params->rowtot[Gi] && T1.params->coltot[Gi]) {
			  t_ia = T1.matrix[Gi][i][a];
			  f_ia = fIA.matrix[Gi][i][a];
			}

			if(Dints.params->rowtot[Gjk] && Dints.params->coltot[Gjk]) {
			  D_jkbc = Dints.matrix[Gjk][jk][bc];
			  t_jkbc = T2.matrix[Gjk][jk][bc];
			}

			VABC[Gab][ab][c] += t_ia * D_jkbc + f_ia * t_jkbc;
		      }

		      /* -t_ib * D_jkac - f_ib * t_jkac */
		      if(Gi == Gb && Gjk == Gac) {
			t_ib = D_jkac = f_ib = t_jkac = 0.0;

			if(T1.params->rowtot[Gi] && T1.params->coltot[Gi]) {
			  t_ib = T1.matrix[Gi][i][b];
			  f_ib = fIA.matrix[Gi][i][b];
			}

			if(Dints.params->rowtot[Gjk] && Dints.params->coltot[Gjk]) {
			  D_jkac = Dints.matrix[Gjk][jk][ac];
			  t_jkac = T2.matrix[Gjk][jk][ac];
			}

			VABC[Gab][ab][c] -= t_ib * D_jkac + f_ib * t_jkac;
		      }

		      /* +t_ic * D_jkab + f_ic * t_jkba */
		      if(Gi == Gc && Gjk == Gab) {
			t_ic = D_jkba = f_ic = t_jkba = 0.0;

			if(T1.params->rowtot[Gi] && T1.params->coltot[Gi]) {
			  t_ic = T1.matrix[Gi][i][c];
			  f_ic = fIA.matrix[Gi][i][c];
			}

			if(Dints.params->rowtot[Gjk] && Dints.params->coltot[Gjk]) {
			  D_jkba = Dints.matrix[Gjk][jk][ab];
			  t_jkba = T2.matrix[Gjk][jk][ab];
			}

			VABC[Gab][ab][c] += t_ic * D_jkba + f_ic * t_jkba;
		      }

		      /* -t_ja * D_ikbc - f_ja * t_ikbc*/
		      if(Gj == Ga && Gik == Gbc) {
			t_ja = D_ikbc = f_ja = t_ikbc = 0.0;

			if(T1.params->rowtot[Gj] && T1.params->coltot[Gj]) {
			  t_ja = T1.matrix[Gj][j][a];
			  f_ja = fIA.matrix[Gj][j][a];
			}

			if(Dints.params->rowtot[Gik] && Dints.params->coltot[Gik]) {
			  D_ikbc = Dints.matrix[Gik][ik][bc];
			  t_ikbc = T2.matrix[Gik][ik][bc];
			}

			VABC[Gab][ab][c] -= t_ja * D_ikbc + f_ja * t_ikbc;
		      }

		      /* +t_jb * D_ikac + f_jb * t_ikac */
		      if(Gj == Gb && Gik == Gac) {
			t_jb = D_ikac = f_jb = t_ikac = 0.0;

			if(T1.params->rowtot[Gj] && T1.params->coltot[Gj]) {
			  t_jb = T1.matrix[Gj][j][b];
			  f_jb = fIA.matrix[Gj][j][b];
			}

			if(Dints.params->rowtot[Gik] && Dints.params->coltot[Gik]) {
			  D_ikac = Dints.matrix[Gik][ik][ac];
			  t_ikac = T2.matrix[Gik][ik][ac];
			}

			VABC[Gab][ab][c] += t_jb * D_ikac + f_jb * t_ikac;
		      }

		      /* -t_jc * D_ikba - f_jc * t_ikba */
		      if(Gj == Gc && Gik == Gab) {
			t_jc = D_ikba = f_jc = t_ikba = 0.0;

			if(T1.params->rowtot[Gj] && T1.params->coltot[Gj]) {
			  t_jc = T1.matrix[Gj][j][c];
			  f_jc = fIA.matrix[Gj][j][c];
			}

			if(Dints.params->rowtot[Gik] && Dints.params->coltot[Gik]) {
			  D_ikba = Dints.matrix[Gik][ik][ab];
			  t_ikba = T2.matrix[Gik][ik][ab];
			}

			VABC[Gab][ab][c] -= t_jc * D_ikba + f_jc * t_ikba;
		      }

		      /* -t_ka * D_jibc - f_ka * t_jibc */
		      if(Gk == Ga && Gji == Gbc) {
			t_ka = D_jibc = f_ka = t_jibc = 0.0;

			if(T1.params->rowtot[Gk] && T1.params->coltot[Gk]) {
			  t_ka = T1.matrix[Gk][k][a];
			  f_ka = fIA.matrix[Gk][k][a];
			}

			if(Dints.params->rowtot[Gji] && Dints.params->coltot[Gji]) {
			  D_jibc = Dints.matrix[Gji][ji][bc];
			  t_jibc = T2.matrix[Gji][ji][bc];
			}

			VABC[Gab][ab][c] -= t_ka * D_jibc + f_ka * t_jibc;
		      }

		      /* +t_kb * D_jiac + f_kb * t_jiac */
		      if(Gk == Gb && Gji == Gac) {
			t_kb = D_jiac = f_kb = t_jiac = 0.0;

			if(T1.params->rowtot[Gk] && T1.params->coltot[Gk]) {
			  t_kb = T1.matrix[Gk][k][b];
			  f_kb = fIA.matrix[Gk][k][b];
			}

			if(Dints.params->rowtot[Gji] && Dints.params->coltot[Gji]) {
			  D_jiac = Dints.matrix[Gji][ji][ac];
			  t_jiac = T2.matrix[Gji][ji][ac];
			}

			VABC[Gab][ab][c] += t_kb * D_jiac + f_kb * t_jiac;
		      }

		      /* -t_kc * D_jiab - f_kc * t_jiba*/
		      if(Gk == Gc && Gji == Gab) {
			t_kc = D_jiba = f_kc = t_jiba = 0.0;

			if(T1.params->rowtot[Gk] && T1.params->coltot[Gk]) {
			  t_kc = T1.matrix[Gk][k][c];
			  f_kc = fIA.matrix[Gk][k][c];
			}

			if(Dints.params->rowtot[Gji] && Dints.params->coltot[Gji]) {
			  D_jiba = Dints.matrix[Gji][ji][ab];
			  t_jiba = T2.matrix[Gji][ji][ab];
			}

			VABC[Gab][ab][c] -= t_kc * D_jiba + f_kc * t_jiba;
		      }

		      /*
		      if(fabs(VABC[Gab][ab][c]) > 1e-7)
			fprintf(outfile, "%d %d %d %d %d %d %20.15f\n", I,J,K,A,B,C,VABC[Gab][ab][c]);
		      */

		      /* Sum V and W into V */
		      VABC[Gab][ab][c] += WABC[Gab][ab][c];

		      /* Build the rest of the denominator and divide it into W */
		      denom = dijk;
		      if(fAB.params->rowtot[Ga])
			denom -= fAB.matrix[Ga][a][a];
		      if(fAB.params->rowtot[Gb])
			denom -= fAB.matrix[Gb][b][b];
		      if(fAB.params->rowtot[Gc])
			denom -= fAB.matrix[Gc][c][c];

		      WABC[Gab][ab][c] /= denom;

		    } /* c */
		  } /* ab */
		} /* Gab */

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;
		  ET += dot_block(WABC[Gab], VABC[Gab], Fints.params->coltot[Gab], virtpi[Gc], 1.0/6.0);
		  dpd_free_block(WABC[Gab], Fints.params->coltot[Gab], virtpi[Gc]);
		  dpd_free_block(VABC[Gab], Fints.params->coltot[Gab], virtpi[Gc]);
		}

	      } /* I >= J >= K */

	    } /* k */
	  } /* j */
	} /* i */

      } /* Gk */
    } /* Gj */
  } /* Gi */

  free(WABC);
  free(VABC);
  free(WBCA);
  free(WACB);

  fclose(ijkfile);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&T2, h);
    dpd_buf4_mat_irrep_close(&Eints, h);
    dpd_buf4_mat_irrep_close(&Dints, h);
  }

  dpd_buf4_close(&T2);
  dpd_buf4_close(&Fints);
  dpd_buf4_close(&Eints);
  dpd_buf4_close(&Dints);

  dpd_file2_mat_close(&T1);
  dpd_file2_close(&T1);

  dpd_file2_mat_close(&fIJ);
  dpd_file2_mat_close(&fAB);
  dpd_file2_mat_close(&fIA);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);
  dpd_file2_close(&fIA);

  return ET;
}

}} // namespace psi::CCTRIPLES
