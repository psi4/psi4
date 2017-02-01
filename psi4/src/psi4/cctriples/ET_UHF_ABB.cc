/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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
    \ingroup CCTRIPLES
    \brief Enter brief description of file here
*/
 #include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"
#include "psi4/libparallel/ParallelPrinter.h"
namespace psi { namespace cctriples {

double ET_UHF_ABB(void)
{
  int cnt;
  int h, nirreps;
  int Gi, Gj, Gk, Ga, Gb, Gc, Gd, Gl;
  int Gji, Gij, Gjk, Gkj, Gik, Gki, Gijk;
  int Gab, Gbc, Gac, Gca, Gba;
  int Gid, Gjd, Gkd;
  int Gil, Gjl, Gkl;
  int I, J, K, A, B, C;
  int i, j, k, a, b, c;
  int ij, ji, ik, ki, jk, kj;
  int ab, ba, ac, ca, bc, cb;
  int cd, bd, ad, db, dc;
  int lc, lb, la;
  int id, jd, kd;
  int il, jl, kl;
  int *aoccpi, *avirtpi, *aocc_off, *avir_off;
  int *boccpi, *bvirtpi, *bocc_off, *bvir_off;
  double value_c, value_d, dijk, denom, ET_ABB;
  int nrows, ncols, nlinks;
  double t_ia, t_jb, t_jc, t_kb, t_kc;
  double f_ia, f_jb, f_jc, f_kb, f_kc;
  double D_jkbc, D_ikac, D_ikab, D_ijac, D_ijab;
  double t_jkbc, t_ikac, t_ikab, t_ijac, t_ijab;
  dpdbuf4 T2AB, T2BB, T2BA;
  dpdbuf4 FBBints, FABints, FBAints;
  dpdbuf4 EBBints, EABints, EBAints;
  dpdbuf4 DBBints, DABints;
  dpdfile2 T1A, T1B, fIJ, fij, fAB, fab, fIA, fia;
  double ***WAbc, ***VAbc, ***WAcb, ***WbAc, ***WcAb, ***WbcA;
  int nijk, mijk;

  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi;
  avirtpi = moinfo.avirtpi;
  aocc_off = moinfo.aocc_off;
  avir_off = moinfo.avir_off;
  boccpi = moinfo.boccpi;
  bvirtpi = moinfo.bvirtpi;
  bocc_off = moinfo.bocc_off;
  bvir_off = moinfo.bvir_off;

  global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
  global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");
  global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
  global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
  global_dpd_->file2_mat_init(&fIJ);
  global_dpd_->file2_mat_init(&fij);
  global_dpd_->file2_mat_init(&fAB);
  global_dpd_->file2_mat_init(&fab);
  global_dpd_->file2_mat_init(&fIA);
  global_dpd_->file2_mat_init(&fia);
  global_dpd_->file2_mat_rd(&fIJ);
  global_dpd_->file2_mat_rd(&fij);
  global_dpd_->file2_mat_rd(&fAB);
  global_dpd_->file2_mat_rd(&fab);
  global_dpd_->file2_mat_rd(&fIA);
  global_dpd_->file2_mat_rd(&fia);

  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1A);
  global_dpd_->file2_mat_rd(&T1A);
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->file2_mat_init(&T1B);
  global_dpd_->file2_mat_rd(&T1B);

  global_dpd_->buf4_init(&T2BB, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
  global_dpd_->buf4_init(&T2AB, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  global_dpd_->buf4_init(&T2BA, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");

  global_dpd_->buf4_init(&FBBints, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
  global_dpd_->buf4_init(&FABints, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
  global_dpd_->buf4_init(&FBAints, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");

  global_dpd_->buf4_init(&EBBints, PSIF_CC_EINTS, 0, 10, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
  global_dpd_->buf4_init(&EABints, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
  global_dpd_->buf4_init(&EBAints, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");

  global_dpd_->buf4_init(&DBBints, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
  global_dpd_->buf4_init(&DABints, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&T2BB, h);
    global_dpd_->buf4_mat_irrep_rd(&T2BB, h);

    global_dpd_->buf4_mat_irrep_init(&T2AB, h);
    global_dpd_->buf4_mat_irrep_rd(&T2AB, h);

    global_dpd_->buf4_mat_irrep_init(&T2BA, h);
    global_dpd_->buf4_mat_irrep_rd(&T2BA, h);

    global_dpd_->buf4_mat_irrep_init(&EBBints, h);
    global_dpd_->buf4_mat_irrep_rd(&EBBints, h);

    global_dpd_->buf4_mat_irrep_init(&EABints, h);
    global_dpd_->buf4_mat_irrep_rd(&EABints, h);

    global_dpd_->buf4_mat_irrep_init(&EBAints, h);
    global_dpd_->buf4_mat_irrep_rd(&EBAints, h);

    global_dpd_->buf4_mat_irrep_init(&DBBints, h);
    global_dpd_->buf4_mat_irrep_rd(&DBBints, h);

    global_dpd_->buf4_mat_irrep_init(&DABints, h);
    global_dpd_->buf4_mat_irrep_rd(&DABints, h);
  }

  /* Compute the number of IJK combinations in this spin case */
  nijk = 0;
  for(Gi=0; Gi < nirreps; Gi++)
    for(Gj=0; Gj < nirreps; Gj++)
      for(Gk=0; Gk < nirreps; Gk++)
	for(i=0; i < aoccpi[Gi]; i++) {
	  I = aocc_off[Gi] + i;
	  for(j=0; j < boccpi[Gj]; j++) {
	    J = bocc_off[Gj] + j;
	    for(k=0; k < boccpi[Gk]; k++) {
	      K = bocc_off[Gk] + k;

	      if(J > K) nijk++;
	    }
	  }
	}
  std::shared_ptr<OutFile> printer(new OutFile("ijk.dat",TRUNCATE));
  //ffile(&ijkfile,"ijk.dat",0);
  printer->Printf("Spin Case: ABB\n");
  printer->Printf("Number of IJK combintions: %d\n", nijk);
  printer->Printf("\nCurrent IJK Combination:\n");


  mijk = 0;
  ET_ABB = 0.0;

  WAbc = (double ***) malloc(nirreps * sizeof(double **));
  VAbc = (double ***) malloc(nirreps * sizeof(double **));
  WAcb = (double ***) malloc(nirreps * sizeof(double **));
  WbcA = (double ***) malloc(nirreps * sizeof(double **));
  WcAb = (double ***) malloc(nirreps * sizeof(double **));
  WbAc = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      for(Gk=0; Gk < nirreps; Gk++) {

	Gij = Gji = Gi ^ Gj;
	Gjk = Gkj = Gj ^ Gk;
	Gik = Gki = Gi ^ Gk;

	Gijk = Gi ^ Gj ^ Gk;

	for(i=0; i < aoccpi[Gi]; i++) {
	  I = aocc_off[Gi] + i;
	  for(j=0; j < boccpi[Gj]; j++) {
	    J = bocc_off[Gj] + j;
	    for(k=0; k < boccpi[Gk]; k++) {
	      K = bocc_off[Gk] + k;

	      if(J > K) {

		mijk++;
		printer->Printf("%d\n", mijk);


		ij = EABints.params->rowidx[I][J];
		ji = EBAints.params->rowidx[J][I];
		jk = EBBints.params->rowidx[J][K];
		kj = EBBints.params->rowidx[K][J];
		ik = EABints.params->rowidx[I][K];
		ki = EBAints.params->rowidx[K][I];

		dijk = 0.0;
		if(fIJ.params->rowtot[Gi])
		  dijk += fIJ.matrix[Gi][i][i];
		if(fij.params->rowtot[Gj])
		  dijk += fij.matrix[Gj][j][j];
		if(fij.params->rowtot[Gk])
		  dijk += fij.matrix[Gk][k][k];

		/* Begin connected triples */

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  WAbc[Gab] = global_dpd_->dpd_block_matrix(FABints.params->coltot[Gab], bvirtpi[Gc]);
		}

		for(Gd=0; Gd < nirreps; Gd++) {
		  /* -t_jkcd * F_IdAb */
		  Gab = Gid = Gi ^ Gd;
		  Gc = Gjk ^ Gd;

		  cd = T2BB.col_offset[Gjk][Gc];
		  id = FABints.row_offset[Gid][I];

		  FABints.matrix[Gid] = global_dpd_->dpd_block_matrix(bvirtpi[Gd], FABints.params->coltot[Gid]);
		  global_dpd_->buf4_mat_irrep_rd_block(&FABints, Gid, id, bvirtpi[Gd]);

		  nrows = FABints.params->coltot[Gid];
		  ncols = bvirtpi[Gc];
		  nlinks = bvirtpi[Gd];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
			    &(FABints.matrix[Gid][0][0]), nrows,
			    &(T2BB.matrix[Gjk][jk][cd]), nlinks, 1.0,
			    &(WAbc[Gab][0][0]), ncols);

		  global_dpd_->free_dpd_block(FABints.matrix[Gid], bvirtpi[Gd], FABints.params->coltot[Gid]);

		}

		for(Gl=0; Gl < nirreps; Gl++) {
		  /* -t_IlAb E_jklc */
		  Gab = Gil = Gi ^ Gl;
		  Gc = Gjk ^ Gl;

		  lc = EBBints.col_offset[Gjk][Gl];
		  il = T2AB.row_offset[Gil][I];

		  nrows = T2AB.params->coltot[Gil];
		  ncols = bvirtpi[Gc];
		  nlinks = boccpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			    &(T2AB.matrix[Gil][il][0]), nrows,
			    &(EBBints.matrix[Gjk][jk][lc]), ncols, 1.0,
			    &(WAbc[Gab][0][0]), ncols);
		}

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  WAcb[Gab] = global_dpd_->dpd_block_matrix(FABints.params->coltot[Gab], bvirtpi[Gc]);
		}

		for(Gd=0; Gd < nirreps; Gd++) {
		  /* +t_jkbd * F_IdAc */
		  Gac = Gid = Gi ^ Gd;
		  Gb = Gjk ^ Gd;

		  bd = T2BB.col_offset[Gjk][Gb];
		  id = FABints.row_offset[Gid][I];

		  FABints.matrix[Gid] = global_dpd_->dpd_block_matrix(bvirtpi[Gd], FABints.params->coltot[Gid]);
		  global_dpd_->buf4_mat_irrep_rd_block(&FABints, Gid, id, bvirtpi[Gd]);

		  nrows = FABints.params->coltot[Gid];
		  ncols = bvirtpi[Gb];
		  nlinks = bvirtpi[Gd];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			    &(FABints.matrix[Gid][0][0]), nrows,
			    &(T2BB.matrix[Gjk][jk][bd]), nlinks, 1.0,
			    &(WAcb[Gac][0][0]), ncols);

		  global_dpd_->free_dpd_block(FABints.matrix[Gid], bvirtpi[Gd], FABints.params->coltot[Gid]);

		}

		for(Gl=0; Gl < nirreps; Gl++) {
		  /* +t_IlAc E_jklb */
		  Gac = Gil = Gi ^ Gl;
		  Gb = Gjk ^ Gl;

		  lb = EBBints.col_offset[Gjk][Gl];
		  il = T2AB.row_offset[Gil][I];

		  nrows = T2AB.params->coltot[Gil];
		  ncols = bvirtpi[Gb];
		  nlinks = boccpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
			    &(T2AB.matrix[Gil][il][0]), nrows,
			    &(EBBints.matrix[Gjk][jk][lb]), ncols, 1.0,
			    &(WAcb[Gac][0][0]), ncols);
		}

        global_dpd_->sort_3d(WAcb, WAbc, nirreps, Gijk, FABints.params->coltot, FABints.params->colidx,
		       FABints.params->colorb, FABints.params->rsym, FABints.params->ssym,
		       avir_off, bvir_off, bvirtpi, bvir_off, FABints.params->colidx, acb, 1);

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  global_dpd_->free_dpd_block(WAcb[Gab], FABints.params->coltot[Gab], bvirtpi[Gc]);

		  WbcA[Gab] = global_dpd_->dpd_block_matrix(FBBints.params->coltot[Gab], avirtpi[Gc]);
		}

		for(Gd=0; Gd < nirreps; Gd++) {
		  /* +t_IkAd * F_jdbc */
		  Gbc = Gjd = Gj ^ Gd;
		  Ga = Gik ^ Gd;

		  ad = T2AB.col_offset[Gik][Ga];
		  jd = FBBints.row_offset[Gjd][J];

		  FBBints.matrix[Gjd] = global_dpd_->dpd_block_matrix(bvirtpi[Gd], FBBints.params->coltot[Gjd]);
		  global_dpd_->buf4_mat_irrep_rd_block(&FBBints, Gjd, jd, bvirtpi[Gd]);

		  nrows = FBBints.params->coltot[Gjd];
		  ncols = avirtpi[Ga];
		  nlinks = bvirtpi[Gd];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			    &(FBBints.matrix[Gjd][0][0]), nrows,
			    &(T2AB.matrix[Gik][ik][ad]), nlinks, 1.0,
			    &(WbcA[Gbc][0][0]), ncols);

		  global_dpd_->free_dpd_block(FBBints.matrix[Gjd], bvirtpi[Gd], FBBints.params->coltot[Gjd]);

		  /* -t_IjAd * F_kdbc */
		  Gbc = Gkd = Gk ^ Gd;
		  Ga = Gij ^ Gd;

		  ad = T2AB.col_offset[Gij][Ga];
		  kd = FBBints.row_offset[Gkd][K];

		  FBBints.matrix[Gkd] = global_dpd_->dpd_block_matrix(bvirtpi[Gd], FBBints.params->coltot[Gkd]);
		  global_dpd_->buf4_mat_irrep_rd_block(&FBBints, Gkd, kd, bvirtpi[Gd]);

		  nrows = FBBints.params->coltot[Gkd];
		  ncols = avirtpi[Ga];
		  nlinks = bvirtpi[Gd];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
			    &(FBBints.matrix[Gkd][0][0]), nrows,
			    &(T2AB.matrix[Gij][ij][ad]), nlinks, 1.0,
			    &(WbcA[Gbc][0][0]), ncols);

		  global_dpd_->free_dpd_block(FBBints.matrix[Gkd], bvirtpi[Gd], FBBints.params->coltot[Gkd]);

		}

		for(Gl=0; Gl < nirreps; Gl++) {
		  /* -t_jlbc E_kIlA */
		  Gbc = Gjl = Gj ^ Gl;
		  Ga = Gki ^ Gl;

		  la = EBAints.col_offset[Gki][Gl];
		  jl = T2BB.row_offset[Gjl][J];

		  nrows = T2BB.params->coltot[Gjl];
		  ncols = avirtpi[Ga];
		  nlinks = boccpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			    &(T2BB.matrix[Gjl][jl][0]), nrows,
			    &(EBAints.matrix[Gki][ki][la]), ncols, 1.0,
			    &(WbcA[Gbc][0][0]), ncols);

		  /* +t_klbc E_jIlA */
		  Gbc = Gkl = Gk ^ Gl;
		  Ga = Gji ^ Gl;

		  la = EBAints.col_offset[Gji][Gl];
		  kl = T2BB.row_offset[Gkl][K];

		  nrows = T2BB.params->coltot[Gkl];
		  ncols = avirtpi[Ga];
		  nlinks = boccpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
			    &(T2BB.matrix[Gkl][kl][0]), nrows,
			    &(EBAints.matrix[Gji][ji][la]), ncols, 1.0,
			    &(WbcA[Gbc][0][0]), ncols);
		}

        global_dpd_->sort_3d(WbcA, WAbc, nirreps, Gijk, FBBints.params->coltot, FBBints.params->colidx,
		       FBBints.params->colorb, FBBints.params->rsym, FBBints.params->ssym,
		       bvir_off, bvir_off, avirtpi, avir_off, FABints.params->colidx, cab, 1);

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  global_dpd_->free_dpd_block(WbcA[Gab], FBBints.params->coltot[Gab], avirtpi[Gc]);

		  WcAb[Gab] = global_dpd_->dpd_block_matrix(FBAints.params->coltot[Gab], bvirtpi[Gc]);
		}

		for(Gd=0; Gd < nirreps; Gd++) {
		  /* -t_IkDb * F_jDcA */
		  Gca = Gjd = Gj ^ Gd;
		  Gb = Gik ^ Gd;

		  db = T2AB.col_offset[Gik][Gd];
		  jd = FBAints.row_offset[Gjd][J];

		  FBAints.matrix[Gjd] = global_dpd_->dpd_block_matrix(avirtpi[Gd], FBAints.params->coltot[Gjd]);
		  global_dpd_->buf4_mat_irrep_rd_block(&FBAints, Gjd, jd, avirtpi[Gd]);

		  nrows = FBAints.params->coltot[Gjd];
		  ncols = bvirtpi[Gb];
		  nlinks = avirtpi[Gd];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			    &(FBAints.matrix[Gjd][0][0]), nrows,
			    &(T2AB.matrix[Gik][ik][db]), ncols, 1.0,
			    &(WcAb[Gca][0][0]), ncols);

		  global_dpd_->free_dpd_block(FBAints.matrix[Gjd], avirtpi[Gd], FBAints.params->coltot[Gjd]);

		  /* +t_IjDb * F_kDcA */
		  Gca = Gkd = Gk ^ Gd;
		  Gb = Gij ^ Gd;

		  db = T2AB.col_offset[Gij][Gd];
		  kd = FBAints.row_offset[Gkd][K];

		  FBAints.matrix[Gkd] = global_dpd_->dpd_block_matrix(avirtpi[Gd], FBAints.params->coltot[Gkd]);
		  global_dpd_->buf4_mat_irrep_rd_block(&FBAints, Gkd, kd, avirtpi[Gd]);

		  nrows = FBAints.params->coltot[Gkd];
		  ncols = bvirtpi[Gb];
		  nlinks = avirtpi[Gd];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
			    &(FBAints.matrix[Gkd][0][0]), nrows,
			    &(T2AB.matrix[Gij][ij][db]), ncols, 1.0,
			    &(WcAb[Gca][0][0]), ncols);

		  global_dpd_->free_dpd_block(FBAints.matrix[Gkd], avirtpi[Gd], FBAints.params->coltot[Gkd]);

		}

		for(Gl=0; Gl < nirreps; Gl++) {
		  /* +t_jLcA E_IkLb */
		  Gca = Gjl = Gj ^ Gl;
		  Gb = Gik ^ Gl;

		  lb = EABints.col_offset[Gik][Gl];
		  jl = T2BA.row_offset[Gjl][J];

		  nrows = T2BA.params->coltot[Gjl];
		  ncols = bvirtpi[Gb];
		  nlinks = aoccpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
			    &(T2BA.matrix[Gjl][jl][0]), nrows,
			    &(EABints.matrix[Gik][ik][lb]), ncols, 1.0,
			    &(WcAb[Gca][0][0]), ncols);

		  /* -t_kLcA E_IjLb */
		  Gca = Gkl = Gk ^ Gl;
		  Gb = Gij ^ Gl;

		  lb = EABints.col_offset[Gij][Gl];
		  kl = T2BA.row_offset[Gkl][K];

		  nrows = T2BA.params->coltot[Gkl];
		  ncols = bvirtpi[Gb];
		  nlinks = aoccpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			    &(T2BA.matrix[Gkl][kl][0]), nrows,
			    &(EABints.matrix[Gij][ij][lb]), ncols, 1.0,
			    &(WcAb[Gca][0][0]), ncols);
		}

        global_dpd_->sort_3d(WcAb, WAbc, nirreps, Gijk, FBAints.params->coltot, FBAints.params->colidx,
		       FBAints.params->colorb, FBAints.params->rsym, FBAints.params->ssym,
		       bvir_off, avir_off, bvirtpi, bvir_off, FABints.params->colidx, bca, 1);

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  global_dpd_->free_dpd_block(WcAb[Gab], FBAints.params->coltot[Gab], bvirtpi[Gc]);

		  WbAc[Gab] = global_dpd_->dpd_block_matrix(FBAints.params->coltot[Gab], bvirtpi[Gc]);
		}

		for(Gd=0; Gd < nirreps; Gd++) {
		  /* +t_IkDc * F_jDbA */
		  Gba = Gjd = Gj ^ Gd;
		  Gc = Gik ^ Gd;

		  dc = T2AB.col_offset[Gik][Gd];
		  jd = FBAints.row_offset[Gjd][J];

		  FBAints.matrix[Gjd] = global_dpd_->dpd_block_matrix(avirtpi[Gd], FBAints.params->coltot[Gjd]);
		  global_dpd_->buf4_mat_irrep_rd_block(&FBAints, Gjd, jd, avirtpi[Gd]);

		  nrows = FBAints.params->coltot[Gjd];
		  ncols = bvirtpi[Gc];
		  nlinks = avirtpi[Gd];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
			    &(FBAints.matrix[Gjd][0][0]), nrows,
			    &(T2AB.matrix[Gik][ik][dc]), ncols, 1.0,
			    &(WbAc[Gba][0][0]), ncols);

		  global_dpd_->free_dpd_block(FBAints.matrix[Gjd], avirtpi[Gd], FBAints.params->coltot[Gjd]);

		  /* -t_IjDc * F_kDbA */
		  Gba = Gkd = Gk ^ Gd;
		  Gc = Gij ^ Gd;

		  dc = T2AB.col_offset[Gij][Gd];
		  kd = FBAints.row_offset[Gkd][K];

		  FBAints.matrix[Gkd] = global_dpd_->dpd_block_matrix(avirtpi[Gd], FBAints.params->coltot[Gkd]);
		  global_dpd_->buf4_mat_irrep_rd_block(&FBAints, Gkd, kd, avirtpi[Gd]);

		  nrows = FBAints.params->coltot[Gkd];
		  ncols = bvirtpi[Gc];
		  nlinks = avirtpi[Gd];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			    &(FBAints.matrix[Gkd][0][0]), nrows,
			    &(T2AB.matrix[Gij][ij][dc]), ncols, 1.0,
			    &(WbAc[Gba][0][0]), ncols);

		  global_dpd_->free_dpd_block(FBAints.matrix[Gkd], avirtpi[Gd], FBAints.params->coltot[Gkd]);

		}

		for(Gl=0; Gl < nirreps; Gl++) {
		  /* -t_jLbA * E_IkLc */
		  Gba = Gjl = Gj ^ Gl;
		  Gc = Gik ^ Gl;

		  lc = EABints.col_offset[Gik][Gl];
		  jl = T2BA.row_offset[Gjl][J];

		  nrows = T2BA.params->coltot[Gjl];
		  ncols = bvirtpi[Gc];
		  nlinks = aoccpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			    &(T2BA.matrix[Gjl][jl][0]), nrows,
			    &(EABints.matrix[Gik][ik][lc]), ncols, 1.0,
			    &(WbAc[Gba][0][0]), ncols);

		  /* +t_kLbA * E_IjLc */
		  Gba = Gkl = Gk ^ Gl;
		  Gc = Gij ^ Gl;

		  lc = EABints.col_offset[Gij][Gl];
		  kl = T2BA.row_offset[Gkl][K];

		  nrows = T2BA.params->coltot[Gkl];
		  ncols = bvirtpi[Gc];
		  nlinks = aoccpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
			    &(T2BA.matrix[Gkl][kl][0]), nrows,
			    &(EABints.matrix[Gij][ij][lc]), ncols, 1.0,
			    &(WbAc[Gba][0][0]), ncols);
		}

        global_dpd_->sort_3d(WbAc, WAbc, nirreps, Gijk, FBAints.params->coltot, FBAints.params->colidx,
		       FBAints.params->colorb, FBAints.params->rsym, FBAints.params->ssym,
		       bvir_off, avir_off, bvirtpi, bvir_off, FABints.params->colidx, bac, 1);

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  global_dpd_->free_dpd_block(WbAc[Gab], FBAints.params->coltot[Gab], bvirtpi[Gc]);
		}

		/* Add disconnected triples and finish W and V */
		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  VAbc[Gab] = global_dpd_->dpd_block_matrix(FABints.params->coltot[Gab], bvirtpi[Gc]);

		  for(ab=0; ab < FABints.params->coltot[Gab]; ab++) {
		    A = FABints.params->colorb[Gab][ab][0];
		    Ga = FABints.params->rsym[A];
		    a = A - avir_off[Ga];
		    B = FABints.params->colorb[Gab][ab][1];
		    Gb = FABints.params->ssym[B];
		    b = B - bvir_off[Gb];

		    Gbc = Gb ^ Gc;
		    Gac = Ga ^ Gc;

		    for(c=0; c < bvirtpi[Gc]; c++) {
		      C = bvir_off[Gc] + c;

		      bc = DBBints.params->colidx[B][C];
		      ac = DABints.params->colidx[A][C];

		      /* +t_IA * D_jkbc + f_IA * t_jkbc */
		      if(Gi == Ga && Gjk == Gbc) {
			t_ia = D_jkbc = f_ia = t_jkbc = 0.0;

			if(T1A.params->rowtot[Gi] && T1A.params->coltot[Gi]) {
			  t_ia = T1A.matrix[Gi][i][a];
			  f_ia = fIA.matrix[Gi][i][a];
			}

			if(DBBints.params->rowtot[Gjk] && DBBints.params->coltot[Gjk]) {
			  D_jkbc = DBBints.matrix[Gjk][jk][bc];
			  t_jkbc = T2BB.matrix[Gjk][jk][bc];
			}

			VAbc[Gab][ab][c] += t_ia * D_jkbc + f_ia * t_jkbc;
		      }

		      /* +t_jb * D_IkAc + f_jb * t_IkAc */
		      if(Gj == Gb && Gik == Gac) {
			t_jb = D_ikac = f_jb = t_ikac = 0.0;

			if(T1B.params->rowtot[Gj] && T1B.params->coltot[Gj]) {
			  t_jb = T1B.matrix[Gj][j][b];
			  f_jb = fia.matrix[Gj][j][b];
			}

			if(DABints.params->rowtot[Gik] && DABints.params->coltot[Gik]) {
			  D_ikac = DABints.matrix[Gik][ik][ac];
			  t_ikac = T2AB.matrix[Gik][ik][ac];
			}

			VAbc[Gab][ab][c] += t_jb * D_ikac + f_jb * t_ikac;
		      }

		      /* -t_jc * D_IkAb - f_jc * t_IkAb */
		      if(Gj == Gc && Gik == Gab) {
			t_jc = D_ikab = f_jc = t_ikab = 0.0;

			if(T1B.params->rowtot[Gj] && T1B.params->coltot[Gj]) {
			  t_jc = T1B.matrix[Gj][j][c];
			  f_jc = fia.matrix[Gj][j][c];
			}

			if(DABints.params->rowtot[Gik] && DABints.params->coltot[Gik]) {
			  D_ikab = DABints.matrix[Gik][ik][ab];
			  t_ikab = T2AB.matrix[Gik][ik][ab];
			}

			VAbc[Gab][ab][c] -= t_jc * D_ikab + f_jc * t_ikab;
		      }

		      /* -t_kb * D_IjAc - f_kb * t_IjAc */
		      if(Gk == Gb && Gji == Gac) {
			t_kb = D_ijac = f_kb = t_ijac = 0.0;

			if(T1B.params->rowtot[Gk] && T1B.params->coltot[Gk]) {
			  t_kb = T1B.matrix[Gk][k][b];
			  f_kb = fia.matrix[Gk][k][b];
			}

			if(DABints.params->rowtot[Gji] && DABints.params->coltot[Gji]) {
			  D_ijac = DABints.matrix[Gji][ij][ac];
			  t_ijac = T2AB.matrix[Gji][ij][ac];
			}

			VAbc[Gab][ab][c] -= t_kb * D_ijac + f_kb * t_ijac;
		      }

		      /* +t_kc * D_IjAb + f_kc * t_IjAb */
		      if(Gk == Gc && Gji == Gab) {
			t_kc = D_ijab = f_kc = t_ijab = 0.0;

			if(T1B.params->rowtot[Gk] && T1B.params->coltot[Gk]) {
			  t_kc = T1B.matrix[Gk][k][c];
			  f_kc = fia.matrix[Gk][k][c];
			}

			if(DABints.params->rowtot[Gji] && DABints.params->coltot[Gji]) {
			  D_ijab = DABints.matrix[Gji][ij][ab];
			  t_ijab = T2AB.matrix[Gji][ij][ab];
			}

			VAbc[Gab][ab][c] += t_kc * D_ijab + f_kc * t_ijab;
		      }

		      /* Sum V and W into V */
		      VAbc[Gab][ab][c] += WAbc[Gab][ab][c];

		      /* Build the rest of the denominator and divide it into W */
		      denom = dijk;
		      if(fAB.params->rowtot[Ga])
			denom -= fAB.matrix[Ga][a][a];
		      if(fab.params->rowtot[Gb])
			denom -= fab.matrix[Gb][b][b];
		      if(fab.params->rowtot[Gc])
			denom -= fab.matrix[Gc][c][c];

		      WAbc[Gab][ab][c] /= denom;

		    } /* c */
		  } /* ab */
		} /* Gab */

		/* 1/2 Dot product of final V and W is the energy for this ijk triple */
		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;
		  ET_ABB += dot_block(WAbc[Gab], VAbc[Gab], FABints.params->coltot[Gab], bvirtpi[Gc], 0.5);
		  global_dpd_->free_dpd_block(WAbc[Gab], FABints.params->coltot[Gab], bvirtpi[Gc]);
		  global_dpd_->free_dpd_block(VAbc[Gab], FABints.params->coltot[Gab], bvirtpi[Gc]);
		}

	      } /* J >= K */

	    } /* k */
	  } /* j */
	} /* i */

      } /* Gk */
    } /* Gj */
  } /* Gi */

  free(WAbc);
  free(VAbc);
  free(WAcb);
  free(WbcA);
  free(WcAb);
  free(WbAc);

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_close(&T2BB, h);
    global_dpd_->buf4_mat_irrep_close(&T2AB, h);
    global_dpd_->buf4_mat_irrep_close(&T2BA, h);
    global_dpd_->buf4_mat_irrep_close(&EBBints, h);
    global_dpd_->buf4_mat_irrep_close(&EABints, h);
    global_dpd_->buf4_mat_irrep_close(&EBAints, h);
    global_dpd_->buf4_mat_irrep_close(&DBBints, h);
    global_dpd_->buf4_mat_irrep_close(&DABints, h);
  }

  global_dpd_->buf4_close(&T2BB);
  global_dpd_->buf4_close(&T2AB);
  global_dpd_->buf4_close(&T2BA);
  global_dpd_->buf4_close(&FBBints);
  global_dpd_->buf4_close(&FABints);
  global_dpd_->buf4_close(&FBAints);
  global_dpd_->buf4_close(&EBBints);
  global_dpd_->buf4_close(&EABints);
  global_dpd_->buf4_close(&EBAints);
  global_dpd_->buf4_close(&DBBints);
  global_dpd_->buf4_close(&DABints);

  global_dpd_->file2_mat_close(&T1A);
  global_dpd_->file2_close(&T1A);
  global_dpd_->file2_mat_close(&T1B);
  global_dpd_->file2_close(&T1B);

  global_dpd_->file2_mat_close(&fIJ);
  global_dpd_->file2_mat_close(&fij);
  global_dpd_->file2_mat_close(&fAB);
  global_dpd_->file2_mat_close(&fab);
  global_dpd_->file2_mat_close(&fIA);
  global_dpd_->file2_mat_close(&fia);
  global_dpd_->file2_close(&fIJ);
  global_dpd_->file2_close(&fij);
  global_dpd_->file2_close(&fAB);
  global_dpd_->file2_close(&fab);
  global_dpd_->file2_close(&fIA);
  global_dpd_->file2_close(&fia);

  return ET_ABB;
}

}} // namespace psi::CCTRIPLES
