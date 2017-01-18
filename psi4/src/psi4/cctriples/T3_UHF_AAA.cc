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
    \ingroup DPD
    \brief Enter brief description of file here
*/

/* T3_UHF_AAA(): Computes all connected and disconnected T3(IJK,ABC)
** amplitudes for a given I, J, K combination for input C2, F, and E
** intermediates.  This function will work for AAA or BBB spin cases,
** with either RHF/ROHF or UHF orbitals.
**
** Arguments:
**
**   double ***W: The target connected triples amplitudes in an
**   nirreps x AB x C array.  The memory for this must be allocated
**   externally.
**
**   double ***V: The target disconnected triples amplitudes (if
**   requested) in an nirreps x AB x C array.  The memory for this
**   must be allocated externally.
**
**   int disc: Boolean: 1 == computed disconnected triples; 0 == don't
**
**   int nirreps: Number of irreps.
**
**   int I: Absolute index of orbital I.
**
**   int Gi: Irrep of I.
**
**   int J: Absolute index of orbital J.
**
**   int Gj: Irrep of J.
**
**   int K: Absolute index of orbital K.
**
**   int Gk: Irrep of K.
**
**   dpdbuf4 *C2: Pointer to dpd buffer for double excitation amps,
**   ordered (IJ,AB).
**
**   dpdbuf4 *F: Pointer to dpd buffer for three-virtual-index
**   intermediate, ordered (IA,BC).
**
**   dpdbuf4 *E: Pointer to dpd buffer for three-occupied-index
**   intermediate, ordered (IJ,KA).
**
**   dpdfile *C1: If disconnected T3's are requested, pointer to dpd
**   buffer for single-excitation amps.
**
**   dpdbuf4 *D: If disconnected T3's are requested, pointer to dpd
**   buffer for <IJ||ab> integrals.
**
**   dpdfile2 *fIA: Pointer to the dpd file2 for the occ-vir block of
**   the Fock matrix (or other appropriate one-electron operator).
**
**   dpdfile2 *fIJ: Pointer to the dpd file2 for the occ-occ block of
**   the Fock matrix (or other appropriate one-electron operator).
**
**   dpdfile2 *fAB: Pointer to the dpd file2 for the vir-vir block of
**   the Fock matrix (or other appropriate one-electron operator).
**
**   int *occpi: Number of occupied orbitals per irrep lookup array.
**
**   int *occ_off: Offset lookup for translating between absolute and
**   relative orbital indices for occupied space.
**
**   int *virtpi: Number of virtual orbitals per irrep lookup array.
**
**   int *vir_off: Offset lookup for translating between absolute and
**   relative orbital indices for virtual space.
**
**   double omega: a constant to add to the final denominators -
**   needed for CC3 EOM
**
** TDC, July 2004
** Modified to return disconnected triples, TDC, Feburary 2008
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "psi4/libqt/qt.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"

namespace psi { namespace cctriples {

    void T3_UHF_AAA(double ***W, double ***V, int disc, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk,
		    dpdbuf4 *C2, dpdbuf4 *F, dpdbuf4 *E, dpdfile2 *C1, dpdbuf4 *D, dpdfile2 *fIA, dpdfile2 *fIJ, dpdfile2 *fAB,
		    int *occpi, int *occ_off, int *virtpi, int *vir_off, double omega)
    {
      int h;
      int i, j, k;
      int ij, ji, ik, ki, jk, kj;
      int Gij, Gji, Gik, Gki, Gjk, Gkj, Gijk;
      int Ga, Gb, Gc;
      int Gd, Gl;
      int Gid, Gjd, Gkd;
      int Gab, Gcb, Gca, Gbc, Gac;
      int Gla, Glb, Glc;
      int Gil, Gjl, Gkl;
      int a, b, c, A, B, C;
      int ab, bc, ac;
      int cd, bd, ad;
      int id, jd, kd;
      int la, lb, lc;
      int il, jl, kl;
      int nrows, ncols, nlinks;
      double dijk, denom;
      double ***W2;
      int GE, GF, GC, GX3;
      double t_ia, t_ib, t_ic, t_ja, t_jb, t_jc, t_ka, t_kb, t_kc;
      double f_ia, f_ib, f_ic, f_ja, f_jb, f_jc, f_ka, f_kb, f_kc;
      double t_jkbc, t_jkac, t_jkba, t_ikbc, t_ikac, t_ikba, t_jibc, t_jiac, t_jiba;
      double D_jkbc, D_jkac, D_jkba, D_ikbc, D_ikac, D_ikba, D_jibc, D_jiac, D_jiba;

      GC = C2->file.my_irrep;
      /* F and E are assumed to have same irrep */
      GF = GE =  F->file.my_irrep;
      GX3 = GC^GF;

      global_dpd_->file2_mat_init(C1);
      global_dpd_->file2_mat_rd(C1);

      global_dpd_->file2_mat_init(fIJ);
      global_dpd_->file2_mat_init(fAB);
      global_dpd_->file2_mat_init(fIA);
      global_dpd_->file2_mat_rd(fIJ);
      global_dpd_->file2_mat_rd(fAB);
      global_dpd_->file2_mat_rd(fIA);

      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(C2, h);
	global_dpd_->buf4_mat_irrep_rd(C2, h);

	if(disc) {
	  global_dpd_->buf4_mat_irrep_init(D, h);
	  global_dpd_->buf4_mat_irrep_rd(D, h);
	}

	global_dpd_->buf4_mat_irrep_init(E, h);
	global_dpd_->buf4_mat_irrep_rd(E, h);
      }

      i = I - occ_off[Gi];
      j = J - occ_off[Gj];
      k = K - occ_off[Gk];

      Gij = Gji = Gi ^ Gj;
      Gik = Gki = Gi ^ Gk;
      Gjk = Gkj = Gj ^ Gk;
      Gijk = Gi ^ Gj ^ Gk;

      ij = C2->params->rowidx[I][J];
      ji = C2->params->rowidx[J][I];
      jk = C2->params->rowidx[J][K];
      kj = C2->params->rowidx[K][J];
      ik = C2->params->rowidx[I][K];
      ki = C2->params->rowidx[K][I];

      dijk = 0.0;
      if(fIJ->params->rowtot[Gi]) dijk += fIJ->matrix[Gi][i][i];
      if(fIJ->params->rowtot[Gj]) dijk += fIJ->matrix[Gj][j][j];
      if(fIJ->params->rowtot[Gk]) dijk += fIJ->matrix[Gk][k][k];

      W2 = (double ***) malloc(nirreps * sizeof(double **));

      for(Gab=0; Gab < nirreps; Gab++) {
	Gc = Gab ^ Gijk ^ GX3; /* changed */

	W2[Gab] = global_dpd_->dpd_block_matrix(F->params->coltot[Gab], virtpi[Gc]);

	if(F->params->coltot[Gab] && virtpi[Gc]) {
	  memset(W[Gab][0], 0, F->params->coltot[Gab]*virtpi[Gc]*sizeof(double));
	  if(disc) memset(V[Gab][0], 0, F->params->coltot[Gab]*virtpi[Gc]*sizeof(double));
	}
      }

      for(Gd=0; Gd < nirreps; Gd++) {

	/* +t_kjcd * F_idab */
	Gid = Gi ^ Gd;
	Gab = Gid ^ GF; /* changed */

	Gc = Gjk ^ Gd ^ GC; /* changed */

	cd = C2->col_offset[Gjk][Gc];
	id = F->row_offset[Gid][I];

	F->matrix[Gid] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gid^GF]);
	global_dpd_->buf4_mat_irrep_rd_block(F, Gid, id, virtpi[Gd]);

	nrows = F->params->coltot[Gid^GF];
	ncols = virtpi[Gc];
	nlinks = virtpi[Gd];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gid][0], nrows,
		  &(C2->matrix[Gjk][kj][cd]), nlinks, 1.0, W[Gab][0], ncols);

	global_dpd_->free_dpd_block(F->matrix[Gid], virtpi[Gd], F->params->coltot[Gid^GF]);

	/* +t_ikcd * F_jdab */
	Gjd = Gj ^ Gd;
	Gab = Gjd ^ GF; /* changed */
	Gc = Gik ^ Gd ^ GC;  /* changed */

	cd = C2->col_offset[Gik][Gc];
	jd = F->row_offset[Gjd][J];

	F->matrix[Gjd] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gjd^GF]);
	global_dpd_->buf4_mat_irrep_rd_block(F, Gjd, jd, virtpi[Gd]);

	nrows = F->params->coltot[Gjd^GF];
	ncols = virtpi[Gc];
	nlinks = virtpi[Gd];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gjd][0], nrows,
		  &(C2->matrix[Gik][ik][cd]), nlinks, 1.0, W[Gab][0], ncols);

	global_dpd_->free_dpd_block(F->matrix[Gjd], virtpi[Gd], F->params->coltot[Gjd^GF]);

	/* -t_ijcd * F_kdab */
	Gkd = Gk ^ Gd; /*changed */
	Gab = Gkd ^ GF;
	Gc = Gij ^ Gd ^ GC;

	cd = C2->col_offset[Gij][Gc];
	kd = F->row_offset[Gkd][K];

	F->matrix[Gkd] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gkd^GF]);
	global_dpd_->buf4_mat_irrep_rd_block(F, Gkd, kd, virtpi[Gd]);

	nrows = F->params->coltot[Gkd^GF];
	ncols = virtpi[Gc];
	nlinks = virtpi[Gd];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, F->matrix[Gkd][0], nrows,
		  &(C2->matrix[Gij][ij][cd]), nlinks, 1.0, W[Gab][0], ncols);

	global_dpd_->free_dpd_block(F->matrix[Gkd], virtpi[Gd], F->params->coltot[Gkd^GF]);

      }

      for(Gl=0; Gl < nirreps; Gl++) {

	/* -t_ilab * E_jklc */
	Gil = Gi ^ Gl; /* changed */
	Gab = Gil ^ GC;
	Gc = Gjk ^ Gl ^ GE;

	lc = E->col_offset[Gjk][Gl];
	il = C2->row_offset[Gil][I];

	nrows = C2->params->coltot[Gil^GC];
	ncols = virtpi[Gc];
	nlinks = occpi[Gl];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, C2->matrix[Gil][il], nrows,
		  &(E->matrix[Gjk][jk][lc]), ncols, 1.0, W[Gab][0], ncols);

	/* +t_jlab * E_iklc */
	Gjl = Gj ^ Gl; /* changed */
	Gab = Gjl ^ GC;
	Gc = Gik ^ Gl ^ GE;

	lc = E->col_offset[Gik][Gl];
	jl = C2->row_offset[Gjl][J];

	nrows = C2->params->coltot[Gjl^GC];
	ncols = virtpi[Gc];
	nlinks = occpi[Gl];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, C2->matrix[Gjl][jl], nrows,
		  &(E->matrix[Gik][ik][lc]), ncols, 1.0, W[Gab][0], ncols);

	/* +t_klab * E_jilc */
	Gkl = Gk ^ Gl; /* changed! */
	Gab = Gkl ^ GC;
	Gc = Gji ^ Gl ^ GE;

	lc = E->col_offset[Gji][Gl];
	kl = C2->row_offset[Gkl][K];

	nrows = C2->params->coltot[Gkl^GC];
	ncols = virtpi[Gc];
	nlinks = occpi[Gl];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, C2->matrix[Gkl][kl], nrows,
		  &(E->matrix[Gji][ji][lc]), ncols, 1.0, W[Gab][0], ncols);

      }

      for(Gd=0; Gd < nirreps; Gd++) {
	/* +t_kjbd * F_idca */
	Gid = Gi ^ Gd; /* changed */
	Gca = Gid ^ GF;
	Gb = Gjk ^ Gd ^ GC;

	bd = C2->col_offset[Gjk][Gb];
	id = F->row_offset[Gid][I];

	F->matrix[Gid] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gid^GF]);
	global_dpd_->buf4_mat_irrep_rd_block(F, Gid, id, virtpi[Gd]);

	nrows = F->params->coltot[Gid^GF];
	ncols = virtpi[Gb];
	nlinks = virtpi[Gd];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gid][0], nrows,
		  &(C2->matrix[Gjk][kj][bd]), nlinks, 1.0, W2[Gca][0], ncols);

	global_dpd_->free_dpd_block(F->matrix[Gid], virtpi[Gd], F->params->coltot[Gid^GF]);

	/* +t_ikbd * F_jdca */
	Gjd = Gj ^ Gd;
	Gca = Gjd ^ GF ;
	Gb = Gik ^ Gd ^ GC;

	bd = C2->col_offset[Gik][Gb];
	jd = F->row_offset[Gjd][J];

	F->matrix[Gjd] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gjd^GF]);
	global_dpd_->buf4_mat_irrep_rd_block(F, Gjd, jd, virtpi[Gd]);

	nrows = F->params->coltot[Gjd^GF];
	ncols = virtpi[Gb];
	nlinks = virtpi[Gd];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gjd][0], nrows,
		  &(C2->matrix[Gik][ik][bd]), nlinks, 1.0, W2[Gca][0], ncols);

	global_dpd_->free_dpd_block(F->matrix[Gjd], virtpi[Gd], F->params->coltot[Gjd^GF]);

	/* -t_ijbd * F_kdca */
	Gkd = Gk ^ Gd;
	Gca = Gkd ^ GF;
	Gb = Gij ^ Gd ^ GC;

	bd = C2->col_offset[Gij][Gb];
	kd = F->row_offset[Gkd][K];

	F->matrix[Gkd] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gkd^GF]);
	global_dpd_->buf4_mat_irrep_rd_block(F, Gkd, kd, virtpi[Gd]);

	nrows = F->params->coltot[Gkd^GF];
	ncols = virtpi[Gb];
	nlinks = virtpi[Gd];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t','t',nrows, ncols, nlinks, -1.0, F->matrix[Gkd][0], nrows,
		  &(C2->matrix[Gij][ij][bd]), nlinks, 1.0, W2[Gca][0], ncols);

	global_dpd_->free_dpd_block(F->matrix[Gkd], virtpi[Gd], F->params->coltot[Gkd^GF]);
      }

      for(Gl=0; Gl < nirreps; Gl++) {
	/* -t_ilca * E_jklb */
	Gil = Gi ^ Gl;
	Gca = Gil ^ GC;
	Gb = Gjk ^ Gl ^ GE;

	lb = E->col_offset[Gjk][Gl];
	il = C2->row_offset[Gil][I];

	nrows = C2->params->coltot[Gil^GC];
	ncols = virtpi[Gb];
	nlinks = occpi[Gl];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, C2->matrix[Gil][il], nrows,
		  &(E->matrix[Gjk][jk][lb]), ncols, 1.0, W2[Gca][0], ncols);

	/* +t_jlca * E_iklb */
	Gjl = Gj ^ Gl;
	Gca = Gjl ^ GC;
	Gb = Gik ^ Gl^ GE;

	lb = E->col_offset[Gik][Gl];
	jl = C2->row_offset[Gjl][J];

	nrows = C2->params->coltot[Gjl^GC];
	ncols = virtpi[Gb];
	nlinks = occpi[Gl];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, C2->matrix[Gjl][jl], nrows,
		  &(E->matrix[Gik][ik][lb]), ncols, 1.0, W2[Gca][0], ncols);

	/* +t_klca * E_jilb */
	Gkl = Gk ^ Gl;
	Gca = Gkl ^ GC;
	Gb = Gji ^ Gl ^ GE;

	lb = E->col_offset[Gji][Gl];
	kl = C2->row_offset[Gkl][K];

	nrows = C2->params->coltot[Gkl^GC];
	ncols = virtpi[Gb];
	nlinks = occpi[Gl];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, C2->matrix[Gkl][kl], nrows,
		  &(E->matrix[Gji][ji][lb]), ncols, 1.0, W2[Gca][0], ncols);
      }

      global_dpd_->sort_3d(W2, W, nirreps, Gijk^GX3, F->params->coltot, F->params->colidx,
		  F->params->colorb, F->params->rsym, F->params->ssym, vir_off,
		  vir_off, virtpi, vir_off, F->params->colidx, bca, 1);

      for(Gab=0; Gab < nirreps; Gab++) {
	Gc = Gab ^ Gijk ^ GX3; /* changed */
	if(F->params->coltot[Gab] && virtpi[Gc]) {
	  memset(W2[Gab][0], 0, F->params->coltot[Gab]*virtpi[Gc]*sizeof(double));
	}
      }

      for(Gd=0; Gd < nirreps; Gd++) {
	/* -t_kjad * F_idcb */
	Gid = Gi ^ Gd;
	Gcb = Gid ^ GF;
	Ga = Gkj ^ Gd ^ GC;

	ad = C2->col_offset[Gkj][Ga];
	id = F->row_offset[Gid][I];

	F->matrix[Gid] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gid^GF]);
	global_dpd_->buf4_mat_irrep_rd_block(F, Gid, id, virtpi[Gd]);

	nrows = F->params->coltot[Gid^GF];
	ncols = virtpi[Ga];
	nlinks = virtpi[Gd];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t','t',nrows, ncols, nlinks, -1.0, F->matrix[Gid][0], nrows,
		  &(C2->matrix[Gkj][kj][ad]), nlinks, 1.0, W2[Gcb][0], ncols);

	global_dpd_->free_dpd_block(F->matrix[Gid], virtpi[Gd], F->params->coltot[Gid^GF]);

	/* -t_ikad * F_jdcb */
	Gjd = Gj ^ Gd;
	Gcb = Gjd ^ GF;
	Ga = Gik ^ Gd ^ GC;

	ad = C2->col_offset[Gik][Ga];
	jd = F->row_offset[Gjd][J];

	F->matrix[Gjd] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gjd^GF]);
	global_dpd_->buf4_mat_irrep_rd_block(F, Gjd, jd, virtpi[Gd]);

	nrows = F->params->coltot[Gjd^GF];
	ncols = virtpi[Ga];
	nlinks = virtpi[Gd];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t','t',nrows, ncols, nlinks, -1.0, F->matrix[Gjd][0], nrows,
		  &(C2->matrix[Gik][ik][ad]), nlinks, 1.0, W2[Gcb][0], ncols);

	global_dpd_->free_dpd_block(F->matrix[Gjd], virtpi[Gd], F->params->coltot[Gjd^GF]);

	/* +t_ijad * F_kdcb */
	Gkd = Gk ^ Gd;
	Gcb = Gkd ^ GF;
	Ga = Gij ^ Gd ^ GC;

	ad = C2->col_offset[Gij][Ga];
	kd = F->row_offset[Gkd][K];

	F->matrix[Gkd] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gkd^GF]);
	global_dpd_->buf4_mat_irrep_rd_block(F, Gkd, kd, virtpi[Gd]);

	nrows = F->params->coltot[Gkd^GF];
	ncols = virtpi[Ga];
	nlinks = virtpi[Gd];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gkd][0], nrows,
		  &(C2->matrix[Gij][ij][ad]), nlinks, 1.0, W2[Gcb][0], ncols);

	global_dpd_->free_dpd_block(F->matrix[Gkd], virtpi[Gd], F->params->coltot[Gkd^GF]);

      }

      for(Gl=0; Gl < nirreps; Gl++) {
	/* +t_ilcb * E_jkla */
	Gil = Gi ^ Gl;
	Gcb  = Gil ^ GC;
	Ga = Gjk ^ Gl ^ GE;

	la = E->col_offset[Gjk][Gl];
	il = C2->row_offset[Gil][I];

	nrows = C2->params->coltot[Gil^GC];
	ncols = virtpi[Ga];
	nlinks = occpi[Gl];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, C2->matrix[Gil][il], nrows,
		  &(E->matrix[Gjk][jk][la]), ncols, 1.0, W2[Gcb][0], ncols);

	/* -t_jlcb * E_ikla */
	Gjl = Gj ^ Gl;
	Gcb = Gjl ^ GC;
	Ga = Gik ^ Gl ^ GE;

	la = E->col_offset[Gik][Gl];
	jl = C2->row_offset[Gjl][J];

	nrows = C2->params->coltot[Gjl^GC];
	ncols = virtpi[Ga];
	nlinks = occpi[Gl];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, C2->matrix[Gjl][jl], nrows,
		  &(E->matrix[Gik][ik][la]), ncols, 1.0, W2[Gcb][0], ncols);

	/* -t_klcb * E_jila */
	Gkl = Gk ^ Gl;
	Gcb = Gkl ^ GC;
	Ga = Gji ^ Gl ^ GE;

	la = E->col_offset[Gji][Gl];
	kl = C2->row_offset[Gkl][K];

	nrows = C2->params->coltot[Gkl^GC];
	ncols = virtpi[Ga];
	nlinks = occpi[Gl];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, C2->matrix[Gkl][kl], nrows,
		  &(E->matrix[Gji][ji][la]), ncols, 1.0, W2[Gcb][0], ncols);

      }

      global_dpd_->sort_3d(W2, W, nirreps, Gijk^GX3, F->params->coltot, F->params->colidx,
		  F->params->colorb, F->params->rsym, F->params->ssym, vir_off,
		  vir_off, virtpi, vir_off, F->params->colidx, cba, 1);

      /**** Compute disconnected T3s for given ijk ****/
      if(disc) {
	for(Gab=0; Gab < nirreps; Gab++) {
	  Gc = Gab ^ Gijk;

	  for(ab=0; ab < F->params->coltot[Gab]; ab++) {
	    A = F->params->colorb[Gab][ab][0];
	    Ga = F->params->rsym[A];
	    a = A - vir_off[Ga];
	    B = F->params->colorb[Gab][ab][1];
	    Gb = F->params->ssym[B];
	    b = B - vir_off[Gb];

	    Gbc = Gb ^ Gc;
	    Gac = Ga ^ Gc;

	    for(c=0; c < virtpi[Gc]; c++) {
	      C = vir_off[Gc] + c;

	      bc = D->params->colidx[B][C];
	      ac = D->params->colidx[A][C];

	      /* +t_ia * D_jkbc + f_ia * t_jkbc */
	      if(Gi == Ga && Gjk == Gbc) {
		t_ia = D_jkbc = f_ia = t_jkbc = 0.0;

		if(C1->params->rowtot[Gi] && C1->params->coltot[Gi]) {
		  t_ia = C1->matrix[Gi][i][a];
		  f_ia = fIA->matrix[Gi][i][a];
		}

		if(D->params->rowtot[Gjk] && D->params->coltot[Gjk]) {
		  D_jkbc = D->matrix[Gjk][jk][bc];
		  t_jkbc = C2->matrix[Gjk][jk][bc];
		}

		V[Gab][ab][c] += t_ia * D_jkbc + f_ia * t_jkbc;
	      }

	      /* -t_ib * D_jkac - f_ib * t_jkac */
	      if(Gi == Gb && Gjk == Gac) {
		t_ib = D_jkac = f_ib = t_jkac = 0.0;

		if(C1->params->rowtot[Gi] && C1->params->coltot[Gi]) {
		  t_ib = C1->matrix[Gi][i][b];
		  f_ib = fIA->matrix[Gi][i][b];
		}

		if(D->params->rowtot[Gjk] && D->params->coltot[Gjk]) {
		  D_jkac = D->matrix[Gjk][jk][ac];
		  t_jkac = C2->matrix[Gjk][jk][ac];
		}

		V[Gab][ab][c] -= t_ib * D_jkac + f_ib * t_jkac;
	      }

	      /* +t_ic * D_jkab + f_ic * t_jkba */
	      if(Gi == Gc && Gjk == Gab) {
		t_ic = D_jkba = f_ic = t_jkba = 0.0;

		if(C1->params->rowtot[Gi] && C1->params->coltot[Gi]) {
		  t_ic = C1->matrix[Gi][i][c];
		  f_ic = fIA->matrix[Gi][i][c];
		}

		if(D->params->rowtot[Gjk] && D->params->coltot[Gjk]) {
		  D_jkba = D->matrix[Gjk][jk][ab];
		  t_jkba = C2->matrix[Gjk][jk][ab];
		}

		V[Gab][ab][c] += t_ic * D_jkba + f_ic * t_jkba;
	      }

	      /* -t_ja * D_ikbc - f_ja * t_ikbc*/
	      if(Gj == Ga && Gik == Gbc) {
		t_ja = D_ikbc = f_ja = t_ikbc = 0.0;

		if(C1->params->rowtot[Gj] && C1->params->coltot[Gj]) {
		  t_ja = C1->matrix[Gj][j][a];
		  f_ja = fIA->matrix[Gj][j][a];
		}

		if(D->params->rowtot[Gik] && D->params->coltot[Gik]) {
		  D_ikbc = D->matrix[Gik][ik][bc];
		  t_ikbc = C2->matrix[Gik][ik][bc];
		}

		V[Gab][ab][c] -= t_ja * D_ikbc + f_ja * t_ikbc;
	      }

	      /* +t_jb * D_ikac + f_jb * t_ikac */
	      if(Gj == Gb && Gik == Gac) {
		t_jb = D_ikac = f_jb = t_ikac = 0.0;

		if(C1->params->rowtot[Gj] && C1->params->coltot[Gj]) {
		  t_jb = C1->matrix[Gj][j][b];
		  f_jb = fIA->matrix[Gj][j][b];
		}

		if(D->params->rowtot[Gik] && D->params->coltot[Gik]) {
		  D_ikac = D->matrix[Gik][ik][ac];
		  t_ikac = C2->matrix[Gik][ik][ac];
		}

		V[Gab][ab][c] += t_jb * D_ikac + f_jb * t_ikac;
	      }

	      /* -t_jc * D_ikba - f_jc * t_ikba */
	      if(Gj == Gc && Gik == Gab) {
		t_jc = D_ikba = f_jc = t_ikba = 0.0;

		if(C1->params->rowtot[Gj] && C1->params->coltot[Gj]) {
		  t_jc = C1->matrix[Gj][j][c];
		  f_jc = fIA->matrix[Gj][j][c];
		}

		if(D->params->rowtot[Gik] && D->params->coltot[Gik]) {
		  D_ikba = D->matrix[Gik][ik][ab];
		  t_ikba = C2->matrix[Gik][ik][ab];
		}

		V[Gab][ab][c] -= t_jc * D_ikba + f_jc * t_ikba;
	      }

	      /* -t_ka * D_jibc - f_ka * t_jibc */
	      if(Gk == Ga && Gji == Gbc) {
		t_ka = D_jibc = f_ka = t_jibc = 0.0;

		if(C1->params->rowtot[Gk] && C1->params->coltot[Gk]) {
		  t_ka = C1->matrix[Gk][k][a];
		  f_ka = fIA->matrix[Gk][k][a];
		}

		if(D->params->rowtot[Gji] && D->params->coltot[Gji]) {
		  D_jibc = D->matrix[Gji][ji][bc];
		  t_jibc = C2->matrix[Gji][ji][bc];
		}

		V[Gab][ab][c] -= t_ka * D_jibc + f_ka * t_jibc;
	      }

	      /* +t_kb * D_jiac + f_kb * t_jiac */
	      if(Gk == Gb && Gji == Gac) {
		t_kb = D_jiac = f_kb = t_jiac = 0.0;

		if(C1->params->rowtot[Gk] && C1->params->coltot[Gk]) {
		  t_kb = C1->matrix[Gk][k][b];
		  f_kb = fIA->matrix[Gk][k][b];
		}

		if(D->params->rowtot[Gji] && D->params->coltot[Gji]) {
		  D_jiac = D->matrix[Gji][ji][ac];
		  t_jiac = C2->matrix[Gji][ji][ac];
		}

		V[Gab][ab][c] += t_kb * D_jiac + f_kb * t_jiac;
	      }

	      /* -t_kc * D_jiab - f_kc * t_jiba*/
	      if(Gk == Gc && Gji == Gab) {
		t_kc = D_jiba = f_kc = t_jiba = 0.0;

		if(C1->params->rowtot[Gk] && C1->params->coltot[Gk]) {
		  t_kc = C1->matrix[Gk][k][c];
		  f_kc = fIA->matrix[Gk][k][c];
		}

		if(D->params->rowtot[Gji] && D->params->coltot[Gji]) {
		  D_jiba = D->matrix[Gji][ji][ab];
		  t_jiba = C2->matrix[Gji][ji][ab];
		}

		V[Gab][ab][c] -= t_kc * D_jiba + f_kc * t_jiba;
	      }

	    } /* c */
	  } /* ab */
	} /* Gab */


      } /**** Disconnected T3 complete ****/

      for(Gab=0; Gab < nirreps; Gab++) {
	Gc = Gab ^ Gijk ^ GX3; /* assumes totally symmetric! */
	for(ab=0; ab < F->params->coltot[Gab]; ab++) {
	  A = F->params->colorb[Gab][ab][0];
	  Ga = F->params->rsym[A];
	  a = A - vir_off[Ga];
	  B = F->params->colorb[Gab][ab][1];
	  Gb = F->params->ssym[B];
	  b = B - vir_off[Gb];

	  for(c=0; c < virtpi[Gc]; c++) {
	    denom = dijk;
	    if(fAB->params->rowtot[Ga]) denom -= fAB->matrix[Ga][a][a];
	    if(fAB->params->rowtot[Gb]) denom -= fAB->matrix[Gb][b][b];
	    if(fAB->params->rowtot[Gc]) denom -= fAB->matrix[Gc][c][c];

	    W[Gab][ab][c] /= (omega + denom);
	    if(disc) V[Gab][ab][c] /= (omega + denom);

	  } /* c */
	} /* ab */
      } /* Gab */

      for(Gab=0; Gab < nirreps; Gab++) {
	Gc = Gab ^ Gijk ^ GX3; /* changed */
	global_dpd_->free_dpd_block(W2[Gab], F->params->coltot[Gab], virtpi[Gc]);
      }
      free(W2);

      global_dpd_->file2_mat_close(fIJ);
      global_dpd_->file2_mat_close(fAB);
      if(disc) {
	global_dpd_->file2_mat_close(fIA);
	global_dpd_->file2_mat_close(C1);
      }

      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_close(C2, h);
	if(disc) global_dpd_->buf4_mat_irrep_close(D, h);
	global_dpd_->buf4_mat_irrep_close(E, h);
      }
    }

  }}
