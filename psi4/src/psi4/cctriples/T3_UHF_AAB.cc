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

/* T3_UHF_AAB(): Computes all T3(IJK,ABC) amplitudes for a given I, J,
** K combination for input T2, F, and E intermediates.  This function
** will work for AAB or BBA spin cases with either RHF/ROHF or UHF
** orbitals.
**
** Arguments:
**
**   double ***W: The target triples amplitudes in an nirreps x AB x C
**   array.  The memory for this must be allocated externally.
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
**   dpdbuf4 *T2: Pointer to dpd buffer for double excitation amps,
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
**   double omega: constant to add to denominators - needed for
**   CC3 EOM
**
** TDC, July 2004
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "psi4/libqt/qt.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"

namespace psi { namespace cctriples {

    void T3_UHF_AAB(double ***W, double ***V, int disc, int nirreps,
		    int I, int Gi, int J, int Gj, int K, int Gk,
		    dpdbuf4 *T2AA, dpdbuf4 *T2AB, dpdbuf4 *T2BA, dpdbuf4 *FAA, dpdbuf4 *FAB, dpdbuf4 *FBA,
		    dpdbuf4 *EAA, dpdbuf4 *EAB, dpdbuf4 *EBA, dpdfile2 *T1A, dpdfile2 *T1B,
		    dpdbuf4 *DAA, dpdbuf4 *DAB, dpdfile2 *fIA, dpdfile2 *fia,
		    dpdfile2 *fIJ, dpdfile2 *fij,dpdfile2 *fAB, dpdfile2 *fab,
		    int *aoccpi, int *aocc_off, int *boccpi, int *bocc_off,
		    int *avirtpi, int *avir_off, int *bvirtpi, int *bvir_off, double omega)
    {
      int h;
      int i, j, k;
      int ij, ji, ik, ki, jk, kj;
      int Gij, Gji, Gik, Gki, Gjk, Gkj, Gijk;
      int Ga, Gb, Gc;
      int Gd, Gl;
      int Gid, Gjd, Gkd;
      int Gab, Gcb, Gca, Gac, Gbc;
      int Gla, Glb, Glc;
      int Gil, Gjl, Gkl;
      int a, b, c, A, B, C;
      int ab, bc, ac;
      int dc, bd, ad;
      int id, jd, kd;
      int la, lb, lc;
      int il, jl, kl;
      int nrows, ncols, nlinks;
      double dijk, denom;
      double ***W2;
      int GE, GF, GC, GX3;
      double t_ia, t_ib, t_ja, t_jb, t_kc;
      double f_ia, f_ib, f_ja, f_jb, f_kc;
      double D_jkbc, D_jkac, D_ikbc, D_ikac, D_jiab;
      double t_jkbc, t_jkac, t_ikbc, t_ikac, t_jiab;

      GC = T2AA->file.my_irrep;
      /* F and E are assumed to have same irrep */
      GF = GE =  FAA->file.my_irrep;
      GX3 = GC^GF;

      global_dpd_->file2_mat_init(fIJ);
      global_dpd_->file2_mat_init(fAB);
      global_dpd_->file2_mat_init(fij);
      global_dpd_->file2_mat_init(fab);
      if(disc) {
	global_dpd_->file2_mat_init(fIA);
	global_dpd_->file2_mat_init(fia);
	global_dpd_->file2_mat_init(T1A);
	global_dpd_->file2_mat_init(T1B);
      }

      global_dpd_->file2_mat_rd(fIJ);
      global_dpd_->file2_mat_rd(fAB);
      global_dpd_->file2_mat_rd(fij);
      global_dpd_->file2_mat_rd(fab);
      if(disc) {
	global_dpd_->file2_mat_rd(fIA);
	global_dpd_->file2_mat_rd(fia);
	global_dpd_->file2_mat_rd(T1A);
	global_dpd_->file2_mat_rd(T1B);
      }

      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(T2AA, h);
	global_dpd_->buf4_mat_irrep_rd(T2AA, h);

	global_dpd_->buf4_mat_irrep_init(T2AB, h);
	global_dpd_->buf4_mat_irrep_rd(T2AB, h);

	global_dpd_->buf4_mat_irrep_init(T2BA, h);
	global_dpd_->buf4_mat_irrep_rd(T2BA, h);

	if(disc) {
	  global_dpd_->buf4_mat_irrep_init(DAA, h);
	  global_dpd_->buf4_mat_irrep_rd(DAA, h);

	  global_dpd_->buf4_mat_irrep_init(DAB, h);
	  global_dpd_->buf4_mat_irrep_rd(DAB, h);
	}

	global_dpd_->buf4_mat_irrep_init(EAA, h);
	global_dpd_->buf4_mat_irrep_rd(EAA, h);

	global_dpd_->buf4_mat_irrep_init(EAB, h);
	global_dpd_->buf4_mat_irrep_rd(EAB, h);

	global_dpd_->buf4_mat_irrep_init(EBA, h);
	global_dpd_->buf4_mat_irrep_rd(EBA, h);
      }

      i = I - aocc_off[Gi];
      j = J - aocc_off[Gj];
      k = K - bocc_off[Gk];

      Gij = Gji = Gi ^ Gj;
      Gik = Gki = Gi ^ Gk;
      Gjk = Gkj = Gj ^ Gk;
      Gijk = Gi ^ Gj ^ Gk;

      ij = T2AA->params->rowidx[I][J];
      ji = T2AA->params->rowidx[J][I];
      jk = T2AB->params->rowidx[J][K];
      kj = T2BA->params->rowidx[K][J];
      ik = T2AB->params->rowidx[I][K];
      ki = T2BA->params->rowidx[K][I];

      dijk = 0.0;
      if(fIJ->params->rowtot[Gi]) dijk += fIJ->matrix[Gi][i][i];
      if(fIJ->params->rowtot[Gj]) dijk += fIJ->matrix[Gj][j][j];
      if(fij->params->rowtot[Gk]) dijk += fij->matrix[Gk][k][k];

      W2 = (double ***) malloc(nirreps * sizeof(double **)); /* alpha-beta-alpha */

      /* clear out the old W */
      for(Gab=0; Gab < nirreps; Gab++) {
	Gc = Gab ^ Gijk ^ GX3; /* assumes totally symmetric! */

	if(FAA->params->coltot[Gab] && bvirtpi[Gc]) {
	  memset(W[Gab][0], 0, FAA->params->coltot[Gab]*bvirtpi[Gc]*sizeof(double));
	  memset(V[Gab][0], 0, FAA->params->coltot[Gab]*bvirtpi[Gc]*sizeof(double));
	}
      }

      for(Gd=0; Gd < nirreps; Gd++) {
	/* +C_JkDc * F_IDAB */
	Gid = Gi ^ Gd;
	Gab = Gid ^ GF;
	Gc = Gjk ^ Gd ^ GC;

	dc = T2AB->col_offset[Gjk][Gd];
	id = FAA->row_offset[Gid][I];

	FAA->matrix[Gid] = global_dpd_->dpd_block_matrix(avirtpi[Gd], FAA->params->coltot[Gid^GF]);
	global_dpd_->buf4_mat_irrep_rd_block(FAA, Gid, id, avirtpi[Gd]);

	nrows = FAA->params->coltot[Gid^GF];
	ncols = bvirtpi[Gc];
	nlinks = avirtpi[Gd];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t','n',nrows, ncols, nlinks, 1.0, FAA->matrix[Gid][0], nrows,
		  &(T2AB->matrix[Gjk][jk][dc]), ncols, 1.0, W[Gab][0], ncols);

	global_dpd_->free_dpd_block(FAA->matrix[Gid], avirtpi[Gd], FAA->params->coltot[Gid^GF]);

	/* -C_IkDc * F_JDAB */
	Gjd = Gj ^ Gd;
	Gab = Gjd ^ GF;
	Gc = Gik ^ Gd ^ GC;

	dc = T2AB->col_offset[Gik][Gd];
	jd = FAA->row_offset[Gjd][J];

	FAA->matrix[Gjd] = global_dpd_->dpd_block_matrix(avirtpi[Gd], FAA->params->coltot[Gjd^GF]);
	global_dpd_->buf4_mat_irrep_rd_block(FAA, Gjd, jd, avirtpi[Gd]);

	nrows = FAA->params->coltot[Gjd^GF];
	ncols = bvirtpi[Gc];
	nlinks = avirtpi[Gd];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t','n',nrows, ncols, nlinks, -1.0, FAA->matrix[Gjd][0], nrows,
		  &(T2AB->matrix[Gik][ik][dc]), ncols, 1.0, W[Gab][0], ncols);

	global_dpd_->free_dpd_block(FAA->matrix[Gjd], avirtpi[Gd], FAA->params->coltot[Gjd^GF]);
      }
      for(Gl=0; Gl < nirreps; Gl++) {

	/* -C_ILAB * E_JkLc */
	Gil = Gi ^ Gl;
	Gab = Gil ^ GC;
	Gc = Gjk ^ Gl ^ GE;

	lc = EAB->col_offset[Gjk][Gl];
	il = T2AA->row_offset[Gil][I];

	nrows = T2AA->params->coltot[Gil^GC];
	ncols = bvirtpi[Gc];
	nlinks = aoccpi[Gl];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, T2AA->matrix[Gil][il], nrows,
		  &(EAB->matrix[Gjk][jk][lc]), ncols, 1.0, W[Gab][0], ncols);

	/* +t_JLAB * E_IkLc */
	Gjl = Gj ^ Gl;
	Gab = Gjl ^ GC;
	Gc = Gik ^ Gl ^ GE;

	lc = EAB->col_offset[Gik][Gl];
	jl = T2AA->row_offset[Gjl][J];

	nrows = T2AA->params->coltot[Gjl^GC];
	ncols = bvirtpi[Gc];
	nlinks = aoccpi[Gl];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, T2AA->matrix[Gjl][jl], nrows,
		  &(EAB->matrix[Gik][ik][lc]), ncols, 1.0, W[Gab][0], ncols);
      }

      /* Open memory for an alpha-beta-alpha array */
      for(Gab=0; Gab < nirreps; Gab++) {
	Gc = Gab ^ Gijk ^ GX3;

	W2[Gab] = global_dpd_->dpd_block_matrix(FAB->params->coltot[Gab], avirtpi[Gc]);
      }

      for(Gd=0; Gd < nirreps; Gd++) {
	/* +t_JkBd * F_IdAc */
	Gid = Gi ^ Gd;
	Gac = Gid ^ GF;
	Gb = Gjk ^ Gd ^ GC;

	bd = T2AB->col_offset[Gjk][Gb];
	id = FAB->row_offset[Gid][I];

	FAB->matrix[Gid] = global_dpd_->dpd_block_matrix(bvirtpi[Gd], FAB->params->coltot[Gid^GF]);
	global_dpd_->buf4_mat_irrep_rd_block(FAB, Gid, id, bvirtpi[Gd]);

	nrows = FAB->params->coltot[Gid^GF];
	ncols = avirtpi[Gb];
	nlinks = bvirtpi[Gd];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, FAB->matrix[Gid][0], nrows,
		  &(T2AB->matrix[Gjk][jk][bd]), nlinks, 1.0, W2[Gac][0], ncols);

	global_dpd_->free_dpd_block(FAB->matrix[Gid], bvirtpi[Gd], FAB->params->coltot[Gid^GF]);

	/* -t_IkBd * F_JdAc */
	Gjd = Gj ^ Gd;
	Gac = Gjd ^ GF;
	Gb = Gik ^ Gd ^ GC;

	bd = T2AB->col_offset[Gik][Gb];
	jd = FAB->row_offset[Gjd][J];

	FAB->matrix[Gjd] = global_dpd_->dpd_block_matrix(bvirtpi[Gd], FAB->params->coltot[Gjd^GF]);
	global_dpd_->buf4_mat_irrep_rd_block(FAB, Gjd, jd, bvirtpi[Gd]);

	nrows = FAB->params->coltot[Gjd^GF];
	ncols = avirtpi[Gb];
	nlinks = bvirtpi[Gd];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t','t',nrows, ncols, nlinks, -1.0, FAB->matrix[Gjd][0], nrows,
		  &(T2AB->matrix[Gik][ik][bd]), nlinks, 1.0, W2[Gac][0], ncols);

	global_dpd_->free_dpd_block(FAB->matrix[Gjd], bvirtpi[Gd], FAB->params->coltot[Gjd^GF]);

      }

      for(Gl=0; Gl < nirreps; Gl++) {
	/* -t_IlAc * E_kJlB */
	Gil = Gi ^ Gl;
	Gac = Gil ^ GC;
	Gb = Gkj ^ Gl ^ GE;

	lb = EBA->col_offset[Gkj][Gl];
	il = T2AB->row_offset[Gil][I];

	nrows = T2AB->params->coltot[Gil^GC];
	ncols = avirtpi[Gb];
	nlinks = boccpi[Gl];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, T2AB->matrix[Gil][il], nrows,
		  &(EBA->matrix[Gkj][kj][lb]), ncols, 1.0, W2[Gac][0], ncols);

	/* +t_JlAc * E_kIlB */
	Gjl = Gj ^ Gl;
	Gac = Gjl ^ GC;
	Gb = Gki ^ Gl ^ GE;

	lb = EBA->col_offset[Gki][Gl];
	jl = T2AB->row_offset[Gjl][J];

	nrows = T2AB->params->coltot[Gjl^GC];
	ncols = avirtpi[Gb];
	nlinks = boccpi[Gl];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, T2AB->matrix[Gjl][jl], nrows,
		  &(EBA->matrix[Gki][ki][lb]), ncols, 1.0, W2[Gac][0], ncols);
      }

      /* W(Ac,B) --> W(AB,c) */
      global_dpd_->sort_3d(W2, W, nirreps, Gijk^GX3, FAB->params->coltot, FAB->params->colidx,
		  FAB->params->colorb, FAB->params->rsym, FAB->params->ssym, avir_off,
		  bvir_off, avirtpi, avir_off, FAA->params->colidx, acb, 1);

      /* clean out the alpha-beta-alpha intermediate for next set of terms */
      for(Gab=0; Gab < nirreps; Gab++) {
	Gc = Gab ^ Gijk ^ GX3;
	if(FAB->params->coltot[Gab] && avirtpi[Gc]) {
	  memset(W2[Gab][0], 0, FAB->params->coltot[Gab]*avirtpi[Gc]*sizeof(double));
	}
      }

      for(Gd=0; Gd < nirreps; Gd++) {
	/* -C_JkAd * F_IdBc */
	Gid = Gi ^ Gd;
	Gbc = Gid ^ GF;
	Ga = Gjk ^ Gd ^ GC;

	ad = T2AB->col_offset[Gjk][Ga];
	id = FAB->row_offset[Gid][I];

	FAB->matrix[Gid] = global_dpd_->dpd_block_matrix(bvirtpi[Gd], FAB->params->coltot[Gid^GF]);
	global_dpd_->buf4_mat_irrep_rd_block(FAB, Gid, id, bvirtpi[Gd]);

	nrows = FAB->params->coltot[Gid^GF];
	ncols = avirtpi[Ga];
	nlinks = bvirtpi[Gd];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t','t',nrows, ncols, nlinks, -1.0, FAB->matrix[Gid][0], nrows,
		  &(T2AB->matrix[Gjk][jk][ad]), nlinks, 1.0, W2[Gbc][0], ncols);

	global_dpd_->free_dpd_block(FAB->matrix[Gid], bvirtpi[Gd], FAB->params->coltot[Gid^GF]);

	/* +t_IkAd * F_JdBc */
	Gjd = Gj ^ Gd;
	Gbc = Gjd ^ GF;
	Ga = Gik ^ Gd ^ GC;

	ad = T2AB->col_offset[Gik][Ga];
	jd = FAB->row_offset[Gjd][J];

	FAB->matrix[Gjd] = global_dpd_->dpd_block_matrix(bvirtpi[Gd], FAB->params->coltot[Gjd^GF]);
	global_dpd_->buf4_mat_irrep_rd_block(FAB, Gjd, jd, bvirtpi[Gd]);

	nrows = FAB->params->coltot[Gjd^GF];
	ncols = avirtpi[Ga];
	nlinks = bvirtpi[Gd];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, FAB->matrix[Gjd][0], nrows,
		  &(T2AB->matrix[Gik][ik][ad]), nlinks, 1.0, W2[Gbc][0], ncols);

	global_dpd_->free_dpd_block(FAB->matrix[Gjd], bvirtpi[Gd], FAB->params->coltot[Gjd^GF]);

      }

      for(Gl=0; Gl < nirreps; Gl++) {
	/* +C_IlBc * E_kJlA */
	Gil = Gi ^ Gl;
	Gbc  = Gil ^ GC;
	Ga = Gkj ^ Gl ^ GE;

	la = EBA->col_offset[Gkj][Gl];
	il = T2AB->row_offset[Gil][I];

	nrows = T2AB->params->coltot[Gil^GC];
	ncols = avirtpi[Ga];
	nlinks = boccpi[Gl];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, T2AB->matrix[Gil][il], nrows,
		  &(EBA->matrix[Gkj][kj][la]), ncols, 1.0, W2[Gbc][0], ncols);

	/* -C_JlBc * E_kIlA */
	Gjl = Gj ^ Gl;
	Gbc = Gjl ^ GC;
	Ga = Gki ^ Gl ^ GE;

	la = EBA->col_offset[Gki][Gl];
	jl = T2AB->row_offset[Gjl][J];

	nrows = T2AB->params->coltot[Gjl^GC];
	ncols = avirtpi[Ga];
	nlinks = boccpi[Gl];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, T2AB->matrix[Gjl][jl], nrows,
		  &(EBA->matrix[Gki][ki][la]), ncols, 1.0, W2[Gbc][0], ncols);
      }

      global_dpd_->sort_3d(W2, W, nirreps, Gijk^GX3, FAB->params->coltot, FAB->params->colidx,
		  FAB->params->colorb, FAB->params->rsym, FAB->params->ssym, avir_off,
		  bvir_off, avirtpi, avir_off, FAA->params->colidx, cab, 1);
      /* Close the alpha-beta-alpha array and open a beta-alpha-alpha array */
      for(Gab=0; Gab < nirreps; Gab++) {
	Gc = Gab ^ Gijk^GX3; /* assumes totally symmetric! */

	global_dpd_->free_dpd_block(W2[Gab], FAB->params->coltot[Gab], avirtpi[Gc]);
	W2[Gab] = global_dpd_->dpd_block_matrix(FBA->params->coltot[Gab], avirtpi[Gc]);
      }

      /* Insert cBA terms */
      for(Gd=0; Gd < nirreps; Gd++) {
	/* -C_JIAD * F_kDcB */
	Gkd = Gk ^ Gd;
	Gcb = Gkd ^ GF;
	Ga = Gji ^ Gd ^ GC;

	ad = T2AA->col_offset[Gji][Ga];
	kd = FBA->row_offset[Gkd][K];

	FBA->matrix[Gkd] = global_dpd_->dpd_block_matrix(avirtpi[Gd], FBA->params->coltot[Gkd^GF]);
	global_dpd_->buf4_mat_irrep_rd_block(FBA, Gkd, kd, avirtpi[Gd]);

	nrows = FBA->params->coltot[Gkd^GF];
	ncols = avirtpi[Ga];
	nlinks = avirtpi[Gd];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, FBA->matrix[Gkd][0], nrows,
		  &(T2AA->matrix[Gji][ji][ad]), nlinks, 1.0, W2[Gcb][0], ncols);

	global_dpd_->free_dpd_block(FBA->matrix[Gkd], avirtpi[Gd], FBA->params->coltot[Gkd^GF]);
      }

      for(Gl=0; Gl < nirreps; Gl++) {

	/* -C_kLcB * E_JILA */
	Gkl = Gk ^ Gl;
	Gcb = Gkl ^ GC;
	Ga = Gji ^ Gl ^ GE;

	la = EAA->col_offset[Gji][Gl];
	kl = T2BA->row_offset[Gkl][K];

	nrows = T2BA->params->coltot[Gkl^GC];
	ncols = avirtpi[Ga];
	nlinks = aoccpi[Gl];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, T2BA->matrix[Gkl][kl], nrows,
		  &(EAA->matrix[Gji][ji][la]), ncols, 1.0, W2[Gcb][0], ncols);
      }

      global_dpd_->sort_3d(W2, W, nirreps, Gijk^GX3, FBA->params->coltot, FBA->params->colidx,
		  FBA->params->colorb, FBA->params->rsym, FBA->params->ssym, bvir_off,
		  avir_off, avirtpi, avir_off, FAA->params->colidx, cba, 1);

      /* clean out the beta-alpha-alpha intermediate for next set of terms */
      for(Gab=0; Gab < nirreps; Gab++) {
	Gc = Gab ^ Gijk ^ GX3;
	if(FBA->params->coltot[Gab] && avirtpi[Gc]) {
	  memset(W2[Gab][0], 0, FBA->params->coltot[Gab]*avirtpi[Gc]*sizeof(double));
	}
      }

      for(Gd=0; Gd < nirreps; Gd++) {

	/* +t_JIBD * F_kDcA */
	Gkd = Gk ^ Gd;
	Gca = Gkd ^ GF;
	Gb = Gji ^ Gd ^ GC;

	bd = T2AA->col_offset[Gji][Gb];
	kd = FBA->row_offset[Gkd][K];

	FBA->matrix[Gkd] = global_dpd_->dpd_block_matrix(avirtpi[Gd], FBA->params->coltot[Gkd^GF]);
	global_dpd_->buf4_mat_irrep_rd_block(FBA, Gkd, kd, avirtpi[Gd]);

	nrows = FBA->params->coltot[Gkd^GF];
	ncols = avirtpi[Gb];
	nlinks = avirtpi[Gd];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, FBA->matrix[Gkd][0], nrows,
		  &(T2AA->matrix[Gji][ji][bd]), nlinks, 1.0, W2[Gca][0], ncols);

	global_dpd_->free_dpd_block(FBA->matrix[Gkd], avirtpi[Gd], FBA->params->coltot[Gkd^GF]);
      }

      for(Gl=0; Gl < nirreps; Gl++) {

	/* +C_kLcA * E_JILB */
	Gkl = Gk ^ Gl;
	Gca = Gkl ^ GC;
	Gb = Gji ^ Gl ^ GE;

	lb = EAA->col_offset[Gji][Gl];
	kl = T2BA->row_offset[Gkl][K];

	nrows = T2BA->params->coltot[Gkl^GC];
	ncols = avirtpi[Gb];
	nlinks = aoccpi[Gl];

	if(nrows && ncols && nlinks)
	  C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, T2BA->matrix[Gkl][kl], nrows,
		  &(EAA->matrix[Gji][ji][lb]), ncols, 1.0, W2[Gca][0], ncols);
      }

      global_dpd_->sort_3d(W2, W, nirreps, Gijk^GX3, FBA->params->coltot, FBA->params->colidx,
		  FBA->params->colorb, FBA->params->rsym, FBA->params->ssym,
		  bvir_off, avir_off, avirtpi, avir_off, FAA->params->colidx, bca, 1);

      /*** compute disconnected triples ***/
      if(disc) {
	for(Gab=0; Gab < nirreps; Gab++) {
	  Gc = Gab ^ Gijk;

	  for(ab=0; ab < FAA->params->coltot[Gab]; ab++) {
	    A = FAA->params->colorb[Gab][ab][0];
	    Ga = FAA->params->rsym[A];
	    a = A - avir_off[Ga];
	    B = FAA->params->colorb[Gab][ab][1];
	    Gb = FAA->params->ssym[B];
	    b = B - avir_off[Gb];

	    Gbc = Gb ^ Gc;
	    Gac = Ga ^ Gc;

	    for(c=0; c < bvirtpi[Gc]; c++) {
	      C = bvir_off[Gc] + c;

	      bc = DAB->params->colidx[B][C];
	      ac = DAB->params->colidx[A][C];

	      /* +t_IA * D_JkBc + f_IA * t_JkBc */
	      if(Gi == Ga && Gjk == Gbc) {
		t_ia = D_jkbc = f_ia = t_jkbc = 0.0;

		if(T1A->params->rowtot[Gi] && T1A->params->coltot[Gi]) {
		  t_ia = T1A->matrix[Gi][i][a];
		  f_ia = fIA->matrix[Gi][i][a];
		}

		if(DAB->params->rowtot[Gjk] && DAB->params->coltot[Gjk]) {
		  D_jkbc = DAB->matrix[Gjk][jk][bc];
		  t_jkbc = T2AB->matrix[Gjk][jk][bc];
		}

		V[Gab][ab][c] += t_ia * D_jkbc + f_ia * t_jkbc;
	      }

	      /* -t_IB * D_JkAc - f_IB * t_JkAc */
	      if(Gi == Gb && Gjk == Gac) {
		t_ib = D_jkac = f_ib = t_jkac = 0.0;

		if(T1A->params->rowtot[Gi] && T1A->params->coltot[Gi]) {
		  t_ib = T1A->matrix[Gi][i][b];
		  f_ib = fIA->matrix[Gi][i][b];
		}

		if(DAB->params->rowtot[Gjk] && DAB->params->coltot[Gjk]) {
		  D_jkac = DAB->matrix[Gjk][jk][ac];
		  t_jkac = T2AB->matrix[Gjk][jk][ac];
		}

		V[Gab][ab][c] -= t_ib * D_jkac + f_ib * t_jkac;
	      }

	      /* -t_JA * D_IkBc - f_JA * t_IkBc */
	      if(Gj == Ga && Gik == Gbc) {
		t_ja = D_ikbc = f_ja = t_ikbc = 0.0;

		if(T1A->params->rowtot[Gj] && T1A->params->coltot[Gj]) {
		  t_ja = T1A->matrix[Gj][j][a];
		  f_ja = fIA->matrix[Gj][j][a];
		}

		if(DAB->params->rowtot[Gik] && DAB->params->coltot[Gik]) {
		  D_ikbc = DAB->matrix[Gik][ik][bc];
		  t_ikbc = T2AB->matrix[Gik][ik][bc];
		}

		V[Gab][ab][c] -= t_ja * D_ikbc + f_ja * t_ikbc;
	      }

	      /* +t_JB * D_IkAc + f_JB * t_IkAc */
	      if(Gj == Gb && Gik == Gac) {
		t_jb = D_ikac = f_jb = t_ikac = 0.0;

		if(T1A->params->rowtot[Gj] && T1A->params->coltot[Gj]) {
		  t_jb = T1A->matrix[Gj][j][b];
		  f_jb = fIA->matrix[Gj][j][b];
		}

		if(DAB->params->rowtot[Gik] && DAB->params->coltot[Gik]) {
		  D_ikac = DAB->matrix[Gik][ik][ac];
		  t_ikac = T2AB->matrix[Gik][ik][ac];
		}

		V[Gab][ab][c] += t_jb * D_ikac + f_jb * t_ikac;
	      }

	      /* -t_kc * D_JIAB - f_kc * t_JIAB */
	      if(Gk == Gc && Gji == Gab) {
		t_kc = D_jiab = f_kc = t_jiab = 0.0;

		if(T1B->params->rowtot[Gk] && T1B->params->coltot[Gk]) {
		  t_kc = T1B->matrix[Gk][k][c];
		  f_kc = fia->matrix[Gk][k][c];
		}

		if(DAA->params->rowtot[Gji] && DAA->params->coltot[Gji]) {
		  D_jiab = DAA->matrix[Gji][ji][ab];
		  t_jiab = T2AA->matrix[Gji][ji][ab];
		}

		V[Gab][ab][c] -= t_kc * D_jiab + f_kc * t_jiab;
	      }

	    } /* c */
	  } /* ab */
	} /* Gab */
      }
      /*** disconnected triples complete ***/

      for(Gab=0; Gab < nirreps; Gab++) {
	Gc = Gab ^ Gijk ^ GX3; /* assumes totally symmetric! */

	for(ab=0; ab < FAA->params->coltot[Gab]; ab++) {
	  A = FAA->params->colorb[Gab][ab][0];
	  B = FAA->params->colorb[Gab][ab][1];
	  Ga = FAA->params->rsym[A];
	  Gb = FAA->params->ssym[B];
	  a = A - avir_off[Ga];
	  b = B - avir_off[Gb];

	  for(c=0; c < bvirtpi[Gc]; c++) {
	    C = bvir_off[Gc] + c;

	    denom = dijk;
	    if(fAB->params->rowtot[Ga]) denom -= fAB->matrix[Ga][a][a];
	    if(fAB->params->rowtot[Gb]) denom -= fAB->matrix[Gb][b][b];
	    if(fab->params->rowtot[Gc]) denom -= fab->matrix[Gc][c][c];

	    W[Gab][ab][c] /= (denom + omega);
	    if(disc) V[Gab][ab][c] /= (denom + omega);

	  } /* c */
	} /* ab */
      } /* Gab */

      for(Gab=0; Gab < nirreps; Gab++) {
	Gc = Gab ^ Gijk ^ GX3;
	global_dpd_->free_dpd_block(W2[Gab], FBA->params->coltot[Gab], avirtpi[Gc]);
      }

      free(W2);

      global_dpd_->file2_mat_close(fIJ);
      global_dpd_->file2_mat_close(fij);
      global_dpd_->file2_mat_close(fAB);
      global_dpd_->file2_mat_close(fab);
      if(disc) {
	global_dpd_->file2_mat_close(fIA);
	global_dpd_->file2_mat_close(fia);
	global_dpd_->file2_mat_close(T1A);
	global_dpd_->file2_mat_close(T1B);
      }

      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_close(T2AA, h);
	global_dpd_->buf4_mat_irrep_close(T2AB, h);
	global_dpd_->buf4_mat_irrep_close(T2BA, h);

	global_dpd_->buf4_mat_irrep_close(EAA, h);
	global_dpd_->buf4_mat_irrep_close(EAB, h);
	global_dpd_->buf4_mat_irrep_close(EBA, h);

	if(disc) {
	  global_dpd_->buf4_mat_irrep_close(DAA, h);
	  global_dpd_->buf4_mat_irrep_close(DAB, h);
	}
      }
    }

  }}
