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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/

/* L3_RHF_AAA(): Computes all T3(IJK,ABC) amplitudes for a given I, J,
** K combination for input T2, F, and E intermediates.  This function
** will work for AAA or BBB spin cases, with either RHF/ROHF or UHF
** orbitals.
**
** Arguments:
**
**   double ***W1: The target triples amplitudes in an nirreps x AB x
**   C array.  The memory for this must be allocated externally.
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
** TDC, July 2004
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"

namespace psi { namespace cclambda {

void L3_AAA(double ***W1, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk,
	    dpdbuf4 *T2, dpdbuf4 *F, dpdbuf4 *E, dpdfile2 *fIJ, dpdfile2 *fAB,
	    dpdbuf4 *D, dpdbuf4 *LIJAB, dpdfile2 *LIA, dpdfile2 *FME,
	    int *occpi, int *occ_off, int *virtpi, int *vir_off)
{
  int h;
  int i, j, k;
  int Ga, Gb, Gc;
  int Gij, Gji, Gik, Gki, Gjk, Gkj, Gijk;
  int Gab, Gba, Gbc, Gcb, Gac, Gca;
  int Gd, Gl;
  int Gid, Gjd, Gkd;
  int Gla, Glb, Glc;
  int Gil, Gjl, Gkl;
  int a, b, c, A, B, C;
  int ij, ji, ik, ki, jk, kj;
  int ab, ba, ac, ca, bc, cb;
  int cd, bd, ad;
  int id, jd, kd;
  int la, lb, lc;
  int il, jl, kl;
  int nrows, ncols, nlinks;
  double dijk, denom;
  double ***W2;
  double L1, F1, D2, L2;

  global_dpd_->file2_mat_init(fIJ);
  global_dpd_->file2_mat_init(fAB);
  global_dpd_->file2_mat_rd(fIJ);
  global_dpd_->file2_mat_rd(fAB);

  global_dpd_->file2_mat_init(LIA);
  global_dpd_->file2_mat_rd(LIA);
  global_dpd_->file2_mat_init(FME);
  global_dpd_->file2_mat_rd(FME);

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(T2, h);
    global_dpd_->buf4_mat_irrep_rd(T2, h);

    global_dpd_->buf4_mat_irrep_init(E, h);
    global_dpd_->buf4_mat_irrep_rd(E, h);

    global_dpd_->buf4_mat_irrep_init(LIJAB, h);
    global_dpd_->buf4_mat_irrep_rd(LIJAB, h);

    global_dpd_->buf4_mat_irrep_init(D, h);
    global_dpd_->buf4_mat_irrep_rd(D, h);
  }

  i = I - occ_off[Gi];
  j = J - occ_off[Gj];
  k = K - occ_off[Gk];

  Gij = Gji = Gi ^ Gj;
  Gik = Gki = Gi ^ Gk;
  Gjk = Gkj = Gj ^ Gk;
  Gijk = Gi ^ Gj ^ Gk;

  ij = T2->params->rowidx[I][J];
  ji = T2->params->rowidx[J][I];
  jk = T2->params->rowidx[J][K];
  kj = T2->params->rowidx[K][J];
  ik = T2->params->rowidx[I][K];
  ki = T2->params->rowidx[K][I];

  dijk = 0.0;
  if(fIJ->params->rowtot[Gi]) dijk += fIJ->matrix[Gi][i][i];
  if(fIJ->params->rowtot[Gj]) dijk += fIJ->matrix[Gj][j][j];
  if(fIJ->params->rowtot[Gk]) dijk += fIJ->matrix[Gk][k][k];

  W2 = (double ***) malloc(nirreps * sizeof(double **));

  for(Gab=0; Gab < nirreps; Gab++) {
    Gc = Gab ^ Gijk; /* assumes totally symmetric! */

    W2[Gab] = global_dpd_->dpd_block_matrix(F->params->coltot[Gab], virtpi[Gc]);

    if(F->params->coltot[Gab] && virtpi[Gc]) {
      memset(W1[Gab][0], 0, F->params->coltot[Gab]*virtpi[Gc]*sizeof(double));
    }
  }

  for(Gd=0; Gd < nirreps; Gd++) {

    /* +t_kjcd * F_idab */
    Gab = Gid = Gi ^ Gd; /* assumes totally symmetric! */
    Gc = Gjk ^ Gd;       /* assumes totally symmetric! */

    cd = T2->col_offset[Gjk][Gc];
    id = F->row_offset[Gid][I];

    F->matrix[Gid] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gid]);
    global_dpd_->buf4_mat_irrep_rd_block(F, Gid, id, virtpi[Gd]);

    nrows = F->params->coltot[Gid];
    ncols = virtpi[Gc];
    nlinks = virtpi[Gd];

    if(nrows && ncols && nlinks)
      C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gid][0], nrows,
	      &(T2->matrix[Gjk][kj][cd]), nlinks, 1.0, W1[Gab][0], ncols);

    global_dpd_->free_dpd_block(F->matrix[Gid], virtpi[Gd], F->params->coltot[Gid]);

    /* +t_ikcd * F_jdab */
    Gab = Gjd = Gj ^ Gd; /* assumes totally symmetric! */
    Gc = Gik ^ Gd;       /* assumes totally symmetric! */

    cd = T2->col_offset[Gik][Gc];
    jd = F->row_offset[Gjd][J];

    F->matrix[Gjd] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gjd]);
    global_dpd_->buf4_mat_irrep_rd_block(F, Gjd, jd, virtpi[Gd]);

    nrows = F->params->coltot[Gjd];
    ncols = virtpi[Gc];
    nlinks = virtpi[Gd];

    if(nrows && ncols && nlinks)
      C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gjd][0], nrows,
	      &(T2->matrix[Gik][ik][cd]), nlinks, 1.0, W1[Gab][0], ncols);

    global_dpd_->free_dpd_block(F->matrix[Gjd], virtpi[Gd], F->params->coltot[Gjd]);

    /* -t_ijcd * F_kdab */
    Gab = Gkd = Gk ^ Gd; /* assumes totally symmetric! */
    Gc = Gij ^ Gd;       /* assumes totally symmetric! */

    cd = T2->col_offset[Gij][Gc];
    kd = F->row_offset[Gkd][K];

    F->matrix[Gkd] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gkd]);
    global_dpd_->buf4_mat_irrep_rd_block(F, Gkd, kd, virtpi[Gd]);

    nrows = F->params->coltot[Gkd];
    ncols = virtpi[Gc];
    nlinks = virtpi[Gd];

    if(nrows && ncols && nlinks)
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, F->matrix[Gkd][0], nrows,
	      &(T2->matrix[Gij][ij][cd]), nlinks, 1.0, W1[Gab][0], ncols);

    global_dpd_->free_dpd_block(F->matrix[Gkd], virtpi[Gd], F->params->coltot[Gkd]);

  }

  for(Gl=0; Gl < nirreps; Gl++) {

    /* -t_ilab * E_jklc */
    Gab = Gil = Gi ^ Gl; /* assumes totally symmetric! */
    Gc = Gjk ^ Gl;       /* assumes totally symmetric! */

    lc = E->col_offset[Gjk][Gl];
    il = T2->row_offset[Gil][I];

    nrows = T2->params->coltot[Gil];
    ncols = virtpi[Gc];
    nlinks = occpi[Gl];

    if(nrows && ncols && nlinks)
      C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, T2->matrix[Gil][il], nrows,
	      &(E->matrix[Gjk][jk][lc]), ncols, 1.0, W1[Gab][0], ncols);

    /* +t_jlab * E_iklc */
    Gab = Gjl = Gj ^ Gl; /* assumes totally symmetric! */
    Gc = Gik ^ Gl;       /* assumes totally symmetric! */

    lc = E->col_offset[Gik][Gl];
    jl = T2->row_offset[Gjl][J];

    nrows = T2->params->coltot[Gjl];
    ncols = virtpi[Gc];
    nlinks = occpi[Gl];

    if(nrows && ncols && nlinks)
      C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, T2->matrix[Gjl][jl], nrows,
	      &(E->matrix[Gik][ik][lc]), ncols, 1.0, W1[Gab][0], ncols);

    /* +t_klab * E_jilc */
    Gab = Gkl = Gk ^ Gl; /* assumes totally symmetric! */
    Gc = Gji ^ Gl;       /* assumes totally symmetric! */

    lc = E->col_offset[Gji][Gl];
    kl = T2->row_offset[Gkl][K];

    nrows = T2->params->coltot[Gkl];
    ncols = virtpi[Gc];
    nlinks = occpi[Gl];

    if(nrows && ncols && nlinks)
      C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, T2->matrix[Gkl][kl], nrows,
	      &(E->matrix[Gji][ji][lc]), ncols, 1.0, W1[Gab][0], ncols);

  }

  for(Gd=0; Gd < nirreps; Gd++) {
    /* +t_kjbd * F_idca */
    Gca = Gid = Gi ^ Gd; /* assumes totally symmetric! */
    Gb = Gjk ^ Gd;       /* assumes totally symmetric! */

    bd = T2->col_offset[Gjk][Gb];
    id = F->row_offset[Gid][I];

    F->matrix[Gid] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gid]);
    global_dpd_->buf4_mat_irrep_rd_block(F, Gid, id, virtpi[Gd]);

    nrows = F->params->coltot[Gid];
    ncols = virtpi[Gb];
    nlinks = virtpi[Gd];

    if(nrows && ncols && nlinks)
      C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gid][0], nrows,
	      &(T2->matrix[Gjk][kj][bd]), nlinks, 1.0, W2[Gca][0], ncols);

    global_dpd_->free_dpd_block(F->matrix[Gid], virtpi[Gd], F->params->coltot[Gid]);

    /* +t_ikbd * F_jdca */
    Gca = Gjd = Gj ^ Gd; /* assumes totally symmetric! */
    Gb = Gik ^ Gd;       /* assumes totally symmetric! */

    bd = T2->col_offset[Gik][Gb];
    jd = F->row_offset[Gjd][J];

    F->matrix[Gjd] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gjd]);
    global_dpd_->buf4_mat_irrep_rd_block(F, Gjd, jd, virtpi[Gd]);

    nrows = F->params->coltot[Gjd];
    ncols = virtpi[Gb];
    nlinks = virtpi[Gd];

    if(nrows && ncols && nlinks)
      C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gjd][0], nrows,
	      &(T2->matrix[Gik][ik][bd]), nlinks, 1.0, W2[Gca][0], ncols);

    global_dpd_->free_dpd_block(F->matrix[Gjd], virtpi[Gd], F->params->coltot[Gjd]);

    /* -t_ijbd * F_kdca */
    Gca = Gkd = Gk ^ Gd; /* assumes totally symmetric! */
    Gb = Gij ^ Gd;       /* assumes totally symmetric! */

    bd = T2->col_offset[Gij][Gb];
    kd = F->row_offset[Gkd][K];

    F->matrix[Gkd] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gkd]);
    global_dpd_->buf4_mat_irrep_rd_block(F, Gkd, kd, virtpi[Gd]);

    nrows = F->params->coltot[Gkd];
    ncols = virtpi[Gb];
    nlinks = virtpi[Gd];

    if(nrows && ncols && nlinks)
      C_DGEMM('t','t',nrows, ncols, nlinks, -1.0, F->matrix[Gkd][0], nrows,
	      &(T2->matrix[Gij][ij][bd]), nlinks, 1.0, W2[Gca][0], ncols);

    global_dpd_->free_dpd_block(F->matrix[Gkd], virtpi[Gd], F->params->coltot[Gkd]);
  }

  for(Gl=0; Gl < nirreps; Gl++) {
    /* -t_ilca * E_jklb */
    Gca = Gil = Gi ^ Gl; /* assumes totally symmetric! */
    Gb = Gjk ^ Gl;       /* assumes totally symmetric! */

    lb = E->col_offset[Gjk][Gl];
    il = T2->row_offset[Gil][I];

    nrows = T2->params->coltot[Gil];
    ncols = virtpi[Gb];
    nlinks = occpi[Gl];

    if(nrows && ncols && nlinks)
      C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, T2->matrix[Gil][il], nrows,
	      &(E->matrix[Gjk][jk][lb]), ncols, 1.0, W2[Gca][0], ncols);

    /* +t_jlca * E_iklb */
    Gca = Gjl = Gj ^ Gl; /* assumes totally symmetric! */
    Gb = Gik ^ Gl;       /* assumes totally symmetric! */

    lb = E->col_offset[Gik][Gl];
    jl = T2->row_offset[Gjl][J];

    nrows = T2->params->coltot[Gjl];
    ncols = virtpi[Gb];
    nlinks = occpi[Gl];

    if(nrows && ncols && nlinks)
      C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, T2->matrix[Gjl][jl], nrows,
	      &(E->matrix[Gik][ik][lb]), ncols, 1.0, W2[Gca][0], ncols);

    /* +t_klca * E_jilb */
    Gca = Gkl = Gk ^ Gl; /* assumes totally symmetric! */
    Gb = Gji ^ Gl;       /* assumes totally symmetric! */

    lb = E->col_offset[Gji][Gl];
    kl = T2->row_offset[Gkl][K];

    nrows = T2->params->coltot[Gkl];
    ncols = virtpi[Gb];
    nlinks = occpi[Gl];

    if(nrows && ncols && nlinks)
      C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, T2->matrix[Gkl][kl], nrows,
	      &(E->matrix[Gji][ji][lb]), ncols, 1.0, W2[Gca][0], ncols);
  }

  global_dpd_->sort_3d(W2, W1, nirreps, Gijk, F->params->coltot, F->params->colidx,
	      F->params->colorb, F->params->rsym, F->params->ssym, vir_off,
	      vir_off, virtpi, vir_off, F->params->colidx, bca, 1);

  for(Gab=0; Gab < nirreps; Gab++) {
    Gc = Gab ^ Gijk; /* assumes totally symmetric! */
    if(F->params->coltot[Gab] && virtpi[Gc]) {
      memset(W2[Gab][0], 0, F->params->coltot[Gab]*virtpi[Gc]*sizeof(double));
    }
  }

  for(Gd=0; Gd < nirreps; Gd++) {
    /* -t_kjad * F_idcb */
    Gcb = Gid = Gi ^ Gd; /* assumes totally symmetric! */
    Ga = Gkj ^ Gd;       /* assumes totally symmetric! */

    ad = T2->col_offset[Gkj][Ga];
    id = F->row_offset[Gid][I];

    F->matrix[Gid] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gid]);
    global_dpd_->buf4_mat_irrep_rd_block(F, Gid, id, virtpi[Gd]);

    nrows = F->params->coltot[Gid];
    ncols = virtpi[Ga];
    nlinks = virtpi[Gd];

    if(nrows && ncols && nlinks)
      C_DGEMM('t','t',nrows, ncols, nlinks, -1.0, F->matrix[Gid][0], nrows,
	      &(T2->matrix[Gkj][kj][ad]), nlinks, 1.0, W2[Gcb][0], ncols);

    global_dpd_->free_dpd_block(F->matrix[Gid], virtpi[Gd], F->params->coltot[Gid]);

    /* -t_ikad * F_jdcb */
    Gcb = Gjd = Gj ^ Gd; /* assumes totally symmetric! */
    Ga = Gik ^ Gd;       /* assumes totally symmetric! */

    ad = T2->col_offset[Gik][Ga];
    jd = F->row_offset[Gjd][J];

    F->matrix[Gjd] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gjd]);
    global_dpd_->buf4_mat_irrep_rd_block(F, Gjd, jd, virtpi[Gd]);

    nrows = F->params->coltot[Gjd];
    ncols = virtpi[Ga];
    nlinks = virtpi[Gd];

    if(nrows && ncols && nlinks)
      C_DGEMM('t','t',nrows, ncols, nlinks, -1.0, F->matrix[Gjd][0], nrows,
	      &(T2->matrix[Gik][ik][ad]), nlinks, 1.0, W2[Gcb][0], ncols);

    global_dpd_->free_dpd_block(F->matrix[Gjd], virtpi[Gd], F->params->coltot[Gjd]);

    /* +t_ijad * F_kdcb */
    Gcb = Gkd = Gk ^ Gd; /* assumes totally symmetric! */
    Ga = Gij ^ Gd;       /* assumes totally symmetric! */

    ad = T2->col_offset[Gij][Ga];
    kd = F->row_offset[Gkd][K];

    F->matrix[Gkd] = global_dpd_->dpd_block_matrix(virtpi[Gd], F->params->coltot[Gkd]);
    global_dpd_->buf4_mat_irrep_rd_block(F, Gkd, kd, virtpi[Gd]);

    nrows = F->params->coltot[Gkd];
    ncols = virtpi[Ga];
    nlinks = virtpi[Gd];

    if(nrows && ncols && nlinks)
      C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gkd][0], nrows,
	      &(T2->matrix[Gij][ij][ad]), nlinks, 1.0, W2[Gcb][0], ncols);

    global_dpd_->free_dpd_block(F->matrix[Gkd], virtpi[Gd], F->params->coltot[Gkd]);

  }

  for(Gl=0; Gl < nirreps; Gl++) {
    /* +t_ilcb * E_jkla */
    Gcb  = Gil = Gi ^ Gl; /* assumes totally symmetric! */
    Ga = Gjk ^ Gl;       /* assumes totally symmetric! */

    la = E->col_offset[Gjk][Gl];
    il = T2->row_offset[Gil][I];

    nrows = T2->params->coltot[Gil];
    ncols = virtpi[Ga];
    nlinks = occpi[Gl];

    if(nrows && ncols && nlinks)
      C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, T2->matrix[Gil][il], nrows,
	      &(E->matrix[Gjk][jk][la]), ncols, 1.0, W2[Gcb][0], ncols);

    /* -t_jlcb * E_ikla */
    Gcb  = Gjl = Gj ^ Gl; /* assumes totally symmetric! */
    Ga = Gik ^ Gl;       /* assumes totally symmetric! */

    la = E->col_offset[Gik][Gl];
    jl = T2->row_offset[Gjl][J];

    nrows = T2->params->coltot[Gjl];
    ncols = virtpi[Ga];
    nlinks = occpi[Gl];

    if(nrows && ncols && nlinks)
      C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, T2->matrix[Gjl][jl], nrows,
	      &(E->matrix[Gik][ik][la]), ncols, 1.0, W2[Gcb][0], ncols);

    /* -t_klcb * E_jila */
    Gcb  = Gkl = Gk ^ Gl; /* assumes totally symmetric! */
    Ga = Gji ^ Gl;       /* assumes totally symmetric! */

    la = E->col_offset[Gji][Gl];
    kl = T2->row_offset[Gkl][K];

    nrows = T2->params->coltot[Gkl];
    ncols = virtpi[Ga];
    nlinks = occpi[Gl];

    if(nrows && ncols && nlinks)
      C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, T2->matrix[Gkl][kl], nrows,
	      &(E->matrix[Gji][ji][la]), ncols, 1.0, W2[Gcb][0], ncols);

  }

  global_dpd_->sort_3d(W2, W1, nirreps, Gijk, F->params->coltot, F->params->colidx,
	      F->params->colorb, F->params->rsym, F->params->ssym, vir_off,
	      vir_off, virtpi, vir_off, F->params->colidx, cba, 1);

  /* Add disconnected terms */
  for(Gab=0; Gab < nirreps; Gab++) {
    Gc = Gab ^ Gijk;
    Gba = Gab;
    for(ab=0; ab < F->params->coltot[Gab]; ab++) {
      A = F->params->colorb[Gab][ab][0];
      B = F->params->colorb[Gab][ab][1];
      Ga = F->params->rsym[A];
      Gb = F->params->ssym[B];
      a = A - vir_off[Ga];
      b = B - vir_off[Gb];

      Gbc = Gcb = Gb ^ Gc;
      Gac = Gca = Ga ^ Gc;

      ba = LIJAB->params->colidx[B][A];

      for(c=0; c < virtpi[Gc]; c++) {
	C = vir_off[Gc] + c;

	bc = LIJAB->params->colidx[B][C];
	cb = LIJAB->params->colidx[C][B];
	ac = LIJAB->params->colidx[A][C];
	ca = LIJAB->params->colidx[C][A];

	/* +L_ia * D_jkbc + F_ia * L_jkbc */
	if(Gi == Ga && Gjk == Gbc) {
	  L1 = D2 = F1 = L2 = 0.0;
	  if(LIA->params->rowtot[Gi] && LIA->params->coltot[Gi]) {
	    L1 = LIA->matrix[Gi][i][a];
	    F1 = FME->matrix[Gi][i][a];
	  }
	  if(D->params->rowtot[Gjk] && D->params->coltot[Gjk]) {
	    D2 = D->matrix[Gjk][jk][bc];
	    L2 = LIJAB->matrix[Gjk][jk][bc];
	  }
	  W1[Gab][ab][c] += L1 * D2 + F1 * L2;
	}
	/* -L_ib * D_jkac - F_ib * L_jkac */
	if(Gi == Gb && Gjk == Gac) {
	  L1 = D2 = F1 = L2 = 0.0;
	  if(LIA->params->rowtot[Gi] && LIA->params->coltot[Gi]) {
	    L1 = LIA->matrix[Gi][i][b];
	    F1 = FME->matrix[Gi][i][b];
	  }
	  if(D->params->rowtot[Gjk] && D->params->coltot[Gjk]) {
	    D2 = D->matrix[Gjk][jk][ac];
	    L2 = LIJAB->matrix[Gjk][jk][ac];
	  }
	  W1[Gab][ab][c] -= L1 * D2 + F1 * L2;
	}
	/* -L_ic * D_jkba - F_ic * L_jkba */
	if(Gi == Gc && Gjk == Gba) {
	  L1 = D2 = F1 = L2 = 0.0;
	  if(LIA->params->rowtot[Gi] && LIA->params->coltot[Gi]) {
	    L1 = LIA->matrix[Gi][i][c];
	    F1 = FME->matrix[Gi][i][c];
	  }
	  if(D->params->rowtot[Gjk] && D->params->coltot[Gjk]) {
	    D2 = D->matrix[Gjk][jk][ba];
	    L2 = LIJAB->matrix[Gjk][jk][ba];
	  }
	  W1[Gab][ab][c] -= L1 * D2 + F1 * L2;
	}
	/* -L_ja * D_ikbc - F_ja * L_ikbc */
	if(Gj == Ga && Gik == Gbc) {
	  L1 = D2 = F1 = L2 = 0.0;
	  if(LIA->params->rowtot[Gj] && LIA->params->coltot[Gj]) {
	    L1 = LIA->matrix[Gj][j][a];
	    F1 = FME->matrix[Gj][j][a];
	  }
	  if(D->params->rowtot[Gik] && D->params->coltot[Gik]) {
	    D2 = D->matrix[Gik][ik][bc];
	    L2 = LIJAB->matrix[Gik][ik][bc];
	  }
	  W1[Gab][ab][c] -= L1 * D2 + F1 * L2;
	}
	/* +L_jb * D_ikac + F_jb * L_ikac */
	if(Gj == Gb && Gik == Gac) {
	  L1 = D2 = F1 = L2 = 0.0;
	  if(LIA->params->rowtot[Gj] && LIA->params->coltot[Gj]) {
	    L1 = LIA->matrix[Gj][j][b];
	    F1 = FME->matrix[Gj][j][b];
	  }
	  if(D->params->rowtot[Gik] && D->params->coltot[Gik]) {
	    D2 = D->matrix[Gik][ik][ac];
	    L2 = LIJAB->matrix[Gik][ik][ac];
	  }
	  W1[Gab][ab][c] += L1 * D2 + F1 * L2;
	}
	/* +L_jc * D_ikba + F_jc * L_ikba */
	if(Gj == Gc && Gik == Gba) {
	  L1 = D2 = F1 = L2 = 0.0;
	  if(LIA->params->rowtot[Gj] && LIA->params->coltot[Gj]) {
	    L1 = LIA->matrix[Gj][j][c];
	    F1 = FME->matrix[Gj][j][c];
	  }
	  if(D->params->rowtot[Gik] && D->params->coltot[Gik]) {
	    D2 = D->matrix[Gik][ik][ba];
	    L2 = LIJAB->matrix[Gik][ik][ba];
	  }
	  W1[Gab][ab][c] += L1 * D2 + F1 * L2;
	}
	/* -L_ka * D_jibc - F_ka * L_jibc */
	if(Gk == Ga && Gji == Gbc) {
	  L1 = D2 = F1 = L2 = 0.0;
	  if(LIA->params->rowtot[Gk] && LIA->params->coltot[Gk]) {
	    L1 = LIA->matrix[Gk][k][a];
	    F1 = FME->matrix[Gk][k][a];
	  }
	  if(D->params->rowtot[Gji] && D->params->coltot[Gji]) {
	    D2 = D->matrix[Gji][ji][bc];
	    L2 = LIJAB->matrix[Gji][ji][bc];
	  }
	  W1[Gab][ab][c] -= L1 * D2 + F1 * L2;
	}
	/* +L_kb * D_jiac + F_kb * L_jiac */
	if(Gk == Gb && Gji == Gac) {
	  L1 = D2 = F1 = L2 = 0.0;
	  if(LIA->params->rowtot[Gk] && LIA->params->coltot[Gk]) {
	    L1 = LIA->matrix[Gk][k][b];
	    F1 = FME->matrix[Gk][k][b];
	  }
	  if(D->params->rowtot[Gji] && D->params->coltot[Gji]) {
	    D2 = D->matrix[Gji][ji][ac];
	    L2 = LIJAB->matrix[Gji][ji][ac];
	  }
	  W1[Gab][ab][c] += L1 * D2 + F1 * L2;
	}
	/* +L_kc * D_jiba + F_kc * L_jiba */
	if(Gk == Gc && Gji == Gba) {
	  L1 = D2 = F1 = L2 = 0.0;
	  if(LIA->params->rowtot[Gk] && LIA->params->coltot[Gk]) {
	    L1 = LIA->matrix[Gk][k][c];
	    F1 = FME->matrix[Gk][k][c];
	  }
	  if(D->params->rowtot[Gji] && D->params->coltot[Gji]) {
	    D2 = D->matrix[Gji][ji][ba];
	    L2 = LIJAB->matrix[Gji][ji][ba];
	  }
	  W1[Gab][ab][c] += L1 * D2 + F1 * L2;
	}
      }
    }
  }

  for(Gab=0; Gab < nirreps; Gab++) {
    Gc = Gab ^ Gijk; /* assumes totally symmetric! */

    for(ab=0; ab < F->params->coltot[Gab]; ab++) {
      A = F->params->colorb[Gab][ab][0];
      B = F->params->colorb[Gab][ab][1];
      Ga = F->params->rsym[A];
      Gb = F->params->ssym[B];
      a = A - vir_off[Ga];
      b = B - vir_off[Gb];

      for(c=0; c < virtpi[Gc]; c++) {
	C = vir_off[Gc] + c;

	denom = dijk;
	if(fAB->params->rowtot[Ga]) denom -= fAB->matrix[Ga][a][a];
	if(fAB->params->rowtot[Gb]) denom -= fAB->matrix[Gb][b][b];
	if(fAB->params->rowtot[Gc]) denom -= fAB->matrix[Gc][c][c];

	W1[Gab][ab][c] /= denom;

      } /* c */
    } /* ab */
  } /* Gab */

  for(Gab=0; Gab < nirreps; Gab++) {
    Gc = Gab ^ Gijk; /* assumes totally symmetric! */
    global_dpd_->free_dpd_block(W2[Gab], F->params->coltot[Gab], virtpi[Gc]);
  }
  free(W2);

  global_dpd_->file2_mat_close(fIJ);
  global_dpd_->file2_mat_close(fAB);

  global_dpd_->file2_mat_close(FME);
  global_dpd_->file2_mat_close(LIA);

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_close(T2, h);
    global_dpd_->buf4_mat_irrep_close(E, h);
    global_dpd_->buf4_mat_irrep_close(D, h);
    global_dpd_->buf4_mat_irrep_close(LIJAB, h);
  }
}

}} // namespace psi::cclambda
