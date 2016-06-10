/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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
    \ingroup TRANSQT2
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <exception.h>
#define EXTERN
#include "globals.h"

#include <libmints/wavefunction.h>
#include <libtrans/mospace.h>
#include <libmints/matrix.h>
#include <libmints/vector.h>

// This code is no longer used.  -TDC 12/2015

/*
** semicanonical_fock(): Compute the alpha- and beta-spin Fock
** matrices using the ROHF orbitals from the wfn and
** diagonalize the occ-occ and vir-vir blocks to form semicanonical
** orbitals, which are written in the checkpoint file for the
** subsequent transformation.  These orbitals are required for
** perturbation-based methods (e.g., (T), CC2, CC3, etc.) built upon
** an ROHF reference.
**
** This code is based on MLA's code in the old transqt, which was
** based on TDC's Brueckner code in ccenergy.
**
** TDC, 7/06
*/

namespace psi {

  namespace transqt2 {

void uhf_fock_build(double **fock_a, double **fock_b, double **D_a, double **D_b);

void semicanonical_fock(void)
{
  int i, j;
  double *alpha_evals, *beta_evals;
  double **C_a, **C_b;
  int nirreps, nmo, nso, ij, p, q;

  nirreps = moinfo.nirreps;
  nso = moinfo.nso;
  nmo = moinfo.nmo;

  /* Write Semicanonical Alpha and Beta Fock Matrix Eigenvectors
     and Eigenvalues to the Wfn */

  throw PSIEXCEPTION("SemiCanonical transform does not work at the moment");
  //Process::environment.legacy_wavefunction()->semicanonicalize();

  C_a = block_matrix(nmo,nmo);
  C_b = block_matrix(nmo,nmo);
  C_a = Process::environment.legacy_wavefunction()->Ca()->to_block_matrix();
  C_b = Process::environment.legacy_wavefunction()->Cb()->to_block_matrix();
  alpha_evals = init_array(nmo);
  beta_evals = init_array(nmo);
  alpha_evals = Process::environment.legacy_wavefunction()->epsilon_a()->to_block_vector();
  beta_evals = Process::environment.legacy_wavefunction()->epsilon_b()->to_block_vector();

  if(params.print_lvl > 2) {
    outfile->Printf( "\nAlpha Eigenvalues\n");
    for (i=0; i<nmo; i++)
      outfile->Printf( "%10.7lf\n", alpha_evals[i]);
    outfile->Printf( "\nBeta Eigenvalues\n");
    for (i=0; i<nmo; i++)
      outfile->Printf( "%10.7lf\n", beta_evals[i]);
    

    outfile->Printf( "\nAlpha Eigenvectors\n");
    print_mat(C_a, nso, nmo, "outfile");
    outfile->Printf( "\nBeta Eigenvectors\n");
    print_mat(C_b, nso, nmo, "outfile");
  }

  free_block(C_a);
  free_block(C_b);
  free(alpha_evals);
  free(beta_evals);
}


void uhf_fock_build(double **fock_a, double **fock_b, double **D_a, double **D_b)
{
  throw PSIEXCEPTION("UHF Fock build is deprecated");
  int i, j, ij;
  int nso=0;
  int ntri=0;
  int stat=0;
  double *scratch;
  int lastbuf, idx;
  int p=0; int q=0; int r=0; int s=0;
  int pq=0; int rs=0;
  double value;
  Value *valptr;
  Label *lblptr;
  struct iwlbuf InBuf;
  double **Dt;

  nso = moinfo.nso;
  ntri = nso*(nso+1)/2;

  Dt = block_matrix(nso, nso);
  for(p=0; p < nso; p++)
    for(q=0; q < nso; q++)
      Dt[p][q] = D_a[p][q] + D_b[p][q];

  /* one-electron contributions */
  double **H = Process::environment.legacy_wavefunction()->H()->to_block_matrix();

  for(i=0; i < nso; i++) {
      for(j=0; j <= i; j++) {
          fock_a[i][j] = fock_a[j][i] = H[i][j];
          fock_b[i][j] = fock_b[j][i] = H[i][j];
      }
  }

  iwl_buf_init(&InBuf, PSIF_SO_TEI, 0.0, 1, 1);
  do {

    lastbuf = InBuf.lastbuf;
    lblptr = InBuf.labels;
    valptr = InBuf.values;

    for(idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
      p = abs((int) lblptr[idx++]);
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];
      value = (double) valptr[InBuf.idx];

      pq = INDEX(p,q);
      rs = INDEX(r,s);

      /* outfile->Printf( "%d %d %d %d [%d] [%d] %20.15f\n", p, q, r, s, pq, rs, value); */

      /* (pq|rs) */
      fock_a[p][q] += Dt[r][s] * value;
      fock_a[p][r] -= D_a[q][s] * value;
      fock_b[p][q] += Dt[r][s] * value;
      fock_b[p][r] -= D_b[q][s] * value;

      if(p!=q && r!=s && pq != rs) {

        /* (pq|sr) */
        fock_a[p][q] += Dt[s][r] * value;
        fock_a[p][s] -= D_a[q][r] * value;
        fock_b[p][q] += Dt[s][r] * value;
        fock_b[p][s] -= D_b[q][r] * value;

        /* (qp|rs) */
        fock_a[q][p] += Dt[r][s] * value;
        fock_a[q][r] -= D_a[p][s] * value;
        fock_b[q][p] += Dt[r][s] * value;
        fock_b[q][r] -= D_b[p][s] * value;

        /* (qp|sr) */
        fock_a[q][p] += Dt[s][r] * value;
        fock_a[q][s] -= D_a[p][r] * value;
        fock_b[q][p] += Dt[s][r] * value;
        fock_b[q][s] -= D_b[p][r] * value;

        /* (rs|pq) */
        fock_a[r][s] += Dt[p][q] * value;
        fock_a[r][p] -= D_a[s][q] * value;
        fock_b[r][s] += Dt[p][q] * value;
        fock_b[r][p] -= D_b[s][q] * value;

        /* (rs|qp) */
        fock_a[r][s] += Dt[q][p] * value;
        fock_a[r][q] -= D_a[s][p] * value;
        fock_b[r][s] += Dt[q][p] * value;
        fock_b[r][q] -= D_b[s][p] * value;

        /* (sr|pq) */
        fock_a[s][r] += Dt[p][q] * value;
        fock_a[s][p] -= D_a[r][q] * value;
        fock_b[s][r] += Dt[p][q] * value;
        fock_b[s][p] -= D_b[r][q] * value;

        /* (sr|qp) */
        fock_a[s][r] += Dt[q][p] * value;
        fock_a[s][q] -= D_a[r][p] * value;
        fock_b[s][r] += Dt[q][p] * value;
        fock_b[s][q] -= D_b[r][p] * value;
      }
      else if(p!=q && r!=s && pq==rs) {

        /* (pq|sr) */
        fock_a[p][q] += Dt[s][r] * value;
        fock_a[p][s] -= D_a[q][r] * value;
        fock_b[p][q] += Dt[s][r] * value;
        fock_b[p][s] -= D_b[q][r] * value;

        /* (qp|rs) */
        fock_a[q][p] += Dt[r][s] * value;
        fock_a[q][r] -= D_a[p][s] * value;
        fock_b[q][p] += Dt[r][s] * value;
        fock_b[q][r] -= D_b[p][s] * value;

        /* (qp|sr) */
        fock_a[q][p] += Dt[s][r] * value;
        fock_a[q][s] -= D_a[p][r] * value;
        fock_b[q][p] += Dt[s][r] * value;
        fock_b[q][s] -= D_b[p][r] * value;

      }
      else if(p!=q && r==s) {

        /* (qp|rs) */
        fock_a[q][p] += Dt[r][s] * value;
        fock_a[q][r] -= D_a[p][s] * value;
        fock_b[q][p] += Dt[r][s] * value;
        fock_b[q][r] -= D_b[p][s] * value;

        /* (rs|pq) */
        fock_a[r][s] += Dt[p][q] * value;
        fock_a[r][p] -= D_a[s][q] * value;
        fock_b[r][s] += Dt[p][q] * value;
        fock_b[r][p] -= D_b[s][q] * value;

        /* (rs|qp) */
        fock_a[r][s] += Dt[q][p] * value;
        fock_a[r][q] -= D_a[s][p] * value;
        fock_b[r][s] += Dt[q][p] * value;
        fock_b[r][q] -= D_b[s][p] * value;

      }
      else if(p==q && r!=s) {

        /* (pq|sr) */
        fock_a[p][q] += Dt[s][r] * value;
        fock_a[p][s] -= D_a[q][r] * value;
        fock_b[p][q] += Dt[s][r] * value;
        fock_b[p][s] -= D_b[q][r] * value;

        /* (rs|pq) */
        fock_a[r][s] += Dt[p][q] * value;
        fock_a[r][p] -= D_a[s][q] * value;
        fock_b[r][s] += Dt[p][q] * value;
        fock_b[r][p] -= D_b[s][q] * value;

        /* (sr|pq) */
        fock_a[s][r] += Dt[p][q] * value;
        fock_a[s][p] -= D_a[r][q] * value;
        fock_b[s][r] += Dt[p][q] * value;
        fock_b[s][p] -= D_b[r][q] * value;

      }
      else if(p==q && r==s && pq!=rs) {

        /* (rs|pq) */
        fock_a[r][s] += Dt[p][q] * value;
        fock_a[r][p] -= D_a[s][q] * value;
        fock_b[r][s] += Dt[p][q] * value;
        fock_b[r][p] -= D_b[s][q] * value;

      }
    }

    if(!lastbuf) iwl_buf_fetch(&InBuf);

  } while (!lastbuf);
  iwl_buf_close(&InBuf, 1);

  free_block(Dt);
}

  } // namespace transqt2
} // namespace psi