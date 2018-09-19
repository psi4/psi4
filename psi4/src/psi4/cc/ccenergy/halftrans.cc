/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

/* halftrans(): Routine to transform the last two indices of a dpdbuf4
** between the MO and SO bases.
**
** dpdbuf4 *Buf1:    Pointer to the MO dpdbuf4 (already initialized)
** dpdbuf4 *Buf2:    Pointer to the SO dpdbuf4 (already initialized).
** double ***C1:     Pointer to the left transformation matrix (symmetry blocked, SO x MO)
** double ***C2:     Pointer to the right transformation matrix (symmetry blocked, SO x MO)
** int nirreps:      The number of irreps in the point group
** int **mo_row:     A lookup array.  For a dpdbuf4 with MO indices (ij,ab),
**                   given the irrep h of ij (= ab) and the irrep of orbital a, the
**                   array returns the offset of the start of the set of b molecular
**                   orbitals.
** int **so_row:     Like mo_row, but for a dpdbuf4 with the last two
**                   indices in the SO basis.
** int *mospi_left:  The number of MO's per irrep for the left upper index.
** int *mospi_right: The number of MO's per irrep for the right upper index.
** int *sospi:       The number of SO's per irrep.
** int type:         0 = MO --> SO; 1 = SO --> MO
** double alpha:     Multiplicative factor for the transformation
** double beta:      Multiplicative factor for the target
*/

void CCEnergyWavefunction::halftrans(dpdbuf4 *Buf1, int dpdnum1, dpdbuf4 *Buf2, int dpdnum2, double ***C1, double ***C2,
               int nirreps, int **mo_row, int **so_row, int *mospi_left, int *mospi_right,
               int *sospi, int type, double alpha, double beta)
{
  int h, Gc, Gd, cd, pq, ij;
  double **X;

  for(h=0; h < nirreps; h++) {

    dpd_set_default(dpdnum1);
    global_dpd_->buf4_mat_irrep_init(Buf1, h);

    dpd_set_default(dpdnum2);
    global_dpd_->buf4_mat_irrep_init(Buf2, h);

    if(type==0) { /* alpha * Buf1 --> beta * Buf2 */
      if(alpha != 0.0) { dpd_set_default(dpdnum1); global_dpd_->buf4_mat_irrep_rd(Buf1, h); }
      if(beta != 0.0) { dpd_set_default(dpdnum2); global_dpd_->buf4_mat_irrep_rd(Buf2, h); }
    }
    if(type==1) { /* alpha * Buf2 --> beta * Buf1 */
      if(alpha != 0.0) { dpd_set_default(dpdnum2); global_dpd_->buf4_mat_irrep_rd(Buf2, h); }
      if(beta != 0.0) { dpd_set_default(dpdnum1); global_dpd_->buf4_mat_irrep_rd(Buf1, h); }
    }

    for(Gc=0; Gc < nirreps; Gc++) {
      Gd = h^Gc;

      cd = mo_row[h][Gc];
      pq = so_row[h][Gc];

      if(mospi_left[Gc] && mospi_right[Gd] && sospi[Gc] && sospi[Gd]) {

	if(type == 0) {
	  X = block_matrix(mospi_left[Gc],sospi[Gd]);

	  for(ij=0; ij < Buf1->params->rowtot[h]; ij++) {

	    C_DGEMM('n','t', mospi_left[Gc], sospi[Gd], mospi_right[Gd], 1.0,
		    &(Buf1->matrix[h][ij][cd]), mospi_right[Gd], &(C2[Gd][0][0]),
                    mospi_right[Gd], 0.0, &(X[0][0]), sospi[Gd]);

	    C_DGEMM('n','n', sospi[Gc], sospi[Gd], mospi_left[Gc], alpha,
		    &(C1[Gc][0][0]), mospi_left[Gc], &(X[0][0]), sospi[Gd],
		    beta, &(Buf2->matrix[h][ij][pq]), sospi[Gd]);
	  }
	}
	else {
	  X = block_matrix(sospi[Gc],mospi_right[Gd]);

	  for(ij=0; ij < Buf1->params->rowtot[h]; ij++) {

	    C_DGEMM('n','n', sospi[Gc], mospi_right[Gd], sospi[Gd], 1.0,
		    &(Buf2->matrix[h][ij][pq]), sospi[Gd], &(C2[Gd][0][0]), mospi_right[Gd],
		    0.0, &(X[0][0]), mospi_right[Gd]);

	    C_DGEMM('t','n', mospi_left[Gc], mospi_right[Gd], sospi[Gc], alpha,
		    &(C1[Gc][0][0]), mospi_left[Gc], &(X[0][0]), mospi_right[Gd],
		    beta, &(Buf1->matrix[h][ij][cd]), mospi_right[Gd]);

	  }
	}

	free_block(X);
      }
    }

    dpd_set_default(dpdnum1);
    if(type==1) global_dpd_->buf4_mat_irrep_wrt(Buf1, h);
    global_dpd_->buf4_mat_irrep_close(Buf1, h);

    dpd_set_default(dpdnum2);
    if(type==0) global_dpd_->buf4_mat_irrep_wrt(Buf2, h);
    global_dpd_->buf4_mat_irrep_close(Buf2, h);

  }

}
}} // namespace psi::ccenergy
