/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>

namespace psi { namespace cclambda {

/* halftrans(): Routine to transform the last two indices of a dpdbuf4
** between the MO and SO bases.
** 
** dpdbuf4 *Buf1: Pointer to the MO dpdbuf4 (already initialized)
** dpdbuf4 *Buf2: Pointer to the SO dpdbuf4 (already initialized).
** double ***C:   Pointer to the transformation matrix (symmetry blocked, SO x MO)
** int nirreps:   The number of irreps in the point group
** int **mo_row:  A lookup array.  For a dpdbuf4 with MO indices (ij,ab),
**                given the irrep h of ij (= ab) and the irrep of orbital a, the
**                array returns the offset of the start of the set of b molecular
**                orbitals.
** int **so_row:  Like mo_row, but for a dpdbuf4 with the last two
**                indices in the SO basis.
** int *mospi:    The number of MO's per irrep.
** int *sospi:    The number of SO's per irrep.
** int type:      0 = MO --> SO; 1 = SO --> MO
** double alpha:  multiplicative factor for the transformation
** double beta:   multiplicative factor for the target
*/

void halftrans(dpdbuf4 *Buf1, int dpdnum1, dpdbuf4 *Buf2, int dpdnum2, double ***C, int nirreps, 
	       int **mo_row, int **so_row, int *mospi, int *sospi, int type, double alpha, double beta)
{
  int h, Gc, Gd, cd, pq, ij;
  double **X;

  for(h=0; h < nirreps; h++) {

    dpd_set_default(dpdnum1);
    dpd_buf4_mat_irrep_init(Buf1, h);

    dpd_set_default(dpdnum2);
    dpd_buf4_mat_irrep_init(Buf2, h);

    if(type==0) { /* alpha * Buf1 --> beta * Buf2 */
      if(alpha != 0.0) { dpd_set_default(dpdnum1); dpd_buf4_mat_irrep_rd(Buf1, h); }
      if(beta != 0.0) { dpd_set_default(dpdnum2); dpd_buf4_mat_irrep_rd(Buf2, h); }
    }
    if(type==1) { /* alpha * Buf2 --> beta * Buf1 */
      if(alpha != 0.0) { dpd_set_default(dpdnum2); dpd_buf4_mat_irrep_rd(Buf2, h); }
      if(beta != 0.0) { dpd_set_default(dpdnum1); dpd_buf4_mat_irrep_rd(Buf1, h); }
    }

    for(Gc=0; Gc < nirreps; Gc++) {
      Gd = h^Gc;

      cd = mo_row[h][Gc];
      pq = so_row[h][Gc];

      if(mospi[Gc] && mospi[Gd] && sospi[Gc] && sospi[Gd]) {

	if(type == 0) {
	  X = block_matrix(mospi[Gc],sospi[Gd]);

	  for(ij=0; ij < Buf1->params->rowtot[h]; ij++) {

	    C_DGEMM('n','t', mospi[Gc], sospi[Gd], mospi[Gd], 1.0,
		    &(Buf1->matrix[h][ij][cd]), mospi[Gd], &(C[Gd][0][0]), mospi[Gd],
		    0.0, &(X[0][0]), sospi[Gd]);

	    C_DGEMM('n','n', sospi[Gc], sospi[Gd], mospi[Gc], alpha, 
		    &(C[Gc][0][0]), mospi[Gc], &(X[0][0]), sospi[Gd],
		    beta, &(Buf2->matrix[h][ij][pq]), sospi[Gd]);
	  }
	}
	else {
	  X = block_matrix(sospi[Gc],mospi[Gd]);

	  for(ij=0; ij < Buf1->params->rowtot[h]; ij++) {

	    C_DGEMM('n','n', sospi[Gc], mospi[Gd], sospi[Gd], 1.0,
		    &(Buf2->matrix[h][ij][pq]), sospi[Gd], &(C[Gd][0][0]), mospi[Gd],
		    0.0, &(X[0][0]), mospi[Gd]);

	    C_DGEMM('t','n', mospi[Gc], mospi[Gd], sospi[Gc], alpha, 
		    &(C[Gc][0][0]), mospi[Gc], &(X[0][0]), mospi[Gd],
		    beta, &(Buf1->matrix[h][ij][cd]), mospi[Gd]);

	  }
	}

	free_block(X);
      }
    }

    dpd_set_default(dpdnum1);
    if(type==1) dpd_buf4_mat_irrep_wrt(Buf1, h);
    dpd_buf4_mat_irrep_close(Buf1, h);

    dpd_set_default(dpdnum2);
    if(type==0) dpd_buf4_mat_irrep_wrt(Buf2, h);
    dpd_buf4_mat_irrep_close(Buf2, h);

  }

}

}} // namespace psi::cclambda
