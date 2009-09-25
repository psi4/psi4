/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <stdlib.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <math.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* build_Z_RHF():  Solve the orbital Z-vector equations for RHF refs:
**
**    sum E,M A(AI,EM) D(orb)(E,M) = -X(A,I)
**
** where A(AI,EM) is the orbital Hessian computed in build_A(), X(A,I)
** is the orbital rotation gradient computed in build_X(), and
** D(orb)(E,M) is the final Z-vector we want. 
**
*/

void build_Z_RHF(void)
{
  dpdbuf4 A;
  dpdfile2 X1, D;
  double *X;
  int h, nirreps, a, i, count;

  nirreps = moinfo.nirreps;

  /* Grab only irrep 0 of the orbital Hessian */
  dpd_buf4_init(&A, CC_MISC, 0, 11, 11, 11, 11, 0, "A(EM,AI)");
  dpd_buf4_mat_irrep_init(&A, 0);
  dpd_buf4_mat_irrep_rd(&A, 0);

  /* Place all the elements of the orbital rotation gradient, X into a
     linear array, Z */
  dpd_file2_init(&X1, CC_OEI, 0, 1, 0, "XAI");
  dpd_file2_mat_init(&X1);
  dpd_file2_mat_rd(&X1);
  X = init_array(A.params->rowtot[0]);
  for(h=0,count=0; h < nirreps; h++)
    for(a=0; a < X1.params->rowtot[h]; a++)
      for(i=0; i < X1.params->coltot[h]; i++) 
	X[count++] = -X1.matrix[h][a][i];

  dpd_file2_mat_close(&X1);
  dpd_file2_close(&X1);

  /* Trying out Matt's Pople code --- way to go, Matt! */
  pople(A.matrix[0], X, A.params->rowtot[0], 1, 1e-12, outfile, 0);

  /* Build the orbital component of Dai */
  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  dpd_file2_mat_init(&D);
  for(h=0,count=0; h < nirreps; h++)
    for(a=0; a < D.params->rowtot[h]; a++)
      for(i=0; i < D.params->coltot[h]; i++)
	D.matrix[h][a][i] = X[count++];
  dpd_file2_mat_wrt(&D);
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  free(X);

  dpd_buf4_mat_irrep_close(&A, 0);
  dpd_buf4_close(&A);
}



}} // namespace psi::ccdensity
