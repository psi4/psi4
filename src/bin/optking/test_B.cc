/*! \file
    \ingroup OPTKING
    \brief TEST_BMAT.CC : compares analytic and numerical B matrices
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

#include <libqt/qt.h>

namespace psi { //namespace optking {

int test_B(const cartesians &carts, simples_class &simples, const salc_set &symm) {
  int i, j, atom, xyz, ncarts, natom, nsalcs;
  double **B, **dq, disp_size=0.01, max_error;
  double *coord, *q_plus, *q_minus;

  natom = carts.get_natom();
  ncarts = 3*natom;
  nsalcs = symm.get_num();

  coord = carts.get_coord();
  simples.compute(coord);
  simples.compute_s(coord);
  B = compute_B(simples,symm);
  free_array(coord);

  if (optinfo.mode == MODE_TEST_BMAT) {
    fprintf(outfile,"\nB Matrix - Analytical, dB_i/(dr angstroms) \n");
    print_mat(B, nsalcs, ncarts, outfile );
  }

  dq = block_matrix(nsalcs,ncarts);

  for (i=0; i<nsalcs; ++i) {
    for (atom=0; atom<natom; ++atom) {
      for (xyz=0; xyz<3; ++xyz) {
        coord = carts.get_coord(); /* coord is in au */

        coord[3*atom+xyz] += disp_size;
        simples.compute(coord);
        //simples.print(outfile, 1);
        q_plus = compute_q(simples,symm); /* q is in Ang and radians */
        //fprintf(outfile,"q +\n");
        //print_mat(&q_plus, 1, nsalcs, outfile);

        coord[3*atom+xyz] -= 2.0*disp_size;
        simples.compute(coord);
        q_minus = compute_q(simples,symm);
 
        dq[i][3*atom+xyz] = (q_plus[i]-q_minus[i]) / (2.0*disp_size*_bohr2angstroms);

        free_array(q_plus);
        free_array(q_minus);
        free_array(coord);
      }
    }
  }
  if (optinfo.mode == MODE_TEST_BMAT) {
    fprintf(outfile,"\nB Matrix - Numerical, disp_size = %lf\n",disp_size);
    print_mat(dq, nsalcs, ncarts, outfile );
  }

  max_error = 0.0;
  for (i=0; i<nsalcs; ++i)
    for (j=0; j<3*natom; ++j)
      if ( fabs(B[i][j]-dq[i][j]) > max_error )
        max_error = fabs(B[i][j]-dq[i][j]);

  fprintf(outfile,"\n\tTesting B-matrix numerically. Maximum difference is %.1e.", max_error);
  if (max_error > 5.0e-3) {
    fprintf(outfile, "\nUh-Oh.  Perhaps a bug or your angular coordinates are at a discontinuity.\n");
    fprintf(outfile, "If the latter, restart your optimization at a new or updated geometry.\n");
    fprintf(outfile, "Remove angular coordinates that are fixed by symmetry\n");
  }
  else {
    fprintf(outfile,"  Looks great.\n");
  }

  free_block(B);
  free_block(dq);
  return 0;
}

}//} /* namespace psi::optking */

