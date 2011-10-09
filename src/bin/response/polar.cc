/*! \file
    \ingroup RESPONSE
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace response {

void transpert(const char *pert);
double dot(double *A, double *B, int n);

void polar(void)
{
  int h, h1, nirreps, row, col, dim;
  int Ga, Gi, i, a, ai, aa, ii;
  double ***R, ***S;
  int alpha, beta;
  char **cartcomp, pert[32];
  double ***C;
  double **polar;

  polar = block_matrix(3,3); /* polarizability tensor */

  cartcomp = (char **) malloc(3 * sizeof(char *));
  cartcomp[0] = strdup("X");
  cartcomp[1] = strdup("Y");
  cartcomp[2] = strdup("Z");

  C = moinfo.RPA_inv;

  R = (double ***) malloc(3 * sizeof(double **));
  S = (double ***) malloc(3 * sizeof(double **));
  for(alpha=0; alpha < 3; alpha++) {
    R[alpha] = (double **) malloc(moinfo.nirreps * sizeof(double *));
    S[alpha] = (double **) malloc(moinfo.nirreps * sizeof(double *));
    for(h=0; h < moinfo.nirreps; h++) {
      R[alpha][h] = init_array(2*moinfo.RPA_dim[h]);
      S[alpha][h] = init_array(2*moinfo.RPA_dim[h]);
    }
  }

  /* prepare the MO-basis dipole moment integrals */
  for(alpha=0; alpha < 3; alpha++) {
    sprintf(pert, "Mu_%1s", cartcomp[alpha]);
    transpert(pert);
  }

  /* Set up the dipole vectors for this irrep */
  for(alpha=0; alpha < 3; alpha++) {
    for(h=0; h < moinfo.nirreps; h++) {
      if(dim = moinfo.RPA_dim[h]) {
        for(Ga=0,ai=0; Ga < moinfo.nirreps; Ga++) {
          Gi = h^Ga;
          for(a=0; a < moinfo.virtpi[Ga]; a++) {
            aa = moinfo.qt2pitzer[moinfo.qt_vir[a] + moinfo.vir_off[Ga]];
            for(i=0; i < moinfo.occpi[Gi]; i++,ai++) {
              ii = moinfo.qt2pitzer[moinfo.qt_occ[i] + moinfo.occ_off[Gi]];
              R[alpha][h][ai] = 2 * moinfo.MU[alpha][aa][ii];
              R[alpha][h][ai+dim] = - 2 * moinfo.MU[alpha][aa][ii];
            }
          }
        }

        C_DGEMV('n', 2*dim, 2*dim, 1, &(C[h][0][0]), 2*dim, &(R[alpha][h][0]), 1, 0, &(S[alpha][h][0]), 1);
      } /* if(dim) */
    }
  }

  for(h=0; h < moinfo.nirreps; h++) {
    if(dim = moinfo.RPA_dim[h]) {
      for(alpha=0; alpha < 3; alpha++) {
        for(beta=0; beta < 3; beta++) {
          polar[alpha][beta] += dot(R[alpha][h],S[beta][h],2*dim);
        }
      }
    }
  }

  for(alpha=0; alpha < 3; alpha++) {
    for(h=0; h < moinfo.nirreps; h++) {
      free(R[alpha][h]); free(S[alpha][h]);
    }
  }
  for(alpha=0; alpha < 3; alpha++) {
    free(R[alpha]); free(S[alpha]);
  }
  free(R);
  free(S);

  fprintf(outfile, "\n\tHartree-Fock Electric Polarizability Tensor  [(e^2 a0^2)/E_h]:\n");
  fprintf(outfile, "\t---------------------------------------------------------------\n");
  mat_print(polar, 3, 3, outfile);
  free_block(polar);

  free(cartcomp[0]);
  free(cartcomp[1]);
  free(cartcomp[2]);
  free(cartcomp);
}

}} // namespace psi::response
