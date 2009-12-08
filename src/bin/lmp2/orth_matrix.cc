/*! \file
    \ingroup LMP2
    \brief Compute the matrix that will transform the residuals into an orthogonal basis
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libint/libint.h>
#include <libmints/basisset.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/factory.h>
#include <libmints/symmetry.h>
#include <libmints/wavefunction.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace lmp2{

void orth_matrix() {

  int a, b;
  double *evals_S;
  double **X;
  double **Xt;
  double **temp;
  double **Fbar;
  double **U;

  evals_S = init_array(mo.nso);
  X = block_matrix(mo.nso,mo.nso);

  sq_rsp(mo.nso, mo.nso, pao.ovlp, evals_S, 1, X, 1.0E-14);

  //   ****  Othogonalize the Overlap Matrix S = UXU (X = overlap^(-1/2))  ****
  Xt = block_matrix(mo.nso,mo.nso);
  for(a=0; a < mo.nso; a++) {
    for(b=0; b < mo.nso; b++) {
      if(fabs(evals_S[a]) > 1e-6) {
        Xt[b][a] = X[b][a]/sqrt(evals_S[a]);
      }
      else Xt[b][a] = 0.0;
    }
  }

  //   ****  Form an Initial Orthogonal Fock Matrix (F = [tran(orths)][Hcore][orths] ****
  temp = block_matrix(mo.nso,mo.nso);
  Fbar = block_matrix(mo.nso,mo.nso);

  C_DGEMM('t', 'n', mo.nso, mo.nso, mo.nso, 1, Xt[0], mo.nso, pao.F[0], mo.nso, 0, temp[0], mo.nso);
  C_DGEMM('n', 'n', mo.nso, mo.nso, mo.nso, 1, temp[0], mo.nso, Xt[0], mo.nso, 0, Fbar[0], mo.nso);

  U = block_matrix(mo.nso,mo.nso);
  sq_rsp(mo.nso, mo.nso, Fbar, ao.evals, 1, U, 1.0E-14);

  lo.W = block_matrix(mo.nso,mo.nso);
  C_DGEMM('n', 'n', mo.nso, mo.nso, mo.nso, 1, Xt[0], mo.nso, U[0], mo.nso, 0, lo.W[0], mo.nso);

}

}} // namespace psi::lmp2
