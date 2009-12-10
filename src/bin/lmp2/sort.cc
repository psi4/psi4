/*! \file
    \ingroup LMP2
    \brief localized the SCF MO's
*/
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
//#include <cmath>
//#include <libipv1/ip_lib.h>
//#include <libciomr/libciomr.h>
//#include <libchkpt/chkpt.h>
//#include <libpsio/psio.h>
//#include <libqt/qt.h>
//#include <libint/libint.h>
//#include <libmints/basisset.h>
//#include <libmints/onebody.h>
//#include <libmints/twobody.h>
//#include <libmints/integral.h>
//#include <libmints/factory.h>
//#include <libmints/symmetry.h>
//#include <libmints/wavefunction.h>
//#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace lmp2{

void LMP2::sort_shell(int **A, int n) {
  int i, j, k, a;
  int *val;

  val = (int *) malloc(4 * sizeof(int));

  a=0;
  for(i=0; i < n-1; i++) {
      val[a] = A[a][k=i];
      val[a+1] = A[a+1][k=i];
      val[a+2] = A[a+2][k=i];
      val[a+3] = A[a+3][k=i];

      for(j=i+1; j < n; j++)
        if(A[a+1][j] >= val[a+1]) {
          val[a] = A[a][k=j];
          val[a+1] = A[a+1][k=j];
          val[a+2] = A[a+2][k=j];
          val[a+3] = A[a+3][k=j];
        }

      if(k != i) {
        A[a][k] = A[a][i];
        A[a][i] = val[a];
        A[a+1][k] = A[a+1][i];
        A[a+1][i] = val[a+1];
        A[a+2][k] = A[a+2][i];
        A[a+2][i] = val[a+2];
        A[a+3][k] = A[a+3][i];
        A[a+3][i] = val[a+3];
      }
  }
  for(i=0; i < n-1; i++) {
      val[a] = A[a][k=i];
      val[a+1] = A[a+1][k=i];
      val[a+2] = A[a+2][k=i];
      val[a+3] = A[a+3][k=i];

      for(j=i+1; j < n; j++)
        if(A[a+1][j] >= val[a+1] && A[a+3][j] >= val[a+3]) {
          val[a] = A[a][k=j];
          val[a+1] = A[a+1][k=j];
          val[a+2] = A[a+2][k=j];
          val[a+3] = A[a+3][k=j];
        }

      if(k != i) {
        A[a][k] = A[a][i];
        A[a][i] = val[a];
        A[a+1][k] = A[a+1][i];
        A[a+1][i] = val[a+1];
        A[a+2][k] = A[a+2][i];
        A[a+2][i] = val[a+2];
        A[a+3][k] = A[a+3][i];
        A[a+3][i] = val[a+3];
      }
  }

  free(val);
}


}}
