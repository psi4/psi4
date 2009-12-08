/*! \file
    \ingroup LMP2
    \brief localized the SCF MO's
*/
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
//#include <cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
//#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace lmp2{

void opdm() {

  int i, j, k;
  int nso, nocc; 

  nso = mo.nso;
  nocc = mo.nocc;

  ao.D = block_matrix(nso, nso);

  /* Compute the density matrix in tha AO basis */
  C_DGEMM('n', 't', nso, nso, nocc, 1, ao.C[0], nso, ao.C[0], nso, 0, ao.D[0], nso);

  free_block(ao.C);

}

}} // namespace psi::lmp2
