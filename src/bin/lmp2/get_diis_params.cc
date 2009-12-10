/*! \file
    \ingroup LMP2
    \brief Enter brief description of file here
*/

//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
//#include <libipv1/ip_lib.h>
//#include <libciomr/libciomr.h>
//#include <libchkpt/chkpt.h>
//#include <libpsio/psio.h>
//#include <libqt/qt.h>
//#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace lmp2{

void LMP2::get_diis_params() {

  if(diis == 1) {
    //  **** Set up variables for DIIS extrapolation ****
    if(iter < ndiis-1) div = iter;
    else div = iter%ndiis;

    if(iter <= ndiis) matsize = iter - 1;
    else matsize = ndiis;

    if(div == 0) dmat1 = ndiis - 1;
    else dmat1 = div - 1;

    if(div == 0) dmat2 = ndiis - 2;
    else if(div == 1) dmat2 = ndiis - 1;
    else dmat2 = div - 2;

    nmat = iter%2;
    if(nmat == 1) omat = 0;
    else omat = 1;
  }
  else {
    div = iter%2;
    if(div == 1) dmat1 = 0;
    else dmat1 = 1;
  }

}

}} /* End namespace */

