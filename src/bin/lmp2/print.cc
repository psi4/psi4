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

void print() {


  if(params.myid == 0) {
    if(lo.iter > 0) {
      fprintf(outfile, "%d\t %20.12f\t %20.12f %20.12f\n", lo.iter, lo.Emp2, lo.DEmp2, lo.Drms);
    }
    else {
      fprintf(outfile, "\nBegin LMP2 Iterations\n");
      fprintf(outfile, "Iter\t           LMP2\t\t\t Delta(LMP2)\t\t RMS(T)\n");
      fprintf(outfile, "%d\t %20.12f\n", lo.iter, lo.Emp2);
    }
  }

}

}} // namespace psi::lmp2
