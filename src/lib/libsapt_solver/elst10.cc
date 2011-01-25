/* This function calculates the Elst10 energy */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#include "sapt0.h"

namespace psi { namespace sapt {

void SAPT0::elst10()
{
  if (params_.print)
    fprintf(outfile,"Begining Elst10 Calculation\n\n");

  results_.elst10 = 4.0*C_DDOT(calc_info_.nrio,calc_info_.diagAA,1,
    calc_info_.diagBB,1);

  if (params_.print) {
    fprintf(outfile,"E^{(10)}_{elst}     = %18.12lf mH\n\n",
      results_.elst10*1000.0);
    fflush(outfile);
  }
}

}}
