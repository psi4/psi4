/*! \file
    \ingroup LMP2
    \brief localized the SCF MO's
*/

#include <iostream>
#include <fstream>              // file I/O support
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libiwl/iwl.hpp>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.hpp>
#include <libpsio/psio.hpp>
#include <libqt/qt.h>
#include <libint/libint.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace lmp2{

void LMP2::iterate() {

  iter = 0;
  conv = 0;
  while(conv != 1) {
    get_diis_params();
    amplitudes();
    if(iter > 1 && diis == 1)
      diis_ext();
    energy();
      if(iter > 0)
        check_conv();
    if(myid == 0) 
      print_iteration();
    iter++;
  }

}

}}
