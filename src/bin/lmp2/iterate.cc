/*! \file
    \ingroup LMP2
    \brief localized the SCF MO's
*/

//#include <iostream>
//#include <fstream>              // file I/O support
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
//#include <cmath>
//#include <libiwl/iwl.hpp>
//#include <libipv1/ip_lib.h>
//#include <libciomr/libciomr.h>
//#include <libchkpt/chkpt.hpp>
//#include <libpsio/psio.hpp>
//#include <libqt/qt.h>
//#include <libint/libint.h>
//#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{

extern int myid;
extern int nprocs;

namespace lmp2{

extern int myid_lmp2;
extern int nprocs_lmp2;


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

  if(RILMP2){
    fprintf(outfile,"\tRI-LMP2 correlation energy         = %20.15f\n",Emp2);
    fprintf(outfile,"      * RI-LMP2 total energy               = %20.15f\n\n",
      Escf + Emp2);
    fprintf(outfile,"\tOpposite-Spin correlation energy  = %20.15f\n",E_OS);
    fprintf(outfile,"\tSame-Spin correlation energy      = %20.15f\n\n",E_SS);
    fprintf(outfile,"      * SCS-RI-LMP2 total energy           = %20.15f\n\n",
      Escf + scs_scale_os*E_OS + scs_scale_ss*E_SS);
  }

}
}}
