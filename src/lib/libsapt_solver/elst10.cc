#include "sapt0.h"

namespace psi { namespace sapt {

void SAPT0::elst10()
{
  e_elst10_ = 4.0*C_DDOT(ndf_+3,diagAA_,1,diagBB_,1);
  
  if (print_) {
    fprintf(outfile,"    Elst10,r            = %18.12lf H\n",e_elst10_);
    fflush(outfile);
  }
}

}}
