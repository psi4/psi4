#include <libchkpt/chkpt.hpp>
#include "scf.h"

namespace psi{ namespace mcscf{

void SCF::density_matrix()
{
  // Closed-shell density
  O->zero();
  for(int h = 0; h < nirreps; ++h){
    for(int i = 0; i < docc[h]; ++i){
      O->set(h, i, i, 1.0);
    }
  }
  transform(O,Dc,C_T);

  // ROHF density
  if(reference == rohf){
    O->zero();
    for(int h = 0; h < nirreps; ++h){
      for(int i = docc[h]; i < docc[h] + actv[h]; ++i){
        O->set(h, i, i, 1.0);
      }
    }
    transform(O,Do,C_T);
  }

  // TWOCON density
  if(reference == tcscf){
    for(int I = 0 ; I < nci; ++I){
      O->zero();
      O->set(tcscf_sym[I], tcscf_mos[I], tcscf_mos[I], 1.0);
      transform(O,Dtc[I],C_T);
    }
  }
}

}} /* End Namespaces */
