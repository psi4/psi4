#include <cmath>

#include <libchkpt/chkpt.hpp>
#include "scf.h"

namespace psi{ namespace mcscf{

void SCF::construct_S_inverse_sqrt()
{
  SBlockVector lambda("lambda",nirreps,sopi);
  SBlockMatrix L("L",nirreps,sopi,sopi);
  SBlockMatrix Lambda("Lambda",nirreps,sopi,sopi);

  S.diagonalize(L,lambda);

//   lambda->print();
//   L->print();

  for(int h = 0; h < nirreps; ++h){
    for(int i = 0; i < sopi[h]; ++i){
      Lambda->set(h, i, i, 1.0 / sqrt(lambda->get(h,i)) );
    }
  }


  T.multiply(false,true,Lambda,L);
  S_sqrt_inv.multiply(false,false,L,T);

  for(int h = 0; h < nirreps; ++h){
    for(int i = 0; i < sopi[h]; ++i){
      Lambda->set(h, i, i, sqrt(lambda->get(h,i)) );
    }
  } 

  T.multiply(false,true,Lambda,L);
  S_sqrt.multiply(false,false,L,T);
}

}} /* End Namespaces */
