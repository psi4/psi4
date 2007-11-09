/*! \file symmetrize.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<stdio.h>
#include<cmath>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
namespace psi { namespace CINTS {
void symmetrize_hessian(double **hess)
{
  int i, j;

  for(i=0; i < Molecule.num_atoms*3; i++) {
    for(j=0; j < i; j++) {
      hess[i][j] = hess[j][i] = (hess[i][j] + hess[j][i]);
    }
  }
  
  return;
}
};};
