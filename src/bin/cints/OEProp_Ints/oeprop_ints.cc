/*! \file oeprop_ints.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

#include"moment_ints.h"
#include"moment_deriv1.h"
#include"angmom_ints.h"

namespace psi { namespace CINTS {
void oeprop_ints()
{
  moment_ints();
  moment_deriv1();
  angmom_ints();
  return;
}
};};
