/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

#include"moinfo.h"
#include"moinfo_corr.h"
#include"rmp2_energy.h"


namespace psi { namespace CINTS {
void mp2()
{
  
  init_moinfo();
  init_moinfo_corr();
  switch(UserOptions.reftype) {
  case rhf:
      rmp2_energy();
      break;
  default:
      throw std::domain_error("MP2 energy with specified REFERENCE not implemented");
  }
  cleanup_moinfo_corr();
  cleanup_moinfo();

  return;
}
}}
