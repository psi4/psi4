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

#include"ccinfo.h"


namespace psi { namespace cints {
  extern void cc_bt2();
  
  void direct_cc()
  {
    
    init_ccinfo();
    
    if (UserOptions.make_cc_bt2)
      cc_bt2();
    
    cleanup_ccinfo();
    
    return;
  }
};};
