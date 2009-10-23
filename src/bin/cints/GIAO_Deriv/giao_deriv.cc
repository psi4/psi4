/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/

#include "giao_oe_deriv.h"
#include "giao_te_deriv.h"

#define DO_OEI 1
#define DO_TEI 1

namespace psi { namespace cints {
void giao_deriv()
{
#if DO_OEI
  giao_oe_deriv();
#endif
#if DO_TEI
  giao_te_deriv();
#endif
}

};};
