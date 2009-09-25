/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void FSD(int i, int C_irr);
void WamefSD(int i, int C_irr);
void WmnieSD(int i, int C_irr);

/* This function computes the H-bar singles-doubles block contribution
to a Sigma vector stored at Sigma plus 'i' */

void sigmaSD(int i, int C_irr) {

#ifdef TIME_CCEOM
  timer_on("FSD");     FSD(i, C_irr);     timer_off("FSD");
  timer_on("WamefSD"); WamefSD(i, C_irr); timer_off("WamefSD");
  timer_on("WmnieSD"); WmnieSD(i, C_irr); timer_off("WmnieSD");
#else
  FSD(i, C_irr);
  WamefSD(i, C_irr);
  WmnieSD(i, C_irr);
#endif

  return;
}

}} // namespace psi::cceom
