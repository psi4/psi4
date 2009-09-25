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

void FDD(int i, int C_irr);
void WabefDD(int i, int C_irr);
void WmnijDD(int i, int C_irr);
void WmbejDD(int i, int C_irr);
void WmnefDD(int i, int C_irr);

/* This function computes the H-bar doubles-doubles block contribution
to a Sigma vector stored at Sigma plus 'i' */

void sigmaDD(int i, int C_irr) {

#ifdef TIME_CCEOM
  timer_on("FDD");       FDD(i, C_irr);     timer_off("FDD");
  timer_on("WmnijDD");   WmnijDD(i, C_irr); timer_off("WmnijDD");
  timer_on("WabefDD");   WabefDD(i, C_irr); timer_off("WabefDD");
  timer_on("WmbejDD");   WmbejDD(i, C_irr); timer_off("WmbejDD");
  timer_on("WmnefDD");   WmnefDD(i, C_irr); timer_off("WmnefDD");
#else
  FDD(i, C_irr);
  WmnijDD(i, C_irr);
  WabefDD(i, C_irr);
  WmbejDD(i, C_irr);
  WmnefDD(i, C_irr);
#endif

  return;
}

}} // namespace psi::cceom
