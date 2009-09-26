/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libchkpt/chkpt.h>

#include "defines.h"
#define EXTERN
#include "global.h"

#include "float.h"

namespace psi { namespace input {

void OI_OSrecurs_float(FLOAT** OIX, FLOAT** OIY, FLOAT** OIZ, FLOAT PA[3], FLOAT PB[3],
		 FLOAT gamma, int lmaxi, int lmaxj)
{
  int i,j,k;
  FLOAT pp = 1/(2*gamma);

  OIX[0][0] = OIY[0][0] = OIZ[0][0] = 1.0;

  if (!lmaxi && !lmaxj)
    return;

	/* Upward recursion in j for i=0 */
  if (lmaxj) {
    OIX[0][1] = PB[0];
    OIY[0][1] = PB[1];
    OIZ[0][1] = PB[2];
  }

  for(j=1;j<lmaxj;j++) {
    OIX[0][j+1] = PB[0]*OIX[0][j];
    OIY[0][j+1] = PB[1]*OIY[0][j];
    OIZ[0][j+1] = PB[2]*OIZ[0][j];
    OIX[0][j+1] += j*pp*OIX[0][j-1];
    OIY[0][j+1] += j*pp*OIY[0][j-1];
    OIZ[0][j+1] += j*pp*OIZ[0][j-1];
  }

	/* Upward recursion in i for all j's */

  if (lmaxi) {
    OIX[1][0] = PA[0];
    OIY[1][0] = PA[1];
    OIZ[1][0] = PA[2];
  }
  for(j=1;j<=lmaxj;j++) {
    OIX[1][j] = PA[0]*OIX[0][j];
    OIY[1][j] = PA[1]*OIY[0][j];
    OIZ[1][j] = PA[2]*OIZ[0][j];
    OIX[1][j] += j*pp*OIX[0][j-1];
    OIY[1][j] += j*pp*OIY[0][j-1];
    OIZ[1][j] += j*pp*OIZ[0][j-1];
  }
  for(i=1;i<lmaxi;i++) {
    OIX[i+1][0] = PA[0]*OIX[i][0];
    OIY[i+1][0] = PA[1]*OIY[i][0];
    OIZ[i+1][0] = PA[2]*OIZ[i][0];
    OIX[i+1][0] += i*pp*OIX[i-1][0];
    OIY[i+1][0] += i*pp*OIY[i-1][0];
    OIZ[i+1][0] += i*pp*OIZ[i-1][0];
    for(j=1;j<=lmaxj;j++) {
      OIX[i+1][j] = PA[0]*OIX[i][j];
      OIY[i+1][j] = PA[1]*OIY[i][j];
      OIZ[i+1][j] = PA[2]*OIZ[i][j];
      OIX[i+1][j] += i*pp*OIX[i-1][j];
      OIY[i+1][j] += i*pp*OIY[i-1][j];
      OIZ[i+1][j] += i*pp*OIZ[i-1][j];
      OIX[i+1][j] += j*pp*OIX[i][j-1];
      OIY[i+1][j] += j*pp*OIY[i][j-1];
      OIZ[i+1][j] += j*pp*OIZ[i][j-1];
    }
  }

  return;
}

}} // namespace psi::input
