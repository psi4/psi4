/*! \file
    \ingroup BASIS
    \brief Enter brief description of file here 
*/

#include <stdexcept>
#include <libciomr/libciomr.h>
#include "osrecur.h"

namespace psi {

OI_OSRecursor::OI_OSRecursor(int maxam1, int maxam2) :
  maxam1_(maxam1), maxam2_(maxam2)
{
  if (maxam1 < 0)
    throw std::runtime_error("ERROR: OI_OSRecursor::OI_OSRecursor -- maxam1 must be nonnegative");
  if (maxam2 < 0)
    throw std::runtime_error("ERROR: OI_OSRecursor::OI_OSRecursor -- maxam2 must be nonnegative");

  OIX_ = block_matrix(maxam1_+1,maxam2_+1);
  OIY_ = block_matrix(maxam1_+1,maxam2_+1);
  OIZ_ = block_matrix(maxam1_+1,maxam2_+1);
}

OI_OSRecursor::~OI_OSRecursor()
{
  free_block(OIX_);
  free_block(OIY_);
  free_block(OIZ_);
}

void OI_OSRecursor::compute(PSI_FLOAT PA[3], PSI_FLOAT PB[3], PSI_FLOAT gamma, int am1, int am2)
{
  if (am1 < 0 || am1 > maxam1_)
    throw std::runtime_error("ERROR: OI_OSRecursor::compute -- am1 out of bounds");
  if (am2 < 0 || am2 > maxam2_)
    throw std::runtime_error("ERROR: OI_OSRecursor::compute -- am2 out of bounds");

  int i,j;
  PSI_FLOAT pp = 1/(2*gamma);
  int lmaxi = am1;
  int lmaxj = am2;

  OIX_[0][0] = OIY_[0][0] = OIZ_[0][0] = 1.0;

	/* Upward recursion in j for i=0 */

  OIX_[0][1] = PB[0];
  OIY_[0][1] = PB[1];
  OIZ_[0][1] = PB[2];

  for(j=1;j<lmaxj;j++) {
    OIX_[0][j+1] = PB[0]*OIX_[0][j];
    OIY_[0][j+1] = PB[1]*OIY_[0][j];
    OIZ_[0][j+1] = PB[2]*OIZ_[0][j];
    OIX_[0][j+1] += j*pp*OIX_[0][j-1];
    OIY_[0][j+1] += j*pp*OIY_[0][j-1];
    OIZ_[0][j+1] += j*pp*OIZ_[0][j-1];
  }

  /* Upward recursion in i for all j's */
  if (lmaxi > 0) {
    OIX_[1][0] = PA[0];
    OIY_[1][0] = PA[1];
    OIZ_[1][0] = PA[2];
    for(j=1;j<=lmaxj;j++) {
      OIX_[1][j] = PA[0]*OIX_[0][j];
      OIY_[1][j] = PA[1]*OIY_[0][j];
      OIZ_[1][j] = PA[2]*OIZ_[0][j];
      OIX_[1][j] += j*pp*OIX_[0][j-1];
      OIY_[1][j] += j*pp*OIY_[0][j-1];
      OIZ_[1][j] += j*pp*OIZ_[0][j-1];
    }
    for(i=1;i<lmaxi;i++) {
      OIX_[i+1][0] = PA[0]*OIX_[i][0];
      OIY_[i+1][0] = PA[1]*OIY_[i][0];
      OIZ_[i+1][0] = PA[2]*OIZ_[i][0];
      OIX_[i+1][0] += i*pp*OIX_[i-1][0];
      OIY_[i+1][0] += i*pp*OIY_[i-1][0];
      OIZ_[i+1][0] += i*pp*OIZ_[i-1][0];
      for(j=1;j<=lmaxj;j++) {
        OIX_[i+1][j] = PA[0]*OIX_[i][j];
        OIY_[i+1][j] = PA[1]*OIY_[i][j];
        OIZ_[i+1][j] = PA[2]*OIZ_[i][j];
        OIX_[i+1][j] += i*pp*OIX_[i-1][j];
        OIY_[i+1][j] += i*pp*OIY_[i-1][j];
        OIZ_[i+1][j] += i*pp*OIZ_[i-1][j];
        OIX_[i+1][j] += j*pp*OIX_[i][j-1];
        OIY_[i+1][j] += j*pp*OIY_[i][j-1];
        OIZ_[i+1][j] += j*pp*OIZ_[i][j-1];
      }
    }
  }
}

}

