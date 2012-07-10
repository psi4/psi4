/*! \file
    \ingroup QT
    \brief xyz coordinates for three points and R, theta, and phi, returns the
     coordinates a fourth point; angles should enter function in degrees */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libqt/qt.h>
#include <physconst.h>
#include <psifiles.h>
#include <psi4-dec.h>

#define ZMAT_LINEAR_CUTOFF (1.0e-14)

namespace psi {

void zmat_point(double *A, double *B, double *C, double R_CD, double theta_BCD,
  double phi_ABCD, double *D)
{
  double eAB[3],eBC[3],eX[3],eY[3], cosABC, sinABC;
  int xyz;

  theta_BCD *= _pi/180.0;
  phi_ABCD *= _pi/180.0;

  unit_vec(B,A,eAB); /* vector B->A */
  unit_vec(C,B,eBC); /* vector C->B */
  cosABC = -dot_prod(eBC,eAB);

  sinABC = sqrt(1 - (cosABC * cosABC) );
  if ( (sinABC - ZMAT_LINEAR_CUTOFF) < 0.0 ) {
    throw PsiException("Reference points cannot be colinear.",__FILE__,__LINE__);
  }

  cross_prod(eAB,eBC,eY);
  for(xyz=0;xyz<3;xyz++)
    eY[xyz] /= sinABC;
  cross_prod(eY,eBC,eX);
  for (xyz=0;xyz<3;xyz++)
    D[xyz] = C[xyz] + R_CD * ( - eBC[xyz] * cos(theta_BCD) +
                                 eX[xyz] * sin(theta_BCD) * cos(phi_ABCD) +
                                 eY[xyz] * sin(theta_BCD) * sin(phi_ABCD) );
  return;
}

}

