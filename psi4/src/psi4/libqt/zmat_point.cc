/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup QT
    \brief xyz coordinates for three points and R, theta, and phi, returns the
     coordinates a fourth point; angles should enter function in degrees */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/psi4-dec.h"
#include "psi4/physconst.h"

#define ZMAT_LINEAR_CUTOFF (1.0e-14)

namespace psi {

void zmat_point(double *A, double *B, double *C, double R_CD, double theta_BCD,
  double phi_ABCD, double *D)
{
  double eAB[3],eBC[3],eX[3],eY[3], cosABC, sinABC;
  int xyz;

  theta_BCD *= pc_pi/180.0;
  phi_ABCD *= pc_pi/180.0;

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
