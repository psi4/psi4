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
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

double LCX(const char *pert_c, int irrep_c,
	   const char *pert_x, int irrep_x, double omega);
double HXY(const char *pert_x, int irrep_x, double omega_x,
	   const char *pert_y, int irrep_y, double omega_y);
double LHX1Y1(const char *pert_x, int irrep_x, double omega_x,
	      const char *pert_y, int irrep_y, double omega_y);
double LHX2Y2(const char *pert_x, int irrep_x, double omega_x,
	      const char *pert_y, int irrep_y, double omega_y);
double LHX1Y2(const char *pert_x, int irrep_x, double omega_x,
	      const char *pert_y, int irrep_y, double omega_y);
double cc2_LHX1Y1(const char *pert_x, int irrep_x, double omega_x,
		  const char *pert_y, int irrep_y, double omega_y);
double cc2_LHX1Y2(const char *pert_x, int irrep_x, double omega_x,
		  const char *pert_y, int irrep_y, double omega_y);

void linresp(double *tensor, double A, double B,
	     const char *pert_x, int x_irrep, double omega_x,
	     const char *pert_y, int y_irrep, double omega_y)
{
  double polar, polar_LCX, polar_HXY, polar_LHX1Y1, polar_LHX1Y2, polar_LHX2Y2;

  /* clear out scratch space */
  for(int j=PSIF_CC_TMP; j <= PSIF_CC_TMP11; j++) {
     psio_close(j,0); psio_open(j,0);
  }

  polar_LCX = 0.0;
  polar_HXY = 0.0;
  polar_LHX1Y1 = 0.0;
  polar_LHX2Y2 = 0.0;
  polar_LHX1Y2 = 0.0;

  if((x_irrep^y_irrep)==0) {

    if(omega_y != 0.0) {  /* we assume omega_x = -omega_y */
      timer_on("linear terms");
      polar_LCX = LCX(pert_x, x_irrep, pert_y, y_irrep, omega_y);
      polar_LCX += LCX(pert_y, y_irrep, pert_x, x_irrep, omega_x);
      timer_off("linear terms");

      if(!params.sekino && !params.linear) {
        if(params.wfn == "CC2") {
          timer_on("quad terms");
          polar_HXY = HXY(pert_x, x_irrep, omega_x, pert_y, y_irrep, omega_y);
          polar_LHX1Y1 = cc2_LHX1Y1(pert_x, x_irrep, omega_x, pert_y, y_irrep, omega_y);
          polar_LHX1Y2 = cc2_LHX1Y2(pert_x, x_irrep, omega_x, pert_y, y_irrep, omega_y);
          polar_LHX1Y2 += cc2_LHX1Y2(pert_y, y_irrep, omega_y, pert_x, x_irrep, omega_x);
          timer_off("quad terms");
        }
        else {
          timer_on("quad terms");
          polar_LHX1Y1 = LHX1Y1(pert_x, x_irrep, omega_x, pert_y, y_irrep, omega_y);
          polar_LHX2Y2 = LHX2Y2(pert_x, x_irrep, omega_x, pert_y, y_irrep, omega_y);
          polar_LHX1Y2 = LHX1Y2(pert_x, x_irrep, omega_x, pert_y, y_irrep, omega_y);
          polar_LHX1Y2 += LHX1Y2(pert_y, y_irrep, omega_y, pert_x, x_irrep, omega_x);
          timer_off("quad terms");
        }
      }
    }
    else {
      timer_on("linear terms");
      polar_LCX = LCX(pert_x, x_irrep, pert_y, y_irrep, 0.0);
      polar_LCX += LCX(pert_y, y_irrep, pert_x, x_irrep, 0.0);
      timer_off("linear terms");
      if(!params.sekino && !params.linear) {
        if(params.wfn == "CC2") {
          timer_on("quad terms");
          polar_HXY = HXY(pert_x, x_irrep, 0.0, pert_y, y_irrep, 0.0);
          polar_LHX1Y1 = cc2_LHX1Y1(pert_x, x_irrep, 0.0, pert_y, y_irrep, 0.0);
          polar_LHX1Y2 = cc2_LHX1Y2(pert_x, x_irrep, 0.0, pert_y, y_irrep, 0.0);
          polar_LHX1Y2 += cc2_LHX1Y2(pert_y, y_irrep, 0.0, pert_x, x_irrep,0.0);
          timer_off("quad terms");
        }
        else {
          timer_on("quad terms");
          polar_LHX1Y1 = LHX1Y1(pert_x, x_irrep, 0.0, pert_y, y_irrep, 0.0);
          polar_LHX2Y2 = LHX2Y2(pert_x, x_irrep, 0.0, pert_y, y_irrep, 0.0);
          polar_LHX1Y2 = LHX1Y2(pert_x, x_irrep, 0.0, pert_y, y_irrep, 0.0);
          polar_LHX1Y2 += LHX1Y2(pert_y, y_irrep, 0.0, pert_x, x_irrep, 0.0);
          timer_off("quad terms");
        }
      }
    }

    if(params.sekino || params.linear) /* only linear term needed in Sekino-Bartlett model III */
      polar = polar_LCX;
    else
      polar = polar_LCX+polar_HXY+polar_LHX1Y1+polar_LHX2Y2+polar_LHX1Y2;

    if(params.print & 2) {
      outfile->Printf( "\n\tLinresp tensor <<%s;%s>>\n", pert_x, pert_y);
      outfile->Printf( "\tpolar_LCX    = %20.12f\n", polar_LCX);
      if(params.wfn == "CC2")
        outfile->Printf( "\tpolar_HXY    = %20.12f\n", polar_HXY);
      outfile->Printf( "\tpolar_LHX1Y1 = %20.12f\n", polar_LHX1Y1);
      outfile->Printf( "\tpolar_LHX1Y2 = %20.12f\n", polar_LHX1Y2);
      outfile->Printf( "\tpolar_LHX2Y2 = %20.12f\n", polar_LHX2Y2);

    }

    *tensor = A * polar + B * (*tensor);
  }

}

}} // namespace psi::ccresponse
