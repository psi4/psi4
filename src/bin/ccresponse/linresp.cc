/*! \file
    \ingroup CCRESPONSE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
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
  for(int j=CC_TMP; j <= CC_TMP11; j++) {
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
        if(!strcmp(params.wfn,"CC2")) {
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
        if(!strcmp(params.wfn,"CC2")) {
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
      fprintf(outfile, "\n\tLinresp tensor <<%s;%s>>\n", pert_x, pert_y);
      fprintf(outfile, "\tpolar_LCX    = %20.12f\n", polar_LCX);
      if(!strcmp(params.wfn,"CC2"))
        fprintf(outfile, "\tpolar_HXY    = %20.12f\n", polar_HXY);
      fprintf(outfile, "\tpolar_LHX1Y1 = %20.12f\n", polar_LHX1Y1);
      fprintf(outfile, "\tpolar_LHX1Y2 = %20.12f\n", polar_LHX1Y2);
      fprintf(outfile, "\tpolar_LHX2Y2 = %20.12f\n", polar_LHX2Y2);
      fflush(outfile);
    }

    *tensor = A * polar + B * (*tensor);
  }

}

}} // namespace psi::ccresponse
