/*! \file rad_extent.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include <libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

#define MAX_ITER 100

namespace psi { namespace CINTS {

/*!------------------------------------------------------
  Use Newton-Raphson method to find maximum r for which
  the radial parts of basis functions drop below
  some threshold.
 ------------------------------------------------------*/

void init_rad_extent(double thresh)
{
  int i;
  int iter;
  int shell, prim, first_prim, last_prim, am;
  const double r0 = 4.0;        /* Start at 4.0 bohr ---*/
  double func, dfuncdr, sum, dsumdr;
  double r, r_new, tmp;

  BasisSet.thresh = thresh;
  
  for (shell=0;shell<BasisSet.num_shells;shell++) {
    r_new = r = r0;
    first_prim = BasisSet.shells[shell].fprim - 1;
    last_prim = first_prim + BasisSet.shells[shell].n_prims - 1;
    am = BasisSet.shells[shell].am-1;
    iter = 0;
    do {
      iter++;
      r = r_new;
      /*--------------------------------------
	evaluate the radial part of the shell
	and its first derivative
       --------------------------------------*/
      func = -thresh;
      dfuncdr = 0;
      for(prim=first_prim;prim<=last_prim;prim++) {
	tmp = BasisSet.cgtos[prim].ccoeff[am]*exp(-BasisSet.cgtos[prim].exp*(r*r));
	switch(am) {
	case 0:
	  func += tmp;
	  dfuncdr += -2.0*r*BasisSet.cgtos[prim].exp*tmp;
	  break;
	case 1:
	  func += r*tmp;
	  dfuncdr += (1.0-2.0*r*r*BasisSet.cgtos[prim].exp)*tmp;
	  break;
	default:
	  for(i=1;i<am;i++)
	    tmp *= r;
	  func += r*tmp;
	  dfuncdr += (1.0-2.0*r*r*BasisSet.cgtos[prim].exp)*tmp;
	}
      }

      /*--- We actually want to look at the
	absolute value of the basis function ---*/
      if (func+thresh < 0.0) {
	func = -(func+thresh)-thresh;
	dfuncdr *= -1.0;
      }

      /*--- If the function is too tight - reduce r ---*/
      if (func+thresh == 0.0 && dfuncdr == 0.0) {
	r_new /= 2.5;
	continue;
      }

      /*--- If before or at the maximum (for p-functions and higher)
	    step hopefully past it ---*/
      if (dfuncdr >= 0.0) {
	r_new *= 2.5;
	continue;
      }
	
      r_new = r - func/dfuncdr;
      if ( r_new <= 0.0 ) {
	r_new = r / 2.0;
      }

      if (iter > MAX_ITER)
	throw std::domain_error("Too many iterations while computing radial extents");
      
    } while (fabs(func/thresh) >= SOFT_ZERO);

    BasisSet.shells[shell].rad_extent = r_new;
    if (UserOptions.print_lvl > PRINT_DEBUG)
      fprintf(outfile,"Shell# = %d    Radial extent = %lf\n",shell,r_new);
      
  }

  return;
}


};};
