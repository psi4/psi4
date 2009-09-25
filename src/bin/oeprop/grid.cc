/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "includes.h"
#include "globals.h"

namespace psi { namespace oeprop {

extern void compute_grid_oeprops();
extern void compute_grid_dens_2d();
extern void compute_grid_dens_3d();
extern void compute_grid_mos();


void compute_grid() {

  switch (grid) {
  case 1:
      compute_grid_oeprops();
      break;
      
  case 2:
  case 3:
  case 4:
      compute_grid_dens_2d();
      break;

  case 5:
      compute_grid_mos();
      break;

  case 6:
      compute_grid_dens_3d();
      break;
  }

  return;
}

}} // namespace psi::oeprop
