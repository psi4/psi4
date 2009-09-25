/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "includes.h"
#include "prototypes.h"
#include "globals.h"

#include <vdw_radii.h>

#define VDW_SCALE 0.5

namespace psi { namespace oeprop {

void compute_connectivity()
{
  int i, j;
  int zi, zj;
  double radius_i, radius_j, rij;

  connectivity = init_int_matrix(natom,natom);
  
  for(i=0;i<natom;i++) {
    zi = (int)zvals[i];
    radius_i = (zi <= LAST_VDW_RADII_INDEX) ? atomic_vdw_radii[zi] : atomic_vdw_radii[0];
    radius_i /= _bohr2angstroms;
    for(j=0;j<i;j++) {
      zj = (int)zvals[j];
      radius_j = (zj <= LAST_VDW_RADII_INDEX) ? atomic_vdw_radii[zj] : atomic_vdw_radii[0];
      radius_j /= _bohr2angstroms;

      rij = sqrt((geom[i][0] - geom[j][0])*(geom[i][0] - geom[j][0]) +
		 (geom[i][1] - geom[j][1])*(geom[i][1] - geom[j][1]) +
		 (geom[i][2] - geom[j][2])*(geom[i][2] - geom[j][2]));

      if (rij <= VDW_SCALE*(radius_i + radius_j))
	connectivity[i][j] = connectivity[j][i] = 1;
    }
  }

  return;
}

}} // namespace psi::oeprop
