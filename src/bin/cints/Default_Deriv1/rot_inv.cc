/*! \file
    \ingroup CINTS
    \brief Routines for determining the rotational invariance of gradients.
*/
#include <cstdio>
#include <cstdlib>
#include <libipv1/ip_lib.h>
#include <cmath>
#include <masses.h>
#include <libint/libint.h>
#include <psifiles.h>

#include "defines.h"
#define EXTERN
#include "global.h"

namespace psi { namespace CINTS {
void check_rot_inv()
{
  int atom;
  double cross[] = {0.0, 0.0, 0.0} ;
  double mod_cross;

  for(atom=0; atom<Molecule.num_atoms; atom++) {
    cross[0] += (Molecule.centers[atom].y*Grad[atom][2] - Molecule.centers[atom].z*Grad[atom][1]);
    cross[1] += (Molecule.centers[atom].z*Grad[atom][0] - Molecule.centers[atom].x*Grad[atom][2]);
    cross[2] += (Molecule.centers[atom].x*Grad[atom][1] - Molecule.centers[atom].y*Grad[atom][0]);
  }

  mod_cross = sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);

  if (mod_cross > ROT_INV_TOLER) {
    printf("  Rotational invariance of the energy derivative is violated!\n\n");
    fprintf(outfile,"  Rotational invariance of the energy derivative is violated!\n");
    fprintf(outfile,"  |X cross Grad| = %15.12lf\n\n",mod_cross);
    abort();
  }
  else {
    fprintf(outfile,"  Rotational invariance condition satisfied.\n");
    fprintf(outfile,"  |X cross Grad| = %15.12lf   (it is the accuracy of the computed forces)\n",mod_cross);
    fprintf(outfile,"  So long..\n\n");
  }
  
  return;
}
}}
