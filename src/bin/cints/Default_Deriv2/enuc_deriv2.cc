/*! \file enuc_deriv2.cc
    \ingroup (CINTS)
    \brief Second derivatives of the nuclear repulsion
*/
#include<cstdio>

#include<libipv1/ip_lib.h>
#include<cmath>
#include<libciomr/libciomr.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"symmetrize.h"
#include<Tools/prints.h>

namespace psi { namespace CINTS {
  //! Second derivatives of the nuclear repulsion
void enuc_deriv2()
{
  int i, j;
  int coord_ix, coord_iy, coord_iz, coord_jx, coord_jy, coord_jz;
  double x, y, z, x2, y2, z2;
  double r, r2, r5, pfac;
  double **hess_enuc;

  if (Molecule.num_atoms == 0)
    return;

  hess_enuc = block_matrix(Molecule.num_atoms*3,Molecule.num_atoms*3);
  for(i=1; i<Molecule.num_atoms; i++) {
    coord_ix = i*3;
    coord_iy = coord_ix+1;
    coord_iz = coord_ix+2;
    for(j=0; j<i; j++) {
      coord_jx = j*3;
      coord_jy = coord_jx+1;
      coord_jz = coord_jx+2;

      x = (Molecule.centers[i].x-Molecule.centers[j].x);
      y = (Molecule.centers[i].y-Molecule.centers[j].y);
      z = (Molecule.centers[i].z-Molecule.centers[j].z);
      x2 = x*x; y2 = y*y; z2 = z*z;
      r2 = x2 + y2 + z2;
      r = sqrt(r2);
      r5 = r2*r2*r;
      pfac = Molecule.centers[i].Z_nuc*Molecule.centers[j].Z_nuc/r5;

      hess_enuc[coord_ix][coord_ix] += pfac*(3*x2 - r2);
      hess_enuc[coord_iy][coord_iy] += pfac*(3*y2 - r2);
      hess_enuc[coord_iz][coord_iz] += pfac*(3*z2 - r2);
      hess_enuc[coord_ix][coord_iy] += pfac*3*x*y;
      hess_enuc[coord_ix][coord_iz] += pfac*3*x*z;
      hess_enuc[coord_iy][coord_iz] += pfac*3*y*z;

      hess_enuc[coord_jx][coord_jx] += pfac*(3*x2 - r2);
      hess_enuc[coord_jy][coord_jy] += pfac*(3*y2 - r2);
      hess_enuc[coord_jz][coord_jz] += pfac*(3*z2 - r2);
      hess_enuc[coord_jx][coord_jy] += pfac*3*x*y;
      hess_enuc[coord_jx][coord_jz] += pfac*3*x*z;
      hess_enuc[coord_jy][coord_jz] += pfac*3*y*z;

      hess_enuc[coord_ix][coord_jx] -= pfac*(3*x*x-r2);
      hess_enuc[coord_ix][coord_jy] -= pfac*3*x*y;
      hess_enuc[coord_ix][coord_jz] -= pfac*3*x*z;
      hess_enuc[coord_iy][coord_jx] -= pfac*3*y*x;
      hess_enuc[coord_iy][coord_jy] -= pfac*(3*y*y-r2);
      hess_enuc[coord_iy][coord_jz] -= pfac*3*y*z;
      hess_enuc[coord_iz][coord_jx] -= pfac*3*z*x;
      hess_enuc[coord_iz][coord_jy] -= pfac*3*z*y;
      hess_enuc[coord_iz][coord_jz] -= pfac*(3*z*z-r2);
    }
  }

  symmetrize_hessian(hess_enuc);
  if (UserOptions.print_lvl >= PRINT_OEDERIV)
    print_atommat("Nuclear repulsion component of the molecular Hessian (a.u.)",hess_enuc);

  add_mat(Hess,hess_enuc,Hess,Molecule.num_atoms*3,Molecule.num_atoms*3);
  free_block(hess_enuc);
  
  return;
}
};};
