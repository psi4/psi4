/*! \file enuc_deriv1.cc
    \ingroup (CINTS)
    \brief Compute the energy derivative of the nuclear repulsion.
*/
#include<stdio.h>
#include<libipv1/ip_lib.h>
#include<cmath>
#include<libciomr/libciomr.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include <Tools/prints.h>

namespace psi { namespace CINTS {
  //! Energy derivative of the nuclear repulsion.
void enuc_deriv1()
{
  int i, j;
  double rij2, tmp;
  double **grad_enuc;

  if (Molecule.num_atoms == 0)
    return;

  grad_enuc = block_matrix(Molecule.num_atoms,3);
  for(i=0; i<Molecule.num_atoms; i++)
    for(j=0; j<Molecule.num_atoms; j++)
      if (j != i) {
	rij2 =  (Molecule.centers[i].x-Molecule.centers[j].x)*(Molecule.centers[i].x-Molecule.centers[j].x);
	rij2 += (Molecule.centers[i].y-Molecule.centers[j].y)*(Molecule.centers[i].y-Molecule.centers[j].y);
	rij2 += (Molecule.centers[i].z-Molecule.centers[j].z)*(Molecule.centers[i].z-Molecule.centers[j].z);
	tmp = Molecule.centers[i].Z_nuc*Molecule.centers[j].Z_nuc/(rij2*sqrt(rij2));
	grad_enuc[i][0] -= (Molecule.centers[i].x-Molecule.centers[j].x)*tmp;
	grad_enuc[i][1] -= (Molecule.centers[i].y-Molecule.centers[j].y)*tmp;
	grad_enuc[i][2] -= (Molecule.centers[i].z-Molecule.centers[j].z)*tmp;
      }

  if (UserOptions.print_lvl >= PRINT_OEDERIV)
    print_atomvec("Nuclear repulsion component of the forces (a.u.)",grad_enuc);

  add_mat(Grad,grad_enuc,Grad,Molecule.num_atoms,3);
  free_block(grad_enuc);
  
  return;
}
};};
