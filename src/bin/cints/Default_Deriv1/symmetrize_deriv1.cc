/*! \file symmetrize_deriv1.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/

#include <cstdio>

#include <cmath>
#include <vector>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libint/libint.h>

#include "defines.h"
#define EXTERN
#include "global.h"
using namespace std;

namespace psi { namespace CINTS {

void symmetrize_deriv1()
{
  int atom, atom2;
  int symop;
  double grad[3];
  vector<vector<double> > symm_grad(Molecule.num_atoms);

  for(atom=0; atom<Molecule.num_atoms; atom++) {
    vector<double> vdummy(3);
    symm_grad[atom]=vdummy;
    grad[0] = grad[1] = grad[2] = 0.0;
    for(symop=0;symop<Symmetry.nirreps;symop++) {
      atom2 = Symmetry.ict[symop][atom]-1;
      grad[0] += Symmetry.cartrep[symop][0]*Grad[atom2][0] +
	Symmetry.cartrep[symop][1]*Grad[atom2][1] +
	Symmetry.cartrep[symop][2]*Grad[atom2][2];
      grad[1] += Symmetry.cartrep[symop][3]*Grad[atom2][0] +
	Symmetry.cartrep[symop][4]*Grad[atom2][1] +
	Symmetry.cartrep[symop][5]*Grad[atom2][2];
      grad[2] += Symmetry.cartrep[symop][6]*Grad[atom2][0] +
	Symmetry.cartrep[symop][7]*Grad[atom2][1] +
	Symmetry.cartrep[symop][8]*Grad[atom2][2];
    }
    symm_grad[atom][0] = grad[0]/Symmetry.nirreps;
    symm_grad[atom][1] = grad[1]/Symmetry.nirreps;
    symm_grad[atom][2] = grad[2]/Symmetry.nirreps;
  }

  for(atom=0;atom<Molecule.num_atoms;atom++) {
    Grad[atom][0] = symm_grad[atom][0];
    Grad[atom][1] = symm_grad[atom][1];
    Grad[atom][2] = symm_grad[atom][2];
  }

  return;
}
};};
