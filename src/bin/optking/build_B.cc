/*! \file
    \ingroup OPTKING
    \brief Build B matrix for salcs
*/

#include <libciomr/libciomr.h>
#include "intcos.h"
#include "salc.h"

namespace psi { namespace optking {

extern double **unit_matrix(int dim);
extern double **mass_matrix(int dim, double *masses);

/*
double **old_build_B(const Salc_set & salcs, const Intcos & simples) {
  
  int ncarts = 3*salcs.get_natom();
  int nsalcs = salcs.size();
  int nsimples = simples.size();

  double **B = block_matrix(nsalcs,ncarts);

  mmult(salcs.matrix,0,simples.Bsimp,0,B,0,nsalcs,nsimples,ncarts,0);

  return B;
}
*/

/*! \file
    \ingroup OPTKING
    \brief This function computes G = BuB^t from B
      if use_masses is true u is a diagonal matrix of inverse masses 
      otherwise u is the unit matrix
*/

double **build_G(double **B, int nrows, int ncols, bool use_masses, double *masses) {
  double **u;
  double **G = block_matrix(nrows,nrows);
  double **tmat = block_matrix(ncols,nrows);

  if (use_masses)
    u = mass_matrix(ncols, masses);
  else
    u = unit_matrix(ncols);

  mmult(u,0,B,1,tmat,0,ncols,ncols,nrows,0);
  mmult(B,0,tmat,0,G,0,nrows,ncols,nrows,0);

  free_block(u);
  free_block(tmat);
  return G;
}

}}
