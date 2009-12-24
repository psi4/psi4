/*! \file
    \ingroup OPTKING
    \brief This function computes G via BuB^t where u is a diagonal matrix
    of inverse masses.
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

namespace psi { namespace optking {

double **compute_G(double **B, int num_intcos, const cartesians & carts) {
  double **u, **G, **temp_mat, *masses;
  int i, dim_carts;

  dim_carts = 3*carts.get_natom();
  masses = carts.get_mass();
  u = mass_mat(masses);
  free_array(masses);
  //u = unit_matrix(dim_carts);

  G = init_matrix(num_intcos,num_intcos);
  temp_mat = init_matrix(dim_carts,num_intcos);

  opt_mmult(u,0,B,1,temp_mat,0,dim_carts,dim_carts,num_intcos,0);
  opt_mmult(B,0,temp_mat,0,G,0,num_intcos,dim_carts,num_intcos,0);

  free_matrix(u);
  free_matrix(temp_mat);

  return G;
}

}}
