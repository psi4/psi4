/*! \file molecule.cc
    \ingroup optking
    \brief molecule class (really, molecular system class)
*/

#include "molecule.h"

#include <cmath>
#include <iostream>
#include <sstream>

#include "linear_algebra.h"
#include "print.h"
#include "atom_data.h"
#include "physconst.h"

#define EXTERN
#include "globals.h"

namespace opt {

// compute change in energy according to RFO approximation
inline double DE_rfo_energy(double rfo_t, double rfo_g, double rfo_h) {
  return (rfo_t * rfo_g + 0.5 * rfo_t * rfo_t * rfo_h)/(1 + rfo_t*rfo_t);
}

// Take Rational Function Optimization step
void MOLECULE::rfo_step(void) {
  int i, j;
  int dim = g_nintco();
  double tval, tval2;
  double *fq = p_Opt_data->g_forces_pointer();
  double **H = p_Opt_data->g_H_pointer();
  double *dq = p_Opt_data->g_dq_pointer();

  // build (lower-triangle of) RFO matrix and diagonalize
  double **rfo_mat = init_matrix(dim+1, dim+1);
  for (i=0; i<dim; ++i)
    for (j=0; j<=i; ++j)
      rfo_mat[i][j] = H[i][j];

  for (i=0; i<dim; ++i)
    rfo_mat[dim][i] = - fq[i];

  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile,"RFO mat\n");
    print_matrix(outfile, rfo_mat, dim+1, dim+1);
  }

  double *lambda = init_array(dim+1);
  opt_symm_matrix_eig(rfo_mat, dim+1, lambda);
  if (Opt_params.print_lvl >= 2) {
    fprintf(outfile,"RFO eigenvalues/lambdas\n");
    print_matrix(outfile, &(lambda), 1, dim+1);
  }

  // Do intermediate normalization.  
  // RFO paper says to scale eigenvector to make the last element equal to 1.
  // During the course of an optimization some evects may appear that are bogus leads
  // - the root following can avoid them. 
  for (i=0; i<dim+1; ++i) {
    tval = rfo_mat[i][dim];
    if (fabs(tval) > Opt_params.rfo_normalization_min) {
      for (j=0;j<dim+1;++j)
        rfo_mat[i][j] /= rfo_mat[i][dim];
    }
  }
  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile,"RFO eigenvectors (rows)\n");
    print_matrix(outfile, rfo_mat, dim+1, dim+1);
  }

  int rfo_root, f;
  double rfo_eval;
  double *rfo_u;      // unit vector in step direction
  double rfo_dqnorm;   // norm of step
  double rfo_g;         // gradient in step direction
  double rfo_h;         // hessian in step direction
  double DE_projected; // projected energy change by quadratic approximation

  // *** choose which RFO eigenvector to use
  // if not root following, then use rfo_root'th lowest eigenvalue; default is 0 (lowest)
  if ( (!Opt_params.rfo_follow_root) || (p_Opt_data->g_iteration() == 1)) {
    rfo_root = Opt_params.rfo_root;
    fprintf(outfile,"\tFollowing RFO solution %d.\n", rfo_root);
  }
  else { // do root following
    double * rfo_old_evect = p_Opt_data->g_rfo_eigenvector_pointer();
    tval = 0;
    for (i=0; i<dim; ++i) { // dot only within H block, excluding rows and columns approaching 0
      tval2 = array_dot(rfo_mat[i], rfo_old_evect,dim);
      if (tval2 > tval) {
        tval = tval2;
        rfo_root = i;
      }
    }
    fprintf(outfile,"RFO vector %d has maximal overlap with previous step\n", rfo_root+1);
  }
  p_Opt_data->set_rfo_eigenvector(rfo_mat[rfo_root]);

  // print out lowest energy evects
  if (Opt_params.print_lvl >= 2) {
    for (i=0; i<dim+1; ++i) {
      if ((lambda[i] < 0.0) || (i <rfo_root)) {
        fprintf(outfile,"RFO eigenvalue %d: %15.10lf (or 2*%-15.10lf)\n", i+1, lambda[i],lambda[i]/2);
        fprintf(outfile,"eigenvector:\n");
        print_matrix(outfile, &(rfo_mat[i]), 1, dim+1);
      }
    }
  }
  free_array(lambda);

  for (j=0; j<dim; ++j)
    dq[j] = rfo_mat[rfo_root][j]; // leave out last column

  // Zero steps for frozen fragment
  for (f=0; f<fragments.size(); ++f) {
    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
      fprintf(outfile,"\tZero'ing out displacements for frozen fragment %d\n", f+1);
      for (i=0; i<fragments[f]->g_nintco(); ++i)
        dq[ g_intco_offset(f) + i ] = 0.0;
    }
  }

  apply_intrafragment_step_limit(dq);
  //check_intrafragment_zero_angles(dq);

  // get norm |dq| and unit vector in the step direction
  rfo_dqnorm = sqrt( array_dot(dq, dq, dim) );
  rfo_u = init_array(dim);
  array_copy(rfo_mat[rfo_root], rfo_u, dim);
  array_normalize(rfo_u, dim);
  free_matrix(rfo_mat);

  // get gradient and hessian in step direction
  rfo_g = -1 * array_dot(fq, rfo_u, dim);
  rfo_h = 0;
  for (i=0; i<dim; ++i)
    rfo_h += rfo_u[i] * array_dot(H[i], rfo_u, dim);

  DE_projected = DE_rfo_energy(rfo_dqnorm, rfo_g, rfo_h);
  fprintf(outfile,"\tProjected energy change by RFO approximation: %20.10lf\n", DE_projected);

/* Test step sizes
  double *x_before = g_geom_array();
  double **G = compute_G(true);
  double **G_inv = symm_matrix_inv(G, dim, 1);
  free_matrix(G);
  double *Gdq = init_array(dim);
  opt_matrix_mult(G_inv, 0, &dq, 1, &Gdq, 1, dim, dim, 1, 0);
  free_matrix(G_inv);
  double N = array_dot(dq, Gdq, dim);
  N = sqrt(N);
  free_array(Gdq);
  fprintf(outfile,"Step-size in mass-weighted internal coordinates: %20.10lf\n", N);
*/

  // do displacements for each fragment separately
  for (f=0; f<fragments.size(); ++f) {
    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
      fprintf(outfile,"\tDisplacements for frozen fragment %d skipped.\n", f+1);
      continue;
    }
    fragments[f]->displace(&(dq[g_intco_offset(f)]), &(fq[g_intco_offset(f)]), g_atom_offset(f));
  }

  // do displacements for interfragment coordinates
  for (int I=0; I<interfragments.size(); ++I) {
    if (interfragments[I]->is_frozen() || Opt_params.freeze_interfragment) {
      fprintf(outfile,"\tDisplacements for frozen interfragment %d skipped.\n", I+1);
      continue;
    }
    interfragments[I]->orient_fragment( &(dq[g_interfragment_intco_offset(I)]),
                                        &(fq[g_interfragment_intco_offset(I)]) );
  }

#if defined(OPTKING_PACKAGE_QCHEM)
  // fix rotation matrix for rotations in QCHEM EFP code
  for (int I=0; I<efp_fragments.size(); ++I)
    efp_fragments[I]->displace( I, &(dq[g_efp_fragment_intco_offset(I)]) );
#endif

  symmetrize_geom(); // now symmetrize the geometry for next step

/* Test step sizes
  double *x_after = g_geom_array();
  double *masses = g_masses();

  double sum = 0.0;
  for (int i=0; i<g_natom(); ++i)
    for (int xyz=0; xyz<3; ++xyz)
      sum += (x_before[3*i+xyz] - x_after[3*i+xyz]) * (x_before[3*i+xyz] - x_after[3*i+xyz])
               * masses[i];

  sum = sqrt(sum);
  fprintf(outfile,"Step-size in mass-weighted cartesian coordinates [bohr (amu)^1/2] : %20.10lf\n", sum);
  free_array(x_after);
  free_array(x_before);
  free_array(masses);
*/

  // save values in step data
  p_Opt_data->save_step_info(DE_projected, rfo_u, rfo_dqnorm, rfo_g, rfo_h);

  free_array(rfo_u);
  fflush(outfile);

} // end take RFO step

}

