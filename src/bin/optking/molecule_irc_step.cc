/*! \file molecule_irc_step.cc
    \ingroup optking
    \Gonzalez Schlegel (1990) Reaction Path Following Algorithm
*/

#include "molecule.h"

#include <cmath>
#include <iostream>
#include <sstream>

#include "linear_algebra.h"
#include "print.h"
#include "atom_data.h"
#include "physconst.h"
#include "IRC_data.h"

#define EXTERN
#include "globals.h"

namespace opt {

// compute square root and inverse square root of matrix
void matrix_root(double **A, int dim, bool inverse);
// compute normalization factor for step
double step_N_factor(double **G, double *g, int Nintco);
// return the lowest eigenvector of the matrix (here, the Hessian)
double *lowest_evector(double **H, int Nintco);

// compute change in energy according to quadratic approximation
inline double DE_nr_energy(double step, double grad, double hess)
{
  return (step * grad + 0.5 * step * step * hess);
}

void IRC_DATA::point_converged(opt::MOLECULE &mol)
{
//  p_Opt_data->reset_iteration_to_size();
  sphere_step = 0;
  mol.irc_step();
}

// Gonzalez and Schlegel JCP 2154 (1989).
void MOLECULE::irc_step(void)
{
  bool at_TS = !(p_irc_data->size());          //are we starting from the transition state?
  bool at_FS = !(p_irc_data->sphere_step);     //FS for "first step" - are we starting a new
                                               //constrained optimization?
fprintf(outfile, "\ncoord,sphere: %i, %i\n", p_irc_data->size(), p_irc_data->sphere_step);
//Basic Information
  double s    = Opt_params.IRC_step_size;      //step size
  double *dq  = p_Opt_data->g_dq_pointer();    //internal coordinate change
  double *f_q = p_Opt_data->g_forces_pointer();//internal coordinate gradient
  double **H  = p_Opt_data->g_H_pointer();     //internal coordinate Hessian

//Dimensions
  int Nintco = g_nintco();
  int Natom = g_natom();
  int Ncart = 3 * Natom;

double *x_initial = g_geom_array();

//The G matrix, its square root, and its inverse square root
  double **G         = compute_G(1);                           //G = BUB^t
  double **rootG_reg = matrix_return_copy(G, Nintco, Nintco); //G^1/2
  double **rootG_inv = matrix_return_copy(G, Nintco, Nintco); //G^-1/2
  matrix_root(rootG_reg, Nintco, 0);
  matrix_root(rootG_inv, Nintco, 1);
/*
fprintf(outfile, "G matrix:\n");
print_matrix(outfile, G, Nintco, Nintco);
fprintf(outfile, "G^1/2   :\n");
print_matrix(outfile, rootG_reg, Nintco, Nintco);
fprintf(outfile, "G^-1/2  :\n");
print_matrix(outfile, rootG_inv, Nintco, Nintco);
*/
//mass-weighted Hessian matrix:
  double **H_m = init_matrix(Nintco, Nintco);
  double **T = init_matrix(Nintco, Nintco);
  // T = HG^1/2
  opt_matrix_mult(H, 0, rootG_reg, 0, T, 0, Nintco, Nintco, Nintco, 0);
  // H_m = (G^1/2)H(G^1/2) = (G^1/2)T
  opt_matrix_mult(rootG_reg, 0, T, 0, H_m, 0, Nintco, Nintco, Nintco, 0);
  free_matrix(T);

fprintf(outfile,"Mass-weighted Hessian\n");
print_matrix(outfile, H_m, Nintco, Nintco);

//information to calculate predicted energy change (DE_predicted)
  double DE_projected;
  double dq_n;    //step norm
  double *dq_u;   //unit vector in step direction
  double dq_g;    //gradient in step direction
  double dq_h;    //hessian in step direction

//Pivot Point and Initial Working Geometry:
//step along along normalized, mass-weighted v
  if(at_FS) {
cout << "First point of constrained optimization.\n";
    fprintf(outfile, "\tFirst point of constrained optimization.\n");
    //if starting from TS, follow lowest eigenvalued eigenvalued eigenvector
    //otherwise, follow the gradient (negative of the force vector)
    double *v;
    if (at_TS) {
cout << "stepping from TS\n";
      v = lowest_evector(H, Nintco);

      if(Opt_params.IRC_direction == OPT_PARAMS::FORWARD)
        fprintf(outfile, "Stepping in forward direction from TS.\n");
      else if (Opt_params.IRC_direction == OPT_PARAMS::BACKWARD) {
        fprintf(outfile, "Stepping in backward direction from TS.\n");
        for(int i=0; i<Nintco; i++)
          v[i] *= -1;
      }
    }
    else {
      fprintf(outfile, "Stepping along IRC using gradient.\n");
cout << "Stepping from previous point using gradient\n";
      v = init_array(Nintco);
      for(int i=0; i<Nintco; i++)
        v[i] = -f_q[i];
    }

    // calculate pivot point and next point
    double N = step_N_factor(G, v, Nintco);
    double *dq_pivot = init_array(Nintco);

    for(int i=0; i<Nintco; i++) {
      dq[i] = - array_dot(G[i], v, Nintco) * N * s;
      dq_pivot[i] = dq[i] / 2;
    }

fprintf(outfile, "\nVector to follow: \n");
print_array(outfile, v, Nintco);
fprintf(outfile, "\nDq to pivot point: \n");
print_array(outfile, dq_pivot, Nintco);
fprintf(outfile, "\nDq to next geometry: \n");
print_array(outfile, dq, Nintco);

// Check calculated step size
double **G_inv = symm_matrix_inv(G, Nintco, Nintco);
double *dq_G_inv = init_array(Nintco);

opt_matrix_mult(G_inv, 0, &dq, 1, &dq_G_inv, 1, Nintco, Nintco, 1, 0);
double dq_norm = sqrt( array_dot(dq_G_inv, dq, Nintco) );

fprintf(outfile, "\ndq_norm bf sqrt:  %20.15e\n", array_dot(dq_G_inv, dq, Nintco) );
fprintf(outfile, "\nChecking dq_norm: %20.10f\n", dq_norm);
free_array(dq_G_inv);
free_matrix(G_inv);

    double *q = intco_values();
    double *q_pivot = init_array(Nintco);
    for(int i=0; i<Nintco; i++)
      q_pivot[i] = q[i] + dq_pivot[i];
    double *x = g_geom_array();
    double *f_x = g_grad_array();
    double *f_q_copy = init_array(Nintco);
    array_copy(f_q, f_q_copy, Nintco);

    // q_pivot, q, x,f_q, f_x should NOT be freed; pointers are assigned in IRC_data
    p_irc_data->add_irc_point(p_irc_data->g_next_coord_step(), q_pivot, q, x, f_q_copy, f_x, g_energy());

    free_array(dq_pivot);

    // Do displacements for each fragment separately.
    for (int f=0; f<fragments.size(); ++f) {
      if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
        fprintf(outfile,"\tDisplacements for frozen fragment %d skipped.\n", f+1);
        continue;
      }
      fragments[f]->displace(&(dq[g_intco_offset(f)]), true, g_intco_offset(f));
    }
    // Do displacements for interfragment coordinates.
    double *q_target;
    for (int I=0; I<interfragments.size(); ++I) {
      q_target = interfragments[I]->intco_values();
      for (int i=0; i<interfragments[I]->g_nintco(); ++i)
        q_target[i] += dq[g_interfragment_intco_offset(I) + i];
       interfragments[I]->orient_fragment(q_target);
      free_array(q_target);
    }

    // Compute norm, unit vector, gradient and Hessian in the step direction
    // Save results to Opt_data.
    dq_n = sqrt( array_dot(dq, dq, Nintco) );

    dq_u = init_array(Nintco);
    array_copy(dq, dq_u, Nintco);
    array_normalize(dq_u, Nintco);

    dq_g = -1 * array_dot(f_q, dq_u, Nintco);

    dq_h = 0;
    for(int i=0; i<Nintco; i++)
      dq_h += dq_u[i] * array_dot(H[i], dq_u, Nintco);

    DE_projected = DE_nr_energy(dq_n, dq_g, dq_h);
fprintf(outfile,"DE_projected  %20.15lf\n", DE_projected);

    p_Opt_data->save_step_info(DE_projected, dq_u, dq_n, dq_g, dq_h);
    free_array(dq_u);

    // Increment to keep track of number of constrained optimization steps for this point.
    p_irc_data->sphere_step++;

    symmetrize_geom();
    free_array(v);
    return;
  }

// Constrained Optimization On Hypersphere:
// Step along hypersphere (radius s/2) to find succeeding point with the following properties:
//  |p| = |dq_k - dq_pivot| = s/2
//  (instantaneous acceleration vector Gv at dq_k) is parallel to (radial vector p)

//1. diagonalize H_m:
  double **V = matrix_return_copy(H_m, Nintco, Nintco);
  double *h = init_array(Nintco);
  opt_symm_matrix_eig(V, Nintco, h);

  if (Opt_params.print_lvl > 2) {
    for(int i=0; i<Nintco; i++) {
      fprintf(outfile, "Eigenvector %i of mass-weighted Hessian:\n", i);
      print_array(outfile, V[i], Nintco);
    }
  }
  fprintf(outfile, "Eigenvalues of the mass-weighted Hessian:\n");
  print_array(outfile, h, Nintco);

//2. Express p and g, in mass-weighted coordinates and in the eigenbasis of Hm, the mass-weighted Hessian.
  double *q = intco_values();
  double *q_pivot = p_irc_data->steps.back()->g_q_pivot();

fprintf(outfile, "\nRetrieved pivot point from IRC_data: \n");
print_array(outfile, q_pivot, Nintco);

  double *p = init_array(Nintco);  //vector step from pivot point
  double *p_m = init_array(Nintco);//mass-weighted
  double *p_h = init_array(Nintco);//in the basis of H_m
  for(int i=0; i<Nintco; i++)
    p[i] = q[i] - q_pivot[i];
  for(int i=0; i<Nintco; i++)
    p_m[i] = array_dot(rootG_inv[i], p, Nintco);
  for(int i=0; i<Nintco; i++)
    p_h[i] = array_dot(p_m, V[i], Nintco);

  double *g_m = init_array(Nintco);//mass-weighted gradient
  double *g_h = init_array(Nintco);//in the basis of H_m
  for(int i=0; i<Nintco; i++)
    g_m[i] = -array_dot(rootG_reg[i], f_q, Nintco);
  for(int i=0; i<Nintco; i++)
    g_h[i] = array_dot(g_m, V[i], Nintco);

fprintf(outfile, "\np:\n");
print_array(outfile, p, Nintco);
fprintf(outfile, "\np_m:\n");
print_array(outfile, p_m, Nintco);
fprintf(outfile, "\np_h:\n");
print_array(outfile, p_h, Nintco);
fprintf(outfile, "\ng_m:\n");
print_array(outfile, g_m, Nintco);
fprintf(outfile, "\ng_h:\n");
print_array(outfile, g_h, Nintco);

//3. solve equation 26 in Gonzalez & Schlegel (1990) for lambda based on current p-vector
  double lagrangian=1; // value of lagrangian expression
  double df_1; //1st derivative of lagrangian
  double df_2; //2nd derivative of lagrangian
  int j = 0;
  double old_lambda = -999;
  int lag_iter=0;

  // lagrange multiplier ; use previous value as first guess at next point
  static double lambda = 1.5 * h[0]; // Make 50% lower than the lowest, negative eval

  fprintf(outfile, "\tDetermining lagrangian multiplier for constrained minimization.\n");
  fprintf(outfile, "\t Iter     Multiplier     Lagrangian\n");

  while (fabs(lambda - old_lambda) > 1e-14)
  {
    double tval;
    double D;
    lagrangian = df_1 = df_2 = 0;
    old_lambda = lambda;

    for(int i=0; i<Nintco; i++) {
      D = h[i] - lambda; 
      tval = (h[i]*p_h[i] - g_h[i]) * (h[i]*p_h[i] - g_h[i]) / (D * D);
      lagrangian += tval;
      df_1 += 2 * tval / D ;
      df_2 += 6 * tval / (D * D);
    }
    lagrangian -= (s / 2) * (s / 2);

    fprintf(outfile, "\t%5d%15.3e%15.3e\n", lag_iter, lambda, lagrangian);

    //lambda += 2 * (lagrangian * df_1) / (lagrangian * df_2 - 2 * (df_1 * df_1));
    tval = - lagrangian / df_1 ;
    lambda += tval / (1 + 0.5 * tval * df_2 / df_1); // 'Halley's method'
    //lambda -= lagrangian / df_1; 'Newton's method'

    if (++lag_iter > 50)
      throw("Could not converge lagrangian for constrained minimization");
  }

//4. find dq_m based on equation 24 in Gonzalez & Schlegel (1990)
  double *dq_m = init_array(Nintco);

  // Compute (H_m - lambda * I)^(-1)
  for(int i=0; i<Nintco; i++)
    H_m[i][i] -= lambda;
  double **H_m_lag_inv = symm_matrix_inv(H_m, Nintco, 1);

  // Compute vector (g_m - lambda * p_m) and dot with rows of H_m_lag_inv to get -dq_m
  double *g_lag_p = init_array(Nintco);
  for(int i=0; i<Nintco; i++)
    g_lag_p[i] = g_m[i] - lambda * p_m[i];

  opt_matrix_mult(H_m_lag_inv, 0, &g_lag_p, 1, &dq_m, 1, Nintco, Nintco, 1, 0);
  free_array(g_lag_p);
  free_matrix(H_m_lag_inv);

//5. find dq ( = (G^1/2)dq_m) and do displacements
  opt_matrix_mult(rootG_reg, 0, &dq_m, 1, &dq, 1, Nintco, Nintco, 1, 0);

  // Do displacements for each fragment seperately.
  for (int f=0; f<fragments.size(); ++f) {
    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
      fprintf(outfile,"\tDisplacements for frozen fragment %d skipped.\n", f+1);
      continue;
    }
    fragments[f]->displace(&(dq[g_intco_offset(f)]), true, g_intco_offset(f));
  }
  // Do displacements for interfragment coordinates.
  double *q_target;
  for (int I=0; I<interfragments.size(); ++I) {
    q_target = interfragments[I]->intco_values();
    for (int i=0; i<interfragments[I]->g_nintco(); ++i)
      q_target[i] += dq[g_interfragment_intco_offset(I) + i];
     interfragments[I]->orient_fragment(q_target);
    free_array(q_target);
  }

  // Compute norm, unit vector, gradient and Hessian in the step direction
  dq_n = sqrt( array_dot(dq, dq, Nintco) );

  dq_u = init_array(Nintco);
  array_copy(dq, dq_u, Nintco);
  array_normalize(dq_u, Nintco);

  dq_g = -1 * array_dot(f_q, dq_u, Nintco);

  dq_h = 0;
  for(int i=0; i<Nintco; i++)
    dq_h += dq_u[i] * array_dot(H[i], dq_u, Nintco);

  p_Opt_data->save_step_info(DE_projected, dq_u, dq_n, dq_g, dq_h);
  p_irc_data->sphere_step++;
  symmetrize_geom();

  free_matrix(rootG_reg);
  free_matrix(rootG_inv);
  free_matrix(H_m);
  free_matrix(V);
  free_array(h);
  free_array(p);
  free_array(p_m);
  free_array(p_h);
  free_array(g_m);
  free_array(g_h);
  free_array(dq_m);
  free_array(dq_u);
}

void matrix_root(double **A, int dim, bool inverse)
{
  double **V = matrix_return_copy(A, dim, dim);
  double *A_evals = init_array(dim);

  opt_symm_matrix_eig(V, dim, A_evals);

  /* for(int i=0; i<dim; i++) {
    fprintf(outfile, "G_evect %i:\n", i);
    print_array(outfile, V[i], dim);
  }
  fprintf(outfile, "G_evals:\n");
  print_array(outfile, A_evals, dim); */

  if (inverse) {
    for(int k=0; k<dim; k++)
      if(fabs(A_evals[k]) > Opt_params.redundant_eval_tol)
        A_evals[k] = 1.0/A_evals[k];
  }

  for(int k=0; k<dim; k++) {
    if(fabs(A_evals[k]) > 0)
      A_evals[k] = sqrt(A_evals[k]);
    else
      throw("matrix_root() : cannot take square root of negative eigenvalue.\n");
  }

  zero_matrix(A, dim, dim);

  for(int i=0; i<dim; i++)
    for(int j=0; j<dim; j++)
      for(int k=0; k<dim; k++)
        A[i][j] += V[k][i] * A_evals[k] * V[k][j];

  free_matrix(V);
  free_array(A_evals);
}

double step_N_factor(double **G, double *g, int Nintco)
{
  double N = 0;
  double Ggsum;

  for(int i=0; i<Nintco; i++) {

    Ggsum = 0;
    for(int j = i+1; j<Nintco; j++)
      Ggsum += G[i][j] * g[j];

    N += g[i] * (G[i][i] * g[i] + 2 * Ggsum);
  }

  N = 1 / sqrt(N);
  return N;
}

double *lowest_evector(double **H, int Nintco)
{
  double **V = matrix_return_copy(H, Nintco, Nintco);
  double *V_evals = init_array(Nintco);
  opt_symm_matrix_eig(V, Nintco, V_evals);

  double *min_evector = init_array(Nintco);
  double max_element = 0;

  for(int i=0; i<Nintco; i++) {
    min_evector[i] = V[0][i];
    if( fabs(min_evector[i]) > fabs(max_element) )
      max_element = min_evector[i];
  }

  for(int i=0; i<Nintco; i++)
    min_evector[i] *= max_element / fabs(max_element);

cout << "making largest element of eigenvector positive\n";

  free_matrix(V);
  free_array(V_evals);

  return min_evector;
}

}
