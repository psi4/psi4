/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*! \file molecule_irc_step.cc
    \ingroup optking
    \Gonzalez Schlegel (1990) Reaction Path Following Algorithm
*/

#include "molecule.h"

#include <iostream>
#include <sstream>

#include "linear_algebra.h"
#include "print.h"
#include "atom_data.h"
#include "physconst.h"
#include "IRC_data.h"

#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
#endif

namespace opt {

// compute normalization factor for step
double step_N_factor(double **G, double *g, int Nintco);
// return the lowest eigenvector of the matrix (here, the Hessian)
double *lowest_evector(double **H, int Nintco);

void GS_interpolation(double *p, double *p_0, double *g, double *g_0, double s, int dim);
void interpolation(double *p, double *p_0, double *g, double *g_0, double s, int dim);

double lag_function(double l, double *f, double *h, double *p, double *g, int dim, double s);
// compute change in energy according to quadratic approximation
inline double DE_nr_energy(double step, double grad, double hess)
{
  return (step * grad + 0.5 * step * step * hess);
}

void IRC_DATA::point_converged(opt::MOLECULE &mol)
{
  if(!go)
cout << "we made it.";
  psi::outfile->Printf("\tPoint is converged. Setting sphere_step to 0, and calling irc_step().\n\n");
  if(steps.size() > 1)
  {
    double f_dot = array_dot(steps[steps.size()-1]->g_f_q(), steps[steps.size()-2]->g_f_q(), mol.g_nintco());
    psi::outfile->Printf("\nforce vector - current dotted with previous: %f\n", f_dot);
  }

  sphere_step = 0;

  mol.irc_step();
}

void IRC_DATA::progress_report(opt::MOLECULE &mol) 
{
  double DE;
  int dim = mol.g_nintco();
  int blocks = 4;
  int sign = 1;

  if(Opt_params.IRC_direction == OPT_PARAMS::BACKWARD)
    sign = -1;

//Printing Energies and Energy Changes for Each Step
  psi::outfile->Printf(  "\t----------------------------------------------");
  psi::outfile->Printf("\n\t            ****      IRC Report      ****\n");
  psi::outfile->Printf(  "\t----------------------------------------------\n");
  psi::outfile->Printf(  "\t Step    Energy              Change in Energy");
  psi::outfile->Printf("\n");
  psi::outfile->Printf(  "\t----------------------------------------------\n");
  for (int i=0; i<steps.size(); ++i)
  {
    if (i == 0) DE = g_step(i).g_energy();
    else DE = g_step(i).g_energy() - g_step(i-1).g_energy();

    psi::outfile->Printf("\t %3d %18.12lf  %18.12lf\n", i, g_step(i).g_energy(), DE);
  }
  psi::outfile->Printf(  "\t----------------------------------------------\n");
  psi::outfile->Printf("\n");

//Printing Internal Coordinates for Each step
  psi::outfile->Printf("\t--------------------------------------");
  for(int i=0; i<(dim/blocks)*blocks; i++)
  {
    psi::outfile->Printf("-------------");
  }
  psi::outfile->Printf("\n");
  psi::outfile->Printf("\t              ****     IRC Steps     ****\n");
  psi::outfile->Printf("\t--------------------------------------");
  for(int i=0; i<(dim/blocks)*blocks; i++)
  {
    psi::outfile->Printf("-------------");
  }

  for(int j=0; j < dim/blocks; j++)
  {
    psi::outfile->Printf("\n\t        |          Distance         |\n");
    psi::outfile->Printf("  \t Step   | Step    Arc       Line    |");
    for(int i = (j*blocks); i < ((j+1)* blocks); i++)
    {
      psi::outfile->Printf("    Coord %3d", i);
    }
    psi::outfile->Printf("\n");
    psi::outfile->Printf("\t--------------------------------------");
    for(int i = (j*blocks); i < ((j+1)* blocks); i++)
    {
      psi::outfile->Printf("-------------");
    }
    psi::outfile->Printf("\n");
    for (int i=0; i<steps.size(); ++i)
    {
      psi::outfile->Printf("\t %3d %9.2lf %9.5lf  %9.5lf   ", i, sign*g_step(i).g_step_dist(), sign*g_step(i).g_arc_dist(), sign*g_step(i).g_line_dist());
      for(int k = (j*blocks); k < ((j+1)*blocks); k++)
        psi::outfile->Printf("%13.8f",g_step(i).g_q()[k]);
      psi::outfile->Printf("\n");
    }
    psi::outfile->Printf("\t--------------------------------------");
    for(int i = (j*blocks); i < ((j+1)* blocks); i++)
    {
      psi::outfile->Printf("-------------");
    }
  }
  if(dim % blocks != 0)
  {
    psi::outfile->Printf("\n\t        |          Distance         |\n");
    psi::outfile->Printf("  \t Step   | Step    Arc       Line    |");
    for(int i = (dim - (dim % blocks)); i < dim; i++)
    {
      psi::outfile->Printf("    Coord %3d", i);
    }
    psi::outfile->Printf("\n");
    psi::outfile->Printf("\t--------------------------------------");
    for(int i = (dim - (dim % blocks)); i < dim; i++)
    {
      psi::outfile->Printf("-------------");
    }
    psi::outfile->Printf("\n");
    for (int i=0; i<steps.size(); ++i)
    {
      psi::outfile->Printf("\t %3d %9.2lf %9.5lf  %9.5lf   ", i, sign*g_step(i).g_step_dist(), sign*g_step(i).g_arc_dist(), sign*g_step(i).g_line_dist());
      for(int k = (dim - (dim % blocks)); k < dim; k++)
        psi::outfile->Printf("%13.8f",g_step(i).g_q()[k]);
      psi::outfile->Printf("\n");
    }
    psi::outfile->Printf("\t--------------------------------------");
    for(int i = (dim - (dim % blocks)); i < dim; i++)
    {
      psi::outfile->Printf("-------------");
    }
  }

  psi::outfile->Printf("\n");
  psi::outfile->Printf("\n");

  mol.print_intcos("outfile");
}


// See Gonzalez and Schlegel JCP 2154 (1989).
void MOLECULE::irc_step(void)
{
  // Are we at the TS?  at_TS
  bool at_TS = !(p_irc_data->size());        
  if (at_TS) psi::outfile->Printf("\n\tIRC_DATA is empty, so we are at the transition state.\n");

  // Is this one the first step toward a new path point?  at_FS 
  bool at_FS = !(p_irc_data->sphere_step);    
  if (at_FS) psi::outfile->Printf("\tIRC_DATA->sphere_step == 0, so now is time \
for a first step toward new point on path.\n");

  psi::outfile->Printf("\n\tRxn path step %d, constrained step %d\n", p_irc_data->size(),
    p_irc_data->sphere_step);

  double s    = Opt_params.IRC_step_size;      //step size
  p_irc_data->step_length = s;
  double *dq  = p_Opt_data->g_dq_pointer();    //internal coordinate change
  double *f_q = p_Opt_data->g_forces_pointer();//internal coordinate gradient
  int Nintco = g_nintco();
  int Natom = g_natom();
  int Ncart = 3 * Natom;
  bool answer = 1;

  int opt_iter = p_Opt_data->g_iteration() - 1;
  if(opt_iter > 2 && p_irc_data->in_min_range)
  {
cout << "DE01:    " << p_Opt_data->g_energy(opt_iter) - p_Opt_data->g_energy(opt_iter - 2) << "\n";
cout << "DE02:    " << p_Opt_data->g_energy(opt_iter - 1) - p_Opt_data->g_energy(opt_iter - 3) << "\n";
    if(    fabs(p_Opt_data->g_energy(opt_iter)     - p_Opt_data->g_energy(opt_iter - 2)) < 0.01e-03
        && fabs(p_Opt_data->g_energy(opt_iter - 1) - p_Opt_data->g_energy(opt_iter - 3)) < 0.01e-03 )
    {
      p_irc_data->go = 0;
    }
  }

    if (p_irc_data->sphere_step == 1)
    {

      double *u_f_q = init_array(Nintco);
      double *u_f_q_0 = init_array(Nintco);
      array_copy(f_q, u_f_q, Nintco);
      array_copy(p_irc_data->g_f_q(), u_f_q_0, Nintco);
      array_normalize(u_f_q, Nintco);
      array_normalize(u_f_q_0, Nintco);
      double u_f_q_dot = array_dot(u_f_q, u_f_q_0, Nintco);
      psi::outfile->Printf("\ninternal force vector dot - current with previous: %20.15f\n", u_f_q_dot);
cout << "u_f_q_dot: " << u_f_q_dot << "\n";
cout << "line_dist: " << p_irc_data->g_line_dist(p_irc_data->size()-1);

      if(u_f_q_dot < -0.5)
      {
        p_irc_data->in_min_range = 1;
      }

      if(p_irc_data->size() > 3)
      {
        if (u_f_q_dot < -0.7 ||
            (
              fabs(p_irc_data->g_line_dist(p_irc_data->size()-1) - p_irc_data->g_line_dist(p_irc_data->size()-2)) < s*10e-03
            )
           )
        {
          cout << "\nHouston, we've found a minimum!\n";

          if(Opt_params.IRC_stop == OPT_PARAMS::ASK)
          {
            cout << "Would you like to proceed? (1=yes, 0=no):";
            cin >> answer;
          }
          if(Opt_params.IRC_stop == OPT_PARAMS::STOP || !answer)
          {
            p_irc_data->go = 0;
          }
        }
      }

      free_array(u_f_q);
      free_array(u_f_q_0);
    }
/*    if(!p_irc_data->go)
    {
      double *q = intco_values();
      double *q_pivot = init_array(Nintco);

      double *x = g_geom_array();
      double *f_x = g_grad_array();
      array_scm(f_x, -1, Ncart);
      double *f_q_copy = init_array(Nintco);
      array_copy(f_q, f_q_copy, Nintco);

      p_irc_data->add_irc_point(p_irc_data->g_next_coord_step(), q_pivot, q, x, f_q_copy, f_x, g_energy(),
                                p_irc_data->step_length, p_irc_data->arc_length, p_irc_data->line_length);

      return;
    }
*/


  //The G matrix, its square root, and its inverse square root
  double **G         = compute_G(1);                          //G = BUB^t

  double **rootG_reg = matrix_return_copy(G, Nintco, Nintco); //G^1/2
  matrix_root(rootG_reg, Nintco, 0);

  if (Opt_params.print_lvl > 2) {
    psi::outfile->Printf( "\nrootG matrix:\n");
    print_matrix("outfile", rootG_reg, Nintco, Nintco);
  }

  double **rootG_inv = matrix_return_copy(G, Nintco, Nintco); //G^-1/2
  matrix_root(rootG_inv, Nintco, 1);

  if (Opt_params.print_lvl > 2) {
    psi::outfile->Printf( "G matrix:\n");
    print_matrix("outfile", G, Nintco, Nintco);
  }

  // Compute mass-weighted Hessian matrix:
  double **H  = p_Opt_data->g_H_pointer();     //internal coordinate Hessian
  double **H_m = init_matrix(Nintco, Nintco);
  double **T = init_matrix(Nintco, Nintco);
  // T = HG^1/2
  opt_matrix_mult(H, 0, rootG_reg, 0, T, 0, Nintco, Nintco, Nintco, 0);
  // H_m = (G^1/2)H(G^1/2) = (G^1/2)T
  opt_matrix_mult(rootG_reg, 0, T, 0, H_m, 0, Nintco, Nintco, Nintco, 0);
  free_matrix(T);

  if (Opt_params.print_lvl > 2) {
    psi::outfile->Printf("Mass-weighted Hessian:\n");
    print_matrix("outfile", H_m, Nintco, Nintco);
  }

  // Variables to calculate predicted energy change (DE_predicted)
  double DE_projected;
  double dq_n;    //step norm
  double *dq_u;   //unit vector in step direction
  double dq_g;    //gradient in step direction
  double dq_h;    //hessian in step direction

  //Pivot Point and Initial Working Geometry:
  //step along along normalized, mass-weighted v
  if(at_FS) {
cout << "First point of constrained optimization.\n";
    psi::outfile->Printf( "\tFirst point of constrained optimization.\n");
    // If starting from TS, follow lowest-eigenvalued eigenvector.
    // Otherwise, follow the gradient (negative of the force vector).
    double *v;
    if (at_TS) {
      //v = lowest_evector(H_m, Nintco);
      v = lowest_evector(H, Nintco);

      if(Opt_params.IRC_direction == OPT_PARAMS::FORWARD)
        psi::outfile->Printf( "\tStepping in forward direction from TS.\n");
      else if (Opt_params.IRC_direction == OPT_PARAMS::BACKWARD) {
        psi::outfile->Printf( "\tStepping in backward direction from TS.\n");
        array_scm(v, -1, Nintco);
      }
    }
    else {
      psi::outfile->Printf( "\tStepping along IRC using gradient.\n");
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

    psi::outfile->Printf( "\n\tVector to follow: \n");
    print_array("outfile", v, Nintco);

    if (Opt_params.print_lvl > 2) {
      psi::outfile->Printf( "\nDq to pivot point: \n");
      print_array("outfile", dq_pivot, Nintco);
      psi::outfile->Printf( "\nDq to next geometry: \n");
      print_array("outfile", dq, Nintco);
    }

    free_array(v);

/* // Check calculated step size.
double **G_inv = symm_matrix_inv(G, Nintco, Nintco);
double *G_inv_dq = init_array(Nintco);
opt_matrix_mult(G_inv, 0, &dq, 1, &G_inv_dq, 1, Nintco, Nintco, 1, 0);
double dq_norm = sqrt( array_dot(dq, G_inv_dq, Nintco) );
psi::outfile->Printf( "\nCheck dq_norm to first point (dq G^-1 dq^t)^1/2: %20.15f\n", dq_norm);
opt_matrix_mult(G_inv, 0, &dq_pivot, 1, &G_inv_dq, 1, Nintco, Nintco, 1, 0);
dq_norm = sqrt( array_dot(dq_pivot, G_inv_dq, Nintco) );
psi::outfile->Printf( "\nCheck dq_norm to pivot point (dq G^-1 dq^t)^1/2: %20.15f\n", dq_norm);
free_array(G_inv_dq);
free_matrix(G_inv);
// */

    double *q = intco_values();
    double *q_pivot = init_array(Nintco);
    for(int i=0; i<Nintco; i++)
      q_pivot[i] = q[i] + dq_pivot[i];
    free_array(dq_pivot);

    double *x = g_geom_array();
    double *f_x = g_grad_array();
    array_scm(f_x, -1, Ncart);
    double *f_q_copy = init_array(Nintco);
    array_copy(f_q, f_q_copy, Nintco);

    // q_pivot, q, x,f_q, f_x should NOT be freed; pointers are assigned in IRC_data
    p_irc_data->add_irc_point(p_irc_data->g_next_coord_step(), q_pivot, q, x, f_q_copy, f_x, g_energy(),
                              p_irc_data->step_length, p_irc_data->arc_length, p_irc_data->line_length);

    // Do displacements for each fragment separately.
    for (int f=0; f<fragments.size(); ++f) {
      if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
        psi::outfile->Printf("\tDisplacements for frozen fragment %d skipped.\n", f+1);
        continue;
      }
      fragments[f]->displace(&(dq[g_intco_offset(f)]), &(f_q[g_intco_offset(f)]), g_atom_offset(f));
    }
    // Do displacements for interfragment coordinates.
    for (int I=0; I<interfragments.size(); ++I) {
      if (interfragments[I]->is_frozen() || Opt_params.freeze_interfragment) {
        psi::outfile->Printf("\tDisplacements for frozen interfragment %d skipped.\n", I+1);
        continue;
      }
      interfragments[I]->orient_fragment( &(dq[g_interfragment_intco_offset(I)]),
                                         &(f_q[g_interfragment_intco_offset(I)]) );
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

    psi::outfile->Printf("\tGradient in step direction: %15.10lf\n", dq_g);
    psi::outfile->Printf("\tHessian in step direction : %15.10lf\n", dq_h);

    DE_projected = DE_nr_energy(dq_n, dq_g, dq_h);
    psi::outfile->Printf("\tProjected energy change for next step: %20.15lf\n", DE_projected);

    p_Opt_data->save_step_info(DE_projected, dq_u, dq_n, dq_g, dq_h);
    free_array(dq_u);

    // Increment to keep track of number of constrained optimization steps for this point.
    p_irc_data->sphere_step++;

    // g2D = g_geom_2D();
    // psi::outfile->Printf("Geometry after symmetrization\n");
    // print_matrix("outfile", g2D, Natom, 3);
    // free_matrix(g2D);
    symmetrize_geom();

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
      psi::outfile->Printf( "Eigenvector %i of mass-weighted Hessian:\n", i);
      print_array("outfile", V[i], Nintco);
    }
  }
  psi::outfile->Printf( "\n\tEigenvalues of the mass-weighted Hessian:\n");
  print_array("outfile", h, Nintco);

//2. Express p and g, in mass-weighted coordinates and in the eigenbasis of Hm, the mass-weighted Hessian.
  double *q = intco_values();
  double *q_pivot = p_irc_data->steps.back()->g_q_pivot();

  double *dq_0 = p_Opt_data->g_dq_pointer(p_Opt_data->nsteps() - 2);
  double *f_q_0 = p_Opt_data->g_forces_pointer(p_Opt_data->nsteps() - 2);

  //psi::outfile->Printf( "\nRetrieved pivot point from IRC_data: \n");
  //print_array("outfile", q_pivot, Nintco);

  double *p = init_array(Nintco);  //vector step from pivot point
  double *p_0 = init_array(Nintco);
  double *p_m = init_array(Nintco);//mass-weighted step
  double *p_m0 = init_array(Nintco);
  for(int i=0; i<Nintco; i++)
    p[i] = q[i] - q_pivot[i];
  for(int i=0; i<Nintco; i++)
    p_0[i] = p[i] - dq_0[i];
  for(int i=0; i<Nintco; i++)
    p_m[i] = array_dot(rootG_inv[i], p, Nintco);
  for(int i=0; i<Nintco; i++)
    p_m0[i] = array_dot(rootG_inv[i], p_0, Nintco);

  double *g = init_array(Nintco);
  double *g_0 = init_array(Nintco);
  double *g_m = init_array(Nintco);//mass-weighted gradient
  double *g_m0 = init_array(Nintco);
  for(int i=0; i<Nintco; i++)
    g[i] = -f_q[i];
  for(int i=0; i<Nintco; i++)
    g_0[i] = -f_q_0[i];
  for(int i=0; i<Nintco; i++)
    g_m[i] = array_dot(rootG_reg[i], g, Nintco);
  for(int i=0; i<Nintco; i++)
    g_m0[i] = array_dot(rootG_reg[i], g_0, Nintco);

psi::outfile->Printf( "\np_m before linear interpolation: ");
print_array("outfile", p_m, Nintco);
psi::outfile->Printf( "\ng_m before linear interpolation: ");
print_array("outfile", g_m, Nintco);

if(0 && p_irc_data->sphere_step > 1)
{
  if(0)
    GS_interpolation(p_m, p_m0, g_m, g_m0, s, Nintco);
  else
    interpolation(p_m, p_m0, g_m, g_m0, s, Nintco);
}

psi::outfile->Printf( "\np_m after linear interpolation:  ");
print_array("outfile", p_m, Nintco);
psi::outfile->Printf( "\ng_m after linear interpolation:  ");
print_array("outfile", g_m, Nintco);

  double *p_h = init_array(Nintco);//in the basis of H_m
  for(int i=0; i<Nintco; i++)
    p_h[i] = array_dot(p_m, V[i], Nintco);
  double *g_h = init_array(Nintco);//in the basis of H_m
  for(int i=0; i<Nintco; i++)
    g_h[i] = array_dot(g_m, V[i], Nintco);

  if (Opt_params.print_lvl > 2) {
    psi::outfile->Printf( "\np (q-q_pivot):\n");
    print_array("outfile", p, Nintco);
    psi::outfile->Printf( "\np_m:\n");
    print_array("outfile", p_m, Nintco);
    psi::outfile->Printf( "\np_h:\n");
    print_array("outfile", p_h, Nintco);
    psi::outfile->Printf( "\ng_m:\n");
    print_array("outfile", g_m, Nintco);
    psi::outfile->Printf( "\ng_h:\n");
    print_array("outfile", g_h, Nintco);
  }

//3. solve equation 26 in Gonzalez & Schlegel (1990) for lambda based on current p-vector
  static double lagrangian; // value of lagrangian expression
  double *df = init_array(5); //derivatives of lagrangian
  double h_f;

  static double lambda = 1.5 * h[0]; // Make 50% lower than the lowest (negative) eval

  int lag_iter=0;
  int iter=0;
  int coarse_iter=0;
  double temp_lambda;
  double old_lambda = -999;
  double old_lagrangian;
  double lb_lagrangian = -100;
  double ub_lagrangian = 100;
  double lb_lambda;
  double ub_lambda;
  int order = 4;

  old_lagrangian = lag_function(lambda, df, h, p_h, g_h, Nintco, s);
  lagrangian = lag_function(lambda, df, h, p_h, g_h, Nintco, s);
  while(lagrangian*old_lagrangian > 0 && coarse_iter < 1000)
  {
    if(lagrangian < 0 && fabs(lagrangian) < fabs(lb_lagrangian))
    {
      lb_lagrangian = lagrangian;
      lb_lambda = lambda;
    }

    if(lagrangian > 0 && fabs(lagrangian) < fabs(ub_lagrangian))
    {
      ub_lagrangian = lagrangian;
      ub_lambda = lambda;
    }

    coarse_iter++;
    old_lagrangian = lagrangian;
    lambda -= 1;
    lagrangian = lag_function(lambda, df, h, p_h, g_h, Nintco, s);
  }

  psi::outfile->Printf( "\n\tDetermining lagrangian multiplier for constrained minimization.\n");
  psi::outfile->Printf( "\t   Iter     Multiplier     Lagrangian\n");

  while (fabs(lambda - old_lambda) > 1e-16)
  {
    old_lagrangian = lagrangian;
    lagrangian = lag_function(lambda, df, h, p_h, g_h, Nintco, s);
    h_f = -df[0] / df[1];

    if(lagrangian < 0 && fabs(lagrangian) < fabs(lb_lagrangian))
    {
      lb_lagrangian = lagrangian;
      lb_lambda = lambda;
    }

    if(lagrangian > 0 && fabs(lagrangian) < fabs(ub_lagrangian))
    {
      ub_lagrangian = lagrangian;
      ub_lambda = lambda;
    }

    if(lagrangian*old_lagrangian < 0)
    {
      temp_lambda = lambda;
      lambda = (temp_lambda + old_lambda) / 2;
      old_lambda = temp_lambda;
    }
    else if(order == 1)
    {
      old_lambda = lambda;
      lambda += h_f;
    }
    else if(order == 2)
    {
      old_lambda = lambda;
      lambda += h_f / (1 + 0.5 * h_f * df[2] / df[1]);
    }
    else if(order == 3)
    {
      old_lambda = lambda;
      lambda += h_f * (1 + 0.5 * h_f * df[2] / df[1]) / ( 1 + (df[2]/df[1]) * h_f + (df[3]/(6*df[1])) * h_f * h_f );
    }
    else
    {
      old_lambda = lambda;
      lambda += h_f * (24*df[1] + 24*df[2] * h_f + 4*df[3] * h_f * h_f)
                    / (24*df[1] + 36*h_f * df[2] + 6*(df[2]*df[2]/df[1]) * h_f * h_f
                             + 8*df[3] * h_f * h_f + df[4] * h_f * h_f * h_f);
    }

    psi::outfile->Printf( "\t  %5d%15.3e%15.3e\n", lag_iter, lambda, lagrangian);

    ++lag_iter;

    if (++iter > 50)
    {
      old_lambda = lambda;
      lambda = (lb_lambda + ub_lambda) / 2;
    }
   
    if (lag_iter > 200)
      throw(INTCO_EXCEPT("Could not converge lagrangian for constrained minimization"));
  }


//4. Find dq_m from equation 24 in Gonzalez & Schlegel (1990)
// dq_m = (H_m - lambda I)^-1 [lambda * p_m - g_m]
  double *dq_m = init_array(Nintco);

  // Compute (H_m - lambda * I)^(-1)
  for(int i=0; i<Nintco; i++)
    H_m[i][i] -= lambda;
  double **H_m_lag_inv = symm_matrix_inv(H_m, Nintco, 1);

  double *g_lag_p = init_array(Nintco);
  for(int i=0; i<Nintco; i++)
    g_lag_p[i] = lambda * p_m[i] - g_m[i];

  opt_matrix_mult(H_m_lag_inv, 0, &g_lag_p, 1, &dq_m, 1, Nintco, Nintco, 1, 0);

  free_matrix(H_m_lag_inv);

//5. find dq ( = (G^1/2)dq_m) and do displacements
  opt_matrix_mult(rootG_reg, 0, &dq_m, 1, &dq, 1, Nintco, Nintco, 1, 0);

// Check equations 22 and 23 
/* double ss = 0.0;
  for (int i=0; i<Nintco; ++i)
    ss += (p_m[i] + dq_m[i]) * (p_m[i] + dq_m[i]) ;
  ss -= (0.5 * s) * (0.5 * s);
  psi::outfile->Printf("Eqn. 22 step size check %20.15lf\n", ss);

  double *lhs1 = init_array(Nintco);
  double *lhs2 = init_array(Nintco);
  for (int i=0; i<Nintco; ++i) {
    lhs1[i] = g_m[i] - lambda * p_m[i];
    for (int j=0; j<Nintco; ++j)
      lhs2[i] += H_m[i][j] * dq_m[j];
  }
  psi::outfile->Printf("Eqn. 23 check\n");
  for (int i=0; i<Nintco; ++i)
    psi::outfile->Printf( "%20.15lf\n", lhs1[i] + lhs2[i]);

  free_array(lhs1);
  free_array(lhs2); */

  // Do displacements for each fragment separately.
  for (int f=0; f<fragments.size(); ++f) {
    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
      psi::outfile->Printf("\tDisplacements for frozen fragment %d skipped.\n", f+1);
      continue;
    }
    fragments[f]->displace(&(dq[g_intco_offset(f)]), &(f_q[g_intco_offset(f)]), g_atom_offset(f));
  }
  // Do displacements for interfragment coordinates.
  for (int I=0; I<interfragments.size(); ++I) {
    if (interfragments[I]->is_frozen() || Opt_params.freeze_interfragment) {
      psi::outfile->Printf("\tDisplacements for frozen interfragment %d skipped.\n", I+1);
      continue;
    }
    interfragments[I]->orient_fragment( &(dq[g_interfragment_intco_offset(I)]),
                                       &(f_q[g_interfragment_intco_offset(I)]) );
  }

//find arclength of step and pass to irc_data to be stored on convergence
  double *new_q = intco_values();
  double *old_q = p_irc_data->g_q();
  double *new_dq = init_array(Nintco);
  for(int i=0; i<Nintco; i++)
  {
    p[i] = new_q[i] - q_pivot[i];
    new_dq[i] = new_q[i] - old_q[i];
  }
  for(int i=0; i<Nintco; i++)
  {
    p_m[i] = array_dot(rootG_inv[i], p, Nintco);
    dq_m[i] = array_dot(rootG_inv[i], new_dq, Nintco);
  }
  double *u_dq_m = init_array(Nintco);
  double *u_p_m = init_array(Nintco);
  array_copy(dq_m, u_dq_m, Nintco);
  array_copy(p_m, u_p_m, Nintco);
  array_normalize(u_dq_m, Nintco);
  array_normalize(u_p_m, Nintco);

//beta is the angle between the vector p_m, pointing radially from the hypersphere parallel to the
//gradient at the converged point, and the vector dq_m, connecting the starting point and the new
//point converged on the hypersphere
  double u_dqm_pm = array_dot(u_dq_m, u_p_m, Nintco);
  double beta = acos( fabs( u_dqm_pm ) );
  p_irc_data->arc_length = fabs( s * beta / tan( beta ) );
  p_irc_data->line_length = sqrt( array_dot(new_dq, new_dq, Nintco) );

  free_array(u_dq_m);
  free_array(u_p_m);
  free_array(g_lag_p);
  free_array(new_dq);

  free_array(dq_m);
  // Compute norm, unit vector, gradient and Hessian in the step direction
  dq_n = sqrt( array_dot(dq, dq, Nintco) );

  dq_u = init_array(Nintco);
  array_copy(dq, dq_u, Nintco);
  array_normalize(dq_u, Nintco);

  dq_g = -1 * array_dot(f_q, dq_u, Nintco);

  dq_h = 0;
  for(int i=0; i<Nintco; i++)
    dq_h += dq_u[i] * array_dot(H[i], dq_u, Nintco);

  psi::outfile->Printf("\tGradient in step direction: %15.10lf\n", dq_g);
  psi::outfile->Printf("\tHessian in step direction : %15.10lf\n", dq_h);

  DE_projected = DE_nr_energy(dq_n, dq_g, dq_h);
  psi::outfile->Printf("\tProjected energy change for next step: %20.15lf\n", DE_projected);

  p_Opt_data->save_step_info(DE_projected, dq_u, dq_n, dq_g, dq_h);
  free_array(dq_u);
  p_irc_data->sphere_step++;

  //double **g2D = g_geom_2D();
  //psi::outfile->Printf("Geometry before symmetrization\n");
  //print_matrix("outfile", g2D, Natom, 3);
  //free_matrix(g2D);
  symmetrize_geom();

  free_matrix(rootG_reg);
  free_matrix(rootG_inv);
  free_matrix(H_m);
  free_matrix(V);
  free_array(h);
  free_array(p);
  free_array(p_0);
  free_array(p_m);
  free_array(p_m0);
  free_array(p_h);
  free_array(g_m);
  free_array(g_m0);
  free_array(g_h);
  free_array(df);
}

// Compute step distance normalization factor N= (gGg^t)^-1/2
double step_N_factor(double **G, double *g, int Nintco) {
  double N = 0;

  for(int i=0; i<Nintco; i++) {
    double Ggsum = 0;
    for(int j = i+1; j<Nintco; j++)
      Ggsum += G[i][j] * g[j];
    N += g[i] * (G[i][i] * g[i] + 2 * Ggsum);
  }
  N = 1 / sqrt(N);
  if (Opt_params.print_lvl > 2)
    psi::outfile->Printf("\tNormalizing factor N: %15.10lf\n", N);

  return N;
}

double *lowest_evector(double **H, int Nintco)
{
  double **V = matrix_return_copy(H, Nintco, Nintco);
  double *V_evals = init_array(Nintco);
  opt_symm_matrix_eig(V, Nintco, V_evals);

  // find largest element
  double max_element = -1;
  for(int i=0; i<Nintco; i++)
    if( fabs(V[0][i]) > fabs(max_element) )
      max_element = V[0][i];

  int sign; 
  (max_element == fabs(max_element)) ? sign = 1 : sign = -1;

  double *min_evector = init_array(Nintco);
  for(int i=0; i<Nintco; i++)
    min_evector[i] = sign * V[0][i];

  cout << "making largest element of eigenvector positive\n";

  free_matrix(V);
  free_array(V_evals);

  return min_evector;
}

double lag_function(double l, double *f, double *h, double *p, double *g, int dim, double s)
{
  double tval;
  double D;

  for(int i=0; i<5; i++)
    f[i] = 0;

  for(int i=0; i<dim; i++)
  {
    D = h[i] - l;
    tval = (h[i]*p[i] - g[i]) * (h[i]*p[i] - g[i]) / (D*D);

    f[0] +=       tval;
    f[1] += 2   * tval / (D);
    f[2] += 6   * tval / (D*D);
    f[3] += 24  * tval / (D*D*D);
    f[4] += 120 * tval / (D*D*D*D);
  }
  f[0] -= (s/2) * (s/2);

  return f[0];
}

void GS_interpolation(double *p, double *p_0, double *g, double *g_0, double s, int dim)
{
  double p_p = array_dot(p, p, dim);
  double p0_p0 = array_dot(p_0, p_0, dim);
  double cosTh = array_dot(p, p_0, dim) / sqrt(p_p * p0_p0);
  double Th = acos( cosTh );
  double g_p = array_dot(g, p, dim);
  double g0_p0 = array_dot(g_0, p_0, dim);

  double proj;
  double proj_0;
  double gPer = 0;
  double gPer_0 = 0;
  for(int i=0; i<dim; i++)
  {
    proj = g[i] - (g_p/p_p) * p[i];
    proj_0 = g_0[i] - (g0_p0/p0_p0) * p_0[i];

    gPer += proj * proj;
    gPer_0 += proj_0 * proj_0;
  }

  gPer = sqrt(gPer);
  gPer_0 = sqrt(gPer_0);

  double Th_i = (gPer_0 * Th) / (gPer_0 - gPer);

  double a = sin(Th_i) / sin(Th);
  double b = cos(Th_i) - (sin(Th_i)*cosTh/sin(Th));

  for(int i=0; i<dim; i++)
  {
    p[i] = a*p[i] + b*p_0[i];
    g[i] = (Th_i/Th)*g[i] + (1 - Th_i/Th)*g_0[i];
  }
}

void interpolation(double *p, double *p_0, double *g, double *g_0, double s, int dim)
{
  double p_p = array_dot(p, p, dim);
  double p0_p0 = array_dot(p_0, p_0, dim);
psi::outfile->Printf( "p_norm is %f\n", sqrt(p_p));
psi::outfile->Printf( "p0_norm is %f\n", sqrt(p0_p0));

  double cosTh = array_dot(p, p_0, dim) / sqrt(p_p * p0_p0);
  double Th = acos( cosTh );

print_array("outfile", p, dim);
print_array("outfile", p_0, dim);
print_array("outfile", g, dim);
print_array("outfile", g_0, dim);
psi::outfile->Printf( "step size is %f\n", s);
psi::outfile->Printf( "p_p is %f\n", p_p);
psi::outfile->Printf( "cosTh is %f\n", cosTh);
psi::outfile->Printf( "Th is %f\n", Th);

  double g_p = array_dot(g, p, dim);
  double g0_p0 = array_dot(g_0, p_0, dim);
psi::outfile->Printf( "g_p is %f\n", g_p);
psi::outfile->Printf( "g0_p0 is %f\n", g0_p0);

  double proj;
  double proj_0;
  double gPer = 0;
  double gPer_0 = 0;
  for(int i=0; i<dim; i++)
  {
    //projection of g perpendicular to p found by subtracting from g
    //the projection of g onto p
    proj = g[i] - (g_p/p_p) * p[i];
    proj_0 = g_0[i] - (g0_p0/p0_p0) * p_0[i];

    //we are looking for the norm of the component of g perpendicular to g
    //hence, we dot the elements (e0*e0 + ... + edim*edim)
    gPer += proj * proj;
    gPer_0 += proj_0 * proj_0;
  }
  //having dotted, we take the square root to get the norm
  gPer = sqrt(gPer);
  gPer_0 = sqrt(gPer_0);

psi::outfile->Printf( "gPer is %f\n", gPer);
psi::outfile->Printf( "gPer_0 is %f\n", gPer_0);

  double Th_i = (Th * gPer) / (gPer - gPer_0);
  double cosTh_i = cos( Th_i );
  double sinTh_i = sin( Th_i );
psi::outfile->Printf( "Th_i is %f\n", Th_i);

  double *pPer = init_array(dim);
  for(int i=0; i<dim; i++)
    pPer[i] = p_0[i] - cosTh * p[i];
  array_normalize(pPer, dim);
  array_scm(pPer, sqrt(p0_p0), dim);
print_array("outfile", pPer, dim);

  for(int i=0; i<dim; i++)
  {
    p[i] = cosTh_i*p[i] - sinTh_i*pPer[i];
    g[i] = (1 + Th_i/Th)*g[i] - (Th_i/Th)*g_0[i];
  }

  free_array(pPer);
}

}

