/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

/*! \file molecule_irc_step.cc
    \ingroup optking
    \Gonzalez Schlegel (1990) Reaction Path Following Algorithm
*/

#include "molecule.h"

#include <iostream>
#include <sstream>

#include "linear_algebra.h"
#include "atom_data.h"
#include "psi4/optking/physconst.h"
#include "IRC_data.h"

#include "print.h"
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

//void GS_interpolation(double *p, double *p_0, double *g, double *g_0, double s, int dim);
void GS_interpolation(double *p, double *p_0, double *g, double *g_0, int dim);
//void interpolation(double *p, double *p_0, double *g, double *g_0, double s, int dim);
void interpolation(double *p, double *p_0, double *g, double *g_0, int dim);

double lag_function(double l, double *f, double *h, double *p, double *g, int dim, double s);
// compute change in energy according to quadratic approximation
inline double DE_nr_energy(double step, double grad, double hess)
{
  return (step * grad + 0.5 * step * step * hess);
}

void IRC_DATA::point_converged(opt::MOLECULE &mol)
{

  // Why do the dot below?
  //if(steps.size() > 1)
    //double f_dot = array_dot(steps[steps.size()-1]->g_f_q(), steps[steps.size()-2]->g_f_q(), mol.Ncoord());

  sphere_step = 0;


  mol.irc_step();
}

void IRC_DATA::progress_report(opt::MOLECULE &mol)
{
  double DE;
  int dim = mol.Ncoord();
  int blocks = 4;
  int sign = 1;

  if(Opt_params.IRC_direction == OPT_PARAMS::BACKWARD)
    sign = -1;

//Printing Energies and Energy Changes for Each Step
  oprintf_out("@IRC ----------------------------------------------\n");
  oprintf_out("@IRC            ****      IRC Report      ****\n");
  oprintf_out("@IRC ----------------------------------------------\n");
  oprintf_out("@IRC  Step    Energy              Change in Energy \n");
  oprintf_out("@IRC ----------------------------------------------\n");
  for (std::size_t i=0; i<steps.size(); ++i)
  {
    if (i == 0) DE = g_step(i).g_energy();
    else DE = g_step(i).g_energy() - g_step(i-1).g_energy();

    oprintf_out("@IRC  %3d %18.12lf  %18.12lf\n", i, g_step(i).g_energy(), DE);
  }
  oprintf_out("@IRC ----------------------------------------------\n\n");

//Printing Internal Coordinates for Each step
  oprintf_out("@IRC -----------------------------------------------------\n");
  oprintf_out("@IRC              ****     IRC Steps     ****             \n");
  oprintf_out("@IRC -----------------------------------------------------");
  for(int j=0; j < dim/blocks; j++)
  {
    oprintf_out("\n@IRC        |          Distance         |\n");
    oprintf_out(  "@IRC Step   | Step    Arc       Line    |");
    for(int i = (j*blocks); i < ((j+1)* blocks); i++)
    {
      oprintf_out("    Coord %3d", i);
    }
    oprintf_out("\n");
    oprintf_out("@IRC --------------------------------------");
    for(int i = (j*blocks); i < ((j+1)* blocks); i++)
    {
      oprintf_out("-------------");
    }
    oprintf_out("\n");
    for (std::size_t i=0; i<steps.size(); ++i)
    {
      oprintf_out("@IRC  %3d %9.2lf %9.5lf  %9.5lf   ", i, sign*g_step(i).g_step_dist(), sign*g_step(i).g_arc_dist(), sign*g_step(i).g_line_dist());
      for(int k = (j*blocks); k < ((j+1)*blocks); k++)
        oprintf_out("%13.8f",g_step(i).g_q()[k]);
      oprintf_out("\n");
    }
    oprintf_out("@IRC --------------------------------------");
    for(int i = (j*blocks); i < ((j+1)* blocks); i++)
    {
      oprintf_out("-------------");
    }
  }
  if(dim % blocks != 0)
  {
    oprintf_out("\n@IRC         |          Distance         |\n");
    oprintf_out(  "@IRC  Step   | Step    Arc       Line    |");
    for(int i = (dim - (dim % blocks)); i < dim; i++)
    {
      oprintf_out("    Coord %3d", i);
    }
    oprintf_out("\n");
    oprintf_out("@IRC --------------------------------------");
    for(int i = (dim - (dim % blocks)); i < dim; i++)
    {
      oprintf_out("-------------");
    }
    oprintf_out("\n");
    for (std::size_t i=0; i<steps.size(); ++i)
    {
      oprintf_out("@IRC  %3d %9.2lf %9.5lf  %9.5lf   ", i, sign*g_step(i).g_step_dist(), sign*g_step(i).g_arc_dist(), sign*g_step(i).g_line_dist());
      for(int k = (dim - (dim % blocks)); k < dim; k++)
        oprintf_out("%13.8f",g_step(i).g_q()[k]);
      oprintf_out("\n");
    }
    oprintf_out("@IRC --------------------------------------");
    for(int i = (dim - (dim % blocks)); i < dim; i++)
    {
      oprintf_out("-------------");
    }
  }

  oprintf_out("\n");
  oprintf_out("\n");

  mol.print_coords(psi_outfile, qc_outfile);
  mol.print_simples(psi_outfile, qc_outfile);
}


// See Gonzalez and Schlegel JCP 2154 (1989).
void MOLECULE::irc_step(void)
{
  // Are we at the TS?  at_TS
  bool at_TS = !(p_irc_data->size());
  if (at_TS) oprintf_out("\n  IRC_DATA is empty, so we are at the transition state.\n");

  // Is this one the first step toward a new path point?  at_FS
  bool at_FS = !(p_irc_data->sphere_step);

  double s    = Opt_params.IRC_step_size;      //step size
  p_irc_data->step_length = s;
  double *dq  = p_Opt_data->g_dq_pointer();    //internal coordinate change
  double *f_q = p_Opt_data->g_forces_pointer();//internal coordinate gradient
  int Nintco = Ncoord();
  int Natom = g_natom();
  int Ncart = 3 * Natom;
  bool answer = 1;


  int opt_iter = p_Opt_data->g_iteration() - 1;
  if(opt_iter > 2 && p_irc_data->in_min_range) {
    if(    fabs(p_Opt_data->g_energy(opt_iter)     - p_Opt_data->g_energy(opt_iter - 2)) < 0.01e-03
        && fabs(p_Opt_data->g_energy(opt_iter - 1) - p_Opt_data->g_energy(opt_iter - 3)) < 0.01e-03 ) {
      p_irc_data->go = 0;
    }
  }

    if (p_irc_data->sphere_step == 1) {

      double *u_f_q = init_array(Nintco);
      double *u_f_q_0 = init_array(Nintco);
      array_copy(f_q, u_f_q, Nintco);
      array_copy(p_irc_data->g_f_q(), u_f_q_0, Nintco);
      array_normalize(u_f_q, Nintco);
      array_normalize(u_f_q_0, Nintco);
      double u_f_q_dot = array_dot(u_f_q, u_f_q_0, Nintco);

      if(u_f_q_dot < -0.5) {
        p_irc_data->in_min_range = 1;
      }

      if(p_irc_data->size() > 3) {
        if (u_f_q_dot < -0.7 ||
            (
              fabs(p_irc_data->g_line_dist(p_irc_data->size()-1) - p_irc_data->g_line_dist(p_irc_data->size()-2)) < s*10e-03
            )
           ) {
          oprintf_out("\n@IRC\n@IRC Houston, we've found a minimum!\n@IRC\n");

          if(Opt_params.IRC_stop == OPT_PARAMS::ASK) {
            cout << "Would you like to proceed? (1=yes, 0=no):";
            cin >> answer;
          }
          if(Opt_params.IRC_stop == OPT_PARAMS::STOP || !answer) {
            p_irc_data->go = 0;
          }
        }
      }

      free_array(u_f_q);
      free_array(u_f_q_0);
    }

  //The G matrix, its square root, and its inverse square root
  double **G         = compute_G(1);                          //G = BUB^t

  double **rootG_reg = matrix_return_copy(G, Nintco, Nintco); //G^1/2
  matrix_root(rootG_reg, Nintco, 0);

  if (Opt_params.print_lvl > 2) {
    oprintf_out( "\n@IRC rootG matrix:\n");
    oprint_matrix_out(rootG_reg, Nintco, Nintco);
  }

  double **rootG_inv = matrix_return_copy(G, Nintco, Nintco); //G^-1/2
  matrix_root(rootG_inv, Nintco, 1);

  if (Opt_params.print_lvl > 2) {
    oprintf_out( "@IRC G matrix:\n");
    oprint_matrix_out(G, Nintco, Nintco);
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
    oprintf_out("  Mass-weighted Hessian:\n");
    oprint_matrix_out(H_m, Nintco, Nintco);
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
    // If starting from TS, follow lowest-eigenvalued eigenvector.
    // Otherwise, follow the gradient (negative of the force vector).
    double *v;
    if (at_TS) {
      //v = lowest_evector(H_m, Nintco);
      v = lowest_evector(H, Nintco);

      if(Opt_params.IRC_direction == OPT_PARAMS::FORWARD)
        oprintf_out( "  Stepping in forward direction from TS.\n");
      else if (Opt_params.IRC_direction == OPT_PARAMS::BACKWARD) {
        oprintf_out( "  Stepping in backward direction from TS.\n");
        array_scm(v, -1, Nintco);
      }
    }
    else {
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

    if (Opt_params.print_lvl > 2) {
      oprintf_out( "\n@IRC Dq to pivot point: \n");
      oprint_array_out(dq_pivot, Nintco);
      oprintf_out( "\n@IRC Dq to next geometry: \n");
      oprint_array_out(dq, Nintco);
    }

    free_array(v);

    double *x = g_geom_array();
    double *f_x = g_grad_array();
    array_scm(f_x, -1, Ncart);
    double *f_q_copy = init_array(Nintco);
    array_copy(f_q, f_q_copy, Nintco);

    // RAK 11-14 Try to calculate xyz coordinates of the pivot point
    for (std::size_t f=0; f<fragments.size(); ++f) {
      if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
        oprintf_out("    Displacements for frozen fragment %d skipped.\n", f+1);
        continue;
      }
      fragments[f]->displace(&(dq_pivot[g_coord_offset(f)]), &(f_q[g_coord_offset(f)]), g_atom_offset(f));
    }
    for (std::size_t I=0; I<interfragments.size(); ++I) {
      if (interfragments[I]->is_frozen() || Opt_params.freeze_interfragment) {
        oprintf_out("    Displacements for frozen interfragment %d skipped.\n", I+1);
        continue;
      }
      interfragments[I]->orient_fragment( &(dq_pivot[g_interfragment_coord_offset(I)]), &(f_q[g_interfragment_coord_offset(I)]) );
    }
    double *x_pivot = g_geom_array();

    // new we can calculate q and q_pivot consistently (wrt the 180 discontinuity)
    set_geom_array(x);
    fix_tors_near_180();
    double *q = coord_values();
    set_geom_array(x_pivot);
    double *q_pivot = coord_values();

    // q_pivot, q, x,f_q, f_x should NOT be freed; pointers are assigned in IRC_data
    //p_irc_data->add_irc_point(p_irc_data->g_next_coord_step(), q_pivot, x_pivot, q, x, f_q_copy, f_x, g_energy(),
    //                          p_irc_data->step_length, p_irc_data->arc_length, p_irc_data->line_length);
    p_irc_data->add_irc_point(p_irc_data->g_next_coord_step(), q_pivot, x_pivot, q, x, f_q_copy, f_x, g_energy());

    /* START PRINT CONVERGED POINT INFO */
    int point_number = Opt_params.IRC_direction == OPT_PARAMS::FORWARD ? p_irc_data->size() - 1 : - (p_irc_data->size() - 1);
    set_geom_array(x);
    oprintf_out("\n@IRC\n");
    if(p_irc_data->size() == 1) {
      oprintf_out("@IRC  **** Point %2d on IRC path ****\n", point_number);
    }
    else {
      oprintf_out("@IRC  **** Point %2d on IRC path is optimized ****\n", point_number);
    }
    oprintf_out("@IRC    Final energy:         %20.13lf\n", p_Opt_data->g_energy());
    oprintf_out("@IRC    Arc path distance:    %20.13lf\n", p_irc_data->g_arc_dist());
    oprintf_out("@IRC    Linear path distance: %20.13lf\n", p_irc_data->g_line_dist());
    oprintf_out("@IRC\n");
    print_geom_out_irc();
    if (Opt_params.print_trajectory_xyz_file)
      print_xyz_irc(point_number, Opt_params.IRC_direction == OPT_PARAMS::FORWARD);
    oprintf_out("@IRC\n");
    oprintf_out("@IRC\n\n");
    set_geom_array(x_pivot);
    /* END PRINT CONVERGED POINT INFO */

    // Take step of dq_pivot length again since dq is 2*dq_pivot
    for (std::size_t f=0; f<fragments.size(); ++f) {
      if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
        oprintf_out("    Displacements for frozen fragment %d skipped.\n", f+1);
        continue;
      }
      fragments[f]->displace(&(dq_pivot[g_coord_offset(f)]), &(f_q[g_coord_offset(f)]), g_atom_offset(f));
    }
    // Do displacements for interfragment coordinates.
    for (std::size_t I=0; I<interfragments.size(); ++I) {
      if (interfragments[I]->is_frozen() || Opt_params.freeze_interfragment) {
        oprintf_out("    Displacements for frozen interfragment %d skipped.\n", I+1);
        continue;
      }
      interfragments[I]->orient_fragment( &(dq_pivot[g_interfragment_coord_offset(I)]),
                                         &(f_q[g_interfragment_coord_offset(I)]) );
    }
    // end RAK 11-14
    free_array(dq_pivot);

    // Compute norm, unit vector, gradient and Hessian in the step direction
    // Save results to Opt_data.
    dq_n = sqrt( array_dot(dq, dq, Nintco) );

    oprintf_out("     Norm of target step-size %10.5lf\n", dq_n);

    dq_u = init_array(Nintco);
    array_copy(dq, dq_u, Nintco);
    array_normalize(dq_u, Nintco);

    dq_g = -1 * array_dot(f_q, dq_u, Nintco);

    dq_h = 0;
    for(int i=0; i<Nintco; i++)
      dq_h += dq_u[i] * array_dot(H[i], dq_u, Nintco);

    DE_projected = DE_nr_energy(dq_n, dq_g, dq_h);

    p_Opt_data->save_step_info(DE_projected, dq_u, dq_n, dq_g, dq_h);
    free_array(dq_u);

    // Increment to keep track of number of constrained optimization steps for this point.
    p_irc_data->sphere_step++;

    // g2D = g_geom_2D();
    // oprintf_out("Geometry after symmetrization\n");
    // oprint_matrix_out(g2D, Natom, 3);
    // free_matrix(g2D);
    symmetrize_geom();

    return;
  }  // end first/pivot point

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
      oprintf_out( "  Eigenvector %i of mass-weighted Hessian:\n", i);
      oprint_array_out(V[i], Nintco);
    }
  }

//2. Express p and g, in mass-weighted coordinates and in the eigenbasis of Hm, the mass-weighted Hessian.

  /* RAK 11-14 Recompute pivot point internal coordinates; this is so that the
     torsional convention can be applied and it makes sense to subtract differences. */
  //double *q_pivot = p_irc_data->steps.back()->g_q_pivot();
  //double *q = coord_values();
  // save xyz, note configuration of torsions, compute q's, put xyz back.
  double *x_orig = g_geom_array();
  fix_tors_near_180();
  double *q = coord_values();
  double *x_pivot = p_irc_data->g_x_pivot();
  set_geom_array(x_pivot);
  double *q_pivot = coord_values();
  set_geom_array(x_orig);
  free_array(x_orig);
  // end RAK 11-14

  double *dq_0 = p_Opt_data->g_dq_pointer(p_Opt_data->nsteps() - 2);
  double *f_q_0 = p_Opt_data->g_forces_pointer(p_Opt_data->nsteps() - 2);

  //oprintf_out( "\nRetrieved pivot point from IRC_data: \n");
  //oprint_array_out(q_pivot, Nintco);

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

  if(0 && p_irc_data->sphere_step > 1) {
    if(0)
      //GS_interpolation(p_m, p_m0, g_m, g_m0, s, Nintco);
      GS_interpolation(p_m, p_m0, g_m, g_m0, Nintco);
    else
      //interpolation(p_m, p_m0, g_m, g_m0, s, Nintco);
      interpolation(p_m, p_m0, g_m, g_m0, Nintco);
  }

  double *p_h = init_array(Nintco);//in the basis of H_m
  for(int i=0; i<Nintco; i++)
    p_h[i] = array_dot(p_m, V[i], Nintco);
  double *g_h = init_array(Nintco);//in the basis of H_m
  for(int i=0; i<Nintco; i++)
    g_h[i] = array_dot(g_m, V[i], Nintco);

  if (Opt_params.print_lvl > 2) {
    oprintf_out( "\np (q-q_pivot):\n");
    oprint_array_out(p, Nintco);
    oprintf_out( "\np_m:\n");
    oprint_array_out(p_m, Nintco);
    oprintf_out( "\np_h:\n");
    oprint_array_out(p_h, Nintco);
    oprintf_out( "\ng_m:\n");
    oprint_array_out(g_m, Nintco);
    oprintf_out( "\ng_h:\n");
    oprint_array_out(g_h, Nintco);
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

  oprintf_out( "\n    Determining lagrangian multiplier for constrained minimization.\n");

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

    ++lag_iter;

    if (++iter > 50)
    {
      old_lambda = lambda;
      lambda = (lb_lambda + ub_lambda) / 2;
    }

    if (lag_iter > 200)
      throw(INTCO_EXCEPT("Could not converge lagrangian for constrained minimization"));
  }

  oprintf_out( "\n    Lagrangian multiplier is converged.\n");


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


  // Do displacements for each fragment separately.
  for (std::size_t f=0; f<fragments.size(); ++f) {
    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
      oprintf_out("  Displacements for frozen fragment %d skipped.\n", f+1);
      continue;
    }
    fragments[f]->displace(&(dq[g_coord_offset(f)]), &(f_q[g_coord_offset(f)]), g_atom_offset(f));
  }
  // Do displacements for interfragment coordinates.
  for (std::size_t I=0; I<interfragments.size(); ++I) {
    if (interfragments[I]->is_frozen() || Opt_params.freeze_interfragment) {
      oprintf_out("  Displacements for frozen interfragment %d skipped.\n", I+1);
      continue;
    }
    interfragments[I]->orient_fragment( &(dq[g_interfragment_coord_offset(I)]),
                                       &(f_q[g_interfragment_coord_offset(I)]) );
  }

//find arclength of step and pass to irc_data to be stored on convergence
  // RAK 11-14
  //double *new_q = coord_values();
  //double *old_q = p_irc_data->g_q();
  // save xyz, note configuration of torsions, compute q's, put old xyz back
  x_orig = g_geom_array();
  fix_tors_near_180();
  double *new_q = coord_values();
  double *old_x = p_irc_data->g_x();
  set_geom_array(old_x);
  double *old_q = coord_values();
  set_geom_array(x_orig);
  free_array(x_orig);
  // end RAK 11-14

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

  oprintf_out("  Gradient in step direction: %15.10lf\n", dq_g);
  oprintf_out("  Hessian in step direction : %15.10lf\n", dq_h);

  DE_projected = DE_nr_energy(dq_n, dq_g, dq_h);
  oprintf_out("  Projected energy change for next step: %20.15lf\n", DE_projected);

  p_Opt_data->save_step_info(DE_projected, dq_u, dq_n, dq_g, dq_h);
  free_array(dq_u);

  /* START PRINTING */
  {
    // for IRC only consider forces tangent to the hypersphere search surface
    double **Ginv = symm_matrix_inv(G, Nintco, Nintco);

    // gradient perpendicular to p and tangent to hypersphere is:
    // g_m' = g_m - (g_m^t p_m / p_m^t p_m) p_m, or
    // g'   = g   - (g^t p / (p^t G^-1 p)) G^-1 p
    double *Ginv_p = init_array(Nintco);
    for(int i=0; i<Nintco; i++)
      Ginv_p[i] = array_dot(Ginv[i], p, Nintco);
    free_matrix(Ginv);

    double overlap = array_dot(f_q, p, Nintco) / array_dot(p, Ginv_p, Nintco);

    double* f_tan = init_array(Nintco);
    array_copy(f_q, f_tan, Nintco);

    for (int i=0; i<Nintco; ++i)
      f_tan[i] -= overlap * Ginv_p[i];
    free_array(Ginv_p);

    double E  = p_Opt_data->g_energy();
    double E0 = p_Opt_data->g_last_energy();
    double DE = (p_Opt_data->g_iteration() > 1) ? E - E0 : E;
    double max_force = array_abs_max(f_tan, Nintco);
    double rms_force = array_rms(f_tan, Nintco);
    double max_disp = array_abs_max(dq, Nintco);
    double rms_disp = array_rms(dq, Nintco);
    //double p_Opt_data->g_forces_pointer()

    int point = at_FS ? p_irc_data->size() + 1 : p_irc_data->size();
    if (p_irc_data->go) {
      if (p_irc_data->sphere_step == 1) {
        oprintf_out("\n@IRC");
        oprintf_out("\n@IRC Point   Sphere Step       Energy          DE         MAX Force     RMS Force      MAX Disp      RMS Disp   \n");
        oprintf_out(  "@IRC -----------------------------------------------------------------------------------------------------------\n");
      }
      else {
        oprintf_out("\n     Point   Sphere Step       Energy          DE         MAX Force     RMS Force      MAX Disp      RMS Disp   \n");
        oprintf_out(  "     -----------------------------------------------------------------------------------------------------------\n");
      }
      oprintf_out("@IRC %3d   %8d      %16.8f %10.2e %1s  %10.2e %1s  %10.2e %1s  %10.2e %1s  %10.2e %1s  ~\n",
          point, p_irc_data->sphere_step, E,
          DE, (Opt_params.i_max_DE ? ((fabs(DE) < Opt_params.conv_max_DE) ? "*" : "") : "o"),
          max_force, (Opt_params.i_max_force ? ((fabs(max_force) < Opt_params.conv_max_force) ? "*" : "") : "o"),
          rms_force, (Opt_params.i_rms_force ? ((fabs(rms_force) < Opt_params.conv_rms_force) ? "*" : "") : "o"),
          max_disp, (Opt_params.i_max_disp ? ((fabs(max_disp) < Opt_params.conv_max_disp) ? "*" : "") : "o"),
          rms_disp, (Opt_params.i_rms_disp ? ((fabs(rms_disp) < Opt_params.conv_rms_disp) ? "*" : "") : "o"));
      oprintf_out("     -----------------------------------------------------------------------------------------------------------\n\n");
    }
  }
  /* END PRINTING */



  p_irc_data->sphere_step++;

  //double **g2D = g_geom_2D();
  //oprintf_out("Geometry before symmetrization\n");
  //oprint_matrix_out(g2D, Natom, 3);
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

//void GS_interpolation(double *p, double *p_0, double *g, double *g_0, double s, int dim)
void GS_interpolation(double *p, double *p_0, double *g, double *g_0, int dim)
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

//void interpolation(double *p, double *p_0, double *g, double *g_0, double s, int dim)
void interpolation(double *p, double *p_0, double *g, double *g_0, int dim)
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

  double Th_i = (Th * gPer) / (gPer - gPer_0);
  double cosTh_i = cos( Th_i );
  double sinTh_i = sin( Th_i );

  double *pPer = init_array(dim);
  for(int i=0; i<dim; i++)
    pPer[i] = p_0[i] - cosTh * p[i];
  array_normalize(pPer, dim);
  array_scm(pPer, sqrt(p0_p0), dim);

  for(int i=0; i<dim; i++)
  {
    p[i] = cosTh_i*p[i] - sinTh_i*pPer[i];
    g[i] = (1 + Th_i/Th)*g[i] - (Th_i/Th)*g_0[i];
  }

  free_array(pPer);
}

}
