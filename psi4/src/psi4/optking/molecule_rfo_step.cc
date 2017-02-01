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

/*! \file molecule.cc
    \ingroup optking
    \brief molecule class (really, molecular system class)
*/

#include "molecule.h"

#include <iostream>
#include <sstream>

#include "linear_algebra.h"
#include "atom_data.h"
#include "psi4/optking/physconst.h"

#include "print.h"
#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
#endif

namespace opt {

// compute change in energy according to RFO approximation
inline double DE_rfo_energy(double rfo_t, double rfo_g, double rfo_h) {
  return (rfo_t * rfo_g + 0.5 * rfo_t * rfo_t * rfo_h)/(1 + rfo_t*rfo_t);
}

// Take Rational Function Optimization step
void MOLECULE::rfo_step(void) {
  int i, j;
  int dim = Ncoord();
  double tval, tval2;
  double *fq = p_Opt_data->g_forces_pointer();
  double **H = p_Opt_data->g_H_pointer();
  double *dq = p_Opt_data->g_dq_pointer();
  const int max_projected_rfo_iter = 25;

  oprintf_out("\tTaking RFO optimization step.\n");

  // Determine the eigenvectors/eigenvalues of H.  Used in RS-RFO
  double **Hevects = matrix_return_copy(H,dim,dim);
  double *h = init_array(dim);
  opt_symm_matrix_eig(Hevects, dim, h);

  // Build the original, unscaled RFO matrix.
  // For restricted-step RFO, the matrix is not symmetric, so we must copy the whole matrix.
  double **RFO = init_matrix(dim+1,dim+1);
  for (i=0; i<dim; ++i)
    for (j=0; j<dim; ++j)
      RFO[i][j] = H[i][j];

  for (i=0; i<dim; ++i)
    RFO[i][dim]= RFO[dim][i] = -fq[i];

  if (Opt_params.print_lvl >= 4) {
    oprintf_out("Original, unscaled RFO mat\n");
    oprint_matrix_out(RFO, dim+1, dim+1);
  }

  int rfo_root;     // root to follow
  double rfo_g;         // gradient in step direction
  double rfo_h;         // hessian in step direction
  double DE_projected; // projected energy change by quadratic approximation
  bool symm_rfo_step = false;
  double * rfo_old_evect = p_Opt_data->g_rfo_eigenvector_pointer();
  double lambda, sum;
  double analyticDerivative;  // d(norm step squared) / d(alpha)
  double trust = Opt_params.intrafragment_step_limit;
  double dqtdq = 10;     // square of norm of step
  double alpha = 1;      // scaling factor for RS-RFO, scaling matrix is sI

  double **SRFO = init_matrix(dim+1,dim+1);
  double *SRFOevals = init_array(dim+1);
  double *rfo_u = init_array(dim); // unit vector in step direction
  bool converged = false;

  bool rfo_follow_root = Opt_params.rfo_follow_root;
  double *last_iter_evect = init_array(dim);
  if (p_Opt_data->g_iteration() > 1 && rfo_follow_root)
    array_copy(rfo_old_evect, last_iter_evect, dim); // start with vector from previous iter

  //Iterative sequence to find alpha; we'll give it max_projected_rfo_iter tries
  int iter = -1;
  while (!converged && iter<max_projected_rfo_iter) {

    ++iter;

    // If we get to iteration 14, and we haven't not yet converged, than
    // bail out on the restricted-step algorithm.  Set alpha=1 and apply the crude
    // intrafragment_step_limit below.

    if (iter == max_projected_rfo_iter) {
      oprintf_out("\tFailed to converge alpha.  Doing simple step scaling instead.\n");

      alpha = 1;
    }
    else if (Opt_params.simple_step_scaling) // If simple_step_scaling is on, not an iterative method.
      iter = max_projected_rfo_iter;

    // Scale the RFO matrix.
    for (i=0; i<=dim; i++) {
      for (j=0; j<dim; j++)
        SRFO[j][i] = RFO[j][i] / alpha;
      SRFO[dim][i] = RFO[dim][i];
    }
    if (Opt_params.print_lvl >= 4) {
      oprintf_out("\nScaled RFO matrix.\n");
      oprint_matrix_out( SRFO, dim+1, dim+1);

    }

    //Find the eigenvectors and eigenvalues of RFO matrix.
    opt_asymm_matrix_eig(SRFO, dim+1, SRFOevals);
    if (Opt_params.print_lvl >= 4) {
      oprintf_out("Eigenvectors of scaled RFO matrix.\n");
      oprint_matrix_out(SRFO, dim+1,dim+1);
    }
    if (Opt_params.print_lvl >= 2) {
      oprintf_out("Eigenvalues of scaled RFO matrix.\n");
      int cnt2=0;
      for (i=0; i<dim+1; ++i) {
        oprintf_out("%10.6lf", SRFOevals[i]);
        if (++cnt2 == 6) { oprintf_out("\n"); cnt2 = 0; }
      }
      oprintf_out("\n");

      oprintf_out("First eigenvector (unnormalized) of scaled RFO matrix.\n");
      oprint_array_out(SRFO[0], dim+1);
    }

    // Do intermediate normalization.
    // RFO paper says to scale eigenvector to make the last element equal to 1.
    // During the course of an optimization some evects may appear that are bogus leads
    // - the root following can avoid them.
    for (i=0; i<dim+1; ++i) {
      tval = abs( array_abs_max(SRFO[i], dim)/ SRFO[i][dim] ); // how big is dividing going to make it?
      if (fabs(tval) < Opt_params.rfo_normalization_max) { // same check occurs below for acceptability
        for (j=0;j<dim+1;++j)
          SRFO[i][j] /= SRFO[i][dim];
      }
    }
    if (Opt_params.print_lvl >= 4) {
      oprintf_out("All scaled RFO eigenvectors (rows).\n");
      oprint_matrix_out(SRFO, dim+1, dim+1);
    }

    // Choose which RFO eigenvector to use.
    if ((!rfo_follow_root) || (p_Opt_data->g_iteration() == 1)) {

      rfo_root = Opt_params.rfo_root;
      if (iter == 0)
        oprintf_out("\tGoing to follow RFO solution %d.\n", rfo_root+1);

      while (!symm_rfo_step) {
        symm_rfo_step = coord_combo_is_symmetric(SRFO[rfo_root], dim);

        if (!symm_rfo_step) {
          if (iter == 0)
            oprintf_out("\tRejecting RFO root %d because it breaks the molecular point group.\n", rfo_root+1);
          ++rfo_root;
        }
        else {
          tval = abs( array_abs_max(SRFO[rfo_root], dim)/ SRFO[rfo_root][dim] );
          if (fabs(tval) > Opt_params.rfo_normalization_max) { // matching test in code above
            if (iter == 0)
              oprintf_out("\tRejecting RFO root %d because normalization gives large value.\n", rfo_root+1);
            symm_rfo_step = false;
            ++rfo_root;
          }
        }
        if (rfo_root == dim+1) {
          rfo_root = Opt_params.rfo_root; // no good one found, use the default
          break;
        }
      }
      rfo_follow_root = true;
      array_copy(SRFO[rfo_root], last_iter_evect, dim);
    }
    else { // do root following
      tval = 0;
      for (i=0; i<dim; ++i) { // dot only within H block
        tval2 = array_dot(SRFO[i], last_iter_evect, dim);
        if (tval2 > tval) {
          tval = tval2;
          rfo_root = i;
        }
      }
      array_copy(SRFO[rfo_root], last_iter_evect, dim);
    }
    if (iter == 0)
      oprintf_out("\tUsing RFO vector %d.\n", rfo_root+1);

    // Print only the lowest eigenvalues/eigenvectors
    if (Opt_params.print_lvl >= 2) {
      oprintf_out("\trfo_root is %d\n", rfo_root+1);
      for (i=0; i<dim+1; ++i) {
        if ((SRFOevals[i] < -0.000001) || (i <rfo_root)) {
          oprintf_out("\nScaled RFO eigenvalue %d: %15.10lf (or 2*%-15.10lf)\n", i+1, SRFOevals[i], SRFOevals[i]/2);
          oprintf_out("eigenvector:\n");
          oprint_matrix_out(&(SRFO[i]), 1, dim+1);
          oprintf_out("\n");
        }
      }
    }

    for (j=0; j<dim; ++j)
      dq[j] = SRFO[rfo_root][j]; // leave out last column

    // Project out redundancies in steps.
    // RAK added this projection in 2014, but doesn't seem to be necessary. f,H are already projected.
    project_dq(dq);

    // Zero steps for frozen fragment.  If user has specified.
    for (std::size_t f=0; f<fragments.size(); ++f) {
      if (fragments[f]->is_frozen()) {
        if (iter == 0)
            oprintf_out("\tZero'ing out displacements for frozen fragment %d\n", f+1);
        for (i=0; i<fragments[f]->Ncoord(); ++i)
          dq[ g_coord_offset(f) + i ] = 0.0;
      }
    }

    for (std::size_t I=0; I<interfragments.size(); ++I) {
      if (interfragments[I]->is_frozen()) {
        if (iter == 0)
          oprintf_out("\tZero'ing out displacements for frozen interfragment %d\n", I+1);
        for (int i=0; i<interfragments[I]->Ncoord(); ++i)
          dq[ g_interfragment_coord_offset(I)+i] = 0;
      }
    }
    // Set dq of frozen interfragment coordinates to hard zero.
    for (std::size_t I=0; I<interfragments.size(); ++I) {
      for (int i=0; i< interfragments[I]->Ncoord(); ++i)
        if (interfragments[I]->is_frozen(i))
          dq[g_interfragment_coord_offset(I)+i] = 0.0;
    }

    // Perhaps not necessary:
    // check_intrafragment_zero_angles(dq);

    // get norm |dq| and unit vector in the step direction
    dqtdq = array_dot(dq, dq, dim);

    if (fabs(alpha) > Opt_params.rsrfo_alpha_max) { // don't call it converged if alpha explodes, and give up
      converged = false;
      iter = max_projected_rfo_iter - 1;
    }
    else if (sqrt(dqtdq) < (trust+1e-5))
      converged = true;

    if (iter == 0 && !converged) {

      oprintf_out("\n\tDetermining step-restricting scale parameter for RS-RFO.\n");
      oprintf_out("\tMaximum step size allowed %10.5lf\n\n", trust);
      oprintf_out("\t Iter      |step|        alpha        rfo_root  \n");
      oprintf_out("\t------------------------------------------------\n");
      oprintf_out("\t%5d%12.5lf%14.5lf%12d\n", iter, sqrt(dqtdq), alpha, rfo_root+1);
    }
    else if ( (iter > 0) && !Opt_params.simple_step_scaling)
      oprintf_out("\t%5d%12.5lf%14.5lf%12d\n", iter, sqrt(dqtdq), alpha, rfo_root+1);


    // Find the analytical derivative.
    lambda = -1 * array_dot(fq, dq, dim);
    if (Opt_params.print_lvl >= 2) {
      oprintf_out("dq:\n"); oprint_array_out(dq, dim);
      oprintf_out("fq:\n"); oprint_array_out(fq,dim);
      oprintf_out("\tlambda calculated by (dq^t).(-f)     = %20.10lf\n", lambda);
    }

    // Calculate derivative of step size wrt alpha.
    // Equation 20, Besalu and Bofill, Theor. Chem. Acc., 1999, 100:265-274
    sum = 0;
    for (i=0; i<dim; i++)
      sum += ( pow(array_dot( Hevects[i], fq, dim),2) ) / ( pow(( h[i]-lambda*alpha ),3) );

    analyticDerivative = 2*lambda / (1+alpha*dqtdq ) * sum;
    if (Opt_params.print_lvl >= 2)
      oprintf_out("\tAnalytic derivative d(norm)/d(alpha) = %20.10lf\n", analyticDerivative);

    // Calculate new scaling alpha value.
    // Equation 20, Besalu and Bofill, Theor. Chem. Acc., 1999, 100:265-274
    alpha += 2*(trust * sqrt(dqtdq) - dqtdq) / analyticDerivative;
  }

  if ((iter > 0) && !Opt_params.simple_step_scaling)
    oprintf_out("\t------------------------------------------------\n");

  // Crude/old way to limit step size if restricted step algorithm failed.
  if (!converged)
    apply_intrafragment_step_limit(dq);

  if (Opt_params.print_lvl >= 3) {
    oprintf_out("\tFinal scaled step dq:\n");
    oprint_matrix_out(&dq, 1, dim);
  }

  // Save and print out RFO eigenvector(s)
  p_Opt_data->set_rfo_eigenvector(dq);
  dqtdq = array_dot(dq, dq, dim);
  double rfo_dqnorm = sqrt( dqtdq );
  array_copy(dq, rfo_u, dim);
  array_normalize(rfo_u, dim);

  oprintf_out("\tNorm of target step-size %10.5lf\n", rfo_dqnorm);

 // get gradient and hessian in step direction
  rfo_g = -1 * array_dot(fq, rfo_u, dim);
  rfo_h = 0;
  for (i=0; i<dim; ++i)
    rfo_h += rfo_u[i] * array_dot(Hevects[i], rfo_u, dim);

  DE_projected = DE_rfo_energy(rfo_dqnorm, rfo_g, rfo_h);
  oprintf_out("\tProjected energy change by RFO approximation: %20.10lf\n", DE_projected);

  free_matrix(Hevects);
  free_array(h);
  free_matrix(RFO);
  free_matrix(SRFO);
  free_array(SRFOevals);

  std::vector<int> lin_angles = validate_angles(dq);
  if (!lin_angles.empty())
    throw(INTCO_EXCEPT("New linear angles", lin_angles));

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
  oprintf_out("Step-size in mass-weighted internal coordinates: %20.10lf\n", N);
*/
  for (std::size_t f=0; f<fragments.size(); ++f) {
    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
      oprintf_out("\tDisplacements for frozen fragment %d skipped.\n", f+1);
      continue;
    }
    fragments[f]->displace(&(dq[g_coord_offset(f)]), &(fq[g_coord_offset(f)]), g_atom_offset(f));
  }

  // do displacements for interfragment coordinates
  for (std::size_t I=0; I<interfragments.size(); ++I) {
    if (interfragments[I]->is_frozen() || Opt_params.freeze_interfragment) {
      oprintf_out("\tDisplacements for frozen interfragment %d skipped.\n", I+1);
      continue;
    }
    interfragments[I]->orient_fragment( &(dq[g_interfragment_coord_offset(I)]),
                                        &(fq[g_interfragment_coord_offset(I)]) );
  }

  // fix rotation matrix for rotations in QCHEM EFP code
  for (std::size_t I=0; I<fb_fragments.size(); ++I)
    fb_fragments[I]->displace( I, &(dq[g_fb_fragment_coord_offset(I)]) );

  symmetrize_geom(true); // now symmetrize the geometry for next step

  if (Opt_params.print_lvl > 1) {
    oprintf_out("Symmetrized geometry\n");
    double **g = g_geom_2D();
    oprint_matrix_out(g, g_natom(), 3);
    free_matrix(g);
  }

/* Test step sizes
  double *x_after = g_geom_array();
  double *masses = g_masses();

  double sum = 0.0;
  for (int i=0; i<g_natom(); ++i)
    for (int xyz=0; xyz<3; ++xyz)
      sum += (x_before[3*i+xyz] - x_after[3*i+xyz]) * (x_before[3*i+xyz] - x_after[3*i+xyz])
               * masses[i];

  sum = sqrt(sum);
  oprintf_out("Step-size in mass-weighted cartesian coordinates [bohr (amu)^1/2] : %20.10lf\n", sum);
  free_array(x_after);
  free_array(x_before);
  free_array(masses);
*/

  // save values in step data
  p_Opt_data->save_step_info(DE_projected, rfo_u, rfo_dqnorm, rfo_g, rfo_h);

  free_array(rfo_u);

  // Before quitting, make sure step is reasonable.  It should only be screwball if we are using the
  // "First Guess" after the back-transformation failed.
  double norm = sqrt(array_dot(dq, dq, dim));

  if (norm > 10 * trust) {
    p_Opt_data->restore_previous_consecutive_backsteps(); // it has already been reset to 0
    throw(BAD_STEP_EXCEPT("Step is far too large.\n"));
  }

} // end take RFO step

}
