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

/*! \file molecule.cc
    \ingroup optking
    \brief molecule class (really, molecular system class)
*/

#include "molecule.h"

#include <iostream>
#include <sstream>

#include "linear_algebra.h"
#include "print.h"
#include "atom_data.h"
#include "physconst.h"

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
  int dim = g_nintco();
  int natom = g_natom();
  double tval, tval2;
  double *fq = p_Opt_data->g_forces_pointer();
  double **H = p_Opt_data->g_H_pointer();
  double *dq = p_Opt_data->g_dq_pointer();

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

  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile,"Original, unscaled RFO mat\n");
    print_matrix(outfile, RFO, dim+1, dim+1);
  }

  int rfo_root, f;     // root to follow
  double rfo_eval;
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

  //Iterative sequence to find alpha; we'll give it 15 tries
  int iter = -1;
  while (!converged && iter<16) {

    ++iter;

    // If we get to iteration 14, and we haven't not yet converged, than
    // bail out on the restricted-step algorithm.  Set alpha=1 and apply the crude
    // intrafragment_step_limit below.
    if (iter == 16) {
      fprintf(outfile, "\tFailed to converge alpha.  Doing simple step scaling instead.\n");
      alpha = 1;
    }
    else if (Opt_params.simple_step_scaling) // If simple_step_scaling is on, not an iterative method.
      iter = 16;

    // Scale the RFO matrix.
    for (i=0; i<=dim; i++) {
      for (j=0; j<dim; j++) {
        SRFO[j][i] = RFO[j][i] / alpha;
      };
      SRFO[dim][i] = RFO[dim][i];
    }
    if (Opt_params.print_lvl >= 3) {
      fprintf(outfile,"Scaled RFO matrix.\n");
      print_matrix(outfile,  SRFO, dim+1, dim+1);
      fflush(outfile);
    }

    //Find the eigenvectors and eigenvalues of RFO matrix.
    opt_asymm_matrix_eig(SRFO, dim+1, SRFOevals);
    if (Opt_params.print_lvl >= 3) {
      fprintf(outfile,"Eigenvectors of scaled RFO matrix.\n");
      print_matrix(outfile, SRFO, dim+1,dim+1);
      fprintf(outfile,"Eigenvalues of scaled RFO matrix.\n");
      for (i=0; i<dim+1; ++i)
        fprintf(outfile,"%10.6lf", SRFOevals[i]);
      fprintf(outfile,"\n");
    }

    // Do intermediate normalization.  
    // RFO paper says to scale eigenvector to make the last element equal to 1.
    // During the course of an optimization some evects may appear that are bogus leads
    // - the root following can avoid them. 
    for (i=0; i<dim+1; ++i) {
      tval = SRFO[i][dim];
      if (fabs(tval) > Opt_params.rfo_normalization_min) {
        for (j=0;j<dim+1;++j)
          SRFO[i][j] /= SRFO[i][dim];
      }
    }
    if (Opt_params.print_lvl >= 3) {
      fprintf(outfile,"Scaled RFO eigenvectors (rows).\n");
      print_matrix(outfile, SRFO, dim+1, dim+1);
    }

    // Choose which RFO eigenvector to use.
    // if not root following, then use rfo_root'th lowest eigenvalue; default is 0 (lowest)
    if ( (!Opt_params.rfo_follow_root) || (p_Opt_data->g_iteration() == 1)) {
      rfo_root = Opt_params.rfo_root;
      if (iter == 0)
        fprintf(outfile,"\tGoing to follow RFO solution %d.\n", rfo_root+1);

      // Now test RFO eigenvector and make sure that it is totally symmetric.
      while (!symm_rfo_step) {
        symm_rfo_step = intco_combo_is_symmetric(SRFO[rfo_root], dim);

        if (!symm_rfo_step) {
          fprintf(outfile,"\tRejecting RFO root %d because it breaks the molecular point group.\n", rfo_root+1);
          fprintf(outfile,"\tIf you are doing an energy minimization, there may exist a lower-energy, ");
          fprintf(outfile,"structure with less symmetry.\n");
          ++rfo_root;
        }

        if (rfo_root == dim+1) // quit in the unlikely event we've checked them all
        break;
      }
    }
    else { // do root following
      tval = 0;
      for (i=0; i<dim; ++i) { // dot only within H block, excluding rows and columns approaching 0
        tval2 = array_dot(SRFO[i], rfo_old_evect, dim);
        if (tval2 > tval) {
          tval = tval2;
          rfo_root = i;
        }
      }
      if (iter == 0)
        fprintf(outfile,"RFO vector %d overlaps most with earlier step.\n", rfo_root+1);
    }

    // Print only the lowest eigenvalues/eigenvectors
    if (Opt_params.print_lvl >= 2) {
      for (i=0; i<dim+1; ++i) {
        if ((SRFOevals[i] < 0.0) || (i <rfo_root)) {
          fprintf(outfile,"Scaled RFO eigenvalue %d: %15.10lf (or 2*%-15.10lf)\n", i+1, SRFOevals[i], SRFOevals[i]/2);
          fprintf(outfile,"eigenvector:\n");
          print_matrix(outfile, &(SRFO[i]), 1, dim+1);
        }
      }
    }

    for (j=0; j<dim; ++j)
      dq[j] = SRFO[rfo_root][j]; // leave out last column

    // Zero steps for frozen fragment.  If user has specified.
    for (f=0; f<fragments.size(); ++f) {
      if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
        fprintf(outfile,"\tZero'ing out displacements for frozen fragment %d\n", f+1);
        for (i=0; i<fragments[f]->g_nintco(); ++i)
          dq[ g_intco_offset(f) + i ] = 0.0;
      }
    }

    // Perhaps not necessary:
    // check_intrafragment_zero_angles(dq);
	
    // get norm |dq| and unit vector in the step direction
    dqtdq = array_dot(dq, dq, dim);

    if (sqrt(dqtdq) < (trust+1e-5))
      converged = true;

    if (iter == 0 && !converged) {
      fprintf(outfile,"\n\tDetermining step-restricting scale parameter for RS-RFO.\n");
      fprintf(outfile,"\tMaximum step size allowed %10.5lf\n\n", trust);
      fprintf(outfile,"\t Iter      |step|      alpha \n");
      fprintf(outfile,"\t-------------------------------\n");
      fprintf(outfile,"\t%5d%12.5lf%12.5lf\n", iter, sqrt(dqtdq), alpha);
    }
    else if ( (iter > 0) && !Opt_params.simple_step_scaling)
      fprintf(outfile,"\t%5d%12.5lf%12.5lf\n", iter, sqrt(dqtdq), alpha);
    fflush(outfile);

    // Find the analytical derivative.
    lambda = -1 * array_dot(fq, dq, dim);
    if (Opt_params.print_lvl >= 3)
      fprintf(outfile,"\tlambda calculated by (dq^t).(-f) = %20.10lf\n", lambda);

    // Calculate derivative of step size wrt alpha.
    // Equation 20, Besalu and Bofill, Theor. Chem. Acc., 1999, 100:265-274
    sum = 0;
    for (i=0; i<dim; i++)
      sum += ( pow(array_dot( Hevects[i], fq, dim),2) ) / ( pow(( h[i]-lambda*alpha ),3) );

    analyticDerivative = 2*lambda / (1+alpha*dqtdq ) * sum;
    if (Opt_params.print_lvl >= 3)
      fprintf(outfile,"\tAnalytic derivative d(norm)/d(alpha) = %20.10lf\n", analyticDerivative);

    // Calculate new scaling alpha value.
    // Equation 20, Besalu and Bofill, Theor. Chem. Acc., 1999, 100:265-274
    alpha += 2*(trust * sqrt(dqtdq) - dqtdq) / analyticDerivative;
    if (Opt_params.print_lvl >= 3)
      fprintf(outfile,"\tNew alpha value = %20.10lf\n", alpha);
  }

  if ((iter > 0) && !Opt_params.simple_step_scaling)
    fprintf(outfile,"\t-------------------------------\n");

  // Crude/old way to limit step size if restricted step algorithm failed.
  if (!converged)
    apply_intrafragment_step_limit(dq);

  // Save and print out RFO eigenvector(s) 
  p_Opt_data->set_rfo_eigenvector(dq);
  dqtdq = array_dot(dq, dq, dim);
  double rfo_dqnorm = sqrt( dqtdq );
  array_copy(dq, rfo_u, dim);
  array_normalize(rfo_u, dim);

 // get gradient and hessian in step direction
  rfo_g = -1 * array_dot(fq, rfo_u, dim);
  rfo_h = 0;
  for (i=0; i<dim; ++i)
    rfo_h += rfo_u[i] * array_dot(Hevects[i], rfo_u, dim);

  DE_projected = DE_rfo_energy(rfo_dqnorm, rfo_g, rfo_h);
  fprintf(outfile,"\tProjected energy change by RFO approximation: %20.10lf\n", DE_projected);

  free_matrix(Hevects);
  free_array(h);
  free_matrix(RFO);
  free_matrix(SRFO);
  free_array(SRFOevals);

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

