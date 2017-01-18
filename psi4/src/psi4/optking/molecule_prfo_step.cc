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

// compute change in energy according to P-RFO approximation
// see J. Phys. Chem. 1985, 89, 52-57
inline double DE_rfo_energy(double rfo_t, double rfo_g, double rfo_h) {
  return (rfo_t * rfo_g + 0.5 * rfo_t * rfo_t * rfo_h)/(1 + rfo_t*rfo_t);
}

void MOLECULE::prfo_step(void) {
  double **Horig = p_Opt_data->g_H_pointer();
  double *fq = p_Opt_data->g_forces_pointer();
  double *dq = p_Opt_data->g_dq_pointer();
  int Nintco = Ncoord();
  int rfo_root; // ultimately, should be array of roots to maximize
  int cnt_i;

  oprintf_out("\tTaking PRFO optimization step.\n");

  // don't use Horig anymore -it's the pointer to the good, original Hessian
  double **H = matrix_return_copy(Horig, Nintco, Nintco);

  if (Opt_params.print_lvl >= 3) {
    oprintf_out("\nHessian matrix\n");
    oprint_matrix_out(H, Nintco, Nintco);
  }

  // diagonalize H (technically, only have to semi-diagonalize)
  double * lambda = init_array(Nintco);
  opt_symm_matrix_eig(H, Nintco, lambda);
  double **H_evects = H; // rename for clarity

  if (Opt_params.print_lvl >= 3) {
    oprintf_out("\n\tEigenvalues of Hessian \n");
    oprint_matrix_out(&lambda, 1, Nintco);
    oprintf_out( "\n\tEigenvectors of Hessian (rows) \n");
    oprint_matrix_out(H_evects, Nintco, Nintco);
  }

  // construct diagonalized Hessian with evals on diagonal
  double **H_diag = init_matrix(Nintco,Nintco);
  for (int i=0; i<Nintco; ++i)
    H_diag[i][i] = lambda[i];

  if (Opt_params.print_lvl >= 3) {
    oprintf_out("\n\tH_diag\n");
    oprint_matrix_out(H_diag, Nintco, Nintco);
  }

  // the number of degrees along which to MAXIMIZE; assume 1 for now
  int mu = 1;

  // For now, use vector stored in rfo_root to choose which modes to maximize
  // in future will change to store more than one of these
  if (p_Opt_data->g_iteration() == 1 || !Opt_params.rfo_follow_root) {
    rfo_root = Opt_params.rfo_root;
    oprintf_out("\tMaximizing along %d lowest eigenvalue of Hessian.\n", rfo_root+1);
  }
  else { // do dynamic root-following
    double * rfo_old_evect = p_Opt_data->g_rfo_eigenvector_pointer();
    double max_overlap = 0;
    for (int i=0; i<Nintco; ++i) {
      double tval = fabs (array_dot(H_diag[i], rfo_old_evect, Nintco));
      if (tval > max_overlap) {
        max_overlap = tval;
        rfo_root = i;
      }
    }
    oprintf_out("\tMaximizing along Hessian eigenvalue %d whose \
      eigenvector has maximal overlap with previous step.\n", rfo_root+1);
  }
  p_Opt_data->set_rfo_eigenvector(H_diag[rfo_root]);

  // transform gradient
  double *f_q_Hevect_basis = init_array(Nintco);
  opt_matrix_mult(H_evects, 0, &fq, 1, &f_q_Hevect_basis, 1, Nintco, Nintco, 1, 0);
  if (Opt_params.print_lvl >= 2) {
    oprintf_out("\tInternal forces in au,\n");
    oprint_array_out(fq, Nintco);
    oprintf_out("\tInternal forces in au, in Hevect basis.\n");
    oprint_array_out(f_q_Hevect_basis, Nintco);
  }

  // Build RFO-max.
  double **rfo_max = init_matrix(mu+1, mu+1);

  rfo_max[0][0] = H_diag[rfo_root][rfo_root];

  rfo_max[0][1] = -f_q_Hevect_basis[rfo_root];
  rfo_max[1][0] = -f_q_Hevect_basis[rfo_root];

  if (Opt_params.print_lvl >= 2) {
    oprintf_out("\n RFO max \n");
    oprint_matrix_out(rfo_max,mu+1,mu+1);
  }

  // Build RFO-min.
  double **rfo_min = init_matrix(Nintco-mu+1,Nintco-mu+1);
  cnt_i = 0;
  for (int i=0; i<Nintco; ++i) {
    if (i != rfo_root) {
      rfo_min[cnt_i][cnt_i] = H_diag[i][i];
      ++cnt_i;
    }
  }

  cnt_i = 0;
  for (int i=0; i<Nintco; ++i) {
    if (i != rfo_root) {
      rfo_min[Nintco-mu][cnt_i] = -f_q_Hevect_basis[i];
      rfo_min[cnt_i][Nintco-mu] = rfo_min[Nintco-mu][cnt_i];
      ++cnt_i;
    }
  }

  if (Opt_params.print_lvl >= 2) {
    oprintf_out("\n RFO min \n");
    oprint_matrix_out(rfo_min, Nintco-mu+1 , Nintco-mu+1);
  }

  double* max_evals = init_array(mu+1);
  double* min_evals = init_array(Nintco-mu+1);

  //find eigenvectors and eigenvalues of rfo_max and rfo_min
  //rfo_max and rfo_min now contain eigenvectors as rows.  min/max_evals are eigenvalues
  opt_symm_matrix_eig(rfo_max, mu+1,max_evals);
  opt_symm_matrix_eig(rfo_min, Nintco-mu+1,min_evals);

  if (Opt_params.print_lvl >= 3) {
    oprintf_out("\n RFO min eigenvectors (rows) before normalization\n");
    oprint_matrix_out(rfo_min, Nintco-mu+1, Nintco-mu+1);

    oprintf_out("\n RFO max eigenvectors (rows) before normalization\n");
    oprint_matrix_out(rfo_max, mu+1, mu+1);
  }

  if (Opt_params.print_lvl >= 1) {
    oprintf_out("\n RFO min eigenvalues\n");
    oprint_matrix_out(&min_evals,1, Nintco-mu+1);

    oprintf_out("\n RFO max eigenvalues\n");
    oprint_matrix_out(&max_evals,1, mu+1);
  }

  //Normalize all eigenvectors.
  for (int i=0; i<mu+1; ++i) {
    // how big is dividing going to make it?
    double tval = abs( array_abs_max(rfo_max[i], mu) / rfo_max[i][mu] );
    if (fabs(tval) < Opt_params.rfo_normalization_max) {
      for (int j=0; j<mu+1; ++j)
        rfo_max[i][j] /= rfo_max[i][mu];
    }
  }
  if (Opt_params.print_lvl >= 3) {
    oprintf_out("\n RFO_max eigenvectors (rows)\n");
    oprint_matrix_out(rfo_max, mu+1, mu+1);
  }

  //rfo_min contains normalized eigenvectors as rows
  for (int i=0; i<Nintco-mu+1; ++i) {
    double tval = abs( array_abs_max(rfo_min[i], Nintco-mu) / rfo_min[i][Nintco-mu] );
    if (fabs(tval) < Opt_params.rfo_normalization_max) {
      for (int j=0;j<Nintco-mu+1;++j)
        rfo_min[i][j] /= rfo_min[i][Nintco-mu];
    }
  }
  if (Opt_params.print_lvl >= 3) {
    oprintf_out("\nRFO_min eigenvectors (rows)\n");
    oprint_matrix_out(rfo_min, Nintco-mu+1, Nintco-mu+1);
  }

  double * rfo_step_Hevect_basis = init_array(Nintco);

  // extract step with highest eigenvalue?
  rfo_step_Hevect_basis[rfo_root] = rfo_max[mu][0]; // drop last (2nd) entry

  // extract step with lowest eigenvalue
  cnt_i = 0;
  for (int i=0; i<Nintco; ++i) {
    if (i!= rfo_root) {
      rfo_step_Hevect_basis[i] = rfo_min[0][cnt_i];
      ++cnt_i;
    }
  }

  if (Opt_params.print_lvl >= 2) {
    oprintf_out( "\nRFO step in Hevect basis\n");
    oprint_matrix_out(&rfo_step_Hevect_basis, 1, Nintco);
  }

  // transform back into original basis.
  // write to old dq pointer ?
  opt_matrix_mult(H_evects, 1, &rfo_step_Hevect_basis, 1, &dq, 1, Nintco, Nintco, 1, 0);

  if (Opt_params.print_lvl >= 2) {
    oprintf_out( "\nRFO step in original basis\n");
    oprint_matrix_out(&dq, 1, Nintco);
  }

  apply_intrafragment_step_limit(dq);

  // try to get but with a single extrapolated energy change

  double rfo_dqnorm = sqrt( array_dot(dq, dq, Nintco) );
  double rfo_g = -1 * array_dot(fq, dq, Nintco);

  oprintf_out("\tNorm of target step-size %10.5lf\n", rfo_dqnorm);

  double rfo_h = 0;
  for (int i=0; i<Nintco; ++i)
    rfo_h += dq[i] * array_dot(Horig[i], dq, Nintco);

  double DE_projected = DE_rfo_energy(rfo_dqnorm, rfo_g, rfo_h);

//calculating the projected energy change
/*
double rfo_g_max;
double rfo_h_max;
double DE_projected_max;
double rfo_g_min;
double rfo_h_min;
double DE_projected_min;
double DE_projected;
double rfo_dqnorm_max;
double rfo_dqnorm_min;

//get gradient and hessian in step direction for min and max
  rfo_dqnorm_max = sqrt( array_dot(dq_max, dq_max, Nintco) );
  rfo_g_max = -1 * array_dot(fq, dq_max, Nintco);
  rfo_h_max = 0;
  for (i=0; i<Nintco; ++i) {
    rfo_h_max += dq_max[i] * array_dot(H[i], dq_max, Nintco);
  }
  DE_projected_max = DE_rfo_energy(rfo_dqnorm_max, rfo_g_max, rfo_h_max);

  rfo_dqnorm_min = sqrt( array_dot(dq_min, dq_min, Nintco) );
  rfo_g_min = -1 * array_dot(fq, dq_min, Nintco);
  rfo_h_min = 0;
  for (i=0; i<Nintco; ++i) {
    rfo_h_min += dq_min[i] * array_dot(H[i], dq_min, Nintco);
  }
  DE_projected_min = DE_rfo_energy(rfo_dqnorm_min, rfo_g_min, rfo_h_min);

  DE_projected = DE_projected_min+DE_projected_max;
  oprintf_out("\tProjected energy change by P-RFO approximation: %20.10lf\n", DE_projected);
*/


  // do displacements for each fragment separately
 for (std::size_t f=0; f<fragments.size(); ++f) {
    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
      oprintf_out("\tDisplacements for frozen fragment %d skipped.\n", f+1);
      continue;
    }
    fragments[f]->displace(&(dq[g_coord_offset(f)]), &(fq[g_coord_offset(f)]), g_atom_offset(f));
  }

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

  symmetrize_geom(); // now symmetrize the geometry for next step

  // save values in step data
  p_Opt_data->save_step_info(DE_projected, dq, rfo_dqnorm, rfo_g, rfo_h);


  return;
}

}
