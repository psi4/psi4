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

// compute change in energy according to P-RFO approximation
// see J. Phys. Chem. 1985, 89, 52-57
inline double DE_rfo_energy(double rfo_t, double rfo_g, double rfo_h) {
  return (rfo_t * rfo_g + 0.5 * rfo_t * rfo_t * rfo_h)/(1 + rfo_t*rfo_t);
}

void MOLECULE::prfo_step(void) {
  double **Horig = p_Opt_data->g_H_pointer();
  double *fq = p_Opt_data->g_forces_pointer();
  double *dq = p_Opt_data->g_dq_pointer();
  int Nintco = g_nintco();
  int rfo_root; // ultimately, should be array of roots to maximize
  int cnt_i, cnt_j;

  // don't use Horig anymore -it's the pointer to the good, original Hessian
  double **H = matrix_return_copy(Horig, Nintco, Nintco);

  fprintf(outfile,"\nHessian matrix\n");  
  print_matrix(outfile, H, Nintco, Nintco);

  // diagonalize H (technically, only have to semi-diagonalize)
  double * lambda = init_array(Nintco);
  opt_symm_matrix_eig(H, Nintco, lambda);
  double **H_evects = H; // rename for clarity
 
  fprintf(outfile,"\n\tEigenvalues of Hessian \n");
  print_matrix(outfile, &lambda, 1, Nintco);
  fprintf(outfile, "\n\tEigenvectors of Hessian (rows) \n");
  print_matrix(outfile, H_evects, Nintco, Nintco);

  // construct diagonalized Hessian with evals on diagonal
  double **H_diag = init_matrix(Nintco,Nintco);
  for (int i=0; i<Nintco; ++i)
    H_diag[i][i] = lambda[i];

  fprintf(outfile,"\n\tH_diag\n");
  print_matrix(outfile, H_diag, Nintco, Nintco);

  // the number of degrees along which to MAXIMIZE; assume 1 for now
  int mu = 1;

  // For now, use vector stored in rfo_root to choose which modes to maximize
  // in future will change to store more than one of these
  if (p_Opt_data->g_iteration() == 1 || !Opt_params.rfo_follow_root) {
    rfo_root = Opt_params.rfo_root;
    fprintf(outfile,"\tMaximizing along %d lowest eigenvalue of Hessian.\n", rfo_root+1);
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
    fprintf(outfile,"\tMaximizing along Hessian eigenvalue %d whose \
      eigenvector has maximal overlap with previous step.\n", rfo_root+1);
  }
  p_Opt_data->set_rfo_eigenvector(H_diag[rfo_root]);

  // transform gradient
  double *f_q_Hevect_basis = init_array(Nintco);
  opt_matrix_mult(H_evects, 0, &fq, 1, &f_q_Hevect_basis, 1, Nintco, Nintco, 1, 0);

  // Build RFO-max.
  double **rfo_max = init_matrix(mu+1, mu+1);

  rfo_max[0][0] = H_diag[rfo_root][rfo_root];

  rfo_max[0][1] = -f_q_Hevect_basis[rfo_root];
  rfo_max[1][0] = -f_q_Hevect_basis[rfo_root];

  fprintf(outfile,"\n RFO max \n");
  print_matrix(outfile,rfo_max,mu+1,mu+1);

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
   
  fprintf(outfile,"\n RFO min \n");
  print_matrix(outfile, rfo_min, Nintco-mu+1 , Nintco-mu+1);
 
  double* max_evals = init_array(mu+1);
  double* min_evals = init_array(Nintco-mu+1);
  
  //find eigenvectors and eigenvalues of rfo_max and rfo_min
  //rfo_max and rfo_min now contain eigenvectors as rows.  min/max_evals are eigenvalues
  opt_symm_matrix_eig(rfo_max, mu+1,max_evals);
  opt_symm_matrix_eig(rfo_min, Nintco-mu+1,min_evals);
 
  fprintf(outfile,"\n RFO min eigenvectors (rows) before normalization\n");
  print_matrix(outfile, rfo_min, Nintco-mu+1, Nintco-mu+1);

  fprintf(outfile,"\n RFO max eigenvectors (rows) before normalization\n");
  print_matrix(outfile, rfo_max, mu+1, mu+1);

  fprintf(outfile,"\n RFO min eigenvalues\n");
  print_matrix(outfile, &min_evals,1, Nintco-mu+1);

  fprintf(outfile,"\n RFO max eigenvalues\n");
  print_matrix(outfile, &max_evals,1, mu+1);

  //Normalize all eigenvectors.
  for (int i=0; i<mu+1; ++i) {
    double tval = rfo_max[i][mu];
    if (fabs(tval) > Opt_params.rfo_normalization_min) {
      for (int j=0; j<mu+1; ++j)
        rfo_max[i][j] /= rfo_max[i][mu];
    }
  }
  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile,"\n RFO_max eigenvectors (rows)\n");
    print_matrix(outfile, rfo_max, mu+1, mu+1);
  }

  //rfo_min contains normalized eigenvectors as rows
  for (int i=0; i<Nintco-mu+1; ++i) {
    double tval = rfo_min[i][Nintco-mu];
    if (fabs(tval) > Opt_params.rfo_normalization_min) {
      for (int j=0;j<Nintco-mu+1;++j)
        rfo_min[i][j] /= rfo_min[i][Nintco-mu];
    }
  }
  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile,"\nRFO_min eigenvectors (rows)\n");
    print_matrix(outfile, rfo_min, Nintco-mu+1, Nintco-mu+1);
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
  
  fprintf(outfile, "\nRFO step in Hevect basis\n"); 
  print_matrix(outfile, &rfo_step_Hevect_basis, 1, Nintco);
 
  // transform back into original basis.
  // write to old dq pointer ? 
  opt_matrix_mult(H_evects, 1, &rfo_step_Hevect_basis, 1, &dq, 1, Nintco, Nintco, 1, 0);

  fprintf(outfile, "\nRFO step in original basis\n"); 
  print_matrix(outfile, &dq, 1, Nintco);
 
  apply_intrafragment_step_limit(dq);

  // try to get but with a single extrapolated energy change
 
  double rfo_dqnorm = sqrt( array_dot(dq, dq, Nintco) );
  double rfo_g = -1 * array_dot(fq, dq, Nintco);

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
  fprintf(outfile,"\tProjected energy change by P-RFO approximation: %20.10lf\n", DE_projected);
*/


  // do displacements for each fragment separately
 for (int f=0; f<fragments.size(); ++f) {
    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
      fprintf(outfile,"\tDisplacements for frozen fragment %d skipped.\n", f+1);
      continue;
    }
    fragments[f]->displace(&(dq[g_intco_offset(f)]), &(fq[g_intco_offset(f)]), g_atom_offset(f));
  }

  for (int I=0; I<interfragments.size(); ++I) {
    if (interfragments[I]->is_frozen() || Opt_params.freeze_interfragment) {
      fprintf(outfile,"\tDisplacements for frozen interfragment %d skipped.\n", I+1);
      continue;
    }
    interfragments[I]->orient_fragment( &(dq[g_interfragment_intco_offset(I)]),
                                        &(fq[g_interfragment_intco_offset(I)]) );
  }

  symmetrize_geom(); // now symmetrize the geometry for next step

  // save values in step data
  p_Opt_data->save_step_info(DE_projected, dq, rfo_dqnorm, rfo_g, rfo_h);

  fflush(outfile);
  return;
}

}

