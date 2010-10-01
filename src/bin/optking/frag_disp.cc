/*! \file    frag_disp.cc
    \ingroup OPT10
    \brief   displace fragment geometry only dq changes to values of coordinates
*/

#include "frag.h"
#include "linear_algebra.h"
#include "opt_data.h"

#include <cmath>

#define EXTERN
#include "globals.h"

using psi::outfile;

namespace opt {

// dq - displacements to be performed
void FRAG::displace(double *dq, bool print_disp) {
  int i,j;
  int Ncarts = 3 * natom;
  int Nints = intcos.size();
  double dq_max, dq_rms, **G_inv, *new_q;

  fprintf(outfile,"\n\tBack-transformation to cartesian coordinates...\n");
  fprintf(outfile,"\t---------------------------------------------\n");
  fprintf(outfile,"\t Iter     RMS difference    MAX difference   \n");
  fprintf(outfile,"\t---------------------------------------------\n");

  //compute_intco_values(); // lock in orientation of torsions
  fix_tors_near_180(); // subsequent computes will modify torsional values
  double * q_orig = intco_values();
  if (Opt_params.print_lvl > 2) {
    fprintf(outfile,"\nOriginal q internal coordinates\n");
    for (i=0; i<Nints; ++i) fprintf(outfile, "\t%15.10lf\n", q_orig[i]);
  }

  double * q_target = init_array(Nints);
  for (i=0; i<Nints; ++i)
    q_target[i] = q_orig[i] + dq[i];
  if (Opt_params.print_lvl > 2) {
    fprintf(outfile,"\nTarget q internal coordinates\n");
    for (i=0; i<Nints; ++i) fprintf(outfile, "\t%15.10lf\n", q_target[i]);
  }

  double * new_geom   = g_geom_array(); // cart geometry to start each iter
  double * first_geom = init_array(Ncarts); // first try at back-transformation
  double * dx = init_array(Ncarts);
  double * tmp_v_Nints = init_array(Nints);
  double **B = init_matrix(Nints, Ncarts);
  double **G = init_matrix(Nints, Nints);

  bool bmat_iter_done = false;
  int bmat_iter_cnt = 0;

  while ( !bmat_iter_done && (bmat_iter_cnt < Opt_params.bt_max_iter) ) {

    // G = BuBt ; u = unit matrix; G = BBt
    compute_B(B);
    opt_matrix_mult(B, 0, B, 1, G, 0, Nints, Ncarts, Nints, 0);

    // u B^t (G_inv dq) = dx
    G_inv = symm_matrix_inv(G, Nints, true);
    opt_matrix_mult(G_inv, 0, &dq, 1, &tmp_v_Nints, 1, Nints, Nints, 1, 0);
    opt_matrix_mult(B, 1, &tmp_v_Nints, 1, &dx, 1, Ncarts, Nints, 1, 0);
    free_matrix(G_inv);

    /*if (Opt_params.print_lvl > 2) {
      fprintf(outfile,"\ndx increments\n");
      for (i=0; i<Ncarts; ++i) fprintf(outfile,"%15.10lf\n", dx[i]); } */

    for (i=0; i<Ncarts; ++i)
      new_geom[i] += dx[i];

    if (bmat_iter_cnt == 0) { // save first try in case doesn't converge
      for (i=0; i<Ncarts; ++i)
        first_geom[i] = new_geom[i];
    }

    set_geom_array(new_geom);
    //compute_intco_values();
    new_q = intco_values();
    /* if (Opt_params.print_lvl > 2) {
      fprintf(outfile,"\nObtained q internal coordinates\n");
      for (i=0; i<Nints; ++i) fprintf(outfile, "\t%15.10lf\n", new_q[i]);}*/
    for (i=0; i< Nints; ++i)
      dq[i] = q_target[i] - new_q[i]; // calculate dq from _target_
    free_array(new_q);

    // Test for convergence of iterations
    dq_max = dq_rms = 0.0;
    for (i=0; i<Nints; ++i) {
      if (fabs(dq[i]) > dq_max) dq_max = fabs(dq[i]);
      dq_rms += dq[i]*dq[i];
    }
    dq_rms = sqrt(dq_rms / ((double) Nints));

    if ( dq_rms < Opt_params.bt_dq_conv_rms && dq_max < Opt_params.bt_dq_conv_max )
      bmat_iter_done = true;

    fprintf (outfile,"\t%5d   %12.1e     %12.1e\n", bmat_iter_cnt+1, dq_rms, dq_max);

    ++bmat_iter_cnt;
  }

  fprintf(outfile,"\t---------------------------------------------\n");

  if (bmat_iter_done) {
    fprintf(outfile, "\n\tSuccessfully converged to displaced geometry.\n");
  }
  else {
    fprintf(outfile,"Could not converge backtransformation in %d iterations.\n", bmat_iter_cnt);
    fprintf(outfile,"Using first guess instead.\n");
    set_geom_array(first_geom);
  }

  /* Set dq to final, total displacement achieved */
  //compute_intco_values();
  new_q = intco_values();
  for (i=0; i< Nints; ++i)
    dq[i] = new_q[i] - q_orig[i]; // calculate dq from _target_

  if (print_disp) {

    double *f_q = p_Opt_data->g_forces_pointer(); // for printing

    fprintf(outfile,"\n\t---Internal Coordinate Step in ANG or DEG, aJ/ANG or AJ/DEG (%d) ---\n", intcos.size());
    fprintf(outfile,"\t ----------------------------------------------------------------------\n");
    fprintf(outfile,"\t Coordinate             Previous        Force       Change         New \n");
    fprintf(outfile,"\t ----------             --------       ------       ------       ------\n");
    for (i=0; i<intcos.size(); ++i)
      intcos.at(i)->print_disp(outfile, q_orig[i], f_q[i], dq[i], new_q[i]);
    fprintf(outfile,"\t ----------------------------------------------------------------------\n");
  }
  free_array(new_q);

  free_matrix(G);
  free_array(new_geom);
  free_array(first_geom);
  free_array(dx);
  free_array(tmp_v_Nints);
  free_matrix(B);

  free_array(q_target);
  free_array(q_orig);

  return;
}

}

