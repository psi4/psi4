/*! \file    frag_disp.cc
    \ingroup optking
    \brief   displace fragment geometry only dq changes to values of coordinates
*/

#include "frag.h"
#include "linear_algebra.h"
#include "opt_data.h"

#include <cmath>

#define EXTERN
#include "globals.h"

namespace opt {

// dq - displacements to be performed
void FRAG::displace(double *dq, bool print_disp, int atom_offset) {
  int i,j;
  int Ncarts = 3 * natom;
  int Nints = intcos.size();
  double **G_inv, *new_q, dx_max, dx_rms, dq_rms, first_dq_rms;
  double dx_rms_last = -1;

  fprintf(outfile,"\n\tBack-transformation to cartesian coordinates...\n");
  if (Opt_params.print_lvl >= 2) {
    fprintf(outfile,"\t---------------------------------------------------\n");
    fprintf(outfile,"\t Iter        RMS(dx)        Max(dx)        RMS(dq) \n");
    fprintf(outfile,"\t---------------------------------------------------\n");
  }

  fix_tors_near_180(); // subsequent computes will modify torsional values
  double * q_orig = intco_values();
  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile,"\nOriginal q internal coordinates\n");
    for (i=0; i<Nints; ++i) fprintf(outfile, "\t%15.10lf\n", q_orig[i]);
  }

  double * q_target = init_array(Nints);
  for (i=0; i<Nints; ++i)
    q_target[i] = q_orig[i] + dq[i];
  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile,"\nTarget q internal coordinates\n");
    for (i=0; i<Nints; ++i) fprintf(outfile, "\t%15.10lf\n", q_target[i]);
  }

  double * new_geom   = g_geom_array(); // cart geometry to start each iter
  double * first_geom = init_array(Ncarts); // first try at back-transformation
  double * dx = init_array(Ncarts);
  double * tmp_v_Nints = init_array(Nints);
  double **B = init_matrix(Nints, Ncarts);
  double **G = init_matrix(Nints, Nints);

  bool bt_iter_done = false;
  bool bt_converged = true;
  int bmat_iter_cnt = 0;

  while ( !bt_iter_done ) {

    // G = BuBt ; u = unit matrix; G = BBt
    compute_B(B);
    opt_matrix_mult(B, 0, B, 1, G, 0, Nints, Ncarts, Nints, 0);
    //compute_G(G, true); experimenting with mass-weighting

    // u B^t (G_inv dq) = dx
    G_inv = symm_matrix_inv(G, Nints, true);
    opt_matrix_mult(G_inv, 0, &dq, 1, &tmp_v_Nints, 1, Nints, Nints, 1, 0);
    opt_matrix_mult(B, 1, &tmp_v_Nints, 1, &dx, 1, Ncarts, Nints, 1, 0);

    // experimenting with mass-weighting
    //for (i=0; i<Ncarts; ++i)
    //  dx[i] /= mass[i/3];

    free_matrix(G_inv);

    for (i=0; i<Ncarts; ++i)
      new_geom[i] += dx[i];

    // Test for convergence of iterations
    dx_rms = array_rms(dx, Ncarts);
    dx_max = array_abs_max(dx, Ncarts);

    // maximum change and rms change in xyz coordinates < 10^-6
    if ( dx_rms < Opt_params.bt_dx_conv && dx_max < Opt_params.bt_dx_conv)
      bt_iter_done = true;
    else if (fabs(dx_rms - dx_rms_last) < Opt_params.bt_dx_conv_rms_change)
      bt_iter_done = true;                // change in rms change < 10^-12
    else if ( bmat_iter_cnt >= Opt_params.bt_max_iter) {
      bt_iter_done = true;
      bt_converged = false;
    }

    dx_rms_last = dx_rms;

    //compute_intco_values();
    set_geom_array(new_geom);
    new_q = intco_values();
    for (i=0; i< Nints; ++i)
      dq[i] = q_target[i] - new_q[i]; // calculate dq from _target_
    free_array(new_q);
    //fprintf(outfile,"dq from target\n");
    //print_array(outfile, dq, Nints);

    // save first try in case doesn't converge
    if (bmat_iter_cnt == 0) {
      for (i=0; i<Ncarts; ++i)
        first_geom[i] = new_geom[i];
      first_dq_rms = array_rms(dq, Nints);
    }
    dq_rms = array_rms(dq, Nints);

    if (Opt_params.print_lvl >= 2)
      fprintf (outfile,"\t%5d %14.1e %14.1e %14.1e\n", bmat_iter_cnt+1, dx_rms, dx_max, dq_rms);

    ++bmat_iter_cnt;
  }

  if (Opt_params.print_lvl >= 2)
    fprintf(outfile,"\t---------------------------------------------\n");

  if (bt_converged) {
    fprintf(outfile, "\tSuccessfully converged to displaced geometry.\n");
    if (dq_rms > first_dq_rms) {
      fprintf(outfile,"\tFirst geometry is closer to target in internal coordinates, so am using that one.\n");
      set_geom_array(first_geom);
    }
  }
  else {
    fprintf(outfile,"\tCould not converge backtransformation.\n");
    fprintf(outfile,"\tUsing first guess instead.\n");
    set_geom_array(first_geom);
  }

  // Set dq to final, total displacement ACHIEVED
  new_q = intco_values();
  for (i=0; i<Nints; ++i)
    dq[i] = new_q[i] - q_orig[i]; // calculate dq from _target_

  for (i=0; i<Nints; ++i) {
    if (intcos[i]->g_type() == tors_type) { // passed through 180
      if (dq[i] > _pi)
        dq[i] = dq[i] - (2 * _pi);
      else if (dq[i] < (-2 * _pi)) 
        dq[i] = dq[i] + (2 * _pi);
    }
  }

  if (print_disp) {

    double *f_q = p_Opt_data->g_forces_pointer(); // for printing

    fprintf(outfile,"\n\t---Internal Coordinate Step in ANG or DEG, aJ/ANG or AJ/DEG ---\n");
    fprintf(outfile,  "\t ----------------------------------------------------------------------\n");
    fprintf(outfile,  "\t Coordinate             Previous        Force       Change         New \n");
    fprintf(outfile,  "\t ----------             --------       ------       ------       ------\n");
    for (i=0; i<intcos.size(); ++i)
      intcos.at(i)->print_disp(outfile, q_orig[i], f_q[i], dq[i], new_q[i], atom_offset);
    fprintf(outfile,  "\t ----------------------------------------------------------------------\n");
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

  fflush(outfile);
  return;
}

}

