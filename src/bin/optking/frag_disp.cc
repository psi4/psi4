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

/*! \file    frag_disp.cc
    \ingroup optking
    \brief   displace fragment geometry only dq changes to values of coordinates
*/

#include "frag.h"
#include "linear_algebra.h"
#include "opt_data.h"

#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
#endif

namespace opt {

// dq - displacements in intrafragment internal coordinates to be performed; overridden
//      to actual displacements performed
// fq - internal coordinate forces (used for printing)
// atom_offset - number within the molecule of first atom in this fragment (used for printing)
void FRAG::displace(double *dq, double *fq, int atom_offset) {
  int Nints = intcos.size();

  fix_tors_near_180(); // subsequent computes will modify torsional values for phase
  double * q_orig = intco_values();

  // Do your best to backtransform all internal coordinate displacments.
  fprintf(outfile,"\n\tBack-transformation to cartesian coordinates...\n");
  displace_util(dq, false);

  /* Algorithms that compute DQ, and the backtransformation above may
     allow frozen coordinates to drift a bit.  Fix them here. */
  // See if there are any frozen coordinates.
  bool frag_constraints_present = false;
  for (int i=0; i<Nints; ++i)
    if (intcos[i]->is_frozen())
      frag_constraints_present = true;

  if (frag_constraints_present) {
    double *q_before_adjustment = intco_values();

    double *dq_adjust_frozen = init_array(Nints);
    for (int i=0; i<Nints; ++i)
      if (intcos[i]->is_frozen())
        dq_adjust_frozen[i] = q_orig[i] - q_before_adjustment[i];

    fprintf(outfile,"\n\tBack-transformation to cartesian coordinates to adjust frozen coordinates...\n");
    displace_util(dq_adjust_frozen, true);

    free_array(q_before_adjustment);
    free_array(dq_adjust_frozen);
  }

  // Set dq to final, total displacement ACHIEVED
  double *q_final = intco_values();
  for (int i=0; i<Nints; ++i)
    dq[i] = q_final[i] - q_orig[i]; // calculate dq from _target_

  for (int i=0; i<Nints; ++i) {
    if (intcos[i]->g_type() == tors_type) { // passed through 180, but don't think this code is necessary
      if (dq[i] > _pi)
        dq[i] = dq[i] - (2 * _pi);
      else if (dq[i] < (-2 * _pi))
        dq[i] = dq[i] + (2 * _pi);
    }
  }

  fprintf(outfile,"\n\t---Internal Coordinate Step in ANG or DEG, aJ/ANG or AJ/DEG ---\n");
  fprintf(outfile,  "\t ----------------------------------------------------------------------\n");
  fprintf(outfile,  "\t Coordinate             Previous        Force       Change         New \n");
  fprintf(outfile,  "\t ----------             --------       ------       ------       ------\n");
  for (int i=0; i<intcos.size(); ++i)
    intcos.at(i)->print_disp(outfile, q_orig[i], fq[i], dq[i], q_final[i], atom_offset);
  fprintf(outfile,  "\t ----------------------------------------------------------------------\n");

  free_array(q_orig);
  free_array(q_final);
}

void FRAG::displace_util(double *dq, bool focus_on_constraints) {
  int i,j;
  int Ncarts = 3 * natom;
  int Nints = intcos.size();
  double **G_inv, *new_q, dx_max, dx_rms, dq_rms, first_dq_rms;
  double dx_rms_last = -1;

  double bt_dx_conv            = Opt_params.bt_dx_conv;
  double bt_dx_conv_rms_change = Opt_params.bt_dx_conv_rms_change;
  double bt_max_iter           = Opt_params.bt_max_iter;
  if (focus_on_constraints) {
    bt_dx_conv            = 1.0e-12;
    bt_dx_conv_rms_change = 1.0e-12;
    bt_max_iter           = 100;
  }

  if (Opt_params.print_lvl >= 2) {
    fprintf(outfile,"\t---------------------------------------------------\n");
    fprintf(outfile,"\t Iter        RMS(dx)        Max(dx)        RMS(dq) \n");
    fprintf(outfile,"\t---------------------------------------------------\n");
  }

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
    free_matrix(G_inv);

    for (i=0; i<Ncarts; ++i)
      new_geom[i] += dx[i];

    // Test for convergence of iterations
    dx_rms = array_rms(dx, Ncarts);
    dx_max = array_abs_max(dx, Ncarts);

    // maximum change and rms change in xyz coordinates < 10^-6
    if ( dx_rms < bt_dx_conv && dx_max < bt_dx_conv)
      bt_iter_done = true;
    else if (fabs(dx_rms - dx_rms_last) < bt_dx_conv_rms_change)
      bt_iter_done = true;
    else if ( bmat_iter_cnt >= bt_max_iter) {
      bt_iter_done = true;
      bt_converged = false;
    }
    dx_rms_last = dx_rms;

    set_geom_array(new_geom);
    new_q = intco_values();

    if (focus_on_constraints) {
      // We allow the non-constrained coordinates to change slightly, to allow
      // the frozen ones to be converged tightly.  So pretend the others are ok.
      for (i=0; i< Nints; ++i)
        if (!intcos[i]->is_frozen())
          q_target[i] = new_q[i];
    }

    for (i=0; i< Nints; ++i)
      dq[i] = q_target[i] - new_q[i];
    free_array(new_q);

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
  else if (!focus_on_constraints) { // if we are fixing constraints, we'll keep the best we got.
    fprintf(outfile,"\tCould not converge backtransformation.\n");
    fprintf(outfile,"\tUsing first guess instead.\n");
    if (Opt_params.opt_type == OPT_PARAMS::IRC)
      throw(INTCO_EXCEPT("Could not take constrained step in an IRC computation."));
    set_geom_array(first_geom);
  }

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

