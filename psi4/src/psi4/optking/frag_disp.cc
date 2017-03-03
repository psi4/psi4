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

/*! \file    frag_disp.cc
    \ingroup optking
    \brief   displace fragment geometry only dq changes to values of coordinates
*/

#include "frag.h"
#include "linear_algebra.h"
#include "opt_data.h"

#include "print.h"
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
  int Nints = Ncoord();
  // If can't converge back-transformation, then reduce step size as necessary.
  bool ensure_convergence = Opt_params.ensure_bt_convergence;

  fix_tors_near_180(); // subsequent computes will modify torsional values for phase
  fix_oofp_near_180(); // subsequent computes will modify torsional values for phase
  double * q_orig = coord_values();

  // Do your best to backtransform all internal coordinate displacments.
  oprintf_out("\n\tBack-transformation to cartesian coordinates...\n");

  if (ensure_convergence) { // change step as necessary to get convergence
    double *dq_orig = init_array(Nints);
    array_copy(dq, dq_orig, Nints);

    double *orig_geom = g_geom_array();
    bool conv = false;
    int cnt = -1;

    while ( !conv ) {
      ++cnt;

      if (cnt > 0) {
        oprintf_out("Reducing step-size by a factor of %d.\n", 2*cnt);
        for (int i=0; i<Nints; ++i)
          dq[i] = dq_orig[i] / (2*cnt);
      }

      conv = displace_util(dq, false);

      if (!conv) {
        if (cnt == 5) {
          oprintf_out("\tUnable to back-transform even 1/10th of the desired step rigorously.\n");
          break; // we'll stick with the unconverged best guess from the smallest step try
        }
        else
          set_geom_array(orig_geom); // put original geometry back for next try
      }
    }

    if (conv && cnt > 0) {  // We were able to take a modest step.  Try to complete it.
      oprintf_out("%d partial back-transformations left to do.\n", 2*cnt-1);

      for (int j=1; j<2*cnt; ++j) {
        oprintf_out("Mini-step %d of %d.\n", j+1, 2*cnt);

        for (int i=0; i<Nints; ++i)
          dq[i] = dq_orig[i] / (2*cnt);

      // Project to get physically possible step again.  Doesn't seem to help.
      /*
      int Ncart = 3*g_natom();
        oprintf_out("\tProjecting next mini-step to physically possible one.\n");
        double **B = compute_B();
        double **G = init_matrix(Ncart, Ncart);
        opt_matrix_mult(B, 1, B, 0, G, 0, Ncart, Nints, Ncart, 0);
    
        double **G_inv = symm_matrix_inv(G, Ncart, 1);
        free_matrix(G);

        double **B_inv = init_matrix(Ncart, Nints);
        opt_matrix_mult(G_inv, 0, B, 1, B_inv, 0, Ncart, Ncart, Nints, 0);
        free_matrix(G_inv);

        double **P = init_matrix(Nints, Nints);
        opt_matrix_mult(B, 0, B_inv, 0, P, 0, Nints, Ncart, Nints, 0);
        free_matrix(B);
        free_matrix(B_inv);

        double * temp_arr = init_array(Nints);
        opt_matrix_mult(P, 0, &dq, 1, &temp_arr, 1, Nints, Nints, 1, 0);
        array_copy(temp_arr, dq, Nints);
        free_array(temp_arr);
        free_matrix(P);
      */

        // save cartesian geometry and put it back in if back-transformation fails
        double *g = g_geom_array();
        array_copy(g, orig_geom, 3*g_natom());
        free_array(g);

        fix_bend_axes();
        conv = displace_util(dq, false);
        unfix_bend_axes();
        if (!conv) {
          oprintf_out("\tCouldn't converge this mini-step, so quitting with previous geometry.\n");
          set_geom_array(orig_geom);
          break;
        }
      }
    }
    free_array(orig_geom);
  }
  else { // try to back-transform, but continue either way
    fix_bend_axes();
    displace_util(dq, false);
    unfix_bend_axes();
  }

  /* Algorithms that compute DQ, and the backtransformation above may
     allow frozen coordinates to drift a bit.  Fix them here. */
  // See if there are any frozen coordinates.
  bool frag_constraints_present = false;
  for (int i=0; i<Nints; ++i)
    if (coords.simples[i]->is_frozen())
      frag_constraints_present = true;

  if (frag_constraints_present) {
    double *q_before_adjustment = coord_values();

    double *dq_adjust_frozen = init_array(Nints);
    for (int i=0; i<Nints; ++i)
      if (coords.simples[i]->is_frozen())
        dq_adjust_frozen[i] = q_orig[i] - q_before_adjustment[i];

    oprintf_out("\n\tBack-transformation to cartesian coordinates to adjust frozen coordinates...\n");
    fix_bend_axes();
    displace_util(dq_adjust_frozen, true);
    unfix_bend_axes();

    free_array(q_before_adjustment);
    free_array(dq_adjust_frozen);
  }

  // Set dq to final, total displacement ACHIEVED
  double *q_final = coord_values();
  for (int i=0; i<Nints; ++i)
    dq[i] = q_final[i] - q_orig[i]; // calculate dq from _target_

  for (int i=0; i<Nints; ++i) {
    // passed through pi, but don't think this code is necessary; given the way values are computed
    if (coords.simples[i]->g_type() == tors_type ||
        coords.simples[i]->g_type() == oofp_type) {
      if (dq[i] > _pi) {
        dq[i] = dq[i] - (2 * _pi);
        oprintf_out("\tTorsion changed more than pi.  Fixing in displace().\n");
      }
      else if (dq[i] < (-2 * _pi)) {
        dq[i] = dq[i] + (2 * _pi);
        oprintf_out("\tTorsion changed more than pi.  Fixing in displace().\n");
      }
    }
  }

  oprintf_out("\n\t--- Internal Coordinate Step in ANG or DEG, aJ/ANG or AJ/DEG ---\n");
  oprintf_out(  "\t ---------------------------------------------------------------------------\n");
  oprintf_out(  "\t   Coordinate                Previous        Force       Change         New \n");
  oprintf_out(  "\t   ----------                --------       ------       ------       ------\n");
  for (int i=0; i<Ncoord(); ++i) {
    oprintf_out("\t %4d ",i+1);
    coords.print_disp(psi_outfile, qc_outfile, i, q_orig[i], fq[i], dq[i], q_final[i], atom_offset);
  }
  oprintf_out(  "\t ---------------------------------------------------------------------------\n");

  free_array(q_orig);
  free_array(q_final);
}

bool FRAG::displace_util(double *dq, bool focus_on_constraints) {
  int i;
  int Ncarts = 3 * natom;
  int Nints = Ncoord();
  double **G_inv, *new_q, dx_max, dx_rms, dq_rms, first_dq_rms;
  double dx_rms_last = -1;
  bool rval;

  double bt_dx_conv            = Opt_params.bt_dx_conv;
  double bt_dx_conv_rms_change = Opt_params.bt_dx_conv_rms_change;
  double bt_max_iter           = Opt_params.bt_max_iter;
  if (focus_on_constraints) {
    bt_dx_conv            = 1.0e-12;
    bt_dx_conv_rms_change = 1.0e-12;
    bt_max_iter           = 100;
  }

  double * q_orig = coord_values();

  double * q_target = init_array(Nints);
  for (i=0; i<Nints; ++i)
    q_target[i] = q_orig[i] + dq[i];


  if (Opt_params.print_lvl >= 2) {
    oprintf_out("\t In displace_util \n");
    oprintf_out("\t       Original         Target           Dq\n");
    for (i=0; i<Nints; ++i)
      oprintf_out("\t%15.10lf%15.10lf%15.10lf\n", q_orig[i], q_target[i], dq[i]);
  }

  if (Opt_params.print_lvl >= 2) {
    oprintf_out("\t---------------------------------------------------\n");
    oprintf_out("\t Iter        RMS(dx)        Max(dx)        RMS(dq) \n");
    oprintf_out("\t---------------------------------------------------\n");
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

/*
    // B dx    = dq
    // Bt B dx = Bt dq
    // dx      = (B^t B)^-1 B^t dq
    // dx      = G^-1 B^t dq, where G = B^t B.
    // Tried in 2014.  Will it give different results than the code below if there are redundancies?
    // In this form, G is cart x cart, instead of int by int.
    // Disadvantage is that G has always rotation and translations in it (i.e., 0 evals when diagonalized).
    double * tmp_v_Ncarts = init_array(Ncarts);
    double **Gx_inv;
    double **Gx = init_matrix(Ncarts, Ncarts);
    compute_B(B);
    opt_matrix_mult(B, 1, B, 0, Gx, 0, Ncarts, Nints, Ncarts, 0);
    Gx_inv = symm_matrix_inv(Gx, Ncarts, true);
    opt_matrix_mult(B, 1, &dq, 1, &tmp_v_Ncarts, 1, Ncarts, Nints, 1, 0);
    opt_matrix_mult(Gx_inv, 0, &tmp_v_Ncarts, 1, &dx, 1, Ncarts, Ncarts, 1, 0);
    free_matrix(Gx_inv);
*/

    // B dx = dq
    // B dx = (B Bt)(B Bt)^-1 dq
    // B dx = B * (Bt (B Bt)^-1) dq
    //   dx = Bt (B Bt)^-1 dq
    //   dx = Bt G^-1 dq, where G = B B^t.
    compute_B(B,0,0);
    opt_matrix_mult(B, 0, B, 1, G, 0, Nints, Ncarts, Nints, 0);

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
    else if ( dx_rms > 100.0 ) { // Give up.
      bt_iter_done = true;
      bt_converged = false;
    }
    dx_rms_last = dx_rms;

    set_geom_array(new_geom);
    new_q = coord_values();
    //oprintf("%d new_x[0]=%15.10lf  new_q[0]=%15.10lf\n", bmat_iter_cnt, new_geom[0], new_q[0]);

    if (focus_on_constraints) {
      // We allow the non-constrained coordinates to change slightly, to allow
      // the frozen ones to be converged tightly.  So pretend the others are ok.
      for (i=0; i< Nints; ++i)
        if (!coords.simples[i]->is_frozen())
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
      oprintf_out("\t%5d %14.1e %14.1e %14.1e\n", bmat_iter_cnt+1, dx_rms, dx_max, dq_rms);

    ++bmat_iter_cnt;
  }

  if (Opt_params.print_lvl >= 2)
    oprintf_out("\t---------------------------------------------------\n");

  if (Opt_params.print_lvl >= 2) {
    oprintf_out("\n\tReport of back-transformation:\n");
    oprintf_out("\t  int       q_target          Error\n");
    oprintf_out("\t-----------------------------------\n");
    for (i=0; i<Nints; ++i)
      oprintf_out("\t%5d%15.10lf%15.10lf\n", i+1, q_target[i], -dq[i]);
    oprintf_out("\n");
  }

  if (bt_converged) {
    oprintf_out("\tSuccessfully converged to displaced geometry.\n");
    rval = true;
    if (dq_rms > first_dq_rms) {
      oprintf_out("\tFirst geometry is closer to target in internal coordinates, so am using that one.\n");
      oprintf_out("\tFirst geometry has RMS(Delta(q)) = %8.2e\n", first_dq_rms);
      set_geom_array(first_geom);
      rval = false;
    }
  }
  else if (!focus_on_constraints) { // if we are fixing constraints, we'll keep the best we got.
    rval = false;
    oprintf_out("\tCould not converge backtransformation.\n");
    oprintf_out("\tUsing first guess instead.\n");
    if (Opt_params.opt_type == OPT_PARAMS::IRC)
      throw(INTCO_EXCEPT("Could not take constrained step in an IRC computation."));
    set_geom_array(first_geom);
  }
  else rval = true; // not converged and only for constraint fixing

  free_matrix(G);
  free_array(new_geom);
  free_array(first_geom);
  free_array(dx);
  free_array(tmp_v_Nints);
  free_matrix(B);

  free_array(q_target);
  free_array(q_orig);

  return rval;
}

}
