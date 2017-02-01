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

/*! \file    opt_data.cc
    \ingroup optking
    \brief   OPT_DATA associated functions that do not do i/o
*/

#include "opt_data.h"

#include "print.h"
#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include <qchem.h> // used for rem_write
 #include "qcmath.h"
#endif

namespace opt {

//minimal constructor - just allocates memory
STEP_DATA::STEP_DATA(int Nintco, int Ncart) {
  f_q  = init_array(Nintco);
  geom = init_array(Ncart);
  energy = 0.0;
  DE_predicted = 0.0;
  unit_step = init_array(Nintco);
  dq_norm = 0.0;
  dq_gradient = 0.0;
  dq_hessian = 0.0;
  dq = init_array(Nintco);
}

//free memory
STEP_DATA::~STEP_DATA() {
  free_array(f_q);
  free_array(geom);
  free_array(unit_step);
  free_array(dq);
}

// save geometry and energy
void STEP_DATA::save_geom_energy(double *geom_in, double energy_in, int Ncart) {
  array_copy(geom_in, geom, Ncart);
  energy = energy_in;
}

// save rest of stuff
void STEP_DATA::save_step_info(double DE_predicted_in, double *unit_step_in, double dq_norm_in,
    double dq_gradient_in, double dq_hessian_in, int Nintco) {
  DE_predicted = DE_predicted_in;
  array_copy(unit_step_in, unit_step, Nintco);
  dq_norm = dq_norm_in;
  dq_gradient = dq_gradient_in;
  dq_hessian = dq_hessian_in;
}

OPT_DATA::~OPT_DATA() {
  free_matrix(H);
  free_array(rfo_eigenvector);
  for (std::size_t i=0; i<steps.size(); ++i)
    delete steps[i];
  steps.clear();
}

// Check convergence criteria and print status to output file
// return true, if geometry is optimized
// Checks maximum force && (Delta(E) || maximum disp)
bool OPT_DATA::conv_check(opt::MOLECULE &mol) const {
  double *dq =  g_dq_pointer();
  double max_disp = array_abs_max(dq, Nintco);
  double rms_disp = array_rms(dq, Nintco);

  double DE;
  if (g_iteration() > 1)
    DE = g_energy() - g_last_energy();
  else
    DE = g_energy();

  double *f =  g_forces_pointer();

  if (Opt_params.opt_type == OPT_PARAMS::IRC)
    if (!p_irc_data->go) { return 1; }

  // Save forces and put back in below.
  double *f_backup;
  if (Opt_params.opt_type == OPT_PARAMS::IRC || mol.has_fixed_eq_vals()) {
    f_backup = init_array(Nintco);
    array_copy(f, f_backup, Nintco);
  }

  // for IRC only consider forces tangent to the hypersphere search surface
  if (Opt_params.opt_type == OPT_PARAMS::IRC) {
    double **G = mol.compute_G(1);
    double **Ginv = symm_matrix_inv(G, Nintco, Nintco); 
    free_matrix(G);

    // compute p_m, mass-weighted hypersphere vector
    const double *q_pivot = p_irc_data->g_q_pivot();
    double *q = mol.coord_values();
    double *p = init_array(Nintco);
    for (int i=0; i<Nintco; ++i)
      p[i] = q[i] - q_pivot[i];
    free_array(q);

    // gradient perpendicular to p and tangent to hypersphere is:
    // g_m' = g_m - (g_m^t p_m / p_m^t p_m) p_m, or
    // g'   = g   - (g^t p / (p^t G^-1 p)) G^-1 p
    double *Ginv_p = init_array(Nintco);
    for(int i=0; i<Nintco; i++)
      Ginv_p[i] = array_dot(Ginv[i], p, Nintco);
    free_matrix(Ginv);

    double overlap = array_dot(f, p, Nintco) / array_dot(p, Ginv_p, Nintco);

    for (int i=0; i<Nintco; ++i)
      f[i] -= overlap * Ginv_p[i];
    free_array(Ginv_p);

    if (Opt_params.print_lvl >= 2) {
      oprintf_out("\tForces perpendicular to hypersphere.\n");
      oprint_array_out(f, Nintco);
    }

  }

  // Remove arbitrary forces for user-specified equilibrium values. 
  if (mol.has_fixed_eq_vals()) {
    array_copy(f, f_backup, Nintco); // save forces and put back in below
    oprintf_out("\t Forces used to impose fixed constraints are not included.\n");
    oprintf_out("\t  Forces zeroed: ");

    for (int i=0; i<mol.Ncoord(); ++i) {
      if (mol.is_coord_fixed(i)) {
        oprintf_out("%d ", i+1);
        f[i] = 0.0;
      }
    }
    oprintf_out("\n");
  }

  double max_force = array_abs_max(f, Nintco);
  double rms_force = array_rms(f, Nintco);

  if (Opt_params.opt_type != OPT_PARAMS::IRC) {
      oprintf_out( "\n  ==> Convergence Check <==\n\n");
      oprintf_out( "  Measures of convergence in internal coordinates in au.\n");
      oprintf_out( "  Criteria marked as inactive (o), active & met (*), and active & unmet ( ).\n");
  
      oprintf_out( "  ---------------------------------------------------------------------------------------------");
      if (g_iteration() == 1) oprintf_out( " ~");
      oprintf_out( "\n");
  
      oprintf_out( "   Step     Total Energy     Delta E     MAX Force     RMS Force      MAX Disp      RMS Disp   ");
      if (g_iteration() == 1) oprintf_out( " ~");
      oprintf_out( "\n");
  
      oprintf_out( "  ---------------------------------------------------------------------------------------------");
      if (g_iteration() == 1) oprintf_out( " ~");
      oprintf_out( "\n");
  
      oprintf_out( "    Convergence Criteria");
      Opt_params.i_max_DE    ? oprintf_out( "  %10.2e %1s", Opt_params.conv_max_DE,    "*") : oprintf_out("             %1s", "o");
      Opt_params.i_max_force ? oprintf_out( "  %10.2e %1s", Opt_params.conv_max_force, "*") : oprintf_out("             %1s", "o");
      Opt_params.i_rms_force ? oprintf_out( "  %10.2e %1s", Opt_params.conv_rms_force, "*") : oprintf_out("             %1s", "o");
      Opt_params.i_max_disp  ? oprintf_out( "  %10.2e %1s", Opt_params.conv_max_disp,  "*") : oprintf_out("             %1s", "o");
      Opt_params.i_rms_disp  ? oprintf_out( "  %10.2e %1s", Opt_params.conv_rms_disp,  "*") : oprintf_out("             %1s", "o");
  
      if (g_iteration() == 1) oprintf_out( "  ~");
      oprintf_out( "\n");
  
      oprintf_out( "  ---------------------------------------------------------------------------------------------");
      if (g_iteration() == 1) oprintf_out( " ~");
      oprintf_out( "\n");
  
      oprintf_out( "   %4d %16.8f  %10.2e %1s  %10.2e %1s  %10.2e %1s  %10.2e %1s  %10.2e %1s  ~\n", iteration, g_energy(),
        DE, (Opt_params.i_max_DE ? ((fabs(DE) < Opt_params.conv_max_DE) ? "*" : "") : "o"), 
        max_force, (Opt_params.i_max_force ? ((fabs(max_force) < Opt_params.conv_max_force) ? "*" : "") : "o"),
        rms_force, (Opt_params.i_rms_force ? ((fabs(rms_force) < Opt_params.conv_rms_force) ? "*" : "") : "o"),
        max_disp, (Opt_params.i_max_disp ? ((fabs(max_disp) < Opt_params.conv_max_disp) ? "*" : "") : "o"),
        rms_disp, (Opt_params.i_rms_disp ? ((fabs(rms_disp) < Opt_params.conv_rms_disp) ? "*" : "") : "o"));
  
      oprintf_out( "  ---------------------------------------------------------------------------------------------\n\n");
  }

  // Return forces to what they were when conv_check was called
  if (Opt_params.opt_type == OPT_PARAMS::IRC || mol.has_fixed_eq_vals()) {
    array_copy(f_backup, f, Nintco);
    free_array(f_backup);
  }

// Test for convergence
#if defined(OPTKING_PACKAGE_PSI)

  // Q-Chem and Gaussian have convergence tests involving interplay among criteria.
  //   The requirement of i_untampered means that if a user explicitly adds any of the
  //   5 indiv. criteria on top of G_CONVERGENCE, it is required to be met.
  if ( 
         (  // For un-modified Q-Chem or Molpro (Baker) criteria, if forces and either energy change or displacement met, convergence!
               (Opt_params.i_untampered) 
            && ( (Opt_params.general_conv == "QCHEM") || (Opt_params.general_conv == "MOLPRO") )
            && ( (max_force < Opt_params.conv_max_force) && ((fabs(DE) < Opt_params.conv_max_DE) || (max_disp < Opt_params.conv_max_disp)) )
         )
      || (  // For un-modified Gaussian criteria, if max/rms forces/disp met or flat potential forces met, convergence!
               (Opt_params.i_untampered) 
            && ( (Opt_params.general_conv == "GAU") || (Opt_params.general_conv == "GAU_TIGHT") ||
                 (Opt_params.general_conv == "GAU_VERYTIGHT") || (Opt_params.general_conv == "GAU_LOOSE") )
            && (    ( (max_force < Opt_params.conv_max_force) && (rms_force < Opt_params.conv_rms_force) &&
                      (max_disp  < Opt_params.conv_max_disp)  && (rms_disp  < Opt_params.conv_rms_disp) )
                 || (rms_force * 100 < Opt_params.conv_rms_force) )
         )
      || (  // Otherwise, for all criteria, if criterion not active or criterion met, convergence!
            (!(Opt_params.i_max_DE) || (fabs(DE) < Opt_params.conv_max_DE)) &&
            (!(Opt_params.i_max_force) || (max_force < Opt_params.conv_max_force)) &&
            (!(Opt_params.i_rms_force) || (rms_force < Opt_params.conv_rms_force)) &&
            (!(Opt_params.i_max_disp) || (max_disp < Opt_params.conv_max_disp)) &&
            (!(Opt_params.i_rms_disp) || (rms_disp < Opt_params.conv_rms_disp))
         )
      || (
            (Opt_params.opt_type == OPT_PARAMS::IRC) &&
            (!(Opt_params.i_max_DE) || (fabs(DE) < Opt_params.conv_max_DE)) &&
            (!(Opt_params.i_max_disp) || (max_disp < Opt_params.conv_max_disp)) &&
            (!(Opt_params.i_rms_disp) || (rms_disp < Opt_params.conv_rms_disp))
         )
     ) {

    // This environment variable will store the number of iterations required 
    // for convergence; it allows the db() python utilities to collect this information
    psi::Process::environment.globals["OPTIMIZATION ITERATIONS"] = g_iteration();
    return true; // structure is optimized!
  }
  else 
    return false;

#elif defined(OPTKING_PACKAGE_QCHEM)
  // convergence test is forces and either energy change or displacement
  if ((max_force < Opt_params.conv_max_force) &&
      ((fabs(DE) < Opt_params.conv_max_DE) || (max_disp < Opt_params.conv_max_disp)))  {
    return true; // structure is optimized!
  }
  else 
    return false;
#endif
}

void OPT_DATA::summary(void) const {
  double DE, *f, *dq, max_force, max_disp, rms_force, rms_disp;

  oprintf_out( "\n  ==> Optimization Summary <==\n\n");
  oprintf_out( "  Measures of convergence in internal coordinates in au.\n");

  oprintf_out( "  --------------------------------------------------------------------------------------------------------------- ~\n");
  oprintf_out( "   Step         Total Energy             Delta E       MAX Force       RMS Force        MAX Disp        RMS Disp  ~\n");
  oprintf_out( "  --------------------------------------------------------------------------------------------------------------- ~\n");
  for (int i=0; i<iteration; ++i) {

    if (i == 0) DE = g_energy(i);
    else DE = g_energy(i) - g_energy(i-1); 

    f = g_forces_pointer(i);
    max_force = array_abs_max(f, Nintco);
    rms_force = array_rms(f, Nintco);

    dq = g_dq_pointer(i);
    max_disp = array_abs_max(dq, Nintco);
    rms_disp = array_rms(dq, Nintco);

    oprintf_out( "   %4d %20.12lf  %18.12lf    %12.8lf    %12.8lf    %12.8lf    %12.8lf  ~\n", i + 1, g_energy(i),
      DE, max_force, rms_force, max_disp, rms_disp); 
  }
  oprintf_out( "  --------------------------------------------------------------------------------------------------------------- ~\n\n");
}

inline int sign_of_double(double d) {
  if (d>0) return 1;
  else if (d<0) return -1;
  else return 0; 
}

// do hessian update
void OPT_DATA::H_update(opt::MOLECULE & mol) {

  if (Opt_params.H_update == OPT_PARAMS::BFGS)
    oprintf_out("\n\tPerforming BFGS update.\n");
  else if (Opt_params.H_update == OPT_PARAMS::MS)
    oprintf_out("\n\tPerforming Murtagh/Sargent update.\n");
  else if (Opt_params.H_update == OPT_PARAMS::POWELL)
    oprintf_out("\n\tPerforming Powell update.\n");
  else if (Opt_params.H_update == OPT_PARAMS::BOFILL)
    oprintf_out("\n\tPerforming Bofill update.\n");
  else if (Opt_params.H_update == OPT_PARAMS::NONE) {
    oprintf_out("\n\tNo Hessian update performed.\n");
    return;
  }

  int step_this = steps.size()-1;

  /*** Read/compute current internals and forces ***/
  double *f, *x, *q;
  f = steps[step_this]->g_forces_pointer();
  x = steps[step_this]->g_geom_const_pointer();

  mol.set_geom_array(x);
  q = mol.coord_values();
  mol.fix_tors_near_180(); // Fix configurations of torsions.
  mol.fix_oofp_near_180();

  double *f_old, *x_old, *q_old, *dq, *dg;
  double gq, qq, qz, zz, phi;
  dq = init_array(Nintco);
  dg = init_array(Nintco);

  // Don't go further back than the last Hessian calculation 
  int check_start = step_this - steps_since_last_H;
  if (check_start < 0) check_start = 0;
  if (steps_since_last_H)
    oprintf_out("\tPrevious computed or guess Hessian on step %d.\n", check_start+1);

  // Make list of old geometries to update with.  Check each one to see if it is too close for stability.
  std::vector<int> use_steps;

  for (int i_step= ((int) steps.size()) - 2; i_step>=check_start; --i_step) {

    // Read/compute old internals and forces
    f_old = g_forces_pointer(i_step);
    x_old = g_geom_const_pointer(i_step);
  
    mol.set_geom_array(x_old);
    q_old = mol.coord_values();
  
    for (int i=0;i<Nintco;++i) {
      dq[i] = q[i] - q_old[i];
      dg[i] = (-1.0) * (f[i] - f_old[i]); // gradients -- not forces!
    }
  
    gq = array_dot(dq, dg, Nintco);
    qq = array_dot(dq, dq, Nintco);

    // If there is only one left, take it no matter what.
    if (use_steps.empty() && i_step == check_start) {
      use_steps.push_back(i_step);
      break;
    }
 
    if ( (fabs(gq) < Opt_params.H_update_den_tol) ||(fabs(qq) < Opt_params.H_update_den_tol) ) {
      oprintf_out("\tDenominators (dg)(dq) or (dq)(dq) are very small.\n");
      oprintf_out("\t Skipping Hessian update for step %d.\n", i_step+1);
      continue;
    }
  
    double max_change = array_abs_max(dq, Nintco);
    if ( max_change > Opt_params.H_update_dq_tol ) {
      oprintf_out("\tChange in internal coordinate of %5.2e exceeds limit of %5.2e.\n", max_change, Opt_params.H_update_dq_tol );
      oprintf_out("\t Skipping Hessian update for step %d.\n", i_step+1);
      continue;
    }
    use_steps.push_back(i_step);
    if ((int) use_steps.size() == Opt_params.H_update_use_last)
      break;
  }

  oprintf_out("\tSteps to be used in Hessian update:");
  for (std::size_t i=0; i<use_steps.size(); ++i)
    oprintf_out(" %d", use_steps[i]+1);
  oprintf_out("\n");

  double **H = g_H_pointer();
  double **H_new = init_matrix(Nintco, Nintco);

  for (std::size_t s=0; s<use_steps.size(); ++s) {

    int i_step = use_steps[s];

    // Read/compute old internals and forces
    f_old = g_forces_pointer(i_step);
    x_old = g_geom_const_pointer(i_step);

    mol.set_geom_array(x_old);
    q_old = mol.coord_values();

    for (int i=0;i<Nintco;++i) {
      // Turns out you don't have to correct for torsional changes through 180 in the forces.
      dq[i] = q[i] - q_old[i];
      dg[i] = (-1.0) * (f[i] - f_old[i]); // gradients -- not forces!
    }

    gq = array_dot(dq, dg, Nintco);
    qq = array_dot(dq, dq, Nintco);

    // See  J. M. Bofill, J. Comp. Chem., Vol. 15, pages 1-11 (1994)
    //  and Helgaker, JCP 2002 for formula.
    if (Opt_params.H_update == OPT_PARAMS::BFGS) {
      for (int i=0; i<Nintco; ++i)
        for (int j=0; j<Nintco; ++j)
          H_new[i][j] = H[i][j] + dg[i] * dg[j] / gq ;

      double *Hdq = init_array(Nintco);
      opt_matrix_mult(H, 0, &dq, 1, &Hdq, 1, Nintco, Nintco, 1, 0);

      double qHq = array_dot(dq, Hdq, Nintco);

      for (int i=0; i<Nintco; ++i)
        for (int j=0; j<Nintco; ++j)
          H_new[i][j] -=  Hdq[i] * Hdq[j] / qHq ;

      free_array(Hdq);
    }
    else if (Opt_params.H_update == OPT_PARAMS::MS) {
      double *Z = init_array(Nintco);
      opt_matrix_mult(H, 0, &dq, 1, &Z, 1, Nintco, Nintco, 1, 0);

      for (int i=0; i<Nintco; ++i)
        Z[i] = dg[i] - Z[i];

      qz = array_dot(dq, Z, Nintco);

      for (int i=0; i<Nintco; ++i)
        for (int j=0; j<Nintco; ++j)
          H_new[i][j] = H[i][j] + Z[i] * Z[j] / qz ;

      free_array(Z);
    }
    else if (Opt_params.H_update == OPT_PARAMS::POWELL) {
      double * Z = init_array(Nintco);
      opt_matrix_mult(H, 0, &dq, 1, &Z, 1, Nintco, Nintco, 1, 0);

      for (int i=0; i<Nintco; ++i)
        Z[i] = dg[i] - Z[i];

      qz = array_dot(dq, Z, Nintco);

      for (int i=0; i<Nintco; ++i)
        for (int j=0; j<Nintco; ++j)
          H_new[i][j] = H[i][j] - qz/(qq*qq)*dq[i]*dq[j] + (Z[i]*dq[j] + dq[i]*Z[j])/qq;

      free_array(Z);
    }
    else if (Opt_params.H_update == OPT_PARAMS::BOFILL) {
      // Bofill = (1-phi) * MS + phi * Powell
      double *Z = init_array(Nintco);
      opt_matrix_mult(H, 0, &dq, 1, &Z, 1, Nintco, Nintco, 1, 0);

      for (int i=0; i<Nintco; ++i)
        Z[i] = dg[i] - Z[i];

      qz = array_dot(dq, Z, Nintco);
      zz = array_dot(Z, Z, Nintco);

      phi = 1.0 - qz*qz/(qq*zz);
      if (phi < 0.0) phi = 0.0;
      if (phi > 1.0) phi = 1.0;

      for (int i=0; i<Nintco; ++i) // (1-phi)*MS
        for (int j=0; j<Nintco; ++j)
          H_new[i][j] = H[i][j] + (1.0-phi) * Z[i] * Z[j] / qz ;

      for (int i=0; i<Nintco; ++i)
        for (int j=0; j<Nintco; ++j) // (phi * Powell)
          H_new[i][j] += phi * (-1.0*qz/(qq*qq)*dq[i]*dq[j] + (Z[i]*dq[j] + dq[i]*Z[j])/qq);

      free_array(Z);
    } // end BOFILL

    if (Opt_params.H_update_limit) { // limit changes in H
      // Changes to the Hessian from the update scheme are limited to the larger of
      // (H_update_limit_scale)*(the previous value) and H_update_limit_max.
      double max_limit   = Opt_params.H_update_limit_max;
      double scale_limit = Opt_params.H_update_limit_scale;
      double max;
  
      // compute change in Hessian
      for (int i=0; i<Nintco; ++i)
        for (int j=0; j<Nintco; ++j)
          H_new[i][j] -= H[i][j];
  
      for (int i=0; i<Nintco; ++i) {
        for (int j=0; j<Nintco; ++j) {
          double val = fabs(scale_limit*H[i][j]);
          max = ((val > max_limit) ? val : max_limit);
  
        if (fabs(H_new[i][j]) < max)
          H[i][j] += H_new[i][j];
        else // limit change to max
          H[i][j] += max * sign_of_double(H_new[i][j]);
        }
      }
    }
    else { // copy H_new into H
      for (int i=0; i<Nintco; ++i)
        for (int j=0; j<Nintco; ++j)
          H[i][j] = H_new[i][j];
    }

    free_array(q_old);
    zero_matrix(H_new, Nintco, Nintco);
  } // end loop over steps to use in update
  free_array(q);
  free_array(dq);
  free_array(dg);

  // put original geometry back into molecule object
  mol.set_geom_array(x);


  free_matrix(H_new);
  if (Opt_params.print_lvl >= 2) {
    oprintf_out("\nUpdated Hessian (in au)\n");
    oprint_matrix_out(H, Nintco, Nintco);
  }
  return;
}

// read entry from binary file ; file pointer must be in right place for qchem code
void STEP_DATA::read(int istep, int Nintco, int Ncart) {
  char lbl[80];
  sprintf(lbl, "f_q %d", istep);
  opt_io_read_entry(lbl, (char *) f_q,    Nintco*sizeof(double));
  sprintf(lbl, "geom %d", istep);
  opt_io_read_entry(lbl, (char *) geom,    Ncart*sizeof(double));
  sprintf(lbl, "energy %d", istep);
  opt_io_read_entry(lbl, (char *) &energy,       sizeof(double));
  sprintf(lbl, "DE_predicted %d", istep);
  opt_io_read_entry(lbl, (char *) &DE_predicted, sizeof(double));
  sprintf(lbl, "unit_step %d", istep);
  opt_io_read_entry(lbl, (char *) unit_step, Nintco*sizeof(double));
  sprintf(lbl, "dq_norm %d", istep);
  opt_io_read_entry(lbl, (char *) &dq_norm,      sizeof(double));
  sprintf(lbl, "dq_gradient %d", istep);
  opt_io_read_entry(lbl, (char *) &dq_gradient,  sizeof(double));
  sprintf(lbl, "dq_hessian %d", istep);
  opt_io_read_entry(lbl, (char *) &dq_hessian,   sizeof(double));
  sprintf(lbl, "dq %d", istep);
  opt_io_read_entry(lbl, (char *) dq,       Nintco*sizeof(double));
}

// read entry from binary file ; file pointer must be in right place for qchem code
//write entry to binary file
void STEP_DATA::write(int istep, int Nintco, int Ncart) {
  char lbl[80];
  sprintf(lbl, "f_q %d", istep);
  opt_io_write_entry(lbl, (char *) f_q, Nintco*sizeof(double));
  sprintf(lbl, "geom %d", istep);
  opt_io_write_entry(lbl, (char *) geom, Ncart*sizeof(double));
  sprintf(lbl, "energy %d", istep);
  opt_io_write_entry(lbl, (char *) &energy, sizeof(double));
  sprintf(lbl, "DE_predicted %d", istep);
  opt_io_write_entry(lbl, (char *) &DE_predicted, sizeof(double));
  sprintf(lbl, "unit_step %d", istep);
  opt_io_write_entry(lbl, (char *) unit_step, Nintco*sizeof(double));
  sprintf(lbl, "dq_norm %d", istep);
  opt_io_write_entry(lbl, (char *) &dq_norm, sizeof(double));
  sprintf(lbl, "dq_gradient %d", istep);
  opt_io_write_entry(lbl, (char *) &dq_gradient, sizeof(double));
  sprintf(lbl, "dq_hessian %d", istep);
  opt_io_write_entry(lbl, (char *) &dq_hessian, sizeof(double));
  sprintf(lbl, "dq %d", istep);
  opt_io_write_entry(lbl, (char *) dq, Nintco*sizeof(double));
}

// constructor function reads available data and allocates memory for current step
OPT_DATA::OPT_DATA(int Nintco_in, int Ncart_in) {
  Nintco = Nintco_in;
  Ncart = Ncart_in;
  H = init_matrix(Nintco, Nintco);
  rfo_eigenvector = init_array(Nintco);

  bool data_file_present = opt_io_is_present(); // determine if old data file is present

  if (!data_file_present) {
    oprintf_out( "\tPrevious optimization step data not found.  Starting new optimization.\n\n");
    iteration = 0;
    steps_since_last_H = 0;
    consecutive_backsteps = 0;
  }
  else {
    int Nintco_old, Ncart_old;
    opt_io_open(opt::OPT_IO_OPEN_OLD);
    opt_io_read_entry("Nintco", (char *) &Nintco_old, sizeof(int));
    opt_io_read_entry("Ncart",  (char *) &Ncart_old, sizeof(int));

    if (Nintco_old != Nintco)
      oprintf_out( "\tThe number of coordinates has changed.  Ignoring old data.\n");
    if (Ncart_old != Ncart)
      oprintf_out( "\tThe number of atoms has changed.  Ignoring old data.\n");

    if ( (Nintco_old != Nintco) || (Ncart_old != Ncart) ) {
      iteration = 0;
      steps_since_last_H = 0;
      consecutive_backsteps = 0;
      opt_io_close(0); // close and delete
    }
    else { // read in old optimization data
      opt_io_read_entry("H", (char *) H[0], Nintco * Nintco * sizeof(double) );
      opt_io_read_entry("iteration", (char *) &iteration, sizeof(int));
      opt_io_read_entry("steps_since_last_H", (char *) &steps_since_last_H, sizeof(int));
      opt_io_read_entry("consecutive_backsteps", (char *) &consecutive_backsteps, sizeof(int));
      opt_io_read_entry("rfo_eigenvector", (char *) rfo_eigenvector, Nintco*sizeof(double));
      for (int i=0; i<iteration; ++i) {
        STEP_DATA *one_step = new STEP_DATA(Nintco, Ncart);
        one_step->read(i+1, Nintco, Ncart);
        steps.push_back(one_step);
      }
      opt_io_close(1);
    }
  }

  ++iteration; // increment for current step
  ++steps_since_last_H;
  // create memory for this, current step
  STEP_DATA *one_step = new STEP_DATA(Nintco, Ncart);
  steps.push_back(one_step);
}

// write data to binary file
void OPT_DATA::write(void) {
  opt_io_open(opt::OPT_IO_OPEN_OLD);

  oprintf_out("\tWriting optimization data to binary file.\n");
  opt_io_write_entry("Nintco", (char *) &Nintco, sizeof(int));
  opt_io_write_entry("Ncart" , (char *) &Ncart , sizeof(int));
  opt_io_write_entry("H", (char *) H[0], Nintco * Nintco * sizeof(double) );
  opt_io_write_entry("iteration", (char *) &iteration, sizeof(int));
  opt_io_write_entry("steps_since_last_H", (char *) &steps_since_last_H, sizeof(int));
  opt_io_write_entry("consecutive_backsteps", (char *) &consecutive_backsteps, sizeof(int));
  opt_io_write_entry("rfo_eigenvector", (char *) rfo_eigenvector, Nintco*sizeof(double));

  for (std::size_t i=0; i<steps.size(); ++i)
    steps[i]->write(i+1, Nintco, Ncart);
  opt_io_close(1);

  
  return;
}

// Report on performance of last step
// Eventually might have this function return false to reject a step
bool OPT_DATA::previous_step_report(void) const {
  oprintf_out(  "\tCurrent energy   : %20.10lf\n\n", p_Opt_data->g_energy());

  if (steps.size() == 1) {
    // optking changes the intrafragment_step_limit keyword value when it changes the trust
    // radius.  So on the first step, save the initial user-given step size so that it can be
    // reset to the user-defined value after one molecule is optimized.
    Opt_params.intrafragment_step_limit_orig = Opt_params.intrafragment_step_limit;
    return true;
  }
  bool dont_backup_yet = ((steps.size() < 5) ? true : false);

  oprintf_out("\tEnergy change for the previous step:\n");
  oprintf_out("\t\tProjected    : %20.10lf\n", p_Opt_data->g_last_DE_predicted());
  oprintf_out("\t\tActual       : %20.10lf\n",
      p_Opt_data->g_energy() - p_Opt_data->g_last_energy());

  double Energy_ratio = (p_Opt_data->g_energy() - p_Opt_data->g_last_energy()) / g_last_DE_predicted();
  double Energy_change = p_Opt_data->g_energy() - p_Opt_data->g_last_energy();

  if (Opt_params.print_lvl >= 2)
    oprintf_out("\tEnergy ratio = %10.5lf\n", Energy_ratio);

  // Minimum search
  if (Opt_params.opt_type == OPT_PARAMS::MIN) {

    // Predicted up. Actual down.  OK.  Do nothing.
    if (p_Opt_data->g_last_DE_predicted() > 0 && Energy_ratio < 0.0) {
        return true;
    }
    // Actual step is  up.
    else if (Energy_change > 0) {

      // Always throw exception for bad step.
      if (Opt_params.dynamic && !dont_backup_yet)
        throw(BAD_STEP_EXCEPT("Energy has increased in a minimization.\n"));
      // Not dynamic.  Do limited backsteps only upon request.  Otherwise, keep going.
      else if (consecutive_backsteps < Opt_params.consecutive_backsteps_allowed)
        throw(BAD_STEP_EXCEPT("Energy has increased in a minimization.\n"));
    }
    // Predicted down.  Actual down.
    else if (Energy_ratio < 0.25) {
      decrease_trust_radius();
    }
    else if (Energy_ratio > 0.75) {
      increase_trust_radius();
    }
  }

  return true;
}

// These functions are a bit out of place, given that the trust radius is not
// stored inside opt_data.
void OPT_DATA::increase_trust_radius(void) const {
  std::string module = "OPTKING";
  std::string key = "INTRAFRAG_STEP_LIMIT";
  double max = Opt_params.intrafragment_step_limit_max;
  if (Opt_params.intrafragment_step_limit != max) {
    //double new_val = Opt_params.intrafragment_step_limit * 2;
    double new_val = Opt_params.intrafragment_step_limit * 3;
    Opt_params.intrafragment_step_limit = ((new_val > max) ? max : new_val);
    oprintf_out("\tEnergy ratio indicates good step: Trust radius increased to %6.3e.\n\n",
        Opt_params.intrafragment_step_limit);

    // Save new trust radii in environment
#if defined(OPTKING_PACKAGE_PSI)
    psi::Process::environment.options.set_double(module, key, Opt_params.intrafragment_step_limit);
#elif defined(OPTKING_PACKAGE_QCHEM)
    int qchem_limit_val = (int) (Opt_params.intrafragment_step_limit * 1000);
    rem_write(qchem_limit_val, REM_GEOM_OPT2_INTRAFRAG_STEP_LIMIT);
#endif
  }
  return;
}

void OPT_DATA::decrease_trust_radius(void) const {
  std::string module = "OPTKING";
  std::string key = "INTRAFRAG_STEP_LIMIT";
  double min = Opt_params.intrafragment_step_limit_min;
  if (Opt_params.intrafragment_step_limit != min) {
    double new_val = Opt_params.intrafragment_step_limit / 4;
    Opt_params.intrafragment_step_limit = ((new_val < min) ? min : new_val);
    oprintf_out("\tEnergy ratio indicates iffy step: Trust radius decreased to %6.3e.\n\n",
            Opt_params.intrafragment_step_limit);
#if defined(OPTKING_PACKAGE_PSI)
    psi::Process::environment.options.set_double(module, key, Opt_params.intrafragment_step_limit);
#elif defined(OPTKING_PACKAGE_QCHEM)
    int qchem_limit_val = (int) (Opt_params.intrafragment_step_limit * 1000);
    rem_write(qchem_limit_val, REM_GEOM_OPT2_INTRAFRAG_STEP_LIMIT);
#endif
  }
  return;
}

void OPT_DATA::reset_trust_radius(void) const {
  std::string module = "OPTKING";
  std::string key = "INTRAFRAG_STEP_LIMIT";

  Opt_params.intrafragment_step_limit = Opt_params.intrafragment_step_limit_orig;

#if defined(OPTKING_PACKAGE_PSI)
    psi::Process::environment.options.set_double(module, key, Opt_params.intrafragment_step_limit);
#elif defined(OPTKING_PACKAGE_QCHEM)
    int qchem_limit_val = (int) (Opt_params.intrafragment_step_limit * 1000);
    rem_write(qchem_limit_val, REM_GEOM_OPT2_INTRAFRAG_STEP_LIMIT);
#endif
  return;
}

} // end ::opt
