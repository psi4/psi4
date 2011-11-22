/*! \file    opt_data.cc
    \ingroup optking
    \brief   OPT_DATA associated functions that do not do i/o
*/

#include "opt_data.h"

#include <cmath>

#define EXTERN
#include "globals.h"

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
  for (int i=0; i<steps.size(); ++i)
    delete steps[i];
  steps.clear();
}

// Check convergence criteria and print status to output file
// return true, if geometry is optimized
// Checks maximum force && (Delta(E) || maximum disp)
bool OPT_DATA::conv_check(opt::MOLECULE &mol) const {
  double *dq =  g_dq_pointer();
  double max_disp = array_abs_max(dq, Nintco);

  double DE;
  if (g_iteration() > 1) DE = g_energy() - g_last_energy();
  else DE = g_energy();

  double *f =  g_forces_pointer();

  // for IRC only consider forces tangent to the hypersphere search surface
  double *f_backup;
  if (Opt_params.opt_type == OPT_PARAMS::IRC) {
    // save forces and put back in below
    f_backup = init_array(Nintco);
    array_copy(f, f_backup, Nintco);

    double **G = mol.compute_G(1);
    double **Ginv = symm_matrix_inv(G, Nintco, Nintco); 
    free_matrix(G);

    // compute p_m, mass-weighted hypersphere vector
    const double *q_pivot = p_irc_data->g_q_pivot();
    double *q = mol.intco_values();
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

    if (Opt_params.print_lvl >= 1) {
      fprintf(outfile,"\tForces perpendicular to hypersphere.\n");
      print_array(outfile, f, Nintco);
    }

    fprintf(outfile,"\tFor IRC computations, the forces perpendicular to the rxnpath are tested.\n");
  }

  double max_force = array_abs_max(f, Nintco);

  fprintf(outfile, "\n\tConvergence Check Cycle %4d: (using internal coordinates in au)\n", iteration);
  fprintf(outfile, "\t                    Actual        Tolerance     Converged?\n");
  fprintf(outfile, "\t----------------------------------------------------------\n");

  if ( fabs(Opt_params.conv_max_force) < 1.0e-15 ) fprintf(outfile, "\tMAX Force        %10.1e\n", max_force);
  else fprintf(outfile, "\tMAX Force        %10.1e %14.1e %11s\n", max_force, Opt_params.conv_max_force,
       ((max_force < Opt_params.conv_max_force) ? "yes" : "no"));

  if ( fabs(Opt_params.conv_max_DE) < 1.0e-15 ) fprintf(outfile, "\tEnergy Change    %10.1e\n", fabs(DE));
  else fprintf(outfile, "\tEnergy Change    %10.1e %14.1e %11s\n", DE, Opt_params.conv_max_DE,
       ((fabs(DE) < Opt_params.conv_max_DE) ? "yes" : "no"));

  if ( fabs(Opt_params.conv_max_disp) < 1.0e-15 ) fprintf(outfile, "\tMAX Displacement %10.1e\n", max_disp);
  else fprintf(outfile, "\tMAX Displacement %10.1e %14.1e %11s\n", max_disp, Opt_params.conv_max_disp,
       ((max_disp < Opt_params.conv_max_disp) ? "yes" : "no"));

  fprintf(outfile, "\t----------------------------------------------------------\n");
  printf("\tMAX Force %10.1e : Energy Change %10.1e : MAX Displacement %10.1e\n", max_force, DE, max_disp);

  // return all forces to canonical place
  if (Opt_params.opt_type == OPT_PARAMS::IRC) {
    array_copy(f_backup, f, Nintco);
    free_array(f_backup);
  }

  // convergence test is forces and either energy change or displacement
  if ((max_force < Opt_params.conv_max_force) &&
      ((fabs(DE) < Opt_params.conv_max_DE) || (max_disp < Opt_params.conv_max_disp)))  {

#if defined(OPTKING_PACKAGE_PSI)
    // This environment variable will store the number of iterations required 
    // for convergence; it allows the db() python utilities to collect this information
    psi::Process::environment.globals["OPTIMIZATION ITERATIONS"] = g_iteration();
#endif

    return true; // structure is optimized!
  }
  else 
    return false;
}

void OPT_DATA::summary(void) const {
  double DE, *f, *dq, max_force, max_disp;

  fprintf(outfile,"\n\t            ****  Optimization Summary  ****\n");
  fprintf(outfile,"\t----------------------------------------------------------------------------\n");
  fprintf(outfile,"\t Step         Energy             Delta(E)      MAX force    MAX Delta(q)  \n");
  fprintf(outfile,"\t----------------------------------------------------------------------------\n");

  for (int i=0; i<iteration; ++i) {

    if (i == 0) DE = g_energy(i);
    else DE = g_energy(i) - g_energy(i-1); 

    f =  g_forces_pointer(i);
    max_force = array_abs_max(f, Nintco);

    dq =  g_dq_pointer(i);
    max_disp = array_abs_max(f, Nintco);

    fprintf(outfile,"\t %3d  %18.12lf  %18.12lf  %10.2e   %10.2e\n", i+1, g_energy(i),
      DE, max_force, max_disp);

  }
  fprintf(outfile,"\t----------------------------------------------------------------------\n");
  fprintf(outfile,"\n");
}

inline int sign_of_double(double d) {
  if (d>0) return 1;
  else if (d<0) return -1;
  else return 0; 
}

// do hessian update
void OPT_DATA::H_update(opt::MOLECULE & mol) {

  if (Opt_params.H_update == OPT_PARAMS::BFGS)
    fprintf(outfile,"\n\tPerforming BFGS update");
  else if (Opt_params.H_update == OPT_PARAMS::MS)
    fprintf(outfile,"\n\tPerforming Murtagh/Sargent update");
  else if (Opt_params.H_update == OPT_PARAMS::POWELL)
    fprintf(outfile,"\n\tPerforming Powell update");
  else if (Opt_params.H_update == OPT_PARAMS::BOFILL)
    fprintf(outfile,"\n\tPerforming Bofill update");
  else if (Opt_params.H_update == OPT_PARAMS::NONE) {
    fprintf(outfile,"\n\tNo Hessian update performed.\n");
    return;
  }

  int i_step, step_start, i, j;
  int step_this = steps.size()-1;

  /*** Read/compute current internals and forces ***/
  double *f, *x, *q;
  f = steps[step_this]->g_forces_pointer();
  x = steps[step_this]->g_geom_const_pointer();

  mol.set_geom_array(x);
  q = mol.intco_values();

  mol.fix_tors_near_180(); // fix configuration for torsions

  if (Opt_params.H_update_use_last == 0) { // use all available old gradients
    step_start = 0;
  }
  else {
    step_start = steps.size() - 1 - Opt_params.H_update_use_last;
    if (step_start < 0) step_start = 0;
  }
  fprintf(outfile," with previous %d gradient(s).\n", step_this-step_start);

  double *f_old, *x_old, *q_old, *dq, *dg;
  double gq, qq, qz, zz, phi;

  double **H = g_H_pointer();
  double **H_new = init_matrix(Nintco, Nintco);

  dq = init_array(Nintco);
  dg = init_array(Nintco);

  for (i_step=step_start; i_step < step_this; ++i_step) {

    // Read/compute old internals and forces
    f_old = g_forces_pointer(i_step);
    x_old = g_geom_const_pointer(i_step);

    mol.set_geom_array(x_old);
    q_old = mol.intco_values();

    for (i=0;i<Nintco;++i) {
      dq[i] = q[i] - q_old[i];
      dg[i] = (-1.0) * (f[i] - f_old[i]); // gradients -- not forces!
    }
    gq = array_dot(dq, dg, Nintco);
    qq = array_dot(dq, dq, Nintco);

    // skip Hessian updates with very small denominators (dq)(dq) and (dq)(dg_q)
    if ( (fabs(gq) < Opt_params.H_update_den_tol) ||(fabs(qq) < Opt_params.H_update_den_tol) ) {
      fprintf(outfile,"\tDenominators (dg)(dq) or (dq)(dq) are very small.\n");
      fprintf(outfile,"\t Skipping Hessian update for step %d.\n", i_step + 1);
      continue;
    }

    // skip Hessian updates if a dq element is too large ; helps to avoid use of inapplicable data
    // as well as torsional angles that stray too far from one side to the other of 180
    if ( array_abs_max(dq, Nintco) > Opt_params.H_update_dq_tol ) {
      fprintf(outfile,"\tChange in internal coordinate exceeds limit.\n");
      fprintf(outfile,"\t Skipping Hessian update for step %d.\n", i_step + 1);
      continue;
    }

    // Schlegel 1987 Ab Initio Methods in Quantum Chemistry 
    // To make formulas work for Hessian, i.e., the 2nd derivatives switch dx and dg
    if (Opt_params.H_update == OPT_PARAMS::BFGS) {
      double **temp_mat, **X;
      // Let a = dg^T.dq and X = (I - dg*dq^T) / a
      // Then H = X * H_old * X^T + dg*dg^T/a .
      X = unit_matrix(Nintco);
      for (i=0; i<Nintco; ++i)
        for (j=0; j<Nintco; ++j)
          X[i][j] -= (dg[i] * dq[j]) / gq ;

      temp_mat = init_matrix(Nintco,Nintco);
      opt_matrix_mult(X, 0, H, 0, temp_mat, 0, Nintco, Nintco, Nintco, 0);
      opt_matrix_mult(temp_mat, 0, X, 1, H_new, 0, Nintco, Nintco, Nintco, 0);

      for (i=0; i<Nintco; ++i)
        for (j=0; j<Nintco; ++j)
          H_new[i][j] += dg[i] * dg[j] / gq ;

      free_matrix(temp_mat);
      free_matrix(X);
    }
    else if (Opt_params.H_update == OPT_PARAMS::MS) {
      // Equations taken from Bofill article below
      double *Z = init_array(Nintco);
      opt_matrix_mult(H, 0, &dq, 1, &Z, 1, Nintco, Nintco, 1, 0);
      for (i=0; i<Nintco; ++i)
        Z[i] = dg[i] - Z[i];

      qz = array_dot(dq, Z, Nintco);

      for (i=0; i<Nintco; ++i)
        for (j=0; j<Nintco; ++j)
          H_new[i][j] = H[i][j] + Z[i] * Z[j] / qz ;

      free_array(Z);
    }
    else if (Opt_params.H_update == OPT_PARAMS::POWELL) {
      // Equations taken from Bofill article below
      double * Z = init_array(Nintco);
      opt_matrix_mult(H, 0, &dq, 1, &Z, 1, Nintco, Nintco, 1, 0);
      for (i=0; i<Nintco; ++i)
        Z[i] = dg[i] - Z[i];

      qz = array_dot(dq, Z, Nintco);

      for (i=0; i<Nintco; ++i)
        for (j=0; j<Nintco; ++j)
          H_new[i][j] = H[i][j] - qz/(qq*qq)*dq[i]*dq[j] + (Z[i]*dq[j] + dq[i]*Z[j])/qq;

      free_array(Z);
    }
    else if (Opt_params.H_update == OPT_PARAMS::BOFILL) {
      /* This functions performs a Bofill update on the Hessian according to
      J. M. Bofill, J. Comp. Chem., Vol. 15, pages 1-11 (1994). */
      // Bofill = (1-phi) * MS + phi * Powell
      double *Z = init_array(Nintco);
      opt_matrix_mult(H, 0, &dq, 1, &Z, 1, Nintco, Nintco, 1, 0);
      for (i=0; i<Nintco; ++i)
        Z[i] = dg[i] - Z[i];

      qz = array_dot(dq, Z, Nintco);
      zz = array_dot(Z, Z, Nintco);

      phi = 1.0 - qz*qz/(qq*zz);
      if (phi < 0.0) phi = 0.0;
      if (phi > 1.0) phi = 1.0;

      for (i=0; i<Nintco; ++i) // (1-phi)*MS
        for (j=0; j<Nintco; ++j)
          H_new[i][j] = H[i][j] + (1.0-phi) * Z[i] * Z[j] / qz ;

      for (i=0; i<Nintco; ++i)
        for (j=0; j<Nintco; ++j) // (phi * Powell)
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
      for (i=0; i<Nintco; ++i)
        for (j=0; j<Nintco; ++j)
          H_new[i][j] -= H[i][j];
  
      for (i=0; i<Nintco; ++i) {
        for (j=0; j<Nintco; ++j) {
          double val = H[i][j];
          max = ((scale_limit*val) > max_limit) ? (scale_limit*val) : max_limit;
  
        if (fabs(H_new[i][j]) < max)
          H[i][j] += H_new[i][j];
        else // limit change to max
          H[i][j] += max * sign_of_double(H_new[i][j]);
        }
      }
    }
    else { // copy H_new into H
      for (i=0; i<Nintco; ++i)
        for (j=0; j<Nintco; ++j)
          H[i][j] = H_new[i][j];
    }

    free_array(q_old);
  } // end loop over steps to use in update
  free_array(q);
  free_array(dq);
  free_array(dg);

  // put original geometry back into molecule object
  mol.set_geom_array(x);


  free_matrix(H_new);
  if (Opt_params.print_lvl >= 2) {
    fprintf(outfile, "Updated Hessian (in au)\n");
    print_matrix(outfile, H, Nintco, Nintco);
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
  rfo_eigenvector = init_array(Nintco+1);

  bool data_file_present = opt_io_is_present(); // determine if old data file is present

  if (!data_file_present) {
    fprintf(outfile, "\tPrevious optimization step data not found.  Starting new optimization.\n\n");
    iteration = 0;
    consecutive_backsteps = 0;
  }
  else {
    int Nintco_old, Ncart_old;
    opt_io_open(opt::OPT_IO_OPEN_OLD);
    opt_io_read_entry("Nintco", (char *) &Nintco_old, sizeof(int));
    opt_io_read_entry("Ncart",  (char *) &Ncart_old, sizeof(int));

    if (Nintco_old != Nintco)
      fprintf(outfile, "\tThe number of coordinates has changed.  Ignoring old data.\n");
    if (Ncart_old != Ncart)
      fprintf(outfile, "\tThe number of atoms has changed.  Ignoring old data.\n");

    if ( (Nintco_old != Nintco) || (Ncart_old != Ncart) ) {
      iteration = 0;
      consecutive_backsteps = 0;
      opt_io_close(0); // close and delete
    }
    else { // read in old optimization data
      opt_io_read_entry("H", (char *) H[0], Nintco * Nintco * sizeof(double) );
      opt_io_read_entry("iteration", (char *) &iteration, sizeof(int));
      opt_io_read_entry("consecutive_backsteps", (char *) &consecutive_backsteps, sizeof(int));
      opt_io_read_entry("rfo_eigenvector", (char *) rfo_eigenvector, (Nintco+1)*sizeof(double));
      for (int i=0; i<iteration; ++i) {
        STEP_DATA *one_step = new STEP_DATA(Nintco, Ncart);
        one_step->read(i+1, Nintco, Ncart);
        steps.push_back(one_step);
      }
      opt_io_close(1);
    }
  }

  ++iteration; // increment for current step
  // create memory for this, current step
  STEP_DATA *one_step = new STEP_DATA(Nintco, Ncart);
  steps.push_back(one_step);
}

// write data to binary file
void OPT_DATA::write(void) {
  opt_io_open(opt::OPT_IO_OPEN_OLD);

  fprintf(outfile,"\tWriting optimization data to binary file.\n");
  opt_io_write_entry("Nintco", (char *) &Nintco, sizeof(int));
  opt_io_write_entry("Ncart" , (char *) &Ncart , sizeof(int));
  opt_io_write_entry("H", (char *) H[0], Nintco * Nintco * sizeof(double) );
  opt_io_write_entry("iteration", (char *) &iteration, sizeof(int));
  opt_io_write_entry("consecutive_backsteps", (char *) &consecutive_backsteps, sizeof(int));
  opt_io_write_entry("rfo_eigenvector", (char *) rfo_eigenvector, (Nintco+1)*sizeof(double));

  for (int i=0; i<steps.size(); ++i)
    steps[i]->write(i+1, Nintco, Ncart);
  opt_io_close(1);

  fflush(outfile);
  return;
}

// Report on performance of last step
// Eventually might have this function return false to reject a step
bool OPT_DATA::previous_step_report(void) const {
  fprintf(outfile,  "\tCurrent energy   : %20.10lf\n\n", p_Opt_data->g_energy());

  if (steps.size() == 1) 
    return true;

  fprintf(outfile,"\tEnergy change for the previous step:\n");
  fprintf(outfile,"\t\tProjected    : %20.10lf\n", p_Opt_data->g_last_DE_predicted());
  fprintf(outfile,"\t\tActual       : %20.10lf\n",
      p_Opt_data->g_energy() - p_Opt_data->g_last_energy());

  double Energy_ratio = (p_Opt_data->g_energy() - p_Opt_data->g_last_energy()) / g_last_DE_predicted();

  // Minimum search
  if (Opt_params.opt_type == OPT_PARAMS::MIN) {
    if (Energy_ratio < 0.0 && consecutive_backsteps < Opt_params.consecutive_backsteps_allowed) {
      throw(BAD_STEP_EXCEPT("Energy has increased in a minimization.\n"));
    }
    else if (Energy_ratio < 0.25)
    {
      decrease_trust_radius();
    }
    else if (Energy_ratio > 0.75)
    {
      increase_trust_radius();
    }
  }

  return true;
}

// These functions are a bit out of place, given that the trust radius is not
// stored inside opt_data.
void OPT_DATA::increase_trust_radius(void) const {
  // don't let step_limit get larger than 1.0 au
  std::string module = "OPTKING";
  std::string key = "INTRAFRAGMENT_STEP_LIMIT";
  if (Opt_params.intrafragment_step_limit < 1.0) {
    Opt_params.intrafragment_step_limit *= 2;
    fprintf(outfile,"\tEnergy ratio indicates good step: Trust radius increased to %5.2e.\n\n",
        Opt_params.intrafragment_step_limit);
#if defined(OPTKING_PACKAGE_PSI)
    psi::Process::environment.options.set_double(module, key, Opt_params.intrafragment_step_limit);
#endif
  }
  return;
}

void OPT_DATA::decrease_trust_radius(void) const {
  std::string module = "OPTKING";
  std::string key = "INTRAFRAGMENT_STEP_LIMIT";
  if (Opt_params.intrafragment_step_limit > 0.0005) {
    Opt_params.intrafragment_step_limit /= 4;
    fprintf(outfile,"\tEnergy ratio indicates iffy step: Trust radius decreased to %5.2e.\n\n",
          Opt_params.intrafragment_step_limit);
#if defined(OPTKING_PACKAGE_PSI)
    psi::Process::environment.options.set_double(module, key, Opt_params.intrafragment_step_limit);
#endif
  }
  return;
}

} // end ::opt


