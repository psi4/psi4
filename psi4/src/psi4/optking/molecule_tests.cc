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

using namespace std;

// test the analytic B matrix (and displacement code) by comparing
// analytic DqDx to finite-difference DqDx
void MOLECULE::test_B(void) {
  int Natom = g_natom();
  int Nintco = Ncoord();
  const double disp_size = 0.01;
  // 5-point formula should be good to h^4; a few will be slightly worse
  const double MAX_ERROR = 50*disp_size*disp_size*disp_size*disp_size;

  oprintf_out("\n\tTesting B-matrix numerically...\n");

  double **B_analytic = compute_B();

  if (Opt_params.print_lvl >= 3) {
    oprintf_out( "Analytic B matrix in au\n");
    oprint_matrix_out(B_analytic, Nintco, 3*Natom);
  }

  double **coord, *q_p, *q_m, **B_fd;
  double *q_p2, *q_m2;

  B_fd = init_matrix(Nintco, 3*Natom);
  coord = g_geom_2D(); // in au

  // account for changes in sign of dihedral
  fix_tors_near_180();
  fix_oofp_near_180();
  fix_bend_axes();

  try {

  for (int atom=0; atom<Natom; ++atom) {
    for (int xyz=0; xyz<3; ++xyz) {

      coord[atom][xyz] -= disp_size;
      q_m   = coord_values(coord);
      coord[atom][xyz] -= disp_size;
      q_m2  = coord_values(coord);
      coord[atom][xyz] += 3*disp_size;
      q_p  = coord_values(coord);
      coord[atom][xyz] += disp_size;
      q_p2 = coord_values(coord);
      coord[atom][xyz] -= 2*disp_size; // restore to original
      for (int i=0; i<Nintco; ++i)
        B_fd[i][3*atom+xyz] = (q_m2[i]-8*q_m[i]+8*q_p[i]-q_p2[i]) / (12.0*disp_size);
      free_array(q_p);  free_array(q_m);
      free_array(q_p2); free_array(q_m2);

    }
  }

  } catch (const char *s) {
    oprintf_out("Unable to compute all internal coordinate values at displaced geometries.\n");
    oprintf_out("%s\n",s);
    free_matrix(coord);
    free_matrix(B_analytic);
    free_matrix(B_fd);
    return;
  }

  free_matrix(coord);
  if (Opt_params.print_lvl >= 3) {
    oprintf_out("\nNumerical B matrix in au, disp_size = %lf\n",disp_size);
    oprint_matrix_out(B_fd, Nintco, 3*Natom);

  }

  double max_error = -1.0;
  int max_error_intco = -1;
  for (int i=0; i<Nintco; ++i)
    for (int j=0; j<3*Natom; ++j)
      if ( fabs(B_analytic[i][j] - B_fd[i][j]) > max_error ) {
        max_error = fabs(B_analytic[i][j] - B_fd[i][j]);
        max_error_intco = i;
      }

  oprintf_out("\t\tMaximum difference is %.1e for internal coordinate %d.\n",
    max_error, max_error_intco+1);
  string coord_def = get_coord_definition_from_global_index(max_error_intco);
  oprintf_out("\t\tThis coordinate is %s\n", coord_def.c_str() );
  if (max_error > MAX_ERROR) {
    oprintf_out("\tB-matrix could be in error. However, numerical test will slightly\n");
    oprintf_out("\tfail for linear bond angles.  This is OK.\n");
  }
  else {
    oprintf_out("\t...Passed.\n");
  }

  unfix_bend_axes();
  free_matrix(B_analytic);
  free_matrix(B_fd);

  return;
}

void MOLECULE::test_derivative_B(void) {
  int Natom = g_natom();
  int Nintco = Ncoord();
  const double disp_size = 0.01;
  // 5-point formula should be good to h^4; a few will be slightly worse
  const double MAX_ERROR = 10*disp_size*disp_size*disp_size*disp_size;

  bool warn;
  double **coord, *q;
  double **dq2dx2_analytic, **dq2dx2_fd;

  dq2dx2_fd = init_matrix(3*Natom, 3*Natom);
  coord = g_geom_2D();     // in au

  q = coord_values(coord); // necesessary to set torsional near-180 variables?

  oprintf_out("\n\tTesting Derivative B-matrix numerically...\n");
  for (int i=0; i<Nintco; ++i) {
    warn = false;
    oprintf_out("\t\tInternal coordinate %d : ", i+1);
    dq2dx2_analytic = compute_derivative_B(i);
    if (Opt_params.print_lvl >= 3) {
      oprintf_out( "Analytic B' (Dq2Dx2) matrix in au\n");
      oprint_matrix_out(dq2dx2_analytic, 3*Natom, 3*Natom);
    }

    // compute B' matrix from B matrices
    for (int atom_a=0; atom_a<Natom; ++atom_a) {
      for (int xyz_a=0; xyz_a<3; ++xyz_a) {
        double **B_p, **B_p2, **B_m, **B_m2;

        coord[atom_a][xyz_a] += disp_size;
        set_geom_array(coord[0]);
        B_p  = compute_B();

        coord[atom_a][xyz_a] += disp_size;
        set_geom_array(coord[0]);
        B_p2 = compute_B();

        coord[atom_a][xyz_a] -= 3.0*disp_size;
        set_geom_array(coord[0]);
        B_m  = compute_B();

        coord[atom_a][xyz_a] -= disp_size;
        set_geom_array(coord[0]);
        B_m2 = compute_B();

        coord[atom_a][xyz_a] += 2*disp_size; // restore coord
        set_geom_array(coord[0]);

        for (int atom_b=0; atom_b<Natom; ++atom_b)
          for (int xyz_b=0; xyz_b<3; ++xyz_b)
            dq2dx2_fd[3*atom_a+xyz_a][3*atom_b+xyz_b] =
             (B_m2[i][3*atom_b+xyz_b]-8*B_m[i][3*atom_b+xyz_b]+8*B_p[i][3*atom_b+xyz_b]-B_p2[i][3*atom_b+xyz_b])
               / (12.0*disp_size);

        free_matrix(B_p); free_matrix(B_m);
        free_matrix(B_p2); free_matrix(B_m2);
      }
    } // atom_a

    if (Opt_params.print_lvl >= 3) {
      oprintf_out("\nNumerical B' matrix by values in au, disp_size = %lf\n",disp_size);
      oprint_matrix_out(dq2dx2_fd, 3*Natom, 3*Natom);
    }

    double max_error = 0.0;
    for (int ii=0; ii<3*Natom; ++ii)
      for (int j=0; j<3*Natom; ++j)
        if ( fabs(dq2dx2_analytic[ii][j] - dq2dx2_fd[ii][j]) > max_error )
          max_error = fabs(dq2dx2_analytic[ii][j] - dq2dx2_fd[ii][j]);

    oprintf_out("Maximum difference is %.1e. ", max_error);
    if (max_error > MAX_ERROR) {
      oprintf_out( "Uh-Oh.  See below\n");
      warn = true;
    }
    else { oprintf_out(" Passed.\n"); }

    if (warn) {
      oprintf_out( "\nWarning: Perhaps a bug or your angular coordinates are at a discontinuity.\n");
      oprintf_out( "Try restarting your optimization at a new or updated geometry.\n");
      oprintf_out( "Also, remove angular coordinates that are fixed by symmetry.\n");
    }

    free_matrix(dq2dx2_analytic);

  }
  oprintf_out("\n");
  free_matrix(dq2dx2_fd);

  return;
}

}
