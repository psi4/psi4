/*! \file molecule.cc
    \ingroup optking
    \brief molecule class (really, molecular system class)
*/

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

using namespace std;

// test the analytic B matrix (and displacement code) by comparing
// analytic DqDx to finite-difference DqDx
void MOLECULE::test_B(void) {
  int Natom = g_natom();
  int Nintco = g_nintco();
  const double disp_size = 0.01;
  // 5-point formula should be good to h^4; a few will be slightly worse
  const double MAX_ERROR = 50*disp_size*disp_size*disp_size*disp_size;

  fprintf(outfile,"\n\tTesting B-matrix numerically...\n");

  double **B_analytic = compute_B();

  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile, "Analytic B matrix in au\n");
    print_matrix(outfile, B_analytic, Nintco, 3*Natom);
    fflush(outfile);
  }

  double **coord, *q_p, *q_m, **B_fd;
  double *q_p2, *q_m2;

  B_fd = init_matrix(Nintco, 3*Natom);
  coord = g_geom_2D(); // in au

  // account for changes in sign of dihedral
  fix_tors_near_180();

  try {

  for (int atom=0; atom<Natom; ++atom) {
    for (int xyz=0; xyz<3; ++xyz) {

      coord[atom][xyz] -= disp_size;
      q_m   = intco_values(coord);
      coord[atom][xyz] -= disp_size;
      q_m2  = intco_values(coord);
      coord[atom][xyz] += 3*disp_size;
      q_p  = intco_values(coord);
      coord[atom][xyz] += disp_size;
      q_p2 = intco_values(coord);
      coord[atom][xyz] -= 2*disp_size; // restore to original
      for (int i=0; i<Nintco; ++i)
        B_fd[i][3*atom+xyz] = (q_m2[i]-8*q_m[i]+8*q_p[i]-q_p2[i]) / (12.0*disp_size);
      free_array(q_p);  free_array(q_m);
      free_array(q_p2); free_array(q_m2);

    }
  }

  } catch (const char *s) {
    fprintf(outfile,"Unable to compute all internal coordinate values at displaced geometries.\n");
    fprintf(outfile,"%s\n",s);
    free_matrix(coord);
    free_matrix(B_analytic);
    free_matrix(B_fd);
    return;
  }

  free_matrix(coord);
  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile,"\nNumerical B matrix in au, disp_size = %lf\n",disp_size);
    print_matrix(outfile, B_fd, Nintco, 3*Natom);
    fflush(outfile);
  }

  double max_error = -1.0;
  int max_error_intco = -1;
  for (int i=0; i<Nintco; ++i)
    for (int j=0; j<3*Natom; ++j)
      if ( fabs(B_analytic[i][j] - B_fd[i][j]) > max_error ) {
        max_error = fabs(B_analytic[i][j] - B_fd[i][j]);
        max_error_intco = i;
      }

  fprintf(outfile,"\t\tMaximum difference is %.1e for internal coordinate %d.\n",
    max_error, max_error_intco+1);
  string coord_def = get_intco_definition_from_global_index(max_error_intco);
fprintf(outfile,"a, \n");
  fprintf(outfile,"\t\tThis coordinate is %s\n", coord_def.c_str() );
fprintf(outfile,"b, \n");
  if (max_error > MAX_ERROR) {
    fprintf(outfile, "\t\tB-matrix could be in error.  However, numerical test will fail for ");
    fprintf(outfile, "linear bond angles.  This is OK.\n");
  }
  else {
    fprintf(outfile,"\t...Passed.\n");
  }

  free_matrix(B_analytic);
  free_matrix(B_fd);
  fflush(outfile);
  return;
}

void MOLECULE::test_derivative_B(void) {
  int Natom = g_natom();
  int Nintco = g_nintco();
  const double disp_size = 0.01;
  // 5-point formula should be good to h^4; a few will be slightly worse
  const double MAX_ERROR = 10*disp_size*disp_size*disp_size*disp_size;

  bool warn;
  double **coord, *q;
  double **dq2dx2_analytic, **dq2dx2_fd;

  dq2dx2_fd = init_matrix(3*Natom, 3*Natom);
  coord = g_geom_2D();     // in au

  q = intco_values(coord); // necesessary to set torsional near-180 variables?

  fprintf(outfile,"\n\tTesting Derivative B-matrix numerically...\n");
  for (int i=0; i<Nintco; ++i) {
    warn = false;
    fprintf(outfile,"\t\tInternal coordinate %d : ", i+1); fflush(outfile);
    dq2dx2_analytic = compute_derivative_B(i);
    if (Opt_params.print_lvl >= 3) {
      fprintf(outfile, "Analytic B' (Dq2Dx2) matrix in au\n");
      print_matrix(outfile, dq2dx2_analytic, 3*Natom, 3*Natom); fflush(outfile);
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
      fprintf(outfile,"\nNumerical B' matrix by values in au, disp_size = %lf\n",disp_size);
      print_matrix(outfile, dq2dx2_fd, 3*Natom, 3*Natom);
    }

    double max_error = 0.0;
    for (int ii=0; ii<3*Natom; ++ii)
      for (int j=0; j<3*Natom; ++j)
        if ( fabs(dq2dx2_analytic[ii][j] - dq2dx2_fd[ii][j]) > max_error )
          max_error = fabs(dq2dx2_analytic[ii][j] - dq2dx2_fd[ii][j]);

    fprintf(outfile,"Maximum difference is %.1e. ", max_error);
    if (max_error > MAX_ERROR) {
      fprintf(outfile, "Uh-Oh.  See below\n");
      warn = true;
    }
    else { fprintf(outfile," Passed.\n"); }

    if (warn) {
      fprintf(outfile, "\nWarning: Perhaps a bug or your angular coordinates are at a discontinuity.\n");
      fprintf(outfile, "Try restarting your optimization at a new or updated geometry.\n");
      fprintf(outfile, "Also, remove angular coordinates that are fixed by symmetry.\n");
    }

    free_matrix(dq2dx2_analytic);
    fflush(outfile);
  }
  fprintf(outfile,"\n");
  free_matrix(dq2dx2_fd);
  fflush(outfile);
  return;
}

}

