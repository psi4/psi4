/*! \file molecule.cc
    \ingroup OPT10
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
  const double disp_size = 0.001;

  fprintf(outfile,"\n\tTesting B-matrix numerically...\n");

  double **B_analytic = compute_B();
  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile, "Analytic B matrix in au\n");
    print_matrix(outfile, B_analytic, Nintco, 3*Natom);
  }

  double **coord, *q_plus, *q_minus, **B_fd;

  coord = g_geom_2D(); // in au
  B_fd = init_matrix(Nintco, 3*Natom);

  for (int atom=0; atom<Natom; ++atom) {
    for (int xyz=0; xyz<3; ++xyz) {
      coord[atom][xyz] += disp_size;
      q_plus  = intco_values(coord);
      coord[atom][xyz] -= 2.0*disp_size;
      q_minus = intco_values(coord);
      coord[atom][xyz] += disp_size; // restore

      for (int i=0; i<Nintco; ++i)
        B_fd[i][3*atom+xyz] = (q_plus[i]-q_minus[i]) / (2.0*disp_size);

      free_array(q_plus);
      free_array(q_minus);
    }
  }
  free_matrix(coord);
  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile,"\nNumerical B matrix in au, disp_size = %lf\n",disp_size);
    print_matrix(outfile, B_fd, Nintco, 3*Natom);
  }

  double max_error = 0.0;
  for (int i=0; i<Nintco; ++i)
    for (int j=0; j<3*Natom; ++j)
      if ( fabs(B_analytic[i][j] - B_fd[i][j]) > max_error )
        max_error = fabs(B_analytic[i][j] - B_fd[i][j]);

  fprintf(outfile,"\t\tMaximum difference is %.1e.", max_error);
  if (max_error > 1.0e-3) {
    fprintf(outfile, "\nUh-Oh.  Perhaps a bug or your angular coordinates are at a discontinuity.\n");
    fprintf(outfile, "If the latter, restart your optimization at a new or updated geometry.\n");
    fprintf(outfile, "Remove angular coordinates that are fixed by symmetry\n");
  }
  else {
    fprintf(outfile,"  Passed.\n");
  }

  free_matrix(B_analytic);
  free_matrix(B_fd);
  return;
}

void MOLECULE::test_derivative_B(void) {
  int Natom = g_natom();
  int Nintco = g_nintco();
  const double disp_size = 0.001;
  bool warn = false;

  double **coord, *q;
  double **dq2dx2_analytic, **dq2dx2_fd;

  coord = g_geom_2D(); // in au
  dq2dx2_fd = init_matrix(3*Natom, 3*Natom);

  q = intco_values(coord);

  fprintf(outfile,"\n\tTesting Derivative B-matrix by finite differences...\n");
  for (int i=0; i<Nintco; ++i) {
    fprintf(outfile,"\t\tInternal coordinate %d : ", i+1);
    dq2dx2_analytic = compute_derivative_B(i);
    if (Opt_params.print_lvl >= 3) {
      fprintf(outfile, "Analytic B' (Dq2Dx2) matrix in au\n");
      print_matrix(outfile, dq2dx2_analytic, 3*Natom, 3*Natom);
    }

    for (int atom_a=0; atom_a<Natom; ++atom_a) {
      for (int xyz_a=0; xyz_a<3; ++xyz_a) {

        for (int atom_b=0; atom_b<Natom; ++atom_b) {
          for (int xyz_b=0; xyz_b<3; ++xyz_b) {

            if (atom_a == atom_b && xyz_a == xyz_b) { // diagonal elements
              double *q_plus, *q_minus;
              coord[atom_a][xyz_a] += disp_size;
              q_plus  = intco_values(coord);
              coord[atom_a][xyz_a] -= 2.0*disp_size;
              q_minus = intco_values(coord);
              coord[atom_a][xyz_a] += disp_size; // restore coord

              dq2dx2_fd[3*atom_a+xyz_a][3*atom_a+xyz_a] = (q_plus[i]+q_minus[i]-2.0*q[i])
                  / (disp_size*disp_size);

              free_array(q_plus); free_array(q_minus);
            }
            else { // off-diagonal
              double *q_pp, *q_pm, *q_mp, *q_mm;
              coord[atom_a][xyz_a] += disp_size;
              coord[atom_b][xyz_b] += disp_size;
              q_pp  = intco_values(coord);
              coord[atom_a][xyz_a] -= 2.0*disp_size;
              coord[atom_b][xyz_b] -= 2.0*disp_size;
              q_mm  = intco_values(coord);
              coord[atom_a][xyz_a] += 2.0*disp_size;
              q_pm  = intco_values(coord);
              coord[atom_a][xyz_a] -= 2.0*disp_size;
              coord[atom_b][xyz_b] += 2.0*disp_size;
              q_mp  = intco_values(coord);
              coord[atom_a][xyz_a] += disp_size; //restore coord
              coord[atom_b][xyz_b] -= disp_size;

              dq2dx2_fd[3*atom_a+xyz_a][3*atom_b+xyz_b] = (q_pp[i]+q_mm[i]-q_pm[i]-q_mp[i])
                  / (4.0*disp_size*disp_size);

              free_array(q_pp); free_array(q_mm); free_array(q_pm); free_array(q_mp);
            }
          }
        } // atom_b
      }
    } // atom_a

    if (Opt_params.print_lvl >= 3) {
      fprintf(outfile,"\nNumerical B' matrix in au, disp_size = %lf\n",disp_size);
      print_matrix(outfile, dq2dx2_fd, 3*Natom, 3*Natom);
    }

    double max_error = 0.0;
    for (int i=0; i<3*Natom; ++i)
      for (int j=0; j<3*Natom; ++j)
        if ( fabs(dq2dx2_analytic[i][j] - dq2dx2_fd[i][j]) > max_error )
          max_error = fabs(dq2dx2_analytic[i][j] - dq2dx2_fd[i][j]);

    fprintf(outfile,"Maximum difference is %.1e. ", max_error);
    if (max_error > 5.0e-3) {
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
  return;
}

}

