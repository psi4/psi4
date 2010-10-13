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

bool is_integer(const char *check);

using namespace std;

// constructor reads geometries of fragments and creates fragment objects
// Currently, the assumed format is :
// header row
// natom_i natom_j ... [# in each fragment]  total energy
// atomic_symbol/atomic_number  x y z
// ...
// dE/dx  dE/dy   dE/dz
//
MOLECULE::MOLECULE(ifstream & fin) {
  string lbl;
  stringstream sline (stringstream::in | stringstream::out);
  stringstream sline2 (stringstream::in | stringstream::out);
  char cline[256];

  try {
    fin.clear();
    fin.exceptions(ios_base::failbit);

    fin.getline(cline, 256);   // header row
    fin.getline(cline, 256);   // natom_i and energy row
    sline << cline;

    // cnt number of fragments (integer entries)
    sline >> lbl;
    int nfrag = 0;
    while (is_integer(lbl.c_str())) {
      nfrag++;
      sline >> lbl; 
    }

    // record natoms for each fragment
    sline2 << cline;
    int * natom_i = init_int_array(nfrag);
    int i, frag, nallatom=0;
    for (frag=0; frag<nfrag; ++frag) {
      sline2 >> natom_i[frag];
      nallatom += natom_i[frag];
    }

    // rewind, determine number of geometries
    fin.seekg(0);
    int nlines = 0;
    try { // ignore exceptions only here
      while (fin.peek() != EOF) {
        fin.getline(cline, 256);
        ++nlines;
      }
    }
    catch (std::ios_base::failure & bf) { }
    int ngeom = nlines / (2*nallatom+2);

    fin.clear();
    fin.seekg(0);
    for (i=0; i< (ngeom-1)*(2*nallatom+2); ++i)
      fin.getline(cline, 256);

    fin.getline(cline, 256); // header line

    for (frag=0; frag<nfrag; ++frag)
      fin >> i;
    fin >> energy;

    FRAG * one_frag;

    for (frag=0; frag<nfrag; ++frag) {
      double *Z = init_array(natom_i[frag]);
      double **geom = init_matrix(natom_i[frag], 3);

      for (i=0; i<natom_i[frag]; ++i) {
        // check to see if atomic number or symbol
        if ( isalpha((lbl.c_str())[0]) )
          Z[i] = Element_to_Z(lbl);
        else
          fin >> Z[i];

        fin >> geom[i][0];
        fin >> geom[i][1];
        fin >> geom[i][2];
      }

      one_frag = new FRAG(natom_i[frag], Z, geom);
      fragments.push_back(one_frag);
    } // end loop nfrag

    for (frag=0; frag<nfrag; ++frag) {
      double **grad = init_matrix(natom_i[frag], 3);
      for (i=0; i<natom_i[frag]; ++i) {
        fin >> grad[i][0];
        fin >> grad[i][1];
        fin >> grad[i][2];
      }
      fragments[frag]->set_grad(grad);
      free_matrix(grad);
    }
    free_int_array(natom_i);

  } // end try reading geometries
  catch (std::ios_base::failure & bf) {
    printf("Error reading molecular geometry and gradient\n");
    fprintf(outfile,"Error reading molecular geometry and gradient\n");
    throw "Error reading molecular geometry and gradient\n";
  }

  return;
}


bool is_integer(const char *check) {
  for (int i=0; i<strlen(check); ++i) {
    if (!isdigit(check[i]))
      return false;
  }
  return true;
}


// compute forces in internal coordinates in au
// forces in internal coordinates, f_q = G_inv B u f_x
// if u is unit matrix, f_q = (BB^T)^(-1) * B f_x
void MOLECULE::forces (void) {
  double *f_x, *temp_arr, **B, **G, **G_inv;
  int Ncart = 3*g_natom();
  int Nintco = g_nintco();

  B = compute_B();
  if (Opt_params.print_lvl >= 2) {
    fprintf(outfile, "B matrix\n");
    print_matrix(outfile, B, Nintco, Ncart);
  }

  f_x = g_grad_array(); // Hartree / bohr
  array_scm(f_x, -1, Ncart); // switch gradient -> forces

  temp_arr = init_array(Nintco);
  opt_matrix_mult(B, 0, &f_x, 1, &temp_arr, 1, Nintco, Ncart, 1, 0);
  free_array(f_x);

  G = init_matrix(Nintco, Nintco);
  opt_matrix_mult(B, 0, B, 1, G, 0, Nintco, Ncart, Nintco, 0);
  free_matrix(B);

  G_inv = symm_matrix_inv(G, Nintco, 1);
  free_matrix(G);

  double * f_q = p_Opt_data->g_forces_pointer();
  opt_matrix_mult(G_inv, 0, &temp_arr, 1, &f_q, 1, Nintco, Nintco, 1, 0);
  free_matrix(G_inv);
  free_array(temp_arr);

  if (Opt_params.print_lvl > 3) {
    fprintf(outfile,"Internal forces in au\n");
    print_matrix(outfile, &f_q, 1, Nintco);

    // test by transforming f_q back to cartesian forces and compare
    B = compute_B();
    temp_arr = init_array(Ncart);
    opt_matrix_mult(B, 1, &f_q, 1, &temp_arr, 1, Ncart, Nintco, 1, 0);
    fprintf(outfile,"Recomputed forces in cartesian coordinates\n");
    print_matrix(outfile, &temp_arr, 1, Ncart);
    free_array(temp_arr);
    free_matrix(B);
  }

  return ;
}

// project redundancies (and constraints) out of forces and Hessian matrix
// add constraints here later
void MOLECULE::project_f_and_H(void) {
  int Nintco = g_nintco();
  int Ncart = 3*g_natom();

  double **G = init_matrix(Nintco, Nintco);
  double **B = compute_B();
  opt_matrix_mult(B, 0, B, 1, G, 0, Nintco, Ncart, Nintco, 0);
  free_matrix(B);

  double **G_inv = symm_matrix_inv(G, Nintco, 1);
  double **P = init_matrix(Nintco, Nintco);
  opt_matrix_mult(G, 0, G_inv, 0, P, 0, Nintco, Nintco, Nintco, 0);
  free_matrix(G);
  free_matrix(G_inv);

  double *f_q = p_Opt_data->g_forces_pointer();
  double **H = p_Opt_data->g_H_pointer();

  // f_q~ = P f_q
  double * temp_arr = init_array(Nintco);
  opt_matrix_mult(P, 0, &f_q, 1, &temp_arr, 1, Nintco, Nintco, 1, 0);
  array_copy(temp_arr, f_q, Nintco);
  free_array(temp_arr);

  // H~ = PHP + 1000(1-P)
  double ** temp_mat = init_matrix(Nintco, Nintco);
  opt_matrix_mult(H, 0, P, 0, temp_mat, 0, Nintco, Nintco, Nintco, 0);
  opt_matrix_mult(P, 0, temp_mat, 0, H, 0, Nintco, Nintco, Nintco, 0);
  free_matrix(temp_mat);

  for (int i=0; i<Nintco; ++i)
    P[i][i] = (P[i][i] - 1) * 2000;

  for (int i=0; i<Nintco; ++i)
    for (int j=0; j<Nintco; ++j)
      H[i][j] -= P[i][j];

  free_matrix(P);
}

void MOLECULE::write_geom(void) {
  for (int i=0; i<fragments.size(); ++i)
    fragments[i]->write_geom(outfile);
}

}

#ifdef PSI4
#include <libmints/molecule.h>
#include <libchkpt/chkpt.h>
namespace opt {

void MOLECULE::write_geom_chkpt(void) {
  using namespace psi;
  double *x = g_geom_array();

  double **geom_2D = init_matrix(g_natom(), 3);
  int cnt = 0;
  for (int i=0; i<g_natom(); ++i)
    for (int xyz=0; xyz < 3; ++xyz)
      geom_2D[i][xyz] = x[cnt++];

fprintf(outfile,"Print new geom\n");
print_matrix(outfile, geom_2D, g_natom(), 3);

  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_geom(geom_2D);
  chkpt_close();
}

void MOLECULE::write_geom_to_active_molecule() {
    using namespace psi;
    double *x = g_geom_array();

    double **geom_2D = init_matrix(g_natom(), 3);
    int cnt = 0;
    for (int i=0; i<g_natom(); ++i)
      for (int xyz=0; xyz < 3; ++xyz)
        geom_2D[i][xyz] = x[cnt++];

    // Modify "active molecule" with new geometry
    Process::environment.molecule()->set_geometry(geom_2D);
}

}
#endif

namespace opt {

void MOLECULE::apply_intrafragment_step_limit(double * & dq) {
  int i, f;
  double scale = 1.0;
  double limit = Opt_params.intrafragment_step_limit;

  for (f=0; f<fragments.size(); ++f)
    for (i=0; i<fragments[f]->g_nintco(); ++i)
      if (scale * fabs(dq[g_intco_offset(f)+i]) > limit)
        scale = limit / fabs(dq[g_intco_offset(f)+i]);

  if (scale != 1.0) {
    fprintf(outfile,"\tChange in coordinate exceeds step limit of %10.5lf.\n", limit);
    fprintf(outfile,"\tScaling displacements by %10.5lf\n", scale);

    for (f=0; f<fragments.size(); ++f)
      for (i=0; i<fragments[f]->g_nintco(); ++i)
        dq[g_intco_offset(f)+i] *= scale;
  }
}

// don't let any angles get smaller than 0.0
void MOLECULE::check_intrafragment_zero_angles(double const * const dq) {
  for (int f=0; f<fragments.size(); ++f)
    fragments[f]->check_zero_angles(&(dq[g_intco_offset(f)]));
}

void MOLECULE::H_guess(void) const {
  double **H, **H_frag;
  int f, i, j;

  H = p_Opt_data->g_H_pointer();
  for (f=0; f<fragments.size(); ++f) {
    H_frag = fragments[f]->H_guess();
    for (i=0; i<fragments[f]->g_nintco(); ++i)
      for (j=0; j<fragments[f]->g_nintco(); ++j)
        H[g_intco_offset(f) + i][g_intco_offset(f) + j] = H_frag[i][j];
    free_matrix(H_frag);
  }

  if (Opt_params.print_lvl >= 2) {
    fprintf(outfile,"Initial Hessian guess\n");
    print_matrix(outfile, H, g_nintco(), g_nintco());
  }
  return;
}

// test the analytic B matrix (and displacement code) by comparing
// analytic DqDx to finite-difference DqDx
/*
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
*/

/*
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
*/

double ** MOLECULE::cartesian_H_to_internals(void) const {
  int Nintco = g_nintco();
  int Ncart = 3*g_natom();

  // compute A = u B^t (B u B^t)^-1 where u=unit matrix and -1 is generalized inverse
  double **B = compute_B();
  double **G = init_matrix(Nintco, Nintco);
  opt_matrix_mult(B, 0, B, 1, G, 0, Nintco, Ncart, Nintco, 0);

  double **G_inv = symm_matrix_inv(G, Nintco, true);
  free_matrix(G);

  double **A = init_matrix(Ncart, Nintco);
  opt_matrix_mult(B, 1, G_inv, 0, A, 0, Ncart, Nintco, Nintco, 0);
  free_matrix(G_inv);
  free_matrix(B);

  // compute gradient in internal coordinates, A^t g_x = g_q
  double *grad_x = g_grad_array();
  double *grad_q = init_array(Nintco);
  opt_matrix_mult(A, 1, &grad_x, 1, &grad_q, 1, Nintco, Ncart, 1, 0);
  free_array(grad_x);

  // read in cartesian H
  double **H_cart = p_Opt_data->read_cartesian_H();

  // transform cartesian H to internals; A^t (H_x - K) A
  // K_ij = sum_q ( grad_q[q] d^2(q)/(dxi dxj) )
  double **dq2dx2;

  for (int q=0; q<Nintco; ++q) {
    dq2dx2 = compute_derivative_B(q); // d^2(q)/ dx_i dx_j

    for (int i=0; i<Ncart; ++i)
      for (int j=0; j<Ncart; ++j)
        H_cart[i][j] -= grad_q[q] * dq2dx2[i][j];

    free_matrix(dq2dx2);
  }
  free_array(grad_q);

  double **temp_mat = init_matrix(Ncart, Nintco);
  opt_matrix_mult(H_cart, 0, A, 0, temp_mat, 0, Ncart, Ncart, Nintco, 0);
  free_matrix(H_cart);

  double **H_int = init_matrix(Nintco, Nintco);
  opt_matrix_mult(A, 1, temp_mat, 0, H_int, 0, Nintco, Ncart, Nintco, 0);
  free_matrix(temp_mat);

  free_matrix(A);

  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile, "Hessian transformed to internal coordinates:\n");
    print_matrix(outfile, H_int, Nintco, Nintco); fflush(outfile);
  }

  // Check by transforming internal coordinate Hessian back into cartesian coordinates:
/*
  B = compute_B();
  temp_mat = init_matrix(Ncart, Nintco);
  opt_matrix_mult(B, 1, H_int, 0, temp_mat, 0, Ncart, Nintco, Nintco, 0);
  H_cart = init_matrix(Ncart, Ncart);
  opt_matrix_mult(temp_mat, 0, B, 0, H_cart, 0, Ncart, Nintco, Ncart, 0);
  free_matrix(temp_mat);

  for (int q=0; q<Nintco; ++q) {
    dq2dx2 = compute_derivative_B(q); // d^2(q)/ dx_i dx_j
    for (int i=0; i<Ncart; ++i)
      for (int j=0; j<Ncart; ++j)
        H_cart[i][j] += grad_q[q] * dq2dx2[i][j];
    free_matrix(dq2dx2);
  }
  free_array(grad_q);
  fprintf(outfile, "Hessian transformed back into Cartesian coordinates\n");
  print_matrix(outfile, H_cart, Ncart, Ncart); fflush(outfile);
  free_matrix(B);
*/
  return H_int;
}

}


