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
#include "psi4/libparallel/ParallelPrinter.h"
#include "print.h"
#define EXTERN
#include "globals.h"
#include "linear_algebra.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
 #include "psi4/libmints/molecule.h"
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
 #include "EFP.h"
#endif

namespace opt {

using namespace std;

// if allocate_fragment, then read the number of atoms and allocate memory for
//   a fragment of that size.  Otherwise, this is an empty constructor.

MOLECULE::MOLECULE(int num_atoms) {

  if (num_atoms > 0) {
    FRAG *one_frag = new FRAG(num_atoms);
    fragments.push_back(one_frag);
  }

  return;
}

// return vector of reciprocal masses of dimension = num of cartesians
double * MOLECULE::g_u_vector(void) const {
  double *m = g_masses();
  int Natom = g_natom();

  double *u = init_array(3*Natom);
  for (int a=0; a<Natom; ++a) {
    u[3*a+0] = 1.0 / m[a];
    u[3*a+1] = 1.0 / m[a];
    u[3*a+2] = 1.0 / m[a];
  }
  return u;
}

// compute forces in internal coordinates in au
// forces in internal coordinates, f_q = G_inv B u f_x
// if u is unit matrix, f_q = (BB^T)^(-1) * B f_x
void MOLECULE::forces(void) {
  double *f_x, *temp_arr, **B, **G, **G_inv;
  int Ncart = 3*g_natom();
  int Nintco = Ncoord();

  // f_x
  f_x = g_grad_array(); // Hartree / bohr
  array_scm(f_x, -1, Ncart); // switch gradient -> forces

  if (Opt_params.print_lvl > 3)
    oprint_array_out_precise(f_x, Ncart);

  // B (u f_x)
  B = compute_B();
  if (Opt_params.print_lvl >= 3) {
    oprintf_out( "B matrix\n");
    oprint_matrix_out(B, Nintco, Ncart);
  }
  temp_arr = init_array(Nintco);
  opt_matrix_mult(B, 0, &f_x, 1, &temp_arr, 1, Nintco, Ncart, 1, 0);
  free_array(f_x);

  // G^-1 = (BuBt)^-1
  G = init_matrix(Nintco, Nintco);
  for (int i=0; i<Nintco; ++i)
    for (int k=0; k<Ncart; ++k)
      for (int j=0; j<Nintco; ++j)
        G[i][j] += B[i][k] * /* u[k] * */ B[j][k];
  free_matrix(B);

  G_inv = symm_matrix_inv(G, Nintco, 1);
  free_matrix(G);

  double * f_q = p_Opt_data->g_forces_pointer();
  opt_matrix_mult(G_inv, 0, &temp_arr, 1, &f_q, 1, Nintco, Nintco, 1, 0);
  free_matrix(G_inv);
  free_array(temp_arr);

  // append exernally determined fb forces
  double * fb_force;
  for (std::size_t f=0; f<fb_fragments.size(); ++f) {
    fb_force = fb_fragments[f]->get_forces_pointer();
    for (int i=0; i<fb_fragments[f]->Ncoord(); ++i)
      f_q[ g_fb_fragment_coord_offset(f) + i ] =  fb_force[i] ;
  }

  if (Opt_params.print_lvl >= 3) {
    oprintf_out("Internal forces in au\n");
    oprint_array_out_precise(f_q, Ncoord());
  }

/*
    // test by transforming f_q back to cartesian forces and compare
  if (Opt_params.print_lvl > 3) {
    B = compute_B();
    temp_arr = init_array(Ncart);
    opt_matrix_mult(B, 1, &f_q, 1, &temp_arr, 1, Ncart, Nintco, 1, 0);
    oprintf_out("Recomputed forces in cartesian coordinates\n");
    oprint_matrix_out(&temp_arr, 1, Ncart);
    free_array(temp_arr);
    free_matrix(B);
  }
*/

  return ;
}

// Tell whether there are any fixed equilibrium values
bool MOLECULE::has_fixed_eq_vals(void) {
  for (std::size_t f=0; f<fragments.size(); ++f)
    for (int i=0; i<fragments[f]->Ncoord(); ++i)
      if (fragments[f]->coord_has_fixed_eq_val(i))
        return true;

  return false;
}

// Is any coordinate present that is not a cartesian.
bool MOLECULE::is_noncart_present(void) const {

  if (interfragments.size()) return true;

  for (std::size_t f=0; f<fragments.size(); ++f)
    if (fragments[f]->is_noncart_present())
      return true;

  return false;
}

// Apply extra forces for internal coordinates with user-defined
// equilibrium values.
void MOLECULE::apply_constraint_forces(void) {
  double * f_q = p_Opt_data->g_forces_pointer();
  double **H = p_Opt_data->g_H_pointer();
  int N = Ncoord();
  int iter = p_Opt_data->g_iteration();
  double k;

  int cnt = -1;
  for (std::size_t f=0; f<fragments.size(); ++f) {
    for (int i=0; i<fragments[f]->Ncoord(); ++i) {
      ++cnt;
      if (fragments[f]->coord_has_fixed_eq_val(i)) {
        double eq_val = fragments[f]->coord_fixed_eq_val(i);
        double val = fragments[f]->coord_value(i);

        // Increase force constant by 5% of initial value per iteration
        k = (1 + 0.05 * (iter-1)) * Opt_params.fixed_coord_force_constant;
        H[cnt][cnt] = k;

        double force = (eq_val - val) * k;
        oprintf_out("\tAdding user-defined constraint: Fragment %d; Coordinate %d:\n", f+1, i+1);
        oprintf_out("\t\tValue = %12.6f; Fixed value    = %12.6f\n", val, eq_val);
        oprintf_out("\t\tForce = %12.6f; Force constant = %12.6f\n", force, k);
        f_q[cnt] = force;

        // If user eq. value is specified delete coupling between this coordinate and others.
        oprintf_out("\tRemoving off-diagonal coupling between coordinate %d and others.\n", cnt+1);
        for (int j=0; j<N; ++j)
          if (j != cnt)
            H[j][cnt] = H[cnt][j] = 0;
      }
    }
  }
  return;
}

// project redundancies (and constraints) out of forces and Hessian matrix
// add constraints here later
void MOLECULE::project_f_and_H(void) {
  int Nintco = Ncoord();

  // compute G = B B^t
  double **G = compute_G(false);

  // Put 1's on diagonal for FB coordinates
  for (std::size_t I=0; I<fb_fragments.size(); ++I)
    for (int i=0; i<fb_fragments[i]->Ncoord(); ++i)
      G[g_fb_fragment_coord_offset(I) + i][g_fb_fragment_coord_offset(I) + i] = 1.0;

  // compute P = G G^-1
  double **G_inv = symm_matrix_inv(G, Nintco, 1);
  double **P = init_matrix(Nintco, Nintco);
  opt_matrix_mult(G, 0, G_inv, 0, P, 0, Nintco, Nintco, Nintco, 0);
  free_matrix(G);
  free_matrix(G_inv);

  if (Opt_params.print_lvl >= 3) {
    oprintf_out("\tProjection matrix for redundancies.\n");
    oprint_matrix_out(P, Nintco, Nintco);
  }

  // add constraints to projection matrix
  double **C = compute_constraints();

  bool constraints_present = false;
  for (int i=0; i<Nintco; ++i)
    if (C[i][i] != 0)
      constraints_present = true;

  // P = P' - P' C (CPC)^-1 C P'
  if (constraints_present) {
    double **T = init_matrix(Nintco,Nintco);
    opt_matrix_mult(P, 0, C, 0,  T, 0, Nintco, Nintco, Nintco, 0);
    double **T2 = init_matrix(Nintco,Nintco);
    opt_matrix_mult(C, 0, T, 0, T2, 0, Nintco, Nintco, Nintco, 0);
    double **T3 = symm_matrix_inv(T2, Nintco, 1);

    opt_matrix_mult( C, 0,  P, 0,  T, 0, Nintco, Nintco, Nintco, 0);
    opt_matrix_mult(T3, 0,  T, 0, T2, 0, Nintco, Nintco, Nintco, 0);
    free_matrix(T);
    opt_matrix_mult( C, 0, T2, 0, T3, 0, Nintco, Nintco, Nintco, 0);
    opt_matrix_mult( P, 0, T3, 0, T2, 0, Nintco, Nintco, Nintco, 0);
    free_matrix(T3);
    for (int i=0; i<Nintco; ++i)
      for (int j=0; j<Nintco; ++j)
        P[i][j] -= T2[i][j];
    free_matrix(T2);
  }
  free_matrix(C);

  // Project redundancies and contraints out of forces
  // Unconstrained forces will not need this unless we have interpolated forces,
  // or some other as yet unimplemented method.
  double *f_q = p_Opt_data->g_forces_pointer();
  // f_q~ = P f_q
  double * temp_arr = init_array(Nintco);
  opt_matrix_mult(P, 0, &f_q, 1, &temp_arr, 1, Nintco, Nintco, 1, 0);
  array_copy(temp_arr, f_q, Ncoord());
  free_array(temp_arr);

  if (Opt_params.print_lvl >= 3) {
    oprintf_out("\tInternal forces in au, after projection of redundancies and constraints.\n");
    if (Opt_params.fb_fragments)
      oprintf_out("\tFB external coordinates are not projected.\n");
    oprint_array_out(f_q, Ncoord());
    //oprint_array_out_precise(f_q, Ncoord());
  }

  // Project redundances and constraints out of Hessian matrix
  // Peng, Ayala, Schlegel, JCC 1996 give H -> PHP + 1000(1-P)
  // The second term appears unnecessary and sometimes messes up Hessian updating.

  double **H = p_Opt_data->g_H_pointer();
  double **temp_mat = init_matrix(Nintco, Nintco);
  opt_matrix_mult(H, 0, P, 0, temp_mat, 0, Nintco, Nintco, Nintco, 0);
  opt_matrix_mult(P, 0, temp_mat, 0, H, 0, Nintco, Nintco, Nintco, 0);
  free_matrix(temp_mat);

  /*for (int i=0; i<Nintco;++i)
    H[i][i] += 1000 * (1.0 - P[i][i]);

  for (int i=0; i<Nintco; ++i)
    for (int j=0; j<i; ++j)
      H[j][i] = H[i][j] = H[i][j] + 1000 * (1.0 - P[i][j]);*/

  if (Opt_params.print_lvl >= 3) {
    oprintf_out("Projected (PHP) Hessian matrix\n");
    if (Opt_params.fb_fragments)
      oprintf_out("FB external coordinates are not projected.\n");
    oprint_matrix_out(H, Ncoord(), Ncoord());
  }
  free_matrix(P);

}

// project redundancies out of displacement vector; so far, doesn't seem to make much difference
void MOLECULE::project_dq(double *dq) {
  int Nintco = Ncoord();
  int Ncart = 3*g_natom();

  double *dq_orig; //only for printing
  if (Opt_params.print_lvl >=2) {
    dq_orig = init_array(Nintco);
    array_copy(dq, dq_orig, Ncoord());
  }

  double **B = compute_B();

  //double **G = compute_G(true);
  double **G = init_matrix(Ncart, Ncart);
  opt_matrix_mult(B, 1, B, 0, G, 0, Ncart, Nintco, Ncart, 0);

/*  will need fixed if this function ever helps
#if defined (OPTKING_PACKAGE_QCHEM)
  // Put 1's on diagonal for FB coordinates
  for (std::size_t I=0; I<fb_fragments.size(); ++I)
    for (int i=0; i<fb_fragments[i]->Ncoord(); ++i)
      G[g_fb_fragment_coord_offset(I) + i][g_fb_fragment_coord_offset(I) + i] = 1.0;
#endif
*/

  // B dx = dq
  // B^t B dx = B^t dq
  // dx = (B^t B)^-1 B^t dq
  double **G_inv = symm_matrix_inv(G, Ncart, 1);
  free_matrix(G);

  double **B_inv = init_matrix(Ncart, Nintco);
  opt_matrix_mult(G_inv, 0, B, 1, B_inv, 0, Ncart, Ncart, Nintco, 0);
  free_matrix(G_inv);

  double **P = init_matrix(Nintco, Nintco);
  opt_matrix_mult(B, 0, B_inv, 0, P, 0, Nintco, Ncart, Nintco, 0);
  free_matrix(B);
  free_matrix(B_inv);

  double * temp_arr = init_array(Nintco);
  opt_matrix_mult(P, 0, &dq, 1, &temp_arr, 1, Nintco, Nintco, 1, 0);
  array_copy(temp_arr, dq, Ncoord());
  free_array(temp_arr);
  free_matrix(P);

  if (Opt_params.print_lvl >=2) {
    oprintf_out("Projection of redundancies out of step:\n");
    oprintf_out("\tOriginal dq     Projected dq     Difference\n");
    for (int i=0; i<Nintco; ++i)
      oprintf_out("\t%12.6lf    %12.6lf   %12.6lf\n", dq_orig[i], dq[i], dq[i]-dq_orig[i]);
    free_array(dq_orig);
  }
}

void MOLECULE::apply_intrafragment_step_limit(double * & dq) {
  int dim = Ncoord();
  double scale = 1.0;
  double limit = Opt_params.intrafragment_step_limit;

  for (std::size_t f=0; f<fragments.size(); ++f)
    for (int i=0; i<fragments[f]->Ncoord(); ++i)
      if ((scale * sqrt(array_dot(dq,dq,dim)))> limit)
        scale = limit / sqrt(array_dot(dq,dq,dim));

  if (scale != 1.0) {
    oprintf_out("\tChange in coordinate exceeds step limit of %10.5lf.\n", limit);
    oprintf_out("\tScaling displacements by %10.5lf\n", scale);

    for (std::size_t f=0; f<fragments.size(); ++f)
      for (int i=0; i<fragments[f]->Ncoord(); ++i)
        dq[g_coord_offset(f)+i] *= scale;
  }

}

// Identify if some angles are passing through 0 or going to 180.
std::vector<int> MOLECULE::validate_angles(double const * const dq) {

  std::vector<int> lin_angle;
  std::vector<int> frag_angle;

  for (std::size_t f=0; f<fragments.size(); ++f) {
    frag_angle = fragments[f]->validate_angles(&(dq[g_coord_offset(f)]), g_atom_offset(f));

    for (std::size_t i=0; i<frag_angle.size(); ++i)
      lin_angle.push_back( frag_angle[i] );

    frag_angle.clear();
  }

  if (!lin_angle.empty()) {
    oprintf_out("\tNewly linear bends that need to be incoporated into the internal coordinates:\n");
    for (std::size_t i=0; i<lin_angle.size(); i+=3)
      oprintf_out("\t%5d%5d%5d\n", lin_angle[i]+1, lin_angle[i+1]+1, lin_angle[i+2]+1);
  }
  return lin_angle;
}

void MOLECULE::H_guess(void) const {
  double **H = p_Opt_data->g_H_pointer();

  if (Opt_params.intrafragment_H == OPT_PARAMS::SCHLEGEL)
    oprintf_out("\tGenerating empirical Hessian (Schlegel '84) for each fragment.\n");
  else if (Opt_params.intrafragment_H == OPT_PARAMS::FISCHER)
    oprintf_out("\tGenerating empirical Hessian (Fischer & Almlof '92) for each fragment.\n");
  else if (Opt_params.intrafragment_H == OPT_PARAMS::SIMPLE)
    oprintf_out("\tGenerating simple diagonal Hessian (.5 .2 .1) for each fragment.\n");
  else if (Opt_params.intrafragment_H == OPT_PARAMS::LINDH_SIMPLE)
    oprintf_out("\tGenerating diagonal Hessian from Lindh (1995) for each fragment.\n");
  else if (Opt_params.intrafragment_H == OPT_PARAMS::LINDH)
    oprintf_out("\tUsing model Hessian from Lindh (1995).\n");

  if (Opt_params.intrafragment_H == OPT_PARAMS::SCHLEGEL ||
      Opt_params.intrafragment_H == OPT_PARAMS::FISCHER  ||
      Opt_params.intrafragment_H == OPT_PARAMS::LINDH_SIMPLE  ||
      Opt_params.intrafragment_H == OPT_PARAMS::SIMPLE) {
    for (std::size_t f=0; f<fragments.size(); ++f) {
      double **H_frag = fragments[f]->H_guess();

      for (int i=0; i<fragments[f]->Ncoord(); ++i)
        for (int j=0; j<fragments[f]->Ncoord(); ++j)
          H[g_coord_offset(f) + i][g_coord_offset(f) + j] = H_frag[i][j];

      free_matrix(H_frag);
    }

    for (std::size_t I=0; I<interfragments.size(); ++I) {
      double **H_interfrag = interfragments[I]->H_guess();

      for (int i=0; i<interfragments[I]->Ncoord(); ++i)
        for (int j=0; j<interfragments[I]->Ncoord(); ++j)
          H[g_interfragment_coord_offset(I) + i][g_interfragment_coord_offset(I) + j] =
            H_interfrag[i][j];

      free_matrix(H_interfrag);
    }

    for (std::size_t I=0; I<fb_fragments.size(); ++I) {
      double **H_fb_frag = fb_fragments[I]->H_guess();

      for (int i=0; i<fb_fragments[I]->Ncoord(); ++i)
        for (int j=0; j<fb_fragments[I]->Ncoord(); ++j)
          H[g_fb_fragment_coord_offset(I) + i][g_fb_fragment_coord_offset(I) + j] = H_fb_frag[i][j];

      free_matrix(H_fb_frag);
    }
  }
  else if (Opt_params.intrafragment_H == OPT_PARAMS::LINDH) {
    double **H_xyz = Lindh_guess();      // generate Lindh cartesian Hessian
    //bool read_H_worked = cartesian_H_to_internals(H_xyz); // transform to internals
    cartesian_H_to_internals(H_xyz); // transform to internals
    // if fails, then what?  Fix later. TODO
    free_matrix(H_xyz);
  }

  if (Opt_params.print_lvl >= 2) {
    oprintf_out("\nInitial Hessian guess\n");
    oprint_matrix_out(H, Ncoord(), Ncoord());
    offlush_out();
  }

  return;
}

bool MOLECULE::cartesian_H_to_internals(double **H_cart) const {
  int Nintco = Ncoord();
  int Ncart = 3*g_natom();
  bool success = true; // to be dynamic later

  double **H_int = p_Opt_data->g_H_pointer();

  // If the "internals" are really cartesian, do nothing.
  if (Opt_params.coordinates == OPT_PARAMS::CARTESIAN && !is_noncart_present()) {
    opt_matrix_copy(H_cart, H_int, Ncart, Ncart);
    return true;
  }

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
  //double **H_cart = p_Opt_data->read_cartesian_H();

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
  //free_matrix(H_cart);

  //double **H_int = init_matrix(Nintco, Nintco);
  opt_matrix_mult(A, 1, temp_mat, 0, H_int, 0, Nintco, Ncart, Nintco, 0);
  free_matrix(temp_mat);

  free_matrix(A);

  if (Opt_params.print_lvl >= 3) {
    oprintf_out( "Hessian transformed to internal coordinates:\n");
    oprint_matrix_out(H_int, Nintco, Nintco);
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
  oprintf_out( "Hessian transformed back into Cartesian coordinates\n");
  oprint_matrix_out(H_cart, Ncart, Ncart);
  free_matrix(B);
*/
  return success;
}

double *MOLECULE::g_masses(void) const {
  double *u = init_array(g_natom());
  int cnt = 0;
  for (std::size_t f=0; f<fragments.size(); ++f)
    for (int i=0; i<fragments[f]->g_natom(); ++i)
      u[cnt++] = fragments[f]->g_mass(i);
  return u;
}

double *MOLECULE::g_Z(void) const {
  double *Zs = init_array(g_natom());
  int cnt = 0;
  for (std::size_t f=0; f<fragments.size(); ++f) {
    double *frag_Z = fragments[f]->g_Z_pointer();
    for (int i=0; i<fragments[f]->g_natom(); ++i)
      Zs[cnt++] = frag_Z[i];
  }
  return Zs;
}

// compute B matrix - leave rows for FB coordinates empty
double ** MOLECULE::compute_B(void) const {
  double **B = init_matrix(Ncoord(), 3*g_natom());

  for (std::size_t f=0; f<fragments.size(); ++f)
    fragments[f]->compute_B(B, g_coord_offset(f), g_atom_offset(f));

  for (std::size_t I=0; I<interfragments.size(); ++I) {
    int A_off = g_atom_offset( interfragments[I]->g_A_index());
    int B_off = g_atom_offset( interfragments[I]->g_B_index());
    int coord_off = g_interfragment_coord_offset(I);

    interfragments[I]->compute_B(B, coord_off, A_off, B_off);
  }
  return B;
}

double ** MOLECULE::compute_derivative_B(int intco_index) const {
  int cnt_intcos = 0;
  int fragment_index = -1;
  int coordinate_index = 0;
  bool is_interfragment = true;

  for (std::size_t f=0; f<fragments.size(); ++f) {
    for (int i=0; i<fragments[f]->Ncoord(); ++i) {
      if (cnt_intcos++ == intco_index) {
        fragment_index = f;
        coordinate_index = i;
        is_interfragment = false;
        break;
      }
    }
  }

  if (is_interfragment) {  // intco_index not yet found
    for (std::size_t f=0; f<interfragments.size(); ++f) {
      for (int i=0; i<interfragments[f]->Ncoord(); ++i) {
        if (cnt_intcos++ == intco_index) {
          fragment_index = f;
          coordinate_index = i;
          break;
        }
      }
    }
  }

  if (fragment_index == -1)
    throw(INTCO_EXCEPT("MOLECULE::compute_derivative_B() could not find intco_index"));

  double **dq2dx2_frag;

  double **dq2dx2 = init_matrix(3*g_natom(), 3*g_natom());

  if (!is_interfragment) {
    dq2dx2_frag = fragments[fragment_index]->compute_derivative_B(coordinate_index);
    // dimension is now the number of fragment atoms * 3 x the number of fragment atoms * 3

    int frag_natom = fragments[fragment_index]->g_natom();
    int atom_off = g_atom_offset(fragment_index);

    for (int a=0; a<frag_natom; ++a) {
      int atom_a = atom_off + a;

      for (int b=0; b<frag_natom; ++b) {
        int atom_b = atom_off + b;

      for (int xyz_a=0; xyz_a<3; ++xyz_a)
        for (int xyz_b=0; xyz_b<3; ++xyz_b)
          dq2dx2[3*atom_a + xyz_a][3*atom_b + xyz_b] = dq2dx2_frag[3*a+xyz_a][3*b+xyz_b];

/* old way before more complicated combinations
    // 2x2, 3x3 or 4x4 and given by natom_intco
    int natom_intco = fragments[fragment_index]->g_simple_natom(coordinate_index);
    int atom_a, atom_b;
    for (int a=0; a<natom_intco; ++a) {
      atom_a = g_atom_offset(fragment_index) + fragments[fragment_index]->g_simple_atom(coordinate_index, a);
      for (int b=0; b<natom_intco; ++b) {
        atom_b = g_atom_offset(fragment_index) + fragments[fragment_index]->g_simple_atom(coordinate_index, b);
      for (int xyz_a=0; xyz_a<3; ++xyz_a)
        for (int xyz_b=0; xyz_b<3; ++xyz_b)
          dq2dx2[3*atom_a + xyz_a][3*atom_b + xyz_b] = dq2dx2_frag[3*a+xyz_a][3*b+xyz_b];
      }
    }
*/
      }
    }
    free_matrix(dq2dx2_frag);
  }
/* TODO fix derivative B matrices for interfragment coordinates
  else { // interfragment cordinate
    dq2dx2_frag = interfragments[fragment_index]->compute_derivative_B(coordinate_index);

    // dimension is (3*natomA + 3*natomB) x (3*natomA + 3*natomB) -> 3*natom X 3*natom
    int nA = interfragments[fragment_index]->g_natom_A();
    int nB = interfragments[fragment_index]->g_natom_B();
    int iA = 3*g_atom_offset(interfragments[fragment_index]->g_A_index());
    int iB = 3*g_atom_offset(interfragments[fragment_index]->g_B_index());

    //oprint_matrix_out(dq2dx2_frag, 3*(nA+nB), 3*(nA+nB));
    for (int a=0; a<3*nA; ++a) // A-A block
      for (int aa=0; aa<3*nA; ++aa)
        dq2dx2[iA + a][iA + aa] = dq2dx2_frag[a][aa];

    for (int a=0; a<3*nA; ++a) // A-B block
      for (int bb=0; bb<3*nB; ++bb)
        dq2dx2[iA + a][iB + bb] = dq2dx2_frag[a][3*nA + bb];

    for (int b=0; b<3*nB; ++b) // B-A block
      for (int aa=0; aa<3*nA; ++aa)
        dq2dx2[iB + b][iA + aa] = dq2dx2_frag[3*nA + b][aa];

    for (int b=0; b<3*nB; ++b) // B-B block
      for (int bb=0; bb<3*nB; ++bb)
        dq2dx2[iB + b][iB + bb] = dq2dx2_frag[3*nA + b][3*nA + bb];

    free_matrix(dq2dx2_frag);
  }
*/
  return dq2dx2;
}

double ** MOLECULE::compute_G(bool use_masses) const {
  int Nintco = Ncoord();
  int Ncart = 3*g_natom();

  double **B = compute_B();
  double **G = init_matrix(Nintco, Nintco);

  if (use_masses) {
    double *u = g_masses();

    for (int i=0; i<Nintco; ++i)
      for (int a=0; a<g_natom(); ++a)
        for(int xyz=0; xyz<3; ++xyz)
          B[i][3*a+xyz] /= sqrt(u[a]);

    free_array(u);
  }

  opt_matrix_mult(B, 0, B, 1, G, 0, Nintco, Ncart, Nintco, 0);
  free_matrix(B);

  //oprintf_out("G matrix\n");
  //oprint_matrix_out(G, Ncoord(), Ncoord());

  return G;
}

// Apply strings of atoms for frozen and fixed coordinates;
bool MOLECULE::apply_input_constraints(void) {
  bool frozen_present = false;
  bool fixed_present = false;

  if (   !Opt_params.frozen_distance_str.empty()
      || !Opt_params.frozen_bend_str.empty()
      || !Opt_params.frozen_dihedral_str.empty()
      || !Opt_params.frozen_cartesian_str.empty() ) {
    oprintf_out("\tAssuming in current code that numbering for constraints corresponds to unified fragment.\n");
    frozen_present = fragments[0]->apply_frozen_constraints(Opt_params.frozen_distance_str,
      Opt_params.frozen_bend_str, Opt_params.frozen_dihedral_str, Opt_params.frozen_cartesian_str);
  }

  if (   !Opt_params.fixed_distance_str.empty()
      || !Opt_params.fixed_bend_str.empty()
      || !Opt_params.fixed_dihedral_str.empty() ) {
    oprintf_out("\tAssuming in current code that numbering for constraints corresponds to unified fragment.\n");
    fixed_present = fragments[0]->apply_fixed_constraints(Opt_params.fixed_distance_str,
      Opt_params.fixed_bend_str, Opt_params.fixed_dihedral_str);
  }

  return (fixed_present || frozen_present);
}

// Add cartesian coordinates
int MOLECULE::add_cartesians(void) {
  int nadded = 0;
  for (std::size_t f=0; f<fragments.size(); ++f)
    nadded += fragments[f]->add_cartesians();
  return nadded;
}

// freeze all fragments in molecule
void MOLECULE::freeze_intrafragments(void) {
  oprintf_out("\tSetting all fragments to frozen.\n");
  for (std::size_t f=0; f<fragments.size(); ++f)
    fragments[f]->freeze();
}

// Freeze all coordinates within fragments
void MOLECULE::freeze_intrafragment_coords(void) {
  oprintf_out("\tSetting all coordinates within each fragment to frozen.\n");
  for (std::size_t f=0; f<fragments.size(); ++f)
    fragments[f]->freeze_coords();
}

// Determine trivial coordinate combinations, i.e., don't combine..
int MOLECULE::form_trivial_coord_combinations(void) {
  int nadded = 0;
  for (std::size_t f=0; f<fragments.size(); ++f)
    nadded += fragments[f]->form_trivial_coord_combinations();
  for (std::size_t I=0; I<interfragments.size(); ++I)
    nadded += interfragments[I]->form_trivial_coord_combinations();
  return nadded;
}

// Determine initial delocalized coordinate coefficients.
int MOLECULE::form_delocalized_coord_combinations(void) {
  int nadded = 0;
  for (std::size_t f=0; f<fragments.size(); ++f)
    nadded += fragments[f]->form_delocalized_coord_combinations();

  // We can throw out coordinates which are asymmetric wrt to the ENTIRE system.
  if (g_nfragment() == 2 && g_nfb_fragment() == 0) { // there is only 1 fragment

    // Determine asymmetric combinations; Check each coordinate = row of B.
    int Natom = g_natom();
    double **Bs = fragments[0]->compute_B();
    double **orig_geom = g_geom_2D();
    std::vector<int> asymm_coord;

    // Construct each displaced geometry, then test it.
    for (int cc=0; cc<Ncoord(); ++cc) {
      double **displaced_geom = matrix_return_copy(orig_geom, Natom, 3);
      for (int a=0; a<Natom; ++a)
        for (int xyz=0; xyz<3; ++xyz)
          displaced_geom[a][xyz] += 0.1 * Bs[cc][3*a+xyz];

      bool symm_rfo_step = false;
#if defined(OPTKING_PACKAGE_PSI)
      psi::Process::environment.legacy_molecule()->set_geometry(displaced_geom);
      symm_rfo_step = psi::Process::environment.legacy_molecule()->valid_atom_map(Opt_params.symm_tol);
      psi::Process::environment.legacy_molecule()->set_geometry(orig_geom);
#elif defined(OPTKING_PACKAGE_QCHEM)
      // TODO QCHEM
      symm_rfo_step = true;
#endif
      free_matrix(displaced_geom);

      if (!symm_rfo_step)
        asymm_coord.push_back(cc);
    }

    // Remove asymmetric combinations
    nadded -= asymm_coord.size();
    if (asymm_coord.size())
      oprintf_out("\tRemoving the following coordinates because they break molecular point group:\n\t");
    for (std::size_t i=0; i<asymm_coord.size(); ++i)
      oprintf_out(" %d", asymm_coord[i]+1);
    oprintf_out("\n");

    for (std::size_t i=0; i<asymm_coord.size(); ++i) {
      fragments[0]->erase_combo_coord(asymm_coord[i]);
      for (std::size_t j=i; j<asymm_coord.size(); ++j)
        asymm_coord[j] -= 1;
    }
    free_matrix(Bs);
    free_matrix(orig_geom);
    asymm_coord.clear();
  }
  oprintf_out("\tA total of %d delocalized coordinates added.\n\n", nadded);
  return nadded;
}


// Determine Pulay natural coordinate combinations.
int MOLECULE::form_natural_coord_combinations(void) {
  int nadded = 0;
  for (std::size_t f=0; f<fragments.size(); ++f)
    nadded += fragments[f]->form_natural_coord_combinations();
  return nadded;
}

// The values of the FB coordinates will be taken to be the total change
// since the beginning of the optimization.  These values are determined

// Compute constraint matrix.
double ** MOLECULE::compute_constraints(void) {
  double **C, **C_frag, **C_inter;
  int i, j;

  C = init_matrix(Ncoord(), Ncoord());

  for (std::size_t f=0; f<fragments.size(); ++f) {
    C_frag = fragments[f]->compute_constraints();

    for (i=0; i<fragments[f]->Ncoord(); ++i)
      for (j=0; j<fragments[f]->Ncoord(); ++j)
        C[g_coord_offset(f)+i][g_coord_offset(f)+j] = C_frag[i][j];

    free_matrix(C_frag);
  }

  for (std::size_t I=0; I<interfragments.size(); ++I) {
    C_inter = interfragments[I]->compute_constraints(); // Ncoord() X Ncoord

    for (i=0; i<interfragments[I]->Ncoord(); ++i)
      for (j=0; j<interfragments[I]->Ncoord(); ++j)
        C[g_interfragment_coord_offset(I)+i][g_interfragment_coord_offset(I)+j] = C_inter[i][j];

    free_matrix(C_inter);
  }

  if (Opt_params.print_lvl >= 3) {
    oprintf_out("Constraint matrix\n");
    oprint_matrix_out(C, Ncoord(), Ncoord());
  }

  return C;
}

// Test a linear combination of internal coordinate displacements to see if it
// breaks the molecular point group.
bool MOLECULE::coord_combo_is_symmetric(double *intco_combo, int dim) {
  int natom = g_natom();
  double norm = array_norm(intco_combo, dim);

  double **B = compute_B();
  double **orig_geom = g_geom_2D();

  double **displaced_geom = matrix_return_copy(orig_geom, natom, 3);
  for (int xyz=0; xyz<3; ++xyz)
    for (int atom=0; atom<natom; ++atom)
      for (int intco=0; intco<dim; ++intco)
        displaced_geom[atom][xyz] += 0.1/norm * intco_combo[intco] * B[intco][3*atom+xyz];

  //oprintf_out("Displaced geometry in coord_combo_is_symmetric:\n");
  //oprint_matrix_out(displaced_geom, natom, 3);

  bool symm_rfo_step = false;
#if defined(OPTKING_PACKAGE_PSI)
  psi::Process::environment.legacy_molecule()->set_geometry(displaced_geom);
  symm_rfo_step = psi::Process::environment.legacy_molecule()->valid_atom_map(Opt_params.symm_tol);
  psi::Process::environment.legacy_molecule()->set_geometry(orig_geom);
#elif defined(OPTKING_PACKAGE_QCHEM)
  // TODO QCHEM
  symm_rfo_step = true;
#endif
  free_matrix(displaced_geom);

  if (symm_rfo_step)
    return true;
  else
    return false;
}

bool MOLECULE::is_coord_fixed(int coord_index) {
  int cnt = 0;
  for (std::size_t f=0; f<fragments.size(); ++f)
    for  (int i=0; i<fragments[f]->Ncoord(); ++i) {
      if (cnt == coord_index)
        return fragments[f]->coord_has_fixed_eq_val(i);
      ++cnt;
    }
    return false;
}

// Add dummy FB fragment which contains no atoms.  Read in the energy
// and the forces from QChem.  This will only work (maybe:) for QChem
void MOLECULE::add_fb_fragments(void) {

#if defined(OPTKING_PACKAGE_QCHEM)
  // get number of FB fragments
  int num_fb_frags = ::EFP::GetInstance()->NFragments();
  oprintf_out("\tAdding %d EFP fragments.\n", num_fb_frags);

  // get energy
  energy = ::EFP::GetInstance()->GetEnergy();

  FB_FRAG *one_frag;

  for (int i=0; i<num_fb_frags; ++i) {
    one_frag = new FB_FRAG();
    // add 6 simples just to act as placeholders ;
    // we read external forces and don't compute B matrix for these
    one_frag->add_dummy_coords(6);

    // get gradient
    double *g = init_array(6);
    ::EFP::GetInstance()->GetGrad(i,g);
    one_frag->set_forces(g);
    free_array(g);

    // See note below on EFP values.

    fb_fragments.push_back(one_frag);
  }
#elif defined(OPTKING_PACKAGE_PSI)

  return; // TODO

#endif

}

// from the data in opt_data after that data is read.
void MOLECULE::update_fb_values(void) {

  for (std::size_t i=0; i<fb_fragments.size(); ++i) {
    double *vals = init_array(6);

    double *dq;
    for (int iter=0; iter<p_Opt_data->g_iteration(); ++iter) {
      dq = p_Opt_data->g_dq_pointer(iter);
      for (int coord=0; coord < 6; ++coord)
        vals[coord] += dq[g_fb_fragment_coord_offset(i)+coord] ;
    }
    fb_fragments[i]->set_values(vals);
    free_array(vals);
  }
}

}
