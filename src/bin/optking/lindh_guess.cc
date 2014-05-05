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

/*! \file Lindh_guess.cc
    \ingroup optking
    \brief Function to generate a model cartesian Hessian according to
     R. Lindh, A. Bernhardsson, G. Karlstrom, P.-A. Malmqvist, CPL, 241, 423, 1995.
*/

#include "molecule.h"

#include "v3d.h"

#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
#endif

namespace opt {

using namespace v3d;

// Returns cartesian Lindh guess Hessian for whole system
double **MOLECULE::Lindh_guess(void) const {

  // Build one temporary fragment that contains ALL the atoms.
  int natom = g_natom();
  double **coord_xyz = g_geom_2D();
  double *atomic_numbers = g_Z();

  FRAG * frag = new FRAG(natom, atomic_numbers, coord_xyz);

  double **H_xyz = frag->Lindh_guess();

  delete frag;
  return H_xyz;
}

// return period from atomic number
inline int period(int Z) {
  if      (Z <=  2) return 1;
  else if (Z <= 10) return 2;
  else if (Z <= 18) return 3;
  else if (Z <= 36) return 4;
  else              return 5;
}

// return Lindh alpha value from two periods
inline double alpha_table(int perA, int perB) {
  if (perA == 1) {
    if (perB == 1)
      return 1.000;
    else
      return 0.3949;
  }
  else {
    if (perB == 1)
      return 0.3949;
    else
      return 0.2800;
  }
}

inline double r_ref_table_orig(int perA, int perB) {
  if (perA == 1) {
    if (perB == 1) return 1.35;
    else if (perB == 2) return 2.10;
    else return 2.53;
  }
  else if (perA == 2) {
    if (perB == 1) return 2.10;
    else if (perB == 2) return 2.87;
    else return 3.40;
  }
  else {
    if (perB == 1) return 2.53;
    else return 3.40;
  }
}

// rho_ij = e^(alpha (r^2,ref - r^2))
double FRAG::Lindh_rho_orig(int A, int B, double RAB) const {

  int perA = period((int) Z[A]);
  int perB = period((int) Z[B]);

  double alpha = alpha_table(perA, perB);
  double r_ref = r_ref_table_orig(perA, perB);

  return exp(-alpha * (RAB*RAB - r_ref*r_ref));
}

// Build cartesian hessian according to model in  Lindh paper.
double ** FRAG::Lindh_guess(void) {

  // Build distance matrix
  double **R = init_matrix(natom, natom);
  for (int A=0; A<natom; ++A)
    for (int B=0; B<=A; ++B)
      R[B][A] = R[A][B] = v3d_dist(geom[A], geom[B]);

  // Define "close" atoms;
  double const dist_limit = 5.5;
  bool **close = init_bool_matrix(natom, natom);
  for (int A=0; A<natom; ++A)
    for (int B=0; B<natom; ++B)
      close[A][B] = true;
      //if (R[A][B] < dist_limit)
// make all close for now

    const double k_r   = 0.45;
    const double k_phi = 0.15;
    const double k_tau = 0.005;
    double Lindh_k;

    // Model is defined including ALL possible bends, angles, and torsions..
    double **Hx = init_matrix(3*natom,3*natom);

    // loop over all possible stretches, between i j in fragment
    for (int i=0; i<natom; ++i)
      for (int j=i+1; j<natom; ++j) {

        Lindh_k = k_r * Lindh_rho_orig(i, j, R[i][j]);

        STRE *s1 = new STRE(i,j);
        double **sB = s1->DqDx(geom);

        // H_x1_x2 = dq_i/dq_x1 d2E/dq_i^2 dq_i/dq_x2 
        for (int a=0; a < s1->g_natom(); ++a)
          for (int xyz_a=0; xyz_a<3; ++xyz_a)
            for (int b=0; b < s1->g_natom(); ++b)
              for (int xyz_b=0; xyz_b<3; ++xyz_b)
                Hx[3*s1->g_atom(a)+xyz_a][3*s1->g_atom(b)+xyz_b] +=
                  Lindh_k * sB[a][xyz_a] * sB[b][xyz_b];

        delete s1;
      }

    // loop over all possible i j k
    for (int i=0; i<natom; ++i)
      for (int j=0; j<natom; ++j)
        if (j != i)
          for (int k=i+1; k<natom; ++k)
            if (k != j) {

              Lindh_k = k_phi * Lindh_rho_orig(i, j, R[i][j]) *
                                Lindh_rho_orig(j, k, R[j][k]);

              BEND *b1 = new BEND(i,j,k);
              double **bB = b1->DqDx(geom);

              // H_x1_x2 = dq_i/dq_x1 d2E/dq_i^2 dq_i/dq_x2 
              for (int a=0; a<b1->g_natom(); ++a)
                for (int xyz_a=0; xyz_a<3; ++xyz_a)
                  for (int b=0; b<b1->g_natom(); ++b)
                    for (int xyz_b=0; xyz_b<3; ++xyz_b)
                      Hx[3*b1->g_atom(a)+xyz_a][3*b1->g_atom(b)+xyz_b] +=
                        Lindh_k * bB[a][xyz_a] * bB[b][xyz_b];

              delete b1;
            }

    for (int i=0; i<natom; ++i)
      for (int j=0; j<natom; ++j)
        if (j != i)
          for (int k=0; k<natom; ++k)
            if ( k!=i && k!=j)
              for (int l=i+1; l<natom; ++l)
                if ( l!=j && l!=k) {

                  Lindh_k = k_tau * Lindh_rho_orig(i, j, R[i][j]) *
                                    Lindh_rho_orig(j, k, R[j][k]) * 
                                    Lindh_rho_orig(k, l, R[k][l]);

                  TORS *t1 = new TORS(i,j,k,l);
                  double **tB = t1->DqDx(geom);

                  // H_x1_x2 = dq_i/dq_x1 d2E/dq_i^2 dq_i/dq_x2 
                  for (int a=0; a<t1->g_natom(); ++a)
                    for (int xyz_a=0; xyz_a<3; ++xyz_a)
                      for (int b=0; b<t1->g_natom(); ++b)
                        for (int xyz_b=0; xyz_b<3; ++xyz_b)
                          Hx[3*t1->g_atom(a)+xyz_a][3*t1->g_atom(b)+xyz_b] +=
                            Lindh_k * tB[a][xyz_a] * tB[b][xyz_b];

                  delete t1;
            }

/* Use this code to build all internals and compute forces;  otherwise just do one intco at a time
  // Generate coordinates between ALL atoms, but not 
  // complete reversals or using duplicate i,j,k,l.
  // Create stretch coordinates.
  for (int i=0; i<natom; ++i) {
    for (int j=i+1; j<natom; ++j) {
      if (close[i][j]) {
        STRE *one_stre = new STRE(i, j);
        intcos.push_back(one_stre);
      }

  // Create bend coodinates.
  for (int i=0; i<natom; ++i)
    for (int j=0; j<natom; ++j)
      if (close[j][i] && (j != i))
        for (int k=i+1; k<natom; ++k)
          if (close[k][j] && close[k][i] && (k != j)) {
            BEND *one_bend = new BEND(i, j, k);
            intcos.push_back(one_bend);
          }

  // Create torsion coodinates.
  for (int i=0; i<natom; ++i)
    for (int j=0; j<natom; ++j)
      if (close[j][i] && (j != i))
        for (int k=0; k<natom; ++k)
          if (close[k][j] && close[k][i] && (k!=i) && (k!=j))
            for (int l=i+1; l<natom; ++l)
              if (close[l][k] && close[l][j] && close[l][i] && l!=j && l!=k) {
                TORS *one_tors = new TORS(i, j, k, l);
                intcos.push_back(one_tors);
              }

  free_bool_matrix(close);

  const double k_r   = 0.45;
  const double k_phi = 0.15;
  const double k_tau = 0.005;
  double k;

  double **Hx = init_matrix(3*natom, 3*natom);


  for (int i=0; i<intcos.size(); ++i) {  // loop over intcos
    SIMPLE * q = intcos.at(i);

    double **Bintco = q->DqDx(geom); // dq_i / da_xyz
    int natom_intco = q->g_natom();

    if (q->g_type() == stre_type) {
      int a = q->g_atom(0);
      int b = q->g_atom(1);
      k = k_r * Lindh_rho_orig(a, b, R[a][b]);
    }
    else if (intcos.at(i)->g_type() == bend_type) {
      int a = q->g_atom(0);
      int b = q->g_atom(1);
      int c = q->g_atom(2);
      k = k_phi * Lindh_rho_orig(a, b, R[a][b])
                * Lindh_rho_orig(b, c, R[b][c]);
    }
    else if (intcos.at(i)->g_type() == tors_type) {
      int a = q->g_atom(0);
      int b = q->g_atom(1);
      int c = q->g_atom(2);
      int d = q->g_atom(3);
      k = k_tau * Lindh_rho_orig(a, b, R[a][b])
                * Lindh_rho_orig(b, c, R[b][c])
                * Lindh_rho_orig(c, d, R[c][d]);
    }

    // Hxy += k dq/dx dq/dy
    for (int a=0; a < natom_intco; ++a) {
      int x1 = q->g_atom(a);

      for (int b=0; b < natom_intco; ++b) {
        int x2 = q->g_atom(b);

        for (int xyz1=0; xyz1<3; ++xyz1)
          for (int xyz2=0; xyz2<3; ++xyz2)
             Hx[3*x1 + xyz1][3*x2 + xyz2] += k * Bintco[a][xyz1] * Bintco[b][xyz2];
      }
    }
    free_matrix(Bintco);

// If we include energy gradient term, then need gradient in internals.  // Is it worth it?

  // Compute g_q = (BB^t)^-1 B g_x
  long int Nintco = intcos.size();
  double **B = compute_B();
  double *g_x = g_grad_array();
  double *temp_arr = init_array(Nintco);
  opt_matrix_mult(B, 0, &g_x, 1, &temp_arr, 1, Nintco, 3*natom, 1, 0);
  free_array(g_x);

  double **G = init_matrix(Nintco, Nintco);
  for (int i=0; i<Nintco; ++i)
    for (int k=0; k<3*natom; ++k)
      for (int j=0; j<Nintco; ++j)
        G[i][j] += B[i][k] * B[j][k];
  free_matrix(B);
  double **G_inv = symm_matrix_inv(G, Nintco, 1);
  free_matrix(G);

  double *g_q = init_array(Nintco);
  opt_matrix_mult(G_inv, 0, &temp_arr, 1, &g_q, 1, Nintco, Nintco, 1, 0);
  free_matrix(G_inv);
  free_array(temp_arr);
  // Done computing g_q
    // Hxy += dE/dq_i d2q_i/dxdy
    double **Dq2 = q->Dq2Dx2(geom);
    for (int a=0; a < natom_intco; ++a) {
      int x1 = q->g_atom(a);
      for (int b=0; b < natom_intco; ++b) {
        int x2 = q->g_atom(b);
        for (int xyz1=0; xyz1<3; ++xyz1)
          for (int xyz2=0; xyz2<3; ++xyz2)
             Hx[3*x1 + xyz1][3*x2 + xyz2] += k * g_q[i] * Dq2[3*a+xyz1][3*b+xyz2];
      }
    }
    free_matrix(Dq2);
  //free_array(g_q);
*/

  free_matrix(R);

  if (Opt_params.print_lvl >= 2) {
    fprintf(outfile,"Lindh cartesian Hessian guess\n");
    print_matrix(outfile, Hx, 3*natom, 3*natom);
    fflush(outfile);
  }
  return Hx;
}

}
