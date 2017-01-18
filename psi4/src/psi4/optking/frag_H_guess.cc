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

/*! \file    H_guess.cc
    \ingroup optking
    \brief   generates diagonal empirical Hessians in a.u. such as
      Schlegel, Theor. Chim. Acta, 66, 333 (1984) and
      Fischer and Almlof, J. Phys. Chem., 96, 9770 (1992).
*/

#include "frag.h"
#include "cov_radii.h"
#include "psi4/optking/physconst.h"
#include "v3d.h"
#include "psi4/psi4-dec.h"
#include "print.h"
#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
#endif

namespace opt {

using namespace v3d;

// return period from atomic number
static inline int period(int Z) {
  if      (Z <=  2) return 1;
  else if (Z <= 10) return 2;
  else if (Z <= 18) return 3;
  else if (Z <= 36) return 4;
  else              return 5;
}

// return generic distance from two periods;
// based on DZP RHF, I have suggested: 1.38 1.9 2.53
// my values certainly work better for water
static inline double r_ref_table(int perA, int perB) {
  if (perA == 1) {
    if (perB == 1) return 1.35;     // Lindh 1.35
    else if (perB == 2) return 2.1; // Lindh 2.1
    else return 2.53;
  }
  else if (perA == 2) {
// based on DZP RHF, I have suggested: 1.9 2.87 3.40
    if (perB == 1) return 2.1;      // Lindh 2.1
    else if (perB == 2) return 2.87;
    else return 3.40;
  }
  else {
    if (perB == 1) return 2.53;
    else return 3.40;
  }
}

// return Lindh alpha value from two periods
static inline double alpha_table(int perA, int perB) {
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

// rho_ij = e^(alpha (r^2,ref - r^2))
double FRAG::Lindh_rho(int A, int B, double RAB) const {

  int perA = period((int) Z[A]);
  int perB = period((int) Z[B]);

  double alpha = alpha_table(perA, perB);
  double r_ref = r_ref_table(perA, perB);

  return exp(-alpha * (RAB*RAB - r_ref*r_ref));
}

// covalent bond length in bohr from atomic numbers
static inline double Rcov(double ZA, double ZB) {
  return (cov_radii[(int) ZA] + cov_radii[(int) ZB]) / _bohr2angstroms;
}

// This function generates various diagonal Hessian guesses.
double ** FRAG::H_guess(void) {

  double **R = init_matrix(natom, natom);
  for (int i=0; i<natom; ++i)
    for (int j=i+1; j<natom; ++j)
      R[i][j] = R[j][i] = v3d_dist(geom[i], geom[j]);

  // to hold diagonal force constant guesses for simples
  double *f = init_array(coords.simples.size());

  int cnt = 0;

  // Form diagonal Hessian in simple internals
  if (Opt_params.intrafragment_H == OPT_PARAMS::SCHLEGEL) {
    for (std::size_t i=0; i<coords.simples.size(); ++i) {
      SIMPLE_COORDINATE *q = coords.simples.at(i);

      int a,b,c;
      double A,B,rBCcov;

      switch (q->g_type()) {

        case (stre_type) :
          if (q->is_hbond()) f[cnt++] = 0.03;
          else {
            a = q->g_atom(0);
            b = q->g_atom(1);
            int perA = period((int) Z[a]);
            int perB = period((int) Z[b]);

            if ( perA==1 && perB==1) B = -0.244;
            else if ((perA==1 && perB==2) || (perB==1 && perA==2)) B = 0.352;
            else if (perA==2 && perB==2)                           B = 1.085;
            else if ((perA==1 && perB==3) || (perB==1 && perA==3)) B = 0.660;
            else if ((perA==2 && perB==3) || (perB==2 && perA==3)) B = 1.522;
            else B = 2.068;

            A = 1.734;
            // force constants in au / bohr^2
            f[cnt++] = A/((R[a][b]-B)*(R[a][b]-B)*(R[a][b]-B));
            // for aJ/Ang^2, * _hartree2aJ / SQR(_bohr2angstroms);
          }
        break;

        case (bend_type) :
          a = q->g_atom(0);
          b = q->g_atom(1);
          c = q->g_atom(2);

          f[cnt++] = ( ( (((int) Z[a]) == 1) || (((int) Z[c]) == 1) ) ? 0.160 : 0.250);
        break;

        case (tors_type) :
          b = q->g_atom(1);
          c = q->g_atom(2);
          A = 0.0023;
          B = 0.07;
          rBCcov = Rcov(Z[b], Z[(c)]);
          if (R[b][c] > (rBCcov + A/B)) B = 0.0; // keep > 0
          f[cnt++] = (A - (B*(R[b][c] - rBCcov)));
        break;

        case (oofp_type) :
          f[cnt++] = 0.045;
        break;

        case (cart_type) :
          f[cnt++] = 0.1;
        break;

        default:
          oprintf_out("H_guess encountered unknown internal type.\n");
          f[cnt++] = 1.0;
      } // end switch coordinate type
    } // loop over coords.simples
  } // end Schlegel
  else if (Opt_params.intrafragment_H == OPT_PARAMS::FISCHER) {
    for (std::size_t i=0; i<coords.simples.size(); ++i) {
      SIMPLE_COORDINATE *q = coords.simples.at(i);

      int a,b,c,d,L;
      double A,B,C,D,E,rABcov,rBCcov,rBDcov;

      switch (q->g_type()) {

        case (stre_type) :
          if (q->is_hbond()) f[cnt++] = 0.03;
          else {
            a = q->g_atom(0);
            b = q->g_atom(1);

            rABcov = Rcov(Z[a], Z[b]);

            A = 0.3601;
            B = 1.944;
            f[cnt++] = A * exp(-B*(R[a][b] - rABcov));
          }
        break;

        case (bend_type) :
          a = q->g_atom(0);
          b = q->g_atom(1);
          c = q->g_atom(2);

          rABcov = Rcov(Z[a], Z[b]);
          rBCcov = Rcov(Z[b], Z[c]);

          A = 0.089; B =  0.11; C =  0.44; D = -0.42;
          f[cnt++] = A + B/(pow(rABcov*rBCcov, D)) *
            exp(-C*( R[a][b] + R[b][c] - rABcov - rBCcov));
        break;

        case (tors_type) :
          b = q->g_atom(1);
          c = q->g_atom(2);

          rBCcov = Rcov(Z[b], Z[c]);

          // count number of additional bonds on central atoms
          L = 0;
          for (int j=0; j<natom; ++j) {
            if (j == c) continue;
            if (connectivity[b][j]) ++L;
          }
          for (int j=0; j<natom; ++j) {
            if (j == b) continue;
            if (connectivity[c][j]) ++L;
          }
          A = 0.0015; B = 14.0; C = 2.85; D = 0.57; E = 4.00;
          f[cnt++] = A + B * pow(L,D) / pow(R[b][c] * rBCcov, E) * exp(-C * (R[b][c] - rBCcov));
        break;

        case (oofp_type) :
          a = q->g_atom(0);
          b = q->g_atom(1);
          c = q->g_atom(2);
          d = q->g_atom(3);

          rABcov = Rcov(Z[a], Z[b]);
          rBCcov = Rcov(Z[b], Z[c]);
          rBDcov = Rcov(Z[b], Z[d]);

          double phi; // compute value of angle
          if (!v3d_oofp(geom[a], geom[b], geom[c], geom[d], phi))
            phi = _pi/4;

          A = 0.0025; B = 0.0061; C = 3.0; D = 4.0; E = 0.8;
          f[cnt++] = A + B * pow(rBCcov*rBDcov, E) * pow(cos(phi), D) * exp(-C * (R[a][b] - rABcov));
        break;

        case (cart_type) :
          f[cnt++] = 0.1;
        break;

        default:
          oprintf_out("H_guess encountered unknown internal type.\n");
          f[cnt++] = 1.0;
      } // end switch simples
    } // end loop simples
  } // end Fischer
  else if (Opt_params.intrafragment_H == OPT_PARAMS::SIMPLE) {
    for (std::size_t i=0; i<coords.simples.size(); ++i) {
      SIMPLE_COORDINATE *q = coords.simples.at(i);
      switch (q->g_type()) {
        case (stre_type) :
          f[cnt++] = 0.5;
        break;

        case (bend_type) :
          f[cnt++] = 0.2;
        break;

        case (tors_type) :
          f[cnt++] = 0.1;
        break;

        case (oofp_type) :
          f[cnt++] = 0.1;
        break;

        case (cart_type) :
          f[cnt++] = 0.1;
        break;

        default:
          oprintf_out("H_guess encountered unknown internal type.\n");
          f[cnt++] = 1.0;
      }
    }
  }
  else if (Opt_params.intrafragment_H == OPT_PARAMS::LINDH_SIMPLE) {

    const double k_r   = 0.45;
    const double k_phi = 0.15;
    const double k_tau = 0.005;

    for (std::size_t i=0; i<coords.simples.size(); ++i) {
      SIMPLE_COORDINATE * q = coords.simples.at(i);

      int a,b,c,d;

      switch (q->g_type()) {

        case (stre_type) :
          a = q->g_atom(0);
          b = q->g_atom(1);
          f[cnt++] = k_r * Lindh_rho(a, b, R[a][b]);
        break;

        case (bend_type) :
          a = q->g_atom(0);
          b = q->g_atom(1);
          c = q->g_atom(2);
          f[cnt++] = k_phi * Lindh_rho(a, b, R[a][b]) * Lindh_rho(b, c, R[b][c]);
        break;

        case (tors_type) :
          a = q->g_atom(0);
          b = q->g_atom(1);
          c = q->g_atom(2);
          d = q->g_atom(3);
          f[cnt++] = k_tau * Lindh_rho(a, b, R[a][b]) * Lindh_rho(b, c, R[b][c]) * Lindh_rho(c, d, R[c][d]);
        break;

        case (oofp_type) :
          f[cnt++] = 0.1;
        break;

        case (cart_type) :
          f[cnt++] = 0.1;
        break;

        default:
          oprintf_out("H_guess encountered unknown internal type.\n");
          f[cnt++] = 1.0;
      }
    }
  }
  else {
    oprintf_out("FRAG::H_guess(): Unknown Hessian guess type.\n");
  }

  free_matrix(R);

  if (Opt_params.print_lvl > 1) {
    oprintf_out("diagonal Hessian values for simple coordinates.\n");
    oprint_array_out(f, coords.simples.size());
  }

  // all off-diagonal entries are zero, so this seem silly
  double **H_simple = init_matrix(coords.simples.size(), coords.simples.size());
  for (std::size_t i=0; i<coords.simples.size(); ++i)
    H_simple[i][i] = f[i];
  free_array(f);

  double **H = coords.transform_simples_to_combo(H_simple);
  free_matrix(H_simple);

  return H;
}

}
