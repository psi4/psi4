/*! \file    H_guess.cc
    \ingroup optking
    \brief   generates empirical Hessian according to
      Schlegel, Theor. Chim. Acta, 66, 333 (1984) or
      Fischer and Almlof, J. Phys. Chem., 96, 9770 (1992).
      currently returned in atomic units
*/

#include "frag.h"
#include "cov_radii.h"
#include "physconst.h"
#include <cmath>
#include "v3d.h"
#define EXTERN
#include "globals.h"

namespace opt {

using namespace v3d;

// return period from atomic number
int period(int Z) { 
  if      (Z <=  2) return 1;
  else if (Z <= 10) return 2;
  else if (Z <= 18) return 3;
  else if (Z <= 36) return 4;
  else              return 5;
}

// return generic distance from two periods; I modifed based on DZP RHF
inline double r_ref_table(int perA, int perB) {
  if (perA == 1) {
    if (perB == 1) return 1.38;     // Lindh used 1.35;
    else if (perB == 2) return 1.9; // Lindh used 2.10; O-H=1.78; C-H=2.05
    else return 2.53;
  }
  else if (perA == 2) {
    if (perB == 1) return 1.9;      // Lindh used 2.10; O-H=1.78; C-H=2.05
    else if (perB == 2) return 2.87;
    else return 3.40;
  }
  else {
    if (perB == 1) return 2.53;
    else return 3.40;
  }
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

// rho_ij = e^(alpha (r^2,ref - r^2))
double FRAG::Lindh_rho(int A, int B, double RAB) const {

  int perA = period((int) Z[A]);
  int perB = period((int) Z[B]);

  double alpha = alpha_table(perA, perB);
  double r_ref = r_ref_table(perA, perB);

  return exp(-alpha * (RAB*RAB - r_ref*r_ref));
}

// covalent bond length in bohr from atomic numbers
inline double Rcov(double ZA, double ZB) {
  return (cov_radii[(int) ZA] + cov_radii[(int) ZB]) / _bohr2angstroms;
}

double ** FRAG::H_guess(void) {
  int i, j, a, b, c, d, cnt = 0, perA, perB;
  double rABcov, rBCcov, rBDcov;
  double A, B, C, D, L, E, val;

  double **R = init_matrix(natom, natom);
  for (int A=0; A<natom; ++A)
    for (int B=0; B<natom; ++B)
      R[A][B] = v3d_dist(geom[A], geom[B]);

  double *f = init_array(intcos.size());

  // Form diagonal Hessian in simple internals
  if (Opt_params.intrafragment_H == OPT_PARAMS::SCHLEGEL) {
    for (i=0; i<intcos.size(); ++i) {
      SIMPLE *q = intcos.at(i);

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

          if ( (((int) Z[a]) == 1) || (((int) Z[c]) == 1) )
            val = 0.160;
          else
            val = 0.250;
          f[cnt++] = val;
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

      } // end switch coordinate type
    } // loop over intcos
  } // end Schlegel
  else if (Opt_params.intrafragment_H == OPT_PARAMS::FISCHER) {
    for (i=0; i<intcos.size(); ++i) {
      SIMPLE *q = intcos.at(i);

      switch (q->g_type()) { 

        case (stre_type) :
          if (q->is_hbond()) f[cnt++] = 0.03;
          else {
            a = q->g_atom(0);
            b = q->g_atom(1);
  
            rABcov = Rcov(Z[a], Z[b]);
  
            A = 0.3601; B = 1.944;
            f[cnt++] = A * exp(-B*(R[a][b] - rABcov)); 
          }
        break;

        case (bend_type) :
          a = q->g_atom(0);
          b = q->g_atom(1);
          c = q->g_atom(2);

          rABcov = Rcov(Z[a], Z[b]);
          rBCcov = Rcov(Z[b], Z[c]);

          A = 0.089; B = 0.11; C = 0.44; D = -0.42;
          f[cnt++] = A + B/(pow(rABcov*rBCcov, D)) *
            exp(-C*( R[a][b] + R[b][c] - rABcov - rBCcov));
        break;

        case (tors_type) :
          b = q->g_atom(1);
          c = q->g_atom(2);

          rBCcov = Rcov(Z[b], Z[c]);

          // count number of additional bonds on central atoms
          L = 0;
          for (j=0; j<natom; ++j) {
            if (j == c) continue;
            if (connectivity[b][j]) ++L;
          }
          for (j=0; j<natom; ++j) {
            if (j == b) continue;
            if (connectivity[c][j]) ++L;
          }
          A = 0.0015; B = 14.0; C = 2.85; D = 0.57; E = 4.00;
          f[cnt++] = A + B * pow(L,D) / pow(R[b][c] * rBCcov, E) * exp(-C * (R[b][c] - rBCcov));
        break;
      } // end switch intcos
    } // end loop intcos
  } // end Fischer
  else if (Opt_params.intrafragment_H == OPT_PARAMS::SIMPLE) {
    for (i=0; i<intcos.size(); ++i) {
      SIMPLE *q = intcos.at(i);
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
      }
    }
  }
  else if (Opt_params.intrafragment_H == OPT_PARAMS::LINDH) {

    const double k_r   = 0.45;
    const double k_phi = 0.15;
    const double k_tau = 0.005;
    double k;

    for (i=0; i<intcos.size(); ++i) {
      SIMPLE * q = intcos.at(i);

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
          f[cnt++] = k_phi * Lindh_rho(a, b, R[a][b])
                    * Lindh_rho(b, c, R[b][c]);
        break;

        case (tors_type) :
          a = q->g_atom(0);
          b = q->g_atom(1);
          c = q->g_atom(2);
          d = q->g_atom(3);
          f[cnt++] = k_tau * Lindh_rho(a, b, R[a][b])
                    * Lindh_rho(b, c, R[b][c])
                    * Lindh_rho(c, d, R[c][d]);
        break;
      }
    }
  }

  free_matrix(R);

  double **H = init_matrix(intcos.size(), intcos.size());
  for (i=0; i<intcos.size(); ++i)
    H[i][i] = f[i];
  free_array(f);
  return H;
}

}

