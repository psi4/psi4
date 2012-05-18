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

// covalent bond length in bohr from atomic numbers
inline double Rcov(double ZA, double ZB) {
  return (cov_radii[(int) ZA] + cov_radii[(int) ZB]) / _bohr2angstroms;
}

// period from atomic number
int period(double ZA) { 
  int iZA;
  iZA = (int) ZA;
  if      (iZA <=  2) return 1;
  else if (iZA <= 10) return 2;
  else if (iZA <= 18) return 3;
  else if (iZA <= 36) return 4;
  else                return 5;
}

double ** FRAG::H_guess(void) {
  int i, j, atomA, atomB, atomC, atomD, cnt = 0, perA, perB;
  double rAB, rBC, rABcov, rBCcov, rBDcov;
  double A, B, C, D, L, E, val;
  double *f;

  f = init_array(intcos.size());

  // Form diagonal Hessian in simple internals
  if (Opt_params.intrafragment_H == OPT_PARAMS::SCHLEGEL) {
    for (i=0; i<intcos.size(); ++i) {

      switch (intcos[i]->g_type()) { 

        case (stre_type) :
          if (intcos[i]->is_hbond()) f[cnt++] = 0.03;
          else {
            atomA = intcos[i]->g_atom(0);
            atomB = intcos[i]->g_atom(1);
            perA = period(Z[atomA]);
            perB = period(Z[atomB]);
  
            if ( perA==1 && perB==1) B = -0.244;
            else if ((perA==1 && perB==2) || (perB==1 && perA==2)) B = 0.352;
            else if (perA==2 && perB==2)                           B = 1.085;
            else if ((perA==1 && perB==3) || (perB==1 && perA==3)) B = 0.660;
            else if ((perA==2 && perB==3) || (perB==2 && perA==3)) B = 1.522;
            else B = 2.068;
    
            rAB = v3d_dist(geom[atomA], geom[atomB]);
  
            A = 1.734;
            // force constants in au / bohr^2
            f[cnt++] = A/((rAB-B)*(rAB-B)*(rAB-B));
            // for aJ/Ang^2, * _hartree2aJ / SQR(_bohr2angstroms);
          }
        break;

        case (bend_type) :
          atomA = intcos[i]->g_atom(0);
          atomB = intcos[i]->g_atom(1);
          atomC = intcos[i]->g_atom(2);

          if ( (((int) Z[atomA]) == 1) || (((int) Z[atomC]) == 1) )
            val = 0.160;
          else
            val = 0.250;
          f[cnt++] = val;
        break;

        case (tors_type) :
          atomB = intcos[i]->g_atom(1);
          atomC = intcos[i]->g_atom(2);
          A = 0.0023;
          B = 0.07;
          rBCcov = Rcov(Z[atomB], Z[(atomC)]);
          rBC = v3d_dist(geom[atomB], geom[atomC]);
          if (rBC > (rBCcov + A/B)) B = 0.0; // keep > 0
          f[cnt++] = (A - (B*(rBC - rBCcov)));
        break;

      } // end switch coordinate type
    } // loop over intcos
  } // end Schlegel
  else if (Opt_params.intrafragment_H == OPT_PARAMS::FISCHER) {
    for (i=0; i<intcos.size(); ++i) {

      switch (intcos[i]->g_type()) { 

        case (stre_type) :
          if (intcos[i]->is_hbond()) f[cnt++] = 0.03;
          else {
            atomA = intcos[i]->g_atom(0);
            atomB = intcos[i]->g_atom(1);
  
            rAB   = v3d_dist(geom[atomA], geom[atomB]);
            rABcov = Rcov(Z[atomA], Z[atomB]);
  
            A = 0.3601; B = 1.944;
            f[cnt++] = A * exp(-B*(rAB - rABcov)); 
          }
        break;

        case (bend_type) :
          atomA = intcos[i]->g_atom(0);
          atomB = intcos[i]->g_atom(1);
          atomC = intcos[i]->g_atom(2);

          rABcov = Rcov(Z[atomA], Z[atomB]);
          rBCcov = Rcov(Z[atomB], Z[atomC]);
          rAB = v3d_dist(geom[atomA], geom[atomB]);
          rBC = v3d_dist(geom[atomB], geom[atomC]);

          A = 0.089; B = 0.11; C = 0.44; D = -0.42;
          f[cnt++] = A + B/(pow(rABcov*rBCcov, D)) *
            exp(-C*( rAB + rBC - rABcov - rBCcov));
        break;

        case (tors_type) :
          atomB = intcos[i]->g_atom(1);
          atomC = intcos[i]->g_atom(2);

          rBCcov = Rcov(Z[atomB], Z[atomC]);
          rBC = v3d_dist(geom[atomB], geom[atomC]);

          // count number of additional bonds on central atoms
          L = 0;
          for (j=0; j<natom; ++j) {
            if (j == atomC) continue;
            if (connectivity[atomB][j]) ++L;
          }
          for (j=0; j<natom; ++j) {
            if (j == atomB) continue;
            if (connectivity[atomC][j]) ++L;
          }
          A = 0.0015; B = 14.0; C = 2.85; D = 0.57; E = 4.00;
          f[cnt++] = A + B * pow(L,D) / pow(rBC * rBCcov, E) * exp(-C * (rBC - rBCcov));
        break;
      } // end switch intcos
    } // end loop intcos
  } // end Fischer
  else if (Opt_params.intrafragment_H == OPT_PARAMS::SIMPLE) {
    for (i=0; i<intcos.size(); ++i) {
      switch (intcos[i]->g_type()) {
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


  double **H = init_matrix(intcos.size(), intcos.size());
  for (i=0; i<intcos.size(); ++i)
    H[i][i] = f[i];
  free_array(f);
  return H;
}

}

