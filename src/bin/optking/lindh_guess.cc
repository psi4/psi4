/*! \file Lindh_guess.cc
    \ingroup optking
    \brief Function to generate a model cartesian Hessian according to
     R. Lindh, A. Bernhardsson, G. Karlstrom, P.-A. Malmqvist, CPL, 241, 423, 1995.
*/

#include "molecule.h"
#include <cmath>
#include "v3d.h"

#define EXTERN
#include "globals.h"

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

inline int period(int Z) {
  if      (Z <=  2) return 1;
  else if (Z <= 10) return 2;
  else if (Z <= 18) return 3;
  else if (Z <= 36) return 4;
  else              return 5;
}

inline double alpha_table(int perA, int perB) {
  if (perA == 1) {
    if (perB == 1) return 1000;
    else return 0.3949;
  }
  else {
    if (perB == 1)
      return 0.3949;
    else
      return 0.2800;
  }
}

inline double r_ref_table(int perA, int perB) {
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
double FRAG::Lindh_rho(int A, int B) const {

  int perA = period((int) Z[A]);
  int perB = period((int) Z[B]);

  double alpha = alpha_table(perA, perB);
  double r_ref = r_ref_table(perA, perB);

  double r = v3d_dist(geom[A], geom[B]);

  return (exp(alpha * (r_ref*r_ref - r*r)));
}

double ** FRAG::Lindh_guess(void) {
  int N_xyz = 3*natom;

  // Generate coordinates between ALL atoms, but not 
  // complete reversals or using duplicate i,j,k,l.
  // Create stretch coordinates.
  for (int i=0; i<natom; ++i) {
    for (int j=i+1; j<natom; ++j) {
      STRE *one_stre = new STRE(i, j);
      intcos.push_back(one_stre);
printf("stre %d %d \n", i, j);
    }
  }

  // Create bend coodinates.
  for (int i=0; i<natom; ++i) {
    for (int j=0; j<natom; ++j) {
      if (j != i) {
        for (int k=i+1; k<natom; ++k) {
          if (k != j) {
            BEND *one_bend = new BEND(i, j, k);
            intcos.push_back(one_bend);
printf("bend %d %d %d\n", i, j, k);
          }
        }
      }
    }
  }

  // Create torsion coodinates.
  for (int i=0; i<natom; ++i) {
    for (int j=0; j<natom; ++j) {
      if (j != i) {
        for (int k=0; k<natom; ++k) {
          if (k!=i && k!=j) {
            for (int l=i+1; l<natom; ++l) {
              if (l!=j && l!=k) {
                TORS *one_tors = new TORS(i, j, k, l);
                intcos.push_back(one_tors);
printf("tors %d %d %d %d\n", i, j, k, l);
              }
            }
          }
        }
      }
    }
  }

// We are ignoring gradient term for now.  It's not clear if Helgaker used it.
  const double k_r   = 0.45;
  const double k_phi = 0.15;
  const double k_tau = 0.005;

  double **H_xyz = init_matrix(N_xyz, N_xyz);
  double k;

  for (int i=0; i<intcos.size(); ++i) {  // loop over intcos
    SIMPLE * q = intcos.at(i);

    double **Bintco = q->DqDx(geom); // dq_i / da_xyz
    int natom_intco = q->g_natom();

    if (q->g_type() == stre_type) {
      k = k_r * Lindh_rho(q->g_atom(0), q->g_atom(1));
    }
    else if (intcos.at(i)->g_type() == bend_type) {
      k = k_phi * Lindh_rho(q->g_atom(0), q->g_atom(1))
                * Lindh_rho(q->g_atom(1), q->g_atom(2));
    }
    else if (intcos.at(i)->g_type() == tors_type) {
      k = k_tau * Lindh_rho(q->g_atom(0), q->g_atom(1))
                * Lindh_rho(q->g_atom(1), q->g_atom(2))
                * Lindh_rho(q->g_atom(2), q->g_atom(3));
    }

    for (int a=0; a < natom_intco; ++a) {
      int x1 = q->g_atom(a);

      for (int b=0; b < natom_intco; ++b) {
        int x2 = q->g_atom(b);

        for (int xyz1=0; xyz1<3; ++xyz1)
          for (int xyz2=0; xyz2<3; ++xyz2)
             H_xyz[3*x1 + xyz1][3*x2 + xyz2] += k * Bintco[a][xyz1] * Bintco[b][xyz2];
      }
    }
    free_matrix(Bintco);
  }

  if (Opt_params.print_lvl >= 2) {
    fprintf(outfile,"Lindh cartesian Hessian guess\n");
    print_matrix(outfile, H_xyz, natom, natom);
    fflush(outfile);
  }
  return H_xyz;
}

}
