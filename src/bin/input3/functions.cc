/*! \file
 \ingroup INPUT
 \brief Enter brief description of file here 
 */
/* This file contains functions such as Rotation, distance, and bond angle
 calculations */

#define EXTERN
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libchkpt/chkpt.h>
#include <cmath>
#include <symmetry.h>
#include "input.h"
#include "global.h"
#include "defines.h"
#include <physconst.h>

namespace psi { namespace input {

using namespace psi;

void punt(const char *mess) {
  fprintf(outfile, "  error: %s\n", mess);
  fprintf(stderr, "  INPUT error: %s\n", mess);
  stop_io();
  abort();
}

/* 
 UNIT_VECTORS(): Function calculates unit vectors between each
 pair of atoms and stores them in E. 
 */
double ***unit_vectors(double **Geom, double **Distances) {
  int i, j;
  double ***E;
  double distance;
  
  E = (double ***)malloc(num_atoms*sizeof(double **));
  for (i=0; i<num_atoms; i++) {
    E[i] = (double **)malloc(num_atoms*sizeof(double *));
    for (j=0; j<num_atoms; j++) {
      E[i][j] = (double *)malloc(3*sizeof(double));
    }
  }
  
  for (i=0; i<num_atoms; i++) {
    for (j=0; j<num_atoms; j++) {
      if (i != j) {
        distance = sqrt((Geom[i][0] - Geom[j][0]) * (Geom[i][0] - Geom[j][0])
            + (Geom[i][1] - Geom[j][1]) * (Geom[i][1] - Geom[j][1])
            + (Geom[i][2] - Geom[j][2]) * (Geom[i][2] - Geom[j][2]));
        E[i][j][0] = -(Geom[i][0] - Geom[j][0]) / distance;
        E[i][j][1] = -(Geom[i][1] - Geom[j][1]) / distance;
        E[i][j][2] = -(Geom[i][2] - Geom[j][2]) / distance;
      } else {
        E[i][j][0] = 0.0;
        E[i][j][1] = 0.0;
        E[i][j][2] = 0.0;
      }
    }
  }
  
  return (E);
}

/*
 ** CALC_BOND_ANGLES(): Function calculates all non-redundant bond angles
 **   between all combinations of three atoms, and writes to file fpo.
 **   Note: the first index of BondAngles is always smaller; i.e.
 **   angle 5-3-1 is stored only as 1-3-5
 */
double ***calc_bond_angles(double ***E, double **Distances) {
  int i, j, k;
  double dotprod=0;
  double angle=0;
  
  double ***BA;
  
  BA = (double ***)malloc(num_atoms*sizeof(double **));
  for (i=0; i<num_atoms; i++) {
    BA[i] = (double **)malloc(num_atoms*sizeof(double *));
    for (j=0; j<num_atoms; j++) {
      BA[i][j] = (double *)malloc(num_atoms*sizeof(double));
    }
  }
  
  fprintf(outfile, "\n  -Bond Angles:\n\n") ;
  for (i=0; i<num_atoms; i++) {
    for (j=0; j<num_atoms; j++) {
      for (k=i; k<num_atoms; k++) {
        if ((i != j) && (i != k)&& (j != k)) {
          dotprod = dot_prod(E[j][i], E[j][k]);
          if (dotprod > 1.00000)
            angle = 0.0000;
          else if (dotprod < -1.00000)
            angle = _pi;
          else
            angle = acos(dotprod);
          BA[i][j][k] = angle;
          if (((Distances[i][j] < 3.000) && (Distances[j][k] < 3.000))) {
            fprintf(outfile, "   %2d -%2d -%2d    %14.8lf\n", i+1, j+1, k+1,
                    (angle*180.00/_pi));
          }
        }
      }
    }
  }
  
  return (BA);
}

void calc_oop_angles(double ***E, double ***BondAngles, double **R) {
  int i, j, k, l;
  double tval1=0;
  double tval2=0;
  double angle=0;
  double crossprod[3];
  
  fprintf(outfile, "\n  -Out of Plane Angles:\n") ;
  fprintf(outfile, "   (Angle between vector 14 and plane 234)\n\n");
  
  for (i=0; i<num_atoms; i++) {
    for (j=0; j<num_atoms; j++) {
      for (k=j; k<num_atoms; k++) {
        for (l=0; l<num_atoms; l++) {
          if ( (i!=j) && (i!=k) & (i!=l)&& (j!=k)&& (j!=l)&& (k!=l)) {
            cross_prod(E[l][j], E[l][k], crossprod) ;
            tval1 = dot_prod(crossprod, E[l][i]) ;
            /* make sure j<k when BondAngles is accessed */
            if (j < k)
              tval2 = BondAngles[j][l][k];
            else
              tval2 = BondAngles[k][l][j];
            if (sin(tval2) > 0.0001)
              tval1 = tval1 / sin(tval2) ;
            else
              tval1 = 0.0;
            if (tval1 > 1.000000)
              angle = (_pi / 2.00000);
            else if (tval1 < -1.000000)
              angle = -(_pi / 2.0000);
            else
              angle = asin(tval1) ;
            if (((R[i][j] < 3.0) && (R[j][k] < 3.0)&& (R[k][l] < 3.0))) {
              fprintf(outfile, "   %2d -%2d -%2d -%2d    %13.8lf\n", i+1, j+1,
                      k+1, l+1, angle * 180.0/(double)_pi) ;
            }
          }
        }
      }
    }
  }
}

void calc_tors_angles(double ***E, double ***BondAngles, double **R) {
  int i, j, k, l, xyz;
  double cross1[3], cross2[3], cross3[3];
  double tval, angle, phi2, phi3;
  double sign, norm3;
  
  fprintf(outfile, "\n  -Torsional Angles:\n\n");
  for (i=0; i<num_atoms; i++) {
    for (j=0; j<num_atoms; j++) {
      for (k=0; k<num_atoms; k++) {
        for (l=i; l<num_atoms; l++) {
          if ( (i!=j) && (i!=k) & (i!=l)&& (j!=k)&& (j!=l)&& (k!=l)) {
            cross_prod(E[i][j], E[j][k], cross1) ;
            cross_prod(E[j][k], E[k][l], cross2) ;
            if (i < k)
              phi2 = BondAngles[i][j][k];
            else
              phi2 = BondAngles[k][j][i];
            if (j < l)
              phi3 = BondAngles[j][k][l];
            else
              phi3 = BondAngles[l][k][j];
            tval = dot_prod(cross1, cross2) ;
            if ((sin(phi2) > 0.00001) && (sin(phi3) > 0.00001)) {
              tval /= sin(phi2) ;
              tval /= sin(phi3) ;
            } else
              tval = 2.0;
            if (tval > 0.99999)
              angle = 0.0000;
            else if (tval < -0.99999)
              angle = _pi;
            else
              angle = acos(tval) ;
            
            /* compute the sign */
            cross_prod(cross1, cross2, cross3);
            norm3 = sqrt(dot_prod(cross3, cross3));
            sign = 1.0;
            if (fabs(norm3) > 0.00001) {
              for (xyz=0; xyz<3; ++xyz)
                cross3[xyz] *= 1.0/norm3;
              tval = dot_prod(cross3, E[j][k]);
              if (tval < 0.0)
                sign = -1.0;
            }
            
            if ( ((R[i][j] < 3.0) && (R[j][k] < 3.0)&&(R[k][l] < 3.0))) {
              fprintf(outfile, "   %2d -%2d -%2d -%2d    %13.8lf\n", i+1, j+1,
                      k+1, l+1, sign * angle * 180.0/_pi) ;
            }
          }
        }
      }
    }
  }
}

/*
 Function to calculate all internuclear distances 
 */
void calc_distance(double **geom, double *A, int num) {
  
  int i = 0;
  int j = 0;
  double dist;
  double temp1=0.0;
  double temp2=0.0;
  double temp3=0.0;
  
  for (i=1; i<num; i++) {
    for (j=0; j<i; j++) {
      temp1 = (geom[i][0]-geom[j][0])*(geom[i][0]-geom[j][0]);
      temp2 = (geom[i][1]-geom[j][1])*(geom[i][1]-geom[j][1]);
      temp3 = (geom[i][2]-geom[j][2])*(geom[i][2]-geom[j][2]);
      dist = sqrt(temp1 + temp2 + temp3);
      if (!expert && dist < ZERO_BOND_DISTANCE&&nuclear_charges[i] != 0.0
          &&nuclear_charges[j] != 0.0) {
        printf("  Atoms %d and %d are separated by only %lf!\n", i+1, j+1, dist);
        punt("Invalid geometry");
      } else
        A[ioff[i]+j] = dist;
    }
  }
}

/*
 Function to calculate the nuclear repulsion energy
 */
void Nuc_repulsion(double *Distance, double *repulsion) {
  
  int i, j;
  *repulsion = 0.0;
  
  for (i=1; i<num_atoms; i++)
    for (j=0; j<i; j++)
      if (nuclear_charges[i] != 0.0&& nuclear_charges[j] != 0.0)
        *repulsion += nuclear_charges[i]*nuclear_charges[j]/Distance[ioff[i]+j];
  return;
}

/*
 Function returns 1 if vectors are equivalent and 0 if they are not
 */

int vect_equiv(double *A, double *B) {
  if ((fabs(A[0] - B[0]) < ZERO) && (fabs(A[1] - B[1]) < ZERO)&& (fabs(A[2]
      - B[2]) < ZERO))
    return 1;
  else
    return 0;
}

/*
 Function returns unit vector - median of vectors A and B
 */

void median_vec(double *A, double *B, double *median) {
  double norm, x, y, z;
  
  x = (A[0] + B[0])/2;
  y = (A[1] + B[1])/2;
  z = (A[2] + B[2])/2;
  
  norm = sqrt(x*x + y*y + z*z);
  median[0] = x/norm;
  median[1] = y/norm;
  median[2] = z/norm;
  return;
}

/*
 Function returns 1 if vectors are inversion symmetry related, 0 otherwise
 */

int inv_related(double *A, double *B) {
  if ((fabs(A[0] + B[0]) < ZERO) && (fabs(A[1] + B[1]) < ZERO)&& (fabs(A[2]
      + B[2]) < ZERO))
    return 1;
  else
    return 0;
}

/*
 Function rotates geometry[][] by new_coord[][] and memorizes it
 */

void rotate_full_geom(double **new_coord) {
  double **new_geom;
  int i;
  
  new_geom = block_matrix(num_allatoms, 3);
  mmult(full_geom, 0, new_coord, 0, new_geom, 0, num_allatoms, 3, 3, 0);
  for (i=0; i<num_allatoms; i++) {
    full_geom[i][0] = new_geom[i][0];
    full_geom[i][1] = new_geom[i][1];
    full_geom[i][2] = new_geom[i][2];
  }
  
  free_block(new_geom);
  
  memorize_rotation(new_coord);
  
  return;
}

/*
 Function updates Rref to remember the effect of the rotation
 described by R
 */

void memorize_rotation(double **R) {
  double **Rref_new;
  int i, j;
  
  Rref_new = block_matrix(3, 3);
  mmult(R, 1, Rref, 0, Rref_new, 0, 3, 3, 3, 0);
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      Rref[i][j] = Rref_new[i][j];
  free_block(Rref_new);
  
  return;
}

/*
 Function forms a (3x3) matrix out of 3 3-vectors (aligned in columns)
 */

void vectors_to_matrix(double *v1, double *v2, double *v3, double **matrix) {
  matrix[0][0] = v1[0];
  matrix[1][0] = v1[1];
  matrix[2][0] = v1[2];
  matrix[0][1] = v2[0];
  matrix[1][1] = v2[1];
  matrix[2][1] = v2[2];
  matrix[0][2] = v3[0];
  matrix[1][2] = v3[1];
  matrix[2][2] = v3[2];
  return;
}

/*Functions to allocate memory for character arrays*/

char *init_char_array(int B) {
  char *A;
  int i;
  A = (char *)malloc((B+1)*sizeof(char));
  for (i=0; i<B; i++) {
    A[i] = ' ';
  }
  A[B] = '\0';
  return A;
}

char **init_char_matrix(int A, int B) {
  char **mat;
  int i, j;
  
  mat = (char **)malloc(A*sizeof(char *));
  
  for (i=0; i<A; i++)
    mat[i] = init_char_array(B);
  return mat;
}

void free_char_matrix(char **mat, int A) {
  int i;
  
  for (i=0; i<A; i++)
    free(mat[i]);
  free(mat);
}

void setup(int maxioff) {
  int i, maxdf;
  
  if (maxioff<100)
    maxioff = 100;
  maxdf = 300;
  
  df = (double *) malloc(sizeof(double)*maxdf);
  ioff = (int *) malloc(sizeof(int)*maxioff);
  
  /* df[i] = (i-1)!! */
  df[0] = 1;
  df[1] = 1;
  df[2] = 1;
  for (i=3; i<maxdf; i++)
    df[i] = (i-1)*df[i-2];
  
  ioff[0] = 0;
  for (i=1; i<maxioff; i++)
    ioff[i] = ioff[i-1]+i;
}

double int_pow(double a, int p)
{
  register int i;
  double b = 1.0;
  
  for (i=0; i<p; i++)
    b = b*a;
  return b;
}

/*
 Make canonical and reference frames equivalent
 */
void canon_eq_ref_frame() {
  int i, j;
  
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      Rref[i][j] = (i == j) ? 1.0 : 0.0;
  
  return;
}

}} // namespace psi::input
