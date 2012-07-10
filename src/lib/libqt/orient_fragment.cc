/*! \file
    \ingroup QT
    \brief function which places one fragment into the coordinate system of another

   int natoms_A, natoms_B - number of atoms in each fragment
   int P_A, P_B - number of ref pts to worry about, less than 3 for linear fragments
   double **geom_A, double **geom_B - geometry of fragments A and B - geom B is changed
   double **ref_coeff_A, double **ref_coeff_B - linear combinations which specify reference atoms
   double R_AB - distance between reference atoms #1 on each fragment
   theta_A, theta_B, tau, chi-A, chi-B - interfragment angles
*/
#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <math.h>
#include <physconst.h>
#include <psifiles.h>
#include <psi4-dec.h>

namespace psi {

void orient_fragment(int natom_A, int natom_B, int P_A, int P_B, double **geom_A, double **geom_B,
  double **ref_coeff_A, double **ref_coeff_B, double R_AB, double theta_A, double theta_B,
  double tau, double phi_A, double phi_B, FILE *outfile)
{
  int i, j, errcod, pts, xyz;
  double tval, norm, B_angle, R_B1B2, R_B2B3, e12[3], e12b[3], e12c[3], e12d[3], erot[3];
  double **ref_A, **ref_B, **ref_B_final;
  double sign, cross1[3], cross2[3], cross3[3], phi2, phi3;

  ref_A = block_matrix(3,3);
  ref_B = block_matrix(3,3);
  ref_B_final = block_matrix(3,3);

  /* stick SOMETHING in for non-specified reference atoms - necessary to prevent complaints
     about collinear reference atoms and to make zmat_point() work in such cases */
  if (P_A < 3) {
      for (xyz=0; xyz<3; ++xyz)
        ref_A[2][xyz] = (xyz+1)/_pi;
  }
  if (P_A < 2) {
      for (xyz=0; xyz<3; ++xyz)
        ref_A[1][xyz] = (xyz+1)/(2*_pi);
  }

  for (pts=0; pts<P_A; ++pts)
    for (xyz=0; xyz<3; ++xyz)
      for (i=0; i<natom_A; ++i)
        ref_A[pts][xyz] += ref_coeff_A[pts][i] * geom_A[i][xyz];

  for (pts=0; pts<P_B; ++pts)
    for (xyz=0; xyz<3; ++xyz)
      for (i=0; i<natom_B;++i)
        ref_B[pts][xyz] += ref_coeff_B[pts][i] * geom_B[i][xyz];

fprintf(outfile,"Coordinates for reference points on fragment A\n");
print_mat(ref_A,P_A,3,outfile);
fprintf(outfile,"Coordinates for reference points on fragment B (original) \n");
print_mat(ref_B,P_B,3,outfile);
fprintf(outfile,"\t(1/)R_AB:%10.5f, theta_A:%10.5f, theta_B:%10.5f\n", R_AB, theta_A, theta_B);
fprintf(outfile,"\t     tau:%10.5f,   phi_A:%10.5f,   phi_B:%10.5f\n", tau, phi_A, phi_B);

  /* compute B1-B2 distance, B2-B3 distance, and B1-B2-B3 angle */
  R_B1B2 = 0.0;
  if (P_B>1) {
    for (xyz=0; xyz<3; ++xyz)
      R_B1B2 += (ref_B[1][xyz]-ref_B[0][xyz])*(ref_B[1][xyz]-ref_B[0][xyz]);
    R_B1B2 = sqrt(R_B1B2);
  }
  R_B2B3 = 0.0;
  B_angle = 0.0;
  if (P_B>2) {
    for (xyz=0; xyz<3; ++xyz)
      R_B2B3 += (ref_B[2][xyz]-ref_B[1][xyz])*(ref_B[2][xyz]-ref_B[1][xyz]);
    R_B2B3 = sqrt(R_B2B3);
    unit_vec(ref_B[1],ref_B[0],e12);
    unit_vec(ref_B[1],ref_B[2],e12b);
    B_angle = acos(dot_prod(e12,e12b))*180.0/_pi;
  }

    /* determine location of reference pts for B in coordinate system of A */
  zmat_point(ref_A[2], ref_A[1], ref_A[0], R_AB, theta_A, phi_A, ref_B_final[0]);
  if (P_B>1)
    zmat_point(ref_A[1], ref_A[0], ref_B_final[0], R_B1B2, theta_B, tau, ref_B_final[1]);
  if (P_B>2)
    zmat_point(ref_A[0], ref_B_final[0], ref_B_final[1], R_B2B3, B_angle, phi_B, ref_B_final[2]);

fprintf(outfile,"Target reference points for fragment B\n");
print_mat(ref_B_final,P_B,3,outfile);

  /* translate geom_B to place B1 in correct location */
  for (xyz=0; xyz<3; ++xyz) {
    tval = ref_B_final[0][xyz] - ref_B[0][xyz];
    for (i=0; i<natom_B; ++i)
      geom_B[i][xyz] += tval;
  }

  for (pts=0; pts<P_B; ++pts)
    for (xyz=0; xyz<3; ++xyz) {
      ref_B[pts][xyz] = 0.0;
      for (i=0; i<natom_B;++i)
        ref_B[pts][xyz] += ref_coeff_B[pts][i] * geom_B[i][xyz];
    }

//fprintf(outfile,"Reference points after translation (to fix point B1):\n");
//print_mat(ref_B,P_B,3,outfile);

  if (P_B>1) { /* move fragment B to place reference point B2 in correct location */
    /* Determine rotational angle and axis */
    unit_vec(ref_B[1],       ref_B[0], e12);  /* v B1->B2 */
    unit_vec(ref_B_final[1], ref_B[0], e12b); /* v B1->B2_final */
    B_angle = acos(dot_prod(e12b,e12));
    fprintf(outfile,"Rotation by %f degrees (to fix point B2)\n", 180.0*B_angle/_pi);
    if (fabs(B_angle) > 1.0e-7) {
      cross_prod(e12,e12b,erot);

      /* Move B to put B1 at origin */
      for (xyz=0; xyz<3; ++xyz)
        for (i=0; i<natom_B;++i)
          geom_B[i][xyz] -= ref_B[0][xyz];

      /* Rotate B */
      rotate_vecs(erot, B_angle, geom_B, natom_B);

      /* Move B back to coordinate system of A */
      for (xyz=0; xyz<3; ++xyz)
        for (i=0; i<natom_B;++i)
          geom_B[i][xyz] += ref_B[0][xyz];

      /* Check location of reference points now */
      for (pts=0; pts<P_B; ++pts)
        for (xyz=0; xyz<3; ++xyz) {
          ref_B[pts][xyz] = 0.0;
          for (i=0; i<natom_B;++i)
            ref_B[pts][xyz] += ref_coeff_B[pts][i] * geom_B[i][xyz];
        }

      //fprintf(outfile,"Reference points after rotation (to fix point B2) \n");
      //print_mat(ref_B,P_B,3,outfile);
    }
  }

  if (P_B==3) { /* move fragment B to place reference point B3 in correct location */
    /* Determine rotational angle and axis */
    unit_vec(ref_B[1], ref_B[0], erot);  /* B1 -> B2 is rotation axis */

    /* Calculate B3-B1-B2-B3' torsion angle */
    unit_vec(ref_B[2], ref_B[0], e12);  /* v B1->B3 */
    unit_vec(ref_B[1], ref_B[0], e12b); /* v B1->B2 */
    phi2 = acos(dot_prod(e12,e12b));
    unit_vec(ref_B[0], ref_B[1], e12c);  /* v B2->B1 */
    unit_vec(ref_B_final[2], ref_B[1], e12d); /* v B2->B3' */
    phi3 = acos(dot_prod(e12c,e12d));

    cross_prod(e12 , e12b, cross1) ; /* B3->B1 x B1->B2 */
    cross_prod(e12c, e12d, cross2) ; /* B1->B2 x B2->B3 */
    tval = dot_prod(cross1, cross2) ;

    if ((sin(phi2) > 0.00001) && (sin(phi3) > 0.00001)) {
      tval /= sin(phi2) ;
      tval /= sin(phi3) ;
    }
    else tval = 2.0;

    if (tval > 0.99999) B_angle = 0.0000;
    else if (tval < -0.99999) B_angle = _pi;
    else B_angle = acos(tval) ;

    sign = 1.0; /* check sign */
    cross_prod(cross1, cross2, cross3);
    norm = sqrt(dot_prod(cross3, cross3));
    if (fabs(norm) > 0.00001) {
      for (xyz=0; xyz<3; ++xyz)
        cross3[xyz] *= 1.0/norm;
      tval = dot_prod(cross3, e12b);
      if (tval < 0.0) sign = -1.0;
    }
    B_angle *= sign;

    if (fabs(B_angle) > 1.0e-7) {
      fprintf(outfile,"Rotation by %f degrees (to fix point B3)\n", 180.0*B_angle/_pi);

      /* Move B to put B2 at origin */
      for (xyz=0; xyz<3; ++xyz)
        for (i=0; i<natom_B;++i)
          geom_B[i][xyz] -= ref_B[1][xyz];

      rotate_vecs(erot, B_angle, geom_B, natom_B);

      /* Translate B1 back to coordinate system of A */
      for (xyz=0; xyz<3; ++xyz)
        for (i=0; i<natom_B;++i)
          geom_B[i][xyz] += ref_B[1][xyz];

      for (pts=0; pts<P_B; ++pts)
        for (xyz=0; xyz<3; ++xyz) {
          ref_B[pts][xyz] = 0.0;
          for (i=0; i<natom_B;++i)
            ref_B[pts][xyz] += ref_coeff_B[pts][i] * geom_B[i][xyz];
        }
        //fprintf(outfile,"Reference points on B after rotation for B2 \n");
        //print_mat(ref_B,P_B,3,outfile);
      }
  }

   /* check to see if desired reference points were obtained */
   tval = 0.0;
   for (i=0; i<P_B; ++i)
     for (xyz=0; xyz<3; ++xyz)
       tval += fabs(ref_B[i][xyz] - ref_B_final[i][xyz]);
   if (tval > 1.0e10) {
     throw PsiException("Unable to construct multi-fragment geometry.",__FILE__,__LINE__);
   }
   else
     fprintf(outfile,"Successfully constructed multifragment geometry.\n");

  free_block(ref_A);
  free_block(ref_B);
  free_block(ref_B_final);
  return;
}

}
