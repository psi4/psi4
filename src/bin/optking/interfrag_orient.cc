/*! \file interfrag_orient.cc
    \ingroup optking
    \brief function moves the geometry of fragment B so that the interfragment coordinates
      have the given values

   ndA = # of ref pts on A to worry about
   ndB = # of ref pts on B to worry about

   Value at least
    ndA   ndB
   ------------
     1     1     R_AB
     2     1     + theta_A
     1     2     + theta_B
     2     2     + theta_A + theta_B + fix tau
     3     2     + phi_A
     2     3     + phi_A + phi_B
   ------------
   
   returns true if successful, false if not
*/

#include "frag.h"
#include "interfrag.h"
#include "print.h"
#include "v3d.h"

#define EXTERN
#include "globals.h"

namespace opt {

using namespace v3d;

void zmat_point(double *A, double *B, double *C, double R_CD, double theta_BCD,
  double phi_ABCD, double *D);

void rotate_vecs(double *w, double phi, double **v, int num_v);

// arguments specify the bond length and angle in radians desired for interfragment coordinates
bool INTERFRAG::orient_fragment(double *dq, double *fq) {

  int pts, i, xyz;
  double tval, norm, B_angle, R_B1B2, R_B2B3, e12[3], e12b[3], e12c[3], e12d[3], erot[3];
  double **ref_A, **ref_B, **ref_B_final;
  double sign, cross1[3], cross2[3], cross3[3], phi2, phi3;

  // fill-in unused values with defaults to make code below work
  double R_AB, theta_A, theta_B, tau, phi_A, phi_B;
  R_AB  = 1.0;
  theta_A = theta_B = tau = phi_A = phi_B = _pi/2;

  double *q_orig = intco_values();
  double *q_target = init_array(g_nintco());

  for (i=0; i<g_nintco(); ++i) {
    if (D_on[i])
      q_target[i] = q_orig[i] + dq[i];
  }

  int cnt = 0;
  if (D_on[0]) R_AB    = q_target[cnt++];
  if (D_on[1]) theta_A = q_target[cnt++];
  if (D_on[2]) theta_B = q_target[cnt++];
  if (D_on[3]) tau     = q_target[cnt++];
  if (D_on[4]) phi_A   = q_target[cnt++];
  if (D_on[5]) phi_B   = q_target[cnt++];

  // Make labels for printing
  std::vector<string> lbl(6);
  for (i=0; i<6; ++i)
    if (inter_frag->intcos[i]->is_frozen())
      lbl[i] = "*";

  if (inter_frag->intcos[0]->is_inverse_stre())
    lbl[0] += "1/R_AB";
  else
    lbl[0] += "R_AB";
  lbl[1] += "theta_A";
  lbl[2] += "theta_B";
  lbl[3] += "tau";
  lbl[4] += "phi_A";
  lbl[5] += "phi_B";

  fprintf(outfile,"\t---Interfragment coordinates between fragments %d and %d\n", A_index+1, B_index+1);
  fprintf(outfile,"\t---Internal Coordinate Step in ANG or DEG, aJ/ANG or AJ/DEG ---\n");
  fprintf(outfile,"\t ----------------------------------------------------------------------\n");
  fprintf(outfile,"\t Coordinate             Previous        Force       Change         New \n");
  fprintf(outfile,"\t ----------             --------       ------       ------       ------\n");

  cnt = 0;
  for (i=0; i<6; ++i) {
    double val, force, change, target;
    if (D_on[i]) {
      val    = q_orig[cnt];
      force  = fq[cnt];
      change = dq[cnt];
      target = q_target[cnt];

      if (i == 0) { // change units for bond length coordinate
        if (inter_frag->intcos[0]->is_inverse_stre()) { // 1/R(AB)
          val    /= _bohr2angstroms;
          force  /= _hartree2aJ/_bohr2angstroms;
          change /= _bohr2angstroms;
          target /= _bohr2angstroms;
        }
        else { // R(AB)
          val    *= _bohr2angstroms;
          force  *= _hartree2aJ/_bohr2angstroms;
          change *= _bohr2angstroms;
          target *= _bohr2angstroms;
        }
      } else { // change units for angle in degrees
        val    *= 180.0/_pi;
        force  *= _hartree2aJ*_pi/180.0;
        change *= 180.0/_pi;
        target *= 180.0/_pi;
      }
      fprintf(outfile,"\t%-20s%12.5f%13.5f%13.5f%13.5f\n", lbl[i].c_str(), val, force, change, target);
      ++cnt;
    }
  }
  fprintf(outfile,  "\t ----------------------------------------------------------------------\n");

  // copy B->geom in case this fails
  double **B_geom = B->g_geom();

  ref_A = init_matrix(3,3);
  ref_B = init_matrix(ndB,3);
  ref_B_final = init_matrix(ndB,3);

  // stick SOMETHING in for non-specified reference atoms for zmat_point() function
  if (ndA < 3)
    for (xyz=0; xyz<3; ++xyz)
      ref_A[2][xyz] = (xyz+1);

  if (ndA < 2)
    for (xyz=0; xyz<3; ++xyz)
      ref_A[1][xyz] = (xyz+2);

  // compute current location of reference points on A and B
  for (pts=0; pts<ndA; ++pts)
    for (i=0; i<A->natom;++i)
      for (xyz=0; xyz<3; ++xyz)
        ref_A[pts][xyz] += weightA[pts][i] * A->geom[i][xyz];

  for (pts=0; pts<ndB; ++pts)
    for (i=0; i<B->natom;++i)
      for (xyz=0; xyz<3; ++xyz)
        ref_B[pts][xyz] += weightB[pts][i] * B_geom[i][xyz];

  // compute B1-B2 distance, B2-B3 distance, and B1-B2-B3 angle
  if (ndB>1)
    R_B1B2 = v3d_dist(ref_B[1], ref_B[0]);

  if (ndB>2) {
    R_B2B3 = v3d_dist(ref_B[2], ref_B[1]);
    v3d_angle(ref_B[0], ref_B[1], ref_B[2], B_angle);
  }

  // determine target location of reference pts for B in coordinate system of A
  zmat_point(ref_A[2], ref_A[1], ref_A[0], R_AB, theta_A, phi_A, ref_B_final[0]);
  if (ndB>1)
    zmat_point(ref_A[1], ref_A[0], ref_B_final[0], R_B1B2, theta_B, tau, ref_B_final[1]);
  if (ndB>2)
    zmat_point(ref_A[0], ref_B_final[0], ref_B_final[1], R_B2B3, B_angle, phi_B, ref_B_final[2]);

  //fprintf(outfile,"ref_B original location\n");
  //print_matrix(outfile, ref_B, ndB, 3);
  //fprintf(outfile,"ref_B_final target\n");
  //print_matrix(outfile, ref_B_final, ndB, 3);

  // translate B->geom to place B1 in correct location
  for (xyz=0; xyz<3; ++xyz) {
    tval = ref_B_final[0][xyz] - ref_B[0][xyz];
    for (i=0; i<B->natom; ++i)
      B_geom[i][xyz] += tval;
  }

  // recompute B reference points
  zero_matrix(ref_B, ndB, 3);
  for (pts=0; pts<ndB; ++pts)
    for (i=0; i<B->natom;++i)
      for (xyz=0; xyz<3; ++xyz)
        ref_B[pts][xyz] += weightB[pts][i] * B_geom[i][xyz];

  //fprintf(outfile,"ref_B with B1 corrected\n");
  //print_matrix(outfile, ref_B, ndB, 3);

  if (ndB>1) { /* move fragment B to place reference point B2 in correct location */
    /* Determine rotational angle and axis */
    v3d_eAB(ref_B[0], ref_B[1], e12);  /* v B1->B2 */
    v3d_eAB(ref_B[0], ref_B_final[1], e12b); /* v B1->B2_final */
    B_angle = acos(v3d_dot(e12b,e12));

    if (fabs(B_angle) > 1.0e-7) {
      v3d_cross_product(e12,e12b,erot);

      /* Move B to put B1 at origin */
      for (xyz=0; xyz<3; ++xyz)
        for (i=0; i<B->natom;++i)
          B_geom[i][xyz] -= ref_B[0][xyz];

      /* Rotate B */
      rotate_vecs(erot, B_angle, B_geom, B->natom);

      /* Move B back to coordinate system of A */
      for (xyz=0; xyz<3; ++xyz)
        for (i=0; i<B->natom;++i)
          B_geom[i][xyz] += ref_B[0][xyz];

      // recompute current B reference points
      zero_matrix(ref_B, ndB, 3);
      for (pts=0; pts<ndB; ++pts) 
        for (xyz=0; xyz<3; ++xyz)
          for (i=0; i<B->natom;++i)
            ref_B[pts][xyz] += weightB[pts][i] * B_geom[i][xyz];
    }
    //fprintf(outfile,"ref_B with B2 corrected\n");
    //print_matrix(outfile, ref_B, ndB, 3);
  }
  if (ndB==3) { // move fragment B to place reference point B3 in correct location
    // Determine rotational angle and axis
    v3d_eAB(ref_B[0], ref_B[1], erot);  /* B1 -> B2 is rotation axis */

    // Calculate B3-B1-B2-B3' torsion angle
    v3d_tors(ref_B[2], ref_B[0], ref_B[1], ref_B_final[2], B_angle);

    //fprintf(outfile,"B_angle: %15.10lf\n",B_angle);
    if (fabs(B_angle) > 1.0e-10) {

      // Move B to put B2 at origin
      for (xyz=0; xyz<3; ++xyz)
        for (i=0; i<B->natom;++i)
          B_geom[i][xyz] -= ref_B[1][xyz];

      rotate_vecs(erot, B_angle, B_geom, B->natom);

      // Translate B1 back to coordinate system of A
      for (xyz=0; xyz<3; ++xyz)
        for (i=0; i<B->natom;++i)
          B_geom[i][xyz] += ref_B[1][xyz];

      // update B reference points
      zero_matrix(ref_B, ndB, 3);
      for (pts=0; pts<ndB; ++pts)
        for (xyz=0; xyz<3; ++xyz)
          for (i=0; i<B->natom;++i)
            ref_B[pts][xyz] += weightB[pts][i] * B_geom[i][xyz];
    }
    //fprintf(outfile,"ref_B with B3 corrected\n");
    //print_matrix(outfile, ref_B, ndB, 3);
  }

  // check to see if desired reference points were obtained
  tval = 0.0;
  for (i=0; i<ndB; ++i)
    for (xyz=0; xyz<3; ++xyz)
      tval += (ref_B[i][xyz] - ref_B_final[i][xyz])*(ref_B[i][xyz] - ref_B_final[i][xyz]);
  tval = sqrt(tval);

  free_matrix(ref_A);
  free_matrix(ref_B);
  free_matrix(ref_B_final);

  fprintf(outfile,"\tDifference from target, |x_target - x_achieved| = %.2e\n",tval);

  if (tval > 1.0e-8) {
    fprintf(outfile,"\tUnsuccessful at orienting fragments!\n");
    fflush(outfile);
    return false;
  }
  else {
    fprintf(outfile,"\tSuccessfully oriented fragments.\n");
    fflush(outfile);
    B->set_geom(B_geom);
    free_matrix(B_geom);
    return true;
  }
}

/* Given the xyz coordinates for three points and R, theta, and phi, returns the
coordinates of a fourth point; angles in radians */
void zmat_point(double *A, double *B, double *C, double R_CD, double theta_BCD,
  double phi_ABCD, double *D) {

  double eAB[3],eBC[3],eX[3],eY[3], cosABC, sinABC;

  v3d_eAB(A,B,eAB); /* vector B->A */
  v3d_eAB(B,C,eBC); /* vector C->B */
  cosABC = -v3d_dot(eBC,eAB);

  sinABC = sqrt(1 - (cosABC * cosABC) );
  if ( (sinABC - 1.0e-14) < 0.0 ) {
    printf("Reference points cannot be colinear.");
    throw(INTCO_EXCEPT("Reference points cannot be colinear.", true));
  }

  v3d_cross_product(eAB,eBC,eY);
  for(int xyz=0;xyz<3;xyz++)
    eY[xyz] /= sinABC;
  v3d_cross_product(eY,eBC,eX);
  for (int xyz=0;xyz<3;xyz++)
    D[xyz] = C[xyz] + R_CD * ( - eBC[xyz] * cos(theta_BCD) +
                                 eX[xyz] * sin(theta_BCD) * cos(phi_ABCD) +
                                 eY[xyz] * sin(theta_BCD) * sin(phi_ABCD) );
  return;
}

/*!
** rotate_vecs(): Rotate a set of vectors around an arbitrary axis
**
** \brief Rotate a set of vectors around an arbitrary axis
** Vectors are rows of input matrix
**
** \param  w     double *  : axis to rotate around (wx, wy, wz) - gets normalized here
** \param  phi   double    : magnitude of rotation
** \param  v   double ** : points to rotate - column dim is 3; overwritten on exit
** \param  num_v  int       :
**
** Returns: none
**
** Rollin King, Feb. 2008
** \ingroup OPT
*/
void rotate_vecs(double *w, double phi, double **v, int num_v) {
  double **R, **v_new, wx, wy, wz, cp, norm;

  norm = sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

  w[0] /= norm; w[1] /= norm; w[2] /= norm;

  wx = w[0]; wy = w[1]; wz = w[2];
  cp = 1.0 - cos(phi);

  R = init_matrix(3,3);

  R[0][0] =     cos(phi) + wx*wx*cp;
  R[0][1] = -wz*sin(phi) + wx*wy*cp;
  R[0][2] =  wy*sin(phi) + wx*wz*cp;
  R[1][0] =  wz*sin(phi) + wx*wy*cp;
  R[1][1] =     cos(phi) + wy*wy*cp;
  R[1][2] = -wx*sin(phi) + wy*wz*cp;
  R[2][0] = -wy*sin(phi) + wx*wz*cp;
  R[2][1] =  wx*sin(phi) + wy*wz*cp;
  R[2][2] =     cos(phi) + wz*wz*cp;

  v_new = init_matrix(num_v,3);
  opt_matrix_mult(R, 0, v, 1, v_new, 1, 3, 3, num_v, 0);

  for (int i=0; i<num_v; ++i)
    for (int j=0; j<3; ++j)
      v[i][j] = v_new[i][j];

  free_matrix(v_new);
  free_matrix(R);
}

} // opt

