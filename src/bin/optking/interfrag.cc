/*! \file stre.h
    \ingroup OPT10
    \brief STRE class declaration
*/

#include "frag.h"
#include "interfrag.h"
#include "print.h"
#include "v3d.h" // for H_guess

#define EXTERN
#include "globals.h"

namespace opt {

INTERFRAG::INTERFRAG(FRAG *A_in, FRAG *B_in, int A_index_in, int B_index_in,
    double **weightA_in, double **weightB_in, int ndA_in, int ndB_in) {

  A = A_in;
  B = B_in;
  A_index = A_index_in;
  B_index = B_index_in;
  weightA = weightA_in;
  weightB = weightB_in;
  ndA = ndA_in;
  ndB = ndB_in;

  double **inter_geom = init_matrix(6,3); // some rows may be unused

  // create pseudo-fragment with atomic numbers at 6
  // the atomic numbers may only affect Hessian guess routines
  double *Z = init_array(6);
  for (int i=0; i<6; ++i) Z[i] = 6;
  inter_frag = new FRAG(6, Z, inter_geom);

  update_reference_points();

  if (ndA == 3 && ndB == 3) {
    for (int i=0; i<6; ++i) D_on[i] = true;

    STRE *one_stre  = new STRE(2, 3);    // RAB
    BEND *one_bend  = new BEND(1, 2, 3); // theta_A
    BEND *one_bend2 = new BEND(2, 3, 4); // theta_B
    TORS *one_tors  = new TORS(1, 2, 3, 4); // tau
    TORS *one_tors2 = new TORS(0, 1, 2, 3); // phi_A
    TORS *one_tors3 = new TORS(2, 3, 4, 5); // phi_B

    inter_frag->intcos.push_back(one_stre);
    inter_frag->intcos.push_back(one_bend);
    inter_frag->intcos.push_back(one_bend2);
    inter_frag->intcos.push_back(one_tors);
    inter_frag->intcos.push_back(one_tors2);
    inter_frag->intcos.push_back(one_tors3);
  }
}

// update location of reference points using given geometries
void INTERFRAG::update_reference_points(GeomType new_geom_A, GeomType new_geom_B) {
  for (int i=0; i<6; ++i)
    for (int xyz=0; xyz<3; ++xyz)
      inter_frag->geom[i][xyz] = 0.0;

  for (int xyz=0; xyz<3; ++xyz) {
    for (int a=0; a<A->g_natom(); ++a) {
      inter_frag->geom[0][xyz] += weightA[2][a] * new_geom_A[a][xyz];
      inter_frag->geom[1][xyz] += weightA[1][a] * new_geom_A[a][xyz];
      inter_frag->geom[2][xyz] += weightA[0][a] * new_geom_A[a][xyz];
    }
    for (int b=0; b<B->g_natom(); ++b) {
      inter_frag->geom[3][xyz] += weightB[0][b] * new_geom_B[b][xyz];
      inter_frag->geom[4][xyz] += weightB[1][b] * new_geom_B[b][xyz];
      inter_frag->geom[5][xyz] += weightB[2][b] * new_geom_B[b][xyz];
    }
  }
}

int INTERFRAG::g_nintco(void) const {
  int dim = 0;
  for (int i=0; i<6; ++i)
    if (D_on[i]) ++dim;
  return dim;
}

// compute and return coordinate values - using given fragment geometries
double * INTERFRAG::intco_values(GeomType new_geom_A, GeomType new_geom_B) {
  update_reference_points(new_geom_A, new_geom_B);
  double *q = init_array(g_nintco());

  for (int i=0; i<g_nintco(); ++i)
    q[i] = inter_frag->intcos.at(i)->value(inter_frag->geom);

  return q;
}

// returns B matrix (internals by 3*natomA + 3*natomB)
double ** INTERFRAG::compute_B(GeomType new_geom_A, GeomType new_geom_B) {

  update_reference_points(new_geom_A, new_geom_B);
  int natomA = A->natom;
  int natomB = B->natom;

  double **B = init_matrix(g_nintco(), 3*(natomA+natomB));

  int cnt=0, xyz;
  double **B_ref; // derivative of interfragment D wrt reference point position

  if (D_on[0]) {
    B_ref = inter_frag->intcos.at(0)->DqDx(inter_frag->geom); // RAB, returns (2,3)
    for (xyz=0; xyz<3; ++xyz) {
      for (int a=0; a<natomA; ++a)
        B[cnt][3*a + xyz] += weightA[0][a] * B_ref[0][xyz];
      for (int b=0; b<natomB; ++b)
        B[cnt][3*natomA + 3*b + xyz] += weightB[0][b] * B_ref[1][xyz];
    }
    free_matrix(B_ref);
    ++cnt;
  }

  if (D_on[1]) {
    B_ref = inter_frag->intcos.at(1)->DqDx(inter_frag->geom); // theta_A, returns (3,3)
    for (xyz=0; xyz<3; ++xyz)  {
      for (int a=0; a<natomA; ++a) {
        B[cnt][3*a + xyz] += weightA[1][a] * B_ref[0][xyz];
        B[cnt][3*a + xyz] += weightA[0][a] * B_ref[1][xyz];
      }
      for (int b=0; b<natomB; ++b)
        B[cnt][3*natomA + 3*b + xyz] += weightB[0][b] * B_ref[2][xyz];
    }
    free_matrix(B_ref);
    ++cnt;
  }

  if (D_on[2]) {
    B_ref = inter_frag->intcos.at(2)->DqDx(inter_frag->geom); // theta_B, returns (3,3)
    for (xyz=0; xyz<3; ++xyz)  {
      for (int a=0; a<natomA; ++a)
        B[cnt][3*a + xyz] += weightA[0][a] * B_ref[0][xyz];
      for (int b=0; b<natomB; ++b) {
        B[cnt][3*natomA + 3*b + xyz] += weightB[0][b] * B_ref[1][xyz];
        B[cnt][3*natomA + 3*b + xyz] += weightB[1][b] * B_ref[2][xyz];
      }
    }
    free_matrix(B_ref);
    ++cnt;
  }

  if (D_on[3]) {
    B_ref = inter_frag->intcos.at(3)->DqDx(inter_frag->geom); // tau, returns (4,3)
    for (xyz=0; xyz<3; ++xyz)  {
      for (int a=0; a<natomA; ++a) {
        B[cnt][3*a + xyz] += weightA[1][a] * B_ref[0][xyz];
        B[cnt][3*a + xyz] += weightA[0][a] * B_ref[1][xyz];
      }
      for (int b=0; b<natomB; ++b) {
        B[cnt][3*natomB + 3*b + xyz] += weightB[0][b] * B_ref[2][xyz];
        B[cnt][3*natomB + 3*b + xyz] += weightB[1][b] * B_ref[3][xyz];
      }
    }
    free_matrix(B_ref);
    ++cnt;
  }

  if (D_on[4]) {
    B_ref = inter_frag->intcos.at(4)->DqDx(inter_frag->geom); // phi_A, returns (4,3)
    for (xyz=0; xyz<3; ++xyz)  {
      for (int a=0; a<natomA; ++a) {
        B[cnt][3*a + xyz] += weightA[2][a] * B_ref[0][xyz];
        B[cnt][3*a + xyz] += weightA[1][a] * B_ref[1][xyz];
        B[cnt][3*a + xyz] += weightA[0][a] * B_ref[2][xyz];
      }
      for (int b=0; b<natomB; ++b)
        B[cnt][3*natomA + 3*b + xyz] += weightB[0][b] * B_ref[3][xyz];
    }
    free_matrix(B_ref);
    ++cnt;
  }

  if (D_on[5]) {
    B_ref = inter_frag->intcos.at(5)->DqDx(inter_frag->geom); // phi_B, returns (4,3)
    for (xyz=0; xyz<3; ++xyz)  {
      for (int a=0; a<natomA; ++a)
        B[cnt][3*a + xyz] += weightA[0][a] * B_ref[0][xyz];
      for (int b=0; b<natomB; ++b) {
        B[cnt][3*natomB + 3*b + xyz] += weightB[0][b] * B_ref[1][xyz];
        B[cnt][3*natomB + 3*b + xyz] += weightB[1][b] * B_ref[2][xyz];
        B[cnt][3*natomB + 3*b + xyz] += weightB[2][b] * B_ref[3][xyz];
      }
    }
    free_matrix(B_ref);
    ++cnt;
  }

  return B;
}

// returns derivative B matrix for one internal, returns (3*natomA + 3*natomB) by same
double **INTERFRAG::compute_derivative_B(int J, GeomType new_geom_A, GeomType new_geom_B) {

  update_reference_points(new_geom_A, new_geom_B);

  int nA = A->natom;
  int nB = B->natom;

  int cnt=0, xyz_a, xyz_b;
  double **B2_ref; // 2nd derivative of interfragment D wrt reference point position

  int P, Pp;             // P and P' are fragments 0 or 1 (A or B)
  int P_nref, Pp_nref;   // # of reference points of P and P'
  int P_natom, Pp_natom; // # of atoms in fragments P and P'
  int K, Kp;             // reference atoms 0-2 on fragment P and P'
  int P_off, Pp_off;     // offset for 1st atom of fragment list for both fragments (0 or 3*nA)
  int P_off_ref, Pp_off_ref; // offset for 1st reference atom for fragment (0 or 3*ndA);
  double **P_weight, **Pp_weight;

  double **B2 = init_matrix(3*(nA+nB) , 3*(nA+nB));  // d^2(D)/d(r)^2 (sides are dA, then dB)
  int i, i_xyz, j, j_xyz;

  B2_ref = inter_frag->intcos.at(J)->Dq2Dx2(inter_frag->geom); // d^2(D)/d(reference point)^2

  for (i_xyz=0; i_xyz<3; ++i_xyz) {
    for (j_xyz=0; j_xyz<3; ++j_xyz) {

      if (J == 0) { // RAB, A0-B0

        for (i=0; i<nA; ++i)
          for (j=0; j<nA; ++j) // d^2(R) / d_Ai d_Aj
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[0][i] * weightA[0][j] * B2_ref[i_xyz][j_xyz];

        for (i=0; i<nB; ++i)
          for (j=0; j<nB; ++j) // d^2(R) / d_Bi d_Bj
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[0][i] * weightB[0][j] * B2_ref[3+i_xyz][3+j_xyz];

        for (i=0; i<nB; ++i)
          for (j=0; j<nA; ++j) // d^2(R) / d_Bi d_Aj
            B2[3*nA+3*i+i_xyz][3*j+j_xyz] += weightB[0][i] * weightA[0][j] * B2_ref[3+i_xyz][j_xyz];

      }
      else if (J == 1) { //theta_A, angle is A1-A0-B0

        for (i=0; i<nA; ++i)
          for (j=0; j<nA; ++j) { // d^2(D) / d_Ai[0,1] d_Aj[0,1]
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[1][i] * weightA[1][j] * B2_ref[i_xyz][j_xyz];
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[0][i] * weightA[0][j] * B2_ref[3+i_xyz][3+j_xyz];
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[1][i] * weightA[0][j] * B2_ref[3+i_xyz][j_xyz];
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[0][i] * weightA[1][j] * B2_ref[i_xyz][3+j_xyz];
          }

        for (i=0; i<nB; ++i)
          for (j=0; j<nB; ++j) // d^2(D) / d_Bi[0] d_Bj[0]
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[0][i] * weightB[0][j] * B2_ref[6+i_xyz][6+j_xyz];

        for (i=0; i<nB; ++i)
          for (j=0; j<nA; ++j) {// d^2(D) / d_Bi[0] d_Aj[1,0]
            B2[3*nA+3*i+i_xyz][3*j+j_xyz] += weightB[0][i] * weightA[0][j] * B2_ref[6+i_xyz][3+j_xyz];
            B2[3*nA+3*i+i_xyz][3*j+j_xyz] += weightB[0][i] * weightA[1][j] * B2_ref[6+i_xyz][j_xyz];
          }

      }
      else if (J == 2) { //theta_B, angle is A0-B0-B1

        for (i=0; i<nA; ++i)
          for (j=0; j<nA; ++j) // d^2(D) / d_Ai[0] d_Aj[0]
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[0][i] * weightA[0][j] * B2_ref[i_xyz][j_xyz];

        for (i=0; i<nB; ++i)
          for (j=0; j<nB; ++j) { // d^2(D) / d_Bi[0,1] d_Bj[0,1]
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[0][i] * weightB[0][j] * B2_ref[3+i_xyz][3+j_xyz];
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[1][i] * weightB[1][j] * B2_ref[6+i_xyz][6+j_xyz];
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[1][i] * weightB[0][j] * B2_ref[6+i_xyz][3+j_xyz];
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[0][i] * weightB[1][j] * B2_ref[3+i_xyz][6+j_xyz];
          }

        for (i=0; i<nB; ++i)
          for (j=0; j<nA; ++j) {// d^2(D) / d_Bi[0,1] d_Ai[0]
            B2[3*nA+3*i+i_xyz][3*j+j_xyz] += weightB[0][i] * weightA[0][j] * B2_ref[3+i_xyz][j_xyz];
            B2[3*nA+3*i+i_xyz][3*j+j_xyz] += weightB[1][i] * weightA[0][j] * B2_ref[6+i_xyz][j_xyz];
          }

      }
      else if (J == 3) { //tau, angle is A1-A0-B0-B1

//fprintf(outfile,"tau Dq2\n");
//print_matrix(outfile,B2_ref, 12, 12);

        for (i=0; i<nA; ++i)
          for (j=0; j<=nA; ++j) { // d^2(D) / d_Ai[0,1] d_Aj[0,1]
            B2[3*i+i_xyz][3*j+j_xyz]  = weightA[1][i] * weightA[1][j] * B2_ref[i_xyz][j_xyz];
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[0][i] * weightA[0][j] * B2_ref[3+i_xyz][3+j_xyz];
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[0][i] * weightA[1][j] * B2_ref[3+i_xyz][j_xyz];
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[1][i] * weightA[0][j] * B2_ref[i_xyz][3+j_xyz];
          }

        for (i=0; i<nB; ++i)
          for (j=0; j<nB; ++j) { // d^2(D) / d_Bi[0,1] d_Bj[0,1]
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[0][i] * weightB[0][j] * B2_ref[6+i_xyz][6+j_xyz];
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[1][i] * weightB[1][j] * B2_ref[9+i_xyz][9+j_xyz];
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[1][i] * weightB[0][j] * B2_ref[9+i_xyz][6+j_xyz];
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[0][i] * weightB[1][j] * B2_ref[6+i_xyz][9+j_xyz];
          }

        for (i=0; i<nB; ++i)
          for (j=0; j<nA; ++j) {// d^2(D) / d_Bi[0,1] d_Aj[0,1]
            B2[3*nA+3*i+i_xyz][3*j+j_xyz] += weightB[0][i] * weightA[0][j] * B2_ref[6+i_xyz][3+j_xyz];
            B2[3*nA+3*i+i_xyz][3*j+j_xyz] += weightB[0][i] * weightA[1][j] * B2_ref[6+i_xyz][j_xyz];
            B2[3*nA+3*i+i_xyz][3*j+j_xyz] += weightB[1][i] * weightA[0][j] * B2_ref[9+i_xyz][3+j_xyz];
            B2[3*nA+3*i+i_xyz][3*j+j_xyz] += weightB[1][i] * weightA[1][j] * B2_ref[9+i_xyz][j_xyz];
          }

      }
      else if (J == 4) { //phi_A, angle is A2-A1-A0-B0

        for (i=0; i<nA; ++i)
          for (j=0; j<nA; ++j) { // d^2(D) / d_Ai[0,1,2] d_Aj[0,1,2]
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[2][i] * weightA[2][j] * B2_ref[i_xyz][j_xyz];
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[1][i] * weightA[1][j] * B2_ref[3+i_xyz][3+j_xyz];
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[0][i] * weightA[0][j] * B2_ref[6+i_xyz][6+j_xyz];
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[1][i] * weightA[0][j] * B2_ref[3+i_xyz][6+j_xyz];
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[0][i] * weightA[1][j] * B2_ref[6+i_xyz][3+j_xyz];
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[2][i] * weightA[0][j] * B2_ref[i_xyz][6+j_xyz];
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[0][i] * weightA[2][j] * B2_ref[6+i_xyz][j_xyz];
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[2][i] * weightA[1][j] * B2_ref[i_xyz][3+j_xyz];
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[1][i] * weightA[2][j] * B2_ref[3+i_xyz][j_xyz];
          }

        for (i=0; i<nB; ++i)
          for (j=0; j<=i; ++j) // d^2(D) / d_Bi[0] d_Bj[0]
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[0][i] * weightB[0][j] * B2_ref[9+i_xyz][9+j_xyz];

        for (i=0; i<nB; ++i)
          for (j=0; j<nA; ++j) {// d^2(D) / d_Bi[0] d_Aj[0,1,2]
            B2[3*nA+3*i+i_xyz][3*j+j_xyz] += weightB[0][i] * weightA[2][j] * B2_ref[9+i_xyz][j_xyz];
            B2[3*nA+3*i+i_xyz][3*j+j_xyz] += weightB[0][i] * weightA[1][j] * B2_ref[9+i_xyz][3+j_xyz];
            B2[3*nA+3*i+i_xyz][3*j+j_xyz] += weightB[0][i] * weightA[0][j] * B2_ref[9+i_xyz][6+j_xyz];
          }

      }
      else if (J == 5) { //phi_B, angle is A0-B0-B1-B2

        for (i=0; i<nA; ++i)
          for (j=0; j<nA; ++j) // d^2(D) / d_Ai[0] d_Aj[0]
            B2[3*i+i_xyz][3*j+j_xyz] += weightA[0][i] * weightA[0][j] * B2_ref[i_xyz][j_xyz];

        for (i=0; i<nB; ++i)
          for (j=0; j<nB; ++j) { // d^2(D) / d_Bi[0,1,2] d_Bj[0,1,2]
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[0][i] * weightB[0][j] * B2_ref[3+i_xyz][3+j_xyz];
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[1][i] * weightB[1][j] * B2_ref[6+i_xyz][6+j_xyz];
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[2][i] * weightB[2][j] * B2_ref[9+i_xyz][9+j_xyz];
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[1][i] * weightB[0][j] * B2_ref[6+i_xyz][3+j_xyz];
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[0][i] * weightB[1][j] * B2_ref[3+i_xyz][6+j_xyz];
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[2][i] * weightB[0][j] * B2_ref[9+i_xyz][3+j_xyz];
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[0][i] * weightB[2][j] * B2_ref[3+i_xyz][9+j_xyz];
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[2][i] * weightB[1][j] * B2_ref[9+i_xyz][6+j_xyz];
            B2[3*nA+3*i+i_xyz][3*nA+3*j+j_xyz] += weightB[1][i] * weightB[2][j] * B2_ref[6+i_xyz][9+j_xyz];
          }

        for (i=0; i<nB; ++i)
          for (j=0; j<nA; ++j) {// d^2(D) / d_Bi[0,1,2] d_Aj[0]
            B2[3*nA+3*i+i_xyz][3*j+j_xyz] += weightB[0][i] * weightA[0][j] * B2_ref[3+i_xyz][j_xyz];
            B2[3*nA+3*i+i_xyz][3*j+j_xyz] += weightB[1][i] * weightA[0][j] * B2_ref[6+i_xyz][j_xyz];
            B2[3*nA+3*i+i_xyz][3*j+j_xyz] += weightB[2][i] * weightA[0][j] * B2_ref[9+i_xyz][j_xyz];
          }
      }
    }
  }

  // for now, fill in symmetric matrix
  for (i=0; i<3*nB; ++i) // BA->AB
    for (j=0; j<3*nA; ++j)
      B2[j][3*nA+i] = B2[3*nA+i][j];

  free_matrix(B2_ref);

  return B2;
}


void INTERFRAG::print_intcos(FILE *fp) const {
  fprintf(fp,"\t---Interfragment Coordinates Between Fragments %d and %d---\n", 
    A_index+1, B_index+1);
  fprintf(fp,"\t * Reference Points *\n");
  int cnt=0;
  for (int i=2; i>=0; --i, ++cnt) {
    if (i<ndA) {
      fprintf(fp,"\t\t %d A%d :", cnt+1, i+1);
      for (int j=0; j<A->g_natom(); ++j)
        if (weightA[i][j] != 0.0) fprintf(fp," %d/%5.3f", j+1, weightA[i][j]);
      fprintf(fp,"\n");
    }
  }
  for (int i=0; i<3; ++i, ++cnt) {
    if (i < ndB) {
      fprintf(fp,"\t\t %d B%d :", cnt+1, i+1);
      for (int j=0; j<B->g_natom(); ++j)
        if (weightB[i][j] != 0.0) fprintf(fp," %d/%5.3f", j+1, weightB[i][j]);
      fprintf(fp,"\n");
    }
  }
  inter_frag->print_intcos(fp);
}

void INTERFRAG::print_intco_dat(FILE *fp, int off_A, int off_B) const {
  for (int i=0; i<ndA; ++i) {
    fprintf(fp,"A%d",i+1);
    for (int j=0; j<A->g_natom(); ++j)
      if (weightA[i][j] != 0.0) fprintf(fp," %d", j+1+off_A);
    fprintf(fp,"\n");
  }
  for (int i=0; i<ndB; ++i) {
    fprintf(fp,"B%d",i+1);
    for (int j=0; j<B->g_natom(); ++j)
      if (weightB[i][j] != 0.0) fprintf(fp," %d", j+1+off_B);
    fprintf(fp,"\n");
  }
}

double ** INTERFRAG::H_guess(void) {
  double **H;
  if (Opt_params.interfragment_H == OPT_PARAMS::FISCHER_LIKE) {
        // H_guess uses Opt_params.intrafragment_H on inter_frag, so set and restore value
        OPT_PARAMS::INTRAFRAGMENT_HESSIAN i = Opt_params.intrafragment_H ;
        Opt_params.intrafragment_H = OPT_PARAMS::FISCHER;
        H = inter_frag->H_guess();
        Opt_params.intrafragment_H = i;
  }
  else {
    H = init_matrix(inter_frag->g_nintco(), inter_frag->g_nintco());
    int cnt = 0;
    if (D_on[0]) { // FISCHER formula for stretch
      double rAB = v3d::v3d_dist(inter_frag->geom[2], inter_frag->geom[3]);
      H[cnt][cnt] = 0.3601 * exp(-1.944*(rAB - 4.0));
      if (H[cnt][cnt] > 3) H[cnt][cnt] = 3;
      if (Opt_params.frag_dist_rho) H[cnt][cnt] *= pow(rAB,4);
      ++cnt;
    }
    if (D_on[1]) { H[cnt][cnt] = 0.001;  ++cnt; }
    if (D_on[2]) { H[cnt][cnt] = 0.001;  ++cnt; }
    if (D_on[3]) { H[cnt][cnt] = 0.0005; ++cnt; }
    if (D_on[4]) { H[cnt][cnt] = 0.0005; ++cnt; }
    if (D_on[5]) { H[cnt][cnt] = 0.0005; ++cnt; }
  }
  return H;
}

} // opt

