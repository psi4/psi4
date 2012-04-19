/*! \file stre.h
    \ingroup optking
    \brief STRE class declaration
*/

#include "frag.h"
#include "interfrag.h"
#include "print.h"
#include "v3d.h" // for H_guess

#include <sstream>

#define EXTERN
#include "globals.h"

namespace opt {

using namespace v3d;
using namespace std;

// constructor for given weight linear combination reference points
INTERFRAG::INTERFRAG(FRAG *A_in, FRAG *B_in, int A_index_in, int B_index_in,
    double **weightA_in, double **weightB_in, int ndA_in, int ndB_in, bool principal_axes) {

  A = A_in;
  B = B_in;
  A_index = A_index_in;
  B_index = B_index_in;
  weightA = weightA_in;
  weightB = weightB_in;
  ndA = ndA_in;
  ndB = ndB_in;

  double **inter_geom = init_matrix(6,3); // some rows may be unused

  // Create pseudo-fragment with atomic numbers at 6 positions.
  // The atomic numbers may affect Hessian guess routines at most.
  double *Z = init_array(6);
  for (int i=0; i<6; ++i) Z[i] = 6; // assume C for now
  inter_frag = new FRAG(6, Z, inter_geom);

  update_reference_points();

  add_coordinates_of_reference_pts();
}

// adds the coordinates connecting A2-A1-A0-B0-B1-B2
// sets D_on to indicate which ones (of the 6) are unusued
void INTERFRAG::add_coordinates_of_reference_pts(void) {
  STRE *one_stre = NULL;  // RAB
  BEND *one_bend = NULL;  // theta_A
  BEND *one_bend2 = NULL; // theta_B
  TORS *one_tors  = NULL; // tau
  TORS *one_tors2 = NULL; // phi_A
  TORS *one_tors3 = NULL; // phi_B

  // turn all coordinates on ; turn off unused ones below
  for (int i=0; i<6; ++i) D_on[i] = true;
 
  if (ndA == 3 && ndB == 3) {
    one_stre  = new STRE(2, 3);    // RAB
    one_bend  = new BEND(1, 2, 3); // theta_A
    one_bend2 = new BEND(2, 3, 4); // theta_B
    one_tors  = new TORS(1, 2, 3, 4); // tau
    one_tors2 = new TORS(0, 1, 2, 3); // phi_A
    one_tors3 = new TORS(2, 3, 4, 5); // phi_B
  }
  else if (ndA == 3 && ndB == 2) {
    D_on[5] = false; // no phi_B
    one_stre  = new STRE(2, 3);    // RAB
    one_bend  = new BEND(1, 2, 3); // theta_A
    one_bend2 = new BEND(2, 3, 4); // theta_B
    one_tors  = new TORS(1, 2, 3, 4); // tau
    one_tors2 = new TORS(0, 1, 2, 3); // phi_A
  }
  else if (ndA == 2 && ndB == 3) {
    D_on[4] = false; // no phi_A
    one_stre  = new STRE(2, 3);    // RAB
    one_bend  = new BEND(1, 2, 3); // theta_A
    one_bend2 = new BEND(2, 3, 4); // theta_B
    one_tors  = new TORS(1, 2, 3, 4); // tau
    one_tors3 = new TORS(2, 3, 4, 5); // phi_B
  }
  else if (ndA == 3 && ndB == 1) {
    D_on[2] = D_on[3] = D_on[5] = false; // no theta_B, tau, phi_B
    STRE *one_stre  = new STRE(2, 3);    // RAB
    BEND *one_bend  = new BEND(1, 2, 3); // theta_A
    TORS *one_tors2 = new TORS(0, 1, 2, 3); // phi_A
  }
  else if (ndA == 1 && ndB == 3) {
    D_on[1] = D_on[3] = D_on[4] = false; // no theta_A, tau, phi_A
    one_stre  = new STRE(2, 3);    // RAB
    one_bend2 = new BEND(2, 3, 4); // theta_B
    one_tors3 = new TORS(2, 3, 4, 5); // phi_B
  }
  else if (ndA == 2 && ndB == 2) {
    D_on[4] = D_on[5] = false; // no phi_A, phi_B
    one_stre  = new STRE(2, 3);    // RAB
    one_bend  = new BEND(1, 2, 3); // theta_A
    one_bend2 = new BEND(2, 3, 4); // theta_B
    one_tors  = new TORS(1, 2, 3, 4); // tau
  }
  else if (ndA == 2 && ndB == 1) {
    D_on[2] = D_on[4] = D_on[5] = false; // no theta_B, phi_A, phi_B
    one_stre  = new STRE(2, 3);    // RAB
    one_bend  = new BEND(1, 2, 3); // theta_A
    one_tors  = new TORS(1, 2, 3, 4); // tau
  }
  else if (ndA == 1 && ndB == 2) {
    D_on[1] = D_on[4] = D_on[5] = false; // no theta_A, phi_A, phi_B
    one_stre  = new STRE(2, 3);    // RAB
    one_bend2 = new BEND(2, 3, 4); // theta_B
    one_tors  = new TORS(1, 2, 3, 4); // tau
  }
  else if (ndA == 1 && ndB == 1) {
    D_on[1] = D_on[2] = D_on[3] = D_on[4] =  D_on[5] = false;
    one_stre  = new STRE(2, 3);    // RAB
  }
  else {
    throw(INTCO_EXCEPT("INTERFRAG::INTERFRAG Num. reference points on each fragment must be at least 1."));
  }

  // check if stretch is a H-bond or includes something H-bond like (remember stretch is
  // in general between linear combinations of atoms
  double ang;

  bool *is_XA = init_bool_array(g_natom());
  for (int a=0; a<A->g_natom(); ++a)
    if (A->Z[a] == 7 || A->Z[a] == 8 || A->Z[a] == 9 || A->Z[a] == 17)
      is_XA[a] = true;

  bool *is_XB = init_bool_array(g_natom());
  for (int b=0; b<B->g_natom(); ++b)
    if (B->Z[b] == 7 || B->Z[b] == 8 || B->Z[b] == 9 || B->Z[b] == 17)
      is_XB[b] = true;

  // Look for A[X]-A[H] ... B[Y]
  for (int h=0; h<A->natom; ++h) {
    if (weightA[0][h] != 0.0 && A->Z[h] == 1.0) { // H atom is (part of) A[0]
      for (int x=0; x<A->natom; ++x) {
        if (A->connectivity[x][h] && is_XA[x]) {  // electronegative X atom is present
          for (int y=0; y<B->natom; ++y) {
            if (weightB[0][y] != 0.0 && is_XB[y]) { // electronegative Y atom is part of B[0]
              if (v3d_angle(A->geom[x], A->geom[h], B->geom[y], ang)) { //check angle
                if (ang < _pi/2)
                  one_stre->make_hbond();
              }
            }
          }
        }
      }
    }
  }

  // Look for A[Y]...B[H]-B[X]
  for (int h=0; h<B->natom; ++h) {
    if (weightB[0][h] != 0.0 && B->Z[h] == 1.0) { // H atom is (part of) B[0]
      for (int x=0; x<B->natom; ++x) {
        if (B->connectivity[x][h] && is_XB[x]) {  // electronegative X atom is present
          for (int y=0; y<A->natom; ++y) {
            if (weightA[0][y] != 0.0 && is_XA[y]) { // electronegative Y atom is part of A[0]
              if (v3d_angle(B->geom[x], B->geom[h], A->geom[y], ang)) { //check angle
                if (ang < _pi/2)
                  one_stre->make_hbond();
              }
            }
          }
        }
      }
    }
  }

  if (Opt_params.interfragment_distance_inverse) {
    one_stre->make_inverse_stre(); 
    fprintf(outfile,"Using interfragment 1/R distance coordinate.\n");
  }
  if (one_stre->is_hbond())
    fprintf(outfile,"Detected H-bonding interfragment coordinate.\n");

  if (one_stre  != NULL) inter_frag->intcos.push_back(one_stre);
  if (one_bend  != NULL) inter_frag->intcos.push_back(one_bend);
  if (one_bend2 != NULL) inter_frag->intcos.push_back(one_bend2);
  if (one_tors  != NULL) inter_frag->intcos.push_back(one_tors);
  if (one_tors2 != NULL) inter_frag->intcos.push_back(one_tors2);
  if (one_tors3 != NULL) inter_frag->intcos.push_back(one_tors3);
}

// update location of reference points using given geometries
void INTERFRAG::update_reference_points(GeomType new_geom_A, GeomType new_geom_B) {

  zero_matrix(inter_frag->geom, 6, 3);

  if (!principal_axes) { // use fixed weights
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
  else { // using principal axes

    // first reference point on each fragment is the COM of each
    double *fragment_com = A->com();
    for (int xyz=0; xyz<3; ++xyz)
      inter_frag->geom[2][xyz] = fragment_com[xyz];

    double **axes, *moi;
    int i = A->principal_axes(new_geom_A, axes, moi);

    fprintf(outfile,"Number of principal axes returned is %d\n", i);

    for (i=0; i<ndA-1; ++i) // i can only be 0 or 1
      for (int xyz=0; xyz<3; ++xyz)
        inter_frag->geom[1-i][xyz] = fragment_com[xyz] + axes[i][xyz];

    free_array(moi);
    free_matrix(axes);
    free_array(fragment_com);

    fragment_com = B->com();
    for (int xyz=0; xyz<3; ++xyz)
      inter_frag->geom[3][xyz] = fragment_com[xyz];

    i = B->principal_axes(new_geom_B, axes, moi);

    fprintf(outfile,"Number of principal axes returned is %d\n", i);

    for (i=0; i<ndB-1; ++i) // i can only be 0 or 1
      for (int xyz=0; xyz<3; ++xyz)
        inter_frag->geom[4+i][xyz] = fragment_com[xyz] + axes[i][xyz];

    free_array(moi);
    free_matrix(axes);
    free_array(fragment_com);

    if (Opt_params.print_lvl >= 3) {
      fprintf(outfile,"\tndA: %d ; ndB: %d\n", ndA, ndB);
      fprintf(outfile,"\tReference points are at the following locations.\n");
      for (int i=2; i>2-ndA; --i)
        fprintf(outfile,"%15.10lf %15.10lf %15.10lf\n",
          inter_frag->geom[i][0], inter_frag->geom[i][1], inter_frag->geom[i][2]);
      for (int i=0; i<ndB; ++i)
        fprintf(outfile,"%15.10lf %15.10lf %15.10lf\n",
          inter_frag->geom[3+i][0], inter_frag->geom[3+i][1], inter_frag->geom[3+i][2]);
    }
  }
}

int INTERFRAG::g_nintco(void) const {
  int dim = 0;
  for (int i=0; i<6; ++i)
    if (D_on[i]) ++dim;
  return dim;
}

// freeze coordinate i if D_freeze[i]; index runs 0->6 as does D_on
/*
void INTERFRAG::freeze(bool *D_freeze) {
  int cnt = -1;
  for (int i=0; i<6; ++i) {
    if (D_on[i]) {
      ++cnt;
      if (D_freeze[i])
        inter_frag->intcos[cnt]->freeze();
    }
  }
}
*/

// freeze coordinate i; index is among the coordinates that are 'on'
void INTERFRAG::freeze(int index_to_freeze) {
  if (index_to_freeze<0 || index_to_freeze>g_nintco()) {
    fprintf(outfile,"INTERFRAG::freeze() : Invalid index %d\n", index_to_freeze);
    return;
  }
  inter_frag->intcos[index_to_freeze]->freeze();
}

void INTERFRAG::freeze(void) {
  inter_frag->freeze();
}

// is coordinate J frozen?  J runs over only active coordinates.
bool INTERFRAG::is_frozen(int J) { 
  if (J < 0 || J >= g_nintco())
    throw(INTCO_EXCEPT("INTERFRAG::is_frozen() index J runs only over active coordinates"));
  return inter_frag->intcos[J]->is_frozen();
}

// are all interfragment modes frozen?
bool INTERFRAG::is_frozen(void) {
  return inter_frag->is_frozen();
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

  if (!principal_axes) {

  int cnt=0, xyz;
  double **B_ref; // derivative of interfragment D wrt reference point position

  if (D_on[0]) {
    B_ref = inter_frag->intcos.at(cnt)->DqDx(inter_frag->geom); // RAB, returns (2,3)
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
    B_ref = inter_frag->intcos.at(cnt)->DqDx(inter_frag->geom); // theta_A, returns (3,3)
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
    B_ref = inter_frag->intcos.at(cnt)->DqDx(inter_frag->geom); // theta_B, returns (3,3)
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
    B_ref = inter_frag->intcos.at(cnt)->DqDx(inter_frag->geom); // tau, returns (4,3)
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
    B_ref = inter_frag->intcos.at(cnt)->DqDx(inter_frag->geom); // phi_A, returns (4,3)
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
    B_ref = inter_frag->intcos.at(cnt)->DqDx(inter_frag->geom); // phi_B, returns (4,3)
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

  }
  else { // principal axis
  //natomA natomB //double **B = init_matrix(g_nintco(), 3*(natomA+natomB));

/*

    double **A_u = init_matrix(3,3);
    double *A_lambda = init_array(3);
    int nA_lambda = Afrag->principal_axes(A, A_u, A_lambda);
    double A_mass = Afrag->masses();

    double **B_u = init_matrix(3,3);
    double *B_lambda = init_array(3);
    int nB_lambda = Bfrag->principal_axes(B, B_u, B_lambda);

    // First reference points are the centers of mass.  These are fixed weights.
    if (D_on[0]) {
      B_ref = inter_frag->intcos.at(cnt)->DqDx(inter_frag->geom); // RAB, returns (2,3)
      for (xyz=0; xyz<3; ++xyz) {
        for (int a=0; a<natomA; ++a)
          B[cnt][3*a + xyz] += A_mass[a]/A_mass_total * B_ref[0][xyz];
        for (int b=0; b<natomB; ++b)
          B[cnt][3*natomA + 3*b + xyz] += B_mass[b]/B_mass_total * B_ref[1][xyz];
      }
      free_matrix(B_ref);
      ++cnt;
    }

TODO

    free_matrix(A_u);
    free_matrix(B_u);
    free_array(A_lambda);
    free_array(B_lambda);
*/
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


void INTERFRAG::print_intcos(FILE *fp, int off_A, int off_B) const {
  fprintf(fp,"\t---Interfragment Coordinates Between Fragments %d and %d---\n", 
    A_index+1, B_index+1);
  fprintf(fp,"\t * Reference Points *\n");
  int cnt=0;
  for (int i=2; i>=0; --i, ++cnt) {
    if (i<ndA) {
      fprintf(fp,"\t\t %d A%d :", cnt+1, i+1);
      for (int j=0; j<A->g_natom(); ++j)
        if (weightA[i][j] != 0.0)
          fprintf(fp," %d/%5.3f", off_A+j+1, weightA[i][j]);
      fprintf(fp,"\n");
    }
  }
  for (int i=0; i<3; ++i, ++cnt) {
    if (i < ndB) {
      fprintf(fp,"\t\t %d B%d :", cnt+1, i+1);
      for (int j=0; j<B->g_natom(); ++j)
        if (weightB[i][j] != 0.0)
          fprintf(fp," %d/%5.3f", off_B+j+1, weightB[i][j]);
      fprintf(fp,"\n");
    }
  }
  fflush(fp);
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
  fflush(fp);
}

// TODO - fix this up later
std::string INTERFRAG::get_intco_definition(int coord_index, int off_A, int off_B) const {
  ostringstream iss;
  for (int i=0; i<ndA; ++i) {
    iss << "A" << i+1;
    for (int j=0; j<A->g_natom(); ++j)
      if (weightA[i][j] != 0.0)
        iss << j+1+off_A;
    iss << "\n";
  }
  for (int i=0; i<ndB; ++i) {
    iss << "B" << i+1;
    for (int j=0; j<B->g_natom(); ++j)
      if (weightB[i][j] != 0.0)
        iss << j+1+off_B;
    iss << "\n";
  }
  return iss.str();
}

// Make the initial Hessian guess for interfragment coordinates
double ** INTERFRAG::H_guess(void) {
  double **H;

  // use formulas from Fischer et al - not designed for interfragment modes
  if (Opt_params.interfragment_H == OPT_PARAMS::FISCHER_LIKE) {
    // H_guess uses intrafragment_H on inter_frag, so set and restore value
    OPT_PARAMS::INTRAFRAGMENT_HESSIAN i = Opt_params.intrafragment_H ;
    Opt_params.intrafragment_H = OPT_PARAMS::FISCHER;
    H = inter_frag->H_guess();
    Opt_params.intrafragment_H = i;
  }
  else { // DEFAULT
    H = init_matrix(inter_frag->g_nintco(), inter_frag->g_nintco());
    int cnt=0;
    double rAB;

    if (Opt_params.interfragment_distance_inverse)
      rAB = inter_frag->intcos[0]->value(inter_frag->geom);

    if (inter_frag->intcos[0]->is_hbond()) {

      H[cnt][cnt] = 0.03;
      if (Opt_params.interfragment_distance_inverse)
        H[cnt][cnt] *= pow(rAB,4);
      ++cnt;

      if (D_on[1]) { H[cnt][cnt] = 0.007; ++cnt; }
      if (D_on[2]) { H[cnt][cnt] = 0.007; ++cnt; }
      if (D_on[3]) { H[cnt][cnt] = 0.002; ++cnt; }
      if (D_on[4]) { H[cnt][cnt] = 0.002; ++cnt; }
      if (D_on[5]) { H[cnt][cnt] = 0.002; ++cnt; }
    }
    else {

      H[cnt][cnt] = 0.007;
      if (Opt_params.interfragment_distance_inverse)
        H[cnt][cnt] *= pow(rAB,4);
      ++cnt;

      if (D_on[1]) { H[cnt][cnt] = 0.003;  ++cnt; }
      if (D_on[2]) { H[cnt][cnt] = 0.003;  ++cnt; }
      if (D_on[3]) { H[cnt][cnt] = 0.001; ++cnt; }
      if (D_on[4]) { H[cnt][cnt] = 0.001; ++cnt; }
      if (D_on[5]) { H[cnt][cnt] = 0.001; ++cnt; }
    }
  
  }
  return H;
}

// return matrix with 1's on diagonal for frozen coordinates
double ** INTERFRAG::compute_constraints(void) const {
  double **C = init_matrix(g_nintco(), g_nintco());
  int cnt = 0;
  for (int i=0; i<6; ++i) {
    if (D_on[i]) {
      if (inter_frag->intcos[cnt++]->is_frozen())
        C[i][i] = 1.0;
    }
  }
  return C;
}

} // opt

