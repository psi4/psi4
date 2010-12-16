/*! \file bend.cc
    \ingroup OPT10
    \brief bend class
*/

#include "bend.h"

#include <sstream>
#include "v3d.h"
#include "physconst.h"
#include <math.h>

#include "print.h"
#define EXTERN
#include "globals.h"

namespace opt {

using namespace v3d;
using namespace std;

// constructor - makes sure A<C
BEND::BEND(int A_in, int B_in, int C_in, bool freeze_in) : SIMPLE(bend_type, 3, freeze_in) {
  //fprintf(stdout,"constructing BEND A_in:%d B_in:%d C_in:%d, frozen %d\n",
  //  A_in, B_in, C_in, freeze_in);

  if (A_in == B_in || B_in == C_in || A_in == C_in)
    throw("BEND::BEND() Atoms defining bend are not unique.");

  s_atom[1] = B_in;
  if (A_in <= C_in) {
    s_atom[0] = A_in;
    s_atom[2] = C_in;
  }
  else {
    s_atom[0] = C_in;
    s_atom[2] = A_in;
  }
}

// compute angle and store value in radians
double BEND::value(GeomType geom) const {
  double tval;
  if ( ! v3d_angle(geom[s_atom[0]], geom[s_atom[1]], geom[s_atom[2]], tval) )
    throw("BEND::compute_val: could not compute angle");
  return tval;
}

inline int zeta(const int a, const int m, const int n) {
  if (a == m) return 1;
  else if (a == n) return -1;
  else return 0;
}

inline int delta(const int i, const int j) {
  if (i == j) return 1;
  else return 0;
}

double ** BEND::DqDx(GeomType geom) const {
  double **dqdx = init_matrix(3,3);

  double u[3], v[3], w[3];
  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[0]], u); // B->A
  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[2]], v); // B->C
  double Lu = v3d_norm(u); // RBA
  double Lv = v3d_norm(v); // RBC
  v3d_scm(1.0/Lu, u); // eBA
  v3d_scm(1.0/Lv, v); // eBC

  // determine w' vector
  if (!v3d_is_parallel(u,v))
    v3d_cross_product(u,v,w);
  else {
    double tvect[3];
    tvect[0] = 1; tvect[1] = -1; tvect[2] = 1;
    if (!v3d_is_parallel(u,tvect) && !v3d_is_parallel(v,tvect))
      v3d_cross_product(u,tvect,w);
    else {
      tvect[0] = -1; tvect[1] = 1; tvect[2] = 1;
      v3d_cross_product(u,tvect,w);
    }
  }
  v3d_normalize(w);

  double uXw[3],wXv[3];
  v3d_cross_product(u, w, uXw);
  v3d_cross_product(w, v, wXv);

  for (int a=0; a<3; ++a)
    for (int i=0; i<3; ++i)
      dqdx[a][i] = zeta(a,0,1)*uXw[i]/Lu + zeta(a,2,1)*wXv[i]/Lv;

  return dqdx;
}

double ** BEND::Dq2Dx2(GeomType geom) const {
  double **dq2dx2 = init_matrix(9,9);

  double u[3], v[3], w[3];
  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[0]], u); // B->A
  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[2]], v); // B->C
  double Lu = v3d_norm(u); // RBA
  double Lv = v3d_norm(v); // RBC
  v3d_scm(1.0/Lu, u);  // eBA
  v3d_scm(1.0/Lv, v);  // eBC

  // determine w' vector
  if (!v3d_is_parallel(u,v)) 
    v3d_cross_product(u,v,w);
  else {
    double tvect[3]; 
    tvect[0] = 1; tvect[1] = -1; tvect[2] = 1;
    if (!v3d_is_parallel(u,tvect) && !v3d_is_parallel(v,tvect)) 
      v3d_cross_product(u,tvect,w);
    else {
      tvect[0] = -1; tvect[1] = 1; tvect[2] = 1;
      v3d_cross_product(u,tvect,w);
    }
  }
  v3d_normalize(w);

  // compute first derivatives
  double uXw[3],wXv[3];
  v3d_cross_product(u, w, uXw);
  v3d_cross_product(w, v, wXv);

  double **dqdx = init_matrix(3,3);
  for (int a=0; a<3; ++a)
    for (int i=0; i<3; ++i)
      dqdx[a][i] = zeta(a,0,1)*uXw[i]/Lu + zeta(a,2,1)*wXv[i]/Lv;

  double cos_q = v3d_dot(u,v);
  if (1.0-cos_q*cos_q <= 1.0e-12) return dqdx; // leave 2nd derivatives empty - sin 0 = 0 in denominator
  double sin_q = sqrt(1.0 - cos_q*cos_q);

  double tval;
  for (int a=0; a<3; ++a)
    for (int i=0; i<3; ++i) //i = a_xyz
      for (int b=0; b<3; ++b)
        for (int j=0; j<3; ++j) { //j=b_xyz
          tval = zeta(a,0,1)*zeta(b,0,1)*(u[i]*v[j]+u[j]*v[i]-3*u[i]*u[j]*cos_q+delta(i,j)*cos_q) / (Lu*Lu*sin_q);

          tval +=zeta(a,2,1)*zeta(b,2,1)*(v[i]*u[j]+v[j]*u[i]-3*v[i]*v[j]*cos_q+delta(i,j)*cos_q) / (Lv*Lv*sin_q);

          tval +=zeta(a,0,1)*zeta(b,2,1)*(u[i]*u[j]+v[j]*v[i]-u[i]*v[j]*cos_q-delta(i,j)) / (Lu*Lv*sin_q);

          tval +=zeta(a,2,1)*zeta(b,0,1)*(v[i]*v[j]+u[j]*u[i]-v[i]*u[j]*cos_q-delta(i,j)) / (Lu*Lv*sin_q);

          tval -= cos_q / sin_q * dqdx[a][i] * dqdx[b][j];

          dq2dx2[3*a+i][3*b+j] = tval;
        }

  free_matrix(dqdx);

  return dq2dx2;
}


void BEND::print(FILE *fp, GeomType geom, int off) const {
  ostringstream iss(ostringstream::out); // create stream; allow output to it
  iss << "A(" << s_atom[0]+1+off << "," << s_atom[1]+1+off << "," << s_atom[2]+1+off << ")" ;
  double val = value(geom);
  if (!s_frozen)
    fprintf(fp,"\t %-15s  =  %15.6lf\t%15.6lf\n", iss.str().c_str(), val, val/_pi*180.0);
  else
    fprintf(fp,"\t*%-15s  =  %15.6lf\t%15.6lf\n", iss.str().c_str(), val, val/_pi*180.0);
}

void BEND::print_disp(FILE *fp, const double q_old, const double f_q,
    const double dq, const double q_new, int atom_offset) const {
  ostringstream iss(ostringstream::out); // create stream; allow output to it
  if (s_frozen) iss << "*";
  iss << "A(" << s_atom[0]+atom_offset+1 << "," << s_atom[1]+atom_offset+1 << "," << s_atom[2]+atom_offset+1 << ")" ;
  fprintf(fp,"\t %-15s = %13.6lf%13.6lf%13.6lf%13.6lf\n",
    iss.str().c_str(), q_old/_pi*180.0, f_q*_hartree2aJ*_pi/180.0,
      dq/_pi*180.0, q_new/_pi*180.0);
}

void BEND::print_intco_dat(FILE *fp, int off) const {
  if (s_frozen)
    fprintf(fp, "B*%6d%6d%6d\n", s_atom[0]+1+off, s_atom[1]+1+off, s_atom[2]+1+off);
  else
    fprintf(fp, "B %6d%6d%6d\n", s_atom[0]+1+off, s_atom[1]+1+off, s_atom[2]+1+off);
}

void BEND::print_s(FILE *fp, GeomType geom) const {
  fprintf(fp,"S vector for bend, A(%d %d %d): \n",
    s_atom[0]+1, s_atom[1]+1, s_atom[2]+1);
  double **dqdx = DqDx(geom);
  fprintf(fp, "Atom 1: %12.8f %12.8f,%12.8f\n", dqdx[0][0], dqdx[0][1], dqdx[0][2]);
  fprintf(fp, "Atom 2: %12.8f %12.8f,%12.8f\n", dqdx[1][0], dqdx[1][1], dqdx[1][2]);
  fprintf(fp, "Atom 3: %12.8f %12.8f,%12.8f\n", dqdx[2][0], dqdx[2][1], dqdx[2][2]);
  free_matrix(dqdx);
}

bool BEND::operator==(const SIMPLE & s2) const {
  fprintf(stdout,"BEND::operator==\n");
  if (bend_type != s2.g_type())
    return false;

  if (this->s_atom[1] != s2.g_atom(1))
    return false;

  if (this->s_atom[0] == s2.g_atom(0) && this->s_atom[2] == s2.g_atom(2))
    return true;
  else if (this->s_atom[0] == s2.g_atom(2) && this->s_atom[2] == s2.g_atom(0))
    return true;
  else
    return false;
}

}

