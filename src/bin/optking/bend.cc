/*! \file    bend.cc : Class for bending coordinate.
     \ingroup OPTKING
*/

#include "bend.h"

#include <sstream>
#include "v3d.h"
#include "physconst.h"
#include <math.h>
#include "linear_algebra.h"

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
  linear_bend = false;

  if (A_in == B_in || B_in == C_in || A_in == C_in)
    throw(INTCO_EXCEPT("BEND::BEND() Atoms defining bend are not unique."));

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
  double phi=0.0, tval;
  if (!linear_bend) {
    if ( ! v3d_angle(geom[s_atom[0]], geom[s_atom[1]], geom[s_atom[2]], phi) )
      throw(INTCO_EXCEPT("BEND::compute_val: could not compute angle",true));
  }
  else { //linear bending complement
    double u[3], v[3], w[3], tvect[3];

    v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[0]], u); // B->A
    v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[2]], v); // B->C
    v3d_normalize(u); // eBA
    v3d_normalize(v); // eBC

    // determine w' vector
    if (!v3d_is_parallel(u,v))
      v3d_cross_product(u,v,w);
    else {
      tvect[0] = 1; tvect[1] = -1; tvect[2] = 1;
      if (!v3d_is_parallel(u,tvect) && !v3d_is_parallel(v,tvect))
        v3d_cross_product(u,tvect,w);
      else {
        tvect[0] = -1; tvect[1] = 1; tvect[2] = 1;
        v3d_cross_product(u,tvect,w);
      }
    }
    v3d_normalize(w);

    double *origin = init_array(3);
    // linear complement is sum of 2 angles, u.w + w.v
    if (!v3d_angle(u, origin, w, phi))
      throw(INTCO_EXCEPT("BEND::value: could not compute linear bend",true));

    if (!v3d_angle(w, origin, v, tval))
      throw(INTCO_EXCEPT("BEND::value: could not compute linear bend",true));
    phi += tval;
  }
  return phi;
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
  double u[3], v[3], w[3];
  double tvect[3];
  double **dqdx = init_matrix(3,3);

  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[0]], u); // B->A
  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[2]], v); // B->C
  double Lu = v3d_norm(u); // RBA
  double Lv = v3d_norm(v); // RBC
  v3d_scm(1.0/Lu, u); // eBA
  v3d_scm(1.0/Lv, v); // eBC

  // determine w vector
  if (!v3d_is_parallel(u,v)) {
    if (linear_bend) { // use (-1)*(angle bisector of u and v) for perpendicular w vector
      v3d_axpy(1.0, geom[s_atom[0]], geom[s_atom[2]], tvect);
      v3d_scm(0.5, tvect);
      v3d_eAB(tvect, geom[s_atom[1]], w);
    }
    else
      v3d_cross_product(u,v,w); // use w = uXv
  }
  else { // 0-1-2 is linear ; arbitrarily choose a direction
    tvect[0] = 1; tvect[1] = -1; tvect[2] = 1;
    if (!v3d_is_parallel(u,tvect) && !v3d_is_parallel(v,tvect)) {
      v3d_cross_product(u,tvect,w);
    }
    else {
      tvect[0] = -1; tvect[1] = 1; tvect[2] = 1;
      v3d_cross_product(u,tvect,w);
    }
    if (linear_bend) { // use the complement vector w = w x u 
      v3d_cross_product(w,u,tvect);
      array_copy(tvect, w, 3);
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
  double u[3], v[3], w[3];
  double tvect[3]; 
  double **dq2dx2 = init_matrix(9,9);

  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[0]], u); // B->A
  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[2]], v); // B->C
  double Lu = v3d_norm(u); // RBA
  double Lv = v3d_norm(v); // RBC
  v3d_scm(1.0/Lu, u);  // eBA
  v3d_scm(1.0/Lv, v);  // eBC

  // determine w' vector
  if (!v3d_is_parallel(u,v)) {
    if (linear_bend) { // use (-1)*(angle bisector of u and v) for perpendicular w vector
      v3d_axpy(1.0, geom[s_atom[0]], geom[s_atom[2]], tvect);
      v3d_scm(0.5, tvect);
      v3d_eAB(tvect, geom[s_atom[1]], w);
    }
    else
      v3d_cross_product(u,v,w); // use w = uXv
  }
  else { // 0-1-2 is linear ; arbitrarily choose a direction
    tvect[0] = 1; tvect[1] = -1; tvect[2] = 1;
    if (!v3d_is_parallel(u,tvect) && !v3d_is_parallel(v,tvect)) 
      v3d_cross_product(u,tvect,w);
    else {
      tvect[0] = -1; tvect[1] = 1; tvect[2] = 1;
      v3d_cross_product(u,tvect,w);
    }
    if (linear_bend) { // use the complement vector w = w x u 
      v3d_cross_product(w,u,tvect);
      array_copy(tvect, w, 3);
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
  iss << get_definition_string(off);

  double val = value(geom);
  if (!s_frozen)
    fprintf(fp,"\t %-15s  =  %15.6lf\t%15.6lf\n", iss.str().c_str(), val, val/_pi*180.0);
  else
    fprintf(fp,"\t*%-15s  =  %15.6lf\t%15.6lf\n", iss.str().c_str(), val, val/_pi*180.0);
  fflush(fp);
}

// function to return string of coordinate definition
std::string BEND::get_definition_string(int off) const { 
  ostringstream iss(ostringstream::out); // create stream; allow output to it

  if (linear_bend)
    iss << "L(";
  else
    iss << "B(";

  iss << s_atom[0]+1+off << "," << s_atom[1]+1+off << "," << s_atom[2]+1+off << ")" << flush ;
  return iss.str();
}

void BEND::print_disp(FILE *fp, const double q_old, const double f_q,
    const double dq, const double q_new, int atom_offset) const {
  ostringstream iss(ostringstream::out); // create stream; allow output to it
  if (s_frozen) iss << "*";

  if (linear_bend)
    iss << "L(";
  else
    iss << "B(";

  iss << s_atom[0]+atom_offset+1 << "," << s_atom[1]+atom_offset+1 << "," << s_atom[2]+atom_offset+1 << ")" << flush ;

  fprintf(fp,"\t %-15s = %13.6lf%13.6lf%13.6lf%13.6lf\n",
    iss.str().c_str(), q_old/_pi*180.0, f_q*_hartree2aJ*_pi/180.0,
      dq/_pi*180.0, q_new/_pi*180.0);
  fflush(fp);
}

void BEND::print_intco_dat(FILE *fp, int off) const {
  if (linear_bend) {
    if (s_frozen)
      fprintf(fp, "L*%6d%6d%6d\n", s_atom[0]+1+off, s_atom[1]+1+off, s_atom[2]+1+off);
    else
      fprintf(fp, "L %6d%6d%6d\n", s_atom[0]+1+off, s_atom[1]+1+off, s_atom[2]+1+off);
  }
  else {
    if (s_frozen)
      fprintf(fp, "B*%6d%6d%6d\n", s_atom[0]+1+off, s_atom[1]+1+off, s_atom[2]+1+off);
    else
      fprintf(fp, "B %6d%6d%6d\n", s_atom[0]+1+off, s_atom[1]+1+off, s_atom[2]+1+off);
  }
  fflush(fp);
}

void BEND::print_s(FILE *fp, GeomType geom) const {

  if (linear_bend)
    fprintf(fp,"S vector for bend, L(%d %d %d): \n", s_atom[0]+1, s_atom[1]+1, s_atom[2]+1);
  else
    fprintf(fp,"S vector for bend, B(%d %d %d): \n", s_atom[0]+1, s_atom[1]+1, s_atom[2]+1);

  double **dqdx = DqDx(geom);
  fprintf(fp, "Atom 1: %12.8f %12.8f,%12.8f\n", dqdx[0][0], dqdx[0][1], dqdx[0][2]);
  fprintf(fp, "Atom 2: %12.8f %12.8f,%12.8f\n", dqdx[1][0], dqdx[1][1], dqdx[1][2]);
  fprintf(fp, "Atom 3: %12.8f %12.8f,%12.8f\n", dqdx[2][0], dqdx[2][1], dqdx[2][2]);
  free_matrix(dqdx);
  fflush(fp);
}

bool BEND::operator==(const SIMPLE & s2) const {
  if (bend_type != s2.g_type())
    return false;

  if (this->s_atom[1] != s2.g_atom(1))
    return false;

  if (this->is_linear_bend() != s2.is_linear_bend())
    return false;

  if (this->s_atom[0] == s2.g_atom(0) && this->s_atom[2] == s2.g_atom(2))
    return true;
  else if (this->s_atom[0] == s2.g_atom(2) && this->s_atom[2] == s2.g_atom(0))
    return true;
  else
    return false;
}

}

