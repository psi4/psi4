/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file    bend.cc : Class for bending coordinate.
     \ingroup OPTKING
*/

#include "bend.h"

#include <sstream>
#include <math.h>

#include "v3d.h"
#include "psi4/optking/physconst.h"
#include "linear_algebra.h"
#include "psi4/libparallel/ParallelPrinter.h"
#include "print.h"
#define EXTERN
#include "globals.h"

namespace opt {

using namespace v3d;
using namespace std;

// constructor - makes sure A<C
BEND::BEND(int A_in, int B_in, int C_in, bool freeze_in) : SIMPLE_COORDINATE(bend_type, 3, freeze_in) {
  //oprintf("constructing BEND A_in:%d B_in:%d C_in:%d, frozen %d\n", A_in, B_in, C_in, freeze_in);
  _bend_type = 0;
  axes_fixed = false;

  if (A_in == B_in || B_in == C_in || A_in == C_in)
    throw(INTCO_EXCEPT("BEND::BEND() Atoms defining bend are not unique.", true)); // bad error, abort

  s_atom[1] = B_in;
  if (A_in <= C_in) {
    s_atom[0] = A_in;
    s_atom[2] = C_in;
  }
  else {
    s_atom[0] = C_in;
    s_atom[2] = A_in;
  }

  x[0] = x[1] = x[2] = 0.0;
  w[0] = w[1] = w[2] = 0.0;
}

void BEND::compute_axes(GeomType geom) const {
  double u[3], v[3], w2[3], tv1[3], tv2[3];
  //tv1[0] =  1; tv1[1] = -1; tv1[2] = 1; // arbitrary search vectors
  //tv2[0] = -1; tv2[1] =  1; tv2[2] = 1;
  tv1[0] = 1; tv1[1] = 0; tv1[2] = 0; // more likely not to create 2 bends
  tv2[0] = 0; tv2[1] = 0; tv2[2] = 1; // that both break a symmetry plane
  v3d_normalize(tv1);
  v3d_normalize(tv2);

  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[0]], u); // B->A
  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[2]], v); // B->C
  v3d_normalize(u); // eBA
  v3d_normalize(v); // eBC

  if (_bend_type == 0) {        // not a linear-bend type
    v3d_cross_product(u, v, w); // orthogonal vector
    v3d_normalize(w);
    v3d_axpy(1, u, v, x);       // angle bisector
    v3d_normalize(x);
  }
  // l-b type, but u and v are not parallel (angle is not ~180)
  else if (!v3d_is_parallel(u,v)) {
    // _bend_type == 1 case
    v3d_cross_product(u, v, w); // orthogonal vector
    v3d_normalize(w);
    v3d_axpy(1.0, u, v, x);     // angle bisector
    v3d_normalize(x);
    if (_bend_type == 2) {
      array_copy(w, w2, 3);  // x_normal -> w_complement
      array_copy(x, w, 3);  // -w_normal -> x_complement
      v3d_scm(-1.0, w);
      array_copy(w2, x, 3);
    }
  }
  // l-b type, u and v are parallel but not to tv1.
  else if (!v3d_is_parallel(u,tv1) && !v3d_is_parallel(v,tv1)) {
    v3d_cross_product(u, tv1, w);
    v3d_normalize(w);
    v3d_cross_product(w,   u, x);
    v3d_normalize(x);
    if (_bend_type == 2) {
      array_copy(w, w2, 3);  // x_normal -> w_complement
      array_copy(x, w, 3);  // -w_normal -> x_complement
      v3d_scm(-1.0, w);
      array_copy(w2, x, 3);
    }
  }
  // l-b type, u and v are parallel but not to tv2.
  else if (!v3d_is_parallel(u,tv2) && !v3d_is_parallel(v,tv2)) {
    v3d_cross_product(u, tv2, w);
    v3d_normalize(w);
    v3d_cross_product(w,   u, x);
    v3d_normalize(x);
    if (_bend_type == 2) {
      array_copy(w, w2, 3);  // x_normal -> w_complement
      array_copy(x, w, 3);  // -w_normal -> x_complement
      v3d_scm(-1.0, w);
      array_copy(w2, x, 3);
    }
  }
  return;
}

// compute angle and store value in radians
double BEND::value(GeomType geom) const {
  double phi=0.0, tval=0.0;

  if (!axes_fixed)
    compute_axes(geom);

  double u[3], v[3];
  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[0]], u); // B->A
  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[2]], v); // B->C
  v3d_normalize(u); // eBA
  v3d_normalize(v); // eBC

/*
  if (_bend_type == 0) { // normal bend; could use axes but don't need to
    if ( ! v3d_angle(geom[s_atom[0]], geom[s_atom[1]], geom[s_atom[2]], phi) )
      throw(INTCO_EXCEPT("BEND::compute_val: could not compute angle"));
  }
  else {
*/
    // linear bend is sum of 2 angles, u.x + v.x
    double *origin = init_array(3);
    if (!v3d_angle(u, origin, x, phi))
      throw(INTCO_EXCEPT("BEND::value: could not compute linear bend",true));

    if (!v3d_angle(x, origin, v, tval))
      throw(INTCO_EXCEPT("BEND::value: could not compute linear bend",true));
    phi += tval;
    free_array(origin);
 // }
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
  double u[3], v[3];
  double **dqdx = init_matrix(3,3);

  if (!axes_fixed)
    compute_axes(geom);

  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[0]], u); // B->A
  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[2]], v); // B->C
  double Lu = v3d_norm(u); // RBA
  double Lv = v3d_norm(v); // RBC
  v3d_scm(1.0/Lu, u); // eBA
  v3d_scm(1.0/Lv, v); // eBC

  /*if (_bend_type == 2) {
    printf("bend_type == 2\n");
    printf("w %15.10f %15.10f %15.10f\n", w[0], w[1], w[2]);
    printf("x %15.10f %15.10f %15.10f\n", x[0], x[1], x[2]);
  }*/

  double uXw[3],wXv[3];
  v3d_cross_product(u, w, uXw);
  v3d_cross_product(w, v, wXv);

  for (int a=0; a<3; ++a)
    for (int i=0; i<3; ++i)
      dqdx[a][i] = zeta(a,0,1)*uXw[i]/Lu + zeta(a,2,1)*wXv[i]/Lv;

  return dqdx;
}

double ** BEND::Dq2Dx2(GeomType geom) const {
  double u[3], v[3];
  double **dq2dx2 = init_matrix(9,9);

  if (!axes_fixed)
    compute_axes(geom);

  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[0]], u); // B->A
  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[2]], v); // B->C
  double Lu = v3d_norm(u); // RBA
  double Lv = v3d_norm(v); // RBC
  v3d_scm(1.0/Lu, u);  // eBA
  v3d_scm(1.0/Lv, v);  // eBC

  // compute first derivatives
  double uXw[3],wXv[3];
  v3d_cross_product(u, w, uXw);
  v3d_cross_product(w, v, wXv);

  double **dqdx = init_matrix(3,3);
  for (int a=0; a<3; ++a)
    for (int i=0; i<3; ++i)
      dqdx[a][i] = zeta(a,0,1)*uXw[i]/Lu + zeta(a,2,1)*wXv[i]/Lv;

  double q = value(geom);
  double cos_q = cos(q); // cos_q = v3d_dot(u,v);

  if (1.0-cos_q*cos_q <= 1.0e-12) return dq2dx2; // leave 2nd derivatives empty - sin 0 = 0 in denominator
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


void BEND::print(std::string psi_fp, FILE *qc_fp, GeomType geom, int off) const {
  ostringstream iss(ostringstream::out); // create stream; allow output to it
  iss << get_definition_string(off);

  double val = value(geom);
  if (!s_frozen)
    oprintf(psi_fp, qc_fp, "\t %-15s  =  %15.6lf\t%15.6lf\n", iss.str().c_str(), val, val/_pi*180.0);
  else
    oprintf(psi_fp, qc_fp, "\t*%-15s  =  %15.6lf\t%15.6lf\n", iss.str().c_str(), val, val/_pi*180.0);
}

// function to return string of coordinate definition
std::string BEND::get_definition_string(int off) const {
  ostringstream iss(ostringstream::out); // create stream; allow output to it

  if (_bend_type == 0)
    iss << "B(";
  else if (_bend_type == 1)
    iss << "L(";
  else
    iss << "l(";

  iss << s_atom[0]+1+off << "," << s_atom[1]+1+off << "," << s_atom[2]+1+off << ")" << flush ;
  return iss.str();
}

void BEND::print_disp(std::string psi_fp, FILE *qc_fp, const double q_old, const double f_q,
    const double dq, const double q_new, int atom_offset) const {
  ostringstream iss(ostringstream::out); // create stream; allow output to it
  if (s_frozen) iss << "*";

  if (_bend_type == 0)
    iss << "B(";
  else if (_bend_type == 1)
    iss << "L(";
  else
    iss << "l(";

  iss << s_atom[0]+atom_offset+1 << "," << s_atom[1]+atom_offset+1 << "," << s_atom[2]+atom_offset+1 << ")" << flush ;

  oprintf(psi_fp, qc_fp, "%-15s = %13.6lf%13.6lf%13.6lf%13.6lf\n",
    iss.str().c_str(), q_old/_pi*180.0, f_q*_hartree2aJ*_pi/180.0,
      dq/_pi*180.0, q_new/_pi*180.0);
}

void BEND::print_intco_dat(std::string psi_fp, FILE *qc_fp, int off) const {
  if (_bend_type == 0) {
    if (s_frozen)
      oprintf(psi_fp, qc_fp,  "B*%6d%6d%6d", s_atom[0]+1+off, s_atom[1]+1+off, s_atom[2]+1+off);
    else
      oprintf(psi_fp, qc_fp,  "B %6d%6d%6d", s_atom[0]+1+off, s_atom[1]+1+off, s_atom[2]+1+off);
  }
  else if (_bend_type == 1) {
    if (s_frozen)
      oprintf(psi_fp, qc_fp,  "L*%6d%6d%6d", s_atom[0]+1+off, s_atom[1]+1+off, s_atom[2]+1+off);
    else
      oprintf(psi_fp, qc_fp,  "L %6d%6d%6d", s_atom[0]+1+off, s_atom[1]+1+off, s_atom[2]+1+off);
  }
  else {
    if (s_frozen)
      oprintf(psi_fp, qc_fp,  "l*%6d%6d%6d", s_atom[0]+1+off, s_atom[1]+1+off, s_atom[2]+1+off);
    else
      oprintf(psi_fp, qc_fp,  "l %6d%6d%6d", s_atom[0]+1+off, s_atom[1]+1+off, s_atom[2]+1+off);
  }
  if (s_has_fixed_eq_val)
    oprintf(psi_fp, qc_fp,  "%10.5lf", s_fixed_eq_val);
  oprintf(psi_fp, qc_fp,  "\n");
}

void BEND::print_s(std::string psi_fp, FILE *qc_fp, GeomType geom) const {
  if (_bend_type == 0)
    oprintf(psi_fp, qc_fp, "S vector for bend, B(%d %d %d): \n", s_atom[0]+1, s_atom[1]+1, s_atom[2]+1);
  else if (_bend_type == 1)
    oprintf(psi_fp, qc_fp, "S vector for bend, L(%d %d %d): \n", s_atom[0]+1, s_atom[1]+1, s_atom[2]+1);
  else
    oprintf(psi_fp, qc_fp, "S vector for bend, l(%d %d %d): \n", s_atom[0]+1, s_atom[1]+1, s_atom[2]+1);

  double **dqdx = DqDx(geom);
  oprintf(psi_fp, qc_fp,  "Atom 1: %12.8f %12.8f,%12.8f\n", dqdx[0][0], dqdx[0][1], dqdx[0][2]);
  oprintf(psi_fp, qc_fp,  "Atom 2: %12.8f %12.8f,%12.8f\n", dqdx[1][0], dqdx[1][1], dqdx[1][2]);
  oprintf(psi_fp, qc_fp,  "Atom 3: %12.8f %12.8f,%12.8f\n", dqdx[2][0], dqdx[2][1], dqdx[2][2]);
  free_matrix(dqdx);
}

bool BEND::operator==(const SIMPLE_COORDINATE & s2) const {
  if (bend_type != s2.g_type()) // warning: bend_type here is an INTCO_TYPE
    return false;

  if (this->s_atom[1] != s2.g_atom(1))
    return false;

  if (this->_bend_type != s2.g_bend_type())
    return false;

  if (this->s_atom[0] == s2.g_atom(0) && this->s_atom[2] == s2.g_atom(2))
    return true;
  else if (this->s_atom[0] == s2.g_atom(2) && this->s_atom[2] == s2.g_atom(0))
    return true;
  else
    return false;
}

}
