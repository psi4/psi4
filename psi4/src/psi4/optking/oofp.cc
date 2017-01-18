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

/*! \file oofp.cc
    \ingroup optking
    \brief oofp class
*/

#include "oofp.h"

#include <sstream>

#include "v3d.h"
#include "psi4/optking/physconst.h"
#include "opt_params.h"
#include "psi4/libparallel/ParallelPrinter.h"
#include "print.h"
#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
#endif

namespace opt {

extern OPT_PARAMS Opt_params;
using namespace v3d;
using std::ostringstream;

// constructor.  Definition (A,B,C,D) means angle AB with the CBD plane; canonical order is C < D
OOFP::OOFP(int A_in, int B_in, int C_in, int D_in, bool freeze_in) : SIMPLE_COORDINATE(oofp_type, 4, freeze_in) {

  if ( A_in==B_in || A_in==C_in || A_in==D_in || B_in==C_in || B_in==D_in || C_in==D_in)
    throw(INTCO_EXCEPT((char *)"OOFP::OOFP() Atoms defining oofp are not unique."));

  s_atom[0] = A_in;
  s_atom[1] = B_in;

  if (C_in < D_in) {
    s_atom[2] = C_in;
    s_atom[3] = D_in;
  }
  else {
    s_atom[2] = D_in;
    s_atom[3] = C_in;
  }
  near_180 = 0;
}

void OOFP::fix_oofp_near_180(GeomType geom) {
  double tval = value(geom);
  if ( tval > Opt_params.fix_tors_near_pi)
    near_180 = +1;
  else if ( tval < -1*Opt_params.fix_tors_near_pi)
    near_180 = -1;
  else
    near_180 = 0;
  return;
}

// compute angle and return value in radians
double OOFP::value(GeomType geom) const {
  double tau;

  if (! v3d_oofp(geom[s_atom[0]], geom[s_atom[1]], geom[s_atom[2]], geom[s_atom[3]], tau) )
    throw(INTCO_EXCEPT((char *)"OOFP::compute_val: unable to compute out-of-plane value",true));

  // Extend domain of out-of-plane angles values by checking past.
  if (near_180 == -1 && tau > Opt_params.fix_tors_near_pi)
    return (tau - 2.0 * _pi);
  else if (near_180 == +1 && tau < -1*Opt_params.fix_tors_near_pi)
    return (tau + 2.0 * _pi);
  else
    return tau;

/* Here is old psi3 code
      for (j=0;j<3;++j) {
        eBA[j] = geom[3*A+j] - geom[3*B+j];
        eBC[j] = geom[3*C+j] - geom[3*B+j];
        eBD[j] = geom[3*D+j] - geom[3*B+j];
      }

      rBA = sqrt( SQR(eBA[0]) + SQR(eBA[1]) + SQR(eBA[2]) );
      rBC = sqrt( SQR(eBC[0]) + SQR(eBC[1]) + SQR(eBC[2]) );
      rBD = sqrt( SQR(eBD[0]) + SQR(eBD[1]) + SQR(eBD[2]) );

      scalar_div(rBA,eBA);
      scalar_div(rBC,eBC);
      scalar_div(rBD,eBD);

      dot_array(eBC,eBD,3,&phi_CBD);

      if (phi_CBD > 1.0) phi_CBD = 0.0;
      else if (phi_CBD < -1.0) phi_CBD = _pi ;
      else phi_CBD = acos(phi_CBD) ;

      cross_product(eBC,eBD,tmp);

      dot_array(tmp,eBA,3,&dotprod);

      if (sin(phi_CBD) > optinfo.sin_phi_denominator_tol) dotprod = dotprod / sin(phi_CBD) ;
      else dotprod = 0.0 ;

      if (dotprod > 1.0) angle = _pi / 2.0;
      else if (dotprod < -1.0) angle = -1.0 * _pi / 2.0000;
      else angle = asin(dotprod) ;
*/
}

/*
Maybe someday I'll need to know if a out-of-plane is being intepreted as having
past through 180, but for now it doesn't seem I do.
// returns 'true' if the out-of-plane value is being corrected due to being past 180 or -180
// relative to when the out-of-plane values were fixed.
bool OOFP::fix_oofp_value_corrected(GeomType geom) const {
  double tau;
  v3d_oofp(geom[s_atom[0]], geom[s_atom[1]], geom[s_atom[2]], geom[s_atom[3]], tau);

  if (near_180 == -1 && tau > Opt_params.fix_tors_near_pi)
    return true;
  else if (near_180 == +1 && tau < -1*Opt_params.fix_tors_near_pi)
    return true;
  else
    return false;
}
*/

inline int zeta(const int a, const int m, const int n) {
  if (a == m) return 1;
  else if (a == n) return -1;
  else return 0;
}

inline int delta(const int i, const int j) {
  if (i == j) return 1;
  else return 0;
}


// out-of-plane is m-o-p-n
double ** OOFP::DqDx(GeomType geom) const {
  double **dqdx = init_matrix(4,3);

  // Get coordinates
  double A[3], B[3], C[3], D[3];
  for (int j=0; j<3; ++j) {
    A[j] = geom[s_atom[0]][j];
    B[j] = geom[s_atom[1]][j];
    C[j] = geom[s_atom[2]][j];
    D[j] = geom[s_atom[3]][j];
  }

  // compute E vectors
  double eBA[3], eBC[3], eBD[3];
  v3d_eAB(B,A,eBA);
  v3d_eAB(B,C,eBC);
  v3d_eAB(B,D,eBD);

  // compute value, C-B-D angle, and AB distance
  double val = value(geom);
  double phi_CBD;
  v3d_angle(C, B, D, phi_CBD);
  double rBA = v3d_dist(B, A);
  double rBC = v3d_dist(B, C);
  double rBD = v3d_dist(B, D);

  // initialize working arrays
  double *tmp = init_array(3);
  double *tmp2 = init_array(3);
  double *tmp3 = init_array(3);
  double *zero = init_array(3);

  // S vector for A

  v3d_cross_product(eBC, eBD, tmp);
  v3d_scm(1.0 / (cos(val) * sin(phi_CBD)), tmp);
  v3d_axpy(tan(val), eBA, zero, tmp2);
  for (int j=0;j<3;++j)
    dqdx[0][j] = (tmp[j] - tmp2[j])/rBA;

  // S vector for C

  v3d_cross_product(eBD, eBA, tmp);
  v3d_scm(1.0 / (cos(val) * sin(phi_CBD)),tmp);

  v3d_axpy(cos(phi_CBD), eBD, zero, tmp2);
  v3d_axpy( -1.0, tmp2, eBC, tmp3);

  v3d_scm( tan(val) / (sin(phi_CBD) * sin(phi_CBD)), tmp3);

  for (int j=0;j<3;++j)
    dqdx[2][j] = (tmp[j] - tmp3[j])/rBC;

  // S vector for D

  v3d_cross_product(eBA, eBC, tmp);
  v3d_scm(1.0 / (cos(val) * sin(phi_CBD)), tmp);

  v3d_axpy(cos(phi_CBD), eBC, zero, tmp2);
  v3d_axpy(-1.0, tmp2, eBD, tmp3);
  v3d_scm(tan(val) / (sin(phi_CBD) * sin(phi_CBD)), tmp3);

  for (int j=0;j<3;++j)
    dqdx[3][j] = (tmp[j] - tmp3[j])/rBD;

  // S vector for D
  for (int j=0;j<3;++j)
    dqdx[1][j]  = -1.0 * dqdx[0][j] - dqdx[2][j] - dqdx[3][j];

/* here is old psi3 code
      for (j=0;j<3;++j) {
        eBA[j] = geom[3*A+j] - geom[3*B+j];
        eBC[j] = geom[3*C+j] - geom[3*B+j];
        eBD[j] = geom[3*D+j] - geom[3*B+j];
      }

      rBA = sqrt( SQR(eBA[0]) + SQR(eBA[1]) + SQR(eBA[2]) );
      rBC = sqrt( SQR(eBC[0]) + SQR(eBC[1]) + SQR(eBC[2]) );
      rBD = sqrt( SQR(eBD[0]) + SQR(eBD[1]) + SQR(eBD[2]) );

      scalar_div(rBA,eBA);
      scalar_div(rBC,eBC);
      scalar_div(rBD,eBD);

      dot_array(eBC,eBD,3,&phi_CBD);

      if (phi_CBD > 1.0) phi_CBD = 0.0;
      else if (phi_CBD < -1.0) phi_CBD = _pi;
      else phi_CBD = acos(phi_CBD);

      cross_product(eBC,eBD,tmp);
      scalar_div(cos(val_rad)*sin(phi_CBD),tmp);
      for (j=0;j<3;++j)
         tmp2[j] = tan(val_rad) * eBA[j];
      for (j=0;j<3;++j)
         tmp3[j] = (tmp[j] - tmp2[j])/rBA;
      s_A[0] = tmp3[0];
      s_A[1] = tmp3[1];
      s_A[2] = tmp3[2];

      cross_product(eBD,eBA,tmp);
      scalar_div(cos(val_rad)*sin(phi_CBD),tmp);
      for (j=0;j<3;++j)
        tmp2[j] = cos(phi_CBD) * eBD[j];
      for (j=0;j<3;++j)
        tmp3[j] = eBC[j] - tmp2[j];
      scalar_mult(tan(val_rad)/SQR(sin(phi_CBD)),tmp3,3);
      for (j=0;j<3;++j)
         tmp2[j] = (tmp[j] - tmp3[j])/rBC;
      s_C[0] = tmp2[0];
      s_C[1] = tmp2[1];
      s_C[2] = tmp2[2];
      cross_product(eBA,eBC,tmp);
      scalar_div(cos(val_rad)*sin(phi_CBD),tmp);
      for (j=0;j<3;++j)
        tmp2[j] = cos(phi_CBD) * eBC[j];
      for (j=0;j<3;++j)
        tmp3[j] = eBD[j] - tmp2[j];
      scalar_mult(tan(val_rad)/SQR(sin(phi_CBD)),tmp3,3);

      for (j=0;j<3;++j)
        s_D[j] = (tmp[j] - tmp3[j])/rBD;

      for (j=0;j<3;++j)
        s_B[j]  = (-1.0) * s_A[j] - s_C[j] - s_D[j];

*/
  return dqdx;
}

// out-of-plane is m-o-p-n
// There are several errors in JCP, 22, 9164, (2002)
// I identified incorrect signs by making the equations invariant to reversing the atom indices
// (0,1,2,3) -> (3,2,1,0) and checking terms against finite differences.  Also, the last terms
// with sin^2 in the denominator are incorrectly given as only sin^1 in the paper.
// -RAK 2010
double ** OOFP::Dq2Dx2(GeomType ) const {
  double **dq2dx2 = init_matrix(12,12);
/* not yet implemented
  double u[3], v[3], w[3];
  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[0]], u); // u=m-o eBA
  v3d_axpy(-1, geom[s_atom[2]], geom[s_atom[3]], v); // v=n-p eCD
  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[2]], w); // w=p-o eBC
  double Lu = v3d_norm(u); // RBA
  double Lv = v3d_norm(v); // RCD
  double Lw = v3d_norm(w); // RBC
  v3d_scm(1.0/Lu, u);  // eBA
  v3d_scm(1.0/Lv, v);  // eCD
  v3d_scm(1.0/Lw, w);  // eBC

  double cos_u = v3d_dot(u,w);
  double cos_v = -1.0 * v3d_dot(v,w);
  if (1.0 - cos_u*cos_u <= 1.0e-12) return dq2dx2; // abort and leave zero if 0 or 180 angle
  if (1.0 - cos_v*cos_v <= 1.0e-12) return dq2dx2; // abort and leave zero if 0 or 180 angle
  double sin_u = sqrt(1.0 - cos_u*cos_u);
  double sin_v = sqrt(1.0 - cos_v*cos_v);

  double uXw[3], vXw[3];
  v3d_cross_product(u, w, uXw);
  v3d_cross_product(v, w, vXw);

  double sinu4 = sin_u*sin_u*sin_u*sin_u;
  double sinv4 = sin_v*sin_v*sin_v*sin_v;
  double cosu3 = cos_u*cos_u*cos_u;
  double cosv3 = cos_v*cos_v*cos_v;

  double tval, tval1, tval1b, tval2, tval2b, tval3, tval3b, tval4, tval4b;

  int k; // cartesian ; not i or j
  for (int a=0; a<4; ++a) {
    for (int b=0; b<=a; ++b) {
      for (int i=0; i<3; ++i) { //i = a_xyz
        for (int j=0; j<3; ++j) { //j=b_xyz
          tval = 0;

          if ((a==0 && b==0) || (a==1 && b==0) || (a==1 && b ==1))
            tval +=  zeta(a,0,1)*zeta(b,0,1)*(uXw[i]*(w[j]*cos_u-u[j]) + uXw[j]*(w[i]*cos_u-u[i])) / (Lu*Lu*sinu4);

          // above under reversal of atom indices, u->v ; w->(-w) ; uXw->(-uXw)
          if ((a==3 && b==3) || (a==3 && b==2) || (a==2 && b==2))
            tval += zeta(a,3,2)*zeta(b,3,2)*(vXw[i]*(w[j]*cos_v+v[j]) + vXw[j]*(w[i]*cos_v+v[i])) / (Lv*Lv*sinv4);

          if ((a==1 && b==1) || (a==2 && b==1) || (a==2 && b==0) || (a==1 && b==0))
            tval +=   (zeta(a,0,1)*zeta(b,1,2)+zeta(a,2,1)*zeta(b,1,0)) *
              (uXw[i] * (w[j] - 2*u[j]*cos_u + w[j]*cos_u*cos_u) +
               uXw[j] * (w[i] - 2*u[i]*cos_u + w[i]*cos_u*cos_u)) / (2*Lu*Lw*sinu4);

          if ((a==3 && b==2) || (a==3 && b==1) || (a==2 && b==2) || (a==2 && b==1))
            tval +=  (zeta(a,3,2)*zeta(b,2,1)+zeta(a,1,2)*zeta(b,2,3)) *
              (vXw[i] * (w[j] + 2*v[j]*cos_v + w[j]*cos_v*cos_v) +
               vXw[j] * (w[i] + 2*v[i]*cos_v + w[i]*cos_v*cos_v)) / (2*Lv*Lw*sinv4);

          if ((a==1 && b==1) || (a==2 && b==2) || (a==2 && b==1))
            tval +=   zeta(a,1,2)*zeta(b,2,1)*
              (uXw[i]*(u[j] + u[j]*cos_u*cos_u - 3*w[j]*cos_u + w[j]*cosu3) +
               uXw[j]*(u[i] + u[i]*cos_u*cos_u - 3*w[i]*cos_u + w[i]*cosu3)) / (2*Lw*Lw*sinu4);

          if ((a==2 && b==1) || (a==2 && b==2) || (a==1 && b==1))
            tval +=  zeta(a,2,1)*zeta(b,1,2)*
              (vXw[i]*(-v[j] - v[j]*cos_v*cos_v - 3*w[j]*cos_v + w[j]*cosv3) +
               vXw[j]*(-v[i] - v[i]*cos_v*cos_v - 3*w[i]*cos_v + w[i]*cosv3)) / (2*Lw*Lw*sinv4);

          if ((a != b) && (i != j)) {
            if ( i!=0 && j!=0 ) k = 0; // k is unique coordinate not i or j
            else if ( i!=1 && j!=1 ) k = 1;
            else k = 2;

            if (a==1 && b==1)
              tval +=  zeta(a,0,1)*zeta(b,1,2) * (j-i) *
                pow(-0.5, fabs(j-i)) * (+w[k]*cos_u - u[k]) / (Lu*Lw*sin_u*sin_u);

            if ((a==3 && b==2) || (a==3 && b==1) || (a==2 && b==2) || (a==2 && b==1))
              tval +=  zeta(a,3,2)*zeta(b,2,1) * (j-i) *
                pow(-0.5, fabs(j-i)) * (-w[k]*cos_v - v[k]) / (Lv*Lw*sin_v*sin_v);

            if ((a==2 && b==1) || (a==2 && b==0) || (a==1 && b==1) || (a==1 && b==0))
              tval +=  zeta(a,2,1)*zeta(b,1,0) * (j-i) *
                pow(-0.5, fabs(j-i)) * (-w[k]*cos_u + u[k]) / (Lu*Lw*sin_u*sin_u);

            if (a==2 && b==2)
              tval +=  zeta(a,1,2)*zeta(b,2,3) * (j-i) *
                pow(-0.5, fabs(j-i)) * (+w[k]*cos_v + v[k]) / (Lv*Lw*sin_v*sin_v);
          }
          dq2dx2[3*a+i][3*b+j] = dq2dx2[3*b+j][3*a+i] = tval;
        } // j
      } // i
    } // atom b
  } // atom a
*/
  return dq2dx2;
}


void OOFP::print(std::string psi_fp, FILE *qc_fp, GeomType geom, int off) const {
  ostringstream iss(ostringstream::out); // create stream; allow output to it
  iss << get_definition_string(off);
  double val = value(geom);
  if (!s_frozen)
    oprintf(psi_fp, qc_fp, "\t %-15s  =  %15.6lf\t%15.6lf\n", iss.str().c_str(), val, val/_pi*180.0);
  else
    oprintf(psi_fp, qc_fp, "\t*%-15s  =  %15.6lf\t%15.6lf\n", iss.str().c_str(), val, val/_pi*180.0);

}

// function to return string of coordinate definition
std::string OOFP::get_definition_string(int off) const {
  ostringstream iss(ostringstream::out); // create stream; allow output to it
  iss << "D(" << s_atom[0]+1+off << "," << s_atom[1]+1+off << "," << s_atom[2]+1+off << ","
    << s_atom[3]+1+off << ")" << std::flush;
  return iss.str();
}

void OOFP::print_disp(std::string psi_fp, FILE *qc_fp, const double q_old, const double f_q,
    const double dq, const double q_new, int atom_offset) const {
  ostringstream iss(ostringstream::out); // create stream; allow output to it
  if (s_frozen) iss << "*";
  iss << "D(" << s_atom[0]+atom_offset+1 << "," << s_atom[1]+atom_offset+1 << ","
    << s_atom[2]+atom_offset+1 << "," << s_atom[3]+atom_offset+1 << ")" << std::flush;
  oprintf(psi_fp, qc_fp, "%-15s = %13.6lf%13.6lf%13.6lf%13.6lf\n",
    iss.str().c_str(), q_old/_pi*180.0, f_q*_pi/180.0,dq/_pi*180.0, q_new/_pi*180.0);

}

void OOFP::print_intco_dat(std::string psi_fp, FILE *qc_fp, int off) const {
  if (s_frozen)
    oprintf(psi_fp, qc_fp, "O*%6d%6d%6d%6d", s_atom[0]+1+off, s_atom[1]+1+off,
      s_atom[2]+1+off, s_atom[3]+1+off);
  else
    oprintf(psi_fp, qc_fp, "O %6d%6d%6d%6d", s_atom[0]+1+off, s_atom[1]+1+off,
      s_atom[2]+1+off, s_atom[3]+1+off);

  if (s_has_fixed_eq_val)
    oprintf(psi_fp, qc_fp, "%10.5lf", s_fixed_eq_val);
  oprintf(psi_fp, qc_fp, "\n");
}

void OOFP::print_s(std::string psi_fp, FILE *qc_fp, GeomType geom) const {
  oprintf(psi_fp, qc_fp, "S vector for oofp, D(%d %d %d %d): \n",
    s_atom[0]+1, s_atom[1]+1, s_atom[2]+1, s_atom[3]+1);
  double **dqdx = DqDx(geom);
  oprintf(psi_fp, qc_fp, "Atom 1: %12.8f %12.8f,%12.8f\n", dqdx[0][0], dqdx[0][1], dqdx[0][2]);
  oprintf(psi_fp, qc_fp, "Atom 2: %12.8f %12.8f,%12.8f\n", dqdx[1][0], dqdx[1][1], dqdx[1][2]);
  oprintf(psi_fp, qc_fp, "Atom 3: %12.8f %12.8f,%12.8f\n", dqdx[2][0], dqdx[2][1], dqdx[2][2]);
  oprintf(psi_fp, qc_fp, "Atom 4: %12.8f %12.8f,%12.8f\n", dqdx[3][0], dqdx[3][1], dqdx[3][2]);
  free_matrix(dqdx);
}

bool OOFP::operator==(const SIMPLE_COORDINATE & s2) const {
  if (oofp_type != s2.g_type())
    return false;

  if (this->s_atom[0] == s2.g_atom(0) && this->s_atom[1] == s2.g_atom(1) &&
      this->s_atom[2] == s2.g_atom(2) && this->s_atom[3] == s2.g_atom(3) )
        return true;
  else if (this->s_atom[0] == s2.g_atom(3) && this->s_atom[1] == s2.g_atom(2) &&
           this->s_atom[2] == s2.g_atom(1) && this->s_atom[3] == s2.g_atom(0) )
        return true;
  else
        return false;
}

/*bool OOFP::check_oofp_for_bad_angles(GeomType geom) const {
  double phi_1=0, phi_2=0;

  v3d_angle(geom[s_atom[0]], geom[s_atom[1]], geom[s_atom[2]], phi_1);
  v3d_angle(geom[s_atom[1]], geom[s_atom[2]], geom[s_atom[3]], phi_2);

  double lim = Opt_params.oofp_angle_lim;

  if (phi_1 < lim || phi_1 > (_pi - lim))
    return false;

  if (phi_2 < lim || phi_2 > (_pi - lim))
    return false;

  return true;
}*/


}
