/*! \file tors.cc
    \ingroup OPT10
    \brief tors class
*/

#include "tors.h"

#include <string>
#include <iostream>
#include <sstream>

#include "print.h"
#include "v3d.h"
#include "physconst.h"
#include "opt_params.h"

namespace opt {

extern OPT_PARAMS Opt_params;
using namespace v3d;
using std::ostringstream;

// constructor - canonical order is A < D
TORS::TORS(int A_in, int B_in, int C_in, int D_in) : SIMPLE(tors_type, 4) {
  //fprintf(stdout,"constructing TORS A_in:%d B_in:%d C_in:%d D_in: %d\n", A_in, B_in, C_in, D_in);
  if ( A_in==B_in || A_in==C_in || A_in==D_in || B_in==C_in || B_in==D_in || C_in==D_in)
    throw("TORS::TORS() Atoms defining tors are not unique.");

  if (A_in < D_in) {
    s_atom[0] = A_in;
    s_atom[1] = B_in;
    s_atom[2] = C_in;
    s_atom[3] = D_in;
  }
  else {
    s_atom[0] = D_in;
    s_atom[1] = C_in;
    s_atom[2] = B_in;
    s_atom[3] = A_in;
  }
  near_180 = 0;
}

void TORS::fix_near_180(void) {
  if ( s_val > Opt_params.fix_tors_near_pi)
    near_180 = +1;
  else if ( s_val < -1*Opt_params.fix_tors_near_pi)
    near_180 = -1;
  else
    near_180 = 0;
  return;
}

// compute angle and store value in radians
void TORS::compute_val(double **geom) {
  double tau;

  if (! v3d_tors(geom[s_atom[0]], geom[s_atom[1]], geom[s_atom[2]], geom[s_atom[3]], tau) )
    throw("TORS::compute_val: bond angles will not permit torsion computation");

  // Extend domain of torsion angles by checking past
  // extend domain of torsions so delta(vals) can be calculated
  if (near_180 == -1 && tau > Opt_params.fix_tors_near_pi)
    s_val = -_pi - (_pi - tau);
  else if (near_180 == +1 && tau < -1*Opt_params.fix_tors_near_pi)
    s_val = +_pi + (_pi + tau);
  else
    s_val = tau;
}

double * TORS::g_s(double **geom, const int iatom) const {
  double rAB,rBC,rCD;
  double eAB[3], eBC[3], eCD[3], tmp[3], tmp2[3];
  double phiABC, phiBCD;
  double *s = init_array(3);

  if ( iatom!=0 && iatom!=1 && iatom!=2 && iatom!=3 )
    throw("TORS::g_s: requested atom not 0, 1, 2, or 3");

  v3d_axpy(-1, geom[s_atom[0]], geom[s_atom[1]], eAB);
  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[2]], eBC);
  v3d_axpy(-1, geom[s_atom[2]], geom[s_atom[3]], eCD);

  rAB = v3d_norm(eAB);
  rBC = v3d_norm(eBC);
  rCD = v3d_norm(eCD);

  v3d_scm(1.0/rAB, eAB);
  v3d_scm(1.0/rBC, eBC);
  v3d_scm(1.0/rCD, eCD);

  v3d_angle(geom[s_atom[0]], geom[s_atom[1]], geom[s_atom[2]], phiABC);
  v3d_angle(geom[s_atom[1]], geom[s_atom[2]], geom[s_atom[3]], phiBCD);

  if (iatom == 0) {
    v3d_cross_product(eAB, eBC, tmp);
    v3d_scm(-1.0 / (rAB * sin(phiABC) * sin(phiABC)), tmp);
    s[0] = tmp[0];
    s[1] = tmp[1];
    s[2] = tmp[2];
  }
  else if (iatom == 1) {
    v3d_cross_product(eAB, eBC, tmp);
    v3d_scm( (rBC-rAB*cos(phiABC)) / (rBC*rAB*sin(phiABC)*sin(phiABC)), tmp);
    v3d_cross_product(eCD, eBC, tmp2);
    v3d_scm( cos(phiBCD) / (rBC*sin(phiBCD)*sin(phiBCD)), tmp2);
    s[0] = tmp[0] + tmp2[0];
    s[1] = tmp[1] + tmp2[1];
    s[2] = tmp[2] + tmp2[2];
  }
  else if (iatom == 2) {
    v3d_cross_product(eCD, eBC, tmp);
    v3d_scm( (rBC-rCD*cos(phiBCD)) / (rBC*rCD*sin(phiBCD)*sin(phiBCD)), tmp);
    v3d_cross_product(eAB, eBC, tmp2);
    v3d_scm( cos(phiABC) / (rBC*sin(phiABC)*sin(phiABC)), tmp2);
    s[0] = tmp[0] + tmp2[0];
    s[1] = tmp[1] + tmp2[1];
    s[2] = tmp[2] + tmp2[2];
  }
  else {
    v3d_cross_product(eCD, eBC, tmp);
    v3d_scm(-1.0 / (rCD*sin(phiBCD)*sin(phiBCD)), tmp);
    s[0] = tmp[0];
    s[1] = tmp[1];
    s[2] = tmp[2];
  }
  return s;
}

void TORS::print (const FILE *fp) const {
  ostringstream iss(ostringstream::out); // create stream; allow output to it
  iss << "T(" << s_atom[0]+1 << "," << s_atom[1]+1 << "," << s_atom[2]+1 << "," << s_atom[3]+1 << ")";
  fprintf(const_cast<FILE *>(fp),"\t %-15s  =  %15.6lf\t%15.6lf\n",
    iss.str().c_str(), s_val, s_val/_pi*180.0);
}

void TORS::print_disp(const FILE *fp, const double q_old, const double f_q,
    const double dq, const double q_new) const {
  ostringstream iss(ostringstream::out); // create stream; allow output to it
  iss << "T(" << s_atom[0]+1 << "," << s_atom[1]+1 << "," << s_atom[2]+1 << "," << s_atom[3]+1 << ")";
  fprintf(const_cast<FILE *>(fp),"\t %-15s = %13.6lf%13.6lf%13.6lf%13.6lf\n",
    iss.str().c_str(), q_old/_pi*180.0, f_q*_pi/180.0,dq/_pi*180.0, q_new/_pi*180.0);
}


void TORS::print_s(const FILE *fp, double **geom) const {
  fprintf(const_cast<FILE *>(fp),"S vector for tors, T(%d %d %d %d): \n",
    s_atom[0]+1, s_atom[1]+1, s_atom[2]+1, s_atom[3]+1);
  double *s0, *s1, *s2, *s3;
  s0 = g_s(geom,0);
  s1 = g_s(geom,1);
  s2 = g_s(geom,2);
  s3 = g_s(geom,3);
  fprintf(const_cast<FILE *>(fp), "Atom 1: %12.8f %12.8f,%12.8f\n", s0[0], s0[1], s0[2]);
  fprintf(const_cast<FILE *>(fp), "Atom 2: %12.8f %12.8f,%12.8f\n", s1[0], s1[1], s1[2]);
  fprintf(const_cast<FILE *>(fp), "Atom 3: %12.8f %12.8f,%12.8f\n", s2[0], s2[1], s2[2]);
  fprintf(const_cast<FILE *>(fp), "Atom 4: %12.8f %12.8f,%12.8f\n", s3[0], s3[1], s3[2]);
  free_array(s0);
  free_array(s1);
  free_array(s2);
}

bool TORS::operator==(const SIMPLE & s2) const {
  fprintf(stdout,"TORS::operator==\n");
  if (tors_type != s2.g_type())
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

}

