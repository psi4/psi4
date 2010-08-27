/*! \file bend.cc
    \ingroup OPT10
    \brief bend class
*/

#include "bend.h"

#include <sstream>
#include "v3d.h"
#include "physconst.h"

namespace opt {

using namespace v3d;
using namespace std;

// constructor - makes sure A<C
BEND::BEND(int A_in, int B_in, int C_in) : SIMPLE(bend_type, 3) {
  //fprintf(stdout,"constructing BEND A_in:%d B_in:%d C_in:%d\n", A_in, B_in, C_in);
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
void BEND::compute_val(double **geom) {
  double tval;
  if ( ! v3d_angle(geom[s_atom[0]], geom[s_atom[1]], geom[s_atom[2]], tval) )
    throw("BEND::compute_val: could not compute angle");
  s_val = tval;
}

double * BEND::g_s(double **geom, const int iatom) const {
  int xyz;
  bool rval;
  double ang, rBA, rBC;
  double eBA[3], eBC[3];
  double *s = init_array(3);

  if ( iatom!=0 && iatom!=1 && iatom!=2 )
    throw("BEND::g_s: requested atom not 0, 1, or 2");

  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[0]], eBA);
  v3d_axpy(-1, geom[s_atom[1]], geom[s_atom[2]], eBC);
  rBA = v3d_norm(eBA);
  rBC = v3d_norm(eBC);
  v3d_scm(1.0/rBA, eBA);
  v3d_scm(1.0/rBC, eBC);

  // compute value of angle in radians
  v3d_angle(geom[s_atom[0]], geom[s_atom[1]], geom[s_atom[2]], ang);

  if (iatom == 0) {
    for (xyz=0;xyz<3;++xyz)
      s[xyz] = (eBA[xyz]*cos(ang) - eBC[xyz]) / (rBA*sin(ang));
  }
  else if (iatom == 1) {
    for (xyz=0;xyz<3;++xyz)
      s[xyz] = ((rBA - rBC*cos(ang))*eBA[xyz] + (rBC-rBA*cos(ang))*eBC[xyz])
                      / (rBA * rBC * sin(ang));
  }
  else if (iatom == 2) {
    for (xyz=0;xyz<3;++xyz)
      s[xyz] = (eBC[xyz]*cos(ang) - eBA[xyz]) / (rBC*sin(ang));
  }
  return s;
}

void BEND::print (const FILE *fp) const {
  ostringstream iss(ostringstream::out); // create stream; allow output to it
  iss << "A(" << s_atom[0]+1 << "," << s_atom[1]+1 << "," << s_atom[2]+1 << ")" ;
  fprintf(const_cast<FILE *>(fp),"\t %-15s  =  %15.6lf\t%15.6lf\n",
    iss.str().c_str(), s_val, s_val/_pi*180.0);
}

void BEND::print_disp(const FILE *fp, const double q_old, const double f_q,
    const double dq, const double q_new) const {
  ostringstream iss(ostringstream::out); // create stream; allow output to it
  iss << "A(" << s_atom[0]+1 << "," << s_atom[1]+1 << "," << s_atom[2]+1 << ")" ;
  fprintf(const_cast<FILE *>(fp),"\t %-15s = %13.6lf%13.6lf%13.6lf%13.6lf\n",
    iss.str().c_str(), q_old/_pi*180.0, f_q*_hartree2aJ*_pi/180.0,
      dq/_pi*180.0, q_new/_pi*180.0);
}


void BEND::print_s(const FILE *fp, double **geom) const {
  fprintf(const_cast<FILE *>(fp),"S vector for bend, A(%d %d %d): \n",
    s_atom[0]+1, s_atom[1]+1, s_atom[2]+1);
  double *s0, *s1, *s2;
  s0 = g_s(geom,0);
  s1 = g_s(geom,1);
  s2 = g_s(geom,2);
  fprintf(const_cast<FILE *>(fp), "Atom 1: %12.8f %12.8f,%12.8f\n", s0[0], s0[1], s0[2]);
  fprintf(const_cast<FILE *>(fp), "Atom 2: %12.8f %12.8f,%12.8f\n", s1[0], s1[1], s1[2]);
  fprintf(const_cast<FILE *>(fp), "Atom 3: %12.8f %12.8f,%12.8f\n", s2[0], s2[1], s2[2]);
  free_array(s0);
  free_array(s1);
  free_array(s2);
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

