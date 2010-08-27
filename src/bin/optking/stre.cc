/*! \file stre.cc
    \ingroup OPT10
    \brief stretch class definition
*/

#include "stre.h"

#include <string>
//#include <iostream>
#include <sstream>

#include "v3d.h"
#include "physconst.h"

namespace opt {

using namespace v3d;
using std::ostringstream;

// constructor - makes sure A<B
STRE::STRE(int A_in, int B_in) : SIMPLE(stre_type, 2) {
  //fprintf(stdout,"constructing STRE A_in:%d B_in:%d\n", A_in, B_in);
  if (A_in == B_in) throw("STRE::STRE() atoms defining strech are not unique.");

  if (A_in < B_in) {
    s_atom[0] = A_in;
    s_atom[1] = B_in;
  }
  else {
    s_atom[0] = B_in;
    s_atom[1] = A_in;
  }
}

// compute value and store in au
void STRE::compute_val(double **geom) {
  s_val = v3d_dist(geom[s_atom[0]], geom[s_atom[1]]);
}

// compute s vector for atom and return
double * STRE::g_s(double **geom, const int iatom) const {
  double *eBA = init_array(3);

  if ( (iatom != 0) && (iatom != 1) )
    throw("STRE::g_s: requested atom not 0 or 1");

  if (! v3d_eAB(geom[s_atom[1]], geom[s_atom[0]], eBA) )
    throw("STRE::g_s: could not normalize s vector.");

  if (iatom == 1) {
      eBA[0] *= -1.0;
      eBA[1] *= -1.0;
      eBA[2] *= -1.0;
  }

  return eBA;
}

// print stretch and value
void STRE::print(const FILE *fp) const {
  ostringstream iss(ostringstream::out); // create stream; allow output to it
  iss << "R(" << s_atom[0]+1 << "," << s_atom[1]+1 << ")" ;
  fprintf(const_cast<FILE *>(fp),"\t %-15s  =  %15.6lf\t%15.6lf\n",
    iss.str().c_str(), s_val, s_val*_bohr2angstroms);
}

// print displacement
void STRE::print_disp(const FILE *fp, const double q_orig, const double f_q,
    const double dq, const double new_q) const {
  ostringstream iss(ostringstream::out);
  iss << "R(" << s_atom[0]+1 << "," << s_atom[1]+1 << ")" ;
  fprintf(const_cast<FILE *>(fp),"\t %-15s = %13.6lf%13.6lf%13.6lf%13.6lf\n",
    iss.str().c_str(), q_orig*_bohr2angstroms, f_q*_hartree2aJ/_bohr2angstroms,
    dq*_bohr2angstroms, new_q*_bohr2angstroms);
}


// print s vectors
void STRE::print_s(const FILE *fp, double **geom) const {
  fprintf(const_cast<FILE *>(fp),"S vector for stretch R(%d %d): \n",
    s_atom[0]+1, s_atom[1]+1);
  double *s0, *s1;
  s0 = g_s(geom, 0);
  s1 = g_s(geom, 1);
  fprintf(const_cast<FILE *>(fp),"Atom 1: %12.8f %12.8f,%12.8f\n", s0[0],s0[1],s0[2]);
  fprintf(const_cast<FILE *>(fp),"Atom 2: %12.8f %12.8f,%12.8f\n", s1[0],s1[1],s1[2]);
  free_array(s0);
  free_array(s1);
}

bool STRE::operator==(const SIMPLE & s2) const {
  fprintf(stdout,"STRE::==(SIMPLE &)\n");

  if (stre_type != s2.g_type())
    return false;

  if ((this->s_atom[0] == s2.g_atom(0)) && (this->s_atom[1] == s2.g_atom(1)))
    return true;
  else if ((this->s_atom[0] == s2.g_atom(1)) && (this->s_atom[1] == s2.g_atom(0)))
    return true;
  else
    return false;
}

}

