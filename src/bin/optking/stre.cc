/*! \file stre.cc
    \ingroup optking
    \brief stretch class definition
*/

#include "stre.h"
#include <sstream>
#include "opt_except.h"

#include "v3d.h"
#include "physconst.h"

namespace opt {

using namespace v3d;
using std::ostringstream;

// constructor - makes sure A<B ; default freeze_in =false
STRE::STRE(int A_in, int B_in, bool freeze_in) : SIMPLE(stre_type, 2, freeze_in) {
  //fprintf(stdout,"constructing STRE A_in:%d B_in:%d, frozen %d\n",
  //  A_in, B_in, freeze_in);
  if (A_in == B_in) throw(INTCO_EXCEPT("STRE::STRE() atoms defining strech are not unique."));

  if (A_in < B_in) {
    s_atom[0] = A_in;
    s_atom[1] = B_in;
  }
  else {
    s_atom[0] = B_in;
    s_atom[1] = A_in;
  }
  hbond = false;
  inverse_stre = false;
}

// compute value and store in au
double STRE::value(GeomType geom) const {
  double tval = v3d_dist(geom[s_atom[0]], geom[s_atom[1]]);
  if (inverse_stre)
    return (1.0/tval);
  else
    return tval;
}

inline int delta(const int i, const int j) {
  if (i == j) return 1;
  else return 0;
}

// Equations for B and B' elements from
// Bakken & Helgaker, JCP, 117, 9160 (2002).

// compute and return array of first derivative (B marix elements)
// return order is:
// atom_i  cartesian
// 0       x y z
// 1       x y z ..
double ** STRE::DqDx(GeomType geom) const {
  double eAB[3], val;
  double **dqdx = init_matrix(2,3);

  if (! v3d_eAB(geom[s_atom[0]], geom[s_atom[1]], eAB) )
    throw(INTCO_EXCEPT("STRE::DqDx: could not normalize s vector", true));

  if (inverse_stre)
    val = value(geom); // val = (1/R)

  for (int a=0; a<2; ++a)
    for (int a_xyz=0; a_xyz<3; ++a_xyz) {
      dqdx[a][a_xyz] = eAB[a_xyz];
      if (a == 0) dqdx[a][a_xyz] *= -1;
      if (inverse_stre)
        dqdx[a][a_xyz] *= -1.0*val*val; // -(1/R)^2 * (dR/da)
    }

  return dqdx;
}

// compute and return array of second derivative (B' matrix elements)
// return order is:
// atom_a atom_b cartesians
// 0      1      xx xy xz yx yy yz zx zy zz
double ** STRE::Dq2Dx2(GeomType geom) const {
  double eAB[3];
  double **dq2dx2 = init_matrix(6,6);
  double tval;

  if (! v3d_eAB(geom[s_atom[0]], geom[s_atom[1]], eAB) )
    throw(INTCO_EXCEPT("STRE::Dq2Dx2: could not normalize s vector", true));

  if (!inverse_stre) {
    double length = value(geom);

    for (int a=0; a<2; ++a)
      for (int a_xyz=0; a_xyz<3; ++a_xyz)
        for (int b=0; b<2; ++b)
          for (int b_xyz=0; b_xyz<3; ++b_xyz) {
            tval = (eAB[a_xyz] * eAB[b_xyz] - delta(a_xyz,b_xyz))/length;
            if (a == b) tval *= -1.0;
            dq2dx2[3*a+a_xyz][3*b+b_xyz] = tval;
          }
  }
  else { // using 1/R
    double val = value(geom);

    double **dqdx = DqDx(geom); // matrix is (2,3)

    for (int a=0; a<2; ++a)
      for (int a_xyz=0; a_xyz<3; ++a_xyz)
        for (int b=0; b<2; ++b)
          for (int b_xyz=0; b_xyz<3; ++b_xyz)
            dq2dx2[3*a+a_xyz][3*b+b_xyz] = 2.0 / val * dqdx[a][a_xyz] * dqdx[b][b_xyz];

    free_matrix(dqdx);
  }
  return dq2dx2;
}

// print stretch and value
void STRE::print(FILE *fp, GeomType geom, int off) const {
  ostringstream iss(ostringstream::out); // create stream; allow output to it
  if (inverse_stre)
    iss << "1/R(" << s_atom[0]+1+off << "," << s_atom[1]+1+off << ")" << std::flush ;
  else
    iss << "R(" << s_atom[0]+1+off << "," << s_atom[1]+1+off << ")" << std::flush ;

  double val = value(geom);
  if (!s_frozen)
    fprintf(fp,"\t %-15s  =  %15.6lf\t%15.6lf\n", iss.str().c_str(), val, val*_bohr2angstroms);
  else
    fprintf(fp,"\t*%-15s  =  %15.6lf\t%15.6lf\n", iss.str().c_str(), val, val*_bohr2angstroms);
  fflush(fp);
}

// function to return string of coordinate definition
std::string STRE::get_definition_string(int off) const {
  ostringstream iss(ostringstream::out); // create stream; allow output to it
  if (inverse_stre)
    iss << "1/R(" << s_atom[0]+1+off << "," << s_atom[1]+1+off << ")" << std::flush ;
  else
    iss << "R(" << s_atom[0]+1+off << "," << s_atom[1]+1+off << ")" << std::flush ;
  return iss.str();
}

void STRE::print_intco_dat(FILE *fp, int off) const {
  if (hbond) {
    if (s_frozen)
      fprintf(fp, "H*%6d%6d\n", s_atom[0]+1+off, s_atom[1]+1+off);
    else
      fprintf(fp, "H %6d%6d\n", s_atom[0]+1+off, s_atom[1]+1+off);
  }
  else {
    if (s_frozen)
      fprintf(fp, "R*%6d%6d\n", s_atom[0]+1+off, s_atom[1]+1+off);
    else
      fprintf(fp, "R %6d%6d\n", s_atom[0]+1+off, s_atom[1]+1+off);
  }
  fflush(fp);
}

// print displacement
void STRE::print_disp(FILE *fp, const double q_orig, const double f_q,
    const double dq, const double new_q, int atom_offset) const {
  ostringstream iss(ostringstream::out);
  if (s_frozen) iss << "*";
  iss << "R(" << s_atom[0]+atom_offset+1 << "," << s_atom[1]+atom_offset+1 << ")" << std::flush ;
  fprintf(fp,"\t %-15s = %13.6lf%13.6lf%13.6lf%13.6lf\n",
    iss.str().c_str(), q_orig*_bohr2angstroms, f_q*_hartree2aJ/_bohr2angstroms,
    dq*_bohr2angstroms, new_q*_bohr2angstroms);
  fflush(fp);
}


// print s vectors
void STRE::print_s(FILE *fp, GeomType geom) const {
  fprintf(fp,"S vector for stretch R(%d %d): \n",
    s_atom[0]+1, s_atom[1]+1);
  double **dqdx = DqDx(geom);
  fprintf(fp,"Atom 1: %12.8f %12.8f,%12.8f\n", dqdx[0][0],dqdx[0][1],dqdx[0][2]);
  fprintf(fp,"Atom 2: %12.8f %12.8f,%12.8f\n", dqdx[1][0],dqdx[1][1],dqdx[1][2]);
  free_matrix(dqdx);
  fflush(fp);
}

bool STRE::operator==(const SIMPLE & s2) const {
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

