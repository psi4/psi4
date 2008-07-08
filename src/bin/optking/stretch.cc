/*! \file stretch.cc
    \ingroup OPTKING
    \brief Member functions for stretch class
*/

#include "stretch.h"

namespace psi { namespace optking {

// constructor - makes sure A<B
Stretch::Stretch(int id_in, int A_in, int B_in) {
  id = id_in;
  if (A_in < B_in) {
    A = A_in;
    B = B_in;
  }
  else {
    B = A_in;
    A = B_in;
  }
  sA[0] = sA[1] = sA[2] = 0;
  sB[0] = sB[1] = sB[2] = 0;
  val = 0;
}

void Stretch::print(int print_flag) const {
  if (print_flag == 0)
    fprintf(outfile,"    (%d %d %d)\n", id, A+1, B+1);
  else if (print_flag == 1)
    fprintf(outfile,"   R%d: (%d %d) =  %.8lf \n", id,A+1,B+1,val);
}

void Stretch::print_s(void) const {
  fprintf(outfile,"S vector for stretch %d, R(%d %d): \n", id, A+1, B+1);
  fprintf(outfile,"Atom A: %12.8f %12.8f,%12.8f\n", sA[0],sA[1],sA[2]);
  fprintf(outfile,"Atom B: %12.8f %12.8f,%12.8f\n", sB[0],sB[1],sB[2]);
}

// sets values of stretch - makes A < B
void Stretch::set(int id_in, int A_in, int B_in) {
  id = id_in;
  if (A_in < B_in) {
    A = A_in;
    B = B_in;
  }
  else {
    A = B_in;
    B = A_in;
  }
}

// compute value - given geom
void Stretch::compute(double **geom) {
  val = v3d_dist(geom[A],geom[B]);
}

// compute s vectors - given geom
void Stretch::compute_s(double **geom) {
  int xyz;
  double eBA[3];

  for (xyz=0;xyz<3;++xyz)
    eBA[xyz] = geom[A][xyz] - geom[B][xyz];

  v3d_normalize(eBA);

  for (xyz=0;xyz<3;++xyz) {
    sA[xyz] = eBA[xyz];
    sB[xyz] = -1.0 * eBA[xyz];
  }
}

// return component of s vector from coordinate index
double Stretch::get_s(int atom_index, int xyz) const {
  if (atom_index == 0) return sA[xyz];
  else if (atom_index == 1) return sB[xyz];
  else {
    fprintf(outfile,"invalid call to Stretch::s_val\n");
    return 0.0;
  }
}

// return atom (A or B) given coordinate index
int Stretch::get_atom(int atom_index) const {
  if (atom_index == 0) return A;
  else if (atom_index == 1) return B;
  else {
    fprintf(outfile,"invalid call to Stretch::atom\n");
    return 0;
  }
}

bool Stretch::operator==(const Stretch & s2) const {
  if ( this->A == s2.A && this->B == s2.B)
    return true;
  else
    return false;
}

}}

