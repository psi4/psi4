/*! \file bend.cc
    \ingroup OPTKING
    \brief member functions for Bend class
*/

#include "bend.h"

namespace psi { namespace optking {

// constructor - makes sure A<C
Bend::Bend(int id_in, int A_in, int B_in, int C_in) {
  id = id_in;
  if (A_in <= C_in) {
    A = A_in;
    B = B_in;
    C = C_in;
  }
  else {
    C = A_in;
    B = B_in;
    A = C_in;
  }
  sA[0] = sA[1] = sA[2] = 0;
  sB[0] = sB[1] = sB[2] = 0;
  sC[0] = sC[1] = sC[2] = 0;
  val = 0;
}

void Bend::print (int print_flag) const {
  if (print_flag == 0)
    fprintf(outfile,"    (%d %d %d %d)\n", id, A+1, B+1, C+1);
  else if (print_flag == 1)
    fprintf(outfile,"   A%d: (%d %d %d) =  %.8lf \n", id,A+1,B+1,C+1,val);
}

void Bend::print_s(void) const {
  fprintf(outfile,"S vector for bend %d, A(%d %d %d): \n", id, A+1, B+1, C+1);
  fprintf(outfile,"Atom A: %12.8f %12.8f,%12.8f\n", sA[0],sA[1],sA[2]);
  fprintf(outfile,"Atom B: %12.8f %12.8f,%12.8f\n", sB[0],sB[1],sB[2]);
  fprintf(outfile,"Atom C: %12.8f %12.8f,%12.8f\n", sC[0],sC[1],sC[2]);
}

// sets values of bend - makes A < C
void Bend::set(int id_in, int A_in, int B_in, int C_in) {
  id = id_in;
  if (A_in <= C_in) {
    A = A_in;
    B = B_in;
    C = C_in;
  }
  else {
    C = A_in;
    B = B_in;
    A = C_in;
  }
}

// compute value
void Bend::compute(double **geom) {
  val = v3d_angle(geom[A],geom[B],geom[C]);
}

// compute s vectors
void Bend::compute_s(double **geom) {
  int xyz;
  double val_rad, rBA, rBC;
  double eBA[3], eBC[3];

  val_rad = val*_pi/180.0;

  for (xyz=0;xyz<3;++xyz) {
    eBA[xyz] = geom[A][xyz] - geom[B][xyz];
    eBC[xyz] = geom[C][xyz] - geom[B][xyz];
  }

  rBA = v3d_norm(eBA);
  rBC = v3d_norm(eBC);
  v3d_scm(1.0/rBA,eBA);
  v3d_scm(1.0/rBC,eBC);

  for (xyz=0;xyz<3;++xyz) {
    sA[xyz] = (eBA[xyz]*cos(val_rad) - eBC[xyz]) / (rBA*sin(val_rad));
    sB[xyz] = ((rBA - rBC*cos(val_rad))*eBA[xyz] + (rBC-rBA*cos(val_rad))*eBC[xyz])
                    / (rBA * rBC * sin(val_rad));
    sC[xyz] = (eBC[xyz]*cos(val_rad) - eBA[xyz]) / (rBC*sin(val_rad));
  }
  return;
}

// return component of s vector from coordinate index
double Bend::get_s(int atom_index, int xyz) const {
  if (atom_index == 0) return sA[xyz];
  else if (atom_index == 1) return sB[xyz];
  else if (atom_index == 2) return sC[xyz];
  else {
    fprintf(outfile,"invalid call to s_val\n");
    return 0.0;
  }
}

// return atom index of s vector from coordinate index
int Bend::get_atom(int atom_index) const {
  if (atom_index == 0) return A;
  else if (atom_index == 1) return B;
  else if (atom_index == 2) return C;
  else {
    fprintf(outfile,"invalid call to atom\n");
    return 0;
  }
}

bool Bend::operator==(const Bend & s2) const {
  if ( this->A == s2.A && this->B == s2.B && this->C == s2.C )
    return true;
  else
    return false;
}

}}

