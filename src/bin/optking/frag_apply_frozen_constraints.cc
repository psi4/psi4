/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file frag_apply_frozen_constraints.cc
    \ingroup optking
    \brief apply constraints specified by given strings

     FRAG::apply_frozen_constraints(string) parses string for the atoms that define bonds, bends, and angles
     to be frozen.  Return false if there are none.

     FRAG::apply_fixed_constraints(string) parses string for the atoms that define bonds, bends, and angles
     to be fixed.  Return false if there are none.
*/

#include "frag.h"

#include <string>
#include <sstream>
#include <vector>

#include "print.h"
#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
#endif

namespace opt {

using std::string;

struct fixed_coord {
  std::vector<int> atoms;
  double eq_val;
};

std::vector<int> split_to_ints(string &s);
std::vector<fixed_coord> split_to_fixed_coord(string &str, int N);

// R_string = integer list of atoms, each 2 of which are frozen
// B_string = integer list of atoms, each 3 of which are frozen
// D_string = integer list of atoms, each 4 of which are frozen

bool FRAG::apply_frozen_constraints(string R_string, string B_string, string D_string)
{
  std::vector<int> R_atoms = split_to_ints(R_string);
  std::vector<int> B_atoms = split_to_ints(B_string);
  std::vector<int> D_atoms = split_to_ints(D_string);

  if (R_atoms.size() % 2 != 0)
    throw(INTCO_EXCEPT("Frozen distance string should contain even number of atoms."));
  if (B_atoms.size() % 3 != 0)
    throw(INTCO_EXCEPT("Frozen bend string should contain 3*(whole number) of atoms."));
  if (D_atoms.size() % 4 != 0)
    throw(INTCO_EXCEPT("Frozen dihedral string should contain 4*(whole number) of atoms."));

  if (R_atoms.size()) {
    oprintf_out("\tFrozen distance atom list: \n");
    for (int i=0; i<R_atoms.size(); i+=2)
      oprintf_out("\t %5d %5d\n", R_atoms[i]+1, R_atoms[i+1]+1);
  }
  if (B_atoms.size()) {
    oprintf_out("\tFrozen bend atom list: \n");
    for (int i=0; i<B_atoms.size(); i+=3)
      oprintf_out("\t %5d %5d %5d\n", B_atoms[i]+1, B_atoms[i+1]+1, B_atoms[i+2]+1);
  }
  if (D_atoms.size()) {
    oprintf_out("\tFrozen dihedral atom list: \n");
    for (int i=0; i<D_atoms.size(); i+=4)
      oprintf_out("\t %5d %5d %5d %5d\n", D_atoms[i]+1, D_atoms[i+1]+1,
        D_atoms[i+2]+1, D_atoms[i+3]+1);
  }

  if (!R_atoms.size() && !B_atoms.size() && !D_atoms.size())
    return false;

  // do frozen distances
  for (int i=0; i<R_atoms.size(); i+=2) {
    int a = R_atoms[i];
    int b = R_atoms[i+1];

    if (a >= natom || b >= natom)
      throw(INTCO_EXCEPT("Impossibly large index for atom in frozen distance string."));

    STRE *one_stre = new STRE(a, b, 1); // create frozen stretch between atoms a and b

    // check if coord is already present; returns 1 past the end if not found
    int index = find(one_stre);

    if (index == coords.simples.size())
      coords.simples.push_back(one_stre);// add it
    else { 
      coords.simples[index]->freeze();   // it's there already; make sure it's frozen
      delete one_stre;
    }
  }
  // do frozen bends
  for (int i=0; i<B_atoms.size(); i+=3) {
    int a = B_atoms[i];
    int b = B_atoms[i+1];
    int c = B_atoms[i+2];

    if (a >= natom || b >= natom || c >= natom)
      throw(INTCO_EXCEPT("Impossibly large index for atom in frozen bend string."));

    BEND *one_bend = new BEND(a, b, c, 1); // create frozen bend between atoms a,b,c

    // check if coord is already present; returns 1 past the end if not found
    int index = find(one_bend);

    if (index == coords.simples.size())
      coords.simples.push_back(one_bend);// add it
    else { 
      coords.simples[index]->freeze();   // it's there already; make sure it's frozen
      delete one_bend;
    }
  }
  // do dihedral angles
  for (int i=0; i<D_atoms.size(); i+=4) {
    int a = D_atoms[i];
    int b = D_atoms[i+1];
    int c = D_atoms[i+2];
    int d = D_atoms[i+3];

    if (a >= natom || b >= natom || c >= natom || d >= natom)
      throw(INTCO_EXCEPT("Impossibly large index for atom in frozen dihedral string."));

    TORS *one_tors = new TORS(a, b, c, d, 1); // create frozen dihedral between a,b,c,d

    // check if coord is already present; returns 1 past the end if not found
    int index = find(one_tors);

    if (index == coords.simples.size())
      coords.simples.push_back(one_tors);// add it
    else { 
      coords.simples[index]->freeze();   // it's there already; make sure it's frozen
      delete one_tors;
    }
  }
  return true;
}

bool FRAG::apply_fixed_constraints(string R_string, string B_string, string D_string)
{
  std::vector<fixed_coord> R = split_to_fixed_coord(R_string, 2);
  std::vector<fixed_coord> B = split_to_fixed_coord(B_string, 3);
  std::vector<fixed_coord> D = split_to_fixed_coord(D_string, 4);

  if (!R.size() && !B.size() && !D.size())
    return false;

  if (R.size()) {
    oprintf_out("\tFixed distance atom list: \n");
    for (int i=0; i<R.size(); ++i)
      oprintf_out("\t %5d %5d\n", R[i].atoms[0]+1, R[i].atoms[1]+1);
  }

  if (B.size()) {
    oprintf_out("\tFixed bend atom list: \n");
    for (int i=0; i<B.size(); ++i)
      oprintf_out("\t %5d %5d %5d\n", B[i].atoms[0]+1, B[i].atoms[1]+1, B[i].atoms[2]+1);
  }

  if (D.size()) {
    oprintf_out("\tFixed dihedral atom list: \n");
    for (int i=0; i<D.size(); ++i)
      oprintf_out("\t %5d %5d %5d %5d\n", D[i].atoms[0]+1, D[i].atoms[i+1]+1, D[i].atoms[2]+1, D[i].atoms[3]+1);
  }
  

  // do fixed distances
  for (int i=0; i<R.size(); ++i) {
    int a = R[i].atoms[0];
    int b = R[i].atoms[1];

    if (a >= natom || b >= natom)
      throw(INTCO_EXCEPT("Impossibly large index for atom in fixed distance string."));

    STRE *one_stre = new STRE(a, b, 0);
    // Insist on user-specified fixed coordinates to be given in Angstroms/radians
    one_stre->set_fixed_eq_val(R[i].eq_val/_bohr2angstroms);

    // check if coord is already present; returns 1 past the end if not found
    int index = find(one_stre);

    if (index == coords.simples.size())
      coords.simples.push_back(one_stre);// add it
    else { 
      coords.simples[index]->set_fixed_eq_val(R[i].eq_val/_bohr2angstroms); // it's there already, add the fixed value
      delete one_stre;
    }
  }
  // do fixed bends
  for (int i=0; i<B.size(); ++i) {
    int a = B[i].atoms[0];
    int b = B[i].atoms[1];
    int c = B[i].atoms[2];
    
    if (a >= natom || b >= natom || c >= natom)
      throw(INTCO_EXCEPT("Impossibly large index for atom in fixed bend string."));
    
    BEND *one_bend = new BEND(a, b, c, 0);
    // Insist on user-specified fixed coordinates to be given in Angstroms/radians
    one_bend->set_fixed_eq_val(B[i].eq_val/180.0*_pi);

    // check if coord is already present; returns 1 past the end if not found
    int index = find(one_bend);
  
    if (index == coords.simples.size())
      coords.simples.push_back(one_bend);// add it
    else {
      coords.simples[index]->set_fixed_eq_val(B[i].eq_val/180.0*_pi); // it's there already, add the fixed value
      delete one_bend;
    }
  }
  // do fixed dihedrals
  for (int i=0; i<D.size(); ++i) {
    int a = D[i].atoms[0];
    int b = D[i].atoms[1];
    int c = D[i].atoms[2];
    int d = D[i].atoms[3];

    if (a >= natom || b >= natom || c >= natom || d >= natom)
      throw(INTCO_EXCEPT("Impossibly large index for atom in fixed dihedral string."));

    TORS *one_tors = new TORS(a, b, c, d, 0);
    // Insist on user-specified fixed coordinates to be given in Angstroms/radians
    one_tors->set_fixed_eq_val(D[i].eq_val/180.0*_pi);

    // check if coord is already present; returns 1 past the end if not found
    int index = find(one_tors);

    if (index == coords.simples.size())
      coords.simples.push_back(one_tors);// add it
    else {
      coords.simples[index]->set_fixed_eq_val(D[i].eq_val/180.0*_pi); // it's there already, add the fixed value
      delete one_tors;
    }
  }

  return true;
}

template <typename T>
T StringToNumber ( const string & Text ) { 
  stringstream ss(Text);
  T result;
  return ss >> result ? result : -1;
}

std::vector<int> split_to_ints(string &str) {

  // Replace commas and ( and ) with spaces so that commas don't break it
  for (int i=0; i<str.size(); ++i) {
    if ( str[i] == ',' || str[i] == '(' || str[i] == ')' )
      str[i] = ' ';
  }

  char delim = ' ';
  std::stringstream ss(str);
  string item;
  std::vector<int> elems;

  while (std::getline(ss, item, delim)) {
    if (item.find_first_not_of(" ") != string::npos) { // Remove any empty entries (like first one)
      int a = StringToNumber<int>(item);
      if (a == -1) // change to int failed
        throw(INTCO_EXCEPT("Frozen atom string includes non-whole number."));
      elems.push_back(a-1); // start internal numbering at 0
    }
  }
  return elems;
}

// N = number of integers before each double/value
std::vector<fixed_coord> split_to_fixed_coord(string &str, int N) {

  // Replace commas and ( and ) with spaces so that commas don't break it
  for (int i=0; i<str.size(); ++i)
    if ( str[i] == ',' || str[i] == '(' || str[i] == ')' )
      str[i] = ' ';

  std::vector<fixed_coord> C;
  fixed_coord one_coord;

  char delim = ' ';
  std::stringstream ss(str);
  string item;
  int atom_cnt = 0;

  while (std::getline(ss, item, delim)) {
    if (item.find_first_not_of(" ") != string::npos) { // Remove any empty entries (like first one)

      if (atom_cnt < N) {
        int a = StringToNumber<int>(item);
        if (a == -1) // change to int failed
          throw(INTCO_EXCEPT("Fixed atoms string includes non-whole number for atom."));
        one_coord.atoms.push_back(a-1); // start internal numbering at 0
        ++atom_cnt;
      }
      else { // read eq val
        double val = StringToNumber<double>(item);
        if (val == -1) // change to float failed
          throw(INTCO_EXCEPT("Fixed atoms string includes non-float for value."));
        one_coord.eq_val = val; // start internal numbering at 0
        atom_cnt = 0;
        C.push_back(one_coord);  // save this coordinate and clear
        one_coord.eq_val = 0;
        one_coord.atoms.clear();
      }

    }
  }
  return C;
}

} // namespace opt

