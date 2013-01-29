/*! \file frag_apply_constraint_list.cc
    \ingroup optking
    \brief read constraints from string for fragment

     FRAG::add_constraint_list(string) parses string for the atoms that define bonds, bends, and angles
     to be frozen.  Return false if there are none.
*/

#include "frag.h"

#include <string>
#include <sstream>
#include <vector>

#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
#endif

namespace opt {

using std::string;

std::vector<int> split_to_ints(string &s);

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
    fprintf(outfile,"\tFrozen distance atom list: \n");
    for (int i=0; i<R_atoms.size(); i+=2)
      fprintf(outfile,"\t %5d %5d\n", R_atoms[i]+1, R_atoms[i+1]+1);
  }
  if (B_atoms.size()) {
    fprintf(outfile,"\tFrozen bend atom list: \n");
    for (int i=0; i<B_atoms.size(); i+=3)
      fprintf(outfile,"\t %5d %5d %5d\n", B_atoms[i]+1, B_atoms[i+1]+1, B_atoms[i+2]+1);
  }
  if (D_atoms.size()) {
    fprintf(outfile,"\tFrozen dihedral atom list: \n");
    for (int i=0; i<D_atoms.size(); i+=4)
      fprintf(outfile,"\t %5d %5d %5d %5d\n", D_atoms[i]+1, D_atoms[i+1]+1,
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

    // check if intco is already present; returns 1 past the end if not found
    int index = find(one_stre);

    if (index == intcos.size())
      intcos.push_back(one_stre);// add it
    else { 
      intcos[index]->freeze();   // it's there already; make sure it's frozen
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

    // check if intco is already present; returns 1 past the end if not found
    int index = find(one_bend);

    if (index == intcos.size())
      intcos.push_back(one_bend);// add it
    else { 
      intcos[index]->freeze();   // it's there already; make sure it's frozen
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

    // check if intco is already present; returns 1 past the end if not found
    int index = find(one_tors);

    if (index == intcos.size())
      intcos.push_back(one_tors);// add it
    else { 
      intcos[index]->freeze();   // it's there already; make sure it's frozen
      delete one_tors;
    }
  }
  return true;
}

std::vector<int> split_to_ints(string &str) {
  // Replace commas and ( and ) with spaces so that commas don't break it
  size_t pos = 0;
  while ((pos = str.find(",", pos)) != string::npos) {
     str.replace(pos, 1, " ");
     pos += 1;
  }
  pos = 0;
  while ((pos = str.find(")", pos)) != string::npos) {
     str.replace(pos, 1, " ");
     pos += 1;
  }
  pos = 0;
  while ((pos = str.find("(", pos)) != string::npos) {
     str.replace(pos, 1, " ");
     pos += 1;
  }

  char delim = ' ';
  std::stringstream ss(str);
  string item;
  std::vector<int> elems;

  while (std::getline(ss, item, delim)) {
    if (item.find_first_not_of(" ") != string::npos) { // Remove any empty entries (like first one)
      int a = atoi(item.c_str());
      if (!a) // change to int failed
        throw(INTCO_EXCEPT("Frozen atom string includes non-whole number."));
      elems.push_back(a-1); // start internal numbering at 0
    }
  }
  return elems;
}

} // namespace opt

