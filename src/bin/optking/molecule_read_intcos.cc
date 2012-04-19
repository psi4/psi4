/*! \file molecule_read_intcos.cc
    \ingroup optking
    \brief read internal coordinates from file

     MOLECULE::read_intcos(void) reads internal coordinates from text file
     return false if coordinates cannot be read - possibly b/c they are not there
*/

#include "molecule.h"

#include <cmath>
#include <iostream>
#include <sstream>

#define EXTERN
#include "globals.h"
// maximum number of atoms that can be used to define a referent point
#define MAX_REF_ATOM_INTCO (20)

namespace opt {

using namespace std;

// convert string to integer
bool stoi(string s, int *a);

// convert string to boolean
bool stob(string s, bool *a);

// convert string to float
bool stof(string s, double *val);

// check string for trailing "*".  If present, remove it and return true
bool has_asterisk(string & s);

bool myline(ifstream & fin, vector<string> & tokens, int & line_num);

// clear tokens
// read line and tokenize it into tokens
// lines beginning with % or empty lines are skipped/ignored
// return false, if no more lines
bool myline(ifstream & fin, vector<string> & tokens, int & line_num) {
  string sline;
  stringstream streamline;
  bool read_next = true;
  bool line_present = false;

  tokens.clear();

  while (read_next && !fin.eof()) {
    line_present = getline(fin, sline);
    //printf("getline read: %s\n", sline.c_str());
    read_next = false;

    if (line_present) {
      ++line_num;
      streamline << sline;
      while(streamline >> sline)
        tokens.push_back(sline);

      if (tokens.empty()) {
        tokens.clear();
        streamline.clear();
        read_next = true;
      }
      else if (tokens[0][0] == '%') {
        tokens.clear();
        streamline.clear();
        read_next = true;
      }
      else
        return true;
    }
  }
  //printf("myline returns false\n");
  return false;
}

bool MOLECULE::read_intcos(std::ifstream & fintco) {
  stringstream error;
  int a, b, natom, line_num=0;
  bool D_on[6];     // interfragment coordinates active
  bool D_frozen[6]; // interfragment coordinates frozen
  FRAG * frag1;
  int first_atom, last_atom;
  int first_frag, second_frag;
  vector<string> vline;
  vector<int> A1; vector<int> A2; vector<int> A3;
  vector<int> B1; vector<int> B2; vector<int> B3;
  int ndA=0, ndB=0;

  bool line_present = false;

  // read in line and tokenize
  line_present = myline(fintco, vline, line_num);

  int cnt =0;

  while (line_present) {

    if ((vline[0] == "F") || (vline[0] == "F*")) { // read first and last atom for fragment
      bool frozen = has_asterisk(vline[0]);
       
      if (vline.size() != 3) {
        error << "Format of fragment line is \"F integer(1st_atom) integer(last_atom)\", line " << line_num << ".\n";
        throw(INTCO_EXCEPT(error.str().c_str()));
      }
      if ( !stoi(vline[1], &first_atom) || !stoi(vline[2], &last_atom)) {
        error << "Format of fragment line is \"F integer(1st_atom) integer(last_atom)\", line " << line_num << ".\n";
        throw(INTCO_EXCEPT(error.str().c_str()));
      }
      --first_atom; // start at 0 internally
      --last_atom;
      if (first_atom > last_atom) {
        error << "Last atom must be greater than first atom, line " << line_num << ".\n";
        throw(INTCO_EXCEPT(error.str().c_str()));
      }

      // create fragment
      frag1 = new FRAG(last_atom - first_atom + 1);
      if (frozen) frag1->freeze();
      fragments.push_back(frag1);

      // now read all internal coordinates for that fragment
      line_present = myline(fintco, vline, line_num);
      bool end_of_fragment = false;
      while(line_present && !end_of_fragment) {
        if (fragments.back()->read_intco(vline, first_atom))
          line_present = myline(fintco, vline, line_num);
        else
          end_of_fragment = true; // break out of fragment reading lines
      }

      if (!line_present) // no more lines
        return true;
    }
    else if (vline[0] == "I" || vline[0] == "I*") { // read interfragment definition
      bool frozen_I = has_asterisk(vline[0]);
      if (vline.size() != 3) {
        error << "Format of interfragment line is \"I integer(1st_frag) integer(2nd_frag)\" in line " << line_num << ".\n";
        throw(INTCO_EXCEPT(error.str().c_str()));
      }
      if ( !stoi(vline[1], &first_frag) || !stoi(vline[2], &second_frag) ) {
        error << "Format of interfragment line is \"I integer(1st_frag) integer(2nd_frag)\" in line " << line_num << ".\n";
        throw(INTCO_EXCEPT(error.str().c_str()));
      }
      --first_frag;
      --second_frag;
      if ( first_frag >= fragments.size() || second_frag >= fragments.size()) {
        error << "Fragments can only be referenced if already defined, error in line " << line_num << ".\n";
        throw(INTCO_EXCEPT(error.str().c_str()));
      }

      // read mandatory line with 6 booleans
      line_present = myline(fintco, vline, line_num);
      if (!line_present) {
        error << "Missing line " << line_num+1 << " to indicate with six 1/0's which coordinates are active.\n";
        throw(INTCO_EXCEPT(error.str().c_str()));
      }
      if (vline.size() != 6 ) {
        error << "Indicate with _six_ 1/0's which coordinates are active, error in line " << line_num << ".\n";
        throw(INTCO_EXCEPT(error.str().c_str()));
      }
      for (int i=0; i<6; ++i) {
        D_frozen[i] = has_asterisk(vline[i]);
        if (!stob(vline[i], &(D_on[i]))) {
          error << "Indicate with six 1/0's which coordinates are active, error in line " << line_num << ".\n";
        throw(INTCO_EXCEPT(error.str().c_str()));
        }
      }
      if (frozen_I) { // freeze fragment; override other
        for (int i=0; i<6; ++i) D_frozen[i] = true;
      }

      A1.clear(); A2.clear(); A3.clear();
      B1.clear(); B2.clear(); B3.clear();

      line_present = myline(fintco, vline, line_num);
      bool end_of_if = false;

      while (line_present && !end_of_if) {

        if (vline[0] == "A1" || vline[0] == "A2" || vline[0] == "A3" ||
            vline[0] == "B1" || vline[0] == "B2" || vline[0] == "B3") {

          if (vline.size() < 2) {
            error << "Missing atoms list for reference atom, error in line" << line_num << ".\n";
          throw(INTCO_EXCEPT(error.str().c_str()));
          }

          for (int i=1; i<vline.size(); ++i) {
            int a;
            if (!stoi(vline[i], &a)) {
              error << "Could not read atom list, error in line" << line_num << ".\n";
            throw(INTCO_EXCEPT(error.str().c_str()));
            }
            a--;

            if (vline[0] == "A1")      A1.push_back(a - g_atom_offset(first_frag));
            else if (vline[0] == "A2") A2.push_back(a - g_atom_offset(first_frag));
            else if (vline[0] == "A3") A3.push_back(a - g_atom_offset(first_frag));
            else if (vline[0] == "B1") B1.push_back(a - g_atom_offset(second_frag));
            else if (vline[0] == "B2") B2.push_back(a - g_atom_offset(second_frag));
            else if (vline[0] == "B3") B3.push_back(a - g_atom_offset(second_frag));
          }

        }
        else { // done reading reference points
          end_of_if = true;
        }

        line_present = myline(fintco, vline, line_num);
      }

      // Done reading reference atoms, now error-check
      // Determine which reference atoms are needed for interfragment coodinates requested
      bool ref_A_needed[3], ref_B_needed[3];
      for (int i=0; i<3; ++i) {
        ref_A_needed[i] = false;
        ref_B_needed[i] = false;
      }

      if (D_on[0]) {
        ref_A_needed[0] = true; // A1 
        ref_B_needed[0] = true; // B1 
      }
      if (D_on[1]) { //theta_A
        ref_A_needed[1] = true; // A2 
        ref_A_needed[0] = true; // A1 
        ref_B_needed[0] = true; // B1 
      }
      if (D_on[2]) { //theta_B
        ref_A_needed[0] = true; // A1 
        ref_B_needed[0] = true; // B1 
        ref_B_needed[1] = true; // B2 
      }
      if (D_on[3]) { //tau
        ref_A_needed[1] = true; // A2 
        ref_A_needed[0] = true; // A1 
        ref_B_needed[0] = true; // B1 
        ref_B_needed[1] = true; // B2 
      }
      if (D_on[4]) { //phi_A
        ref_A_needed[2] = true; // A3 
        ref_A_needed[1] = true; // A2 
        ref_A_needed[0] = true; // A1 
        ref_B_needed[0] = true; // B1 
      }
      if (D_on[5]) { //phi_B
        ref_A_needed[0] = true; // A1 
        ref_B_needed[0] = true; // B1 
        ref_B_needed[1] = true; // B2 
        ref_B_needed[2] = true; // B3 
      }

      if ((ref_A_needed[0] && A1.empty()) || (ref_A_needed[1] && A2.empty()) || (ref_A_needed[2] && A3.empty())
        || (ref_B_needed[0] && B1.empty()) || (ref_B_needed[1] && B2.empty()) || (ref_B_needed[2] && B3.empty())) {
          error << "Not all necessary reference atom supplied, error in line" << line_num << ".\n";
          throw(INTCO_EXCEPT(error.str().c_str()));
      }

      ndA = ndB = 0;
      if (!A1.empty()) ++ndA;
      if (!A2.empty()) ++ndA;
      if (!A3.empty()) ++ndA;
      if (!B1.empty()) ++ndB;
      if (!B2.empty()) ++ndB;
      if (!B3.empty()) ++ndB;

      double **weightA = init_matrix(ndA, fragments[first_frag]->g_natom());
      double **weightB = init_matrix(ndB, fragments[second_frag]->g_natom());

      if (ndA > 0) {
        for (int i=0; i<A1.size(); ++i)
          weightA[0][A1[i]] = 1.0;
      }
      if (ndA > 1) {
        for (int i=0; i<A2.size(); ++i) 
          weightA[1][A2[i]] = 1.0;
      }
      if (ndA > 2) {
        for (int i=0; i<A3.size(); ++i) 
          weightA[2][A3[i]] = 1.0;
      }
      if (ndB > 0) {
        for (int i=0; i<B1.size(); ++i) 
          weightB[0][B1[i]] = 1.0;
      }
      if (ndB > 1) {
        for (int i=0; i<B2.size(); ++i) 
          weightB[1][B2[i]] = 1.0;
      }
      if (ndB > 2) {
        for (int i=0; i<B3.size(); ++i) 
          weightB[2][B3[i]] = 1.0;
      }

      INTERFRAG * one_IF = new INTERFRAG(fragments[first_frag], fragments[second_frag],
      first_frag, second_frag, weightA, weightB, ndA, ndB);

      // freeze coordinates, need to convert to the index of active coordinates
      int if_index = 0;
      for (int i=0; i<6; ++i) {
        if (D_on[i]) {
          if (D_frozen[i])
            one_IF->freeze(if_index++);
        }
      }

      interfragments.push_back(one_IF);

    } // end of if vline[0] == 'I'
    else {
      error << "Unknown initial character on line " << line_num << ".\n";
      throw(INTCO_EXCEPT(error.str().c_str()));
    }
  }
  
  return true;
} // end MOLECULE::read_intcos


// this function belongs to the FRAG class - not MOLECULE but the reading of the intco
// definitions is so linked with that above, I'll put the function here
// reads internal coordinate definition line
// offset is the first atom number in the fragment so fragment can store relative numbering
bool FRAG::read_intco(vector<string> & s, int offset) {
  int a, b, c, d;
  std::string error;
  double eq_val; // for imposing a constraint
  bool has_eq_val=false;

  bool frozen=false;
  frozen = has_asterisk(s[0]); // removes asterisk

  if ((s[0] == "R") || (s[0] == "H")) {
    if ( s.size() != 3 && s.size() != 4)
      throw(INTCO_EXCEPT("Format of stretch entry is \"R atom_1 atom_2\""));
    if ( s.size() == 4 ) {
      if (stof(s[3], &eq_val))
        has_eq_val = true;
      else
        throw(INTCO_EXCEPT("Format of stretch entry is \"R atom_1 atom_2 (eq_val)\""));
    }
    if ( !stoi(s[1], &a) || !stoi(s[2], &b) )
      throw(INTCO_EXCEPT("Format of stretch entry is \"R atom_1 atom_2\""));
    --a; --b;

    STRE *one_stre = new STRE(a-offset, b-offset, frozen);
    if (s[0] == "H") one_stre->make_hbond();
    if (has_eq_val) one_stre->set_fixed_eq_val(eq_val);

    if ( !present(one_stre) )
      intcos.push_back(one_stre);
    else
      delete one_stre;

    return true;
  }
  else if ((s[0] == "B") || (s[0] == "L")) {
    if (s.size() != 4 && s.size() != 5)
      throw(INTCO_EXCEPT("Format of bend entry is \"B atom_1 atom_2 atom_3\""));
    if ( s.size() == 5 ) {
      if (stof(s[4], &eq_val))
        has_eq_val = true;
      else
        throw(INTCO_EXCEPT("Format of bend entry is \"B atom_1 atom_2 atom_3 (eq_val)\""));
    }
    if ( !stoi(s[1], &a) || !stoi(s[2], &b) || !stoi(s[3], &c) )
      throw(INTCO_EXCEPT("Format of bend entry is \"B atom_1 atom_2 atom_3\""));
    --a; --b; --c;

    BEND *one_bend = new BEND(a-offset, b-offset, c-offset, frozen);
    if (s[0] == "L") one_bend->make_linear_bend();
    if (has_eq_val) one_bend->set_fixed_eq_val(eq_val);

    if ( !present(one_bend) )
      intcos.push_back(one_bend);
    else
      delete one_bend;

    return true;
  }
  else if (s[0] == "D") {
    if (s.size() != 5 && s.size() != 6)
      throw(INTCO_EXCEPT("Format of dihedral entry is \"D atom_1 atom_2 atom_3 atom_4\""));
    if ( s.size() == 6 ) {
      if (stof(s[5], &eq_val))
        has_eq_val = true;
      else
        throw(INTCO_EXCEPT("Format of dihedral entry is \"D atom_1 atom_2 atom_3 atom_4 (eq_val)\""));
    }
    if ( !stoi(s[1], &a) || !stoi(s[2], &b) || !stoi(s[3], &c) || !stoi(s[4], &d) )
      throw(INTCO_EXCEPT("Format of dihedral entry is \"D atom_1 atom_2 atom_3 atom_4\""));
    --a; --b; --c; --d;

    TORS *one_tors = new TORS(a-offset, b-offset, c-offset, d-offset, frozen);
    if (has_eq_val) one_tors->set_fixed_eq_val(eq_val);
    
    if ( !present(one_tors) )
      intcos.push_back(one_tors);
    else
      delete one_tors;

    return true;
  }
  //printf("read_intco returning false\n");
  return false;
}

// convert string to integer
bool stoi(string s, int *a) {
  int i = atoi(s.c_str());
  if (i!=0) {
    *a = i;
    return true;
  }
  return false;
}

// convert string to float
bool stof(string s, double *val) {
  double i = atof(s.c_str());
  if (i!=0) {
    *val = i;
    return true;
  }
  return false;
}

// convert string to boolean
bool stob(string s, bool *a) {
  if (s == "1") {
    *a = true;
    return true;
  }
  else if (s == "0") {
    *a = false;
    return true;
  }
  else 
   return false;
}

// removes asterisk and returns true if it was present
bool has_asterisk(string & s) {
  if (s[s.size()-1] == '*') {
    s.erase(s.size()-1);
    return true;
  }
  else
    return false;
}

} // namespace opt

