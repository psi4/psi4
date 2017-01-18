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

/*! \file molecule_read_coords.cc
    \ingroup optking
    \brief read internal coordinates from file

     MOLECULE::read_coords(void) reads internal coordinates from text file
     return false if coordinates cannot be read - possibly b/c they are not there
*/

#include "molecule.h"

#include <iostream>
#include <sstream>

#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
#endif

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
    getline(fin, sline);
    if(sline.size() == 0)
      line_present = false;
    else
      line_present = true;
    // printf("getline read: %s\n", sline.c_str());
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

bool MOLECULE::read_coords(std::ifstream & fintco) {
  stringstream error;
  int line_num=0;
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
        if (fragments.back()->read_coord(vline, first_atom))
          line_present = myline(fintco, vline, line_num);
        else
          end_of_fragment = true; // break out of fragment reading lines
      }

      if (!line_present) // no more lines
        return true;
    }
    else if ((vline[0] == "C") || (vline[0] == "C*")) { // Combination coordinate
      //bool frozen = has_asterisk(vline[0]); // implement later ?
      int combo_length = 0;

      if (vline.size() != 2) {
        error << "Format of combination line header is \"C integer(number of lines) in line" << line_num << ".\n";
        throw(INTCO_EXCEPT(error.str().c_str()));
      }
      if (!stoi(vline[1], &combo_length)) {
        error << "Format of combination line header is \"C integer(number of lines) in line" << line_num << ".\n";
        throw(INTCO_EXCEPT(error.str().c_str()));
      }

      vector<int> cc_index;
      vector<double> cc_coeff;
      for (int i=0; i<combo_length; ++i) { // Read each simple id and coefficient
        line_present = myline(fintco, vline, line_num);
        int simple_id = 0;
        double simple_coeff = 0;

        if (!line_present) {
          error << "Missing line " << line_num << " in combination coordinate list.\n";
          throw(INTCO_EXCEPT(error.str().c_str()));
        }
        if (vline.size() != 2) {
          error << "Format of combination line \"integer(simple index)  double(coefficient)" << line_num << ".\n";
          throw(INTCO_EXCEPT(error.str().c_str()));
        }
        if (!stoi(vline[0], &simple_id)) {
          error << "Format of combination line \"integer(simple index)  double(coefficient)" << line_num << ".\n";
          throw(INTCO_EXCEPT(error.str().c_str()));
        }
        --simple_id;  // internal numbering starts with 0
        if (!stof(vline[1], &simple_coeff)) {
          error << "Format of combination line \"integer(simple index)  double(coefficient)" << line_num << ".\n";
          throw(INTCO_EXCEPT(error.str().c_str()));
        }

        cc_index.push_back(simple_id);
        cc_coeff.push_back(simple_coeff);
      }
      // Normalize
      double tval = 0.0;
      for (std::size_t j=0; j<cc_coeff.size(); ++j)
        tval += cc_coeff[j] * cc_coeff[j];
      tval = 1.0/sqrt(tval);
      for (std::size_t j=0; j<cc_coeff.size(); ++j)
        cc_coeff[j] *= tval;
      
      fragments.back()->add_combination_coord(cc_index, cc_coeff);

// read next line
      line_present = myline(fintco, vline, line_num);
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
      if ( first_frag >= (int) fragments.size() || second_frag >= (int) fragments.size()) {
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

          for (std::size_t i=1; i<vline.size(); ++i) {
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

      if (ndA > 0)
        for (std::size_t i=0; i<A1.size(); ++i)
          weightA[0][A1[i]] = 1.0/A1.size();
       
      if (ndA > 1)  
        for (std::size_t i=0; i<A2.size(); ++i) 
          weightA[1][A2[i]] = 1.0/A2.size();
       
      if (ndA > 2)  
        for (std::size_t i=0; i<A3.size(); ++i) 
          weightA[2][A3[i]] = 1.0/A3.size();
       
      if (ndB > 0)  
        for (std::size_t i=0; i<B1.size(); ++i) 
          weightB[0][B1[i]] = 1.0/B1.size();
       
      if (ndB > 1)  
        for (std::size_t i=0; i<B2.size(); ++i) 
          weightB[1][B2[i]] = 1.0/B2.size();
       
      if (ndB > 2)  
        for (std::size_t i=0; i<B3.size(); ++i) 
          weightB[2][B3[i]] = 1.0/B3.size();
       

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
} // end MOLECULE::read_coords


// this function belongs to the FRAG class - not MOLECULE but the reading of the intco
// definitions is so linked with that above, I'll put the function here
// reads internal coordinate definition line
// offset is the first atom number in the fragment so fragment can store relative numbering
bool FRAG::read_coord(vector<string> & s, int offset) {
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
    if (s[0] == "H") one_stre->set_hbond(true);
    if (has_eq_val) one_stre->set_fixed_eq_val(eq_val);

    if ( !present(one_stre) )
      coords.simples.push_back(one_stre);
    else
      delete one_stre;

    return true;
  }
  if (s[0] == "X") { // read cartesian coordinate
    if ( s.size() != 3 && s.size() != 4)
      throw(INTCO_EXCEPT("Format of cartesian entry is \"X atom_1 xyz\""));
    if ( s.size() == 4 ) {
      if (stof(s[3], &eq_val))
        has_eq_val = true;
      else
        throw(INTCO_EXCEPT("Format of cartesian entry is \"X atom_1 xyz (eq_val)\""));
    }
    if ( !stoi(s[1], &a) ) // Read atom
      throw(INTCO_EXCEPT("Format of cartesian entry is \"X atom_1 xyz\""));
    --a;

    // read x, y or z
    if      ( s[2] == "X" ) b = 0;
    else if ( s[2] == "Y" ) b = 1;
    else if ( s[2] == "Z" ) b = 2;
    else throw(INTCO_EXCEPT("Format of cartesian entry is \"X atom_1 xyz\""));

    CART *one_cart = new CART(a-offset, b, frozen);
    if (has_eq_val) one_cart->set_fixed_eq_val(eq_val);

    if ( !present(one_cart) )
      coords.simples.push_back(one_cart);
    else
      delete one_cart;

    return true;
  }
  else if ((s[0] == "B") || (s[0] == "L") || s[0] == "l") {
    if (s.size() != 4 && s.size() != 5)
      throw(INTCO_EXCEPT("Format of bend entry is \"B[L,l] atom_1 atom_2 atom_3\""));
    if ( s.size() == 5 ) {
      if (stof(s[4], &eq_val))
        has_eq_val = true;
      else
        throw(INTCO_EXCEPT("Format of bend entry is \"B[L,l] atom_1 atom_2 atom_3 (eq_val)\""));
    }
    if ( !stoi(s[1], &a) || !stoi(s[2], &b) || !stoi(s[3], &c) )
      throw(INTCO_EXCEPT("Format of bend entry is \"B[L,l] atom_1 atom_2 atom_3\""));
    --a; --b; --c;

    BEND *one_bend = new BEND(a-offset, b-offset, c-offset, frozen);

    if (s[0] == "L") one_bend->make_lb_normal();
    else if (s[0] == "l") one_bend->make_lb_complement();

    if (has_eq_val) one_bend->set_fixed_eq_val(eq_val);

    if ( !present(one_bend) )
      coords.simples.push_back(one_bend);
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
      coords.simples.push_back(one_tors);
    else
      delete one_tors;

    return true;
  }
  else if (s[0] == "O") {
    if (s.size() != 5 && s.size() != 6)
      throw(INTCO_EXCEPT("Format of out-of-plane entry is \"O atom_1 atom_2 atom_3 atom_4\""));
    if ( s.size() == 6 ) {
      if (stof(s[5], &eq_val))
        has_eq_val = true;
      else
        throw(INTCO_EXCEPT("Format of out_of_plane entry is \"O atom_1 atom_2 atom_3 atom_4 (eq_val)\""));
    }
    if ( !stoi(s[1], &a) || !stoi(s[2], &b) || !stoi(s[3], &c) || !stoi(s[4], &d) )
      throw(INTCO_EXCEPT("Format of out_of_plane entry is \"O atom_1 atom_2 atom_3 atom_4\""));
    --a; --b; --c; --d;

    OOFP *one_oofp = new OOFP(a-offset, b-offset, c-offset, d-offset, frozen);
    if (has_eq_val) one_oofp->set_fixed_eq_val(eq_val);
   
    if ( !present(one_oofp) )
      coords.simples.push_back(one_oofp);
    else
      delete one_oofp;

    return true;
  }

  //printf("read_coord returning false\n");
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
  double d;
  try {
    d = std::stod(s, NULL);
  }
  catch(...) {
    return false;
  }
  *val = d;
  return true;
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
