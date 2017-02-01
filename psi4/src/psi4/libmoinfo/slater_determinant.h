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

#ifndef _psi_src_lib_libmoinfo_slater_determinant_h_
#define _psi_src_lib_libmoinfo_slater_determinant_h_

/*! \file    slater_determinant.h
    \ingroup LIBMOINFO
    \brief   This class stores all the basic info regarding a Slater determinant
*/

// C LIBS
#include <cstdlib>
// STL
#include <vector>
#include <string>

namespace psi {

class SlaterDeterminant{
public:
//  SlaterDeterminant();
//  SlaterDeterminant(SlaterDeterminant& det);
  SlaterDeterminant(int alfa_sym_,int beta_sym_,std::vector<bool> alfa_bits_,std::vector<bool> beta_bits_);
  ~SlaterDeterminant();

  // Get functions
  int get_alfa_sym() const {return alfa_sym;}
  int get_beta_sym() const {return beta_sym;}
  size_t get_alfa_string() const {return alfa_string;}
  size_t get_beta_string() const {return beta_string;}
  std::vector<bool> get_alfa_bits() const {return alfa_bits;}
  std::vector<bool> get_beta_bits() const {return beta_bits;}

  // Set functions
  void set_alfa_bits(std::vector<bool> alfa_bits_) {alfa_bits = alfa_bits_;}
  void set_beta_bits(std::vector<bool> beta_bits_) {beta_bits = beta_bits_;}

  // Properties
  bool is_closed_shell();

  // Print functions
  std::string get_label();
private:
  // Class private functions
  void startup();
  void cleanup();
  char get_occupation_symbol(int i);
  // Class private data
  int    alfa_sym;              // Symmetry of the alfa string
  int    beta_sym;              // Symmetry of the beta string
  size_t alfa_string;           // Address of the alfa string
  size_t beta_string;           // Address of the beta string
  std::vector<bool> alfa_bits;  // Bit representation of the alfa string
  std::vector<bool> beta_bits;  // Bit representation of the beta string
};

}

#endif // _psi_src_lib_libmoinfo_slater_determinant_h_