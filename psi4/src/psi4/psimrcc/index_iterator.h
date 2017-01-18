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

#ifndef _psi_src_bin_psimrcc_index_iterator_h_
#define _psi_src_bin_psimrcc_index_iterator_h_

/*! \file    index_iterator.h
    \ingroup (PSIMRCC)
    \brief   This class is used to iterate over n-tuples of MOs indices (p,q,r,..)
*/

#include <string>

namespace psi{ namespace psimrcc{

class CCIndex;

class CCIndexIterator{
public:
  // Class Constructor and Destructor
  explicit CCIndexIterator(std::string str);
  explicit CCIndexIterator(std::string str,int select_irrep);
  explicit CCIndexIterator(CCIndex* index);
  explicit CCIndexIterator(CCIndex* index,int select_irrep);
  ~CCIndexIterator();

  // Class Public Methods
  bool first();
  void next();
  bool end()   {return(absolute >= max_abs);}

  template <int N>
  short ind_abs() {return tuples[absolute][N];}

  int   sym()     {return symmetry;}
  size_t rel()    {return relative;}
  size_t abs()    {return absolute;}
private:
  // Class private functions
  void        startup(int min_sym,int max_sym);

  // Generica data
  int                               nirreps;

  // Index object
  CCIndex*                          ccindex;

  // Internal iterator
  size_t                            relative;  // Relative address of the current tuple
  size_t                            absolute;      // Absolute address of the current tuple
  size_t                            max_abs;  // Max absolute address of the current tuple
  size_t                            min_abs;  // Min absolute address of the current tuple
  int                               symmetry; // Symmetry of the current tuple
  int                               current_block;

  // Properties of the tuples
  int                               nelements;
  int**                             element_irrep;
  short**                           tuples;
  std::vector<size_t>               block_last;
  std::vector<int>                  block_symmetry;
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_index_iterator_h_