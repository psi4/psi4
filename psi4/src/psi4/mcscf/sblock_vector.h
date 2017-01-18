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

#ifndef _psi_src_lib_libmemtrix_sblock_vector_h_
#define _psi_src_lib_libmemtrix_sblock_vector_h_

#include <string>

#include "block_vector.h"

namespace psi{ namespace mcscf{

// Smart version of BlockVector
class SBlockVector
{
public:
  SBlockVector();
  SBlockVector(std::string label, int nirreps, int*& rows_size);
  SBlockVector(std::string label, int nirreps, vecint& rows_size);
  ~SBlockVector() {if(block_vector_)
                     if(block_vector_->subtract_reference())
                          block_vector_ = 0;}

  // Manual allocation
  void allocate(std::string label, int nirreps, int*& rows_size);
  void allocate(std::string label, int nirreps, vecint& rows_size);

  void subtract_reference(){
      if(block_vector_){
          if(block_vector_->subtract_reference())
              block_vector_ = 0;
      }
  }

  // Copy constructor and assignment operator
  SBlockVector            (const   SBlockVector& src);
  SBlockVector& operator= (const   SBlockVector& src);

  // Allow access to the implementation object
  const BlockVector* operator-> () const {return block_vector_;}
  BlockVector*       operator-> ()       {return block_vector_;}

  // Access the implementation object
  BlockVector* getBlockVector() {return block_vector_;}

  // Checking functions
  bool is_allocated() {return (block_vector_);}
  void check(const char* cstr);

  void copy(SBlockVector& source);
private:
  SBlockVector(BlockVector* block_vector);

  BlockVector* block_vector_;
};

}}

#endif // _psi_src_lib_libmemtrix_sblock_vector_h_