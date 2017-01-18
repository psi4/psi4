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

#ifndef _psi_src_lib_libmemtrix_block_vector_h_
#define _psi_src_lib_libmemtrix_block_vector_h_

#include <string>
#include <vector>

#include "vector_base.h"
#include "psi4/libpsi4util/memory_manager.h"

typedef std::vector<int> vecint;

namespace psi{ namespace mcscf{

class BlockVector
{
public:
  BlockVector();
  BlockVector(std::string label, int nirreps, size_t*& rows_size);
  BlockVector(std::string label, int nirreps, int*& rows_size);
  BlockVector(std::string label, int nirreps, vecint& rows_size);
  ~BlockVector();

  void print();
  void copy(BlockVector& source);

  // Inlines
  void        set(int h, int i, double value) {vector_base_[h]->set(i,value);}
  void        add(int h, int i, double value) {vector_base_[h]->add(i,value);}
  double      get(int h, int i)               {return(vector_base_[h]->get(i));}

  VectorBase* getVectorBase(int h) {return(vector_base_[h]);}

  // Reference counting related
  unsigned int ref ()  const { return ref_;}   // Number of references
  void add_reference      () { ref_++;}
  bool subtract_reference () { if (--ref_ == 0){ delete this; return true;} return false;}
  // Reference count
  unsigned int ref_;
private:
  // Vector label and pointer
  std::string label_;
  VectorBase** vector_base_;

  // Block sizes etc.
  size_t*  rows_size_;
  size_t*  rows_offset_;
  int   nirreps_;

  void startup(std::string label, int nirreps, size_t*& rows_size);
  void startup(std::string label, int nirreps, int*& rows_size);
  void startup(std::string label, int nirreps, vecint& rows_size);
  void cleanup();
};

}}

#endif // _psi_src_lib_libmemtrix_block_vector_h_
