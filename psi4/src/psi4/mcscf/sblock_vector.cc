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

#include <cstdlib>
#include "psi4/psifiles.h"
#include <cstdio>

#include "sblock_vector.h"

#include "psi4/psi4-dec.h"

namespace psi{ namespace mcscf{

SBlockVector::SBlockVector()
 : block_vector_(0)
{
}

SBlockVector::SBlockVector(std::string label, int nirreps, int*& rows_size)
 : block_vector_(0)
{
  block_vector_ = new BlockVector(label,nirreps,rows_size);
  block_vector_->add_reference();
}

SBlockVector::SBlockVector(std::string label, int nirreps, vecint& rows_size)
 : block_vector_(0)
{
  block_vector_ = new BlockVector(label,nirreps,rows_size);
  block_vector_->add_reference();
}

SBlockVector::SBlockVector(BlockVector* block_vector)
 : block_vector_(block_vector)
{
  block_vector_->add_reference();
}

SBlockVector::SBlockVector(const SBlockVector& src)
{
  block_vector_ = src.block_vector_;
  block_vector_->add_reference();
}

void SBlockVector::allocate(std::string label, int nirreps, int*& rows_size)
{
  block_vector_ = new BlockVector(label,nirreps,rows_size);
  block_vector_->add_reference();
}

void SBlockVector::allocate(std::string label, int nirreps, vecint& rows_size)
{
  block_vector_ = new BlockVector(label,nirreps,rows_size);
  block_vector_->add_reference();
}

SBlockVector& SBlockVector::operator= (const SBlockVector& src)
{
  // Make sure we don't copy ourself!
  if (block_vector_ == src.block_vector_) return *this;

  block_vector_->subtract_reference();  // Remove reference from existing object
  block_vector_ = src.block_vector_;
  block_vector_->add_reference();       // Add reference to our new object

  return *this;
}

void SBlockVector::check(const char* cstr)
{
  if(!is_allocated()){
    outfile->Printf("\n\n  Error: SBlockVector operation '%s' is using an uninitialized matrix",cstr);

    exit(PSI_RETURN_FAILURE);
  }
}

void SBlockVector::copy(SBlockVector& source)
{
  block_vector_->copy(*source.block_vector_);
}

}}
