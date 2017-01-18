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
#include <cstdio>

#include "psi4/psifiles.h"

#include "sblock_matrix.h"

#include "psi4/psi4-dec.h"

namespace psi{ namespace mcscf{

SBlockMatrix::SBlockMatrix()
 : block_matrix_(0)
{
}

SBlockMatrix::SBlockMatrix(std::string label, int nirreps, size_t*& rows_size, size_t*& cols_size)
 : block_matrix_(0)
{
  block_matrix_ = new BlockMatrix(label,nirreps,rows_size,cols_size);
  block_matrix_->add_reference();
}

SBlockMatrix::SBlockMatrix(std::string label, int nirreps, int*& rows_size, int*& cols_size)
 : block_matrix_(0)
{
  block_matrix_ = new BlockMatrix(label,nirreps,rows_size,cols_size);
  block_matrix_->add_reference();
}

SBlockMatrix::SBlockMatrix(std::string label, int nirreps, vecint& rows_size, vecint& cols_size)
 : block_matrix_(0)
{
  block_matrix_ = new BlockMatrix(label,nirreps,rows_size,cols_size);
  block_matrix_->add_reference();
}

SBlockMatrix::SBlockMatrix(BlockMatrix* block_matrix)
 : block_matrix_(block_matrix)
{
  block_matrix_->add_reference();
}

SBlockMatrix::SBlockMatrix(SBlockMatrix& src)
{
  block_matrix_ = src.block_matrix_;
  block_matrix_->add_reference();
}

void SBlockMatrix::allocate(std::string label, int nirreps, size_t*& rows_size, size_t*& cols_size)
{
  block_matrix_ = new BlockMatrix(label,nirreps,rows_size,cols_size);
  block_matrix_->add_reference();
}

void SBlockMatrix::allocate(std::string label, int nirreps, int*& rows_size, int*& cols_size)
{
  block_matrix_ = new BlockMatrix(label,nirreps,rows_size,cols_size);
  block_matrix_->add_reference();
}

void SBlockMatrix::allocate(std::string label, int nirreps, vecint& rows_size, vecint& cols_size)
{
  block_matrix_ = new BlockMatrix(label,nirreps,rows_size,cols_size);
  block_matrix_->add_reference();
}

SBlockMatrix& SBlockMatrix::operator+= (SBlockMatrix& src)
{
  check("operator+=");  src.check("operator+=");
  *(block_matrix_) += *(src.block_matrix_);
  return *this;
}

SBlockMatrix& SBlockMatrix::operator-= (SBlockMatrix& src)
{
  check("operator-=");  src.check("operator-=");
  *(block_matrix_) -= *(src.block_matrix_);
  return *this;
}

SBlockMatrix& SBlockMatrix::operator= (SBlockMatrix& src)
{
  check("operator=");  src.check("operator=");
  *(block_matrix_) = *(src.block_matrix_);
  // Make sure we don't copy ourself!
/*  if (block_matrix_ == src.block_matrix_) return *this;

  block_matrix_->subtract_reference();  // Remove reference from existing object
  block_matrix_ = src.block_matrix_;
  block_matrix_->add_reference();       // Add reference to our new object
*/
  return *this;
}

// void SBlockMatrix::copy(SBlockMatrix& src)
// {
//   block_matrix_->subtract_reference();
//   block_matrix_ = new BlockMatrix();
//   *block_matrix_ = *src.block_matrix_;
//   block_matrix_->add_reference();       // Add reference to our new object
// }

void SBlockMatrix::multiply(bool transpose_A, bool transpose_B, SBlockMatrix& A, SBlockMatrix& B)
{
  check("multiply");  A.check("multiply");  B.check("multiply");
  block_matrix_->multiply(transpose_A,transpose_B,A.getBlockMatrix(),B.getBlockMatrix());
}

void SBlockMatrix::diagonalize(SBlockMatrix& eigenmatrix,SBlockVector& eigenvalues)
{
  check("diagonalize"); eigenmatrix.check("diagonalize");  eigenvalues.check("multiply");
  block_matrix_->diagonalize(eigenmatrix.getBlockMatrix(),eigenvalues.getBlockVector());
}

double dot(SBlockMatrix& A,SBlockMatrix& B)
{
  A.check("dot");  B.check("dot");
  return( dot(A.getBlockMatrix(),B.getBlockMatrix()) );
}

void SBlockMatrix::check(const char* cstr)
{
  if(!is_allocated()){
    outfile->Printf("\n\n  Error: SBlockMatrix operation '%s' is using an uninitialized matrix",cstr);

    exit(PSI_RETURN_FAILURE);
  }
}

}}
