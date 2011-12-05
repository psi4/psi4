#include <cstdlib>
#include <cstdio>

#include <psifiles.h>

#include "sblock_matrix.h"

#include <psi4-dec.h>

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
    fprintf(outfile,"\n\n  Error: SBlockMatrix operation '%s' is using an uninitialized matrix",cstr);
    fflush(outfile);
    exit(PSI_RETURN_FAILURE);
  }
}

}}
