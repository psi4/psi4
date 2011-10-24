#include <libutil/libutil.h>
#include <cstdio>

#include <libutil/memory_manager.h>
#include "block_matrix.h"
#include "matrix_base.h"

#include <psi4-dec.h>

namespace psi{ namespace mcscf{

BlockMatrix::BlockMatrix()
 : nirreps_(0) ,ref_(0), matrix_base_(0), rows_size_(0), cols_size_(0), rows_offset_(0), cols_offset_(0)
{
}


BlockMatrix::BlockMatrix(std::string label, int nirreps, size_t*& rows_size, size_t*& cols_size)
 : ref_(0), matrix_base_(0), rows_size_(0), cols_size_(0), rows_offset_(0), cols_offset_(0)
{
  startup(label,nirreps,rows_size,cols_size);
}

BlockMatrix::BlockMatrix(std::string label, int nirreps, int*& rows_size, int*& cols_size)
 : ref_(0), matrix_base_(0), rows_size_(0), cols_size_(0), rows_offset_(0), cols_offset_(0)
{
  startup(label,nirreps,rows_size,cols_size);
}

BlockMatrix::BlockMatrix(std::string label, int nirreps, vecint& rows_size, vecint& cols_size)
 : ref_(0), matrix_base_(0), rows_size_(0), cols_size_(0), rows_offset_(0), cols_offset_(0)
{
  startup(label,nirreps,rows_size,cols_size);
}

BlockMatrix::~BlockMatrix()
{
  cleanup();
}

void BlockMatrix::startup(std::string label, int nirreps, size_t*& rows_size, size_t*& cols_size)
{
  label_   = label;
  nirreps_ = nirreps;

  // Allocate and compute the offsets
  allocate1(size_t,rows_size_,nirreps);
  allocate1(size_t,cols_size_,nirreps);
  for(int h = 0; h < nirreps; ++h){
    rows_size_[h] = rows_size[h];
    cols_size_[h] = cols_size[h];
  }

  // Allocate and compute the offsets
  allocate1(size_t,rows_offset_,nirreps);
  allocate1(size_t,cols_offset_,nirreps);
  rows_offset_[0] = 0;
  cols_offset_[0] = 0;
  for(int h = 1; h < nirreps; ++h){
    rows_offset_[h] = rows_offset_[h-1] + rows_size[h-1];
    cols_offset_[h] = cols_offset_[h-1] + cols_size[h-1];
  }

  // Allocate the blocks
  matrix_base_ = new MatrixBase*[nirreps_];
  for(int h = 0; h < nirreps_; ++h){
    matrix_base_[h] = new MatrixBase(rows_size_[h],cols_size_[h]);
  }
}

void BlockMatrix::startup(std::string label, int nirreps, vecint& rows_size, vecint& cols_size)
{
  label_   = label;
  nirreps_ = nirreps;

  // Allocate and compute the offsets
  allocate1(size_t,rows_size_,nirreps);
  allocate1(size_t,cols_size_,nirreps);
  for(int h = 0; h < nirreps; ++h){
    rows_size_[h] = rows_size[h];
    cols_size_[h] = cols_size[h];
  }

  // Allocate and compute the offsets
  allocate1(size_t,rows_offset_,nirreps);
  allocate1(size_t,cols_offset_,nirreps);
  rows_offset_[0] = 0;
  cols_offset_[0] = 0;
  for(int h = 1; h < nirreps; ++h){
    rows_offset_[h] = rows_offset_[h-1] + rows_size[h-1];
    cols_offset_[h] = cols_offset_[h-1] + cols_size[h-1];
  }

  // Allocate the blocks
  matrix_base_ = new MatrixBase*[nirreps_];
  for(int h = 0; h < nirreps_; ++h){
    matrix_base_[h] = new MatrixBase(rows_size_[h],cols_size_[h]);
  }
}

void BlockMatrix::startup(std::string label, int nirreps, int*& rows_size, int*& cols_size)
{
  label_   = label;
  nirreps_ = nirreps;

  // Allocate and compute the offsets
  allocate1(size_t,rows_size_,nirreps);
  allocate1(size_t,cols_size_,nirreps);
  for(int h = 0; h < nirreps; ++h){
    rows_size_[h] = rows_size[h];
    cols_size_[h] = cols_size[h];
  }

  // Allocate and compute the offsets
  allocate1(size_t,rows_offset_,nirreps);
  allocate1(size_t,cols_offset_,nirreps);
  rows_offset_[0] = 0;
  cols_offset_[0] = 0;
  for(int h = 1; h < nirreps; ++h){
    rows_offset_[h] = rows_offset_[h-1] + rows_size[h-1];
    cols_offset_[h] = cols_offset_[h-1] + cols_size[h-1];
  }

  // Allocate the blocks
  matrix_base_ = new MatrixBase*[nirreps_];
  for(int h = 0; h < nirreps_; ++h){
    matrix_base_[h] = new MatrixBase(rows_size_[h],cols_size_[h]);
  }
}

void BlockMatrix::cleanup()
{
  if(matrix_base_){
    for(int h = 0; h < nirreps_; ++h){
      delete matrix_base_[h];
    }
    delete[] matrix_base_;
    matrix_base_ = 0;
  }
  release1(rows_size_);
  release1(cols_size_);
  release1(rows_offset_);
  release1(cols_offset_);
}

void BlockMatrix::zero()
{
  for(int h = 0; h < nirreps_; ++h)
    matrix_base_[h]->zero();
}

void BlockMatrix::zero_diagonal()
{
  for(int h = 0; h < nirreps_; ++h)
    matrix_base_[h]->zero_diagonal();
}

void BlockMatrix::print()
{
  fprintf(outfile,"\n\n  ## %s ##\n",label_.c_str());
  for(int h = 0; h < nirreps_; ++h){
    fprintf(outfile,"\n[%zu*%zu]\n",rows_size_[h],cols_size_[h]);
    matrix_base_[h]->print();
  }
  fflush(outfile);
}

void BlockMatrix::scale(double factor)
{
  for(int h = 0; h < nirreps_; ++h)
    matrix_base_[h]->scale(factor);
}

void BlockMatrix::transpose()
{
  for(int h = 0; h < nirreps_; ++h)
    matrix_base_[h]->transpose();
}

void BlockMatrix::multiply(bool transpose_A, bool transpose_B, BlockMatrix* A, BlockMatrix* B)
{
  for(int h = 0; h < nirreps_; ++h)
    getMatrixBase(h)->multiply(transpose_A,         transpose_B,
                               A->getMatrixBase(h), B->getMatrixBase(h));
}

void BlockMatrix::diagonalize(BlockMatrix* eigenvectors,BlockVector* eigenvalues)
{
  for(int h = 0; h < nirreps_; ++h)
    getMatrixBase(h)->diagonalize(eigenvectors->getMatrixBase(h),
                            eigenvalues->getVectorBase(h));
}

double dot(BlockMatrix* A,BlockMatrix* B)
{
  double value = 0.0;
  for(int h = 0; h < A->nirreps_; ++h)
    value += dot(A->getMatrixBase(h),B->getMatrixBase(h));
  return(value);
}

BlockMatrix& BlockMatrix::operator=(BlockMatrix& rhs)
{
  if(this == &rhs){
    return(*this);
  }

  for(int h=0; h < nirreps_; ++h){
    if(rows_size_[h] * cols_size_[h]>0){
      for(int i = 0; i < rows_size_[h]; ++i)
        for(int j = 0; j < cols_size_[h]; ++j)
          matrix_base_[h]->set(i,j, rhs.matrix_base_[h]->get(i,j) );
    }
  }
  return(*this);
}

BlockMatrix& BlockMatrix::operator+=(const BlockMatrix& rhs)
{
  for(int h=0; h < nirreps_; ++h)
    *matrix_base_[h] += *rhs.matrix_base_[h];
  return(*this);
}

}}


// double operator^(const BlockMatrix& rhs,const BlockMatrix& lhs)
// {
//   double value = 0.0;
//   for(int h=0; h < nirreps_; ++h){
//     value += dot(rhs->getMatrixBase(h),lhs->getMatrixBase(h));
//   }
//   return(value);
// }


/*


BlockMatrix::BlockMatrix(std::string label_, int nirreps_, int*& block_size_)
: label(label_),nirreps(nirreps_),block_size(block_size_)
{
  // Compute the block_offset
  allocate1(double,block_offset,nirreps);
  block_offset[0] = 0;
  for(int h = 1; h < nirreps; ++h){
    block_offset[h] = block_offset[h-1] + block_size[h-1];
  }

  // Allocate the matrix
  allocate1(double**,matrix,nirreps);
  for(int h = 0; h < nirreps; ++h){
    allocate2(double,matrix[h],block_size[h],block_size[h]);
  }
}

BlockMatrix::~BlockMatrix()
{
  cleanup();
}

void BlockMatrix::cleanup()
{
  for(int h=0;h<nirreps;h++){
    release2(matrix[h]);
  }
  release1(matrix);
  release1(block_offset);
}

void BlockMatrix::minus(BlockMatrix* B)
{
  for(int h=0; h < nirreps; ++h){
    double** A_matrix_block = matrix[h];
    double** B_matrix_block = B->get_block(h);
    if(block_size[h]>0){
      for(int i = 0; i < block_size[h]; ++i)
        for(int j = 0; j < block_size[h]; ++j)
          A_matrix_block[i][j] -= B_matrix_block[i][j];
    }
  }
}






double operator^(const BlockMatrix& rhs,const BlockMatrix& lhs)
{
  double value = 0.0;
  int nirreps = rhs.get_nirreps();
  for(int h=0; h < nirreps; ++h){
    const double** rhs_matrix_block = rhs.get_block(h);
    const double** lhs_matrix_block = lhs.get_block(h);
    int block_size = rhs.get_block_size(h);
    if(block_size>0){
      for(int i = 0; i < block_size; ++i)
        for(int j = 0; j < block_size; ++j)
          value += lhs_matrix_block[i][j] * rhs_matrix_block[i][j];
    }
  }
  return(value);
}



*/

