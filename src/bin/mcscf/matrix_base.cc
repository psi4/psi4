#include <cstring>
#include <iostream>

#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libutil/libutil.h>

#include "matrix_base.h"

#include <psi4-dec.h>

extern FILE* outfile;

namespace psi{ namespace mcscf{

extern MemoryManager* memory_manager;

MatrixBase::MatrixBase(size_t rows, size_t cols) : rows_(rows),cols_(cols),elements_(rows*cols),matrix_(NULL)
{
  allocate2(double,matrix_,rows_,cols_);
}

MatrixBase::~MatrixBase()
{
  release2(matrix_);
}

void MatrixBase::print()
{
  for(size_t i=0 ; i < rows_; i++){
    fprintf(outfile,"\n  ");
    for(size_t j=0 ; j < cols_; j++)
      fprintf(outfile,"%10.6f",matrix_[i][j]);
  }
  fprintf(outfile,"\n");
}

void MatrixBase::scale(double factor)
{
  if(elements_>0)
    C_DSCAL(elements_,
            factor,
            &(matrix_[0][0]),
            1);
}

void MatrixBase::transpose()
{
  if(elements_>0){
    double temp;
    for(size_t i = 0; i < rows_; ++i){
      for(size_t j = i + 1; j < cols_; ++j){
        temp          = matrix_[i][j];
        matrix_[i][j] = matrix_[j][i];
        matrix_[j][i] = temp;
      }
    }
  }
}

void MatrixBase::zero()
{
  if(elements_>0)
    memset(&(matrix_[0][0]),'\0', sizeof(double) * elements_);
}

void MatrixBase::zero_diagonal()
{
  if(elements_>0 && (rows_ == cols_))
    for(size_t i=0 ; i < rows_; i++)
      matrix_[i][i] = 0.0;
}


void MatrixBase::multiply(bool transpose_A, bool transpose_B, MatrixBase* A, MatrixBase* B)
{
  char transa = transpose_A ? 't' : 'n';
  char transb = transpose_B ? 't' : 'n';
  if(elements_>0){
    // Multiply A and B
    size_t m = rows_;       // TODO This only works for square matrices!
    size_t n = rows_;
    size_t k = rows_;
    size_t nca = rows_;
    size_t ncb = rows_;
    size_t ncc = rows_;
    C_DGEMM(transa,
            transb,
            m,
            n,
            k,
            1.0,
            A->get_matrix()[0],nca,
            B->get_matrix()[0],ncb,
            0.0,
            get_matrix()[0],ncb);
  }
}

void MatrixBase::diagonalize(MatrixBase* eigenmatrix, VectorBase* eigenvalues)
{
  // Diagonalize the block
  if(elements_>0 && (rows_ == cols_)){
    sq_rsp(rows_,
           cols_,
           matrix_,
           eigenvalues->get_vector(),
           1,
           eigenmatrix->get_matrix(),
           1.0e-14);
  }
}

double dot(MatrixBase* A, MatrixBase* B)
{
  double value = 0.0;
  if(A->rows_ * A->cols_>0){
    for(size_t i = 0; i < A->rows_; ++i)
      for(size_t j = 0; j < A->cols_; ++j)
        value += A->matrix_[i][j] * B->matrix_[i][j];
  }
  return(value);
}

MatrixBase& MatrixBase::operator+=(const MatrixBase& rhs)
{
  if(elements_>0){
    for(size_t i = 0; i < rows_; ++i)
      for(size_t j = 0; j < cols_; ++j)
        matrix_[i][j] += rhs.matrix_[i][j];
  }
  return(*this);
}

MatrixBase& MatrixBase::operator-=(const MatrixBase& rhs)
{
  if(elements_>0){
    for(size_t i = 0; i < rows_; ++i)
      for(size_t j = 0; j < cols_; ++j)
        matrix_[i][j] -= rhs.matrix_[i][j];
  }
  return(*this);
}

}}
