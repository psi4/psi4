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

// Latest revision on April 38, 2013.
#include <stdio.h>
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libiwl/iwl.hpp"
#include "tensors.h"
#include "psi4/libparallel/ParallelPrinter.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"

using namespace std;

namespace psi{ namespace dfoccwave{


/********************************************************************************************/
/************************** 1d array ********************************************************/
/********************************************************************************************/
Tensor1d::Tensor1d(int d1)
{
  A1d_ = NULL;
  dim1_=d1;
  memalloc();
}//

Tensor1d::Tensor1d(string name, int d1)
{
  A1d_ = NULL;
  dim1_=d1;
  name_=name;
  memalloc();
}//

Tensor1d::Tensor1d()
{
  A1d_ = NULL;
  dim1_ = 0;

}//

Tensor1d::~Tensor1d()
{
  release();
}//

void Tensor1d::memalloc()
{
    if (A1d_) release();
    A1d_ = new double[dim1_];
    zero();
}//

void Tensor1d::release()
{
   if (!A1d_) return;
   delete [] A1d_;
   A1d_ = NULL;
}//

void Tensor1d::init(int d1)
{
    dim1_=d1;
    if (A1d_) release();
    A1d_ = new double[dim1_];
}//

void Tensor1d::init(string name, int d1)
{
    dim1_=d1;
    name_=name;
    if (A1d_) release();
    A1d_ = new double[dim1_];
}//

void Tensor1d::zero()
{
    memset(A1d_, 0, sizeof(double)*dim1_);
}//

void Tensor1d::print()
{
  if (name_.length()) outfile->Printf( "\n ## %s ##\n", name_.c_str());
  for (int p=0; p<dim1_; p++){
    outfile->Printf(" %3d %10.7f \n",p,A1d_[p]);
  }

}//

void Tensor1d::print(std::string OutFileRMR)
{
   std::shared_ptr<psi::PsiOutStream> printer=(OutFileRMR=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(OutFileRMR)));
   if (name_.length()) printer->Printf( "\n ## %s ##\n", name_.c_str());
  for (int p=0; p<dim1_; p++){
    printer->Printf(" %3d %10.7f \n",p,A1d_[p]);
  }
}//

void Tensor1d::set(int i, double value)
{
  A1d_[i]=value;
}//

void Tensor1d::set(double *vec)
{
    for (int i=0; i<dim1_; ++i) A1d_[i] = vec[i];
}//

void Tensor1d::set(const SharedTensor1d &vec)
{
    for (int i=0; i<dim1_; ++i) A1d_[i] = vec->A1d_[i];
}//

double Tensor1d::get(int i)
{
  return A1d_[i];
}//

void Tensor1d::add(const SharedTensor1d& a)
{
    /*
    double *lhs, *rhs;
    size_t size = dim1_;
    if (size) {
        lhs = A1d_;
        rhs = Adum->A1d_;
        for (size_t ij=0; ij<size; ++ij) {
            *lhs += *rhs;
            lhs++; rhs++;
        }
    }
    */
      #pragma omp parallel for
      for (int i=0; i<dim1_; ++i) A1d_[i] += a->A1d_[i];

}//

void Tensor1d::add(int i, double value)
{
  A1d_[i]+=value;
}//

void Tensor1d::subtract(const SharedTensor1d& a)
{
    /*
    double *lhs, *rhs;
    size_t size = dim1_;
    if (size) {
        lhs = A1d_;
        rhs = Adum->A1d_;
        for (size_t ij=0; ij<size; ++ij) {
            *lhs -= *rhs;
            lhs++; rhs++;
        }
    }
    */
      #pragma omp parallel for
      for (int i=0; i<dim1_; ++i) A1d_[i] -= a->A1d_[i];
}//

void Tensor1d::subtract(int i, double value)
{
  A1d_[i]-=value;
}//

void Tensor1d::to_shared_vector(SharedVector A)
{
      #pragma omp parallel for
      for (int i=0; i<dim1_; ++i) {
           A->set(0,i,A1d_[i]);
      }
}//

double Tensor1d::rms()
{
  double summ = 0.0;
  for (int i=0; i<dim1_; ++i) summ += A1d_[i] * A1d_[i];
  summ=sqrt(summ/dim1_);

  return summ;
}//

double Tensor1d::rms(const SharedTensor1d& Atemp)
{
  double summ = 0.0;
  for (int i=0; i<dim1_; ++i) summ += (A1d_[i] - Atemp->A1d_[i])  * (A1d_[i] - Atemp->A1d_[i]);
  summ=sqrt(summ/dim1_);

  return summ;
}//

double Tensor1d::dot(const SharedTensor1d &y)
{
  double value = 0.0;
  int incx = 1;
  int incy = 1;
  if (dim1_ == y->dim1_) value = C_DDOT((ULI)dim1_, A1d_, incx, y->A1d_, incy);
  return value;
}//

void Tensor1d::gbmv(bool transa, const SharedTensor2d& a, const SharedTensor1d& b, double alpha, double beta)
{
    char ta = transa ? 't' : 'n';
    int m, n, k, kl, ku, incx, incy, lda;

    m = a->dim1_;
    n = a->dim2_;
    k = b->dim1_;
    kl = m - 1;// # of subdiagonal of matrix A, at most kl = m - 1
    ku = n - 1;// # of superdiagonal of matrix A, at most ku = n - 1
    lda = kl + ku + 1;
    incx = 1;// increments in elements of b vector
    incy = 1;// increments in elements of A1d_

    // A1d_ = alpha * A * b + beta, where A is a band matrix
    if (m && n) {
       C_DGBMV(ta, m, n, kl, ku, alpha, &(a->A2d_[0][0]), lda, b->A1d_, incx, beta, A1d_, incy);
    }
}//

void Tensor1d::gemv(bool transa, const SharedTensor2d& a, const SharedTensor1d& b, double alpha, double beta)
{
    char ta = transa ? 't' : 'n';
    int m, n, k, incx, incy, lda;

    m = a->dim1();
    n = a->dim2();
    lda = n;
    incx = 1;// increments in elements of b vector
    incy = 1;// increments in elements of A1d_

    // A1d_ = alpha * A * b + beta, where A is a general matrix
    if (m && n) {
       C_DGEMV(ta, m, n, alpha, &(a->A2d_[0][0]), lda, b->A1d_, incx, beta, A1d_, incy);
    }
}//

void Tensor1d::gemv(bool transa, int m, int n, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta)
{
    char ta = transa ? 't' : 'n';
    int incx, incy, lda;

    lda = n;
    incx = 1;// increments in elements of b vector
    incy = 1;// increments in elements of A1d_

    // A1d_ = alpha * A * b + beta, where A is a general matrix
    if (m && n) {
       C_DGEMV(ta, m, n, alpha, &(a->A2d_[0][0]), lda, &(b->A2d_[0][0]), incx, beta, A1d_, incy);
    }
}//

void Tensor1d::gemv(bool transa, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta)
{
    char ta = transa ? 't' : 'n';
    int incx, incy, lda, m, n;

    m = a->dim1();
    n = a->dim2();

    /*
    if (transa) {
        m = a->dim1();
        n = a->dim2();
    }

    else {
        n = a->dim1();
        m = a->dim2();
    }
    */

    lda = n;
    incx = 1;// increments in elements of b vector
    incy = 1;// increments in elements of A1d_

    // A1d_ = alpha * A * b + beta, where A is a general matrix
    if (m && n) {
       C_DGEMV(ta, m, n, alpha, &(a->A2d_[0][0]), lda, &(b->A2d_[0][0]), incx, beta, A1d_, incy);
    }
}//

void Tensor1d::gemv(bool transa, int m, int n, const SharedTensor2d& a, const SharedTensor2d& b, int start_a, int start_b, double alpha, double beta)
{
    char ta = transa ? 't' : 'n';
    int incx, incy, lda;

    lda = n;
    incx = 1;// increments in elements of b vector
    incy = 1;// increments in elements of A1d_

    // A1d_ = alpha * A * b + beta, where A is a general matrix
    if (m && n) {
       C_DGEMV(ta, m, n, alpha, a->A2d_[0]+start_a, lda, b->A2d_[0]+start_b, incx, beta, A1d_, incy);
    }
}//

void Tensor1d::gemv(bool transa, int m, int n, const SharedTensor2d& a, const SharedTensor2d& b, int start_a, int start_b, int start_c, double alpha, double beta)
{
    char ta = transa ? 't' : 'n';
    int incx, incy, lda;

    lda = n;
    incx = 1;// increments in elements of b vector
    incy = 1;// increments in elements of A1d_

    // A1d_ = alpha * A * b + beta, where A is a general matrix
    if (m && n) {
       C_DGEMV(ta, m, n, alpha, a->A2d_[0]+start_a, lda, b->A2d_[0]+start_b, incx, beta, A1d_+start_c, incy);
    }
}//

double Tensor1d::xay(const SharedTensor2d &a, const SharedTensor1d &y)
{
  double value = 0.0;
  SharedTensor1d ay = SharedTensor1d(new Tensor1d(a->dim1_));
  ay->gemv(false, a, y, 1.0, 0.0);
  value = dot(ay);
  return value;
}//

void Tensor1d::axpy(const SharedTensor1d &a, double alpha)
{
    ULI length = (ULI)dim1_;
    C_DAXPY(length, alpha, a->A1d_, 1, A1d_, 1);
}

void Tensor1d::scale(double a)
{
    //size_t size = dim1_ ;
    ULI size = (ULI)dim1_ ;
    if (size) C_DSCAL(size, a, A1d_, 1);
}//

void Tensor1d::copy(double *a)
{
    //size_t size;
    //size = dim1_ * sizeof(double);
    //if (size) memcpy(&(A1d_[0]), &(x[0]), size);
    ULI size = (ULI)dim1_ ;
    C_DCOPY(size, a, 1, A1d_, 1);
}//

void Tensor1d::copy(const SharedTensor1d &a)
{
    //size_t size;
    //size = dim1_ * sizeof(double);
    //if (size) memcpy(&(A1d_[0]), &(x->A1d_[0]), size);
    ULI size = (ULI)dim1_ ;
    C_DCOPY(size, a->A1d_, 1, A1d_, 1);
}//

void Tensor1d::row_vector(SharedTensor2d &A, int n)
{
  int dim = A->dim2();
  for (int i=0; i<dim; i++) A1d_[i] = A->get(n, i);
}//

void Tensor1d::column_vector(SharedTensor2d &A, int n)
{
  int dim = A->dim1();
  for (int i=0; i<dim; i++) A1d_[i] = A->get(i, n);
}//

void Tensor1d::dirprd(SharedTensor1d &a, SharedTensor1d &b)
{
  int dima = a->dim1();
  int dimb = b->dim1();

  if (dima == dimb && dima == dim1_) {
      for (int i=0; i<dim1_; i++) A1d_[i] = a->get(i) * b->get(i);
  }
  else throw SanityCheckError("Vector dimensions do NOT match!", __FILE__, __LINE__);
}//

void Tensor1d::symm_packed(const SharedTensor2d &A)
{
    // Form Lower triangular part
    #pragma omp parallel for
    for (int p = 0; p < A->dim1(); p++) {
          for (int q = 0; q <= p; q++) {
               int pq = index2(p,q);
               double perm = (p == q ? 1.0 : 2.0);
               A1d_[pq] = perm * A->get(p,q);
          }
    }

}//

void Tensor1d::ltm(const SharedTensor2d &A)
{
    // Form Lower triangular part
    #pragma omp parallel for
    for (int p = 0; p < A->dim1(); p++) {
          for (int q = 0; q <= p; q++) {
               int pq = index2(p,q);
               A1d_[pq] = A->get(p,q);
          }
    }

}//

/********************************************************************************************/
/************************** 2d array ********************************************************/
/********************************************************************************************/
Tensor2d::Tensor2d(int d1,int d2)
{
  A2d_ = NULL;
  row_idx_ = NULL;
  col_idx_ = NULL;
  row2d1_ = NULL;
  row2d2_ = NULL;
  col2d1_ = NULL;
  col2d2_ = NULL;
  d1_ = 0;
  d2_ = 0;
  dim1_=d1;
  dim2_=d2;
  memalloc();
}//

Tensor2d::Tensor2d(string name, int d1,int d2)
{
  A2d_ = NULL;
  row_idx_ = NULL;
  col_idx_ = NULL;
  row2d1_ = NULL;
  row2d2_ = NULL;
  col2d1_ = NULL;
  col2d2_ = NULL;
  d1_ = 0;
  d2_ = 0;
  dim1_=d1;
  dim2_=d2;
  name_=name;
  memalloc();
}//

Tensor2d::Tensor2d()
{
  A2d_ = NULL;
  row_idx_ = NULL;
  col_idx_ = NULL;
  row2d1_ = NULL;
  row2d2_ = NULL;
  col2d1_ = NULL;
  col2d2_ = NULL;
  d1_ = 0;
  d2_ = 0;
  dim1_ = 0;
  dim2_ = 0;

}//

Tensor2d::Tensor2d(psi::PSIO* psio, unsigned int fileno, string name, int d1,int d2)
{
  A2d_ = NULL;
  row_idx_ = NULL;
  col_idx_ = NULL;
  row2d1_ = NULL;
  row2d2_ = NULL;
  col2d1_ = NULL;
  col2d2_ = NULL;
  d1_ = 0;
  d2_ = 0;
  dim1_=d1;
  dim2_=d2;
  name_=name;
  memalloc();
  read(psio, fileno);
}

Tensor2d::Tensor2d(std::shared_ptr<psi::PSIO> psio, unsigned int fileno, string name, int d1,int d2)
{
  A2d_ = NULL;
  row_idx_ = NULL;
  col_idx_ = NULL;
  row2d1_ = NULL;
  row2d2_ = NULL;
  col2d1_ = NULL;
  col2d2_ = NULL;
  d1_ = 0;
  d2_ = 0;
  dim1_=d1;
  dim2_=d2;
  name_=name;
  memalloc();
  read(psio, fileno);
}

Tensor2d::Tensor2d(psi::PSIO& psio, unsigned int fileno, string name, int d1,int d2)
{
  A2d_ = NULL;
  row_idx_ = NULL;
  col_idx_ = NULL;
  row2d1_ = NULL;
  row2d2_ = NULL;
  col2d1_ = NULL;
  col2d2_ = NULL;
  d1_ = 0;
  d2_ = 0;
  dim1_=d1;
  dim2_=d2;
  name_=name;
  memalloc();
  read(&psio, fileno);
}//

Tensor2d::Tensor2d(string name, int d1, int d2, int d3, int d4)
{
  A2d_ = NULL;
  row_idx_ = NULL;
  col_idx_ = NULL;
  row2d1_ = NULL;
  row2d2_ = NULL;
  col2d1_ = NULL;
  col2d2_ = NULL;
  d1_ = d1;
  d2_ = d2;
  d3_ = d3;
  d4_ = d4;
  dim1_=d1*d2;
  dim2_=d3*d4;
  name_=name;

  // memalloc
  if (A2d_) release();
  A2d_ = block_matrix(dim1_, dim2_);
  zero();

  // row idx
  row_idx_ = init_int_matrix(d1_, d2_);
  memset(row_idx_[0], 0, sizeof(int)*d1_*d2_);
  row2d1_ = new int[dim1_];
  row2d2_ = new int[dim1_];
  memset(row2d1_, 0, sizeof(int)*dim1_);
  memset(row2d2_, 0, sizeof(int)*dim1_);
  for (int i = 0; i < d1_; i++) {
       for (int a = 0; a < d2_; a++) {
            int ia = a + (i * d2_);
            row_idx_[i][a] = ia;
            row2d1_[ia] = i;
            row2d2_[ia] = a;
       }
  }

  // col idx
  col_idx_ = init_int_matrix(d3_, d4_);
  memset(col_idx_[0], 0, sizeof(int)*d3_*d4_);
  col2d1_ = new int[dim2_];
  col2d2_ = new int[dim2_];
  memset(col2d1_, 0, sizeof(int)*dim2_);
  memset(col2d2_, 0, sizeof(int)*dim2_);
  for (int i = 0; i < d3_; i++) {
       for (int a = 0; a < d4_; a++) {
            int ia = a + (i * d4_);
            col_idx_[i][a] = ia;
            col2d1_[ia] = i;
            col2d2_[ia] = a;
       }
  }

}//

Tensor2d::Tensor2d(string name, int d1, int d2, int d3)
{
  A2d_ = NULL;
  row_idx_ = NULL;
  col_idx_ = NULL;
  row2d1_ = NULL;
  row2d2_ = NULL;
  col2d1_ = NULL;
  col2d2_ = NULL;
  d1_ = d1;
  d2_ = d2;
  d3_ = d3;
  d4_ = 0;
  dim1_=d1;
  dim2_=d2*d3;
  name_=name;

  // memalloc
  if (A2d_) release();
  A2d_ = block_matrix(dim1_, dim2_);
  zero();

  // col idx
  col_idx_ = init_int_matrix(d2_, d3_);
  memset(col_idx_[0], 0, sizeof(int)*d2_*d3_);
  col2d1_ = new int[dim2_];
  col2d2_ = new int[dim2_];
  memset(col2d1_, 0, sizeof(int)*dim2_);
  memset(col2d2_, 0, sizeof(int)*dim2_);
  for (int i = 0; i < d2_; i++) {
       for (int a = 0; a < d3_; a++) {
            int ia = a + (i * d3_);
            col_idx_[i][a] = ia;
            col2d1_[ia] = i;
            col2d2_[ia] = a;
       }
  }

}//

Tensor2d::~Tensor2d()
{
  release();
}//

void Tensor2d::memalloc()
{
    if (A2d_) release();
    A2d_ = block_matrix(dim1_, dim2_);
    zero();
}//

void Tensor2d::release()
{
   //if (!A2d_) return;
   //free_block(A2d_);
   if (A2d_) free_block(A2d_);
   if (row_idx_) free_int_matrix(row_idx_);
   if (col_idx_) free_int_matrix(col_idx_);
   if (row2d1_) delete [] row2d1_;
   if (row2d2_) delete [] row2d2_;
   if (col2d1_) delete [] col2d1_;
   if (col2d2_) delete [] col2d2_;

   A2d_ = NULL;
   row_idx_ = NULL;
   col_idx_ = NULL;
   row2d1_ = NULL;
   row2d2_ = NULL;
   col2d1_ = NULL;
   col2d2_ = NULL;
}//

void Tensor2d::init(int d1,int d2)
{
    dim1_=d1;
    dim2_=d2;
    if (A2d_) release();
    A2d_ = block_matrix(dim1_, dim2_);
}//

void Tensor2d::init(string name, int d1,int d2)
{
    dim1_=d1;
    dim2_=d2;
    name_=name;
    if (A2d_) release();
    A2d_ = block_matrix(dim1_, dim2_);
}//

void Tensor2d::zero()
{
    memset(A2d_[0], 0, sizeof(double)*dim1_*dim2_);
}//

void Tensor2d::zero_diagonal()
{
    if (dim1_ == dim2_ ) {
       for (int i=0; i<dim1_; i++) A2d_[i][i] = 0.0;
   }
}//

void Tensor2d::print()
{
  if (A2d_) {
      if (name_.length()) outfile->Printf( "\n ## %s ##\n", name_.c_str());
      print_mat(A2d_,dim1_,dim2_,"outfile");

  }
}//

void Tensor2d::print(std::string OutFileRMR)
{
   std::shared_ptr<psi::PsiOutStream> printer=(OutFileRMR=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(OutFileRMR)));
   if (A2d_) {
      if (name_.length()) printer->Printf( "\n ## %s ##\n", name_.c_str());
      print_mat(A2d_,dim1_,dim2_,OutFileRMR);
  }
}//

void Tensor2d::set(int i, int j, double value)
{
  A2d_[i][j]=value;
}//

void Tensor2d::set(double **A)
{
      if (A == NULL) return;
      #pragma omp parallel for
      for (int i=0; i<dim1_; ++i) {
        for (int j=0; j<dim2_; ++j) {
          A2d_[i][j] = A[i][j];
        }
      }
}//

void Tensor2d::set(SharedTensor2d &A)
{
      if (A == NULL) return;
      #pragma omp parallel for
      for (int i=0; i<dim1_; ++i) {
        for (int j=0; j<dim2_; ++j) {
          A2d_[i][j] = A->A2d_[i][j];
        }
      }
}//

void Tensor2d::set(SharedMatrix A)
{
      #pragma omp parallel for
      for (int i=0; i<dim1_; ++i) {
        for (int j=0; j<dim2_; ++j) {
          A2d_[i][j] = A->get(0,i,j);
        }
      }
}//

void Tensor2d::set2(SharedMatrix A)
{
      #pragma omp parallel for
      for (int i=0; i<dim1_; ++i) {
        for (int j=0; j<dim2_; ++j) {
          A2d_[i][j] = A->get(i,j);
        }
      }
}//

void Tensor2d::set(SharedTensor1d &A)
{
      #pragma omp parallel for
      for (int i=0; i<dim1_; i++) {
        for (int j=0; j<dim2_; j++) {
             int ij = j + (i*dim2_);
             A2d_[i][j] = A->get(ij);
        }
      }
}//

void Tensor2d::set(double *A)
{
      #pragma omp parallel for
      for (int i=0; i<dim1_; i++) {
        for (int j=0; j<dim2_; j++) {
             int ij = j + (i*dim2_);
             A2d_[i][j] = A[ij];
        }
      }
}//

double Tensor2d::get(int i, int j)
{
  return A2d_[i][j];
}//

void Tensor2d::gemm(bool transa, bool transb, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta)
{
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int m, n, k, nca, ncb, ncc;

    m = dim1_;
    n = dim2_;
    k = transa ? a->dim1_ : a->dim2_;
    nca = transa ? m : k; // lda
    ncb = transb ? k : n; // ldb
    ncc = n; // ldc

    if (m && n && k) {
        C_DGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), nca, &(b->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
    }
}//

void Tensor2d::contract(bool transa, bool transb, int m, int n, int k, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta)
{
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int lda, ldb, ldc;

    lda = transa ? m : k;
    ldb = transb ? k : n;
    ldc = n;

    if (m && n && k) {
        C_DGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), lda, &(b->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
    }
}//

void Tensor2d::contract(bool transa, bool transb, int m, int n, int k, const SharedTensor2d& a, const SharedTensor2d& b, int start_a, int start_b, double alpha, double beta)
{
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int lda, ldb, ldc;

    lda = transa ? m : k;
    ldb = transb ? k : n;
    ldc = n;

    if (m && n && k) {
        C_DGEMM(ta, tb, m, n, k, alpha, a->A2d_[0]+start_a, lda, b->A2d_[0]+start_b, ldb, beta, A2d_[0], ldc);
    }
}//

void Tensor2d::contract(bool transa, bool transb, int m, int n, int k, const SharedTensor2d& a, const SharedTensor2d& b,
                        int start_a, int start_b, int start_c, double alpha, double beta)
{
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int lda, ldb, ldc;

    lda = transa ? m : k;
    ldb = transb ? k : n;
    ldc = n;

    if (m && n && k) {
        C_DGEMM(ta, tb, m, n, k, alpha, a->A2d_[0]+start_a, lda, b->A2d_[0]+start_b, ldb, beta, A2d_[0]+start_c, ldc);
    }
}//

void Tensor2d::contract323(bool transa, bool transb, int m, int n, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta)
{
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int k, nca, ncb, ncc;

    k = transb ? b->dim2_ : b->dim1_;
    nca = transa ? m : k;
    ncb = transb ? k : n;
    ncc = n;

    if (m && n && k) {
        #pragma omp parallel for
        for (int Q = 0; Q < dim1_; Q++) {
             C_DGEMM(ta, tb, m, n, k, alpha, a->A2d_[Q], nca, b->A2d_[0], ncb, beta, A2d_[Q], ncc);
        }
    }
}//

void Tensor2d::contract233(bool transa, bool transb, int m, int n, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta)
{
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int k, lda, ldb, ldc;

    k = transa ? a->dim1_ : a->dim2_;
    lda = transa ? m : k;
    ldb = transb ? k : n;
    ldc = n;

    if (m && n && k) {
        #pragma omp parallel for
        for (int Q = 0; Q < dim1_; Q++) {
             C_DGEMM(ta, tb, m, n, k, alpha, a->A2d_[0], lda, b->A2d_[Q], ldb, beta, A2d_[Q], ldc);
        }
    }
}//

void Tensor2d::contract332(bool transa, bool transb, int k, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta)
{
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int m, n, nca, ncb, ncc;

    m = dim1_;
    n = dim2_;
    nca = transa ? m : k;
    ncb = transb ? k : n;
    ncc = n;

    if (m && n && k) {
        //#pragma omp parallel for
        for (int Q = 0; Q < a->dim1(); Q++) {
             C_DGEMM(ta, tb, m, n, k, alpha, a->A2d_[Q], nca, b->A2d_[Q], ncb, beta, A2d_[0], ncc);
        }
    }
}//

void Tensor2d::contract424(int target_x, int target_y, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta)
{

    char ta;
    char tb;
    int lda, ldb, ldc;
    int m, n, k;

    // C(pq,rs) = \sum_{o} A(oq,rs) B(o,p)
    if (target_x == 1 && target_y == 1) {
        ta = 't';
        tb = 'n';
        m = d1_;
        n = d2_*d3_*d4_;
        k = b->dim1();
        lda = m;
        ldb = n;
        ldc = n;

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, b->A2d_[0], lda, a->A2d_[0], ldb, beta, A2d_[0], ldc);
        }
    }

    // C(pq,rs) = \sum_{o} A(oq,rs) B(p,o)
    else if (target_x == 1 && target_y == 2) {
        ta = 'n';
        tb = 'n';
        m = d1_;
        n = d2_*d3_*d4_;
        k = b->dim2();
        lda = k;
        ldb = n;
        ldc = n;

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, b->A2d_[0], lda, a->A2d_[0], ldb, beta, A2d_[0], ldc);
        }
    }

    // C(pq,rs) = \sum_{o} A(po,rs) B(o,q)
    else if (target_x == 2 && target_y == 1) {
        ta = 'n';
        tb = 'n';
        m = d1_*d3_*d4_;
        n = d2_;
        k = b->dim1();
        lda = k;
        ldb = n;
        ldc = n;

        SharedTensor2d temp1 = SharedTensor2d(new Tensor2d("temp1", a->d1_, a->d3_, a->d4_, a->d2_));
        SharedTensor2d temp2 = SharedTensor2d(new Tensor2d("temp2", d1_, d3_, d4_, d2_));
        temp1->sort(1342, a, 1.0, 0.0);

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, temp1->A2d_[0], lda, b->A2d_[0], ldb, 0.0, temp2->A2d_[0], ldc);
        }
        temp1.reset();
        SharedTensor2d temp3 = SharedTensor2d(new Tensor2d("temp3", d1_, d2_, d3_, d4_));
        temp3->sort(1423, temp2, 1.0, 0.0);
        temp2.reset();
        scale(beta);
        add(temp3);
        temp3.reset();
    }

    // C(pq,rs) = \sum_{o} A(po,rs) B(q,o)
    else if (target_x == 2 && target_y == 2) {
        ta = 'n';
        tb = 't';
        m = d1_*d3_*d4_;
        n = d2_;
        k = b->dim2();
        lda = k;
        ldb = k;
        ldc = n;

        SharedTensor2d temp1 = SharedTensor2d(new Tensor2d("temp1", a->d1_, a->d3_, a->d4_, a->d2_));
        SharedTensor2d temp2 = SharedTensor2d(new Tensor2d("temp2", d1_, d3_, d4_, d2_));
        temp1->sort(1342, a, 1.0, 0.0);

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, temp1->A2d_[0], lda, b->A2d_[0], ldb, 0.0, temp2->A2d_[0], ldc);
        }
        temp1.reset();
        SharedTensor2d temp3 = SharedTensor2d(new Tensor2d("temp3", d1_, d2_, d3_, d4_));
        temp3->sort(1423, temp2, 1.0, 0.0);
        temp2.reset();
        scale(beta);
        add(temp3);
        temp3.reset();
    }

    // C(pq,rs) = \sum_{o} A(pq,os) B(o,r)
    else if (target_x == 3 && target_y == 1) {
        ta = 'n';
        tb = 'n';
        m = d1_*d2_*d4_;
        n = d3_;
        k = b->dim1();
        lda = k;
        ldb = n;
        ldc = n;

        SharedTensor2d temp1 = SharedTensor2d(new Tensor2d("temp1", a->d1_, a->d2_, a->d4_, a->d3_));
        SharedTensor2d temp2 = SharedTensor2d(new Tensor2d("temp2", d1_, d2_, d4_, d3_));
        temp1->sort(1243, a, 1.0, 0.0);

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, temp1->A2d_[0], lda, b->A2d_[0], ldb, 0.0, temp2->A2d_[0], ldc);
        }
        temp1.reset();
        SharedTensor2d temp3 = SharedTensor2d(new Tensor2d("temp3", d1_, d2_, d3_, d4_));
        temp3->sort(1243, temp2, 1.0, 0.0);
        temp2.reset();
        scale(beta);
        add(temp3);
        temp3.reset();
    }

    // C(pq,rs) = \sum_{o} A(pq,os) B(r,o)
    else if (target_x == 3 && target_y == 2) {
        ta = 'n';
        tb = 't';
        m = d1_*d2_*d4_;
        n = d3_;
        k = b->dim2();
        lda = k;
        ldb = k;
        ldc = n;

        SharedTensor2d temp1 = SharedTensor2d(new Tensor2d("temp1", a->d1_, a->d2_, a->d4_, a->d3_));
        SharedTensor2d temp2 = SharedTensor2d(new Tensor2d("temp2", d1_, d2_, d4_, d3_));
        temp1->sort(1243, a, 1.0, 0.0);

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, temp1->A2d_[0], lda, b->A2d_[0], ldb, 0.0, temp2->A2d_[0], ldc);
        }
        temp1.reset();
        SharedTensor2d temp3 = SharedTensor2d(new Tensor2d("temp3", d1_, d2_, d3_, d4_));
        temp3->sort(1243, temp2, 1.0, 0.0);
        temp2.reset();
        scale(beta);
        add(temp3);
        temp3.reset();
    }

    // C(pq,rs) = \sum_{o} A(pq,ro) B(o,s)
    else if (target_x == 4 && target_y == 1) {
        ta = 'n';
        tb = 'n';
        m = d1_ * d2_ * d3_;
        n = d4_;
        k = b->dim1();
        lda = k;
        ldb = n;
        ldc = n;
        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), lda, &(b->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
        }
    }

    // C(pq,rs) = \sum_{o} A(pq,ro) B(s,o)
    else if (target_x == 4 && target_y == 2) {
        ta = 'n';
        tb = 't';
        m = d1_ * d2_ * d3_;
        n = d4_;
        k = b->dim2();
        lda = k;
        ldb = k;
        ldc = n;
        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), lda, &(b->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
        }
    }

    /*
    // C(pq,rs) = \sum_{o} A(pq,ro) B(o,s)
    else if (target_x == 4 && target_y == 1) {
    #pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
         for (int q = 0; q < d2_; q++) {
              int pq = row_idx_[p][q];
              for (int r = 0; r < d3_; r++) {
                   for (int s = 0; s < d4_; s++) {
                        int rs = col_idx_[r][s];
                        double sum = 0.0;
                        for (int o = 0; o < b->dim1(); o++) {
                             int ro = a->col_idx_[r][o];
                             sum += a->get(pq,ro) * b->get(o,s);
                        }
                        A2d_[pq][rs] = (alpha * sum) + (beta * A2d_[pq][rs]);
                   }
              }
         }
    }
    }

    // C(pq,rs) = \sum_{o} A(pq,ro) B(s,o)
    else if (target_x == 4 && target_y == 2) {
    #pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
         for (int q = 0; q < d2_; q++) {
              int pq = row_idx_[p][q];
              for (int r = 0; r < d3_; r++) {
                   for (int s = 0; s < d4_; s++) {
                        int rs = col_idx_[r][s];
                        double sum = 0.0;
                        for (int o = 0; o < b->dim2(); o++) {
                             int ro = a->col_idx_[r][o];
                             sum += a->get(pq,ro) * b->get(s,o);
                        }
                        A2d_[pq][rs] = (alpha * sum) + (beta * A2d_[pq][rs]);
                   }
              }
         }
    }
    }

    // C(pq,rs) = \sum_{o} A(oq,rs) B(o,p)
    if (target_x == 1 && target_y == 1) {
    #pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
         for (int q = 0; q < d2_; q++) {
              int pq = row_idx_[p][q];
              for (int r = 0; r < d3_; r++) {
                   for (int s = 0; s < d4_; s++) {
                        int rs = col_idx_[r][s];
                        double sum = 0.0;
                        for (int o = 0; o < b->dim1(); o++) {
                             int oq = a->row_idx_[o][q];
                             sum += a->get(oq,rs) * b->get(o,p);
                        }
                        A2d_[pq][rs] = (alpha * sum) + (beta * A2d_[pq][rs]);
                   }
              }
         }
    }
    }

    // C(pq,rs) = \sum_{o} A(oq,rs) B(p,o)
    else if (target_x == 1 && target_y == 2) {
    #pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
         for (int q = 0; q < d2_; q++) {
              int pq = row_idx_[p][q];
              for (int r = 0; r < d3_; r++) {
                   for (int s = 0; s < d4_; s++) {
                        int rs = col_idx_[r][s];
                        double sum = 0.0;
                        for (int o = 0; o < b->dim2(); o++) {
                             int oq = a->row_idx_[o][q];
                             sum += a->get(oq,rs) * b->get(p,o);
                        }
                        A2d_[pq][rs] = (alpha * sum) + (beta * A2d_[pq][rs]);
                   }
              }
         }
    }
    }

    // C(pq,rs) = \sum_{o} A(po,rs) B(o,q)
    else if (target_x == 2 && target_y == 1) {
    #pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
         for (int q = 0; q < d2_; q++) {
              int pq = row_idx_[p][q];
              for (int r = 0; r < d3_; r++) {
                   for (int s = 0; s < d4_; s++) {
                        int rs = col_idx_[r][s];
                        double sum = 0.0;
                        for (int o = 0; o < b->dim1(); o++) {
                             int po = a->row_idx_[p][o];
                             sum += a->get(po,rs) * b->get(o,q);
                        }
                        A2d_[pq][rs] = (alpha * sum) + (beta * A2d_[pq][rs]);
                   }
              }
         }
    }
    }

    // C(pq,rs) = \sum_{o} A(po,rs) B(q,o)
    else if (target_x == 2 && target_y == 2) {
    #pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
         for (int q = 0; q < d2_; q++) {
              int pq = row_idx_[p][q];
              for (int r = 0; r < d3_; r++) {
                   for (int s = 0; s < d4_; s++) {
                        int rs = col_idx_[r][s];
                        double sum = 0.0;
                        for (int o = 0; o < b->dim2(); o++) {
                             int po = a->row_idx_[p][o];
                             sum += a->get(po,rs) * b->get(q,o);
                        }
                        A2d_[pq][rs] = (alpha * sum) + (beta * A2d_[pq][rs]);
                   }
              }
         }
    }
    }

    // C(pq,rs) = \sum_{o} A(pq,os) B(o,r)
    else if (target_x == 3 && target_y == 1) {
    #pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
         for (int q = 0; q < d2_; q++) {
              int pq = row_idx_[p][q];
              for (int r = 0; r < d3_; r++) {
                   for (int s = 0; s < d4_; s++) {
                        int rs = col_idx_[r][s];
                        double sum = 0.0;
                        for (int o = 0; o < b->dim1(); o++) {
                             int os = a->col_idx_[o][s];
                             sum += a->get(pq,os) * b->get(o,r);
                        }
                        A2d_[pq][rs] = (alpha * sum) + (beta * A2d_[pq][rs]);
                   }
              }
         }
    }
    }

    // C(pq,rs) = \sum_{o} A(pq,os) B(r,o)
    else if (target_x == 3 && target_y == 2) {
    #pragma omp parallel for
    for (int p = 0; p < d1_; p++) {
         for (int q = 0; q < d2_; q++) {
              int pq = row_idx_[p][q];
              for (int r = 0; r < d3_; r++) {
                   for (int s = 0; s < d4_; s++) {
                        int rs = col_idx_[r][s];
                        double sum = 0.0;
                        for (int o = 0; o < b->dim2(); o++) {
                             int os = a->col_idx_[o][s];
                             sum += a->get(pq,os) * b->get(r,o);
                        }
                        A2d_[pq][rs] = (alpha * sum) + (beta * A2d_[pq][rs]);
                   }
              }
         }
    }
    }
    */

    else {
         outfile->Printf("\tcontract424: Unrecognized targets! \n");

    }

}//

void Tensor2d::contract442(int target_a, int target_b, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta)
{

    char ta;
    char tb;
    int lda, ldb, ldc;
    int m, n, k;

    // C(p,q) = \sum_{rst} A(pr,st) B(qr,st)
    if (target_a == 1 && target_b == 1) {
        ta = 'n';
        tb = 't';
        m = dim1_;
        n = dim2_;
        k = a->d2_ * a->d3_ * a->d4_;
        lda = k;
        ldb = k;
        ldc = n;

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), lda, &(b->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
        }
    }

    // C(p,q) = \sum_{rst} A(pr,st) B(rq,st)
    else if (target_a == 1 && target_b == 2) {
        ta = 'n';
        tb = 't';
        m = dim1_;
        n = dim2_;
        k = a->d2_ * a->d3_ * a->d4_;
        lda = k;
        ldb = k;
        ldc = n;

        SharedTensor2d temp = SharedTensor2d(new Tensor2d("temp", b->d2_, b->d1_, b->d3_, b->d4_));
        temp->sort(2134, b, 1.0, 0.0);
        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), lda, &(temp->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
        }
        temp.reset();
    }

    // C(p,q) = \sum_{rst} A(pr,st) B(rs,qt)
    else if (target_a == 1 && target_b == 3) {
        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = a->d2_ * a->d3_ * a->d4_;
        lda = k;
        ldb = n;
        ldc = n;

        SharedTensor2d temp = SharedTensor2d(new Tensor2d("temp", b->d1_, b->d2_, b->d4_, b->d3_));
        temp->sort(1243, b, 1.0, 0.0);
        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), lda, &(temp->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
        }
        temp.reset();
    }

    // C(p,q) = \sum_{rst} A(pr,st) B(rs,tq)
    else if (target_a == 1 && target_b == 4) {
        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = a->d2_ * a->d3_ * a->d4_;
        lda = k;
        ldb = n;
        ldc = n;
        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, &(a->A2d_[0][0]), lda, &(b->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
        }
    }

    // C(p,q) = \sum_{rst} A(rs,tp) B(rs,tq) = \sum_{rst} X(p,rst) B(rst,q)
    else if (target_a == 4 && target_b == 4) {
        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = a->d1_ * a->d2_ * a->d3_;
        lda = k;
        ldb = n;
        ldc = n;

        SharedTensor2d temp = SharedTensor2d(new Tensor2d("temp", a->d4_, a->d1_, a->d2_, a->d3_));
        temp->sort(4123, a, 1.0, 0.0);

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, &(temp->A2d_[0][0]), lda, &(b->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
        }
        temp.reset();
    }

    // C(p,q) = \sum_{rst} A(rs,tp) B(rs,qt)
    else if (target_a == 4 && target_b == 3) {
        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = a->d2_ * a->d3_ * a->d1_;
        lda = k;
        ldb = n;
        ldc = n;

        SharedTensor2d temp1 = SharedTensor2d(new Tensor2d("temp", a->d4_, a->d1_, a->d2_, a->d3_));
        SharedTensor2d temp2 = SharedTensor2d(new Tensor2d("temp", b->d1_, b->d2_, b->d4_, b->d3_));
        temp1->sort(4123, a, 1.0, 0.0);
        temp2->sort(1243, b, 1.0, 0.0);
        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, &(temp1->A2d_[0][0]), lda, &(temp2->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
        }
        temp1.reset();
        temp2.reset();
    }

    // C(p,q) = \sum_{rst} A(rp,st) B(rq,st)
    else if (target_a == 2 && target_b == 2) {
        ta = 'n';
        tb = 't';
        m = dim1_;
        n = dim2_;
        k = a->d1_ * a->d3_ * a->d4_;
        lda = k;
        ldb = k;
        ldc = n;

        SharedTensor2d temp1 = SharedTensor2d(new Tensor2d("temp1", a->d2_, a->d1_, a->d3_, a->d4_));
        SharedTensor2d temp2 = SharedTensor2d(new Tensor2d("temp2", b->d2_, b->d1_, b->d3_, b->d4_));
        temp1->sort(2134, a, 1.0, 0.0);
        temp2->sort(2134, b, 1.0, 0.0);
        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, &(temp1->A2d_[0][0]), lda, &(temp2->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
        }
        temp1.reset();
        temp2.reset();
    }

    // C(p,q) = \sum_{rst} A(rs,pt) B(rs,qt)
    else if (target_a == 3 && target_b == 3) {
        ta = 't';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = a->d1_ * a->d2_ * a->d4_;
        lda = m;
        ldb = n;
        ldc = n;

        SharedTensor2d temp1 = SharedTensor2d(new Tensor2d("temp1", a->d1_, a->d2_, a->d4_, a->d3_));
        SharedTensor2d temp2 = SharedTensor2d(new Tensor2d("temp2", b->d1_, b->d2_, b->d4_, b->d3_));
        temp1->sort(1243, a, 1.0, 0.0);
        temp2->sort(1243, b, 1.0, 0.0);
        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, &(temp1->A2d_[0][0]), lda, &(temp2->A2d_[0][0]), ldb, beta, &(A2d_[0][0]), ldc);
        }
        temp1.reset();
        temp2.reset();
    }

    else {
         outfile->Printf("contract442: Unrecognized targets!");

    }

}//

void Tensor2d::gemv(bool transa, const SharedTensor2d& a, const SharedTensor1d& b, double alpha, double beta)
{
    char ta = transa ? 't' : 'n';
    int m, n, k, incx, incy, lda;

    m = a->dim1();
    n = a->dim2();
    lda = n;
    incx = 1;// increments in elements of b vector
    incy = 1;// increments in elements of A2d_

    if (m && n) {
       C_DGEMV(ta, m, n, alpha, &(a->A2d_[0][0]), lda, b->A1d_, incx, beta, &(A2d_[0][0]), incy);
    }
}//

void Tensor2d::gemv(bool transa, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta)
{
    char ta = transa ? 't' : 'n';
    int m, n, k, incx, incy, lda;

    m = a->dim1();
    n = a->dim2();
    lda = n;
    incx = 1;// increments in elements of b vector
    incy = 1;// increments in elements of A2d_

    if (m && n) {
       C_DGEMV(ta, m, n, alpha, &(a->A2d_[0][0]), lda, &(b->A2d_[0][0]), incx, beta, &(A2d_[0][0]), incy);
    }
}//

void Tensor2d::davidson(int n_eigval, const SharedTensor2d& eigvectors, const SharedTensor1d& eigvalues, double cutoff, int print)
{

    david(A2d_, dim1_, n_eigval, eigvalues->A1d_, eigvectors->A2d_, cutoff, print);

}//

void Tensor2d::add(const SharedTensor2d &a)
{
    ULI length = (ULI)dim1_ * (ULI)dim2_;
    C_DAXPY(length, 1.0, a->A2d_[0], 1, A2d_[0], 1);
}//

void Tensor2d::add(double **a)
{
    ULI length = (ULI)dim1_ * (ULI)dim2_;
    C_DAXPY(length, 1.0, a[0], 1, A2d_[0], 1);
}//

void Tensor2d::add(double alpha, const SharedTensor2d &Adum)
{
    SharedTensor2d temp = SharedTensor2d(new Tensor2d(Adum->dim1_, Adum->dim2_));
    temp->copy(Adum);
    temp->scale(alpha);
    add(temp);
}//

void Tensor2d::add(int i, int j, double value)
{
  A2d_[i][j]+=value;
}//

void Tensor2d::subtract(const SharedTensor2d &a)
{
    ULI length = (ULI)dim1_ * (ULI)dim2_;
    C_DAXPY(length, -1.0, a->A2d_[0], 1, A2d_[0], 1);
}//

void Tensor2d::subtract(int i, int j, double value)
{
  A2d_[i][j]-=value;
}//

void Tensor2d::axpy(double **a, double alpha)
{
    ULI length = (ULI)dim1_ * (ULI)dim2_;
    C_DAXPY(length, alpha, a[0], 1, A2d_[0], 1);
}//

void Tensor2d::axpy(const SharedTensor2d &a, double alpha)
{
    ULI length = (ULI)dim1_ * (ULI)dim2_;
    C_DAXPY(length, alpha, a->A2d_[0], 1, A2d_[0], 1);
}//

void Tensor2d::axpy(ULI length, int inc_a, const SharedTensor2d &a, int inc_2d, double alpha)
{
    C_DAXPY(length, alpha, a->A2d_[0], inc_a, A2d_[0], inc_2d);
}//

void Tensor2d::axpy(ULI length, int start_a, int inc_a, const SharedTensor2d &A, int start_2d, int inc_2d, double alpha)
{
    C_DAXPY(length, alpha, A->A2d_[0]+start_a, inc_a, A2d_[0]+start_2d, inc_2d);
}//

double Tensor2d::norm()
{
    double value = 0.0;
    ULI length = (ULI)dim1_ * (ULI)dim2_;
    value = C_DNRM2(length, A2d_[0], 1);
    return value;
}//

double **Tensor2d::transpose2()
{
    double** temp = block_matrix(dim2_, dim1_);
    memset(temp[0], 0, sizeof(double)*dim1_*dim2_);
      #pragma omp parallel for
      for (int i=0; i<dim2_; ++i) {
	for (int j=0; j<dim1_; ++j) {
	  temp[i][j] = A2d_[j][i];
	}
      }

    return temp;
}//

SharedTensor2d Tensor2d::transpose()
{
    SharedTensor2d temp = SharedTensor2d(new Tensor2d(dim2_, dim1_));
      #pragma omp parallel for
      for (int i=0; i<dim2_; ++i) {
	for (int j=0; j<dim1_; ++j) {
	  temp->A2d_[i][j] = A2d_[j][i];
	}
      }

    return temp;
}//

void Tensor2d::trans(const SharedTensor2d &A)
{
      #pragma omp parallel for
      for (int i=0; i<dim1_; ++i) {
	for (int j=0; j<dim2_; ++j) {
	     A2d_[i][j] = A->A2d_[j][i];
	}
      }

}//

void Tensor2d::trans(double **A)
{
      #pragma omp parallel for
      for (int i=0; i<dim1_; ++i) {
	for (int j=0; j<dim2_; ++j) {
	     A2d_[i][j] = A[j][i];
	}
      }

}//

void Tensor2d::copy(double **a)
{
    //size_t size = dim1_ * dim2_ * sizeof(double);
    //if (size) memcpy(&(A2d_[0][0]), &(a[0][0]), size);
    ULI length;
    length = (ULI)dim1_ * (ULI)dim2_;
    C_DCOPY(length, a[0], 1, A2d_[0], 1);
}

void Tensor2d::copy(const SharedTensor2d &Adum)
{
   // Make sure that matrices are in the same size
   bool same = true;
   if (dim2_ != Adum->dim2_ || dim1_ != Adum->dim1_)  same = false;

    if (same == false) {
        release();
        dim1_ = Adum->dim1_;
        dim2_ = Adum->dim2_;
        memalloc();
    }

    // If matrices are in the same size
    ULI length;
    length = (ULI)dim1_ * (ULI)dim2_;
    if (dim1_ != 0 && dim2_ != 0) {
        //memcpy(A2d_[0], Adum->A2d_[0], dim1_ * dim2_ * sizeof(double));
        C_DCOPY(length, Adum->A2d_[0], 1, A2d_[0], 1);
    }
}//

void Tensor2d::copy(ULI length, const SharedTensor2d &A, int inc_a, int inc_2d)
{
    C_DCOPY(length, A->A2d_[0], inc_a, A2d_[0], inc_2d);
}//

void Tensor2d::copy(const SharedTensor2d &A, int start)
{
    memcpy(A2d_[0], A->A2d_[0]+start, dim1_ * dim2_ * sizeof(double));
}//

void Tensor2d::pcopy(const SharedTensor2d &A, int dim_copy, int dim_skip)
{
    double *temp = new double[dim_copy];
    int syc = 0;
    // A[m] is getting the pointer to the m-th row of A.
    // A[0]+m is getting the pointer to the m-th element of A.
    for (int i = 0; i < dim1_*dim2_; i+=dim_copy) {
         memcpy(temp, A->A2d_[0]+syc, dim_copy * sizeof(double));
         memcpy(A2d_[0]+i, temp, dim_copy * sizeof(double));
         syc += dim_copy + dim_skip;
    }
    delete [] temp;

}//

void Tensor2d::pcopy(const SharedTensor2d &A, int dim_copy, int dim_skip, int start)
{
    double *temp = new double[dim_copy];
    int syc = 0;
    // A[m] is getting the pointer to the m-th row of A.
    // A[0]+m is getting the pointer to the m-th element of A.
    for (int i = 0; i < dim1_*dim2_; i+=dim_copy) {
         memcpy(temp, A->A2d_[0]+start+syc, dim_copy * sizeof(double));
         memcpy(A2d_[0]+i, temp, dim_copy * sizeof(double));
         syc += dim_copy + dim_skip;
    }
    delete [] temp;

}//

void Tensor2d::diagonalize(const SharedTensor2d &eigvectors, const SharedTensor1d &eigvalues, double cutoff)
{
   sq_rsp(dim1_, dim2_, A2d_, eigvalues->A1d_, 1, eigvectors->A2d_, cutoff);

}//

void Tensor2d::cdsyev(char jobz, char uplo, const SharedTensor2d& eigvectors, const SharedTensor1d& eigvalues)
{
      if (dim1_) {
	int lwork=3*dim2_;
	double **work = block_matrix(dim1_,lwork);
	memset(work[0],0.0,sizeof(double)*dim1_*lwork);
        C_DSYEV(jobz, uplo, dim1_, &(A2d_[0][0]), dim2_, eigvalues->A1d_, &(work[0][0]), lwork);
	free_block(work);
      }
}//

void Tensor2d::cdgesv(const SharedTensor1d& Xvec)
{
      if (dim1_) {
	int errcod;
	int *ipiv = new int[dim1_];
	memset(ipiv,0,sizeof(int)*dim1_);
	errcod=0;
	errcod = C_DGESV(dim1_, 1, &(A2d_[0][0]), dim2_, &(ipiv[0]), Xvec->A1d_, dim2_);
        delete [] ipiv;
      }
}//

void Tensor2d::cdgesv(const SharedTensor1d& Xvec, int errcod)
{
      if (dim1_) {
	int *ipiv = new int[dim1_];
	memset(ipiv,0,sizeof(int)*dim1_);
	errcod=0;
	errcod = C_DGESV(dim1_, 1, &(A2d_[0][0]), dim2_, &(ipiv[0]), Xvec->A1d_, dim2_);
        delete [] ipiv;
      }
}//

void Tensor2d::cdgesv(double* Xvec)
{
      if (dim1_) {
	int errcod;
	int *ipiv = new int[dim1_];
	memset(ipiv,0,sizeof(int)*dim1_);
	errcod=0;
	errcod = C_DGESV(dim1_, 1, &(A2d_[0][0]), dim2_, &(ipiv[0]), Xvec, dim2_);
        delete [] ipiv;
      }
}//

void Tensor2d::cdgesv(double* Xvec, int errcod)
{
      if (dim1_) {
	int *ipiv = new int[dim1_];
	memset(ipiv,0,sizeof(int)*dim1_);
	errcod=0;
	errcod = C_DGESV(dim1_, 1, &(A2d_[0][0]), dim2_, &(ipiv[0]), Xvec, dim2_);
        delete [] ipiv;
      }
}//

void Tensor2d::lineq_flin(const SharedTensor1d& Xvec, double *det)
{
      if (dim1_) {
	flin(A2d_, Xvec->A1d_, dim1_, 1, det);
      }
}//

void Tensor2d::lineq_flin(double* Xvec, double *det)
{
      if (dim1_) {
	flin(A2d_, Xvec, dim1_, 1, det);
      }
}//

void Tensor2d::lineq_pople(const SharedTensor1d& Xvec, int num_vecs, double cutoff)
{
      if (dim1_) {
	pople(A2d_, Xvec->A1d_, dim1_, num_vecs, cutoff, "outfile", 0);
      }
}//

void Tensor2d::lineq_pople(double* Xvec, int num_vecs, double cutoff)
{
      if (dim1_) {
	pople(A2d_, Xvec, dim1_, num_vecs, cutoff, "outfile", 0);
      }
}//

void Tensor2d::level_shift(double value)
{
      #pragma omp parallel for
      for (int i=0; i<dim1_; ++i) {
	subtract(i, i, value);
      }

}//

void Tensor2d::outer_product(const SharedTensor1d &x, const SharedTensor1d &y)
{
  #pragma omp parallel for
  for (int i=0; i < x->dim1_; i++) {
     for (int j=0; j < y->dim1_; j++) {
          A2d_[i][j] = x->A1d_[i] * y->A1d_[j];
     }
  }

}//

// TODO:
// DGER compute the rank-one update of a general matrix: A <-- A + alpha * x * yT
// dger(m, n, alpha, x, incx, y, incy, a, lda);

void Tensor2d::scale(double a)
{
    //size_t size;
    ULI size;
    size = (ULI)dim1_ * (ULI)dim2_;
    if (size) C_DSCAL(size, a, &(A2d_[0][0]), 1);
}//

void Tensor2d::scale_row(int m, double a)
{
    C_DSCAL((ULI)dim1_, a, &(A2d_[m][0]), 1);
}//

void Tensor2d::scale_column(int n, double a)
{
    C_DSCAL((ULI)dim2_, a, &(A2d_[0][n]), dim1_);
}//

void Tensor2d::identity()
{
    zero();
    for (int i=0; i<dim1_; ++i) A2d_[i][i] = 1.0;
}//

double Tensor2d::trace()
{
    double value = 0.0;
    for (int i=0; i<dim1_; ++i) value += A2d_[i][i];
    return value;
}//

void Tensor2d::transform(const SharedTensor2d& a, const SharedTensor2d& transformer)
{
    SharedTensor2d temp = SharedTensor2d(new Tensor2d(a->dim1_, transformer->dim2_));
    temp->gemm(false, false, a, transformer, 1.0, 0.0);
    gemm(true, false, transformer, temp, 1.0, 0.0);
}//

void Tensor2d::back_transform(const SharedTensor2d& a, const SharedTensor2d& transformer)
{
    SharedTensor2d temp = SharedTensor2d(new Tensor2d(a->dim1_, transformer->dim2_));
    temp->gemm(false, true, a, transformer, 1.0, 0.0);
    gemm(false, false, transformer, temp, 1.0, 0.0);
}//

void Tensor2d::back_transform(const SharedTensor2d& a, const SharedTensor2d& transformer, double alpha, double beta)
{
    SharedTensor2d temp = SharedTensor2d(new Tensor2d(a->dim1_, transformer->dim2_));
    temp->gemm(false, true, a, transformer, 1.0, 0.0);
    gemm(false, false, transformer, temp, alpha, beta);
}//

void Tensor2d::pseudo_transform(const SharedTensor2d& a, const SharedTensor2d& transformer)
{
    SharedTensor2d temp = SharedTensor2d(new Tensor2d(a->dim1_, transformer->dim2_));
    temp->gemm(false, false, a, transformer, 1.0, 0.0);
    gemm(false, false, transformer, temp, 1.0, 0.0);
}//

void Tensor2d::triple_gemm(const SharedTensor2d& a, const SharedTensor2d& b, const SharedTensor2d& c)
{
  if (a->dim2_ == b->dim1_ && b->dim2_ == c->dim1_ && a->dim1_ == dim1_ && c->dim2_ == dim2_) {
    SharedTensor2d bc = SharedTensor2d(new Tensor2d(b->dim1_, c->dim2_));
    bc->gemm(false, false, b, c, 1.0, 0.0);
    gemm(false, false, a, bc, 1.0, 0.0);
  }
  else {
    outfile->Printf("\n Warning!!! Matrix dimensions do NOT match in triple_gemm().\n");

  }

}//

double Tensor2d::vector_dot(double **rhs)
{
    double value = 0.0;
    //size_t size = dim1_ * dim2_;
    ULI size;
    size = (ULI)dim1_ * (ULI)dim2_;
    if (size) value += C_DDOT(size, (&A2d_[0][0]), 1, &(rhs[0][0]), 1);
    return value;
}//

double Tensor2d::vector_dot(const SharedTensor2d &rhs)
{
    double value = 0.0;
    //size_t size = dim1_ * dim2_;
    ULI size;
    size = (ULI)dim1_ * (ULI)dim2_;
    if (size) value += C_DDOT(size, (&A2d_[0][0]), 1, &(rhs->A2d_[0][0]), 1);
    return value;
}//

void Tensor2d::write(std::shared_ptr<psi::PSIO> psio, unsigned int fileno)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) already_open = true;
    else psio->open(fileno, PSIO_OPEN_OLD);
    psio->write_entry(fileno, const_cast<char*>(name_.c_str()), (char*)A2d_[0], sizeof(double) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);     // Close and keep
}//

void Tensor2d::write(std::shared_ptr<psi::PSIO> psio, unsigned int fileno, psio_address start, psio_address *end)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) already_open = true;
    else psio->open(fileno, PSIO_OPEN_OLD);
    ULI size_ = (ULI)dim1_*dim2_*sizeof(double);
    psio->write(fileno, const_cast<char*>(name_.c_str()), (char*)A2d_[0], size_, start, end);
    if (!already_open) psio->close(fileno, 1);     // Close and keep
}//

void Tensor2d::write(psi::PSIO* const psio, unsigned int fileno)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) already_open = true;
    else psio->open(fileno, PSIO_OPEN_OLD);
    psio->write_entry(fileno, const_cast<char*>(name_.c_str()), (char*)A2d_[0], sizeof(double) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);     // Close and keep
}//

void Tensor2d::write(psi::PSIO* const psio, unsigned int fileno, psio_address start, psio_address *end)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) already_open = true;
    else psio->open(fileno, PSIO_OPEN_OLD);
    ULI size_ = (ULI)dim1_*dim2_*sizeof(double);
    psio->write(fileno, const_cast<char*>(name_.c_str()), (char*)A2d_[0], size_, start, end);
    if (!already_open) psio->close(fileno, 1);     // Close and keep
}//

void Tensor2d::write(psi::PSIO& psio, unsigned int fileno)
{
    write(&psio, fileno);
}//

void Tensor2d::write(psi::PSIO& psio, unsigned int fileno, psio_address start, psio_address *end)
{
    write(&psio, fileno, start, end);
}//

void Tensor2d::write(std::shared_ptr<psi::PSIO> psio, const string& filename, unsigned int fileno)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) already_open = true;
    else psio->open(fileno, PSIO_OPEN_OLD);
    psio->write_entry(fileno, const_cast<char*>(filename.c_str()), (char*)A2d_[0], sizeof(double) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);     // Close and keep
}//

void Tensor2d::write(std::shared_ptr<psi::PSIO> psio, unsigned int fileno, bool three_index, bool symm)
{
    // Form Lower triangular part
    if (three_index && symm) {
        int ntri_col = 0.5 * d2_ * (d2_ +1);
        SharedTensor2d temp = SharedTensor2d(new Tensor2d("temp", d1_, ntri_col));
        #pragma omp parallel for
        for (int R = 0; R < d1_; R++) {
              for (int p = 0; p < d2_; p++) {
                   for (int q = 0; q < d3_; q++) {
                        int pq = col_idx_[p][q];
                        int pq_sym = index2(p,q);
                        temp->set(R, pq_sym, A2d_[R][pq]);
                   }
              }
        }

        // Check to see if the file is open
        bool already_open = false;
        if (psio->open_check(fileno)) already_open = true;
        else psio->open(fileno, PSIO_OPEN_OLD);
        psio->write_entry(fileno, const_cast<char*>(name_.c_str()), (char*)temp->A2d_[0], sizeof(double) * dim1_ * ntri_col);
        if (!already_open) psio->close(fileno, 1);     // Close and keep
        temp.reset();
    }

    else {
        // Check to see if the file is open
        bool already_open = false;
        if (psio->open_check(fileno)) already_open = true;
        else psio->open(fileno, PSIO_OPEN_OLD);
        psio->write_entry(fileno, const_cast<char*>(name_.c_str()), (char*)A2d_[0], sizeof(double) * dim1_ * dim2_);
        if (!already_open) psio->close(fileno, 1);     // Close and keep
    }

}//

void Tensor2d::write(std::shared_ptr<psi::PSIO> psio, const string& filename, unsigned int fileno, bool three_index, bool symm)
{
    // Form Lower triangular part
    if (three_index && symm) {
        int ntri_col = 0.5 * d2_ * (d2_ +1);
        SharedTensor2d temp = SharedTensor2d(new Tensor2d("temp", d1_, ntri_col));
        #pragma omp parallel for
        for (int R = 0; R < d1_; R++) {
              for (int p = 0; p < d2_; p++) {
                   for (int q = 0; q < d3_; q++) {
                        int pq = col_idx_[p][q];
                        int pq_sym = index2(p,q);
                        temp->set(R, pq_sym, A2d_[R][pq]);
                   }
              }
        }

        // Check to see if the file is open
        bool already_open = false;
        if (psio->open_check(fileno)) already_open = true;
        else psio->open(fileno, PSIO_OPEN_OLD);
        psio->write_entry(fileno, const_cast<char*>(filename.c_str()), (char*)temp->A2d_[0], sizeof(double) * dim1_ * ntri_col);
        if (!already_open) psio->close(fileno, 1);     // Close and keep
        temp.reset();
    }

    else {
        // Check to see if the file is open
        bool already_open = false;
        if (psio->open_check(fileno)) already_open = true;
        else psio->open(fileno, PSIO_OPEN_OLD);
        psio->write_entry(fileno, const_cast<char*>(filename.c_str()), (char*)A2d_[0], sizeof(double) * dim1_ * dim2_);
        if (!already_open) psio->close(fileno, 1);     // Close and keep
    }

}//

void Tensor2d::write_symm(std::shared_ptr<psi::PSIO> psio, unsigned int fileno)
{
        // Form Lower triangular part
        int ntri_col = 0.5 * dim1_ * (dim1_ +1);
        SharedTensor1d temp = SharedTensor1d(new Tensor1d("temp", ntri_col));
        #pragma omp parallel for
        for (int p = 0; p < dim1_; p++) {
              for (int q = 0; q <= p; q++) {
                   int pq = index2(p,q);
                   temp->set(pq, A2d_[p][q]);
              }
        }

        // Check to see if the file is open
        bool already_open = false;
        if (psio->open_check(fileno)) already_open = true;
        else psio->open(fileno, PSIO_OPEN_OLD);
        psio->write_entry(fileno, const_cast<char*>(name_.c_str()), (char*)&(temp->A1d_[0]), sizeof(double) * ntri_col);
        if (!already_open) psio->close(fileno, 1);     // Close and keep
        temp.reset();

}//

void Tensor2d::write_anti_symm(std::shared_ptr<psi::PSIO> psio, unsigned int fileno)
{
        // Form Lower triangular part
	int ntri_row, ntri_col;
	if (dim1_ > 1) {
            ntri_row = 0.5 * d1_ * (d1_ - 1);
	}
	else if (dim1_ == 1) {
            ntri_row = 1;
	}
	if (dim2_ > 1) {
            ntri_col = 0.5 * d3_ * (d3_ - 1);
	}
	else if (dim2_ == 1) {
            ntri_col = 1;
	}
        SharedTensor2d temp = SharedTensor2d(new Tensor2d("temp", ntri_row, ntri_col));
        #pragma omp parallel for
        for (int p = 1; p < d1_; p++) {
              for (int q = 0; q < p; q++) {
                   int pq = row_idx_[p][q];
                   int pq2 = idx_asym(p,q);
                   for (int r = 1; r < d3_; r++) {
                        for (int s = 0; s < r; s++) {
                             int rs = col_idx_[r][s];
                             int rs2 = idx_asym(r,s);
                             temp->set(pq2, rs2, A2d_[pq][rs]);
                        }
                   }
              }
        }

        // Check to see if the file is open
        bool already_open = false;
        if (psio->open_check(fileno)) already_open = true;
        else psio->open(fileno, PSIO_OPEN_OLD);
        psio->write_entry(fileno, const_cast<char*>(name_.c_str()), (char*)temp->A2d_[0], sizeof(double) * ntri_row * ntri_col);
        if (!already_open) psio->close(fileno, 1);     // Close and keep
        temp.reset();

}//

void Tensor2d::read(psi::PSIO* psio, unsigned int fileno)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) already_open = true;
    else psio->open(fileno, PSIO_OPEN_OLD);
    psio->read_entry(fileno, const_cast<char*>(name_.c_str()), (char*)A2d_[0], sizeof(double) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);     // Close and keep
}

void Tensor2d::read(psi::PSIO* psio, unsigned int fileno, psio_address start, psio_address *end)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) already_open = true;
    else psio->open(fileno, PSIO_OPEN_OLD);
    ULI size_ = (ULI)dim1_*dim2_*sizeof(double);
    psio->read(fileno, const_cast<char*>(name_.c_str()), (char*)A2d_[0], size_, start, end);
    if (!already_open) psio->close(fileno, 1);     // Close and keep
}

void Tensor2d::read(std::shared_ptr<psi::PSIO> psio, unsigned int fileno)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) already_open = true;
    else psio->open(fileno, PSIO_OPEN_OLD);
    psio->read_entry(fileno, const_cast<char*>(name_.c_str()), (char*)A2d_[0], sizeof(double) * dim1_ * dim2_);
    if (!already_open) psio->close(fileno, 1);     // Close and keep
}

void Tensor2d::read(std::shared_ptr<psi::PSIO> psio, unsigned int fileno, psio_address start, psio_address *end)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) already_open = true;
    else psio->open(fileno, PSIO_OPEN_OLD);
    ULI size_ = (ULI)dim1_*dim2_*sizeof(double);
    psio->read(fileno, const_cast<char*>(name_.c_str()), (char*)A2d_[0], size_, start, end);
    if (!already_open) psio->close(fileno, 1);     // Close and keep
}

void Tensor2d::read(psi::PSIO& psio, unsigned int fileno)
{
    read(&psio, fileno);
}//

void Tensor2d::read(psi::PSIO& psio, unsigned int fileno, psio_address start, psio_address *end)
{
    read(&psio, fileno, start, end);
}//

void Tensor2d::read(std::shared_ptr<psi::PSIO> psio, unsigned int fileno, bool three_index, bool symm)
{
    // Form Lower triangular part
    if (three_index && symm) {
        int ntri_col = 0.5 * d2_ * (d2_ +1);
        SharedTensor2d temp = SharedTensor2d(new Tensor2d("temp", d1_, ntri_col));

        // Check to see if the file is open
        bool already_open = false;
        if (psio->open_check(fileno)) already_open = true;
        else psio->open(fileno, PSIO_OPEN_OLD);
        psio->read_entry(fileno, const_cast<char*>(name_.c_str()), (char*)temp->A2d_[0], sizeof(double) * dim1_ * ntri_col);
        if (!already_open) psio->close(fileno, 1);     // Close and keep

        #pragma omp parallel for
        for (int R = 0; R < d1_; R++) {
              for (int p = 0; p < d2_; p++) {
                   for (int q = 0; q < d3_; q++) {
                        int pq = col_idx_[p][q];
                        int pq_sym = index2(p,q);
                        A2d_[R][pq] = temp->get(R, pq_sym);
                   }
              }
        }
        temp.reset();
    }

    else {
        // Check to see if the file is open
        bool already_open = false;
        if (psio->open_check(fileno)) already_open = true;
        else psio->open(fileno, PSIO_OPEN_OLD);
        psio->read_entry(fileno, const_cast<char*>(name_.c_str()), (char*)A2d_[0], sizeof(double) * dim1_ * dim2_);
        if (!already_open) psio->close(fileno, 1);     // Close and keep
    }

}//

void Tensor2d::read_symm(std::shared_ptr<psi::PSIO> psio, unsigned int fileno)
{
        // Form Lower triangular part
        int ntri_col = 0.5 * dim1_ * (dim1_ +1);
        SharedTensor1d temp = SharedTensor1d(new Tensor1d("temp", ntri_col));

        // Check to see if the file is open
        bool already_open = false;
        if (psio->open_check(fileno)) already_open = true;
        else psio->open(fileno, PSIO_OPEN_OLD);
        psio->read_entry(fileno, const_cast<char*>(name_.c_str()), (char*)&(temp->A1d_[0]), sizeof(double) * ntri_col);
        if (!already_open) psio->close(fileno, 1);     // Close and keep

        #pragma omp parallel for
        for (int p = 0; p < dim1_; p++) {
              for (int q = 0; q <= p; q++) {
                   int pq = index2(p,q);
                   A2d_[p][q] = temp->get(pq);
                   A2d_[q][p] = temp->get(pq);
              }
        }
        temp.reset();
}//

void Tensor2d::read_anti_symm(std::shared_ptr<psi::PSIO> psio, unsigned int fileno)
{
        // Form Lower triangular part
	int ntri_row, ntri_col;
	if (dim1_ > 1) {
            ntri_row = 0.5 * d1_ * (d1_ - 1);
	}
	else if (dim1_ == 1) {
            ntri_row = 1;
	}
	if (dim2_ > 1) {
            ntri_col = 0.5 * d3_ * (d3_ - 1);
	}
	else if (dim2_ == 1) {
            ntri_col = 1;
	}


        SharedTensor2d temp = SharedTensor2d(new Tensor2d("temp", ntri_row, ntri_col));

        // Check to see if the file is open
        bool already_open = false;
        if (psio->open_check(fileno)) already_open = true;
        else psio->open(fileno, PSIO_OPEN_OLD);
        psio->read_entry(fileno, const_cast<char*>(name_.c_str()), (char*)temp->A2d_[0], sizeof(double) * ntri_row * ntri_col);
        if (!already_open) psio->close(fileno, 1);     // Close and keep

        #pragma omp parallel for
        for (int p = 1; p < d1_; p++) {
              for (int q = 0; q < p; q++) {
                   int pq = row_idx_[p][q];
                   int qp = row_idx_[q][p];
                   int pq2 = idx_asym(p,q);
                   for (int r = 1; r < d3_; r++) {
                        for (int s = 0; s < r; s++) {
                             int rs = col_idx_[r][s];
                             int sr = col_idx_[s][r];
                             int rs2 = idx_asym(r,s);
                             double value = temp->get(pq2, rs2);
                             A2d_[pq][rs] = value;
                             A2d_[pq][sr] = -1.0 * value;
                             A2d_[qp][rs] = -1.0 * value;
                             A2d_[qp][sr] = value;
                        }
                   }
              }
        }
        temp.reset();

}//

bool Tensor2d::read(PSIO* psio, int itap, const char *label, int dim)
{
    int ntri = 0.5 * dim * (dim + 1);
    double *mybuffer = init_array(ntri);
    memset(mybuffer, 0, sizeof(double)*ntri);
    IWL::read_one(psio, itap, label, mybuffer, ntri, 0, 0, "outfile");

    double **Asq = block_matrix(dim, dim);
    memset(Asq[0], 0, sizeof(double)*dim*dim);
    tri_to_sq(mybuffer,Asq,dim);
    free(mybuffer);

    set(Asq);
    free_block(Asq);
    return true;
}//

bool Tensor2d::read(std::shared_ptr<psi::PSIO> psio, int itap, const char *label, int dim)
{
    int ntri = 0.5 * dim * (dim + 1);
    double *mybuffer = init_array(ntri);
    memset(mybuffer, 0, sizeof(double)*ntri);
    IWL::read_one(psio.get(), itap, label, mybuffer, ntri, 0, 0, "outfile");

    double **Asq = block_matrix(dim, dim);
    memset(Asq[0], 0, sizeof(double)*dim*dim);
    tri_to_sq(mybuffer,Asq,dim);
    free(mybuffer);

    set(Asq);
    free_block(Asq);
    return true;
}//

void Tensor2d::save(std::shared_ptr<psi::PSIO> psio, unsigned int fileno)
{
    write(psio, fileno);
    release();
}//

void Tensor2d::save(psi::PSIO* const psio, unsigned int fileno)
{
    write(psio, fileno);
    release();
}//

void Tensor2d::save(psi::PSIO& psio, unsigned int fileno)
{
    write(&psio, fileno);
    release();
}//

void Tensor2d::load(std::shared_ptr<psi::PSIO> psio, unsigned int fileno, string name, int d1,int d2)
{
    init(name,d1,d2);
    read(psio, fileno);
}//

void Tensor2d::load(psi::PSIO* const psio, unsigned int fileno, string name, int d1,int d2)
{
    init(name,d1,d2);
    read(psio, fileno);
}//

void Tensor2d::load(psi::PSIO& psio, unsigned int fileno, string name, int d1,int d2)
{
    init(name,d1,d2);
    read(&psio, fileno);
}//

void Tensor2d::mywrite(const string& filename)
{
      // write binary data
      ofstream OutFile;
      OutFile.open(const_cast<char*>(filename.c_str()), ios::out | ios::binary);
      OutFile.write((char*)A2d_[0], dim1_*dim2_*sizeof(double));
      OutFile.close();
}//

void Tensor2d::mywrite(int fileno)
{
      ostringstream convert;
      convert << fileno;
      std::string scr = PSIOManager::shared_object()->get_default_path();
      std::string pid_ = psio_getpid();
      //std::string fname = scr + "psi_dfocc." + convert.str();
      std::string fname = scr + "psi." + pid_  + "." + convert.str();

      // write binary data
      ofstream OutFile;
      OutFile.open(const_cast<char*>(fname.c_str()), ios::out | ios::binary);
      OutFile.write((char*)A2d_[0], dim1_*dim2_*sizeof(double));
      OutFile.close();
}//

void Tensor2d::mywrite(int fileno, bool append)
{
      ostringstream convert;
      convert << fileno;
      std::string scr = PSIOManager::shared_object()->get_default_path();
      std::string pid_ = psio_getpid();
      //std::string fname = scr + "psi_dfocc." + convert.str();
      std::string fname = scr + "psi." + pid_  + "." + convert.str();

      // write binary data
      ofstream OutFile;
      if (append) OutFile.open(const_cast<char*>(fname.c_str()), ios::out | ios::binary | ios::app);
      else OutFile.open(const_cast<char*>(fname.c_str()), ios::out | ios::binary);
      OutFile.write((char*)A2d_[0], dim1_*dim2_*sizeof(double));
      OutFile.close();
}//

void Tensor2d::myread(const string& filename)
{
      // read binary data
      ifstream InFile;
      InFile.open(const_cast<char*>(filename.c_str()), ios::in | ios::binary);
      InFile.read((char*)A2d_[0], dim1_*dim2_*sizeof(double));
      InFile.close();

}//

void Tensor2d::myread(int fileno)
{
      ostringstream convert;
      convert << fileno;
      std::string scr = PSIOManager::shared_object()->get_default_path();
      std::string pid_ = psio_getpid();
      //std::string fname = scr + "psi_dfocc." + convert.str();
      std::string fname = scr + "psi." + pid_  + "." + convert.str();

      // read binary data
      ifstream InFile;
      InFile.open(const_cast<char*>(fname.c_str()), ios::in | ios::binary);
      InFile.read((char*)A2d_[0], dim1_*dim2_*sizeof(double));
      InFile.close();

}//

void Tensor2d::myread(int fileno, bool append)
{
      ostringstream convert;
      convert << fileno;
      std::string scr = PSIOManager::shared_object()->get_default_path();
      std::string pid_ = psio_getpid();
      //std::string fname = scr + "psi_dfocc." + convert.str();
      std::string fname = scr + "psi." + pid_  + "." + convert.str();

      // read binary data
      ifstream InFile;
      if (append) InFile.open(const_cast<char*>(fname.c_str()), ios::in | ios::binary | ios::app);
      else InFile.open(const_cast<char*>(fname.c_str()), ios::in | ios::binary);
      InFile.read((char*)A2d_[0], dim1_*dim2_*sizeof(double));
      InFile.close();

}//

void Tensor2d::myread(int fileno, ULI start)
{
      ostringstream convert;
      convert << fileno;
      std::string scr = PSIOManager::shared_object()->get_default_path();
      std::string pid_ = psio_getpid();
      std::string fname = scr + "psi." + pid_  + "." + convert.str();
      //std::string fname = scr + "psi_dfocc." + convert.str();

      // read binary data
      ifstream InFile;
      InFile.open(const_cast<char*>(fname.c_str()), ios::in | ios::binary);
      InFile.seekg(start, ios::beg);
      InFile.read((char*)A2d_[0], dim1_*dim2_*sizeof(double));
      InFile.close();

}//

double **Tensor2d::to_block_matrix()
{
    double **temp = block_matrix(dim1_, dim2_);
    //memcpy(&(temp[0][0]), &(A2d_[0][0]), dim1_ * dim2_ * sizeof(double));
    ULI length;
    length = (ULI)dim1_ * (ULI)dim2_;
    C_DCOPY(length, A2d_[0], 1, temp[0], 1);
    return temp;
}//

double *Tensor2d::to_lower_triangle()
{
    if (dim1_ != dim2_) return NULL;
    int ntri = 0.5 * dim1_ * (dim1_ + 1);
    double *tri = new double[ntri];
    double **temp = to_block_matrix();
    sq_to_tri(temp, tri, dim1_);
    free_block(temp);
    return tri;
}//

void Tensor2d::to_shared_matrix(SharedMatrix A)
{
      #pragma omp parallel for
      for (int i=0; i<dim1_; ++i) {
        for (int j=0; j<dim2_; ++j) {
          A->set(0,i,j,A2d_[i][j]);
        }
      }
}//

void Tensor2d::to_matrix(SharedMatrix A)
{
      #pragma omp parallel for
      for (int i=0; i<dim1_; ++i) {
        for (int j=0; j<dim2_; ++j) {
          A->set(i,j,A2d_[i][j]);
        }
      }
}//

void Tensor2d::to_pointer(double *A)
{
      #pragma omp parallel for
      for (int i=0; i<dim1_; ++i) {
        for (int j=0; j<dim2_; ++j) {
             int ij = j + (i*dim2_);
             A[ij] = A2d_[i][j];
        }
      }
}//

void Tensor2d::mgs()
{
    double rmgs1,rmgs2;
      //#pragma omp parallel for
      for (int k=0; k<dim1_; k++) {// loop-1
	rmgs1 = 0.0;

	for (int i=0; i<dim1_;i++) {
	  rmgs1 += A2d_[i][k] * A2d_[i][k];
	}

	rmgs1 = sqrt(rmgs1);

	for (int i=0; i<dim1_;i++) {
	  A2d_[i][k]/=rmgs1;
	}

	for (int j=(k+1); j<dim1_;j++) {// loop-2
	  rmgs2 = 0.0;

	  for (int i=0; i<dim1_;i++) {
	    rmgs2 += A2d_[i][k] * A2d_[i][j];
	  }

	  for (int i=0; i<dim1_;i++) {
	    A2d_[i][j] -= rmgs2 * A2d_[i][k];
	  }
	}// end 2
      }// end 1

}//

void Tensor2d::gs()
{
    if (dim1_ != 0 && dim2_ != 0 ) {
       schmidt(A2d_, dim1_, dim2_, "outfile");
    }
}//

double *Tensor2d::row_vector(int n)
{
  double *temp = new double[dim2_];
  memset(temp, 0, dim2_ * sizeof(double));
  for (int i=0; i<dim2_; i++) temp[i] = A2d_[n][i];
  return temp;
}//

double *Tensor2d::column_vector(int n)
{
  double *temp = new double[dim1_];
  memset(temp, 0, dim1_ * sizeof(double));
  for (int i=0; i<dim1_; i++) temp[i] = A2d_[i][n];
  return temp;
}//

void Tensor2d::sort(int sort_type, const SharedTensor2d &A, double alpha, double beta)
{

  int d1 = A->d1_;
  int d2 = A->d2_;
  int d3 = A->d3_;
  int d4 = A->d4_;

 if (sort_type == 1243) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int sr = col_idx_[s][r];
                        A2d_[pq][sr] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[pq][sr]);
                   }
              }
         }
    }
 }

 else if (sort_type == 1324) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int pr = row_idx_[p][r];
                        int qs = col_idx_[q][s];
                        A2d_[pr][qs] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[pr][qs]);
                   }
              }
         }
    }
 }

 else if (sort_type == 1342) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int pr = row_idx_[p][r];
                        int sq = col_idx_[s][q];
                        A2d_[pr][sq] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[pr][sq]);
                   }
              }
         }
    }
 }

 else if (sort_type == 1423) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int ps = row_idx_[p][s];
                        int qr = col_idx_[q][r];
                        A2d_[ps][qr] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[ps][qr]);
                   }
              }
         }
    }
 }

 else if (sort_type == 1432) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   int rq = col_idx_[r][q];
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int ps = row_idx_[p][s];
                        A2d_[ps][rq] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[ps][rq]);
                   }
              }
         }
    }
 }

 else if (sort_type == 2134) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int qp = row_idx_[q][p];
                        A2d_[qp][rs] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[qp][rs]);
                   }
              }
         }
    }
 }

 else if (sort_type == 2143) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int qp = row_idx_[q][p];
                        int sr = col_idx_[s][r];
                        A2d_[qp][sr] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[qp][sr]);
                   }
              }
         }
    }
 }

 else if (sort_type == 2314) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int qr = row_idx_[q][r];
                        int ps = col_idx_[p][s];
                        A2d_[qr][ps] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[qr][ps]);
                   }
              }
         }
    }
 }

 else if (sort_type == 2341) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int qr = row_idx_[q][r];
                        int sp = col_idx_[s][p];
                        A2d_[qr][sp] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[qr][sp]);
                   }
              }
         }
    }
 }

 else if (sort_type == 2413) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int qs = row_idx_[q][s];
                        int pr = col_idx_[p][r];
                        A2d_[qs][pr] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[qs][pr]);
                   }
              }
         }
    }
 }

 else if (sort_type == 2431) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int qs = row_idx_[q][s];
                        int rp = col_idx_[r][p];
                        A2d_[qs][rp] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[qs][rp]);
                   }
              }
         }
    }
 }

 else if (sort_type == 3124) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int rp = row_idx_[r][p];
                        int qs = col_idx_[q][s];
                        A2d_[rp][qs] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[rp][qs]);
                   }
              }
         }
    }
 }

 else if (sort_type == 3142) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int rp = row_idx_[r][p];
                        int sq = col_idx_[s][q];
                        A2d_[rp][sq] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[rp][sq]);
                   }
              }
         }
    }
 }

 else if (sort_type == 3214) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   int rq = row_idx_[r][q];
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int ps = col_idx_[p][s];
                        A2d_[rq][ps] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[rq][ps]);
                   }
              }
         }
    }
 }

 else if (sort_type == 3241) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int rq = row_idx_[r][q];
                        int sp = col_idx_[s][p];
                        A2d_[rq][sp] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[rq][sp]);
                   }
              }
         }
    }
 }

 else if (sort_type == 3412) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        A2d_[rs][pq] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[rs][pq]);
                   }
              }
         }
    }
 }

 else if (sort_type == 3421) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int qp = col_idx_[q][p];
                        A2d_[rs][qp] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[rs][qp]);
                   }
              }
         }
    }
 }

 else if (sort_type == 4123) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int sp = row_idx_[s][p];
                        int qr = col_idx_[q][r];
                        A2d_[sp][qr] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[sp][qr]);
                   }
              }
         }
    }
 }

 else if (sort_type == 4132) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int sp = row_idx_[s][p];
                        int rq = col_idx_[r][q];
                        A2d_[sp][rq] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[sp][rq]);
                   }
              }
         }
    }
 }

 else if (sort_type == 4213) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int sq = row_idx_[s][q];
                        int pr = col_idx_[p][r];
                        A2d_[sq][pr] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[sq][pr]);
                   }
              }
         }
    }
 }

 else if (sort_type == 4231) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int sq = row_idx_[s][q];
                        int rp = col_idx_[r][p];
                        A2d_[sq][rp] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[sq][rp]);
                   }
              }
         }
    }
 }

 else if (sort_type == 4312) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int sr = row_idx_[s][r];
                        A2d_[sr][pq] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[sr][pq]);
                   }
              }
         }
    }
 }

 else if (sort_type == 4321) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = A->row_idx_[p][q];
              for (int r = 0; r < d3; r++) {
                   for (int s = 0; s < d4; s++) {
                        int rs = A->col_idx_[r][s];
                        int sr = row_idx_[s][r];
                        int qp = col_idx_[q][p];
                        A2d_[sr][qp] = (alpha*A->A2d_[pq][rs]) + (beta*A2d_[sr][qp]);
                   }
              }
         }
    }
 }

 else {
    outfile->Printf("\tUnrecognized sort type!\n");
    throw PSIEXCEPTION("Unrecognized sort type!");

 }

}//

void Tensor2d::sort3a(int sort_type, int d1, int d2, int d3, const SharedTensor2d &A, double alpha, double beta)
{

 if (sort_type == 132) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              for (int r = 0; r < d3; r++) {
                   int rq = q + (r*d2);
                   int qr = r + (q*d3);
                   A2d_[p][rq] = (alpha*A->A2d_[p][qr]) + (beta*A2d_[p][rq]);
              }
         }
    }
 }

 else {
    outfile->Printf("\tUnrecognized sort type!\n");
    throw PSIEXCEPTION("Unrecognized sort type!");

 }

}//

void Tensor2d::sort3b(int sort_type, int d1, int d2, int d3, const SharedTensor2d &A, double alpha, double beta)
{

 if (sort_type == 132) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = q + (p*d2);
              for (int r = 0; r < d3; r++) {
                   int pr = r + (p*d3);
                   A2d_[pr][q] = (alpha*A->A2d_[pq][r]) + (beta*A2d_[pr][q]);
              }
         }
    }
 }

 else if (sort_type == 213) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = q + (p*d2);
              int qp = p + (q*d1);
              for (int r = 0; r < d3; r++) {
                   A2d_[qp][r] = (alpha*A->A2d_[pq][r]) + (beta*A2d_[qp][r]);
              }
         }
    }
 }

 else if (sort_type == 312) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = q + (p*d2);
              for (int r = 0; r < d3; r++) {
                   int rp = p + (r*d1);
                   A2d_[rp][q] = (alpha*A->A2d_[pq][r]) + (beta*A2d_[rp][q]);
              }
         }
    }
 }

 else if (sort_type == 231) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = q + (p*d2);
              for (int r = 0; r < d3; r++) {
                   int qr = r + (q*d3);
                   A2d_[qr][p] = (alpha*A->A2d_[pq][r]) + (beta*A2d_[qr][p]);
              }
         }
    }
 }

 else if (sort_type == 321) {
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
         for (int q = 0; q < d2; q++) {
              int pq = q + (p*d2);
              for (int r = 0; r < d3; r++) {
                   int rq = q + (r*d2);
                   A2d_[rq][p] = (alpha*A->A2d_[pq][r]) + (beta*A2d_[rq][p]);
              }
         }
    }
 }

 else {
    outfile->Printf("\tUnrecognized sort type!\n");
    throw PSIEXCEPTION("Unrecognized sort type!");

 }

}//

void Tensor2d::apply_denom(int frzc, int occ, const SharedTensor2d &fock)
{
    int aocc = d1_;
    int avir = d3_;

    #pragma omp parallel for
    for (int i = 0; i < aocc; i++) {
         double di = fock->A2d_[i + frzc][i + frzc];
         for (int j = 0; j < aocc; j++) {
              double dij = di + fock->A2d_[j + frzc][j + frzc];
              int ij = row_idx_[i][j];
              for (int a = 0; a < avir; a++) {
                   double dija = dij - fock->A2d_[a + occ][a + occ];
                   for (int b = 0; b < avir; b++) {
                        double dijab = dija - fock->A2d_[b + occ][b + occ];
                        int ab = col_idx_[a][b];
                        A2d_[ij][ab] /= dijab;
                   }
              }
         }
    }
}//

void Tensor2d::apply_denom_os(int frzc, int occA, int occB, const SharedTensor2d &fockA, const SharedTensor2d &fockB)
{
    int aoccA = d1_;
    int aoccB = d2_;
    int avirA = d3_;
    int avirB = d4_;

    #pragma omp parallel for
    for (int i = 0; i < aoccA; i++) {
         double di = fockA->A2d_[i + frzc][i + frzc];
         for (int j = 0; j < aoccB; j++) {
              double dij = di + fockB->A2d_[j + frzc][j + frzc];
              int ij = row_idx_[i][j];
              for (int a = 0; a < avirA; a++) {
                   double dija = dij - fockA->A2d_[a + occA][a + occA];
                   for (int b = 0; b < avirB; b++) {
                        double dijab = dija - fockB->A2d_[b + occB][b + occB];
                        int ab = col_idx_[a][b];
                        A2d_[ij][ab] /= dijab;
                   }
              }
         }
    }
}//

void Tensor2d::apply_denom_chem(int frzc, int occ, const SharedTensor2d &fock)
{
    int aocc = d1_;
    int avir = d2_;

    #pragma omp parallel for
    for (int i = 0; i < aocc; i++) {
         double di = fock->A2d_[i + frzc][i + frzc];
         for (int a = 0; a < avir; a++) {
              double dia = di - fock->A2d_[a + occ][a + occ];
              int ia = row_idx_[i][a];
              for (int j = 0; j < aocc; j++) {
                   double diaj = dia + fock->A2d_[j + frzc][j + frzc];
                   for (int b = 0; b < avir; b++) {
                        double diajb = diaj - fock->A2d_[b + occ][b + occ];
                        int jb = col_idx_[j][b];
                        A2d_[ia][jb] /= diajb;
                   }
              }
         }
    }
}//

void Tensor2d::reg_denom(int frzc, int occ, const SharedTensor2d &fock, double reg)
{
    int aocc = d1_;
    int avir = d3_;

    #pragma omp parallel for
    for (int i = 0; i < aocc; i++) {
         double di = fock->A2d_[i + frzc][i + frzc] - reg;
         for (int j = 0; j < aocc; j++) {
              double dij = di + fock->A2d_[j + frzc][j + frzc];
              int ij = row_idx_[i][j];
              for (int a = 0; a < avir; a++) {
                   double dija = dij - fock->A2d_[a + occ][a + occ];
                   for (int b = 0; b < avir; b++) {
                        double dijab = dija - fock->A2d_[b + occ][b + occ];
                        int ab = col_idx_[a][b];
                        A2d_[ij][ab] /= dijab;
                   }
              }
         }
    }
}//

void Tensor2d::reg_denom_os(int frzc, int occA, int occB, const SharedTensor2d &fockA, const SharedTensor2d &fockB, double reg)
{
    int aoccA = d1_;
    int aoccB = d2_;
    int avirA = d3_;
    int avirB = d4_;

    #pragma omp parallel for
    for (int i = 0; i < aoccA; i++) {
         double di = fockA->A2d_[i + frzc][i + frzc] - reg;
         for (int j = 0; j < aoccB; j++) {
              double dij = di + fockB->A2d_[j + frzc][j + frzc];
              int ij = row_idx_[i][j];
              for (int a = 0; a < avirA; a++) {
                   double dija = dij - fockA->A2d_[a + occA][a + occA];
                   for (int b = 0; b < avirB; b++) {
                        double dijab = dija - fockB->A2d_[b + occB][b + occB];
                        int ab = col_idx_[a][b];
                        A2d_[ij][ab] /= dijab;
                   }
              }
         }
    }
}//

void Tensor2d::reg_denom_chem(int frzc, int occ, const SharedTensor2d &fock, double reg)
{
    int aocc = d1_;
    int avir = d2_;

    #pragma omp parallel for
    for (int i = 0; i < aocc; i++) {
         double di = fock->A2d_[i + frzc][i + frzc] - reg;
         for (int a = 0; a < avir; a++) {
              double dia = di - fock->A2d_[a + occ][a + occ];
              int ia = row_idx_[i][a];
              for (int j = 0; j < aocc; j++) {
                   double diaj = dia + fock->A2d_[j + frzc][j + frzc];
                   for (int b = 0; b < avir; b++) {
                        double diajb = diaj - fock->A2d_[b + occ][b + occ];
                        int jb = col_idx_[j][b];
                        A2d_[ia][jb] /= diajb;
                   }
              }
         }
    }
}//

void Tensor2d::dirprd(const SharedTensor2d &a, const SharedTensor2d &b)
{
     #pragma omp parallel for
     for (int i=0; i<dim1_; i++) {
          for (int j=0; j<dim2_; j++) {
              A2d_[i][j] = a->get(i,j) * b->get(i,j);
          }
     }
}//

void Tensor2d::dirprd123(const SharedTensor1d &a, const SharedTensor2d &b, double alpha, double beta)
{
     int d1 = dim1_;
     int d2 = b->dim1();
     int d3 = b->dim2();
     #pragma omp parallel for
     for (int i=0; i<d1; i++) {
          for (int j=0; j<d2; j++) {
               for (int k=0; k<d3; k++) {
                    int jk = k + (j*d3);
                    //A2d_[i][jk] = a->get(i) * b->get(j,k);
                    A2d_[i][jk] = (alpha * a->get(i) * b->get(j,k)) + (beta * A2d_[i][jk]);
               }
          }
     }
}//

void Tensor2d::dirprd123(bool transb, const SharedTensor1d &a, const SharedTensor2d &b, double alpha, double beta)
{
   if (transb) {
     int d1 = dim1_;
     int d2 = b->dim2();
     int d3 = b->dim1();
     #pragma omp parallel for
     for (int i=0; i<d1; i++) {
          for (int j=0; j<d2; j++) {
               for (int k=0; k<d3; k++) {
                    int jk = k + (j*d3);
                    A2d_[i][jk] = (alpha * a->get(i) * b->get(k,j)) + (beta * A2d_[i][jk]);
               }
          }
     }
   }

   else {
     int d1 = dim1_;
     int d2 = b->dim1();
     int d3 = b->dim2();
     #pragma omp parallel for
     for (int i=0; i<d1; i++) {
          for (int j=0; j<d2; j++) {
               for (int k=0; k<d3; k++) {
                    int jk = k + (j*d3);
                    //A2d_[i][jk] = a->get(i) * b->get(j,k);
                    A2d_[i][jk] = (alpha * a->get(i) * b->get(j,k)) + (beta * A2d_[i][jk]);
               }
          }
     }
   }
}//

void Tensor2d::dirprd112(const SharedTensor1d &a, const SharedTensor1d &b)
{
     #pragma omp parallel for
     for (int i=0; i<dim1_; i++) {
          for (int j=0; j<dim2_; j++) {
               A2d_[i][j] = a->get(i) * b->get(j);
          }
     }
}//

void Tensor2d::dirprd112(const SharedTensor1d &a, const SharedTensor1d &b, double alpha, double beta)
{
     #pragma omp parallel for
     for (int i=0; i<dim1_; i++) {
          for (int j=0; j<dim2_; j++) {
               A2d_[i][j] = (alpha * a->get(i) * b->get(j)) + (beta * A2d_[i][j]);
          }
     }
}//

void Tensor2d::dirprd224(const SharedTensor2d &a, const SharedTensor2d &b)
{
     #pragma omp parallel for
     for (int i=0; i<d1_; i++) {
          for (int j=0; j<d2_; j++) {
               int ij = row_idx_[i][j];
               for (int k=0; k<d3_; k++) {
                    for (int l=0; l<d4_; l++) {
                         int kl = col_idx_[k][l];
                         A2d_[ij][kl] = a->get(i,j) * b->get(k,l);
                    }
               }
          }
     }
}//

void Tensor2d::dirprd224(const SharedTensor2d &a, const SharedTensor2d &b, double alpha, double beta)
{
     #pragma omp parallel for
     for (int i=0; i<d1_; i++) {
          for (int j=0; j<d2_; j++) {
               int ij = row_idx_[i][j];
               for (int k=0; k<d3_; k++) {
                    for (int l=0; l<d4_; l++) {
                         int kl = col_idx_[k][l];
                         A2d_[ij][kl] = (alpha * a->get(i,j) * b->get(k,l)) + (beta * A2d_[ij][kl]);
                    }
               }
          }
     }
}//

double* Tensor2d::to_vector(const SharedTensor2i &pair_idx)
{
     double* temp = new double[dim1_ * dim2_];
     #pragma omp parallel for
     for (int i=0; i<dim1_; i++) {
          for (int j=0; j<dim2_; j++) {
               int ij = pair_idx->get(i,j);
               temp[ij] = A2d_[i][j];
          }
     }
     return temp;
}//

double* Tensor2d::to_vector()
{
     double* temp = new double[dim1_ * dim2_];
     #pragma omp parallel for
     for (int i=0; i<dim1_; i++) {
          for (int j=0; j<dim2_; j++) {
               int ij = (i * dim2_) + j;
               temp[ij] = A2d_[i][j];
          }
     }
     return temp;
}//

double Tensor2d::rms()
{
  double summ = 0.0;
  #pragma omp parallel for
  for (int i=0; i<dim1_; ++i) {
       for (int j=0; j<dim2_; ++j) {
            summ += A2d_[i][j] * A2d_[i][j];
       }
  }
  summ=sqrt(summ/(dim1_*dim2_));

  return summ;
}//

double Tensor2d::rms(const SharedTensor2d& a)
{
  double summ = 0.0;
  #pragma omp parallel for
  for (int i=0; i<dim1_; ++i) {
       for (int j=0; j<dim2_; ++j) {
            summ += (A2d_[i][j] - a->A2d_[i][j]) * (A2d_[i][j] - a->A2d_[i][j]);
       }
  }
  summ=sqrt(summ/(dim1_*dim2_));

  return summ;
}//

void Tensor2d::set_act_oo(int aocc, const SharedTensor2d &a)
{
    #pragma omp parallel for
    for (int i = 0; i < aocc; i++) {
         for (int j = 0; j < aocc; j++) {
              A2d_[i][j] = a->get(i,j);
         }
    }
}//

void Tensor2d::set_act_oo(int frzc, int aocc, const SharedTensor2d &a)
{
    #pragma omp parallel for
    for (int i = 0; i < aocc; i++) {
         for (int j = 0; j < aocc; j++) {
              A2d_[i + frzc][j + frzc] = a->get(i,j);
         }
    }
}//

void Tensor2d::set_oo(const SharedTensor2d &a)
{
    int occ = a->dim1();
    #pragma omp parallel for
    for (int i = 0; i < occ; i++) {
         for (int j = 0; j < occ; j++) {
              A2d_[i][j] = a->get(i,j);
         }
    }
}//

void Tensor2d::set_act_ov(int frzc, const SharedTensor2d &A)
{
    int aocc = A->dim1();
    int avir = A->dim2();
    #pragma omp parallel for
    for (int i = 0; i < aocc; i++) {
         for (int a = 0; a < avir; a++) {
              A2d_[i + frzc][a] = A->get(i,a);
         }
    }
}//

void Tensor2d::set_act_vo(int frzc, const SharedTensor2d &A)
{
    int avir = A->dim1();
    int aocc = A->dim2();
    #pragma omp parallel for
    for (int a = 0; a < avir; a++) {
         for (int i = 0; i < aocc; i++) {
              A2d_[a][i + frzc] = A->get(a,i);
         }
    }
}//

void Tensor2d::set_act_vv(int occ, int avir, const SharedTensor2d &A)
{
    #pragma omp parallel for
    for (int a = 0; a < avir; a++) {
         for (int b = 0; b < avir; b++) {
              A2d_[a + occ][b + occ] = A->get(a,b);
         }
    }
}//

void Tensor2d::set_act_vv(const SharedTensor2d &A)
{
    int avir = A->dim1();
    #pragma omp parallel for
    for (int a = 0; a < avir; a++) {
         for (int b = 0; b < avir; b++) {
              A2d_[a][b] = A->get(a,b);
         }
    }
}//

void Tensor2d::set_vv(int occ, const SharedTensor2d &A)
{
    int vir = A->dim1();
    #pragma omp parallel for
    for (int a = 0; a < vir; a++) {
         for (int b = 0; b < vir; b++) {
              A2d_[a + occ][b + occ] = A->get(a,b);
         }
    }
}//

void Tensor2d::set_ov(const SharedTensor2d &A)
{
    int occ = A->dim1();
    int vir = A->dim2();
    #pragma omp parallel for
    for (int i = 0; i < occ; i++) {
         for (int a = 0; a < vir; a++) {
              A2d_[i][a + occ] = A->get(i,a);
         }
    }
}//

void Tensor2d::set_vo(const SharedTensor2d &A)
{
    int vir = A->dim1();
    int occ = A->dim2();
    #pragma omp parallel for
    for (int a = 0; a < vir; a++) {
         for (int i = 0; i < occ; i++) {
              A2d_[a + occ][i] = A->get(a,i);
         }
    }
}//

void Tensor2d::add_oo(const SharedTensor2d &A, double alpha, double beta)
{
    int occ = A->dim1();
    #pragma omp parallel for
    for (int i = 0; i < occ; i++) {
         for (int j = 0; j < occ; j++) {
              A2d_[i][j] = (alpha * A->get(i,j)) + (beta * A2d_[i][j]);
         }
    }
}//

void Tensor2d::add_vv(int occ, const SharedTensor2d &A, double alpha, double beta)
{
    int vir = A->dim1();
    #pragma omp parallel for
    for (int a = 0; a < vir; a++) {
         for (int b = 0; b < vir; b++) {
              A2d_[a + occ][b + occ] = (alpha *A->get(a,b)) + (beta * A2d_[a + occ][b + occ]);
         }
    }
}//

void Tensor2d::add_ov(const SharedTensor2d &A, double alpha, double beta)
{
    int occ = A->dim1();
    int vir = A->dim2();
    #pragma omp parallel for
    for (int i = 0; i < occ; i++) {
         for (int a = 0; a < vir; a++) {
              A2d_[i][a + occ] = (alpha * A->get(i,a)) + (beta * A2d_[i][a + occ]);
         }
    }
}//

void Tensor2d::add_vo(const SharedTensor2d &A, double alpha, double beta)
{
    int vir = A->dim1();
    int occ = A->dim2();
    #pragma omp parallel for
    for (int a = 0; a < vir; a++) {
         for (int i = 0; i < occ; i++) {
              A2d_[a + occ][i] = (alpha * A->get(a,i)) + (beta * A2d_[a + occ][i]);
         }
    }
}//

void Tensor2d::add_aocc_fc(const SharedTensor2d &A, double alpha, double beta)
{
    int aocc = A->dim1();
    int frzc = A->dim2();
    #pragma omp parallel for
    for (int i = 0; i < aocc; i++) {
         for (int j = 0; j < frzc; j++) {
              A2d_[i + frzc][j] = (alpha * A->get(i,j)) + (beta * A2d_[i + frzc][j]);
         }
    }
}//

void Tensor2d::add_fc_aocc(const SharedTensor2d &A, double alpha, double beta)
{
    int frzc = A->dim1();
    int aocc = A->dim2();
    #pragma omp parallel for
    for (int i = 0; i < frzc; i++) {
         for (int j = 0; j < aocc; j++) {
              A2d_[i][j + frzc] = (alpha * A->get(i,j)) + (beta * A2d_[i][j + frzc]);
         }
    }
}//

void Tensor2d::set3_oo(const SharedTensor2d &A)
{
    int naux = A->d1_;
    int occ = A->d2_;
    #pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
         for (int i = 0; i < occ; i++) {
              for (int j = 0; j < occ; j++) {
                   int ij = A->col_idx_[i][j];
                   int oo = col_idx_[i][j];
                   A2d_[Q][oo] = A->get(Q,ij);
              }
         }
    }
}//

void Tensor2d::set3_act_oo(int frzc, const SharedTensor2d &A)
{
    int naux = A->d1_;
    int aoccA = A->d2_;
    int aoccB = A->d3_;
    int occA = d2_;
    int occB = d3_;
    #pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
         for (int i = 0; i < aoccA; i++) {
              for (int j = 0; j < aoccB; j++) {
                   int ij = A->col_idx_[i][j];
                   int oo = ( (i + frzc) * occB) + j + frzc;
                   A2d_[Q][oo] = A->get(Q,ij);
              }
         }
    }
}//

void Tensor2d::add3_oo(const SharedTensor2d &A, double alpha, double beta)
{
    int naux = A->d1_;
    int occ = A->d2_;
    #pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
         for (int i = 0; i < occ; i++) {
              for (int j = 0; j < occ; j++) {
                   int ij = A->col_idx_[i][j];
                   int oo = col_idx_[i][j];
                   A2d_[Q][oo] = (alpha * A->get(Q,ij)) + (beta * A2d_[Q][oo]);
              }
         }
    }
}//

void Tensor2d::set3_act_ov(int frzc, int aocc, int avir, int vir, const SharedTensor2d &A)
{
    int naux = dim1_;
    #pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
         for (int i = 0; i < aocc; i++) {
              for (int a = 0; a < avir; a++) {
                   int ia = (i * avir) + a;
                   int ov = ( (i + frzc) * vir) + a;
                   A2d_[Q][ov] = A->get(Q,ia);
              }
         }
    }
}//

void Tensor2d::set3_ov(const SharedTensor2d &A)
{
    int naux = dim1_;
    int occ = A->d2_;
    int vir = A->d3_;
    #pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
         for (int i = 0; i < occ; i++) {
              for (int a = 0; a < vir; a++) {
                   int ia = A->col_idx_[i][a];
                   int ov = col_idx_[i][a + occ];
                   A2d_[Q][ov] = A->get(Q,ia);
              }
         }
    }
}//

void Tensor2d::set3_vo(const SharedTensor2d &A)
{
    int naux = dim1_;
    int vir = A->d2_;
    int occ = A->d3_;
    #pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
         for (int a = 0; a < vir; a++) {
              for (int i = 0; i < occ; i++) {
                   int ai = A->col_idx_[a][i];
                   int vo = col_idx_[a + occ][i];
                   A2d_[Q][vo] = A->get(Q,ai);
              }
         }
    }
}//

void Tensor2d::set3_vv(const SharedTensor2d &A, int occ)
{
    int naux = dim1_;
    int vir = A->d2_;
    #pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
         for (int a = 0; a < vir; a++) {
              for (int b = 0; b < vir; b++) {
                   int ab = A->col_idx_[a][b];
                   int vv = col_idx_[a + occ][b + occ];
                   A2d_[Q][vv] = A->get(Q,ab);
              }
         }
    }
}//

void Tensor2d::set3_act_vv(const SharedTensor2d &A)
{
    int naux = dim1_;
    int avir = A->d2_;
    #pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
         for (int a = 0; a < avir; a++) {
              for (int b = 0; b < avir; b++) {
                   int ab = A->col_idx_[a][b];
                   int vv = col_idx_[a][b];
                   A2d_[Q][vv] = A->get(Q,ab);
              }
         }
    }
}//

void Tensor2d::swap_3index_col(const SharedTensor2d &A)
{

  int d1 = A->d1_;
  int d2 = A->d2_;
  int d3 = A->d3_;

    #pragma omp parallel for
    for (int Q = 0; Q < d1; Q++) {
         for (int p = 0; p < d2; p++) {
              for (int q = 0; q < d3; q++) {
                   int pq = A->col_idx_[p][q];
                   int qp = col_idx_[q][p];
                   A2d_[Q][qp] = A->A2d_[Q][pq];
              }
         }
    }
}//

void Tensor2d::form_oo(const SharedTensor2d &A)
{
    int occ = dim1_;
    #pragma omp parallel for
    for (int i = 0; i < occ; i++) {
         for (int j = 0; j < occ; j++) {
              A2d_[i][j] = A->get(i,j);
         }
    }
}//

void Tensor2d::form_act_oo(int frzc, const SharedTensor2d &A)
{
    int occ = dim1_;
    #pragma omp parallel for
    for (int i = 0; i < occ; i++) {
         for (int j = 0; j < occ; j++) {
              A2d_[i][j] = A->get(i + frzc, j + frzc);
         }
    }
}//

void Tensor2d::form_vv(int occ, const SharedTensor2d &A)
{
    int vir = dim1_;
    #pragma omp parallel for
    for (int a = 0; a < vir; a++) {
         for (int b = 0; b < vir; b++) {
              A2d_[a][b] = A->get(a + occ, b + occ);
         }
    }
}//

void Tensor2d::form_act_vv(int occ, const SharedTensor2d &A)
{
    int vir = dim1_;
    #pragma omp parallel for
    for (int a = 0; a < vir; a++) {
         for (int b = 0; b < vir; b++) {
              A2d_[a][b] = A->get(a + occ, b + occ);
         }
    }
}//

void Tensor2d::form_vo(const SharedTensor2d &A)
{
    int vir = dim1_;
    int occ = dim2_;
    #pragma omp parallel for
    for (int a = 0; a < vir; a++) {
         for (int i = 0; i < occ; i++) {
              A2d_[a][i] = A->get(a + occ, i);
         }
    }
}//

void Tensor2d::form_vo(int occ, const SharedTensor2d &A)
{
    int vir = dim1_;
    int occ2 = dim2_;
    #pragma omp parallel for
    for (int a = 0; a < vir; a++) {
         for (int i = 0; i < occ2; i++) {
              A2d_[a][i] = A->get(a + occ, i);
         }
    }
}//

void Tensor2d::form_act_vo(int frzc, const SharedTensor2d &A)
{
    int vir = dim1_;
    int occ = dim2_;
    #pragma omp parallel for
    for (int a = 0; a < vir; a++) {
         for (int i = 0; i < occ; i++) {
              A2d_[a][i] = A->get(a + occ, i + frzc);
         }
    }
}//

void Tensor2d::form_act_vo(int frzc, int occ, const SharedTensor2d &A)
{
    int vir = dim1_;
    int occ2 = dim2_;
    #pragma omp parallel for
    for (int a = 0; a < vir; a++) {
         for (int i = 0; i < occ2; i++) {
              A2d_[a][i] = A->get(a + occ, i + frzc);
         }
    }
}//

void Tensor2d::form_ov(const SharedTensor2d &A)
{
    int vir = dim2_;
    int occ = dim1_;
    #pragma omp parallel for
    for (int i = 0; i < occ; i++) {
         for (int a = 0; a < vir; a++) {
              A2d_[i][a] = A->get(i, a + occ);
         }
    }
}//

void Tensor2d::form_ov(int occ, const SharedTensor2d &A)
{
    int vir = dim2_;
    int occ2 = dim1_;
    #pragma omp parallel for
    for (int i = 0; i < occ2; i++) {
         for (int a = 0; a < vir; a++) {
              A2d_[i][a] = A->get(i, a + occ);
         }
    }
}//

void Tensor2d::form_act_ov(int frzc, const SharedTensor2d &A)
{
    int vir = dim2_;
    int occ = dim1_;
    #pragma omp parallel for
    for (int i = 0; i < occ; i++) {
         for (int a = 0; a < vir; a++) {
              A2d_[i][a] = A->get(i + frzc, a + occ);
         }
    }
}//

void Tensor2d::form_act_ov(int frzc, int occ, const SharedTensor2d &A)
{
    int vir = dim2_;
    int occ2 = dim1_;
    #pragma omp parallel for
    for (int i = 0; i < occ2; i++) {
         for (int a = 0; a < vir; a++) {
              A2d_[i][a] = A->get(i + frzc, a + occ);
         }
    }
}//

void Tensor2d::form_ooAB(const SharedTensor2d &A)
{
    int occA = dim1_;
    int occB = dim2_;
    #pragma omp parallel for
    for (int i = 0; i < occA; i++) {
         for (int j = 0; j < occB; j++) {
              A2d_[i][j] = A->get(i,j);
         }
    }
}//

void Tensor2d::form_b_ij(int frzc, const SharedTensor2d &A)
{
    #pragma omp parallel for
    for (int Q = 0; Q < d1_; Q++) {
         for (int i = 0; i < d2_; i++) {
              for (int j = 0; j < d3_; j++) {
                   int ij = col_idx_[i][j];
                   int oo = A->col_idx_[i + frzc][j + frzc];
                   A2d_[Q][ij] = A->get(Q,oo);
              }
         }
    }
}//

void Tensor2d::form_b_ia(int frzc, const SharedTensor2d &A)
{
    #pragma omp parallel for
    for (int Q = 0; Q < d1_; Q++) {
         for (int i = 0; i < d2_; i++) {
              for (int a = 0; a < d3_; a++) {
                   int ia = col_idx_[i][a];
                   int ov = A->col_idx_[i + frzc][a];
                   A2d_[Q][ia] = A->get(Q,ov);
              }
         }
    }
}//

void Tensor2d::form_b_ab(const SharedTensor2d &A)
{
    #pragma omp parallel for
    for (int Q = 0; Q < d1_; Q++) {
         for (int a = 0; a < d2_; a++) {
              for (int b = 0; b < d3_; b++) {
                   int ab = col_idx_[a][b];
                   int vv = A->col_idx_[a][b];
                   A2d_[Q][ab] = A->get(Q,vv);
              }
         }
    }
}//

void Tensor2d::form_b_kl(const SharedTensor2d &A)
{
    int naux = d1_;
    int aocc = d2_;
    int frzc = d3_;
    #pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
         for (int i = 0; i < aocc; i++) {
              for (int j = 0; j < frzc; j++) {
                   int ij = A->col_idx_[i + frzc][j];
                   int oo = col_idx_[i][j];
                   A2d_[Q][oo] = A->get(Q,ij);
              }
         }
    }
}//

void Tensor2d::form_b_ki(const SharedTensor2d &A)
{
    int naux = d1_;
    int aocc = d2_;
    int occ = d3_;
    int frzc = occ - aocc;
    #pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
         for (int i = 0; i < aocc; i++) {
              for (int j = 0; j < occ; j++) {
                   int ij = A->col_idx_[i + frzc][j];
                   int oo = col_idx_[i][j];
                   A2d_[Q][oo] = A->get(Q,ij);
              }
         }
    }
}//

void Tensor2d::form_b_ka(const SharedTensor2d &A)
{
    int naux = d1_;
    int aocc = d2_;
    int vir = d3_;
    int occ = A->d2_;
    int frzc = occ - aocc;
    #pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
         for (int i = 0; i < aocc; i++) {
              for (int j = 0; j < vir; j++) {
                   int ij = A->col_idx_[i + frzc][j];
                   int oo = col_idx_[i][j];
                   A2d_[Q][oo] = A->get(Q,ij);
              }
         }
    }
}//

void Tensor2d::form_b_li(const SharedTensor2d &A)
{
    int naux = d1_;
    int frzc = d2_;
    int occ = d3_;
    #pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
         for (int i = 0; i < frzc; i++) {
              for (int j = 0; j < occ; j++) {
                   int ij = A->col_idx_[i][j];
                   int oo = col_idx_[i][j];
                   A2d_[Q][oo] = A->get(Q,ij);
              }
         }
    }
}//

void Tensor2d::form_b_il(const SharedTensor2d &A)
{
    int naux = d1_;
    int occ = d2_;
    int frzc = d3_;
    #pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
         for (int i = 0; i < occ; i++) {
              for (int j = 0; j < frzc; j++) {
                   int ij = A->col_idx_[i][j];
                   int oo = col_idx_[i][j];
                   A2d_[Q][oo] = A->get(Q,ij);
              }
         }
    }
}//

void Tensor2d::form_b_la(const SharedTensor2d &A)
{
    int naux = d1_;
    int frzc = d2_;
    int vir = d3_;
    #pragma omp parallel for
    for (int Q = 0; Q < naux; Q++) {
         for (int i = 0; i < frzc; i++) {
              for (int j = 0; j < vir; j++) {
                   int ij = A->col_idx_[i][j];
                   int oo = col_idx_[i][j];
                   A2d_[Q][oo] = A->get(Q,ij);
              }
         }
    }
}//

void Tensor2d::symmetrize()
{
    SharedTensor2d temp = SharedTensor2d(new Tensor2d(dim2_, dim1_));
    temp = transpose();
    add(temp);
    scale(0.5);
    temp.reset();
    /*
    #pragma omp parallel for
    for (int i=0; i<dim2_; ++i) {
	for (int j=0; j<dim1_; ++j) {
	     A2d_[i][j] = 0.5 * (A2d_[i][j] + A2d_[j][i]);
	}
    }
    */

}//

void Tensor2d::symmetrize(const SharedTensor2d &A)
{
    #pragma omp parallel for
    for (int i=0; i<dim1_; ++i) {
	for (int j=0; j<dim2_; ++j) {
	     A2d_[i][j] = 0.5 * (A->A2d_[i][j] + A->A2d_[j][i]);
	}
    }

}//

void Tensor2d::symmetrize3(const SharedTensor2d &A)
{
    SharedTensor2d temp = SharedTensor2d(new Tensor2d("temp", d1_, d3_, d2_));
    temp->swap_3index_col(A);
    add(temp);
    scale(0.5);
    temp.reset();
}//

void Tensor2d::symm_packed(const SharedTensor2d &A)
{
    // Form symetric packed 3-index
    #pragma omp parallel for
    for (int R = 0; R < A->d1_; R++) {
          for (int p = 0; p < A->d2_; p++) {
               for (int q = 0; q <= p; q++) {
                    int pq = A->col_idx_[p][q];
                    int pq_sym = index2(p,q);
                    double perm = (p == q ? 1.0 : 2.0);
                    A2d_[R][pq_sym] = perm * A->get(R, pq);
               }
          }
    }

}//

void Tensor2d::ltm(const SharedTensor2d &A)
{
    // Form Lower triangular part
    #pragma omp parallel for
    for (int R = 0; R < A->d1_; R++) {
          for (int p = 0; p < A->d2_; p++) {
               for (int q = 0; q <= p; q++) {
                    int pq = A->col_idx_[p][q];
                    int pq_sym = index2(p,q);
                    A2d_[R][pq_sym] = A->get(R, pq);
               }
          }
    }

}//

void Tensor2d::expand23(int d1, int d2, int d3, const SharedTensor2d &A)
{
    // Convert Lower triangular to full tensor
    #pragma omp parallel for
    for (int p = 0; p < d1; p++) {
          for (int q = 0; q < d2; q++) {
               for (int r = 0; r < d3; r++) {
		    int pq = (p*d2) + q;
                    int qr = index2(q,r);
                    A2d_[pq][r] = A->get(p,qr);
               }
          }
    }

}//

void Tensor2d::set_row(const SharedTensor2d &A, int n)
{
  #pragma omp parallel for
  for (int i=0; i<d3_; i++) {
       for (int j=0; j<d4_; j++) {
           int ij = col_idx_[i][j];
           A2d_[n][ij] = A->get(i, j);
       }
  }
}//

void Tensor2d::set_column(const SharedTensor2d &A, int n)
{
  #pragma omp parallel for
  for (int i=0; i<d1_; i++) {
       for (int j=0; j<d2_; j++) {
           int ij = row_idx_[i][j];
           A2d_[ij][n] = A->get(i, j);
       }
  }

  /*
    ULI length;
    length = (ULI)dim1_ * (ULI)dim2_;
    C_DCOPY(length, A->A2d_[0], 1, &(A2d_[0][n]), 1);
  */
}//

void Tensor2d::get_row(const SharedTensor2d &A, int n)
{
  #pragma omp parallel for
  for (int i=0; i<A->d3_; i++) {
       for (int j=0; j<A->d4_; j++) {
           int ij = A->col_idx_[i][j];
           A2d_[i][j] = A->get(n, ij);
       }
  }
}//

void Tensor2d::get_column(const SharedTensor2d &A, int n)
{
  #pragma omp parallel for
  for (int i=0; i<A->d1_; i++) {
       for (int j=0; j<A->d2_; j++) {
           int ij = A->row_idx_[i][j];
           A2d_[i][j] = A->get(ij, n);
       }
  }

  /*
    ULI length;
    length = (ULI)dim1_ * (ULI)dim2_;
    C_DCOPY(length, &(A->A2d_[0][n]), 1, A2d_[0], 1);
  */
}//

void Tensor2d::add2row(const SharedTensor2d &A, int n)
{
  #pragma omp parallel for
  for (int i=0; i<d3_; i++) {
       for (int j=0; j<d4_; j++) {
           int ij = col_idx_[i][j];
           A2d_[n][ij] += A->get(i, j);
       }
  }
}//

void Tensor2d::add2col(const SharedTensor2d &A, int n)
{
  #pragma omp parallel for
  for (int i=0; i<d1_; i++) {
       for (int j=0; j<d2_; j++) {
           int ij = row_idx_[i][j];
           A2d_[ij][n] += A->get(i, j);
       }
  }
}//

void Tensor2d::symm4(const SharedTensor2d &a)
{
     #pragma omp parallel for
     for (int i=0; i < a->d1_; i++) {
          for (int j=0; j <= i; j++) {
               int ij = a->row_idx_[i][j];
               int ji = a->row_idx_[j][i];
               int ij2 = index2(i,j);
               for (int k=0; k < a->d3_; k++) {
                    for (int l=0; l <= k; l++) {
                         int kl = a->col_idx_[k][l];
                         int kl2 = index2(k,l);
                         A2d_[ij2][kl2] = 0.5 * ( a->get(ij,kl) + a->get(ji,kl) );
                    }
               }
          }
     }
}//

void Tensor2d::symm_col4(const SharedTensor2d &a)
{
     #pragma omp parallel for
     for (int i=0; i < a->d1_; i++) {
          for (int j=0; j <= i; j++) {
               int ij = a->row_idx_[i][j];
               int ij2 = index2(i,j);
               for (int k=0; k < a->d3_; k++) {
                    for (int l=0; l <= k; l++) {
                         int kl = a->col_idx_[k][l];
                         int lk = a->col_idx_[l][k];
                         int kl2 = index2(k,l);
                         A2d_[ij2][kl2] = 0.5 * ( a->get(ij,kl) + a->get(ij,lk) );
                    }
               }
          }
     }
}//

void Tensor2d::antisymm4(const SharedTensor2d &a)
{
     #pragma omp parallel for
     for (int i=0; i < a->d1_; i++) {
          for (int j=0; j <= i; j++) {
               int ij = a->row_idx_[i][j];
               int ji = a->row_idx_[j][i];
               int ij2 = index2(i,j);
               for (int k=0; k < a->d3_; k++) {
                    for (int l=0; l <= k; l++) {
                         int kl = a->col_idx_[k][l];
                         int kl2 = index2(k,l);
                         A2d_[ij2][kl2] = 0.5 * ( a->get(ij,kl) - a->get(ji,kl) );
                    }
               }
          }
     }
}//

void Tensor2d::antisymm_col4(const SharedTensor2d &a)
{
     #pragma omp parallel for
     for (int i=0; i < a->d1_; i++) {
          for (int j=0; j <= i; j++) {
               int ij = a->row_idx_[i][j];
               int ij2 = index2(i,j);
               for (int k=0; k < a->d3_; k++) {
                    for (int l=0; l <= k; l++) {
                         int kl = a->col_idx_[k][l];
                         int lk = a->col_idx_[l][k];
                         int kl2 = index2(k,l);
                         A2d_[ij2][kl2] = 0.5 * ( a->get(ij,kl) - a->get(ij,lk) );
                    }
               }
          }
     }
}//

void Tensor2d::symm_row_packed4(const SharedTensor2d &a)
{
     #pragma omp parallel for
     for (int i=0; i < a->d1_; i++) {
          for (int j=0; j <= i; j++) {
               int ij = a->row_idx_[i][j];
               int ji = a->row_idx_[j][i];
               int ij2 = index2(i,j);
               double perm = (i == j ? 1.0 : 2.0);
               for (int k=0; k < a->d3_; k++) {
                    for (int l=0; l <= k; l++) {
                         int kl = a->col_idx_[k][l];
                         int kl2 = index2(k,l);
                         A2d_[ij2][kl2] = 0.5 * perm * ( a->get(ij,kl) + a->get(ji,kl) );
                    }
               }
          }
     }
}//

void Tensor2d::symm_col_packed4(const SharedTensor2d &a)
{
     #pragma omp parallel for
     for (int i=0; i < a->d1_; i++) {
          for (int j=0; j <= i; j++) {
               int ij = a->row_idx_[i][j];
               int ji = a->row_idx_[j][i];
               int ij2 = index2(i,j);
               for (int k=0; k < a->d3_; k++) {
                    for (int l=0; l <= k; l++) {
                         int kl = a->col_idx_[k][l];
                         int kl2 = index2(k,l);
                         double perm = (k == l ? 1.0 : 2.0);
                         A2d_[ij2][kl2] = 0.5 * perm * ( a->get(ij,kl) + a->get(ji,kl) );
                    }
               }
          }
     }
}//

void Tensor2d::antisymm_row_packed4(const SharedTensor2d &a)
{
     #pragma omp parallel for
     for (int i=0; i < a->d1_; i++) {
          for (int j=0; j <= i; j++) {
               int ij = a->row_idx_[i][j];
               int ji = a->row_idx_[j][i];
               int ij2 = index2(i,j);
               double perm = (i == j ? 1.0 : 2.0);
               for (int k=0; k < a->d3_; k++) {
                    for (int l=0; l <= k; l++) {
                         int kl = a->col_idx_[k][l];
                         int kl2 = index2(k,l);
                         A2d_[ij2][kl2] = 0.5 * perm * ( a->get(ij,kl) - a->get(ji,kl) );
                    }
               }
          }
     }
}//

void Tensor2d::antisymm_col_packed4(const SharedTensor2d &a)
{
     #pragma omp parallel for
     for (int i=0; i < a->d1_; i++) {
          for (int j=0; j <= i; j++) {
               int ij = a->row_idx_[i][j];
               int ji = a->row_idx_[j][i];
               int ij2 = index2(i,j);
               for (int k=0; k < a->d3_; k++) {
                    for (int l=0; l <= k; l++) {
                         int kl = a->col_idx_[k][l];
                         int kl2 = index2(k,l);
                         double perm = (k == l ? 1.0 : 2.0);
                         A2d_[ij2][kl2] = 0.5 * perm * ( a->get(ij,kl) - a->get(ji,kl) );
                    }
               }
          }
     }
}//

void Tensor2d::tei_cs1_anti_symm(const SharedTensor2d &J, const SharedTensor2d &K)
{
    sort(1243, K, -1.0, 0.0);
    axpy(J, 2.0);
}//

void Tensor2d::tei_cs2_anti_symm(const SharedTensor2d &J, const SharedTensor2d &K)
{
    sort(2134, K, -1.0, 0.0);
    axpy(J, 2.0);
}//

void Tensor2d::tei_cs3_anti_symm(const SharedTensor2d &J, const SharedTensor2d &K)
{
    sort(1432, K, -1.0, 0.0);
    axpy(J, 2.0);
}//

void Tensor2d::tei_cs4_anti_symm(const SharedTensor2d &J, const SharedTensor2d &K)
{
    sort(3214, K, -1.0, 0.0);
    axpy(J, 2.0);
}//

void Tensor2d::P_ijab(const SharedTensor2d &A)
{
    // iajb --> ijab
    sort(1324, A, 1.0, 1.0);

    // jaib --> ijab
    sort(3124, A, -1.0, 1.0);

    // ibja --> ijab
    sort(1342, A, -1.0, 1.0);

    // jbia --> ijab
    sort(3142, A, 1.0, 1.0);

}//

void Tensor2d::cont444(int t_a1, int t_a2, int f_a1, int f_a2, const SharedTensor2d& A,
		     int t_b1, int t_b2, int f_b1, int f_b2, const SharedTensor2d& B,
		     double alpha, double beta)
{
    char ta, tb;
    int nca, ncb, ncc;
    int m, n, k;
    int r1, r2, c1, c2;
    int rr1, rr2, cc1, cc2;
    int dim_t, dim_u;

    // C(pq,rs) = \sum_{tu} A(pq,tu) B(tu,rs)
    if (t_a1 == 3 && t_a2 == 4 && t_b1 == 1 && t_b2 == 2) {
        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = A->dim2_;
        nca = k; // lda
        ncb = n; // ldb
        ncc = n; // ldc

        C_DGEMM(ta, tb, m, n, k, alpha, &(A->A2d_[0][0]), nca, &(B->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
    }

    else {

	// r1
	if (t_a1 == 1) {
	    r1 = t_a1;
	    dim_t = A->d1_;
	}
	else if (t_a2 == 1) {
	    r1 = t_a2;
	    dim_u = A->d1_;
	}
	else if (f_a1 == 1) r1 = f_a1;
	else if (f_a2 == 1) r1 = f_a2;

	// r2
	if (t_a1 == 2) {
	    r2 = t_a1;
	    dim_t = A->d2_;
	}
	else if (t_a2 == 2) {
	    r2 = t_a2;
	    dim_u = A->d2_;
	}
	else if (f_a1 == 2) r2 = f_a1;
	else if (f_a2 == 2) r2 = f_a2;

	// c1
	if (t_a1 == 3) {
	    c1 = t_a1;
	    dim_t = A->d3_;
	}
	else if (t_a2 == 3) {
	    c1 = t_a2;
	    dim_u = A->d3_;
	}
	else if (f_a1 == 3) c1 = f_a1;
	else if (f_a2 == 3) c1 = f_a2;

	// c2
	if (t_a1 == 4) {
	    c2 = t_a1;
	    dim_t = A->d4_;
	}
	else if (t_a2 == 4) {
	    c2 = t_a2;
	    dim_u = A->d4_;
	}
	else if (f_a1 == 4) c2 = f_a1;
	else if (f_a2 == 4) c2 = f_a2;

	//outfile->Printf("\tDimensions of A: %2d, %2d, %2d, %2d  \n", r1,r2,c1,c2);

	// Sort A(..,..) to A(pq,tu)
        SharedTensor2d temp1 = SharedTensor2d(new Tensor2d("temp1", d1_, d2_, dim_t, dim_u));
	#pragma omp parallel for
        for (int p = 0; p < d1_; p++) {
             for (int q = 0; q < d2_; q++) {
                  int pq = temp1->row_idx_[p][q];
                  for (int t = 0; t < dim_t; t++) {
                       for (int u = 0; u < dim_u; u++) {
                            int tu = temp1->col_idx_[t][u];

			    if (r1 == f_a1) rr1 = p;
			    else if (r1 == f_a2) rr1 = q;
			    else if (r1 == t_a1) rr1 = t;
			    else if (r1 == t_a2) rr1 = u;

			    if (r2 == f_a1) rr2 = p;
			    else if (r2 == f_a2) rr2 = q;
			    else if (r2 == t_a1) rr2 = t;
			    else if (r2 == t_a2) rr2 = u;

			    if (c1 == f_a1) cc1 = p;
			    else if (c1 == f_a2) cc1 = q;
			    else if (c1 == t_a1) cc1 = t;
			    else if (c1 == t_a2) cc1 = u;

			    if (c2 == f_a1) cc2 = p;
			    else if (c2 == f_a2) cc2 = q;
			    else if (c2 == t_a1) cc2 = t;
			    else if (c2 == t_a2) cc2 = u;

			    int row = rr2 + (rr1 * A->d2_);
			    int col = cc2 + (cc1 * A->d4_);

                            temp1->A2d_[pq][tu] = A->A2d_[row][col];
                       }
                  }
             }
        }
	//temp1->print();


        // r1
        if (t_b1 == 1) r1 = t_b1;
        else if (t_b2 == 1) r1 = t_b2;
        else if (f_b1 == 1) r1 = f_b1;
        else if (f_b2 == 1) r1 = f_b2;

        // r2
	if (t_b1 == 2) r2 = t_b1;
	else if (t_b2 == 2) r2 = t_b2;
	else if (f_b1 == 2) r2 = f_b1;
	else if (f_b2 == 2) r2 = f_b2;

	// c1
	if (t_b1 == 3) c1 = t_b1;
	else if (t_b2 == 3) c1 = t_b2;
	else if (f_b1 == 3) c1 = f_b1;
	else if (f_b2 == 3) c1 = f_b2;

        // c2
        if (t_b1 == 4) c2 = t_b1;
        else if (t_b2 == 4) c2 = t_b2;
        else if (f_b1 == 4) c2 = f_b1;
        else if (f_b2 == 4) c2 = f_b2;

	//outfile->Printf("\tDimensions of B: %2d, %2d, %2d, %2d  \n", r1,r2,c1,c2);

	// Sort B(..,..) to B(tu,rs)
        SharedTensor2d temp2 = SharedTensor2d(new Tensor2d("temp2", dim_t, dim_u, d3_, d4_));
	#pragma omp parallel for
        for (int t = 0; t < dim_t; t++) {
             for (int u = 0; u < dim_u; u++) {
                  int tu = temp2->row_idx_[t][u];
                  for (int r = 0; r < d3_; r++) {
                       for (int s = 0; s < d4_; s++) {
                            int rs = temp2->col_idx_[r][s];

			    if (r1 == f_b1) rr1 = r;
			    else if (r1 == f_b2) rr1 = s;
			    else if (r1 == t_b1) rr1 = t;
			    else if (r1 == t_b2) rr1 = u;

			    if (r2 == f_b1) rr2 = r;
			    else if (r2 == f_b2) rr2 = s;
			    else if (r2 == t_b1) rr2 = t;
			    else if (r2 == t_b2) rr2 = u;

			    if (c1 == f_b1) cc1 = r;
			    else if (c1 == f_b2) cc1 = s;
			    else if (c1 == t_b1) cc1 = t;
			    else if (c1 == t_b2) cc1 = u;

			    if (c2 == f_b1) cc2 = r;
			    else if (c2 == f_b2) cc2 = s;
			    else if (c2 == t_b1) cc2 = t;
			    else if (c2 == t_b2) cc2 = u;

			    int row = rr2 + (rr1 * B->d2_);
			    int col = cc2 + (cc1 * B->d4_);

                            temp2->A2d_[tu][rs] = B->A2d_[row][col];
                       }
                  }
             }
        }
	//temp2->print();

        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = temp1->dim2();
        nca = k; // lda
        ncb = n; // ldb
        ncc = n; // ldc

        C_DGEMM(ta, tb, m, n, k, alpha, &(temp1->A2d_[0][0]), nca, &(temp2->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
	temp1.reset();
	temp2.reset();

    }// end else



}//

void Tensor2d::cont444(bool delete_a, int t_a1, int t_a2, int f_a1, int f_a2, SharedTensor2d& A,
		     bool delete_b, int t_b1, int t_b2, int f_b1, int f_b2, SharedTensor2d& B,
		     double alpha, double beta)
{

    char ta, tb;
    int nca, ncb, ncc;
    int m, n, k;
    int r1, r2, c1, c2;
    int rr1, rr2, cc1, cc2;
    int dim_t, dim_u;

    // C(pq,rs) = \sum_{tu} A(pq,tu) B(tu,rs)
    if (t_a1 == 3 && t_a2 == 4 && t_b1 == 1 && t_b2 == 2) {
        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = A->dim2_;
        nca = k; // lda
        ncb = n; // ldb
        ncc = n; // ldc

        C_DGEMM(ta, tb, m, n, k, alpha, &(A->A2d_[0][0]), nca, &(B->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
    }

    else {

	// r1
	if (t_a1 == 1) {
	    r1 = t_a1;
	    dim_t = A->d1_;
	}
	else if (t_a2 == 1) {
	    r1 = t_a2;
	    dim_u = A->d1_;
	}
	else if (f_a1 == 1) r1 = f_a1;
	else if (f_a2 == 1) r1 = f_a2;

	// r2
	if (t_a1 == 2) {
	    r2 = t_a1;
	    dim_t = A->d2_;
	}
	else if (t_a2 == 2) {
	    r2 = t_a2;
	    dim_u = A->d2_;
	}
	else if (f_a1 == 2) r2 = f_a1;
	else if (f_a2 == 2) r2 = f_a2;

	// c1
	if (t_a1 == 3) {
	    c1 = t_a1;
	    dim_t = A->d3_;
	}
	else if (t_a2 == 3) {
	    c1 = t_a2;
	    dim_u = A->d3_;
	}
	else if (f_a1 == 3) c1 = f_a1;
	else if (f_a2 == 3) c1 = f_a2;

	// c2
	if (t_a1 == 4) {
	    c2 = t_a1;
	    dim_t = A->d4_;
	}
	else if (t_a2 == 4) {
	    c2 = t_a2;
	    dim_u = A->d4_;
	}
	else if (f_a1 == 4) c2 = f_a1;
	else if (f_a2 == 4) c2 = f_a2;

	//outfile->Printf("\tDimensions of A: %2d, %2d, %2d, %2d  \n", r1,r2,c1,c2);

	// Sort A(..,..) to A(pq,tu)
        SharedTensor2d temp1 = SharedTensor2d(new Tensor2d("temp1", d1_, d2_, dim_t, dim_u));
	#pragma omp parallel for
        for (int p = 0; p < d1_; p++) {
             for (int q = 0; q < d2_; q++) {
                  int pq = temp1->row_idx_[p][q];
                  for (int t = 0; t < dim_t; t++) {
                       for (int u = 0; u < dim_u; u++) {
                            int tu = temp1->col_idx_[t][u];

			    if (r1 == f_a1) rr1 = p;
			    else if (r1 == f_a2) rr1 = q;
			    else if (r1 == t_a1) rr1 = t;
			    else if (r1 == t_a2) rr1 = u;

			    if (r2 == f_a1) rr2 = p;
			    else if (r2 == f_a2) rr2 = q;
			    else if (r2 == t_a1) rr2 = t;
			    else if (r2 == t_a2) rr2 = u;

			    if (c1 == f_a1) cc1 = p;
			    else if (c1 == f_a2) cc1 = q;
			    else if (c1 == t_a1) cc1 = t;
			    else if (c1 == t_a2) cc1 = u;

			    if (c2 == f_a1) cc2 = p;
			    else if (c2 == f_a2) cc2 = q;
			    else if (c2 == t_a1) cc2 = t;
			    else if (c2 == t_a2) cc2 = u;

			    int row = rr2 + (rr1 * A->d2_);
			    int col = cc2 + (cc1 * A->d4_);

                            temp1->A2d_[pq][tu] = A->A2d_[row][col];
                       }
                  }
             }
        }
	//temp1->print();
	if (delete_a) A.reset();


        // r1
        if (t_b1 == 1) r1 = t_b1;
        else if (t_b2 == 1) r1 = t_b2;
        else if (f_b1 == 1) r1 = f_b1;
        else if (f_b2 == 1) r1 = f_b2;

        // r2
	if (t_b1 == 2) r2 = t_b1;
	else if (t_b2 == 2) r2 = t_b2;
	else if (f_b1 == 2) r2 = f_b1;
	else if (f_b2 == 2) r2 = f_b2;

	// c1
	if (t_b1 == 3) c1 = t_b1;
	else if (t_b2 == 3) c1 = t_b2;
	else if (f_b1 == 3) c1 = f_b1;
	else if (f_b2 == 3) c1 = f_b2;

        // c2
        if (t_b1 == 4) c2 = t_b1;
        else if (t_b2 == 4) c2 = t_b2;
        else if (f_b1 == 4) c2 = f_b1;
        else if (f_b2 == 4) c2 = f_b2;

	//outfile->Printf("\tDimensions of B: %2d, %2d, %2d, %2d  \n", r1,r2,c1,c2);

	// Sort B(..,..) to B(tu,rs)
        SharedTensor2d temp2 = SharedTensor2d(new Tensor2d("temp2", dim_t, dim_u, d3_, d4_));
	#pragma omp parallel for
        for (int t = 0; t < dim_t; t++) {
             for (int u = 0; u < dim_u; u++) {
                  int tu = temp2->row_idx_[t][u];
                  for (int r = 0; r < d3_; r++) {
                       for (int s = 0; s < d4_; s++) {
                            int rs = temp2->col_idx_[r][s];

			    if (r1 == f_b1) rr1 = r;
			    else if (r1 == f_b2) rr1 = s;
			    else if (r1 == t_b1) rr1 = t;
			    else if (r1 == t_b2) rr1 = u;

			    if (r2 == f_b1) rr2 = r;
			    else if (r2 == f_b2) rr2 = s;
			    else if (r2 == t_b1) rr2 = t;
			    else if (r2 == t_b2) rr2 = u;

			    if (c1 == f_b1) cc1 = r;
			    else if (c1 == f_b2) cc1 = s;
			    else if (c1 == t_b1) cc1 = t;
			    else if (c1 == t_b2) cc1 = u;

			    if (c2 == f_b1) cc2 = r;
			    else if (c2 == f_b2) cc2 = s;
			    else if (c2 == t_b1) cc2 = t;
			    else if (c2 == t_b2) cc2 = u;

			    int row = rr2 + (rr1 * B->d2_);
			    int col = cc2 + (cc1 * B->d4_);

                            temp2->A2d_[tu][rs] = B->A2d_[row][col];
                       }
                  }
             }
        }
	//temp2->print();
	if (delete_b) B.reset();

        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = temp1->dim2();
        nca = k; // lda
        ncb = n; // ldb
        ncc = n; // ldc

        C_DGEMM(ta, tb, m, n, k, alpha, &(temp1->A2d_[0][0]), nca, &(temp2->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
	temp1.reset();
	temp2.reset();

    }// end else

}//

void Tensor2d::cont444(string idx_c, string idx_a, string idx_b, bool delete_a, bool delete_b, SharedTensor2d& A, SharedTensor2d& B, double alpha, double beta)
{

    char ta, tb;
    int nca, ncb, ncc;
    int m, n, k;
    int r1, r2, c1, c2;
    int rr1, rr2, cc1, cc2;
    int dim_t, dim_u;
    int t_a1, t_a2, f_a1, f_a2;
    int t_b1, t_b2, f_b1, f_b2;

        // Expected order: C(pq,rs) = \sum_{tu} A(pq,tu) B(tu,rs)
	/*
        for (int i = 0; i < idx_a.size(); i++) {
	     outfile->Printf("\tIDX_A[%1d]: %c \n", i, idx_a[i]);
	}
	*/

	// Find dummy & free indices for A
	// f_a1
	if (idx_a[0] == idx_c[0]) f_a1 = 1;
	else if (idx_a[1] == idx_c[0]) f_a1 = 2;
	else if (idx_a[2] == idx_c[0]) f_a1 = 3;
	else if (idx_a[3] == idx_c[0]) f_a1 = 4;

	// f_a2
	if (idx_a[0] == idx_c[1]) f_a2 = 1;
	else if (idx_a[1] == idx_c[1]) f_a2 = 2;
	else if (idx_a[2] == idx_c[1]) f_a2 = 3;
	else if (idx_a[3] == idx_c[1]) f_a2 = 4;

	// t_a1
	if (f_a1 == 1 && f_a2 == 2) {
	    t_a1 = 3;
	    t_a2 = 4;
	}
	else if (f_a1 == 1 && f_a2 == 3) {
	    t_a1 = 2;
	    t_a2 = 4;
	}
	else if (f_a1 == 1 && f_a2 == 4) {
	    t_a1 = 2;
	    t_a2 = 3;
	}
	else if (f_a1 == 2 && f_a2 == 3) {
	    t_a1 = 1;
	    t_a2 = 4;
	}
	else if (f_a1 == 2 && f_a2 == 4) {
	    t_a1 = 1;
	    t_a2 = 3;
	}
	else if (f_a1 == 3 && f_a2 == 4) {
	    t_a1 = 1;
	    t_a2 = 2;
	}
	else if (f_a1 == 2 && f_a2 == 1) {
	    t_a1 = 3;
	    t_a2 = 4;
	}
	else if (f_a1 == 3 && f_a2 == 1) {
	    t_a1 = 2;
	    t_a2 = 4;
	}
	else if (f_a1 == 4 && f_a2 == 1) {
	    t_a1 = 2;
	    t_a2 = 3;
	}
	else if (f_a1 == 3 && f_a2 == 2) {
	    t_a1 = 1;
	    t_a2 = 4;
	}
	else if (f_a1 == 4 && f_a2 == 2) {
	    t_a1 = 1;
	    t_a2 = 3;
	}
	else if (f_a1 == 4 && f_a2 == 3) {
	    t_a1 = 1;
	    t_a2 = 2;
	}
	//outfile->Printf("\tf_a1, f_a2, t_a1, t_a2: %1d, %1d, %1d, %1d  \n", f_a1,f_a2,t_a1,t_a2);


	// Find dummy & free indices for B
	// f_b1
	if (idx_b[0] == idx_c[2]) f_b1 = 1;
	else if (idx_b[1] == idx_c[2]) f_b1 = 2;
	else if (idx_b[2] == idx_c[2]) f_b1 = 3;
	else if (idx_b[3] == idx_c[2]) f_b1 = 4;

	// f_b2
	if (idx_b[0] == idx_c[3]) f_b2 = 1;
	else if (idx_b[1] == idx_c[3]) f_b2 = 2;
	else if (idx_b[2] == idx_c[3]) f_b2 = 3;
	else if (idx_b[3] == idx_c[3]) f_b2 = 4;

	// t_b1
	if (idx_b[0] == idx_a[t_a1-1]) t_b1 = 1;
	else if (idx_b[1] == idx_a[t_a1-1]) t_b1 = 2;
	else if (idx_b[2] == idx_a[t_a1-1]) t_b1 = 3;
	else if (idx_b[3] == idx_a[t_a1-1]) t_b1 = 4;

	// t_b2
	if (idx_b[0] == idx_a[t_a2-1]) t_b2 = 1;
	else if (idx_b[1] == idx_a[t_a2-1]) t_b2 = 2;
	else if (idx_b[2] == idx_a[t_a2-1]) t_b2 = 3;
	else if (idx_b[3] == idx_a[t_a2-1]) t_b2 = 4;
	//outfile->Printf("\tf_b1, f_b2, t_b1, t_b2: %1d, %1d, %1d, %1d  \n", f_b1,f_b2,t_b1,t_b2);


        // Figure out A
	// r1
	if (t_a1 == 1) {
	    r1 = t_a1;
	    dim_t = A->d1_;
	}
	else if (t_a2 == 1) {
	    r1 = t_a2;
	    dim_u = A->d1_;
	}
	else if (f_a1 == 1) r1 = f_a1;
	else if (f_a2 == 1) r1 = f_a2;

	// r2
	if (t_a1 == 2) {
	    r2 = t_a1;
	    dim_t = A->d2_;
	}
	else if (t_a2 == 2) {
	    r2 = t_a2;
	    dim_u = A->d2_;
	}
	else if (f_a1 == 2) r2 = f_a1;
	else if (f_a2 == 2) r2 = f_a2;

	// c1
	if (t_a1 == 3) {
	    c1 = t_a1;
	    dim_t = A->d3_;
	}
	else if (t_a2 == 3) {
	    c1 = t_a2;
	    dim_u = A->d3_;
	}
	else if (f_a1 == 3) c1 = f_a1;
	else if (f_a2 == 3) c1 = f_a2;

	// c2
	if (t_a1 == 4) {
	    c2 = t_a1;
	    dim_t = A->d4_;
	}
	else if (t_a2 == 4) {
	    c2 = t_a2;
	    dim_u = A->d4_;
	}
	else if (f_a1 == 4) c2 = f_a1;
	else if (f_a2 == 4) c2 = f_a2;

	//outfile->Printf("\tDimensions of A: %2d, %2d, %2d, %2d  \n", r1,r2,c1,c2);

	// Sort A(..,..) to A(pq,tu)
        SharedTensor2d temp1 = SharedTensor2d(new Tensor2d("temp1", d1_, d2_, dim_t, dim_u));
	#pragma omp parallel for
        for (int p = 0; p < d1_; p++) {
             for (int q = 0; q < d2_; q++) {
                  int pq = temp1->row_idx_[p][q];
                  for (int t = 0; t < dim_t; t++) {
                       for (int u = 0; u < dim_u; u++) {
                            int tu = temp1->col_idx_[t][u];

			    if (r1 == f_a1) rr1 = p;
			    else if (r1 == f_a2) rr1 = q;
			    else if (r1 == t_a1) rr1 = t;
			    else if (r1 == t_a2) rr1 = u;

			    if (r2 == f_a1) rr2 = p;
			    else if (r2 == f_a2) rr2 = q;
			    else if (r2 == t_a1) rr2 = t;
			    else if (r2 == t_a2) rr2 = u;

			    if (c1 == f_a1) cc1 = p;
			    else if (c1 == f_a2) cc1 = q;
			    else if (c1 == t_a1) cc1 = t;
			    else if (c1 == t_a2) cc1 = u;

			    if (c2 == f_a1) cc2 = p;
			    else if (c2 == f_a2) cc2 = q;
			    else if (c2 == t_a1) cc2 = t;
			    else if (c2 == t_a2) cc2 = u;

			    int row = rr2 + (rr1 * A->d2_);
			    int col = cc2 + (cc1 * A->d4_);

                            temp1->A2d_[pq][tu] = A->A2d_[row][col];
                       }
                  }
             }
        }
	//temp1->print();
	if (delete_a) A.reset();

        // Figure out B
        // r1
        if (t_b1 == 1) r1 = t_b1;
        else if (t_b2 == 1) r1 = t_b2;
        else if (f_b1 == 1) r1 = f_b1;
        else if (f_b2 == 1) r1 = f_b2;

        // r2
	if (t_b1 == 2) r2 = t_b1;
	else if (t_b2 == 2) r2 = t_b2;
	else if (f_b1 == 2) r2 = f_b1;
	else if (f_b2 == 2) r2 = f_b2;

	// c1
	if (t_b1 == 3) c1 = t_b1;
	else if (t_b2 == 3) c1 = t_b2;
	else if (f_b1 == 3) c1 = f_b1;
	else if (f_b2 == 3) c1 = f_b2;

        // c2
        if (t_b1 == 4) c2 = t_b1;
        else if (t_b2 == 4) c2 = t_b2;
        else if (f_b1 == 4) c2 = f_b1;
        else if (f_b2 == 4) c2 = f_b2;

	//outfile->Printf("\tDimensions of B: %2d, %2d, %2d, %2d  \n", r1,r2,c1,c2);

	// Sort B(..,..) to B(tu,rs)
        SharedTensor2d temp2 = SharedTensor2d(new Tensor2d("temp2", dim_t, dim_u, d3_, d4_));
	#pragma omp parallel for
        for (int t = 0; t < dim_t; t++) {
             for (int u = 0; u < dim_u; u++) {
                  int tu = temp2->row_idx_[t][u];
                  for (int r = 0; r < d3_; r++) {
                       for (int s = 0; s < d4_; s++) {
                            int rs = temp2->col_idx_[r][s];

			    if (r1 == f_b1) rr1 = r;
			    else if (r1 == f_b2) rr1 = s;
			    else if (r1 == t_b1) rr1 = t;
			    else if (r1 == t_b2) rr1 = u;

			    if (r2 == f_b1) rr2 = r;
			    else if (r2 == f_b2) rr2 = s;
			    else if (r2 == t_b1) rr2 = t;
			    else if (r2 == t_b2) rr2 = u;

			    if (c1 == f_b1) cc1 = r;
			    else if (c1 == f_b2) cc1 = s;
			    else if (c1 == t_b1) cc1 = t;
			    else if (c1 == t_b2) cc1 = u;

			    if (c2 == f_b1) cc2 = r;
			    else if (c2 == f_b2) cc2 = s;
			    else if (c2 == t_b1) cc2 = t;
			    else if (c2 == t_b2) cc2 = u;

			    int row = rr2 + (rr1 * B->d2_);
			    int col = cc2 + (cc1 * B->d4_);

                            temp2->A2d_[tu][rs] = B->A2d_[row][col];
                       }
                  }
             }
        }
	//temp2->print();
	if (delete_b) B.reset();

        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = temp1->dim2();
        nca = k; // lda
        ncb = n; // ldb
        ncc = n; // ldc

        C_DGEMM(ta, tb, m, n, k, alpha, &(temp1->A2d_[0][0]), nca, &(temp2->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
	temp1.reset();
	temp2.reset();


}//

void Tensor2d::cont343(string idx_c, string idx_a, string idx_b, bool delete_b, SharedTensor2d& A, SharedTensor2d& B, double alpha, double beta)
{

    char ta, tb;
    int nca, ncb, ncc;
    int m, n, k;
    int r1, r2, c1, c2;
    int rr1, rr2, cc1, cc2;
    int t_b1, t_b2, f_b1, f_b2;

        // Expected order: C(Q,pq) = \sum_{rs} A(Q,rs) B(rs,pq)

	// Find dummy & free indices for B
	// f_b1
	if (idx_b[0] == idx_c[0]) f_b1 = 1;
	else if (idx_b[1] == idx_c[0]) f_b1 = 2;
	else if (idx_b[2] == idx_c[0]) f_b1 = 3;
	else if (idx_b[3] == idx_c[0]) f_b1 = 4;

	// f_b2
	if (idx_b[0] == idx_c[1]) f_b2 = 1;
	else if (idx_b[1] == idx_c[1]) f_b2 = 2;
	else if (idx_b[2] == idx_c[1]) f_b2 = 3;
	else if (idx_b[3] == idx_c[1]) f_b2 = 4;

	// t_b1
	if (idx_b[0] == idx_a[0]) t_b1 = 1;
	else if (idx_b[1] == idx_a[0]) t_b1 = 2;
	else if (idx_b[2] == idx_a[0]) t_b1 = 3;
	else if (idx_b[3] == idx_a[0]) t_b1 = 4;

	// t_b2
	if (idx_b[0] == idx_a[1]) t_b2 = 1;
	else if (idx_b[1] == idx_a[1]) t_b2 = 2;
	else if (idx_b[2] == idx_a[1]) t_b2 = 3;
	else if (idx_b[3] == idx_a[1]) t_b2 = 4;


        // Figure out B
        // r1
        if (t_b1 == 1) r1 = t_b1;
        else if (t_b2 == 1) r1 = t_b2;
        else if (f_b1 == 1) r1 = f_b1;
        else if (f_b2 == 1) r1 = f_b2;

        // r2
	if (t_b1 == 2) r2 = t_b1;
	else if (t_b2 == 2) r2 = t_b2;
	else if (f_b1 == 2) r2 = f_b1;
	else if (f_b2 == 2) r2 = f_b2;

	// c1
	if (t_b1 == 3) c1 = t_b1;
	else if (t_b2 == 3) c1 = t_b2;
	else if (f_b1 == 3) c1 = f_b1;
	else if (f_b2 == 3) c1 = f_b2;

        // c2
        if (t_b1 == 4) c2 = t_b1;
        else if (t_b2 == 4) c2 = t_b2;
        else if (f_b1 == 4) c2 = f_b1;
        else if (f_b2 == 4) c2 = f_b2;

	//outfile->Printf("\tDimensions of B: %2d, %2d, %2d, %2d  \n", r1,r2,c1,c2);

	// Sort B(..,..) to B(rs,pq)
        SharedTensor2d temp = SharedTensor2d(new Tensor2d("temp", A->d2_, A->d3_, d2_, d3_));
	#pragma omp parallel for
        for (int r = 0; r < A->d2_; r++) {
             for (int s = 0; s < A->d3_; s++) {
                  int rs = temp->row_idx_[r][s];
                  for (int p = 0; p < d2_; p++) {
                       for (int q = 0; q < d3_; q++) {
                            int pq = temp->col_idx_[p][q];

			    if (r1 == f_b1) rr1 = p;
			    else if (r1 == f_b2) rr1 = q;
			    else if (r1 == t_b1) rr1 = r;
			    else if (r1 == t_b2) rr1 = s;

			    if (r2 == f_b1) rr2 = p;
			    else if (r2 == f_b2) rr2 = q;
			    else if (r2 == t_b1) rr2 = r;
			    else if (r2 == t_b2) rr2 = s;

			    if (c1 == f_b1) cc1 = p;
			    else if (c1 == f_b2) cc1 = q;
			    else if (c1 == t_b1) cc1 = r;
			    else if (c1 == t_b2) cc1 = s;

			    if (c2 == f_b1) cc2 = p;
			    else if (c2 == f_b2) cc2 = q;
			    else if (c2 == t_b1) cc2 = r;
			    else if (c2 == t_b2) cc2 = s;

			    int row = rr2 + (rr1 * B->d2_);
			    int col = cc2 + (cc1 * B->d4_);

                            temp->A2d_[rs][pq] = B->A2d_[row][col];
                       }
                  }
             }
        }
	//temp->print();
	if (delete_b) B.reset();

        ta = 'n';
        tb = 'n';
        m = d1_;
        n = d2_ * d3_;
        k = temp->dim1();
        nca = k; // lda
        ncb = n; // ldb
        ncc = n; // ldc

        C_DGEMM(ta, tb, m, n, k, alpha, &(A->A2d_[0][0]), nca, &(temp->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
	temp.reset();

}//

void Tensor2d::cont442(string idx_c, string idx_a, string idx_b, bool delete_a, bool delete_b, SharedTensor2d& A, SharedTensor2d& B, double alpha, double beta)
{

    char ta, tb;
    int nca, ncb, ncc;
    int m, n, k;
    int r1, r2, c1, c2;
    int rr1, rr2, cc1, cc2;
    int dim_r, dim_s, dim_t;
    int t_a1, t_a2, t_a3, f_a1;
    int t_b1, t_b2, t_b3, f_b1;

        // Expected order: C(pq) = \sum_{rst} A(pr,st) B(rs,tq)

	// Find dummy & free indices for A
	// f_a1
	if (idx_a[0] == idx_c[0]) f_a1 = 1;
	else if (idx_a[1] == idx_c[0]) f_a1 = 2;
	else if (idx_a[2] == idx_c[0]) f_a1 = 3;
	else if (idx_a[3] == idx_c[0]) f_a1 = 4;

	// t_a1, t_a2, t_a3
	if (f_a1 == 1) {
	    t_a1 = 2;
	    t_a2 = 3;
	    t_a3 = 4;
	}
	else if (f_a1 == 2) {
	    t_a1 = 1;
	    t_a2 = 3;
	    t_a3 = 4;
	}
	else if (f_a1 == 3) {
	    t_a1 = 1;
	    t_a2 = 2;
	    t_a3 = 4;
	}
	else if (f_a1 == 4) {
	    t_a1 = 1;
	    t_a2 = 2;
	    t_a3 = 3;
	}

	// Find dummy & free indices for B
	// f_b1
	if (idx_b[0] == idx_c[1]) f_b1 = 1;
	else if (idx_b[1] == idx_c[1]) f_b1 = 2;
	else if (idx_b[2] == idx_c[1]) f_b1 = 3;
	else if (idx_b[3] == idx_c[1]) f_b1 = 4;

	// t_b1
	if (idx_b[0] == idx_a[t_a1-1]) t_b1 = 1;
	else if (idx_b[1] == idx_a[t_a1-1]) t_b1 = 2;
	else if (idx_b[2] == idx_a[t_a1-1]) t_b1 = 3;
	else if (idx_b[3] == idx_a[t_a1-1]) t_b1 = 4;

	// t_b2
	if (idx_b[0] == idx_a[t_a2-1]) t_b2 = 1;
	else if (idx_b[1] == idx_a[t_a2-1]) t_b2 = 2;
	else if (idx_b[2] == idx_a[t_a2-1]) t_b2 = 3;
	else if (idx_b[3] == idx_a[t_a2-1]) t_b2 = 4;

	// t_b3
	if (idx_b[0] == idx_a[t_a3-1]) t_b3 = 1;
	else if (idx_b[1] == idx_a[t_a3-1]) t_b3 = 2;
	else if (idx_b[2] == idx_a[t_a3-1]) t_b3 = 3;
	else if (idx_b[3] == idx_a[t_a3-1]) t_b3 = 4;


        // Figure out A
	// r1
	if (t_a1 == 1) {
	    r1 = t_a1;
	    dim_r = A->d1_;
	}
	else if (t_a2 == 1) {
	    r1 = t_a2;
	    dim_s = A->d1_;
	}
	else if (t_a3 == 1) {
	    r1 = t_a3;
	    dim_t = A->d1_;
	}
	else if (f_a1 == 1) r1 = f_a1;

	// r2
	if (t_a1 == 2) {
	    r2 = t_a1;
	    dim_r = A->d2_;
	}
	else if (t_a2 == 2) {
	    r2 = t_a2;
	    dim_s = A->d2_;
	}
	else if (t_a3 == 2) {
	    r2 = t_a3;
	    dim_t = A->d2_;
	}
	else if (f_a1 == 2) r2 = f_a1;

	// c1
	if (t_a1 == 3) {
	    c1 = t_a1;
	    dim_r = A->d3_;
	}
	else if (t_a2 == 3) {
	    c1 = t_a2;
	    dim_s = A->d3_;
	}
	else if (t_a3 == 3) {
	    c1 = t_a3;
	    dim_t = A->d3_;
	}
	else if (f_a1 == 3) c1 = f_a1;

	// c2
	if (t_a1 == 4) {
	    c2 = t_a1;
	    dim_r = A->d4_;
	}
	else if (t_a2 == 4) {
	    c2 = t_a2;
	    dim_s = A->d4_;
	}
	else if (t_a3 == 4) {
	    c2 = t_a3;
	    dim_t = A->d4_;
	}
	else if (f_a1 == 3) c2 = f_a1;


	// Sort A(..,..) to A(pr,st)
        SharedTensor2d temp1 = SharedTensor2d(new Tensor2d("temp1", dim1_, dim_r, dim_s, dim_t));
	#pragma omp parallel for
        for (int p = 0; p < dim1_; p++) {
             for (int r = 0; r < dim_r; r++) {
                  int pr = temp1->row_idx_[p][r];
                  for (int s = 0; s < dim_s; s++) {
                       for (int t = 0; t < dim_t; t++) {
                            int st = temp1->col_idx_[s][t];

			    if (r1 == f_a1) rr1 = p;
			    else if (r1 == t_a1) rr1 = r;
			    else if (r1 == t_a2) rr1 = s;
			    else if (r1 == t_a3) rr1 = t;

			    if (r2 == f_a1) rr2 = p;
			    else if (r2 == t_a1) rr2 = r;
			    else if (r2 == t_a2) rr2 = s;
			    else if (r2 == t_a3) rr2 = t;

			    if (c1 == f_a1) cc1 = p;
			    else if (c1 == t_a1) cc1 = r;
			    else if (c1 == t_a2) cc1 = s;
			    else if (c1 == t_a3) cc1 = t;

			    if (c2 == f_a1) cc2 = p;
			    else if (c2 == t_a1) cc2 = r;
			    else if (c2 == t_a2) cc2 = s;
			    else if (c2 == t_a3) cc2 = t;

			    int row = rr2 + (rr1 * A->d2_);
			    int col = cc2 + (cc1 * A->d4_);

                            temp1->A2d_[pr][st] = A->A2d_[row][col];
                       }
                  }
             }
        }
	//temp1->print();
	if (delete_a) A.reset();

        // Figure out B
        // r1
        if (t_b1 == 1) r1 = t_b1;
        else if (t_b2 == 1) r1 = t_b2;
        else if (t_b3 == 1) r1 = t_b3;
        else if (f_b1 == 1) r1 = f_b1;

        // r2
	if (t_b1 == 2) r2 = t_b1;
        else if (t_b2 == 2) r2 = t_b2;
        else if (t_b3 == 2) r2 = t_b3;
        else if (f_b1 == 2) r2 = f_b1;

	// c1
	if (t_b1 == 3) c1 = t_b1;
        else if (t_b2 == 3) c1 = t_b2;
        else if (t_b3 == 3) c1 = t_b3;
        else if (f_b1 == 3) c1 = f_b1;

        // c2
        if (t_b1 == 4) c2 = t_b1;
        else if (t_b2 == 4) c2 = t_b2;
        else if (t_b3 == 4) c2 = t_b3;
        else if (f_b1 == 4) c2 = f_b1;

	// Sort B(..,..) to B(rs,tq)
        SharedTensor2d temp2 = SharedTensor2d(new Tensor2d("temp2", dim_r, dim_s, dim_t, dim2_));
	#pragma omp parallel for
        for (int r = 0; r < dim_r; r++) {
             for (int s = 0; s < dim_s; s++) {
                  int rs = temp2->row_idx_[r][s];
                  for (int t = 0; t < dim_t; t++) {
                       for (int q = 0; q < dim2_; q++) {
                            int tq = temp2->col_idx_[t][q];

			    if (r1 == t_b1) rr1 = r;
			    else if (r1 == t_b2) rr1 = s;
			    else if (r1 == t_b3) rr1 = t;
			    else if (r1 == f_b1) rr1 = q;

			    if (r2 == t_b1) rr2 = r;
			    else if (r2 == t_b2) rr2 = s;
			    else if (r2 == t_b3) rr2 = t;
			    else if (r2 == f_b1) rr2 = q;

			    if (c1 == t_b1) cc1 = r;
			    else if (c1 == t_b2) cc1 = s;
			    else if (c1 == t_b3) cc1 = t;
			    else if (c1 == f_b1) cc1 = q;

			    if (c2 == t_b1) cc2 = r;
			    else if (c2 == t_b2) cc2 = s;
			    else if (c2 == t_b3) cc2 = t;
			    else if (c2 == f_b1) cc2 = q;

			    int row = rr2 + (rr1 * B->d2_);
			    int col = cc2 + (cc1 * B->d4_);

                            temp2->A2d_[rs][tq] = B->A2d_[row][col];
                       }
                  }
             }
        }
	//temp2->print();
	if (delete_b) B.reset();

        ta = 'n';
        tb = 'n';
        m = dim1_;
        n = dim2_;
        k = dim_r * dim_s * dim_t;
        nca = k; // lda
        ncb = n; // ldb
        ncc = n; // ldc

        C_DGEMM(ta, tb, m, n, k, alpha, &(temp1->A2d_[0][0]), nca, &(temp2->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
	temp1.reset();
	temp2.reset();


}//

void Tensor2d::cont424(string idx_c, string idx_a, string idx_b, bool delete_a, SharedTensor2d& A, SharedTensor2d& B, double alpha, double beta)
{

    char ta, tb;
    int nca, ncb, ncc;
    int m, n, k;
    int r1, r2, c1, c2;
    int rr1, rr2, cc1, cc2;
    int dim_t;
    int t_a1, f_a1, f_a2, f_a3;
    int t_b1, f_b1;

        // Expected order: C(pq,rs) = \sum_{t} A(pq,rt) B(t,s)

	// Find dummy & free indices for A
	// f_a1
	if (idx_a[0] == idx_c[0]) f_a1 = 1;
	else if (idx_a[1] == idx_c[0]) f_a1 = 2;
	else if (idx_a[2] == idx_c[0]) f_a1 = 3;
	else if (idx_a[3] == idx_c[0]) f_a1 = 4;

	// f_a2
	if (idx_a[0] == idx_c[1]) f_a2 = 1;
	else if (idx_a[1] == idx_c[1]) f_a2 = 2;
	else if (idx_a[2] == idx_c[1]) f_a2 = 3;
	else if (idx_a[3] == idx_c[1]) f_a2 = 4;

	// f_a3
	if (idx_a[0] == idx_c[2]) f_a3 = 1;
	else if (idx_a[1] == idx_c[2]) f_a3 = 2;
	else if (idx_a[2] == idx_c[2]) f_a3 = 3;
	else if (idx_a[3] == idx_c[2]) f_a3 = 4;

	// t_a1
	if (idx_a[0] != idx_c[0] && idx_a[0] != idx_c[1] && idx_a[0] != idx_c[2] && idx_a[0] != idx_c[3]) t_a1 = 1;
	else if (idx_a[1] != idx_c[0] && idx_a[1] != idx_c[1] && idx_a[1] != idx_c[2] && idx_a[1] != idx_c[3]) t_a1 = 2;
	else if (idx_a[2] != idx_c[0] && idx_a[2] != idx_c[1] && idx_a[2] != idx_c[2] && idx_a[2] != idx_c[3]) t_a1 = 3;
	else if (idx_a[3] != idx_c[0] && idx_a[3] != idx_c[1] && idx_a[3] != idx_c[2] && idx_a[3] != idx_c[3]) t_a1 = 4;
	//outfile->Printf("\tf_a1, f_a2, f_a3, t_a1: %1d, %1d, %1d, %1d  \n", f_a1,f_a2,f_a3,t_a1);

	// Find dummy & free indices for B
	// f_b1 and t_b1
	if (idx_b[0] == idx_c[3]) {
	    f_b1 = 1;
	    t_b1 = 2;
	    dim_t = B->dim2();
	}
	else if (idx_b[1] == idx_c[3]) {
	    f_b1 = 2;
	    t_b1 = 1;
	    dim_t = B->dim1();
	}
	//outfile->Printf("\tf_b1, t_b1: %1d, %1d \n", f_b1,t_b1);

        // Figure out A
	// r1
	if (t_a1 == 1) r1 = t_a1;
	else if (f_a1 == 1) r1 = f_a1;
	else if (f_a2 == 1) r1 = f_a2;
	else if (f_a3 == 1) r1 = f_a3;

	// r2
	if (t_a1 == 2) r2 = t_a1;
	else if (f_a1 == 2) r2 = f_a1;
	else if (f_a2 == 2) r2 = f_a2;
	else if (f_a3 == 2) r2 = f_a3;

	// c1
	if (t_a1 == 3) c1 = t_a1;
	else if (f_a1 == 3) c1 = f_a1;
	else if (f_a2 == 3) c1 = f_a2;
	else if (f_a3 == 3) c1 = f_a3;

	// c2
	if (t_a1 == 4) c2 = t_a1;
	else if (f_a1 == 4) c2 = f_a1;
	else if (f_a2 == 4) c2 = f_a2;
	else if (f_a3 == 4) c2 = f_a3;

	// Sort A(..,..) to A(pq,rt)
        SharedTensor2d temp = SharedTensor2d(new Tensor2d("temp", d1_, d2_, d3_, dim_t));
	#pragma omp parallel for
        for (int p = 0; p < d1_; p++) {
             for (int q = 0; q < d2_; q++) {
                  int pq = temp->row_idx_[p][q];
                  for (int r = 0; r < d3_; r++) {
                       for (int t = 0; t < dim_t; t++) {
                            int rt = temp->col_idx_[r][t];

			    if (r1 == f_a1) rr1 = p;
			    else if (r1 == f_a2) rr1 = q;
			    else if (r1 == f_a3) rr1 = r;
			    else if (r1 == t_a1) rr1 = t;

			    if (r2 == f_a1) rr2 = p;
			    else if (r2 == f_a2) rr2 = q;
			    else if (r2 == f_a3) rr2 = r;
			    else if (r2 == t_a1) rr2 = t;

			    if (c1 == f_a1) cc1 = p;
			    else if (c1 == f_a2) cc1 = q;
			    else if (c1 == f_a3) cc1 = r;
			    else if (c1 == t_a1) cc1 = t;

			    if (c2 == f_a1) cc2 = p;
			    else if (c2 == f_a2) cc2 = q;
			    else if (c2 == f_a3) cc2 = r;
			    else if (c2 == t_a1) cc2 = t;

			    int row = rr2 + (rr1 * A->d2_);
			    int col = cc2 + (cc1 * A->d4_);

                            temp->A2d_[pq][rt] = A->A2d_[row][col];
                       }
                  }
             }
        }
	if (delete_a) A.reset();

        ta = 'n';
        if (t_b1 == 1) tb = 'n';
	else tb = 't';
        m = d1_ * d2_ * d3_;
        n = d4_;
        k = temp->d4_;
        nca = k; // lda
	// ldb
	if (tb == 't') ncb = k;
	else ncb = n;
        ncc = n; // ldc

        C_DGEMM(ta, tb, m, n, k, alpha, &(temp->A2d_[0][0]), nca, &(B->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
	temp.reset();
}//

void Tensor2d::cont244(string idx_c, string idx_a, string idx_b, bool delete_b, SharedTensor2d& A, SharedTensor2d& B, double alpha, double beta)
{

    char ta, tb;
    int nca, ncb, ncc;
    int m, n, k;
    int r1, r2, c1, c2;
    int rr1, rr2, cc1, cc2;
    int dim_t;
    int t_a1, f_a1;
    int t_b1, f_b1, f_b2, f_b3;

        // Expected order: C(pq,rs) = \sum_{t} A(p,t) B(tq,rs)

	// Find dummy & free indices for B
	// f_b1
	if (idx_b[0] == idx_c[1]) f_b1 = 1;
	else if (idx_b[1] == idx_c[1]) f_b1 = 2;
	else if (idx_b[2] == idx_c[1]) f_b1 = 3;
	else if (idx_b[3] == idx_c[1]) f_b1 = 4;

	// f_b2
	if (idx_b[0] == idx_c[2]) f_b2 = 1;
	else if (idx_b[1] == idx_c[2]) f_b2 = 2;
	else if (idx_b[2] == idx_c[2]) f_b2 = 3;
	else if (idx_b[3] == idx_c[2]) f_b2 = 4;

	// f_b3
	if (idx_b[0] == idx_c[3]) f_b3 = 1;
	else if (idx_b[1] == idx_c[3]) f_b3 = 2;
	else if (idx_b[2] == idx_c[3]) f_b3 = 3;
	else if (idx_b[3] == idx_c[3]) f_b3 = 4;

	// t_b1
	if (idx_b[0] != idx_c[0] && idx_b[0] != idx_c[1] && idx_b[0] != idx_c[2] && idx_b[0] != idx_c[3]) t_b1 = 1;
	else if (idx_b[1] != idx_c[0] && idx_b[1] != idx_c[1] && idx_b[1] != idx_c[2] && idx_b[1] != idx_c[3]) t_b1 = 2;
	else if (idx_b[2] != idx_c[0] && idx_b[2] != idx_c[1] && idx_b[2] != idx_c[2] && idx_b[2] != idx_c[3]) t_b1 = 3;
	else if (idx_b[3] != idx_c[0] && idx_b[3] != idx_c[1] && idx_b[3] != idx_c[2] && idx_b[3] != idx_c[3]) t_b1 = 4;

	// Find dummy & free indices for A
	// f_a1 and t_a1
	if (idx_a[0] == idx_c[0]) {
	    f_a1 = 1;
	    t_a1 = 2;
	    dim_t = A->dim2();
	}
	else if (idx_a[1] == idx_c[0]) {
	    f_a1 = 2;
	    t_a1 = 1;
	    dim_t = A->dim1();
	}

        // Figure out B
	// r1
	if (t_b1 == 1) r1 = t_b1;
	else if (f_b1 == 1) r1 = f_b1;
	else if (f_b2 == 1) r1 = f_b2;
	else if (f_b3 == 1) r1 = f_b3;

	// r2
	if (t_b1 == 2) r2 = t_b1;
	else if (f_b1 == 2) r2 = f_b1;
	else if (f_b2 == 2) r2 = f_b2;
	else if (f_b3 == 2) r2 = f_b3;

	// c1
	if (t_b1 == 3) c1 = t_b1;
	else if (f_b1 == 3) c1 = f_b1;
	else if (f_b2 == 3) c1 = f_b2;
	else if (f_b3 == 3) c1 = f_b3;

	// c2
	if (t_b1 == 4) c2 = t_b1;
	else if (f_b1 == 4) c2 = f_b1;
	else if (f_b2 == 4) c2 = f_b2;
	else if (f_b3 == 4) c2 = f_b3;

	// Sort B(..,..) to B(tq,rs)
        SharedTensor2d temp = SharedTensor2d(new Tensor2d("temp", dim_t, d2_, d3_, d4_));
	#pragma omp parallel for
        for (int t = 0; t < dim_t; t++) {
             for (int q = 0; q < d2_; q++) {
                  int tq = temp->row_idx_[t][q];
                  for (int r = 0; r < d3_; r++) {
                       for (int s = 0; s < d4_; s++) {
                            int rs = temp->col_idx_[r][s];

			    if (r1 == t_b1) rr1 = t;
			    else if (r1 == f_b1) rr1 = q;
			    else if (r1 == f_b2) rr1 = r;
			    else if (r1 == f_b3) rr1 = s;

			    if (r2 == t_b1) rr2 = t;
			    else if (r2 == f_b1) rr2 = q;
			    else if (r2 == f_b2) rr2 = r;
			    else if (r2 == f_b3) rr2 = s;

			    if (c1 == t_b1) cc1 = t;
			    else if (c1 == f_b1) cc1 = q;
			    else if (c1 == f_b2) cc1 = r;
			    else if (c1 == f_b3) cc1 = s;

			    if (c2 == t_b1) cc2 = t;
			    else if (c2 == f_b1) cc2 = q;
			    else if (c2 == f_b2) cc2 = r;
			    else if (c2 == f_b3) cc2 = s;

			    int row = rr2 + (rr1 * B->d2_);
			    int col = cc2 + (cc1 * B->d4_);

                            temp->A2d_[tq][rs] = B->A2d_[row][col];
                       }
                  }
             }
        }
	if (delete_b) B.reset();

        if (t_a1 == 2) ta = 'n';
	else ta = 't';
	tb = 'n';
        m = d1_;
        n = d2_ * d3_ * d4_;
        k = temp->d1_;
        //nca = transa ? m : k; // lda
        //ncb = transb ? k : n; // ldb
	if (ta == 't') nca = m;
	else nca = k;
	ncb = n;
        ncc = n; // ldc

        C_DGEMM(ta, tb, m, n, k, alpha, &(A->A2d_[0][0]), nca, &(temp->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
	temp.reset();

}//

void Tensor2d::cont233(string idx_c, string idx_a, string idx_b, SharedTensor2d& A, SharedTensor2d& B, double alpha, double beta)
{

    char ta, tb;
    int lda, ldb, ldc;
    int m, n, k;
    int dim_r;
    int t_a1, f_a1;
    int t_b1, f_b1;

        // Expected order: C(Q,pq) = \sum_{r} A(p,r) B(Q,rq)

	// Find dummy & free indices for A
	// f_a1 and t_a1
	if (idx_a[0] == idx_c[0]) {
	    f_a1 = 1;
	    t_a1 = 2;
	    dim_r = A->dim2();
	}
	else if (idx_a[1] == idx_c[0]) {
	    f_a1 = 2;
	    t_a1 = 1;
	    dim_r = A->dim1();
	}

	// Find dummy & free indices for B
	// f_b1 and t_b1
	if (idx_b[0] == idx_c[1]) {
	    f_b1 = 1;
	    t_b1 = 2;
	}
	else if (idx_b[1] == idx_c[1]) {
	    f_b1 = 2;
	    t_b1 = 1;
	}


	m = d2_;
	n = d3_;

        if (t_a1 == 2) ta = 'n';
	else ta = 't';

        if (t_b1 == 1) tb = 'n';
	else tb = 't';

	if (ta == 'n') k = A->dim2();
	else k = A->dim1();

        //lda = transa ? m : k;
	if (ta == 't') lda = m;
	else lda = k;

        //ldb = transb ? k : n;
	if (tb == 't') ldb = k;
	else ldb = n;

        ldc = n;

        #pragma omp parallel for
        for (int Q = 0; Q < dim1_; Q++) {
             C_DGEMM(ta, tb, m, n, k, alpha, A->A2d_[0], lda, B->A2d_[Q], ldb, beta, A2d_[Q], ldc);
        }

}//

void Tensor2d::cont323(string idx_c, string idx_a, string idx_b, bool delete_a, SharedTensor2d& A, SharedTensor2d& B, double alpha, double beta)
{

    char ta, tb;
    int nca, ncb, ncc;
    int m, n, k;
    int dim_r;
    int t_a1, f_a1;
    int t_b1, f_b1;
    int r1, c1;
    int rr1, cc1;

        // Expected order: C(Q,pq) = \sum_{r} A(Q,pr) B(r,q)

	// Find dummy & free indices for A
	// f_a1 and t_a1
	if (idx_a[0] == idx_c[0]) {
	    f_a1 = 1;
	    t_a1 = 2;
	}
	else if (idx_a[1] == idx_c[0]) {
	    f_a1 = 2;
	    t_a1 = 1;
	}

	// Find dummy & free indices for B
	// f_b1 and t_b1
	if (idx_b[0] == idx_c[1]) {
	    f_b1 = 1;
	    t_b1 = 2;
	    dim_r = B->dim2();
	}
	else if (idx_b[1] == idx_c[1]) {
	    f_b1 = 2;
	    t_b1 = 1;
	    dim_r = B->dim1();
	}

        // Figure out A
	// r1
	if (t_a1 == 1) r1 = t_a1;
	else if (f_a1 == 1) r1 = f_a1;

	// c1
	if (t_a1 == 2) c1 = t_a1;
	else if (f_a1 == 2) c1 = f_a1;

	// Sort A(Q,..) to A(Q,pr)
        SharedTensor2d temp = SharedTensor2d(new Tensor2d("temp", d1_, d2_, dim_r));
	#pragma omp parallel for
        for (int Q = 0; Q < dim1_; Q++) {
             for (int p = 0; p < d2_; p++) {
                  for (int r = 0; r < dim_r; r++) {
		       int pr = r + (p*dim_r);

		       if (r1 == f_a1) rr1 = p;
		       else if (r1 == t_a1) rr1 = r;

		       if (c1 == f_a1) cc1 = p;
	               else if (c1 == t_a1) cc1 = r;

		       int col = cc1 + (rr1 * A->d3_);

                       temp->A2d_[Q][pr] = A->A2d_[Q][col];
                  }
             }
        }
	if (delete_a) A.reset();

	m = d1_ * d2_;
	n = d3_;
	k = dim_r;
	ta = 'n';
        if (t_b1 == 1) tb = 'n';
	else tb = 't';
	nca = k; // lda
        //ldb = transb ? k : n;
	if (tb == 't') ncb = k;
	else ncb = n;
        ncc = n;

        C_DGEMM(ta, tb, m, n, k, alpha, &(temp->A2d_[0][0]), nca, &(B->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
	temp.reset();
}//

void Tensor2d::cont332(string idx_c, string idx_a, string idx_b, bool delete_a, bool delete_b, SharedTensor2d& A, SharedTensor2d& B, double alpha, double beta)
{

    char ta, tb;
    int nca, ncb, ncc; // number of columns
    int m, n, k;
    int dim_r;
    int t_a1, f_a1;
    int t_b1, f_b1;
    int ra, ca, rb, cb;
    int rra, cca, rrb, ccb;

        // Expected order: C(pq) = \sum_{Qr} A(Q,rp) B(Q,rq)

	// Find dummy & free indices for A
	// f_a1 and t_a1
	if (idx_a[0] == idx_c[0]) {
	    f_a1 = 1;
	    t_a1 = 2;
	    dim_r = A->d3_;
	}
	else if (idx_a[1] == idx_c[0]) {
	    f_a1 = 2;
	    t_a1 = 1;
	    dim_r = A->d2_;
	}

	// Find dummy & free indices for B
	// f_b1 and t_b1
	if (idx_b[0] == idx_c[1]) {
	    f_b1 = 1;
	    t_b1 = 2;
	}
	else if (idx_b[1] == idx_c[1]) {
	    f_b1 = 2;
	    t_b1 = 1;
	}

        // Figure out A
	// ra
	if (t_a1 == 1) ra = t_a1;
	else if (f_a1 == 1) ra = f_a1;

	// ca
	if (t_a1 == 2) ca = t_a1;
	else if (f_a1 == 2) ca = f_a1;

	// Sort A(Q,..) to A(Q,rp)
        SharedTensor2d temp1 = SharedTensor2d(new Tensor2d("temp1", A->d1_, dim_r, dim1_));
	#pragma omp parallel for
        for (int Q = 0; Q < A->d1_; Q++) {
             for (int r = 0; r < dim_r; r++) {
                  for (int p = 0; p < dim1_; p++) {
		       int rp = p + (r*dim1_);

		       if (ra == t_a1) rra = r;
		       else if (ra == f_a1) rra = p;

		       if (ca == t_a1) cca = r;
	               else if (ca == f_a1) cca = p;

		       int col = cca + (rra * A->d3_);

                       temp1->A2d_[Q][rp] = A->A2d_[Q][col];
                  }
             }
        }
	if (delete_a) A.reset();

        // Figure out B
	// rb
	if (t_b1 == 1) rb = t_b1;
	else if (f_b1 == 1) rb = f_b1;

	// cb
	if (t_b1 == 2) cb = t_b1;
	else if (f_b1 == 2) cb = f_b1;

	// Sort B(Q,..) to B(Q,rq)
        SharedTensor2d temp2 = SharedTensor2d(new Tensor2d("temp2", B->d1_, dim_r, dim2_));
	#pragma omp parallel for
        for (int Q = 0; Q < B->d1_; Q++) {
             for (int r = 0; r < dim_r; r++) {
                  for (int q = 0; q < dim2_; q++) {
		       int rq = q + (r*dim2_);

		       if (rb == t_b1) rrb = r;
		       else if (rb == f_b1) rrb = q;

		       if (cb == t_b1) ccb = r;
	               else if (cb == f_b1) ccb = q;

		       int col = ccb + (rrb * B->d3_);

                       temp2->A2d_[Q][rq] = B->A2d_[Q][col];
                  }
             }
        }
	if (delete_b) B.reset();

	m = dim1_;
	n = dim2_;
	k = A->d1_ * dim_r;
	ta = 't';
	tb = 'n';
	nca = m;
	ncb = n;
        ncc = n;

        C_DGEMM(ta, tb, m, n, k, alpha, &(temp1->A2d_[0][0]), nca, &(temp2->A2d_[0][0]), ncb, beta, &(A2d_[0][0]), ncc);
	temp2.reset();
	temp2.reset();
}//

double Tensor2d::get_max_element()
{
  double value = 0.0;
  #pragma omp parallel for
  for (int i=0; i < dim1_; i++) {
       for (int j=0; j < dim2_; j++) {
           if ( fabs(A2d_[i][j]) > value ) value = fabs(A2d_[i][j]);
       }
  }
  return value;
}//



/********************************************************************************************/
/************************** 3d array ********************************************************/
/********************************************************************************************/
Tensor3d::Tensor3d(int d1,int d2, int d3)
{
  A3d_ = NULL;
  dim1_=d1;
  dim2_=d2;
  dim3_=d3;
  memalloc();
}//

Tensor3d::Tensor3d(string name, int d1,int d2, int d3)
{
  A3d_ = NULL;
  dim1_=d1;
  dim2_=d2;
  dim3_=d3;
  name_=name;
  memalloc();
}//

Tensor3d::Tensor3d()
{
  A3d_ = NULL;
  dim1_ = 0;
  dim2_ = 0;
  dim3_ = 0;

}//

Tensor3d::~Tensor3d()
{
  release();
}//

void Tensor3d::memalloc()
{
    if (A3d_) release();
    A3d_ = init_3d_array(dim1_, dim2_, dim3_);
    zero();
}//

void Tensor3d::init(int d1,int d2, int d3)
{
    dim1_=d1;
    dim2_=d2;
    dim3_=d3;
    memalloc();
}//

void Tensor3d::init(string name, int d1,int d2, int d3)
{
    dim1_=d1;
    dim2_=d2;
    dim3_=d3;
    name_=name;
    memalloc();
}//

void Tensor3d::zero()
{
  memset(&(A3d_[0][0][0]), 0, sizeof(double)*dim1_*dim2_*dim3_);
}//

void Tensor3d::print()
{
  if (name_.length()) outfile->Printf( "\n ## %s ##\n", name_.c_str());
  for (int i=0; i<dim1_;i++){
    outfile->Printf( "\n Irrep: %d\n", i+1);
    print_mat(A3d_[i],dim2_,dim3_,"outfile");
  }

}//

void Tensor3d::release()
{
   if (!A3d_) return;
   free_3d_array(A3d_, dim1_, dim2_);
   A3d_ = NULL;
}//

void Tensor3d::set(int h, int i, int j, double value)
{
  A3d_[h][i][j]=value;
}//

double Tensor3d::get(int h, int i, int j)
{
  return A3d_[h][i][j];
}//


/********************************************************************************************/
/************************** 1i array ********************************************************/
/********************************************************************************************/
Tensor1i::Tensor1i(int d1)
{
  A1i_ = NULL;
  dim1_=d1;
  memalloc();
}//

Tensor1i::Tensor1i(string name, int d1)
{
  A1i_ = NULL;
  dim1_=d1;
  name_=name;
  memalloc();
}//

Tensor1i::Tensor1i()
{
  A1i_ = NULL;
  dim1_ = 0;

}//

Tensor1i::~Tensor1i()
{
  release();
}//

void Tensor1i::memalloc()
{
    if (A1i_) release();
    A1i_ = new int[dim1_];
    zero();
}//

void Tensor1i::init(int d1)
{
    dim1_=d1;
    if (A1i_) release();
    A1i_ = new int[dim1_];
}//

void Tensor1i::init(string name, int d1)
{
    dim1_=d1;
    name_=name;
    if (A1i_) release();
    A1i_ = new int[dim1_];
}//

void Tensor1i::zero()
{
    memset(A1i_, 0, sizeof(int)*dim1_);
}//

void Tensor1i::print()
{
  if (name_.length()) outfile->Printf( "\n ## %s ##\n", name_.c_str());
  for (int p=0; p<dim1_; p++){
    outfile->Printf(" %3d %3d \n",p,A1i_[p]);
  }

}//

void Tensor1i::release()
{
   if (!A1i_) return;
   delete [] A1i_;
   A1i_ = NULL;
}//

void Tensor1i::set(int i, int value)
{
  A1i_[i]=value;
}//

int Tensor1i::get(int i)
{
  return A1i_[i];
}//

void Tensor1i::add(const SharedTensor1i& a)
{
    /*
    int *lhs, *rhs;
    size_t size = dim1_;
    if (size) {
        lhs = A1i_;
        rhs = Adum->A1i_;
        for (size_t ij=0; ij<size; ++ij) {
            *lhs += *rhs;
            lhs++; rhs++;
        }
    }
    */
    #pragma omp parallel for
    for (int i=0; i<dim1_; ++i) A1i_[i] += a->A1i_[i];
}//

void Tensor1i::add(int i, int value)
{
  A1i_[i]+=value;
}//

void Tensor1i::subtract(const SharedTensor1i& a)
{
    /*
    int *lhs, *rhs;
    size_t size = dim1_;
    if (size) {
        lhs = A1i_;
        rhs = Adum->A1i_;
        for (size_t ij=0; ij<size; ++ij) {
            *lhs -= *rhs;
            lhs++; rhs++;
        }
    }
    */
    #pragma omp parallel for
    for (int i=0; i<dim1_; ++i) A1i_[i] -= a->A1i_[i];
}//

void Tensor1i::subtract(int i, int value)
{
  A1i_[i]-=value;
}//

/********************************************************************************************/
/************************** 2i array ********************************************************/
/********************************************************************************************/
Tensor2i::Tensor2i(int d1,int d2)
{
  A2i_ = NULL;
  dim1_=d1;
  dim2_=d2;
  memalloc();
}//

Tensor2i::Tensor2i(string name, int d1,int d2)
{
  A2i_ = NULL;
  dim1_=d1;
  dim2_=d2;
  name_=name;
  memalloc();
}//

Tensor2i::Tensor2i()
{
  A2i_ = NULL;
  dim1_ = 0;
  dim2_ = 0;

}//

Tensor2i::~Tensor2i()
{
  release();
}//

void Tensor2i::memalloc()
{
    if (A2i_) release();
    A2i_ = init_int_matrix(dim1_, dim2_);
    zero();
}//

void Tensor2i::init(int d1,int d2)
{
    dim1_=d1;
    dim2_=d2;
    if (A2i_) release();
    A2i_ = init_int_matrix(dim1_, dim2_);
}//

void Tensor2i::init(string name, int d1,int d2)
{
    dim1_=d1;
    dim2_=d2;
    name_=name;
    if (A2i_) release();
    A2i_ = init_int_matrix(dim1_, dim2_);
}//

void Tensor2i::zero()
{
    memset(A2i_[0], 0, sizeof(int)*dim1_*dim2_);
}//

void Tensor2i::zero_diagonal()
{
    if (dim1_ == dim2_ ) {
       for (int i=0; i<dim1_; i++) A2i_[i][i] = 0.0;
   }
}//

void Tensor2i::print()
{
  if (name_.length()) outfile->Printf( "\n ## %s ##\n", name_.c_str());
  print_int_mat(A2i_,dim1_,dim2_,"outfile");

}//

void Tensor2i::print(std::string OutFileRMR)
{
   std::shared_ptr<psi::PsiOutStream> printer=(OutFileRMR=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(OutFileRMR)));
   if (name_.length()) printer->Printf( "\n ## %s ##\n", name_.c_str());
  print_int_mat(A2i_,dim1_,dim2_,OutFileRMR);
}//

void Tensor2i::release()
{
   if (!A2i_) return;
   free_int_matrix(A2i_);
   A2i_ = NULL;
}//

void Tensor2i::set(int i, int j, int value)
{
  A2i_[i][j]=value;
}//

void Tensor2i::set(int **A)
{
      if (A == NULL) return;
      for (int i=0; i<dim1_; ++i) {
        for (int j=0; j<dim2_; ++j) {
          A2i_[i][j] = A[i][j];
        }
      }
}//

double Tensor2i::get(int i, int j)
{
  return A2i_[i][j];
}//

void Tensor2i::add(const SharedTensor2i& a)
{
    /*
    int *lhs, *rhs;
    size_t size = dim1_ * dim2_;
    if (size) {
        lhs = A2i_[0];
        rhs = Adum->A2i_[0];
        for (size_t ij=0; ij<size; ++ij) {
            *lhs += *rhs;
            lhs++; rhs++;
        }
    }
    */
      #pragma omp parallel for
      for (int i=0; i<dim1_; ++i) {
	for (int j=0; j<dim2_; ++j) {
	  A2i_[i][j] += a->A2i_[i][j];
	}
      }
}//

void Tensor2i::add(int i, int j, int value)
{
  A2i_[i][j]+=value;
}//

void Tensor2i::subtract(const SharedTensor2i& a)
{
    /*
    int *lhs, *rhs;
    size_t size = dim1_ * dim2_;
    if (size) {
        lhs = A2i_[0];
        rhs = Adum->A2i_[0];
        for (size_t ij=0; ij<size; ++ij) {
            *lhs -= *rhs;
            lhs++; rhs++;
        }
    }
    */
      #pragma omp parallel for
      for (int i=0; i<dim1_; ++i) {
	for (int j=0; j<dim2_; ++j) {
	  A2i_[i][j] -= a->A2i_[i][j];
	}
      }
}//

void Tensor2i::subtract(int i, int j, int value)
{
  A2i_[i][j]-=value;
}//

SharedTensor2i Tensor2i::transpose()
{
    SharedTensor2i temp = SharedTensor2i(new Tensor2i(dim2_, dim1_));

      for (int i=0; i<dim2_; ++i) {
	for (int j=0; j<dim1_; ++j) {
	  temp->A2i_[i][j] = A2i_[j][i];
	}
      }

    return temp;
}//


void Tensor2i::copy(const SharedTensor2i& Adum)
{
    // Make sure that matrices are in the same size
    bool same = true;
   if (dim2_ != Adum->dim2_ || dim1_ != Adum->dim1_)  same = false;

    if (same == false) {
        release();
        dim1_ = Adum->dim1_;
        dim2_ = Adum->dim2_;
        memalloc();
    }

    // If matrices are in the same size
    ULI length;
    length = (ULI)dim1_ * (ULI)dim2_;
      if (dim1_ != 0 && dim2_ != 0) {
	memcpy(A2i_[0], Adum->A2i_[0], dim1_ * dim2_ * sizeof(int));
      }
}//

void Tensor2i::copy(int **a)
{
    size_t size = dim1_ * dim2_ * sizeof(int);
    if (size) memcpy(&(A2i_[0][0]), &(a[0][0]), size);
}

void Tensor2i::identity()
{
    zero();
    for (int i=0; i<dim1_; ++i) A2i_[i][i] = 1.0;
}//

int Tensor2i::trace()
{
    int value = 0;
    for (int i=0; i<dim1_; ++i) value += A2i_[i][i];
    return value;
}//

int **Tensor2i::to_int_matrix()
{
    int **temp = init_int_matrix(dim1_, dim2_);
    memcpy(&(temp[0][0]), &(A2i_[0][0]), dim1_ * dim2_ * sizeof(int));
    return temp;
}//

/********************************************************************************************/
/************************** 3i array ********************************************************/
/********************************************************************************************/
Tensor3i::Tensor3i(int d1,int d2, int d3)
{
  A3i_ = NULL;
  dim1_=d1;
  dim2_=d2;
  dim3_=d3;
  memalloc();
}//

Tensor3i::Tensor3i(string name, int d1,int d2, int d3)
{
  A3i_ = NULL;
  dim1_=d1;
  dim2_=d2;
  dim3_=d3;
  name_=name;
  memalloc();
}//

Tensor3i::Tensor3i()
{
  A3i_ = NULL;
  dim1_ = 0;
  dim2_ = 0;
  dim3_ = 0;

}//

Tensor3i::~Tensor3i()
{
  release();
}//

void Tensor3i::memalloc()
{
  if (A3i_) release();
  A3i_ = (int ***) malloc(dim1_ * sizeof(int **));
  for(int i=0; i < dim1_; i++) {
    A3i_[i] = (int **) malloc(dim2_ * sizeof(int *));
    for(int j=0; j < dim2_; j++) {
      A3i_[i][j] = (int *) malloc(dim3_ * sizeof(int));
      for(int k=0; k < dim3_; k++) {
        A3i_[i][j][k] = 0.0;
      }
    }
  }
}//

void Tensor3i::init(int d1,int d2, int d3)
{
    dim1_=d1;
    dim2_=d2;
    dim3_=d3;
    if (A3i_) release();
    memalloc();
}//

void Tensor3i::init(string name, int d1, int d2, int d3)
{
    dim1_=d1;
    dim2_=d2;
    dim3_=d3;
    name_=name;
    if (A3i_) release();
    memalloc();
}//

void Tensor3i::zero()
{
  memset(&(A3i_[0][0][0]), 0, sizeof(int)*dim1_*dim2_*dim3_);
}//

void Tensor3i::print()
{
  if (name_.length()) outfile->Printf( "\n ## %s ##\n", name_.c_str());
  for (int i=0; i<dim1_;i++){
    outfile->Printf( "\n Irrep: %d\n", i+1);
    print_int_mat(A3i_[i],dim2_,dim3_,"outfile");
  }

}//

void Tensor3i::release()
{
   if (!A3i_) return;
   for(int i=0; i < dim1_; i++) {
       for(int j=0; j < dim2_; j++) {
           free(A3i_[i][j]);
       }
   }
   for(int i=0; i < dim1_; i++) free(A3i_[i]);
   free(A3i_);
   A3i_ = NULL;
}//

void Tensor3i::set(int h, int i, int j, int value)
{
  A3i_[h][i][j]=value;
}//

int Tensor3i::get(int h, int i, int j)
{
  return A3i_[h][i][j];
}//

/********************************************************************************************/
/********************************************************************************************/
}} // End Namespaces

