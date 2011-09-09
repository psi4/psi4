#include <cstring>
#include <iostream>

#include <libciomr/libciomr.h>
#include <libmoinfo/libmoinfo.h>
#include <libqt/qt.h>
#include <libutil/libutil.h>
#include <libutil/memory_manager.h>

#include "special_matrices.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;

MatrixBase::MatrixBase(size_t nrows_, size_t ncols_) : nrows(nrows_),ncols(ncols_),matrix(0)
{
  allocate2(double,matrix,nrows,ncols);
  zero();
}

MatrixBase::~MatrixBase()
{
  release2(matrix);
}

void MatrixBase::zero()
{
  size_t nbites = static_cast<size_t>(sizeof(double)) * nrows * ncols;
  if(nbites > 0)
    memset(&(matrix[0][0]),'\0', static_cast<size_t>(sizeof(double)) * nrows * ncols);
}

void MatrixBase::print()
{
  if(nrows * ncols > 0){
    for(size_t p = 0; p < nrows; ++p){
      fprintf(outfile,"\n  ");
      for(size_t q = 0; q < ncols; ++q){
        fprintf(outfile,"%10.6f",matrix[p][q]);
      }
    }
  }
}

double MatrixBase::norm()
{
  double norm = 0.0;
  for(size_t p = 0; p < nrows; ++p){
    for(size_t q = 0; q < ncols; ++q){
      norm += matrix[p][q] * matrix[p][q];
    }
  }
  return norm;
}

// (this) = alpha (this) + beta A
void MatrixBase::add(MatrixBase* A, double alpha, double beta)
{
  if(nrows * ncols > 0){
    double** a = A->get_matrix();
    if(alpha != 1.0)
      C_DSCAL(nrows * ncols,alpha,&(matrix[0][0]),1);
    C_DAXPY(nrows * ncols,beta,&(a[0][0]),1,&(matrix[0][0]),1);
  }
//  for(int p = 0; p < nrows; ++p){
//    for(int q = 0; q < ncols; ++q){
//      matrix[p][q] = alpha * matrix[p][q] + beta * a[p][q];
//    }
//  }
}

// (this) =  alpha A * B + beta (this)
void MatrixBase::contract(MatrixBase* A, MatrixBase* B, double const alpha, double const beta)
{
  size_t max_r = A->get_ncols();
  if((max_r != 0) and (nrows != 0) and (ncols != 0))
    C_DGEMM('n','t',nrows,ncols,max_r,alpha,&(A->get_matrix()[0][0]),max_r,&(B->get_matrix()[0][0]),max_r,beta,&(matrix[0][0]),ncols);
  else if(fabs(beta) < 1.0e-9)
    zero();
  // BUG: There was a tricky bug here.  When alpha is 0.0 and you skip DGEMM the matrix doesn't get set to zero.
}

// (this) =  alpha A * B + beta (this)
void MatrixBase::multiply(MatrixBase* A, MatrixBase* B, double alpha, double beta)
{
  double** a = A->get_matrix();
  double** b = B->get_matrix();
  size_t max_r = A->get_ncols();
  for(size_t p = 0; p < nrows; ++p){
    for(size_t q = 0; q < ncols; ++q){
      double sum = 0.0;
      for(size_t r = 0; r < max_r; ++r){
        sum += a[p][r] * b[q][r];
      }
      matrix[p][q] = alpha * matrix[p][q] + beta * sum;
    }
  }
}

//  Arrange the block according to the irrep of the columns (it makes it easier for matrix multiplication)
BlockMatrix::BlockMatrix(int nirreps_,std::vector<size_t>& rows_size_, std::vector<size_t>& cols_size_,int sym_)
{
  // Copy values
  nirreps   = nirreps_;
  sym       = sym_;

  // Allocate and compute the offsets
  rows_size.assign(nirreps,0);
  cols_size.assign(nirreps,0);

  // Make sure the product of the row and cols belong to irrep = sym
  for(int h = 0; h < nirreps; ++h){
    rows_size[h] = rows_size_[h];
    cols_size[h] = cols_size_[h ^ sym];
  }

  // Allocate and compute the offsets
  rows_offset.assign(nirreps,0);
  cols_offset.assign(nirreps,0);

  rows_offset[0] = 0;
  cols_offset[0] = 0;
  for(int h = 1; h < nirreps; ++h){
    rows_offset[h] = rows_offset[h-1] + rows_size[h-1];
    cols_offset[h] = cols_offset[h-1] + cols_size[h-1];
  }

  // Allocate the blocks
  blocks = new MatrixBase*[nirreps];
  for(int h = 0; h < nirreps; ++h){
    blocks[h] = new MatrixBase(rows_size[h],cols_size[h]);
  }
}

BlockMatrix::~BlockMatrix()
{
  // Deallocate the blocks
  for(int h = 0; h < nirreps; ++h){
    delete blocks[h];
  }
  delete[] blocks;
}

void BlockMatrix::print()
{
  // Deallocate the blocks
  for(int h = 0; h < nirreps; ++h){
    fprintf(outfile,"\n    Block %d",h);
    blocks[h]->print();
  }
}

double BlockMatrix::norm()
{
  double norm = 0.0;
  // Deallocate the blocks
  for(int h = 0; h < nirreps; ++h){
    norm += blocks[h]->norm();
  }
  return norm;
}

void BlockMatrix::zero()
{
  // Loop over the irrep of the summation index
  for(int h = 0; h < moinfo->get_nirreps(); ++h){
    blocks[h]->zero();
  }
}

// (this) = alpha (this) + beta A
void BlockMatrix::add(BlockMatrix* A, double alpha, double beta)
{
  // Loop over the irrep of the summation index
  for(int h = 0; h < moinfo->get_nirreps(); ++h){
    blocks[h]->add(A->blocks[h],alpha,beta);
  }
}

// (this) = alpha A * B + beta (this)
void BlockMatrix::contract(BlockMatrix* A, BlockMatrix* B, double alpha, double beta)
{
  // Loop over the irrep of the summation index
  for(int h = 0; h < moinfo->get_nirreps(); ++h){
    int row_sym = h;
    int col_sym = h ^ sym;
    blocks[row_sym]->contract(A->blocks[row_sym],B->blocks[col_sym],alpha,beta);
  }
}

// (this) = alpha A * B + beta (this)
void BlockMatrix::multiply(BlockMatrix* A, BlockMatrix* B, double alpha, double beta)
{
  // Loop over the irrep of the summation index
  for(int h = 0; h < moinfo->get_nirreps(); ++h){
    int row_sym = h;
    int col_sym = h ^ sym;
    blocks[row_sym]->multiply(A->blocks[row_sym],B->blocks[col_sym],alpha,beta);
  }
}

void BlockMatrix::cyclical_permutation_1_2(BlockMatrix* A, CCIndex* pqr_index,CCIndex* p_index,CCIndex* qr_index)
{
  CCIndexIterator pqr(pqr_index,sym);
  for(pqr.first(); !pqr.end(); pqr.next()){
    int p_sym    = p_index->get_tuple_irrep(pqr.ind_abs<0>());
    int q_sym    = p_index->get_tuple_irrep(pqr.ind_abs<1>());
    int r_sym    = p_index->get_tuple_irrep(pqr.ind_abs<2>());

    size_t p_rel = p_index->get_tuple_rel_index(pqr.ind_abs<0>());
    size_t q_rel = p_index->get_tuple_rel_index(pqr.ind_abs<1>());
    size_t r_rel = p_index->get_tuple_rel_index(pqr.ind_abs<2>());

    size_t qr_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<1>(),pqr.ind_abs<2>());
    size_t pr_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<0>(),pqr.ind_abs<2>());
    size_t qp_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<1>(),pqr.ind_abs<0>());

    double value = A->blocks[p_sym]->get(p_rel,qr_rel)
                 - A->blocks[q_sym]->get(q_rel,pr_rel)
                 - A->blocks[r_sym]->get(r_rel,qp_rel);
    blocks[p_sym]->set(p_rel,qr_rel,value);
  }
}

void BlockMatrix::a_b_permutation_1_2(BlockMatrix* A, CCIndex* pqr_index,CCIndex* p_index,CCIndex* qr_index)
{
  CCIndexIterator pqr(pqr_index,sym);
  for(pqr.first(); !pqr.end(); pqr.next()){
    int p_sym    = p_index->get_tuple_irrep(pqr.ind_abs<0>());
    int q_sym    = p_index->get_tuple_irrep(pqr.ind_abs<1>());

    size_t p_rel = p_index->get_tuple_rel_index(pqr.ind_abs<0>());
    size_t q_rel = p_index->get_tuple_rel_index(pqr.ind_abs<1>());

    size_t qr_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<1>(),pqr.ind_abs<2>());
    size_t pr_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<0>(),pqr.ind_abs<2>());

    double value = A->blocks[p_sym]->get(p_rel,qr_rel) - A->blocks[q_sym]->get(q_rel,pr_rel);
    blocks[p_sym]->set(p_rel,qr_rel,value);
  }
}

void BlockMatrix::a_b_permutation(CCIndex* pqr_index,CCIndex* p_index,CCIndex* qr_index)
{
  CCIndexIterator pqr(pqr_index,sym);
  for(pqr.first(); !pqr.end(); pqr.next()){
    if(pqr.ind_abs<0>() >= pqr.ind_abs<1>()){
      int p_sym    = p_index->get_tuple_irrep(pqr.ind_abs<0>());
      int q_sym    = p_index->get_tuple_irrep(pqr.ind_abs<1>());

      size_t p_rel = p_index->get_tuple_rel_index(pqr.ind_abs<0>());
      size_t q_rel = p_index->get_tuple_rel_index(pqr.ind_abs<1>());

      size_t qr_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<1>(),pqr.ind_abs<2>());
      size_t pr_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<0>(),pqr.ind_abs<2>());

      double value = blocks[p_sym]->get(p_rel,qr_rel) - blocks[q_sym]->get(q_rel,pr_rel);
      blocks[p_sym]->set(p_rel,qr_rel,value);
      blocks[q_sym]->set(q_rel,pr_rel,-value);
    }
  }
}


/**
 * Sets the block matrix B to B(pqr) = z * A(pqr) + a * M(pqr) + b * M(prq) + c * M(qpr) + d * M(qrp) + e * M(rpq) + f * M(rqp)
 */
void BlockMatrix::add_permutation_1_2(double z,BlockMatrix* A, CCIndex* pqr_index,CCIndex* p_index,CCIndex* qr_index,
    double a,double b,double c,double d,double e,double f)
{
  CCIndexIterator pqr(pqr_index,sym);
  for(pqr.first(); !pqr.end(); pqr.next()){
    int p_sym    = p_index->get_tuple_irrep(pqr.ind_abs<0>());
    int q_sym    = p_index->get_tuple_irrep(pqr.ind_abs<1>());
    int r_sym    = p_index->get_tuple_irrep(pqr.ind_abs<2>());

    size_t p_rel = p_index->get_tuple_rel_index(pqr.ind_abs<0>());
    size_t q_rel = p_index->get_tuple_rel_index(pqr.ind_abs<1>());
    size_t r_rel = p_index->get_tuple_rel_index(pqr.ind_abs<2>());

    size_t qr_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<1>(),pqr.ind_abs<2>());
    size_t rq_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<2>(),pqr.ind_abs<1>());

    size_t pr_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<0>(),pqr.ind_abs<2>());
    size_t rp_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<2>(),pqr.ind_abs<0>());

    size_t pq_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<0>(),pqr.ind_abs<1>());
    size_t qp_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<1>(),pqr.ind_abs<0>());

    double value = z * blocks[p_sym]->get(p_rel,qr_rel)
                 + a * A->blocks[p_sym]->get(p_rel,qr_rel)
                 + b * A->blocks[p_sym]->get(p_rel,rq_rel)
                 + c * A->blocks[q_sym]->get(q_rel,pr_rel)
                 + d * A->blocks[q_sym]->get(q_rel,rp_rel)
                 + e * A->blocks[r_sym]->get(r_rel,pq_rel)
                 + f * A->blocks[r_sym]->get(r_rel,qp_rel);
    blocks[p_sym]->set(p_rel,qr_rel,value);
  }
}

/**
 * Sets the block matrix B to B(pqr) = z * A(pqr) + a * M(pqr) + b * M(prq) + c * M(qpr) + d * M(qrp) + e * M(rpq) + f * M(rqp)
 */
void BlockMatrix::add_acb(double z,BlockMatrix* A, CCIndex* pqr_index,CCIndex* p_index,CCIndex* qr_index,double a)
{
  CCIndexIterator pqr(pqr_index,sym);
  for(pqr.first(); !pqr.end(); pqr.next()){
    int p_sym    = p_index->get_tuple_irrep(pqr.ind_abs<0>());

    size_t p_rel = p_index->get_tuple_rel_index(pqr.ind_abs<0>());
    size_t rq_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<2>(),pqr.ind_abs<1>());
    size_t qr_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<1>(),pqr.ind_abs<2>());

    double value = z * blocks[p_sym]->get(p_rel,qr_rel) + a * A->blocks[p_sym]->get(p_rel,rq_rel);
    blocks[p_sym]->set(p_rel,qr_rel,value);
  }
}

/**
 * Sets the block matrix B to B(pqr) = z * A(pqr) + a * M(pqr) + b * M(prq) + c * M(qpr) + d * M(qrp) + e * M(rpq) + f * M(rqp)
 */
void BlockMatrix::add_cab(double z,BlockMatrix* A, CCIndex* pqr_index,CCIndex* p_index,CCIndex* qr_index,double a)
{
  CCIndexIterator pqr(pqr_index,sym);
  for(pqr.first(); !pqr.end(); pqr.next()){
    int p_sym    = p_index->get_tuple_irrep(pqr.ind_abs<0>());
    int r_sym    = p_index->get_tuple_irrep(pqr.ind_abs<2>());

    size_t p_rel = p_index->get_tuple_rel_index(pqr.ind_abs<0>());
    size_t r_rel = p_index->get_tuple_rel_index(pqr.ind_abs<2>());
    size_t pq_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<0>(),pqr.ind_abs<1>());
    size_t qr_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<1>(),pqr.ind_abs<2>());

    double value = z * blocks[p_sym]->get(p_rel,qr_rel) + a * A->blocks[r_sym]->get(r_rel,pq_rel);
    blocks[p_sym]->set(p_rel,qr_rel,value);
  }
}
///**
// * Sets the block matrix B to B(pqr) = z * A(pqr) + a * M(pqr) + b * M(prq) + c * M(qpr) + d * M(qrp) + e * M(rpq) + f * M(rqp)
// */
//void BlockMatrix::add_permutation_1_2_abc(double z,BlockMatrix* A, CCIndex* pqr_index,CCIndex* p_index,CCIndex* qr_index,
//    double a,double b,double c,double d,double e,double f)
//{
//  CCIndexIterator pqr(pqr_index,sym);
//  for(pqr.first(); !pqr.end(); pqr.next()){
//    int p_sym    = p_index->get_tuple_irrep(pqr.ind_abs<0>());
//    int q_sym    = p_index->get_tuple_irrep(pqr.ind_abs<1>());
//    int r_sym    = p_index->get_tuple_irrep(pqr.ind_abs<2>());
//
//    size_t p_rel = p_index->get_tuple_rel_index(pqr.ind_abs<0>());
//    size_t q_rel = p_index->get_tuple_rel_index(pqr.ind_abs<1>());
//    size_t r_rel = p_index->get_tuple_rel_index(pqr.ind_abs<2>());
//
//    size_t qr_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<1>(),pqr.ind_abs<2>());
//    size_t rq_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<2>(),pqr.ind_abs<1>());
//
//    size_t pr_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<0>(),pqr.ind_abs<2>());
//    size_t rp_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<2>(),pqr.ind_abs<0>());
//
//    size_t pq_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<0>(),pqr.ind_abs<1>());
//    size_t qp_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<1>(),pqr.ind_abs<0>());
//
//    double value = z * blocks[p_sym]->get(p_rel,qr_rel)
//                 + a * A->blocks[p_sym]->get(p_rel,qr_rel)
//                 + b * A->blocks[p_sym]->get(p_rel,rq_rel)
//                 + c * A->blocks[q_sym]->get(q_rel,pr_rel)
//                 + d * A->blocks[q_sym]->get(q_rel,rp_rel)
//                 + e * A->blocks[r_sym]->get(r_rel,pq_rel)
//                 + f * A->blocks[r_sym]->get(r_rel,qp_rel);
//    blocks[p_sym]->set(p_rel,qr_rel,value);
//  }
//}
//
//void BlockMatrix::add_a_b_permutation_1_2(BlockMatrix* A, CCIndex* pqr_index,CCIndex* p_index,CCIndex* qr_index)
//{
//  CCIndexIterator pqr(pqr_index,sym);
//  while(++pqr){
//    int p_sym    = p_index->get_tuple_irrep(pqr.ind_abs<0>());
//    int q_sym    = p_index->get_tuple_irrep(pqr.ind_abs<1>());
//
//    size_t p_rel = p_index->get_tuple_rel_index(pqr.ind_abs<0>());
//    size_t q_rel = p_index->get_tuple_rel_index(pqr.ind_abs<1>());
//
//    size_t qr_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<1>(),pqr.ind_abs<2>());
//    size_t pr_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<0>(),pqr.ind_abs<2>());
//
//    double value = A->blocks[p_sym]->get(p_rel,qr_rel)
//                 - A->blocks[q_sym]->get(q_rel,pr_rel);
//    blocks[p_sym]->add(p_rel,qr_rel,value);
//  }
//}

void BlockMatrix::add_c_ab_permutation_1_2(BlockMatrix* A, CCIndex* pqr_index,CCIndex* p_index,CCIndex* qr_index)
{
  CCIndexIterator pqr(pqr_index,sym);
  for(pqr.first(); !pqr.end(); pqr.next()){
    int p_sym    = p_index->get_tuple_irrep(pqr.ind_abs<0>());
    int r_sym    = p_index->get_tuple_irrep(pqr.ind_abs<2>());

    size_t p_rel = p_index->get_tuple_rel_index(pqr.ind_abs<0>());
    size_t r_rel = p_index->get_tuple_rel_index(pqr.ind_abs<2>());

    size_t qr_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<1>(),pqr.ind_abs<2>());
    size_t pq_rel = qr_index->get_tuple_rel_index(pqr.ind_abs<0>(),pqr.ind_abs<1>());

    double value = A->blocks[r_sym]->get(r_rel,pq_rel);
    blocks[p_sym]->add(p_rel,qr_rel,value);
  }
}

IndexMatrix::IndexMatrix()
{}

IndexMatrix::~IndexMatrix()
{
  for(BMMap::iterator iter=matrices.begin();iter!=matrices.end();++iter){
    delete iter->second;
  }
}

void IndexMatrix::add_block_matrix(size_t index,int ref,BlockMatrix* block_matrix)
{
    matrices[std::make_pair(index,ref)] = block_matrix;
}

BlockMatrix* IndexMatrix::get_block_matrix(size_t index,int ref)
{
    BMMap::iterator iter = matrices.find(std::make_pair(index,ref));
  if(iter!=matrices.end()){
      return matrices[std::make_pair(index,ref)];
  }

  fprintf(outfile,"\n  Couldn't find element!");
  fflush(outfile);
  abort();
  return 0;
}

void IndexMatrix::print()
{
  for(BMMap::iterator iter = matrices.begin(); iter != matrices.end(); ++iter){
    fprintf(outfile,"\n  Index = %4d Ref = %d",static_cast<int>(iter->first.first),static_cast<int>(iter->first.second));
    iter->second->print();
  }
}

}}
