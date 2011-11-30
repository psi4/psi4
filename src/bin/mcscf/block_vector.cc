#include <cstdio>

#include <libutil/libutil.h>

#include "block_vector.h"
#include "vector_base.h"

#include <psi4-dec.h>

namespace psi{ namespace mcscf{

extern MemoryManager* memory_manager;

BlockVector::BlockVector()
 : nirreps_(0) ,ref_(0), vector_base_(0), rows_size_(0), rows_offset_(0)
{
}


BlockVector::BlockVector(std::string label, int nirreps, size_t*& rows_size)
 : label_(label), nirreps_(nirreps) ,ref_(0), vector_base_(0), rows_size_(0), rows_offset_(0)
{
  startup(label,nirreps,rows_size);
}

BlockVector::BlockVector(std::string label, int nirreps, int*& rows_size)
 : label_(label), nirreps_(nirreps) ,ref_(0), vector_base_(0), rows_size_(0), rows_offset_(0)
{
  startup(label,nirreps,rows_size);
}

BlockVector::BlockVector(std::string label, int nirreps, vecint& rows_size)
 : label_(label), nirreps_(nirreps) ,ref_(0), vector_base_(0), rows_size_(0), rows_offset_(0)
{
  startup(label,nirreps,rows_size);
}

BlockVector::~BlockVector()
{
  cleanup();
}

void BlockVector::startup(std::string label, int nirreps, size_t*& rows_size)
{
  vector_base_ = new VectorBase*[nirreps_];
  // Allocate the blocks
  for(int h = 0; h < nirreps_; ++h){
    vector_base_[h] = new VectorBase(rows_size[h]);
  }

  allocate1(size_t,rows_size_,nirreps);
  allocate1(size_t,rows_offset_,nirreps);
  // Compute the offsets
  rows_offset_[0] = 0;
  for(int h = 1; h < nirreps; ++h){
    rows_size_[h] = rows_size[h];
    rows_offset_[h] = rows_offset_[h-1] + rows_size[h-1];
  }
}

void BlockVector::startup(std::string label, int nirreps, int*& rows_size)
{
  vector_base_ = new VectorBase*[nirreps_];
  // Allocate the blocks
  for(int h = 0; h < nirreps_; ++h){
    vector_base_[h] = new VectorBase(rows_size[h]);
  }

  allocate1(size_t,rows_size_,nirreps);
  allocate1(size_t,rows_offset_,nirreps);
  // Compute the offsets
  rows_offset_[0] = 0;
  for(int h = 1; h < nirreps; ++h){
    rows_size_[h] = rows_size[h];
    rows_offset_[h] = rows_offset_[h-1] + rows_size[h-1];
  }
}

void BlockVector::startup(std::string label, int nirreps, vecint& rows_size)
{
  vector_base_ = new VectorBase*[nirreps_];
  // Allocate the blocks
  for(int h = 0; h < nirreps_; ++h){
    vector_base_[h] = new VectorBase(rows_size[h]);
  }

  allocate1(size_t,rows_size_,nirreps);
  allocate1(size_t,rows_offset_,nirreps);
  // Compute the offsets
  rows_offset_[0] = 0;
  for(int h = 1; h < nirreps; ++h){
    rows_size_[h] = rows_size[h];
    rows_offset_[h] = rows_offset_[h-1] + rows_size[h-1];
  }
}

void BlockVector::cleanup()
{
  if(vector_base_){
    for(int h = 0; h < nirreps_; ++h){
      delete vector_base_[h];
    }
    delete[] vector_base_;
  }
  release1(rows_size_);
  release1(rows_offset_);
}

void BlockVector::print()
{
  fprintf(outfile,"\n\n  ## %s ##\n",label_.c_str());
  for(int h = 0; h < nirreps_; ++h){
    vector_base_[h]->print();
  }
  fflush(outfile);
}


void BlockVector::copy(BlockVector& source)
{
  for(int h = 0; h < nirreps_; ++h){
    vector_base_[h]->copy(*source.vector_base_[h]);
  }
}

}}

// void BlockVector::zero()
// {
//   for(int h = 0; h < nirreps_; ++h)
//     vector_base_[h]->zero();
// }
//

//
// void BlockVector::scale(double factor)
// {
//   for(int h = 0; h < nirreps_; ++h)
//     vector_base_[h]->scale(factor);
// }
//
// void BlockVector::transpose()
// {
//   for(int h = 0; h < nirreps_; ++h)
//     vector_base_[h]->transpose();
// }
//
// void BlockVector::multiply(bool transpose_A, bool transpose_B, BlockVector* A, BlockVector* B)
// {
//   for(int h = 0; h < nirreps_; ++h)
//     getVectorBase(h)->multiply(transpose_A,
//                          transpose_B,
//                          A->getVectorBase(h),
//                          B->getVectorBase(h));
// }
//
// void BlockVector::diagonalize(BlockVector* eigenvectors,BlockVector* eigenvalues)
// {
//   for(int h = 0; h < nirreps_; ++h)
//     getVectorBase(h)->diagonalize(eigenvectors->getVectorBase(h),
//                             eigenvalues->getVectorBase(h));
// }
//
// BlockVector& BlockVector::operator=(BlockVector& rhs)
// {
//   if(this == &rhs){
//     return(*this);
//   }
//
//   cleanup();
//   startup(rhs.label_,
//           rhs.nirreps_,
//           rhs.rows_size_,
//           rhs.cols_size_);
//
//   for(int h=0; h < nirreps_; ++h){
//     if(rows_size_[h] * cols_size_[h]>0){
//       for(int i = 0; i < rows_size_[h]; ++i)
//         for(int j = 0; j < cols_size_[h]; ++j)
//           vector_base_[h]->set(i,j, rhs.vector_base_[h]->get(i,j) );
//     }
//   }
//   return(*this);
// }

