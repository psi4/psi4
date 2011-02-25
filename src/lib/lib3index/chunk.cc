#include "3index.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libmints/mints.h>

//MKL Header
#ifdef HAVE_MKL
#include <mkl.h>
#endif

//OpenMP Header
//_OPENMP is defined by the compiler if it exists
#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;
using namespace psi;

namespace psi { 

TensorChunk::TensorChunk(shared_ptr<PSIO> psio,
                         unsigned int unit,
                         const std::string& name,
                         int slow_size,
                         int fast_size,
                         unsigned long int memory,
                         int min_rows) :
                         psio_(psio), 
                         unit_(unit),
                         name_(name), 
                         slow_size_(slow_size), 
                         fast_size_(fast_size), 
                         memory_(memory),
                         min_rows_(min_rows),
                         touched_(false),
                         block_(-1)
{
    if (slow_size_ % min_rows_ != 0)
        throw PSIEXCEPTION("TensorChunk: Grain size (minimum required rows) is not a factor of the slow index"); 
    
    tensor_size_ = slow_size_*(unsigned long int) fast_size_;
    max_rows_ = memory_ / (unsigned long int) fast_size_;
    if (max_rows_ < min_rows_)
        throw PSIEXCEPTION("TensorChunk: Minimum required rows is too large for memory provided to TensorChunk");
    if (max_rows_ > slow_size_)
        max_rows_ = slow_size_;
    max_rows_ = (max_rows_ / min_rows_) * min_rows_;

    chunk_ = shared_ptr<Matrix>(new Matrix("Tensor Chunk of " + name_, max_rows_, fast_size_));
     
    int nblocks = slow_size_ / max_rows_;
    if (max_rows_ * nblocks < slow_size_)
        nblocks++;
    block_starts_.resize(nblocks);
    block_sizes_.resize(nblocks);
    
    // Naive distribution
    block_starts_[0] = 0;
    block_sizes_[0] = max_rows_;
    int gimp_delta = 0;
    for (int Q = 1; Q < nblocks; Q++) {
        block_starts_[Q] = block_starts_[Q - 1] + max_rows_;
        if (block_starts_[Q] + max_rows_ >= slow_size_) {
            block_sizes_[Q] = slow_size_ - block_starts_[Q];
            gimp_delta = (block_sizes_[Q - 1] - block_sizes_[Q]) / min_rows_; 
        } else {
            block_sizes_[Q] = max_rows_;
        }
    }
   
    // Now Level the gimp out
    for (int Q = 0; Q < gimp_delta - 1; Q++) {
        block_sizes_[nblocks - 2- Q] -= min_rows_;
        block_starts_[nblocks - 2 - Q] -= (gimp_delta - 1 - Q) * min_rows_;
    }
 
}
TensorChunk::~TensorChunk()
{
}
void TensorChunk::read_block(int index, bool cache)
{
    block_ = index;
    if (block_ >= block_sizes_.size() || block_ < 0)
        throw PSIEXCEPTION("TensorChunk::read_block: Invalid block index");
    if (cache && touched_ && block_sizes_.size() == 1)
        return;
    touched_ = true;

    psio_address block_addr = psio_get_address(PSIO_ZERO,(ULI)(block_starts_[block_]*(ULI)fast_size_*sizeof(double))); 
    psio_->read(unit_,name_.c_str(),(char*) chunk_->pointer()[0], block_sizes_[block_]*fast_size_*sizeof(double), block_addr, &block_addr);
}

}

