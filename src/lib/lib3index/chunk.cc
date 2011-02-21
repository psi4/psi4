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
                         unsigned long int memory) :
                         psio_(psio), 
                         unit_(unit),
                         name_(name), 
                         slow_size_(slow_size), 
                         fast_size_(fast_size), 
                         memory_(memory),
                         touched_(false),
                         block_(0)
{
    tensor_size_ = slow_size_*(unsigned long int) fast_size_;
    max_rows_ = memory_ / (unsigned long int) fast_size_;
    if (max_rows_ > slow_size_)
        max_rows_ = slow_size_;
    chunk_ = shared_ptr<Matrix>(new Matrix("Tensor Chunk of " + name_, max_rows_, fast_size_));
     
    int nblocks = slow_size_ / max_rows_;
    if (max_rows_ * nblocks < slow_size_)
        nblocks++;
    block_starts_.resize(nblocks);
    block_sizes_.resize(nblocks);
    
    // TODO spread the gimp out
    block_starts_[0] = 0;
    block_sizes_[0] = max_rows_;
    for (int Q = 1; Q < nblocks; Q++) {
        block_starts_[Q] = block_starts_[Q - 1] + max_rows_;
        if (block_starts_[Q] + max_rows_ >= slow_size_)
            block_sizes_[Q] = slow_size_ - block_starts_[Q];
        else 
            block_sizes_[Q] = max_rows_;
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
    psio_->read(unit_,name_.c_str(),(char*) chunk_->pointer()[0], block_sizes_[block_]*sizeof(double), block_addr, &block_addr);
}

}

