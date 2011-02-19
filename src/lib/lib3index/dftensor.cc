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

ThreeIndexChunk::ThreeIndexChunk(shared_ptr<PSIO> psio, 
                const std::string& name,
                int slow_size,
                int middle_size,
                int fast_size,
                unsigned long int memory 
                ) : psio_(psio), 
                core_tensor_(NULL),
                name_(name),
                slow_size_(slow_size),
                middle_size_(middle_size),
                fast_size_(fast_size_),
                memory_(memory),
                is_core_(false)    
{
    tensor_size_ = slow_size_*middle_size_*(unsigned long int)fast_size_;
    is_cached_ = tensor_size_ <= memory ; 
    unsigned long int max_required;
    if (is_cached_) {
        max_rows_ = slow_size_;
        max_required = tensor_size_;
    } else {
        max_rows_ = (memory / (middle_size*(unsigned long int)fast_size));
        max_required = max_rows_*middle_size*(unsigned long int)fast_size;
    }
    chunk_ = new double[max_required];

    if (is_cached_) {
        address_ = PSIO_ZERO;
        psio_->read(PSIF_DF_TENSOR, name_.c_str(), (char*)(chunk_), tensor_size_*sizeof(double), address_, &address_); 
        current_row_ = 0;
        current_rows_ = slow_size_;
        is_done_ = false;
    } else { 
        current_row_ = 0;
        current_rows_ = max_rows_;
        is_done_ = false;
    }
}
ThreeIndexChunk::ThreeIndexChunk(double* core_tensor, 
                const std::string& name,
                int slow_size,
                int middle_size,
                int fast_size,
                unsigned long int memory 
                ) : core_tensor_(core_tensor),
                name_(name),
                slow_size_(slow_size),
                middle_size_(middle_size),
                fast_size_(fast_size_),
                memory_(memory),    
                is_core_(true),    
                is_cached_(false)    
{
    tensor_size_ = slow_size_*middle_size_*(unsigned long int)fast_size_;
    if (is_cached_) {
        max_rows_ = slow_size_;
    } else {
        max_rows_ = (memory / (middle_size*(unsigned long int)fast_size));
    }
    chunk_ = &core_tensor_[0]; 

    current_row_ = 0;
    current_rows_ = max_rows_;
    is_done_ = false;
}
ThreeIndexChunk::~ThreeIndexChunk()
{
    if (!is_core_)
        delete[] chunk_;
}
void ThreeIndexChunk::reset()
{
    current_row_ = 0L;
    current_rows_ = max_rows_;
    is_done_ = false;

    if (is_core_) {
        chunk_ = core_tensor_;
    } else {
        if (!is_cached_) {
            address_ = PSIO_ZERO;
            psio_->read(PSIF_DF_TENSOR, name_.c_str(), (char*)(chunk_), max_rows_*middle_size_*fast_size_*sizeof(double), address_, &address_); 
        }    
    } 
}
void ThreeIndexChunk::next()
{
    // Check if complete
    if (current_row_ + current_rows_ >= slow_size_) {
        current_row_ = slow_size_;
        current_rows_ = 0;
        is_done_ = true;
        return;
    }

    current_row_ += current_rows_;
    current_rows_ = (current_row_ + max_rows_ >= slow_size_ ? slow_size_ - current_row_ : max_rows_);
        
    if (is_core_) {
        chunk_ = &core_tensor_[current_row_*middle_size_*(unsigned long int)fast_size_];
    } else {
        if (!is_cached_) {
            address_ = PSIO_ZERO;
            psio_->read(PSIF_DF_TENSOR, name_.c_str(), (char*)(chunk_), current_rows_*middle_size_*fast_size_*sizeof(double), address_, &address_); 
        }    
    }
}        
DFTensor::DFTensor(shared_ptr<PSIO> psio, shared_ptr<BasisSet> primary, shared_ptr<BasisSet> aux, double schwarz) :
    psio_(psio), primary_(primary), auxiliary_(aux), schwarz_(new SchwarzSieve(primary, schwarz)), metric_(new FittingMetric(aux))
{
}
DFTensor::~DFTensor()
{
}
void DFTensor::common_init()
{
}
void DFTensor::form_MO_disk(shared_ptr<Matrix> Co, shared_ptr<Matrix> Cv, const std::string& algorithm, double cond)
{
    is_core_ = false;
    nocc_ = Co->colspi()[0];
    nvir_ = Cv->colspi()[0];
}
shared_ptr<ThreeIndexChunk> DFTensor::get_oo_iterator(unsigned long int memory)
{
    return shared_ptr<ThreeIndexChunk>(new ThreeIndexChunk(
        psio_,
        "DF OO Integrals",
        auxiliary_->nbf(),
        nocc_,
        nocc_,
        memory
        ));
}
shared_ptr<ThreeIndexChunk> DFTensor::get_ov_iterator(unsigned long int memory)
{
    return shared_ptr<ThreeIndexChunk>(new ThreeIndexChunk(
        psio_,
        "DF OV Integrals",
        auxiliary_->nbf(),
        nocc_,
        nvir_,
        memory
        ));
}
shared_ptr<ThreeIndexChunk> DFTensor::get_vv_iterator(unsigned long int memory)
{
    return shared_ptr<ThreeIndexChunk>(new ThreeIndexChunk(
        psio_,
        "DF VV Integrals",
        auxiliary_->nbf(),
        nocc_,
        nvir_,
        memory
        ));
}


}

