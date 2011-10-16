#include "permutation.h"
#include "sort.h"
#include "env.h"
#include "index.h"
#include "exception.h"
#include "malloc.h"
#include "filler.h"
#include "data.h"
#include "runtime.h"
#include "gigmatrix.h"

using namespace yeti;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

UnitEstimaterPtr TensorValueEstimater::unit_ = new UnitEstimater;

ThreadedTensorElementComputer::ThreadedTensorElementComputer(const TensorElementComputerPtr& comp)
    :
    mindepth_(0),
    fillers_(YetiRuntime::nthread_compute(), 0) //create empty vector
{
    fillers_[0] = comp;
    for (uli i=1; i < YetiRuntime::nthread_compute(); ++i)
        fillers_[i] = comp->copy();
}

ThreadedTensorElementComputer::~ThreadedTensorElementComputer()
{
}

TensorElementComputer*
ThreadedTensorElementComputer::get_computer(uli threadnum) const
{
    return fillers_[threadnum].get();
}

void
ThreadedTensorElementComputer::set_mindepth(usi depth)
{
    mindepth_ = depth;
}

usi
ThreadedTensorElementComputer::mindepth() const
{
    return mindepth_;
}

TensorElementComputer::TensorElementComputer()
    : buffer_(0), descr_(0)
{
}

TensorElementComputer::~TensorElementComputer()
{
    if (buffer_)
        ::free(buffer_);
}

TensorValueEstimater*
TensorElementComputer::get_estimater(usi depth) const
{
    return TensorValueEstimater::get_unit_estimater();
}

void
TensorElementComputer::set_index_descr(TensorIndexDescr* descr)
{
    descr_ = descr;
}

void
TensorElementComputer::compute(const uli* indices, double *data, uli n)
{
    raise(SanityCheckError, "tile filler not configured to compute type double");
}

void
TensorElementComputer::compute(const uli* indices, quad* data, uli n)
{
    raise(SanityCheckError, "tile filler not configured to compute type quad");
}

void
TensorElementComputer::compute(const uli* indices, int *data, uli n)
{
    raise(SanityCheckError, "tile filler not configured to compute type int");
}

void
TensorElementComputer::compute(const uli* indices, float *data, uli n)
{
    raise(SanityCheckError, "tile filler not configured to compute type float");
}

void
TensorElementComputer::allocate_buffer(uli maxblocksize)
{
    if (buffer_) //clear the old one
        ::free(buffer_);
    buffer_ = ::malloc(maxblocksize);
}

UnitEstimater*
TensorValueEstimater::get_unit_estimater()
{
    return unit_.get();
}

float
UnitEstimater::max_log(const uli *indices) const
{
    return 0;
}

MemsetElementComputer::MemsetElementComputer(TemplateInfo::type_t elem_type)
    : element_type_(elem_type)
{
}

TemplateInfo::type_t
MemsetElementComputer::element_type(const uli *indices, usi depth)
{
    return element_type_;
}

template <typename data_t>
void
MemsetElementComputer::memset(uli n, data_t* data)
{
    ::memset(data, 0, n * sizeof(data_t));
}

void
MemsetElementComputer::compute(const uli* indices, double* data, uli n)
{
    memset<double>(n, data);
}

void
MemsetElementComputer::compute(const uli* indices, int* data, uli n)
{
    memset<int>(n, data);
}

void
MemsetElementComputer::compute(const uli* indices, float *data, uli n)
{
    memset<float>(n, data);
}

void
MemsetElementComputer::compute(const uli* indices, quad* data, uli n)
{
    memset<quad>(n, data);
}

TensorElementComputer*
MemsetElementComputer::copy() const
{
    return new MemsetElementComputer(element_type_);
}

YetiMatrixElementComputer::YetiMatrixElementComputer(
    const MatrixPtr& m
) :
    M_(m),
    rowstart_(0),
    colstart_(0),
    nrow_(0),
    ncol_(0)
{
}

void
YetiMatrixElementComputer::compute(
    const uli* indices,
    double* data,
    uli n
)
{
    rowstart_ = descr_->get(0)->index_start(indices[0]);
    nrow_ = descr_->get(0)->nelements(indices[0]);
    colstart_ = descr_->get(1)->index_start(indices[1]);
    ncol_ = descr_->get(1)->nelements(indices[1]);

    double* dptr = data;
    for (int r=rowstart_; r < rowstart_ + nrow_; ++r)
    {
        for (int c=colstart_; c < colstart_ + ncol_; ++c, ++dptr)
        {
            (*dptr) = M_->get_element(r, c);
        }
    }
}

TensorElementComputer*
YetiMatrixElementComputer::copy() const
{
    TensorElementComputer* computer = new YetiMatrixElementComputer(M_);
    computer->set_index_descr(descr_);
    return computer;
}

DiagonalMatrixElementComputer::DiagonalMatrixElementComputer(
    const double* vals,
    uli n
)
    : values_(new double[n]),
      diagonal_matrix_estimater_(new DiagonalMatrixEstimater),
      nvalues_(n)
{
    if (vals)
        ::memcpy(values_, vals, n * sizeof(double));
    else
    {
        //identity matrix
        for (uli i=0; i < n; ++i)
            values_[i] = 1.0;
    }
}

DiagonalMatrixElementComputer::~DiagonalMatrixElementComputer()
{
    delete[] values_;
}

TensorValueEstimater*
DiagonalMatrixElementComputer::get_estimater(usi depth) const
{
    return diagonal_matrix_estimater_.get();
}

TensorElementComputer*
DiagonalMatrixElementComputer::copy() const
{
    TensorElementComputer* filler = new DiagonalMatrixElementComputer(values_, nvalues_);
    filler->set_index_descr(descr_);
    return filler;
}

float
DiagonalMatrixEstimater::max_log(const uli *indices) const
{
    return LOG_NONZERO;
    if (indices[0] == indices[1])
    {
        return LOG_NONZERO;
    }
    else
    {
        return LOG_ZERO;
    }
}

void
DiagonalMatrixElementComputer::compute(
    const uli* indices,
    double* data,
    uli n
)
{
    uli rowstart = descr_->get(0)->index_start(indices[0]);
    uli rowstop = rowstart + descr_->get(0)->nelements(indices[0]);
    uli colstart = descr_->get(1)->index_start(indices[1]);
    uli colstop = colstart + descr_->get(1)->nelements(indices[1]);

    double* dptr = data;
    for (int r=rowstart; r < rowstop; ++r)
    {
        for (int c=colstart; c < colstop; ++c, ++dptr)
        {
            if (r==c)
                (*dptr) = values_[c];
            else
                (*dptr) = 0.0;
        }
    }
}

IdentityMatrixElementComputer::IdentityMatrixElementComputer(uli nvals)
    : DiagonalMatrixElementComputer(0, nvals)
{
}

TensorElementFilter::TensorElementFilter()
    : descr_(0)
{
}

void
TensorElementFilter::set_index_descr(TensorIndexDescr* descr)
{
    descr_ = descr;
}

bool
TensorSymmetryFilter::is_zero(const uli* indices) const
{
    usi nindex = descr_->nindex();
    usi irrep = 0;
    for (usi i=0; i < nindex; ++i)
    {
        irrep = irrep ^ descr_->get(i)->irrep(indices[i]);
    }
    return irrep != 0;
}

Diagonal_IJIJ_ValueEstimater::Diagonal_IJIJ_ValueEstimater(
    double iiii,
    double ijij,
    double ijji
) :
    iiii_bound_(LOG_ZERO),
    ijij_bound_(LOG_ZERO),
    ijji_bound_(LOG_ZERO)
{
    if (iiii != 0.0) iiii_bound_ = log10(iiii);
    if (ijij != 0.0) ijij_bound_ = log10(ijij);
    if (ijji != 0.0) ijji_bound_ = log10(ijji);
}

float
Diagonal_IJIJ_ValueEstimater::max_log(const uli* indices) const
{
    uli i = indices[0];
    uli j = indices[1];
    uli k = indices[2];
    uli l = indices[3];
    if  (i == k && j==l)
    {
        if (i == j)
            return iiii_bound_;
        else
            return ijij_bound_;
    }
    else if (i == l && j==k)
    {
        return ijji_bound_;
    }
    else
    {
        return LOG_ZERO;
    }
}

Diagonal_IJIJ_ElementComputer::Diagonal_IJIJ_ElementComputer(
    double iiii,
    double ijij,
    double ijji
) :
    iiii_(iiii),
    ijij_(ijij),
    ijji_(ijji),
    estimater_(new Diagonal_IJIJ_ValueEstimater(iiii, ijij, ijji))
{
}

Diagonal_IJIJ_ElementComputer::~Diagonal_IJIJ_ElementComputer()
{
    delete estimater_;
}

void
Diagonal_IJIJ_ElementComputer::compute(const uli* indices, double* data, uli n)
{
    uli istart = descr_->get(0)->index_start(indices[0]);
    uli istop = istart + descr_->get(0)->nelements(indices[0]);
    uli jstart = descr_->get(1)->index_start(indices[1]);
    uli jstop = jstart + descr_->get(1)->nelements(indices[1]);
    uli kstart = descr_->get(2)->index_start(indices[2]);
    uli kstop = kstart + descr_->get(2)->nelements(indices[2]);
    uli lstart = descr_->get(3)->index_start(indices[3]);
    uli lstop = lstart + descr_->get(3)->nelements(indices[3]);


    double* dataptr = data;
    for (uli i=istart; i < istop; ++i)
    {
        for (uli j=jstart; j < jstop; ++j)
        {
            for (uli k=kstart; k < kstop; ++k)
            {

                for (uli l=lstart; l < lstop; ++l, ++dataptr)
                {
                    if (i == k & j==l)
                    {
                        if (i == j)
                            *dataptr = iiii_;
                        else
                            *dataptr = ijij_;
                    }
                    else if (i == l & j==k)
                    {
                        *dataptr = ijji_;
                    }
                    else
                    {
                        *dataptr = 0;
                    }
                }
            }
        }
    }
}

TensorValueEstimater*
Diagonal_IJIJ_ElementComputer::get_estimater(usi depth) const
{
    return estimater_;
}

TensorElementComputer*
Diagonal_IJIJ_ElementComputer::copy() const
{
    TensorElementComputer* filler = new Diagonal_IJIJ_ElementComputer(iiii_, ijij_, ijji_);
    filler->set_index_descr(descr_);
    return filler;
}


TensorElementComputer*
UnitElementComputer::copy() const
{
    return new UnitElementComputer;
}

void
UnitElementComputer::compute(const uli* indices, double* data, uli n)
{
    for (uli i=0; i < n; ++i) 
        data[i] = 1.0;
}

DenominatorElementComputer::DenominatorElementComputer(
    double* evals
) :
    evals_(evals)
{
}

DenominatorElementComputer::~DenominatorElementComputer()
{
}

void
DenominatorElementComputer::compute(const uli* indices, double* data, uli n)
{
    uli istart = descr_->get(0)->index_start(indices[0]);
    uli istop = istart + descr_->get(0)->nelements(indices[0]);
    uli jstart = descr_->get(1)->index_start(indices[1]);
    uli jstop = jstart + descr_->get(1)->nelements(indices[1]);
    uli kstart = descr_->get(2)->index_start(indices[2]);
    uli kstop = kstart + descr_->get(2)->nelements(indices[2]);
    uli lstart = descr_->get(3)->index_start(indices[3]);
    uli lstop = lstart + descr_->get(3)->nelements(indices[3]);


    double* dataptr = data;
    for (uli i=istart; i < istop; ++i)
    {
        double ei = evals_[i];
        for (uli j=jstart; j < jstop; ++j)
        {
            double ej = evals_[j];
            for (uli k=kstart; k < kstop; ++k)
            {
                double ek = evals_[k];
                for (uli l=lstart; l < lstop; ++l, ++dataptr)
                {
                    double el = evals_[l];
                    double eijkl = ei + ej - ek - el;
                    *dataptr = 1.0 / eijkl;
                }
            }
        }
    }
}

TensorElementComputer*
DenominatorElementComputer::copy() const
{
    TensorElementComputer* filler = new DenominatorElementComputer(evals_);
    filler->set_index_descr(descr_);
    return filler;
}

DoubleArrayElementComputer::DoubleArrayElementComputer(const double* data)
    : data_(data), dataptr_(data)
{
}

TensorElementComputer*
DoubleArrayElementComputer::copy() const
{
    return new DoubleArrayElementComputer(data_);
}

void
DoubleArrayElementComputer::compute(const uli* indices, double* data, uli n)
{
    for (uli i=0; i < n; ++i, ++data, ++dataptr_)
        (*data) = (*dataptr_);
}

