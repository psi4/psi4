#include "elementop.h"
#include "opimpl.h"
#include "class.h"
#include "tensor.h"
#include "tensorblock.h"
#include "exception.h"

using namespace yeti;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

ElementOp::~ElementOp()
{
}

bool
ElementOp::do_update_after() const
{
    return true;
}

void
ElementOp::configure(Tensor* tensor)
{
    //do nothing
}

void
ElementOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    double* data
)
{
    cerr << "element op not implemented for data type double" << endl;
    abort();
}

void
ElementOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    int* data
)
{
    cerr << "element op not implemented for data type int" << endl;
    abort();
}


void
ElementOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    float* data
)
{
    cerr << "element op not implemented for data type float" << endl;
    abort();
}

void
ElementOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    quad* data
)
{
    cerr << "element op not implemented for data type quad" << endl;
    abort();
}

DiaOp::DiaOp(
    const double *ei,
    const double *ea
)
    :
    evals_i_(ei),
    evals_a_(ea)
{
}

void
DiaOp::retrieve(TensorBlock* block) const
{
    block->retrieve_verbatim();
}

void
DiaOp::release(TensorBlock* block) const
{
    block->release_verbatim();
}

void
DiaOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    double* data
)
{
    uli istart_ = index_starts[0];
    uli istop_ = istart_ + sizes[0];
    uli astart_ = index_starts[1];
    uli astop_ = astart_ + sizes[1];

    double* dataptr = data;
    for (uli i=istart_; i < istop_; ++i)
    {
        double ei = evals_i_[i];
        for (uli a=astart_; a < astop_; ++a, ++dataptr)
        {
            double ea = evals_a_[a];
            double eia = ei - ea;
            (*dataptr) /= eia;
        }
    }
}

DijabOp::DijabOp(
    const double *ei,
    const double *ej,
    const double *ea,
    const double *eb,
    double ci,
    double cj,
    double ca,
    double cb
)
    :
    evals_i_(ei),
    evals_j_(ej),
    evals_a_(ea),
    evals_b_(eb),
    ci_(ci),
    cj_(cj),
    ca_(ca),
    cb_(cb)
{
}

void
DijabOp::retrieve(TensorBlock* block) const
{
    block->retrieve_verbatim();
}

void
DijabOp::release(TensorBlock* block) const
{
    block->release_verbatim();
}

void
DijabOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    double* data
)
{
    uli istart_ = index_starts[0];
    uli istop_ = istart_ + sizes[0];
    uli jstart_ = index_starts[1];
    uli jstop_ = jstart_ + sizes[1];
    uli astart_ = index_starts[2];
    uli astop_ = astart_ + sizes[2];
    uli bstart_ = index_starts[3];
    uli bstop_ = bstart_ + sizes[3];


    double* dataptr = data;
    for (uli i=istart_; i < istop_; ++i)
    {
        double ei = ci_ * evals_i_[i];
        for (uli j=jstart_; j < jstop_; ++j)
        {
            double ej = cj_ * evals_j_[j];
            for (uli a=astart_; a < astop_; ++a)
            {
                double ea = ca_ * evals_a_[a];
                for (uli b=bstart_; b < bstop_; ++b, ++dataptr)
                {
                    double eb = cb_ * evals_b_[b];
                    double eijab = ei + ej + ea  + eb;
                    (*dataptr) /= eijab;
                }
            }
        }
    }
}

ScaleOp::ScaleOp(double scale)
    : scale_(scale)
{
}

template <typename data_t>
void
ScaleOp::_scale(
    uli nblock,
    data_t* data
)
{
    data_t* dataptr = data;
    data_t scale = scale_;
    for (uli i=0; i < nblock; ++i, ++dataptr)
    {
        (*dataptr) *= scale;
    }
}


void
ScaleOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    double* data
)
{
    _scale<double>(nblock, data);
}

void
ScaleOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    int* data
)
{
    _scale<int>(nblock, data);
}

void
ScaleOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    float* data
)
{
    _scale<float>(nblock, data);
}

void
ScaleOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    quad* data
)
{
    _scale<quad>(nblock, data);
}

void
ScaleOp::retrieve(TensorBlock* block) const
{
    block->retrieve_verbatim();
}

void
ScaleOp::release(TensorBlock* block) const
{
    block->release_verbatim();
}

NormElementOp::NormElementOp()
    : normsq_(0)
{
}

void
NormElementOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    int* data
)
{
    const int* dataptr = data;
    for (uli i=0; i < nblock; ++i, ++dataptr)
    {
        int d = *dataptr;
        normsq_ += d*d;
    }
}

void
NormElementOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    double* data
)
{
    const double* dataptr = data;
    for (uli i=0; i < nblock; ++i, ++dataptr)
    {
        double d = *dataptr;
        normsq_ += d*d;
    }
}

void
NormElementOp::retrieve(TensorBlock* block) const
{
    block->retrieve_read();
}

void
NormElementOp::release(TensorBlock* block) const
{
    block->release_read();
}

double
NormElementOp::norm() const
{
    return sqrt(normsq_);
}

template <typename data_t>
void
ZeroOp::_zero(
    uli nblock,
    data_t* data
)
{
    ::memset(data, 0, nblock * sizeof(data_t));
}


void
ZeroOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    double* data
)
{
    _zero<double>(nblock, data);
}

void
ZeroOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    int* data
)
{
    _zero<int>(nblock, data);
}

void
ZeroOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    float* data
)
{
    _zero<float>(nblock, data);
}

void
ZeroOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    quad* data
)
{
    _zero<quad>(nblock, data);
}

void
ZeroOp::retrieve(TensorBlock* block) const
{
    block->retrieve_write();
}

void
ZeroOp::release(TensorBlock* block) const
{
    block->release_write();
}

bool
ZeroOp::do_update_after() const
{
    return false;
}

void
Diagonalize_IJIJ_Op::retrieve(TensorBlock* block) const
{
    block->retrieve_verbatim();
}

void
Diagonalize_IJIJ_Op::release(TensorBlock* block) const
{
    block->release_verbatim();
}

void
Diagonalize_IJIJ_Op::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    double* data
)
{
    uli istart = index_starts[0];
    uli istop = istart + sizes[0];
    uli jstart = index_starts[1];
    uli jstop = jstart + sizes[1];
    uli kstart = index_starts[2];
    uli kstop = kstart + sizes[2];
    uli lstart = index_starts[3];
    uli lstop = lstart + sizes[3];

    double* dataptr = data;
    for (uli i=istart; i < istop; ++i)
    {
        for (uli j=jstart; j < jstop; ++j)
        {
            for (uli k=kstart; k < kstop; ++k)
            {
                for (uli l=lstart; l < lstop; ++l, ++dataptr)
                {
                    if      (i==k && j==l); //do nothing
                    else if (i==l && j==k); //do nothing
                    else *dataptr = 0;
                }
            }
        }
    }
}


UpperTriangleGetOp::UpperTriangleGetOp()
    : 
    nindex_(0),
    utri_(0),
    offset_p_(0),
    offset_q_(0)
{

}

UpperTriangleGetOp::~UpperTriangleGetOp()
{
    if (utri_)
        delete[] utri_;
}

void
UpperTriangleGetOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    double* data
)
{
    uli pstart = index_starts[0] - offset_p_;
    uli pstop = pstart + sizes[0];
    uli qstart = index_starts[1] - offset_q_;
    uli qstop = qstart + sizes[1];
    double* dataptr = data;
    for (uli p=pstart; p < pstop; ++p)
    {
        for (uli q=qstart; q < qstop; ++q, ++dataptr)
        {
            if (p < q)
                continue;

            uli idx = p*(p+1)/2 + q;
            utri_[idx] = *dataptr;
        }
    }
}

void
UpperTriangleGetOp::configure(Tensor* tensor)
{
    offset_p_ = tensor->get_block_descr()->get(0)->data_index_start();
    offset_p_ = tensor->get_block_descr()->get(1)->data_index_start();
    nindex_ = tensor->get_block_descr()->get(0)->nelements_data();
    uli ntri = nindex_*(nindex_+1)/2;
    utri_ = new double[ntri];
    if (tensor->get_block_descr()->nindex() != 2)
    {
        raise(SanityCheckError, "can only get upper triangle of two-index tensor");
    }
    if (tensor->get_block_descr()->get(0) != tensor->get_block_descr()->get(1))
    {
        raise(SanityCheckError, "can only get upper triangle of symmetric matrices");
    }
    ::memset(utri_, 0, ntri * sizeof(double));
}

void
UpperTriangleGetOp::retrieve(TensorBlock* block) const
{
    block->retrieve_read();
}

void
UpperTriangleGetOp::release(TensorBlock* block) const
{
    block->release_read();
}

const double*
UpperTriangleGetOp::data() const
{
    return utri_;
}

DoubleArrayGetOp::DoubleArrayGetOp()
    : 
    data_(0),
    dataptr_(0),
    offset_p_(0),
    offset_q_(0),
    offset_r_(0),
    np_(0),
    nq_(0),
    nr_(0)
{
}

DoubleArrayGetOp::~DoubleArrayGetOp()
{
    if (data_)
        delete[] data_;
}

void
DoubleArrayGetOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    double* data
)
{
    uli pstart = index_starts[0] - offset_p_;
    uli qstart = index_starts[1] - offset_q_;
    uli rstart = index_starts[2] - offset_r_;
    uli pstop = pstart + sizes[0];
    uli qstop = qstart + sizes[1];
    uli rstop = rstart + sizes[2];
    double* dataptr = data;
    for (uli p=pstart; p < pstop; ++p)
    {
        for (uli q=qstart; q < qstop; ++q)
        {
            for (uli r=rstart; r < rstop; ++r, ++dataptr)
            {
                uli idx = (p*nq_+q)*nr_+r;
                data_[idx] = *dataptr;
            }
        }
    }
}

void
DoubleArrayGetOp::configure(Tensor* tensor)
{
    descr_ = tensor->get_block_descr();
    size_t size = tensor->get_totalsize();
    data_ = new double[size];
    dataptr_ = data_;
    ::memset(data_, 0, size * sizeof(double));
    np_ = descr_->get(0)->nelements_data();
    offset_p_ = descr_->get(0)->data_index_start();
    nq_ = descr_->get(1)->nelements_data();
    offset_q_ = descr_->get(1)->data_index_start();
    nr_ = descr_->get(2)->nelements_data();
    offset_r_ = descr_->get(2)->data_index_start();
}

void
DoubleArrayGetOp::retrieve(TensorBlock* block) const
{
    block->retrieve_read();
}
void
DoubleArrayGetOp::release(TensorBlock* block) const
{
    block->release_read();
}

void
DoubleArrayGetOp::increment_ptr(uli nelements)
{
    dataptr_ += nelements;
}

const double*
DoubleArrayGetOp::data() const
{
    return data_;
}
