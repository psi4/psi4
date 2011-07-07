#include "permutation.h"
#include "sort.h"
#include "env.h"
#include "index.h"
#include "exception.h"
#include "runtime.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

using namespace yeti;
using namespace std;

DECLARE_MALLOC(Sort);

Sort::Sort(Permutation* p)
    :
    p_(p),
    inv_p_(p->inverse()),
    ntot_(1),
    nindex_(p->nindex())
{
    //do not configure
}

void
Sort::configure(Permutation* p)
{
    p_ = p;
}

void
Sort::configure(const uli* sizes)
{
    usi nidxperm = p_->nindex();

    /** JJW 12/15/2010 */
    nindex_ = nidxperm; //reset nindex in case of previous configures

    //compute the lengths
    size_t src_stride = sizes[nidxperm - 1];
    nstrides_[nidxperm - 1] = 1;
    for (int i=nidxperm-2; i >=0 ; --i) //start on element n-1 and loop backward
    {
        nstrides_[i] = src_stride;
        src_stride *= sizes[i];
    }
    ntot_ = src_stride;
    //now do the permutations
    p_->permute(nstrides_, lengths_);
    p_->permute(sizes, nstrides_);

    return;

    bool p_fixes_last = p_->fixes(nidxperm - 1);

    if (p_fixes_last)
    {
        for (int idx = nidxperm - 2; idx >= 0 && p_->fixes(idx); --idx)
        {
            lengths_[idx] = 1;
            nstrides_[idx] *= nstrides_[idx + 1];
            --nindex_;
        }
    }

#if YETI_SANITY_CHECK
    if (ntot_ * sizeof(double) > SORT_BUFFER_SIZE)
    {
        cerr << "sort buffer is not large enough to handle everything!" << endl;
        abort();
    }
#endif
}

Sort::~Sort()
{
}

char*
Sort::data_buffer() const
{
    return const_cast<char*>(data_buffer_);
}

char*
Sort::metadata_buffer(usi depth) const
{
    size_t increment = SORT_BUFFER_SIZE * depth / (YetiRuntime::max_depth() + 1);
    return const_cast<char*>(metadata_buffer_ + increment);
}

const uli*
Sort::lengths() const
{
    return lengths_;
}

const uli*
Sort::nstrides() const
{
    return nstrides_;
}

usi
Sort::nindex() const
{
    return nindex_;
}

uli
Sort::ntot() const
{
    return ntot_;
}

uli
Sort::nelements() const
{
    uli n = 1;
    for (usi i=0; i < nindex_; ++i)
        n *= nstrides_[i];
    return n;
}

uli
Sort::total_stride() const
{
    uli n = 0;
    for (usi i=0; i < nindex_; ++i)
        n += (nstrides_[i] - 1) * lengths_[i];
    return n;
}

void
Sort::print(std::ostream &os) const
{
    os << "Sort for " << p_ << endl;
    os << "Stride Lengths = " << ClassOutput<const uli*>::str(nindex_, lengths_) << endl;
    os << "No. Strides    = " << ClassOutput<const uli*>::str(nindex_, nstrides_) << endl;
}

Permutation*
Sort::get_permutation() const
{
    return p_;
}

Permutation*
Sort::get_inverse_permutation() const
{
    return inv_p_;
}
