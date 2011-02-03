#include "permutation.h"
#include "sort.h"
#include "env.h"
#include "index.h"
#include "exception.h"
#include "runtime.h"


using namespace yeti;
using namespace std;

DECLARE_MALLOC(Sort);

Sort::Sort(
    const PermutationPtr &p,
    const IndexRangeTuplePtr &tuple
)
    : p_(p),
        lengths_(yeti_malloc_indexset()),
        nstrides_(yeti_malloc_indexset()),
        ntot_(1),
        nindex_(p->nindex())
{
    if (p->nindex() != tuple->size())
        raise(SanityCheckError, "tuple and permutation have different number of indices");

    configure(tuple);
}


Sort::Sort(
    const PermutationPtr& p,
    const size_t* sizes
)
    : p_(p),
        lengths_(yeti_malloc_indexset()),
        nstrides_(yeti_malloc_indexset()),
        ntot_(1),
        nindex_(p->nindex())
{
    configure(sizes);
}

void
Sort::configure(
    const PermutationPtr &p,
    const IndexRangeTuplePtr &tuple
)
{
    p_ = p;
    configure(tuple);
}

void
Sort::configure(const IndexRangeTuplePtr &tuple)
{
    size_t* sizes = yeti_malloc_indexset();
    usi i=0;
    IndexRangeTuple::iterator it(tuple->begin());
    IndexRangeTuple::iterator stop(tuple->end());
    for ( ; it != stop; ++it, ++i)
    {
        sizes[i] = (*it)->n();
    }
    configure(sizes);
    yeti_free_indexset(sizes);
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
}

Sort::~Sort()
{
    yeti_free_indexset(nstrides_);
    yeti_free_indexset(lengths_);
}

usi
Sort::n() const
{
    return nindex_;
}

uli
Sort::ntot() const
{
    return ntot_;
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
    return p_.get();
}
