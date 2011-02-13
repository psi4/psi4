#include <map>
#include <libsmartptr/printstream.h>

#include "index.h"
#include "exception.h"
#include "env.h"
#include "class.h"
#include "permutation.h"
#include "tile.h"
#include "malloc.h"
#include "runtime.h"

using namespace yeti;
using namespace std;



DECLARE_MALLOC(IndexRangeTuple);
DECLARE_MALLOC(Indexer);

//IndexRangeTuplePtr IndexRange::scalar_tuple = 0;

IndexRange::IndexRange(
        uli start,
        const SubindexTuplePtr& tuple
) :
    index_(0),
    parent_(0),
    n_(tuple->size()), 
    start_(start),
    subranges_(tuple)
{
    //if (tuple->size() == 0)
    //    raise(SanityCheckError, "Cannot build TileIndex from tuple of size 0");
}

IndexRange::IndexRange(
        size_t start,
        size_t n
)
    : 
    index_(0),
    parent_(0),
    start_(start),
    n_(n),
    subranges_(0)
{
}

IndexRange::IndexRange(
        uli start,
        uli n,
        uli nper
)
    :
    index_(0),
    parent_(0),
    start_(start),
    n_(0), 
    subranges_(0)
{
    if (n == 0) //empty range
        return;

    vector<size_t> sizes;
    form_sizes(n, nper, sizes);
    n_ = sizes.size();
    subranges_ = new SubindexTuple(n_, 0);
    init(sizes);
}

IndexRange::IndexRange(
    size_t start,
    const std::vector<size_t>& subsizes
)
    :
    index_(0),
    parent_(0),
    start_(start),
    n_(subsizes.size()),
    subranges_(new SubindexTuple(n_, 0))
{
    init(subsizes);
}

IndexRange::IndexRange(
    size_t start,
    size_t n,
    const SubindexTuplePtr& tuple
) :
    index_(0),
    parent_(0),
    start_(start),
    n_(n),
    subranges_(tuple->slice(start, start + n))
{
}

IndexRange::IndexRange(
    uli start,
    const SubindexTuplePtr& tuple,
    uli nper
) :
    index_(0),
    parent_(0),
    start_(start),
    n_(0),
    subranges_(0)
{
    vector<size_t> sizes;
    form_sizes(tuple->size(), nper, sizes);
    n_ = sizes.size();

    subranges_ = new SubindexTuple(n_, 0);

    //loop through and create an index range object for each subrange
    vector<size_t>::const_iterator it(sizes.begin());
    uli startidx = start;
    uli idx = 0;
    for ( ; it != sizes.end(); ++it, ++idx)
    {
        size_t n(*it);
        IndexRange* idxrange = new IndexRange(startidx, n, tuple);
        idxrange->set_parent(this);
        idxrange->set_index(idx);
        subranges_->set(idx, idxrange);
        startidx += n;
    }
}

IndexRange::IndexRange(
    uli start,
    uli n,
    IndexRange* subrange
)   :
    index_(0),
    parent_(0),
    start_(start),
    n_(n), 
    subranges_(new SubindexTuple(n_, subrange))
{
}

IndexRange::IndexRange(
    usi depth,
    uli nper,
    uli start
) :
    index_(0),
    parent_(0),
    start_(start),
    n_(nper), 
    subranges_(0)
{
    if (depth == 1) //no more subranges
        return;

    subranges_ = new SubindexTuple(nper, 0);
    usi newdepth = depth - 1;
    uli startidx = start;
    for (uli i=0; i < nper; ++i, ++startidx)
    {
        IndexRange* range = new IndexRange(newdepth, nper, startidx * nper);
        subranges_->set(i, range);
    }
}

IndexRange::IndexRange(const SubindexTuplePtr &tuple)
    :
    index_(0),
    parent_(0),
    start_(0),
    n_(tuple->size()),
    subranges_(0)
{
    //validate that all of the index ranges are at the same depth
    usi maxdepth = tuple->get(0)->depth();

    SubindexTuple::iterator it(tuple->begin());
    ++it; //skip first one

    SubindexTuple::iterator stop(tuple->end());
    bool even_depth = true;
    for ( ; it != stop; ++it)
    {
        usi depth = (*it)->depth();
        if (depth != maxdepth )
        {
            even_depth = false;
            if (depth > maxdepth)
                maxdepth = depth;
        }
    }

    if (!even_depth) //rebuild the tuple
    {
        SubindexTuplePtr newtuple(new SubindexTuple(tuple->size(), 0));
        it = tuple->begin();
        uli idx = 0;
        for ( ; it != stop; ++it, ++idx)
        {
            IndexRange* range = *it;
            if (range->depth() < maxdepth)
                range = new IndexRange(range, maxdepth);

            newtuple->set(idx, range);
        }
        subranges_ = newtuple;
    }
    else
    {
        subranges_ = tuple;
    }
}


IndexRange::IndexRange(
    IndexRange *range,
    usi depth,
    uli start
)
    :
    index_(0),
    parent_(0),
    start_(start),
    n_(1),
    subranges_(new SubindexTuple(1,0))
{
    IndexRange* newrange = range;
    if (range->depth() == depth)
    {
        cerr << "cannot add layers to range of correct depth" << endl;
        abort();
    }

    while (newrange->depth() + 1 != depth)
    {
        newrange = new IndexRange(newrange, start);
    }

    subranges_->set(0, newrange);

}

IndexRange::IndexRange(
    IndexRange* range,
    uli start
)
    :
    index_(0),
    parent_(0),
    start_(start),
    n_(1),
    subranges_(new SubindexTuple(1, range))
{

}

IndexRange::~IndexRange()
{
#if DESTRUCTOR_PRINT
    class_status("destructor");
#endif
}


usi
IndexRange::depth() const
{
    if (subranges_)
        return subranges_->get(0)->depth() + 1;
    else
        return 0;
}

void
IndexRange::expand_subrange_depth(usi depth)
{
    for (usi i=0; i < subranges_->size(); ++i)
    {
        IndexRange* range = subranges_->get(i);
        range = new IndexRange(range, depth, i);
        range->set_parent(this);
        range->set_index(i);
        subranges_->set(i, range);
    }
}

bool
IndexRange::equals(IndexRange* idx) const
{
    if (!(n_ == idx->n()))
        return false;
    if (!(start_ == idx->start()))
        return false;

    //check that they both have subranges
    bool subme = subranges_;
    bool subother = idx->get_subranges();
    if (subme != subother) //incompatible subrange structure
        return false;

    if (subranges_)
    {
        SubindexTuple::iterator itme(subranges_->begin());
        SubindexTuple::iterator itother(idx->get_subranges()->begin());
        for ( ; itme != subranges_->end(); ++itme, ++itother)
        {
            bool check = (*itme)->equals(*itother);
            if (!check)
                return false;
        }
    }

    return true;
}

uli
IndexRange::finish(usi depth) const
{
    usi mydepth = this->depth();
    if (mydepth == depth)
        return start_ + n_;
    else if (mydepth < depth)
        raise(SanityCheckError, "recursion depth exceeded in start");

    //go to next level
    usi lastindex = n_ - 1;
    return subranges_->get(lastindex)->finish(depth);
}


void
IndexRange::form_sizes(
    uli n,
    uli nper,
    vector<size_t>& sizes
)
{
    if (nper > n)
        nper = n;

    size_t nblocks = n / nper;
    size_t nextra = n % nblocks; //these blocks will need one extra
    for (size_t i=0; i < nextra; ++i)
        sizes.push_back(nper + 1);
    for (size_t i=nextra; i < nblocks; ++i)
        sizes.push_back(nper);
}

IndexRange*
IndexRange::get_first_child() const
{
    return subranges_->get(0);
}


SubindexTuplePtr
IndexRange::get_subranges() const
{
    return subranges_;
}

IndexRange*
IndexRange::get_parent() const
{
    return parent_;
}

IndexRange*
IndexRange::get_subindex(uli i) const
{
    uli relative_index = i - start_;
#if YETI_SANITY_CHECK
    if (!subranges_)
    {
        //null! what are you doing?!?!?
        raise(SanityCheckError, "index range has no sub ranges for subindex");
    }

    if (relative_index >= subranges_->size())
    {
        throw LimitExceeded<size_t>(
                "index tuple size",
                subranges_->size() - 1,
                relative_index,
                __FILE__, __LINE__
              );
    }
#endif
    return subranges_->get(relative_index);
}

void
IndexRange::get_subindices(
    std::list<IndexRange*>& ranges,
    usi getdepth
)
{   usi mydepth = depth();
    if (mydepth == getdepth)
    {
        ranges.push_back(this);
    }
    else if (mydepth > getdepth)
    {
        SubindexTuple::iterator it(subranges_->begin());
        SubindexTuple::iterator stop(subranges_->end());
        for ( ; it != stop; ++it)
        {
            (*it)->get_subindices(ranges, getdepth);
        }
    }
    else
    {
        raise(SanityCheckError, "cannot retrieve at depth");
    }
}

void
IndexRange::get_subranges(std::list<IndexRange*>& sublist, usi getdepth)
{
    usi mydepth = depth();
    if (mydepth == getdepth)
    {
        sublist.push_back(this);
        return;
    }

    SubindexTuple::iterator it(subranges_->begin());
    SubindexTuple::iterator stop(subranges_->end());
    for ( ; it != stop; ++it)
    {
        (*it)->get_subranges(sublist, getdepth);
    }
}

usi
IndexRange::get_subdepth_alignment(IndexRange* range)
{
    if (this == range)
        return depth();

    usi mydepth = depth();
    usi depth = range->depth();

    if (depth > mydepth)
        raise(SanityCheckError, "range of higher depth cannot exist as subrange");

    if (depth == 0 || mydepth == 0)
        raise(SanityCheckError, "cannot align index ranges of zero depth");

    IndexRange* next = range;
    while (next->depth() > 0)
    {
        bool found = has_subrange(next);
        if (found)
            return next->depth(); //this is the alignment depth

        next = next->get_subranges()->get(0);
    }

    raise(SanityCheckError, "index ranges have no subindex alignment depth");
}

uli*
IndexRange::get_zero_set()
{
    uli* indices = yeti_malloc_indexset();
    ::memset(indices, 0, YetiRuntime::max_nindex() * sizeof(uli));
    return indices;
}

bool
IndexRange::has_subrange(IndexRange* range)
{
    usi mydepth = depth();
    usi testdepth = range->depth();
    //this cannot be a subrange because ranges are already aligned
    if (mydepth == testdepth)
        return false;

    bool correct_depth = (mydepth == testdepth + 1);

    SubindexTuple::iterator it = subranges_->begin();
    SubindexTuple::iterator stop = subranges_->end();
    for ( ; it != stop; ++it)
    {
        bool found = false;
        IndexRange* test = *it;
        if (correct_depth)
            found = (test == range);
        else
            found = range->has_subrange(range);

        if (found)
            return true;
    }
    return false;
}

void
IndexRange::increment_offsets()
{
    //set the offsets
    SubindexTuple::iterator it = subranges_->begin();
    SubindexTuple::iterator stop = subranges_->end();
    IndexRange* prev = *it;
    ++it;
    for ( ; it != stop; )
    {
        IndexRange* current = *it;
        current->increment_offsets(prev);

        //move to next
        prev = current;
        ++it;
    }
}

void
IndexRange::increment_offsets(IndexRange *range)
{
    usi maxdepth = range->depth();
    for (usi depth=0; depth <= maxdepth; ++depth)
        increment_offset(depth, range->ntot(depth) + range->start(depth));
}

void
IndexRange::increment_offset(usi depth, uli offset)
{
    if ( this->depth() == depth )
        start_ += offset;
    else //pass along to next range
    {
        SubindexTuple::iterator it(subranges_->begin());
        SubindexTuple::iterator stop(subranges_->end());
        for ( ; it != stop; ++it)
        {
            IndexRange* range = *it;
            range->increment_offset(depth, offset);
        }
    }
}


void
IndexRange::init(const vector<size_t>& subsizes)
{
    incref();
    vector<size_t>::const_iterator it(subsizes.begin());
    size_t start = start_;
    size_t idx = 0;
    for ( ; it != subsizes.end(); ++it, ++idx)
    {
        size_t n(*it);
        //create a simple index range
        IndexRange* newrange = new IndexRange(start, n);
        newrange->set_parent(this);
        newrange->set_index(idx);
        subranges_->set(idx, newrange);
        start += n;
    }
    decref();
}

bool
IndexRange::is_parent() const
{
    bool isparent = subranges_;
    return isparent;
}

uli
IndexRange::n() const
{
    return n_;
}

uli
IndexRange::nsubranges(usi depth_level) const
{
    usi dpth = depth();
    if (depth_level > dpth) //depth is too high
        raise(SanityCheckError, "requested depth is too high");

    if (depth_level == dpth)
        return n_;

    //loop subranges and return
    SubindexTuple::iterator it(subranges_->begin());
    SubindexTuple::iterator stop(subranges_->end());
    size_t ntot = 0;
    for ( ; it != stop; ++it)
    {
        IndexRange* range(*it);
        ntot += range->nsubranges(depth_level);
    }
    return ntot;
}

uli
IndexRange::ntot(usi depth) const
{
    usi mydepth = this->depth();
    if (mydepth == depth)
        return n_;
    else if (mydepth < depth)
        raise(SanityCheckError, "recursion depth exceeded in ntot");

    if (subranges_)
    {
        size_t total = 0;
        SubindexTuple::iterator it(subranges_->begin());
        SubindexTuple::iterator stop(subranges_->end());
        for ( ; it != stop; ++it)
        {
            IndexRange* range(*it);
            total += range->ntot(depth);
        }
        return total;
    }
    else
    {
        return n_;
    }
}

uli
IndexRange::nmax(usi depth) const
{
    usi mydepth = this->depth();
    if (mydepth == depth)
        return n_;
    else if (mydepth < depth)
        return 0;

    SubindexTuple::iterator it(subranges_->begin());
    SubindexTuple::iterator stop(subranges_->end());
    size_t max = 0;
    for ( ; it != stop; ++it)
    {
        size_t n = (*it)->nmax(depth);
        if (n > max)
            max = n;
    }
    return max;
}

uli
IndexRange::start(usi depth) const
{
    usi mydepth = this->depth();
    if (mydepth == depth)
        return start_;
    else if (mydepth < depth)
        raise(SanityCheckError, "recursion depth exceeded in start");

    //go to next level
    return subranges_->get(0)->start(depth);
}







void
IndexRange::print(ostream& os) const
{
    os << Env::indent << stream_printf("I(%d-%d)", start_, start_ + n_ - 1); //print inclusive with -1

    if (subranges_)
    {
        ++Env::indent;
        SubindexTuple::iterator it(subranges_->begin());
        for ( ; it != subranges_->end(); ++it)
        {
            os << endl << *it;
        }
        --Env::indent;
    }

}

void
IndexRange::set_index(uli index)
{
    index_ = index;
}

void
IndexRange::set_parent(IndexRange *range)
{
    parent_ = range;
}

void
IndexRange::set_offsets()
{
    //set the offsets
    SubindexTuple::iterator it = subranges_->begin();
    SubindexTuple::iterator stop = subranges_->end();
    IndexRange* prev = *it;
    ++it;
    for ( ; it != stop; )
    {
        IndexRange* current = *it;
        current->set_offsets(prev);

        //move to next
        prev = current;
        ++it;
    }
}

void
IndexRange::set_offsets(IndexRange* prev)
{
    start_ = prev->start() + prev->n();
}


IndexRange*
IndexRange::shift_bottom_range() const
{
    usi mydepth = depth();
    if (mydepth == 1)
    {
        uli nsubranges = subranges_->size();
        SubindexTuplePtr new_subranges = new SubindexTuple(nsubranges);
        for (uli i=0; i < nsubranges; ++i)
        {
            IndexRange* subsubrange = subranges_->get(i);
            uli nsubsubranges = 1;
            SubindexTuplePtr subsubranges = new SubindexTuple(nsubsubranges);
            subsubranges->set(0, subsubrange);
            IndexRange* subrange = new IndexRange(start_ + i, subsubranges);
            new_subranges->set(i, subrange);
        }
        IndexRange* range = new IndexRange(start_, new_subranges);
        return range;
    }
    else if (mydepth == 0)
    {
        cerr << "Should not call shift bottom range on depth 0 index range" << endl;
        abort();
    }
    else
    {
        uli nsubranges = subranges_->size();
        SubindexTuplePtr newranges = new SubindexTuple(nsubranges);
        for (uli i=0; i < nsubranges; ++i)
        {
            newranges->set(i, subranges_->get(i)->shift_bottom_range());
        }
        IndexRange* newrange = new IndexRange(start_, newranges);
        return newrange;
    }
}

void
IndexRange::sizes(std::map<uli, uli>& sizes) const
{
    if (subranges_)
    {
        SubindexTuple::iterator it(subranges_->begin());
        SubindexTuple::iterator stop(subranges_->end());
        for ( ; it != stop; ++it)
        {
            (*it)->sizes(sizes);
        }
        return;
    }
    ++sizes[n_];
}

void
IndexRange::split(uli range)
{
    subranges_ = new SubindexTuple(n_, 0);
    uli start = 0;
    for (uli i=0; i < n_; ++i, start += range)
    {
        IndexRange* newrange = new IndexRange(start, range);
        newrange->set_parent(this);
        newrange->set_index(i);
        subranges_->set(i, newrange);
    }
}


IndexRange*
IndexRange::split_bottom_range() const
{
    if (subranges_)
    {
        uli nsubranges = subranges_->size();
        SubindexTuplePtr newranges = new SubindexTuple(nsubranges);
        for (uli i=0; i < nsubranges; ++i)
        {
            newranges->set(i, subranges_->get(i)->split_bottom_range());
        }
        IndexRange* newrange = new IndexRange(start_, newranges);
        return newrange;
    }
    else
    {
        SubindexTuplePtr subranges = new SubindexTuple(n_);
        uli nper = 1;
        for (uli i=0; i < n_; ++i)
        {
            IndexRange* subrange = new IndexRange(i + start_, nper);
            subranges->set(i, subrange);
        }
        IndexRange* newrange = new IndexRange(start_, subranges);
        return newrange;
    }
}

uli
IndexRange::start() const
{
    return start_;
}

void
IndexRange::validate()
{
    if (!subranges_)
        return;

    if (subranges_->size() == 0)
        raise(SanityCheckError, "no subranges");

    if (n_ != subranges_->size())
        throw ValuesNotEqual<size_t>("no. subranges", n_, subranges_->size(), __FILE__, __LINE__);
}

IndexDescr::IndexDescr(
        const std::string& id,
        const std::string& descr,
        uli n
)
    : id_(id), descr_(descr), n_(n)
{
}

IndexDescr::~IndexDescr()
{
}

uli
IndexDescr::n() const
{
    return n_;
}

string
IndexDescr::id() const
{
    return id_;
}

string
IndexDescr::descr() const
{
    return descr_;
}

bool
IndexDescr::equals(const constIndexDescrPtr& descr) const
{
    return id_ == descr->id_;
}

Indexer::Indexer(
    const IndexRangeTuplePtr& tuple
)
 : 
    cumulsizes_(yeti_malloc_indexset()),
        nindex_(tuple->nindex()),
        nsets_(1),
        permuted_(yeti_malloc_indexset()),
        tmpindices_(yeti_malloc_indexset()),
        indexmap_(yeti_malloc_perm()),
        offsets_(yeti_malloc_indexset())
{

    //no indexmap... just use the identity
    for (usi i=0; i < nindex_; ++i)
        indexmap_[i] = i;

    init(tuple);
}

Indexer::Indexer(
    const IndexRangeTuplePtr& tuple,
    const usi* indexmap,
    usi nindex,
    const PermutationPtr& p
)
 :      
    cumulsizes_(yeti_malloc_indexset()),
        nindex_(nindex),
        nsets_(1),
        permuted_(yeti_malloc_indexset()),
        tmpindices_(yeti_malloc_indexset()),
        indexmap_(yeti_malloc_perm()),
        offsets_(yeti_malloc_indexset())
{
    if (p && !p->is_identity())
        p->image(indexmap, indexmap_, nindex);
    else
        memcpy(indexmap_, indexmap, nindex_ * sizeof(usi));


    init(tuple);
}

Indexer::Indexer(
    const size_t *sizes,
    const size_t *offsets,
    const usi* indexmap,
    usi nindex,
    const PermutationPtr& p
) :
    cumulsizes_(yeti_malloc_indexset()),
        nindex_(nindex),
        nsets_(1),
        permuted_(yeti_malloc_indexset()),
        tmpindices_(yeti_malloc_indexset()),
        indexmap_(yeti_malloc_perm()),
        offsets_(yeti_malloc_indexset())
{
    if (p && !p->is_identity())
        p->image(indexmap, indexmap_, nindex);
    else
        memcpy(indexmap_, indexmap, nindex_ * sizeof(usi));

    init(sizes, offsets, indexmap);
}

Indexer::Indexer(
    const size_t *sizes,
    const size_t *offsets,
    usi nindex,
    const PermutationPtr& p
) :
    cumulsizes_(yeti_malloc_indexset()),
        nindex_(nindex),
        nsets_(1),
        permuted_(yeti_malloc_indexset()),
        tmpindices_(yeti_malloc_indexset()),
        indexmap_(yeti_malloc_perm()),
        offsets_(yeti_malloc_indexset())
{

    usi* indexmap = yeti_malloc_perm();
    for (usi i=0; i < nindex; ++i)
        indexmap[i] = i;

    if (p && !p->is_identity())
        p->image(indexmap, indexmap_, nindex);
    else
        memcpy(indexmap_, indexmap, nindex_ * sizeof(usi));

    init(sizes, offsets, indexmap);

    yeti_free_perm(indexmap);
}

Indexer::Indexer() : //zero indexer
    cumulsizes_(0),
        nindex_(0),
        nsets_(0),
        permuted_(0),
        tmpindices_(0),
        indexmap_(0),
        offsets_(0)
{
}

Indexer::Indexer(const Indexer *parent)
    :
    cumulsizes_(yeti_malloc_indexset()),
    nindex_(parent->nindex_),
    nsets_(parent->nsets_),
    permuted_(yeti_malloc_indexset()),
    tmpindices_(yeti_malloc_indexset()),
    indexmap_(yeti_malloc_perm()),
    offsets_(yeti_malloc_indexset())
{
    for (usi i=0; i < nindex_; ++i)
    {
        cumulsizes_[i] = parent->cumulsizes_[i];
        indexmap_[i] = parent->indexmap_[i];
        offsets_[i] = parent->offsets_[i];
    }
}

Indexer*
Indexer::copy() const
{
    Indexer* cpy = new Indexer(this);
    return cpy;
}

void
Indexer::permute(const PermutationPtr &p)
{
    p->permute(cumulsizes_, tmpindices_);
    ::memcpy(cumulsizes_, tmpindices_, nindex_ * sizeof(uli));

    p->permute(offsets_, tmpindices_);
    ::memcpy(offsets_, tmpindices_, nindex_ * sizeof(uli));

    //p->permute(indexmap_, reinterpret_cast<usi*>(tmpindices_));
    //::memcpy(indexmap_, tmpindices_, nindex_ * sizeof(usi));
}

Indexer::~Indexer()
{
    if (cumulsizes_)
        yeti_free_indexset(cumulsizes_);
    if (permuted_)
        yeti_free_indexset(permuted_);
    if (tmpindices_)
        yeti_free_indexset(tmpindices_);
    if (offsets_)
        yeti_free_indexset(offsets_);
    if (indexmap_)
        yeti_free_perm(indexmap_);
}

void
Indexer::init(
    const uli* sizes,
    const uli* offsets,
    const usi* indexmap
)
{
    for (int i=nindex_ - 1; i >= 0; --i)
    {
        usi imap = indexmap_[i];
        cumulsizes_[i] = nsets_;
        nsets_ *= sizes[imap];
        offsets_[i] = offsets[imap];
    }
}

void
Indexer::init(const IndexRangeTuplePtr& tuple)
{
    IndexRangeTuple::iterator it(tuple->begin());
    for (int i=nindex_ - 1; i >= 0; --i)
    {
        usi imap = indexmap_[i];
        //create an indexer
        IndexRange* range(tuple->get(imap));
        cumulsizes_[i] = nsets_;
        nsets_ *= range->n();
        offsets_[i] = range->start();
    }
}

void
Indexer::print(ostream& os) const
{
    os << Env::indent << "Indexer" << endl;
    os << Env::indent << "Cumulative sizes" << endl;
    os << Env::indent << " ->";
    for (usi i=0; i < nindex_; ++i)
        os << " " << cumulsizes_[i];
    os << endl << Env::indent << "Offsets" << endl;
    os << Env::indent << " ->";
    for (usi i=0; i < nindex_; ++i)
        os << " " << offsets_[i];
    os << endl << Env::indent << "Index Map" << endl;
    os << Env::indent << " ->";
    for (usi i=0; i < nindex_; ++i)
        os << " " << indexmap_[i];
}

size_t
Indexer::index(const size_t* indices) const
{
    for (usi i=0; i < nindex_; ++i)
    {
        tmpindices_[i] = indices[indexmap_[i]];
    }

    size_t idx = 0;
    const size_t* indexptr = tmpindices_;
    const size_t* nptr = cumulsizes_;
    const size_t* offsetptr = offsets_;
    for (usi i=0; i < nindex_; ++i, ++indexptr, ++nptr, ++offsetptr)
    {
        idx += ((*indexptr) - (*offsetptr)) * (*nptr);
    }

    return idx;
}

usi
Indexer::nindex() const
{
    return nindex_;
}

size_t
Indexer::nsets() const
{
    return nsets_;
}

size_t
Indexer::index(const TilePtr& tile) const
{
    return index(tile->indices());
}

void
Indexer::extract(size_t index, size_t** indices) const
{
    size_t remainder = index;
    size_t* idxlist = new size_t[nindex_];
    for (usi i=0; i < nindex_; ++i)
    {
        size_t ni = remainder / cumulsizes_[i];
        idxlist[i] = ni;
        remainder -= ni * cumulsizes_[i];
    }

    if (remainder != 0)
    {
        raise(SanityCheckError, "remainder is not zero");
    }

    (*indices) = idxlist;
}

std::ostream&
yeti::operator <<(std::ostream& os, IndexRange* range)
{
    if (range)
        range->print(os);
    else
        os << "null range" << endl;
    return os;
}

std::ostream&
yeti::operator<<(std::ostream& os, IndexRangeTuple* tuple)
{
    IndexRangeTuple::iterator it(tuple->begin());
    os << "[ ";

    for ( ; it != tuple->end(); ++it)
    {
        IndexRange* range(*it);
        size_t start = range->start();
        size_t n = range->n();
        size_t end = start + n;
        os << stream_printf("I(%d-%d) ", start, end - 1);
        //-1 for inclusive
    }
    os << "]";
    return os;
}

std::ostream&
yeti::operator<<(std::ostream& os, const IndexRangeTuplePtr& tuple)
{
    os << tuple.get();
    return os;
}

IndexRangeTuple::IndexRangeTuple(usi size)
    :
    nindex_(size),
    indices_(reinterpret_cast<IndexRange**>(yeti_malloc_indexptr()) )
{
    ::memset(indices_, 0, nindex_ * sizeof(IndexRange*));
}

IndexRangeTuple::IndexRangeTuple(
    usi size,
    IndexRange** ranges
) :
    nindex_(size),
    indices_(reinterpret_cast<IndexRange**>(yeti_malloc_indexptr()) )
{
    for (usi i=0; i < size; ++i)
    {
        indices_[i] = ranges[i];
    }
}

IndexRangeTuple::IndexRangeTuple(usi size, IndexRange *tmpl)
    :
    nindex_(size),
    indices_(reinterpret_cast<IndexRange**>( yeti_malloc_indexptr() ) )
{
    for (usi i=0; i < nindex_; ++i)
    {
        indices_[i] = tmpl;
    }
}

IndexRangeTuple::~IndexRangeTuple()
{
    yeti_free_indexptr(indices_);
}

IndexRangeTuple::iterator
IndexRangeTuple::begin() const
{
    return indices_;
}

IndexRangeTuple::iterator
IndexRangeTuple::end() const
{
    return indices_ + nindex_;
}

IndexRangeTuple*
IndexRangeTuple::copy(const PermutationPtr &p) const
{
    IndexRangeTuple* tuple = new IndexRangeTuple(nindex_);
    if (p && !p->is_identity())
    {
        const usi* pmap = p->indexmap();
        for (usi i=0; i < nindex_; ++i)
        {
            tuple->set(i, indices_[pmap[i]]);
        }
    }
    else
    {
        for (usi i=0; i < nindex_; ++i)
            tuple->set(i, indices_[i]);
    }
    return tuple;
}

bool
IndexRangeTuple::is_aligned() const
{
    return maxdepth() == mindepth();
}

usi
IndexRangeTuple::mindepth() const
{
    iterator it(begin());
    iterator stop(end());
    usi min = MAX_DEPTH;
    for ( ; it != stop; ++it)
    {
        usi depth = (*it)->depth();
        if (depth < min)
            min = depth;
    }
    return min;
}

usi
IndexRangeTuple::maxdepth() const
{
    iterator it(begin());
    iterator stop(end());
    usi max = 0;
    for ( ; it != stop; ++it)
    {
        usi depth = (*it)->depth();
        if (depth > max)
            max = depth;
    }
    return max;
}

void
IndexRangeTuple::set(usi idx, IndexRange *range)
{
    if (indices_[idx] != 0)
    {
        cerr << "index already set on tuple" << endl;
        abort();
    }

    indices_[idx] = range;
    range->incref();
}

void
IndexRangeTuple::recurse_get_index_tuples(
    usi index,
    IndexRange** ranges,
    std::list<IndexRangeTuplePtr>& tuples,
    usi depth
)
{
    if (index == nindex())
    {
        IndexRangeTuplePtr newtuple
            = new IndexRangeTuple(index, ranges);
        tuples.push_back(newtuple);
        return;
    }

    std::list<IndexRange*> subindices;
    this->get(index)->get_subindices(subindices, depth);
    std::list<IndexRange*>::const_iterator it(subindices.begin());
    std::list<IndexRange*>::const_iterator stop(subindices.end());
    for ( ; it != stop; ++it)
    {
        ranges[index] = *it;
        recurse_get_index_tuples(index + 1, ranges, tuples, depth);
    }
}

void
IndexRangeTuple::get_index_tuples(
    std::list<IndexRangeTuplePtr>& tuples,
    usi depth
)
{
    IndexRange** ranges = reinterpret_cast<IndexRange**>(yeti_malloc_indexptr());
    recurse_get_index_tuples(0, ranges, tuples, depth);
    yeti_free_indexptr(ranges);
}



IndexRange*
IndexRangeTuple::get(usi idx) const
{
    if (indices_[idx] == 0)
    {
        cerr << "all indices not set on tuple" << endl;
        abort();
    }

    return indices_[idx];
}

IndexRangeTuple*
IndexRangeTuple::get_unit_range(IndexRangeTuple* subtuple)
{
    uli start = 0; uli n = 1;

    uli nindex = subtuple->nindex();

    IndexRangeTuple* tuple = new IndexRangeTuple(nindex);
    for (usi i=0; i < nindex; ++i)
    {
        IndexRange* subrange = subtuple->get(i);
        IndexRange* range = new IndexRange(start, n, subrange);
        tuple->set(i, range);
    }
    return tuple;
}

void
IndexRangeTuple::slice_front(usi nindex)
{
    usi nidx_final = nindex_ - nindex;
    for (usi i=0; i < nidx_final; ++i)
    {
        indices_[i] = indices_[i + nindex];
    }
    nindex_ = nidx_final;
}

void
IndexRangeTuple::permute(Permutation *p)
{
    if (p->is_identity())
        return;

    IndexRange** tmp = reinterpret_cast<IndexRange**>(yeti_malloc_indexset());
    p->permute<IndexRange*>(indices_, tmp);
    for (usi i=0; i < nindex_; ++i)
        indices_[i] = tmp[i];
    yeti_free_indexset(reinterpret_cast<uli*>(tmp));
}

usi
IndexRangeTuple::nindex() const
{
    return nindex_;
}

SubindexTuplePtr
SubindexTuple::slice(uli start, uli stop)
{
    uli n = stop - start;
    SubindexTuplePtr newtuple = new SubindexTuple(n, 0);
    iterator it(begin() + start);
    iterator endit(begin() + stop);

    uli idx = 0;
    for ( ; it != endit; ++it, ++idx)
        newtuple->set(idx, *it);
    return newtuple;
}

SubindexTuple::SubindexTuple(uli n, IndexRange *tmpl)
    : CountableArray<IndexRange>(n, tmpl)
{
}

SubindexTuple::SubindexTuple(uli n)
    : CountableArray<IndexRange>(n)
{
}

uli
SubindexTuple::index(IndexRange *subidx) const
{
    iterator it(begin());
    iterator stop(end());
    uli idx = 0;
    for ( ; it != stop; ++it, ++idx)
    {
        IndexRange* subrange = *it;
        if (subrange == subidx)
            return idx;
    }
    raise(SanityCheckError, "no subindex matching parameter subindex");
}

IndexRangeLocation::IndexRangeLocation(const IndexRangeTuplePtr &tuple)
    :
    n_(tuple->nindex()),
    data_(yeti_malloc_indexset())
{
    for (usi i=0; i < n_; ++i)
    {
        data_[i] = reinterpret_cast<uli>(tuple->get(i));
    }
}

void
IndexRangeLocation::print(std::ostream &os) const
{
    IndexRangeTuplePtr tuple = new IndexRangeTuple(n_);
    IndexRange** list = reinterpret_cast<IndexRange**>(data_);
    for (usi i=0; i < n_; ++i)
    {
        os << (void*) list[i] << endl;
        tuple->set(i, list[i]);
    }
    os << tuple;
}

bool
IndexRangeLocation::lt(const IndexRangeLocationPtr &r)
{
    for (usi i=0; i < n_; ++i)
    {
        if (data_[i] != r->data_[i])
            return data_[i] < r->data_[i];
    }
    return false; //exactly equals
}

IndexRangeLocation::~IndexRangeLocation()
{
    yeti_free_indexset(data_);
}

bool
IndexRangeLocationCompare::operator ()(
    const IndexRangeLocationPtr& l,
    const IndexRangeLocationPtr& r
)
{
    return l->lt(r);
}
