#include <map>
#include <libsmartptr/printstream.h>

#include "index.h"
#include "exception.h"
#include "env.h"
#include "class.h"
#include "permutationimpl.h"
#include "malloc.h"
#include "runtime.h"

using namespace yeti;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif


//IndexRangeTuplePtr IndexRange::scalar_tuple = 0;

IndexRange::IndexRange(
        uli start,
        const SubindexTuplePtr& tuple
) :
    index_(0),
    parent_(0),
    nelements_(tuple->size()),
    start_(start),
    subranges_(tuple),
    has_symmetry_(false),
    irrep_(0)
{
    //if (tuple->size() == 0)
    //    yeti_throw(SanityCheckError, "Cannot build TileIndex from tuple of size 0");
}

IndexRange::IndexRange(
        uli start,
        uli n
)
    : 
    index_(0),
    parent_(0),
    start_(start),
    nelements_(n),
    subranges_(0),
    has_symmetry_(false),
    irrep_(0)
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
    nelements_(0),
    subranges_(0),
    has_symmetry_(false),
    irrep_(0)
{
    if (n == 0) //empty range
        return;

    vector<uli> sizes;
    form_sizes(n, nper, sizes);
    nelements_ = sizes.size();
    subranges_ = new SubindexTuple(nelements_);
    init(sizes);
}

IndexRange::IndexRange(
    uli start,
    const std::vector<uli>& subsizes
)
    :
    index_(0),
    parent_(0),
    start_(start),
    nelements_(subsizes.size()),
    subranges_(new SubindexTuple(nelements_)),
    has_symmetry_(false),
    irrep_(0)
{
    init(subsizes);
}

IndexRange::IndexRange(
    uli start,
    uli n,
    const SubindexTuplePtr& tuple
) :
    index_(0),
    parent_(0),
    start_(start),
    nelements_(n),
    subranges_(tuple->slice(start, start + n)),
    has_symmetry_(false),
    irrep_(0)
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
    nelements_(0),
    subranges_(0),
    has_symmetry_(false),
    irrep_(0)
{
    vector<uli> sizes;
    form_sizes(tuple->size(), nper, sizes);
    nelements_ = sizes.size();

    subranges_ = new SubindexTuple(nelements_);

    //loop through and create an index range object for each subrange
    vector<uli>::const_iterator it(sizes.begin());
    uli startidx = start;
    uli idx = 0;
    for ( ; it != sizes.end(); ++it, ++idx)
    {
        uli n(*it);
        IndexRange* idxrange = new IndexRange(startidx, n, tuple);
        idxrange->set_parent(this);
        idxrange->set_index(idx);
        subranges_->set(idx, idxrange);
        startidx += n;
    }
}

IndexRange::IndexRange(
    usi depth,
    uli nper,
    uli start
) :
    index_(0),
    parent_(0),
    start_(start),
    nelements_(nper),
    subranges_(0),
    has_symmetry_(false),
    irrep_(0)
{
    if (depth == 1) //no more subranges
        return;

    subranges_ = new SubindexTuple(nper);
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
    nelements_(tuple->size()),
    subranges_(0),
    has_symmetry_(false),
    irrep_(0)
{
    //validate that all of the index ranges are at the same depth
    usi maxdepth = tuple->get(0)->depth();
    bool even_depth = true;
    SubindexTuple::iterator it(tuple->begin());
    SubindexTuple::iterator stop(tuple->end());
    for ( ; it != stop; ++it)
    {
        IndexRange* range = *it;
        usi depth = range->depth();
        if (depth != maxdepth )
        {
            even_depth = false;
            if (depth > maxdepth)
                maxdepth = depth;
        }
    }

    if (!even_depth) //rebuild the tuple
    {
        cerr << "all subranges must have even depth" << endl;
        abort();
    }

    subranges_ = tuple;
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
    nelements_(1),
    subranges_(new SubindexTuple(1)),
    has_symmetry_(false),
    irrep_(0)
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
    nelements_(1),
    subranges_(new SubindexTuple(1)),
    has_symmetry_(false),
    irrep_(0)
{
    subranges_->set(0, range);
}

IndexRange::~IndexRange()
{
#if DESTRUCTOR_PRINT
    class_status("destructor");
#endif
}

void
IndexRange::acquire_subranges(const SubindexTuplePtr& subtuple)
{
    if (!subranges_)
    {
        yeti_throw(SanityCheckError, "no subranges in acquire_subranges");
    }

    uli start = 0;
    uli compidx = 0;
    for (uli i=0; i < subranges_->size(); ++i)
    {
        IndexRange* subrange = subranges_->get(i);
        SubindexTuplePtr new_subtuple = new SubindexTuple(subrange->nelements());
        for (uli j=0; j < subrange->nelements(); ++j, ++compidx)
        {
            IndexRange* subsubrange = subtuple->get(compidx);
            new_subtuple->set(j,subsubrange);
        }
        IndexRange* new_subrange = new IndexRange(start, new_subtuple);
        subranges_->set(i, new_subrange);
        start += new_subrange->nelements();
    }
}

bool
IndexRange::is_contiguous() const
{
    if (!subranges_)
        return true;

    SubindexTuple::iterator it(subranges_->begin());
    SubindexTuple::iterator stop(subranges_->end());
    uli start = subranges_->get(0)->start();
    for ( ; it != stop; ++it)
    {
        IndexRange* range = *it;
        if (range->start() != start)
        {
            return false;
        }
        start += range->nelements();
    }
    return true;
}

bool
IndexRange::contains(IndexRange* range) const
{
    if (range->depth() < depth())
    {
        return this->has_subrange(range);
    }

    SubindexTuple::iterator it = range->get_subranges()->begin();
    SubindexTuple::iterator stop = range->get_subranges()->end();
    for ( ; it != stop; ++it)
    {
        IndexRange* subrange = *it;
        bool check = this->contains(subrange);
        if (!check)
            return false;
    }

    //everything checks out
    return true;
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
    if (!(nelements_ == idx->nelements()))
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
        return start_ + nelements_;
    else if (mydepth < depth)
        yeti_throw(SanityCheckError, "recursion depth exceeded in start");

    //go to next level
    usi lastindex = nelements_ - 1;
    return subranges_->get(lastindex)->finish(depth);
}


void
IndexRange::form_sizes(
    uli n,
    uli nper,
    vector<uli>& sizes
)
{
    if (nper > n)
        nper = n;

    uli nblocks = n / nper;
    if (n % nper) //extra
        ++nblocks;

    nper = n / nblocks;
    uli nextra = n - nper * nblocks;

    for (uli i=0; i < nextra; ++i)
        sizes.push_back(nper + 1);
    for (uli i=nextra; i < nblocks; ++i)
        sizes.push_back(nper);
}

IndexRange*
IndexRange::get_composite_range(usi depth) const
{
    usi mydepth = this->depth();
    if (depth > mydepth)
        yeti_throw(SanityCheckError, "depth too high for composite range");

    if (depth == mydepth)
    {
        IndexRange* range = new IndexRange(subranges_);
        range->set_irrep(irrep_);
        return range;
    }

    std::list<IndexRange*> sublist;
    get_subranges(sublist, depth);

    SubindexTuplePtr tuple = new SubindexTuple(sublist.size());
    std::list<IndexRange*>::const_iterator it(sublist.begin());
    std::list<IndexRange*>::const_iterator stop(sublist.end());
    uli idx = 0;
    usi irrep = (*it)->get_irrep();
    bool uniform_irrep = true;
    for ( ; it != stop; ++it, ++idx)
    {
        IndexRange* range = *it;
        if (range->get_irrep() != irrep)
            uniform_irrep = false;

        tuple->set(idx, range);
    }

    uli start = this->start(depth + 1);


    IndexRange* range = new IndexRange(start, tuple);
    if (uniform_irrep)
            range->set_irrep(irrep);
    return range;
}

IndexRange*
IndexRange::get_first_child() const
{
    return subranges_->get(0);
}

usi
IndexRange::get_irrep() const
{
    return irrep_;
}

SubindexTuple*
IndexRange::get_subranges() const
{
    return subranges_.get();
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
        yeti_throw(SanityCheckError, "index range has no sub ranges for subindex");
    }

    if (relative_index >= subranges_->size())
    {
        throw LimitExceeded<uli>(
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
        yeti_throw(SanityCheckError, "cannot retrieve at depth");
    }
}

void
IndexRange::get_subranges(std::list<IndexRange*>& sublist, usi getdepth) const
{
    usi mydepth = depth();
    if (mydepth == getdepth)
    {
        sublist.push_back(const_cast<IndexRange*>(this));
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
        yeti_throw(SanityCheckError, "range of higher depth cannot exist as subrange");

    if (depth == 0 || mydepth == 0)
        yeti_throw(SanityCheckError, "cannot align index ranges of zero depth");

    IndexRange* next = range;
    while (next->depth() > 0)
    {
        bool found = has_subrange(next);
        if (found)
            return next->depth(); //this is the alignment depth

        next = next->get_subranges()->get(0);
    }

    yeti_throw(SanityCheckError, "index ranges have no subindex alignment depth");
    return 0;
}

IndexRange*
IndexRange::get_merged_range(IndexRange *r1, IndexRange *r2)
{
    uli nsubranges = r1->nelements() + r2->nelements();
    SubindexTuplePtr merged_tuple = new SubindexTuple(nsubranges);
    uli nr1 = r1->nelements();
    for (uli i=0; i < nr1; ++i)
    {
        merged_tuple->set(i, r1->get_subranges()->get(i));
    }

    uli nr2 = r2->nelements();
    for (uli i=0; i < nr2; ++i)
    {
        merged_tuple->set(i + nr1, r2->get_subranges()->get(i));
    }

    IndexRange* merged = new IndexRange(r1->start(), merged_tuple);
    return merged;
}

bool
IndexRange::has_subrange(IndexRange* range) const
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
IndexRange::init(const vector<uli>& subsizes)
{
    incref();
    vector<uli>::const_iterator it(subsizes.begin());
    uli start = start_;
    uli idx = 0;
    for ( ; it != subsizes.end(); ++it, ++idx)
    {
        uli n(*it);
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
IndexRange::nelements() const
{
    return nelements_;
}

uli
IndexRange::nsubranges(usi depth_level) const
{
    usi mydepth = depth();
    if (depth_level > mydepth) //depth is too high
        yeti_throw(SanityCheckError, "requested depth is too high");

    if (depth_level == mydepth)
        return nelements_;

    //loop subranges and return
    SubindexTuple::iterator it(subranges_->begin());
    SubindexTuple::iterator stop(subranges_->end());
    uli ntot = 0;
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
        return nelements_;
    else if (mydepth < depth)
        yeti_throw(SanityCheckError, "recursion depth exceeded in ntot");

    if (subranges_)
    {
        uli total = 0;
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
        return nelements_;
    }
}

uli
IndexRange::nmax(usi depth) const
{
    usi mydepth = this->depth();
    if (mydepth == depth)
        return nelements_;
    else if (mydepth < depth)
        return 0;

    SubindexTuple::iterator it(subranges_->begin());
    SubindexTuple::iterator stop(subranges_->end());
    uli max = 0;
    for ( ; it != stop; ++it)
    {
        uli n = (*it)->nmax(depth);
        if (n > max)
            max = n;
    }
    return max;
}

void
IndexRange::offset(const IndexRangePtr& range)
{
    usi mydepth = depth();
    for (usi i=0; i < mydepth; ++i)
    {
        uli offset = range->finish(i);
        increment_offset(i, offset);
    }
}

uli
IndexRange::start(usi depth) const
{
    usi mydepth = this->depth();
    if (mydepth == depth)
        return start_;
    else if (mydepth < depth)
        yeti_throw(SanityCheckError, "recursion depth exceeded in start");

    //go to next level
    return subranges_->get(0)->start(depth);
}

void
IndexRange::print(ostream& os) const
{
    os << Env::indent;
    if (nelements_)
        os << stream_printf("I(%lu-%lu)", start_, start_ + nelements_ - 1); //print inclusive with -1
    else
        os << "empty index range";

    if (has_symmetry_)
        os << " symmetry: true";
    else
        os << " symmetry: false";

    os << " irrep: " << irrep_;

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

bool
IndexRange::has_symmetry() const
{
    return has_symmetry_;
}

void
IndexRange::set_irrep(usi irrep)
{
    irrep_ = irrep;
    has_symmetry_ = true;

    if (subranges_)
    {
        SubindexTuple::iterator it = subranges_->begin();
        SubindexTuple::iterator stop = subranges_->end();
        for ( ; it != stop; ++it)
        {
            IndexRange* subrange = *it;
            subrange->set_irrep(irrep);
        }
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
IndexRange::set_offset(uli offset)
{
    start_ = offset;
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
    start_ = prev->start() + prev->nelements();
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
    ++sizes[nelements_];
}

void
IndexRange::split(uli range)
{
    subranges_ = new SubindexTuple(nelements_);
    uli start = 0;
    for (uli i=0; i < nelements_; ++i, start += range)
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
        SubindexTuplePtr subranges = new SubindexTuple(nelements_);
        uli nper = 1;
        for (uli i=0; i < nelements_; ++i)
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

uli
IndexRange::stop() const
{
    return start_ + nelements_;
}

void
IndexRange::subrange_index_location(IndexRange* range, uli* indices) const
{
    bool check = _subrange_index_location(range, indices);

    if (!check)
    {
        yeti_throw(SanityCheckError, "index range is not a subrange");
    }
}

bool
IndexRange::_subrange_index_location(IndexRange* range, uli* indexset) const
{
    usi mydepth = depth();
    usi testdepth = range->depth();

    //this cannot be a subrange because ranges are already aligned
    if (mydepth == testdepth)
    {
        return this == range;
    }

    SubindexTuple::iterator it = subranges_->begin();
    SubindexTuple::iterator stop = subranges_->end();
    uli idx = 0;
    for ( ; it != stop; ++it, ++idx)
    {
        IndexRange* subrange = *it;
        bool found = subrange->_subrange_index_location(range, indexset);
        if (found)
        {
            indexset[mydepth] = idx + start_;
            return true;
        }
    }
    return false;
}

void
IndexRange::validate()
{
    if (!subranges_)
        return;

    if (subranges_->size() == 0)
        yeti_throw(SanityCheckError, "no subranges");

    if (nelements_ != subranges_->size())
        throw ValuesNotEqual<uli>("no. subranges", nelements_, subranges_->size(), __FILE__, __LINE__);
}

IndexDescr::IndexDescr(
    const std::string& descr,
    const IndexRangePtr& range
)
    :
    descr_(descr),
    range_list_(new SubindexTuple(range->depth() + 1)),
    depth_(range->depth()),
    nelements_per_index_(0),
    data_start_indices_(0),
    data_offset_(0),
    nelements_total_(0),
    parent_range_(range),
    index_to_irrep_(0),
    top_index_offset_(range->start()),
    total_data_sizes_(new size_t[range->nelements()]),
    ntop_index_(range->nelements()),
    max_nelements_data_(0),
    average_nelements_data_(0),
    total_data_size_(0),
    total_metadata_size_(0),
    max_nelements_metadata_(0),
    average_nelements_metadata_(0),
    nranges_data_(0),
    subdescr_(0),
    descr_id_(DEFAULT_INDEX_DESCR_ID)
{


    if (!range->is_contiguous())
    {
        cerr << "Cannot build index descr from non-contiguous index range" << endl;
        range->print(cout); cout << endl;
        abort();
    }

    for (usi d=1; d < depth_; ++d)
    {
        IndexRange* composite_range = range->get_composite_range(d-1);
        range_list_->set(d, composite_range);
    }
    range_list_->set(depth_, range.get());

    IndexRange* lowest_metadata_range = this->get_range(1);
    data_offset_ = lowest_metadata_range->start();
    uli nranges = lowest_metadata_range->nelements();
    nelements_per_index_ = new uli[nranges];
    data_start_indices_ = new uli[nranges];
    for (uli i=0; i < nranges; ++i)
    {
        IndexRange* range = lowest_metadata_range
                                ->get_subindex(i + data_offset_);
        uli n = range->nelements();
        nelements_per_index_[i] = n;
        data_start_indices_[i] = range->start();
        nelements_total_ += n;

        if (n > max_nelements_data_)
            max_nelements_data_ = n;

        total_data_size_ += n;
    }

    uli nelems = range->nelements();
    index_to_irrep_ = new usi[nelems];
    for (uli i=0; i < nelems; ++i)
    {
        index_to_irrep_[i] = range->get_subranges()->get(i)->get_irrep();
    }


    usi data_depth = 1;
    usi metadata_depth = 2;
    total_metadata_size_ = 0;
    for (uli i=0; i < ntop_index_; ++i)
    {
        total_data_sizes_[i] = range->get_subranges()->get(i)->ntot();

        uli nelements_metadata = range->get_subranges()->get(i)->ntot(data_depth);
        if (nelements_metadata > max_nelements_metadata_)
            max_nelements_metadata_ = nelements_metadata;
        total_metadata_size_ += nelements_metadata;
    }

    average_nelements_data_ = total_data_size_ / nelements(data_depth);
    average_nelements_metadata_ = total_metadata_size_ / nelements(metadata_depth);
    nranges_data_ = lowest_metadata_range->nelements();

    if (depth_ > 2)
    {
        IndexRange* subrange = range_list_->get(depth_ - 1);
        std::string subdescr_str = "subdescr " + descr_;
        subdescr_ = new IndexDescr(subdescr_str, subrange);
    }
}

IndexDescr::~IndexDescr()
{
    delete[] nelements_per_index_;
    delete[] data_start_indices_;
    delete[] index_to_irrep_;
    delete[] total_data_sizes_;
    range_list_ = 0;
}

IndexDescr*
IndexDescr::get_subdescr() const
{
    return subdescr_.get();
}

uli
IndexDescr::get_descr_id() const
{
    return descr_id_;
}

IndexRange*
IndexDescr::get_top_range() const
{
    return range_list_->get(depth_);
}

IndexRange*
IndexDescr::get_parent_range() const
{
    return parent_range_.get();
}

bool
IndexDescr::has_symmetry() const
{
    return parent_range_->has_symmetry();
}

bool
IndexDescr::is_equivalent(IndexDescr* descr) const
{
    if (this == descr)
        return true;
    else if (descr_id_ == descr->get_descr_id() && descr_id_ != DEFAULT_INDEX_DESCR_ID)
        return true;

    //nope!
    return false;
}

usi
IndexDescr::depth() const
{
    return depth_;
}

usi
IndexDescr::irrep(uli idx) const
{
    return index_to_irrep_[idx - top_index_offset_];
}

void
IndexDescr::print(std::ostream &os) const
{
    os << descr_;
    for (usi i=1; i <= depth_; ++i)
    {
        os << endl << "Depth=" << i << endl;
        os << range_list_->get(i);
    }
    os << endl << "data offset: " << this->data_offset_;
    os << endl << "total data size: " << this->total_data_size_;
    os << endl << "total metadata size: " << this->total_metadata_size_;
}

void
IndexDescr::set_descr_id(uli id)
{
    descr_id_ = id;
}

uli
IndexDescr::start(usi depth) const
{
    return range_list_->get(depth)
            ->get_subranges()->get(0)->start();
}

uli
IndexDescr::index_start(uli index) const
{
    return data_start_indices_[index - data_offset_];
}

uli
IndexDescr::index_increment(uli index) const
{
    return data_start_indices_[index - data_offset_] 
            - data_start_indices_[0];
}

uli
IndexDescr::data_index_start() const
{
    return data_start_indices_[0];
}

uli
IndexDescr::max_nelements_data() const
{
    return max_nelements_data_;
}

uli
IndexDescr::average_nelements_data() const
{
    return average_nelements_data_;
}

uli
IndexDescr::max_nelements_metadata() const
{
    return max_nelements_metadata_;
}

uli
IndexDescr::average_nelements_metadata() const
{
    return average_nelements_metadata_;
}

uli
IndexDescr::nelements(usi depth, uli index) const
{
    return get_range(depth)->get_subindex(index)->nelements();
}

uli
IndexDescr::nelements(usi depth) const
{
    return get_range(depth)->nelements();
}

uli
IndexDescr::nelements(uli index) const
{
    return nelements_per_index_[index - data_offset_];
}

uli
IndexDescr::nelements_total() const
{
    return nelements_total_;
}

uli
IndexDescr::nelements_data() const
{
    return nelements_total_;
}

string
IndexDescr::descr() const
{
    return descr_;
}

IndexRange*
IndexDescr::get_range(usi depth, uli index) const
{
    return range_list_->get(depth)->get_subindex(index);
}

IndexRange*
IndexDescr::get_range(usi depth) const
{
    return range_list_->get(depth);
}

size_t
IndexDescr::total_data_size(uli idx) const
{
    return total_data_sizes_[idx - top_index_offset_];
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

SubindexTuplePtr
SubindexTuple::slice(uli start, uli stop)
{
    uli n = stop - start;
    SubindexTuplePtr newtuple = new SubindexTuple(n);
    iterator it(begin() + start);
    iterator endit(begin() + stop);

    uli idx = 0;
    for ( ; it != endit; ++it, ++idx)
        newtuple->set(idx, *it);
    return newtuple;
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
    yeti_throw(SanityCheckError, "no subindex matching parameter subindex");
    return 0;
}

TensorIndexDescr::TensorIndexDescr(usi nindex)
    :
    indices_(nindex),
    nindex_(nindex),
    has_symmetry_(true)
{
}

TensorIndexDescr*
TensorIndexDescr::copy() const
{
    TensorIndexDescr* cpy = new TensorIndexDescr(nindex_);
    for (usi i=0; i < nindex_; ++i)
        cpy->set(i, indices_.get(i));
    return cpy;
}

usi
TensorIndexDescr::depth() const
{
    return indices_.get(0)->depth();
}

IndexDescr*
TensorIndexDescr::get(usi index) const
{
    return indices_.get(index);
}

bool
TensorIndexDescr::has_symmetry() const
{
    return has_symmetry_;
}

void
TensorIndexDescr::set(
    usi index,
    const IndexDescrPtr& descr
)
{
    has_symmetry_ = has_symmetry_ && descr->has_symmetry();
    indices_.set(index, descr.get());
}

size_t
TensorIndexDescr::max_nelements_data() const
{
    size_t size = 1;
    for (usi i=0; i < nindex_; ++i)
    {
        size *= indices_.get(i)->max_nelements_data();
    }
    return size;
}

size_t
TensorIndexDescr::nblocks_tot_at_depth(usi depth) const
{
    size_t size = 1;
    for (usi i=0; i < nindex_; ++i)
    {
        size *= indices_.get(i)->nelements(depth);
    }
    return size;
}

size_t
TensorIndexDescr::average_nelements_data() const
{
    size_t size = 1;
    for (usi i=0; i < nindex_; ++i)
    {
        size *= indices_.get(i)->average_nelements_data();
    }
    return size;
}

size_t
TensorIndexDescr::max_nblocks_data() const
{
    size_t nblocks = 1;
    for (usi i=0; i < nindex_; ++i)
    {
        nblocks *= indices_.get(i)->max_nelements_metadata();
    }
    return nblocks;
}

size_t
TensorIndexDescr::average_nelements_metadata() const
{
    size_t size = 1;
    for (usi i=0; i < nindex_; ++i)
    {
        size *= indices_.get(i)->average_nelements_metadata();
    }
    return size;
}

size_t
TensorIndexDescr::max_nelements_metadata() const
{
    size_t nblocks = 1;
    for (usi i=0; i < nindex_; ++i)
    {
        nblocks *= indices_.get(i)->max_nelements_metadata();
    }
    return nblocks;
}

usi
TensorIndexDescr::nindex() const
{
    return nindex_;
}

void
TensorIndexDescr::permute(Permutation* p)
{
    IndexDescr* tmp[NINDEX];
    p->permute(indices_.begin(), tmp);
    ::memcpy(indices_.begin(), tmp, nindex_ * sizeof(IndexDescr*));
}

uli
TensorIndexDescr::get_nelements_data(
    usi depth,
    const uli *indices
) const
{
    uli ntot = 1;
    for (usi i=0; i < nindex_; ++i)
    {
        ntot *= get(i)->get_range(depth)->get_subindex(indices[i])->ntot();
    }
    return ntot;
}

uli
TensorIndexDescr::get_nelements_metadata(
    usi depth,
    const uli *indices
) const
{
    uli ntot = 1;
    for (usi i=0; i < nindex_; ++i)
    {
        ntot *= get(i)->get_range(depth)->get_subindex(indices[i])->ntot(1); //metadata depth of 1
    }
    return ntot;
}

size_t
TensorIndexDescr::size(
    const uli *indices,
    usi depth
)
{
    size_t ntot = 1;
    for (usi i=0; i < nindex_; ++i)
    {
        ntot *= get(i)->get_range(depth)
            ->get_subindex(indices[i])->nelements();
    }
    return ntot;
}

void
TensorIndexDescr::get_nelements(
    usi depth,
    const uli* indexset,
    uli* nelements
) const
{
     for (usi i=0; i < nindex_; ++i)
    {
        nelements[i] = get(i)->nelements(depth, indexset[i]);
    }
}

void
TensorIndexDescr::get_nelements(const uli* indexset, uli* nelements) const
{
    for (usi i=0; i < nindex_; ++i)
    {
        nelements[i] = get(i)->nelements(indexset[i]);
    }
}

void
TensorIndexDescr::get_nelements(const uli* indexset, uli* nelements, Permutation* p) const
{
    const usi* indexmap = p->indexmap();
    for (usi i=0; i < nindex_; ++i)
    {
        usi imap = indexmap[i];
        nelements[i] = get(imap)->nelements(indexset[i]);
    }
}


void
TensorIndexDescr::get_index_starts(const uli* indexset, uli* starts) const
{
     for (usi i=0; i < nindex_; ++i)
    {
        starts[i] = get(i)->index_start(indexset[i]);
    }
}

void
TensorIndexDescr::get_index_increments(const uli* indexset, uli* increments) const
{
     for (usi i=0; i < nindex_; ++i)
    {
        increments[i] = get(i)->index_increment(indexset[i]);
    }
}

void
TensorIndexDescr::get_index_starts(uli* starts) const
{
     for (usi i=0; i < nindex_; ++i)
    {
        starts[i] = get(i)->data_index_start();
    }
}

uli
TensorIndexDescr::nelements(const uli* indices)
{
    uli ntot = 1;
    for (usi i=0; i < nindex_; ++i)
    {
        ntot *= get(i)->nelements(indices[i]);
    }
    return ntot;
}


void
TensorIndexDescr::get_nelements_data(uli* nelements) const
{
    for (usi i=0; i < nindex_; ++i)
    {
        nelements[i] = get(i)->nelements_data();
    }
}

void
TensorIndexDescr::print(std::ostream &os) const
{
    os << "Tensor Index Descr " << endl;

    uli ntot = 1;
    for (uli i=0; i < nindex_; ++i)
        ntot *= get(i)->get_top_range()->nelements();

    os << "N=" << ntot << endl;
    for (usi i=0; i < nindex_; ++i)
    {
        os << "Index " << i << ": " << get(i)->descr() << endl;
    }
}

TensorIndexDescr*
TensorIndexDescr::get_subdescr() const
{
    TensorIndexDescr* newdescr = new TensorIndexDescr(nindex_);
    for (usi i=0; i < nindex_; ++i)
    {
        newdescr->set(i, get(i)->get_subdescr());
    }
    return newdescr;
}

ulli
TensorIndexDescr::totalsize() const
{
    ulli size = 1;
    for (usi i=0; i < nindex_; ++i)
        size *= indices_.get(i)->nelements_total();
    return size;
}

size_t
TensorIndexDescr::total_data_size(const uli* indices) const
{
    size_t ntot = 1;
    for (usi i=0; i < nindex_; ++i)
        ntot *= indices_.get(i)->total_data_size(indices[i]);
    return ntot;
}

