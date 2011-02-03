#include "permutationimpl.h"
#include "env.h"
#include "runtime.h"

#include <libyeti/permutation.h>
#include <libyeti/exception.h>
#include <libyeti/index.h>

#include <stdarg.h>
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace yeti;
using namespace smartptr;
using namespace std;

SerialDeclare(PermutationSet);
SerialDeclare(PermutationGroup);

#define PERM_DEBUG 1

bool
less_than(const uli* a, const uli* b, usi n)
{
    for (usi i=0; i < n; ++i)
    {
        if (a[i] > b[i])
            return false;
        else if (a[i] < b[i])
            return true;
    }
    return false; //not less than, exactly equal
}

bool
permutation_less::operator() (const PermutationPtr& p, const PermutationPtr& q) const
{
    //this should not be used to sort permutations
    //when computing the quotient set it is very important that, for example
    //the permutation (3 2 0 1) comes before (3 2 1 0)
#if 0
    usi p_order = p->order();
    usi q_order = q->order();

    if (p_order != q_order)
        return p_order < q_order;
#endif

    usi p_rank = p->rank();
    usi q_rank = q->rank();

    if (p_rank != q_rank)
        return p_rank < q_rank;

    usi p_n = p->nindex();
    usi q_n = p->nindex();
    if (p_n != q_n)
        return p_n < q_n;

    short p_sign = p->sign();
    short q_sign = q->sign();
    if (p_sign != q_sign)
    {
        return p_sign < q_sign;
    }

    const usi* pptr = p->pmap_;
    const usi* qptr = q->pmap_;
    for (usi i=0; i < p_n; ++i, ++pptr, ++qptr)
    {
        usi pelem = *pptr;
        usi qelem = *qptr;
        if (pelem != qelem)
            return pelem < qelem;
    }

    return false; //exactly equal
}

bool
find_match(const usi* list, usi match, usi n)
{
    for (usi i=0; i < n; ++i)
    {
        if (list[i] == match)
            return true;
    }
    return false;
}

Permutation::Permutation(usi n, short sign, usi i, usi j)
    : nindex_(n), pmap_(yeti_malloc_perm()), sign_(sign), rank_(0), order_(0)
{
    pmap_[0] = i;
    pmap_[1] = j;

#if PERM_DEBUG
    validate();
#endif

    init();
}

Permutation::Permutation(usi n, short sign, usi i, usi j, usi k)
    : nindex_(n), pmap_(yeti_malloc_perm()), sign_(sign), rank_(0), order_(0)
{
    pmap_[0] = i;
    pmap_[1] = j;
    pmap_[2] = k;

#if PERM_DEBUG
    validate();
#endif

    init();
}

Permutation::Permutation(usi n, short sign, usi i, usi j, usi k, usi l)
    : nindex_(n), pmap_(yeti_malloc_perm()), sign_(sign), rank_(0), order_(0)
{
    pmap_[0] = i;
    pmap_[1] = j;
    pmap_[2] = k;
    pmap_[3] = l;

#if PERM_DEBUG
    validate();
#endif

    init();
}

Permutation::Permutation(usi n)
    : nindex_(n),
    pmap_(yeti_malloc_perm()),
    sign_(1),
    rank_(0),
    order_(1) //identity has order 1
{
    for (usi i=0; i < n; ++i)
        pmap_[i] = i;

#if PERM_DEBUG
    validate();
#endif

}

Permutation::Permutation(
    usi n,
    short sign,
    const usi* cycle
)
    :
    nindex_(n),
    pmap_(yeti_malloc_perm()),
    sign_(sign),
    rank_(0),
    order_(0)
{
    memcpy(pmap_, cycle, nindex_ * sizeof(usi));

#if PERM_DEBUG
    validate();
#endif

    init();
}

Permutation::Permutation(
    usi* cycle,
    usi n,
    short sign
)
    :
    nindex_(n),
    pmap_(cycle),
    sign_(sign),
    rank_(0),
    order_(0)
{

#if PERM_DEBUG
    validate();
#endif

    init();
}

void
Permutation::validate()
{
    usi* count = yeti_malloc_perm();
    memset(count, 0, sizeof(usi) * nindex_);
    for (usi i=0; i < nindex_; ++i)
    {
        if (pmap_[i] >= nindex_)
            raise(SanityCheckError, "element out of range in pmap");

        ++count[pmap_[i]];
    }

    for (usi i=0; i < nindex_; ++i)
    {
        if (count[i] != 1)
        {
            raise(SanityCheckError, "permutation is not a 1-1 mapping");
        }
    }

    yeti_free_perm(count);
}

Permutation::Permutation(
    usi n,
    const PermutationPtr& p,
    const usi* subset
) :
    nindex_(n),
    pmap_(yeti_malloc_perm()),
    sign_(p->sign_),
    rank_(p->rank_),
    order_(p->order_)
{
    for (usi i=0; i < n; ++i)
        pmap_[i] = i;

    for (usi i=0; i < p->nindex(); ++i)
        pmap_[subset[i]] = subset[p->pmap_[i]];

#if PERM_DEBUG
    validate();
#endif

}

Permutation::Permutation(
    usi n,
    usi i,
    usi j,
    short sign
) : nindex_(n), pmap_(yeti_malloc_perm()), sign_(sign), rank_(2), order_(2)
{
    for (usi idx=0; idx < n; ++idx)
    {
        pmap_[idx] = idx;
    }
    pmap_[i] = j;
    pmap_[j] = i;

#if PERM_DEBUG
    validate();
#endif

}


Permutation::Permutation(
    usi n,
    const PermutationPtr& p,
    expand_t type
) :
    nindex_(n), sign_(p->sign_), rank_(p->rank_), order_(p->order_),
    pmap_(yeti_malloc_perm())
{
    for (usi i=0; i < n; ++i)
        pmap_[i] = i;

    usi offset = 0;
    usi stop = p->nindex();
    if (type == ExpandBackward)
    {
        offset = n - p->nindex();
        stop = n;
    }

    const usi* pmapptr = p->pmap_;
    for (usi i=offset; i < stop; ++i, ++pmapptr)
        pmap_[i] = *pmapptr + offset;

#if PERM_DEBUG
    validate();
#endif

}


Permutation::~Permutation()
{
#if DESTRUCTOR_PRINT
    class_status("destructor");
#endif
    yeti_free_perm(pmap_);
}

void
Permutation::init()
{
    rank_ = compute_rank(pmap_);
    order_ = compute_order();
}

usi
Permutation::compute_order()
{
    usi* image = new usi[nindex_];
    usi* target = new usi[nindex_];
    memcpy(target, pmap_, nindex_ * sizeof(usi));
    usi order = 1;
    while (compute_rank(target) != 0)
    {
        if (order > 25)
        {
            print();
            abort();
        }
        this->permute<usi>(target, image);
        memcpy(target, image, nindex_ * sizeof(usi));
        ++order;
    }
    delete[] target;
    delete[] image;

    return order;
}

usi
Permutation::compute_rank(const usi* arr)
{
    usi rank = 0;
    for (usi i=0; i < nindex(); ++i)
    {
        if (arr[i] != i)
            ++rank;
    }
    return rank;
}

usi
Permutation::rank() const
{
    return rank_;
}

usi
Permutation::order() const
{
    return order_;
}


IndexRangeTuplePtr
Permutation::permute(const IndexRangeTuplePtr &src) const
{
    if (is_identity())
        return src; //nothing new

    make(tuple, IndexRangeTuple, src->size());
    for (usi i=0; i < nindex_; ++i)
    {
        tuple->set(i, src->get(pmap_[i]));
    }

    return tuple;
}

uli*
Permutation::permute(const uli* indices) const
{
    return permute<uli>(indices);
}

void
Permutation::permute(const uli* src, uli* dest) const
{
    permute<uli>(src, dest);
}

int*
Permutation::permute(const int* indices) const
{
    return permute<int>(indices);
}

usi*
Permutation::permute(const usi* indices) const
{
    return permute<usi>(indices);
}

void
Permutation::image(const usi *src, usi *img, usi nindex) const
{
    //change
    for (usi i=0; i < nindex; ++i)
        img[i] = pmap_[src[i]];
}

bool
Permutation::maps_to(const uli *src, const uli *dst) const
{
    for (usi i=0; i < nindex_; ++i)
    {
        if (src[pmap_[i]] != dst[i])
            return false;
    }
    return true;
}

PermutationPtr
Permutation::inverse() const
{
    if (order_ == 2)
    {
        return const_cast<Permutation*>(this);
    }

    raise(SanityCheckError, "not implemented");
}

usi
Permutation::nindex() const
{
    return nindex_;
}

short
Permutation::sign() const
{
    return sign_;
}

const usi*
Permutation::indexmap() const
{
    return pmap_;
}

PermutationPtr
Permutation::product(const constPermutationPtr& p) const
{
    usi* product_map = yeti_malloc_perm();
    permute(p->pmap_, product_map);
    short sign = sign_ * p->sign();
    make(newp, Permutation, product_map, nindex_, sign);
    return newp;
}

bool
Permutation::valid_subgroup(const usi* subset, usi n) const
{
    if (n == nindex_) //there is no subgroup. This is the full group
        return true;

    //figure out which indices are not included in the subset
    usi nomit = nindex_ - n;
    usi* omitted = new usi[nomit];
    usi idx = 0;
    for (usi i=0; i < nindex_; ++i)
    {
        if (!find_match(subset, i, n)) //this index is not in the subset
        {
            omitted[idx] = i; 
            ++idx;
        }
        if (idx == nomit)
            break;
    }

    //the permutation cannot affect any of the ommitted indices
    //if it does, this permutation is not valid for the subgroup
    for (usi i=0; i < nomit; ++i)
    {
        usi idx = omitted[i];
        //the omitted index is not mapped into itself
        if (pmap_[idx] != idx)
        {
            delete[] omitted;
            return false;
        }
    }

    delete[] omitted;
    return true;
}

PermutationPtr
Permutation::compress(const usi* subset, usi nindex) const
{
    usi* vals = yeti_malloc_perm();
    const usi* subsetptr = subset;

    //map the subset indices
    usi* subsetmap = yeti_malloc_perm();
    for (usi i=0; i < nindex; ++i)
        subsetmap[subset[i]] = i;

    for (usi i=0; i < nindex; ++i, ++subsetptr)
        vals[i] = subsetmap[pmap_[subset[i]]];
            
    make(perm, Permutation, nindex, sign_, vals);
    yeti_free_perm(subsetmap);
    yeti_free_perm(vals);

    return perm;
}

bool
Permutation::equals(const constPermutationPtr& p) const
{
    if (nindex_ != p->nindex_)
        return false;

    usi nmin = nindex_ < p->nindex_ ? nindex_ : p->nindex_;

    const usi* test = p->pmap_;
    for (usi i=0; i < nmin; ++i)
    {
        if (pmap_[i] != test[i])
            return false;
    }

    //if the permutations are not aligned, validate
    //that the permutations fix every extra element
    if (nindex_ > p->nindex_)
    {
        for (usi  i=nmin; i < nindex_; ++i)
        {
            if (pmap_[i] != i)
                return false;
        }
    }
    else if (nindex_ < p->nindex_)
    {
        for (usi  i=nmin; i < p->nindex_; ++i)
        {
            if (p->pmap_[i] != i)
                return false;
        }
    }

    return true;
}

void
Permutation::print(std::ostream& os) const
{
    os << Env::indent << "Permutation of " << nindex_ << " elements: ";
    if (sign_ == 1)
        os << "+1";
    else if (sign_ == -1)
        os << "-1";
    os << " (";
    for (usi i=0; i < nindex_; ++i)
        os << " " << pmap_[i];
    os << " )";
}

bool
Permutation::permutes(usi index) const
{
    return pmap_[index] != index;
}

bool
Permutation::fixes(usi index) const
{
    return pmap_[index] == index;
}

usi
Permutation::image(usi index) const
{
    return pmap_[index];
}

bool
Permutation::fixes_arr(const size_t* vals) const
{
    const size_t* valptr = vals;
    const usi* mapptr = pmap_;
    for (usi i=0; i < nindex_; ++i, ++valptr, ++mapptr)
    {
        if (vals[i] != vals[pmap_[i]])
            return false;
        //if (*valptr != vals[*mapptr])
        //    return false;
    }
    return true;
}

bool
Permutation::is_transposition(usi i, usi j) const
{
    for (usi idx =0; idx < nindex_; ++idx)
    {
        if (idx == i || idx == j)
            continue;

        //this is not just a transposition of i and j
        if (permutes(idx))
            return false;
    }

    //if gets permuted to j, this is that transposition
    return permutes(i);
}

bool
Permutation::lt(const PermutationPtr &p, const size_t *vals) const
{
    for (usi i=0; i < nindex_; ++i)
    {
        size_t val_me = vals[pmap_[i]];
        size_t val_p = vals[p->pmap_[i]];
        if (val_me != val_p)
            return val_me < val_p;
    }
    return false; //exactly equal
}

bool
Permutation::improves_sort(const size_t* indices) const
{
    if (rank_ == 0) //identity
        return false;

    const size_t* valptr = indices;
    usi* mapptr = pmap_;
    size_t check, val;
    for (usi i=0; i < nindex_; ++i, ++valptr, ++mapptr)
    {
        val = *valptr;
        check = indices[*mapptr];
        if (val != check)
            return check < val;
    }
    return false; //exactly equal
}



bool
Permutation::is_identity() const
{
    return rank_ == 0;
}

PermutationSet::PermutationSet(usi n)
    : nindex_(n)
{
    SetRuntime(PermutationSet);
}

PermutationSet::PermutationSet(const PermutationSetPtr &set)
    : nindex_(set->nindex()), perms_(set->perms_)
{
    SetRuntime(PermutationSet);
}

PermutationSet::PermutationSet(const XMLArchivePtr& xml)
{
    SetRuntime(PermutationSet);
}

void
PermutationSet::serialize(const XMLArchivePtr& xml) const
{
}

PermutationSet::~PermutationSet()
{
}

PermutationPtr
PermutationSet::get_lowest_permutation(const size_t* vals) const
{
    //if no permutations, just return null
    if (perms_.size() == 0)
        return 0;

    PermutationPtr p(perms_[0]);
    for (iterator it(begin() + 1); it != end(); ++it)
    {
        PermutationPtr next(*it);
        if (next->lt(p, vals))
            p = next;
    }

    return p;
}

usi
PermutationSet::nindex() const
{
    return nindex_;
}

size_t
PermutationSet::order() const
{
    return perms_.size();
}


void
PermutationSet::add(const PermutationPtr& p)
{
    if (p->nindex() != nindex_) //this permutation is not valid for this group
        throw ValuesNotEqual<usi>("permutation group order", nindex_, p->nindex(), __FILE__, __LINE__);

    if (contains(p))
        return;

    perms_.push_back(p);
}

void
PermutationSet::add_and_expand(
    const PermutationSetPtr& p,
    Permutation::expand_t type
)
{
    iterator it(p->begin());
    for ( ; it != p->end(); ++it)
    {
        PermutationPtr p(*it);
        make(expanded_p, Permutation, nindex_, p, type);
        add(expanded_p);
    }
}

void
PermutationSet::add(const PermutationSetPtr& p)
{
    iterator it(p->begin());
    for ( ; it != p->end(); ++it)
        add(*it);
}

bool
PermutationSet::contains(const PermutationPtr& p) const
{
    vector<PermutationPtr>::const_iterator it(perms_.begin());
    for ( ; it != perms_.end(); ++it)
    {
        if ( p->equals(*it) )
            return true;
    }
    return false;
}

bool
PermutationSet::equals(const PermutationSetPtr &set) const
{
    if (set->order() != order())
        return false;

    if (set->nindex() != nindex())
        return false;

    PermutationSet::iterator itme(begin());
    PermutationSet::iterator stop(end());
    PermutationSet::iterator itother(set->begin());
    for ( ; itme != stop; ++itme, ++itother)
    {
        bool eq = (*itme)->equals(*itother);
        if (!eq)
            return false;
    }
    return true;
}

template <class iterator, class element>
element
tmpl_lowest_image(iterator it, iterator stop, const uli *src, uli *dst, usi nindex)
{
    uli size = nindex * sizeof(uli);

    uli* workspace = yeti_malloc_indexset();

    element lowest = *it;
    (*it)->permute(src, dst);
    ++it; //move to second element

    for ( ; it != stop; ++it)
    {
        (*it)->permute(src, workspace);
        if ( less_than(workspace, dst, nindex) )
        {
            ::memcpy(dst, workspace, size);
            lowest = *it;
        }
    }

    yeti_free_indexset(workspace);

    return lowest;
}

Permutation*
PermutationSet::lowest_image(Permutation **perm_list, usi nperms, const uli *src, uli *dst, usi nindex)
{
    Permutation* p =
        tmpl_lowest_image<Permutation**,Permutation*>(perm_list, perm_list + nperms, src, dst, nindex);
}

Permutation*
PermutationSet::lowest_image(const uli* src, uli* dst) const
{
    PermutationPtr p =
        tmpl_lowest_image<PermutationSet::iterator,PermutationPtr>(begin(), end(), src, dst, nindex());
    return p.get();
}

void
PermutationSet::sanity_check(
    const constPermutationSetPtr& grp
) const
{
    if (nindex_ != grp->nindex())
    {
        throw ValuesNotEqual<usi>("permutation sets act on a different indices", nindex_, grp->nindex(), __FILE__, __LINE__);
    }
}

PermutationSetPtr
PermutationSet::intersection_set(const constPermutationSetPtr& grp) const
{
    sanity_check(grp);

    make(intergrp, PermutationSet, nindex_);
    vector<PermutationPtr>::const_iterator it(perms_.begin());
    for ( ; it != perms_.end(); ++it)
    {
        PermutationPtr p(*it);
        if (grp->contains(p))
            intergrp->add(p);
    }
    return intergrp;
}

bool
PermutationSet::improves_sort(const size_t *indices) const
{
    PermutationSet::iterator it(begin());
    PermutationSet::iterator stop(end());
    for ( ; it != stop; ++it)
    {
        if ( (*it)->improves_sort(indices) )
            return true;
    }
    return false; //found none
}

PermutationSetPtr
PermutationSet::union_set(const constPermutationSetPtr& grp) const
{
    sanity_check(grp);

    make(uniongrp, PermutationSet, nindex_);

    vector<PermutationPtr>::const_iterator it(perms_.begin());
    for ( ; it != perms_.end(); ++it) //add all the elements
        uniongrp->add(*it);

    it = grp->perms_.begin();
    for ( ; it != grp->perms_.end(); ++it)
    {
        PermutationPtr p(*it);
        if (!contains(p)) //add any elements we missed
            uniongrp->add(p);
    }
    return uniongrp;
}

PermutationSetPtr
PermutationSet::compressed_subset(const usi* set, usi n) const
{
    make(subset, PermutationSet, n);
    vector<PermutationPtr>::const_iterator it(perms_.begin());
    for ( ; it != perms_.end(); ++it)
    {
        PermutationPtr p(*it);
        if (p->valid_subgroup(set, n))
            subset->add(p->compress(set, n));
    }
    return subset;
}


PermutationSetPtr
PermutationSet::subset(const usi* set, usi n) const
{
    make(subset, PermutationSet, nindex_);
    vector<PermutationPtr>::const_iterator it(perms_.begin());
    for ( ; it != perms_.end(); ++it)
    {
        PermutationPtr p(*it);
        if (p->valid_subgroup(set, n))
            subset->add(p);
    }

    return subset;
}

PermutationSetPtr
PermutationSet::isotropy_set(const size_t* vals) const
{   
    PermutationSet* me = const_cast<PermutationSet*>(this);
    if (perms_.size() == 0)
        return me; //just send me back

    make(subset, PermutationSet, nindex_);
    vector<PermutationPtr>::const_iterator it(perms_.begin());
    for ( ; it != perms_.end(); ++it)
    {
        PermutationPtr p(*it);
        if (p->fixes_arr(vals))
            subset->add(p);
    }

    //to conserve memory, if the subgrp is the same, just return me
    if (subset->order() == (usi) perms_.size())
        return me;

    return subset;
}

PermutationSetPtr
PermutationSet::orbit(const PermutationPtr &p) const
{
    PermutationSet::iterator it(begin());
    PermutationSet::iterator stop(end());
    PermutationSetPtr newset(new PermutationSet(nindex_));
    for ( ; it != stop; ++it)
    {
        PermutationPtr next(*it);
        PermutationPtr product = p->product(next);
        newset->add(product);
    }
    return newset;
}

PermutationSetPtr
PermutationSet::isotropy_set(const IndexSetPtr& idxset) const
{
    return isotropy_set(idxset->data());
}

bool
PermutationSet::contains(const PermutationSetPtr& p) const
{
    PermutationSet::iterator it(p->begin());
    for ( ; it != p->end(); ++it)
    {
        if (!contains(*it))
            return false;
    }
    return true;
}

void
PermutationSet::clear()
{
    perms_.clear();
    make(p, Permutation, nindex());
    perms_.push_back(p);
}

PermutationSet::iterator
PermutationSet::begin() const
{
    return perms_.begin();
}

PermutationSet::iterator
PermutationSet::end() const
{
    return perms_.end();
}

bool
PermutationSet::has_transposition(usi i, usi j) const
{
    iterator it(begin());
    for ( ; it != end(); ++it)
    {
        if ( (*it)->is_transposition(i, j) )
            return true;
    }
    return false;
}

void
PermutationSet::print(ostream& os) const
{
    os << Env::indent << this->rtinfo_->classname() << " of Order " << perms_.size();
    vector<PermutationPtr>::const_iterator it(perms_.begin());
    ++Env::indent;
    for (; it != perms_.end(); ++it)
    {
        os << endl << (*it);
        //(*it)->print(os);
    }
    --Env::indent;
}

void
PermutationSet::sort()
{
    permutation_less p_lt_q;
    ::sort(perms_.begin(), perms_.end(), p_lt_q);
}

PermutationGroup::PermutationGroup(usi n)
    : PermutationSet(n), closed_(true) //only identity, so closed
{
    SetRuntime(PermutationGroup);
    //make the identity permutation
    make(p, Permutation, n);
    perms_.push_back(p);
}

PermutationGroup::PermutationGroup(const PermutationSetPtr& set)
    : PermutationSet(set), closed_(false)
{
    SetRuntime(PermutationGroup);
}

PermutationGroup::PermutationGroup(const XMLArchivePtr& xml)
    : PermutationSet(xml)
{
    SetRuntime(PermutationGroup);
}

PermutationGroup::PermutationGroup(
    usi n,
    const PermutationGroupPtr& grp,
    const usi* subset
) : PermutationSet(n), closed_(false)
{
    SetRuntime(PermutationGroup);
    for ( iterator it(grp->begin()); it != grp->end(); ++it )
    {
        make(p, Permutation, n, *it, subset);
        add(p);
    }
    closed_ = true;
}

void
PermutationGroup::serialize(const XMLArchivePtr& xml) const
{
}

PermutationGroup::~PermutationGroup()
{
}

void
PermutationGroup::sanity_check(
    const constPermutationGroupPtr& grp
) const
{
    PermutationSet::sanity_check(grp);
    if (!grp->closed())
        raise(SanityCheckError, "parameter group is not closed in method call");

    sanity_check();
}

void
PermutationGroup::sanity_check() const
{
    if (!closed_)
        raise(SanityCheckError, "cannot call method on non-closed group");
}

PermutationGroupPtr
PermutationGroup::conjugate(const PermutationPtr &p) const
{
    if (p->order() != 2)
        raise(SanityCheckError, "can only conjugate by an element of order 2");

    PermutationGroupPtr newgroup(new PermutationGroup(this->nindex_));
    PermutationGroup::iterator it(begin());
    PermutationGroup::iterator stop(end());
    for ( ; it != stop; ++it)
    {
        PermutationPtr next(*it);
        PermutationPtr conj = (p->product(next))->product(p);
        newgroup->add(conj);
    }
    newgroup->closed_ = true;
    return newgroup;
}

void
PermutationGroup::close()
{
    if (closed_)
        return;

    unsigned int check = 0;
    while (check != order())
    {
        check = order();
        int size = perms_.size();
        for (int i=0; i < size; ++i)
        {
            PermutationPtr p(perms_[i]);
            for (int j=0; j < size; ++j)
            {
                PermutationPtr q(perms_[j]);
                PermutationPtr product(p->product(q));
                add(product);
            }
        }
    }
    sort();
    closed_ = true;
}

void
PermutationGroup::add(const PermutationPtr& p)
{
    usi size = perms_.size();
    PermutationSet::add(p);

    if (perms_.size() != size) //we added one
        closed_ = false;
}


PermutationGroupPtr
PermutationGroup::intersection_grp(const constPermutationGroupPtr& grp) const
{
    PermutationSetPtr set(this->intersection_set(grp));
    make(pgrp, PermutationGroup, set);
    pgrp->closed_ = true;
    return pgrp;
}

PermutationGroupPtr
PermutationGroup::union_grp(const constPermutationGroupPtr& grp) const
{
    PermutationSetPtr set(this->union_set(grp));
    make(pgrp, PermutationGroup, set);
    pgrp->close();
    return pgrp;
}

PermutationPtr
PermutationGroup::get_identity() const
{
    for ( iterator it(begin()); it != end(); ++it )
    {
        PermutationPtr p(*it);
        if (p->is_identity())
            return p;
    }

    raise(SanityCheckError, "how does the permutation group not contain the identity?");
}

PermutationGroupPtr
PermutationGroup::compressed_subgrp(const usi* subset, usi n) const
{
    PermutationSetPtr set(this->compressed_subset(subset, n));
    make(grp, PermutationGroup, set);
    grp->closed_ = true;
    return grp;
}

PermutationGroupPtr
PermutationGroup::subgrp(const usi* subset, usi n) const
{
    PermutationSetPtr set(this->subset(subset, n));
    make(grp, PermutationGroup, set);
    grp->closed_ = true;
    return grp;
}

bool
PermutationGroup::improves_sort(const size_t *indices) const
{
    if (nindex() == 1)
        return false; //only identity element

    return PermutationSet::improves_sort(indices);
}

PermutationGroupPtr
PermutationGroup::isotropy_grp(const size_t* vals) const
{
    PermutationSetPtr set(this->isotropy_set(vals));
    make(grp, PermutationGroup, set);
    grp->closed_ = true;
    return grp;
}

PermutationGroupPtr
PermutationGroup::isotropy_grp(const IndexSetPtr& vals) const
{
    PermutationSetPtr set(this->isotropy_set(vals));
    make(grp, PermutationGroup, set);
    grp->closed_ = true;
    return grp;
}

bool
PermutationGroup::closed() const
{
    return closed_;
}

PermutationSetPtr
PermutationGroup::quotient_set(const PermutationGroupPtr &pgrp) const
{
    sanity_check(pgrp);

    //validate that the passed in group is actually a subgrp
    if (!contains(pgrp))
        raise(SanityCheckError, "parameter subgroup for quotient is not a subgroup");

    make(qgrp, PermutationSet, nindex());
    PermutationSet::iterator it_numer(begin());

    for ( ; it_numer != end(); ++it_numer)
    {
        PermutationPtr e(*it_numer);

        bool has_match = false;
        PermutationSet::iterator it_denom(pgrp->begin());
        for ( ; it_denom != pgrp->end(); ++it_denom)
        {
            PermutationPtr f(*it_denom);
            PermutationPtr product(f->product(e));
            if (qgrp->contains(product))
            {
                has_match = true;
                break;
            }
        }

        if (!has_match) //not mapped to anything else
            qgrp->add(e);
    }
    return qgrp;
}

void
PermutationGroup::add(const PermutationSetPtr &set)
{
    PermutationSet::add(set);
}
