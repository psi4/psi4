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

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

#define PERM_DEBUG 1

static uli factorial_[NINDEX_PLUS_ONE];

uli Permutation::factorial_[NINDEX_PLUS_ONE];
Permutation* Permutation::plus_perms_[NINDEX_PLUS_ONE][NPERMUTATIONS];
Permutation* Permutation::minus_perms_[NINDEX_PLUS_ONE][NPERMUTATIONS];
Permutation* Permutation::identities_[NINDEX_PLUS_ONE];
uli Permutation::inverses_[NINDEX_PLUS_ONE][NPERMUTATIONS];
uli Permutation::products_[NINDEX_PLUS_ONE][NPERMUTATIONS][NPERMUTATIONS];


uli
factorial(usi n)
{
    uli idx=1;
    for (usi i=2; i <= n; ++i)
    {
        idx *= i;
    }
    return idx;
}

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
permutation_less::operator() (Permutation* p, Permutation* q) const
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

    const usi* pptr = p->pmap_;
    const usi* qptr = q->pmap_;
    for (usi i=0; i < p_n; ++i, ++pptr, ++qptr)
    {
        usi pelem = *pptr;
        usi qelem = *qptr;
        if (pelem != qelem)
            return pelem < qelem;
    }

    short p_sign = p->sign();
    short q_sign = q->sign();
    if (p_sign != q_sign)
    {
        return p_sign < q_sign;
    }

    return false; //exactly equal
}

struct permutation_rank_less 
{
    bool operator() (Permutation* p, Permutation* q) const
    {
        return p->rank() < q->rank();   
    }
};

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

Permutation::Permutation(
    usi n,
    short sign,
    const usi* cycle
)
    :
    nindex_(n),
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

void
Permutation::validate()
{
    usi count[NINDEX];
    memset(count, 0, sizeof(usi) * nindex_);
    for (usi i=0; i < nindex_; ++i)
    {
        if (pmap_[i] >= nindex_)
            yeti_throw(SanityCheckError, "element out of range in pmap");

        ++count[pmap_[i]];
    }

    for (usi i=0; i < nindex_; ++i)
    {
        if (count[i] != 1)
        {
            yeti_throw(SanityCheckError, "permutation is not a 1-1 mapping");
        }
    }

}

Permutation::~Permutation()
{
}

void
Permutation::init()
{
    rank_ = compute_rank(pmap_);
    order_ = compute_order();

    id_index_ = id_index(nindex_, pmap_);
}

void
Permutation::init_permutations()
{
    for (usi i=0; i <= NINDEX; ++i)
        factorial_[i] = factorial(i);

    for (usi i=1; i <= NINDEX; ++i)
    {
        build_permutations(i);
    }

    //do "zero" permutations
    usi null_indices[NINDEX];
    short plus = 1;
    usi zero = 0;
    Permutation* zeroperm = new (Permutation::ValidPermutationMalloc)
            Permutation(zero, plus, null_indices);
    identities_[0] = zeroperm;
    plus_perms_[0][0] = zeroperm;
    minus_perms_[0][0] = zeroperm;
    inverses_[0][0] = 0;

    products_[0][0][0] = 0;
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

bool
Permutation::is_inverse(Permutation* p) const
{
    uli pidx = p->id_index();
    uli invidx = inverses_[nindex_][id_index_];
    return pidx == invidx;
}

Permutation*
Permutation::inverse() const
{
    uli idx = inverses_[nindex_][id_index_];
    return get_permutation(sign_, nindex_, idx);

    usi pmap[NINDEX];
    for (usi i=0; i < nindex_; ++i)
    {
        pmap[pmap_[i]] = i;
    }
    Permutation* p = new Permutation(nindex_, sign_, pmap);
    return p;
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

uli
Permutation::id_index() const
{
    return id_index_;
}

uli
Permutation::id_index(usi nindex, const usi* indexmap)
{
    usi effective_indices[NINDEX];
    usi last_index = 0;
    usi next_index = 0;
    ::memcpy(effective_indices, indexmap, nindex * sizeof(usi));
    for (usi i=0; i < nindex; ++i)
    {
        last_index = effective_indices[i];
        for (usi j=i+1; j < nindex; ++j)
        {
            next_index = effective_indices[j];
            if (last_index < next_index)
                effective_indices[j] = next_index - 1;
        }

    }

    uli idx = 0;
    for (usi i=0; i < nindex-1; ++i)
    {
        idx += effective_indices[i] * factorial_[nindex - i - 1];
    }

    return idx;
}

const usi*
Permutation::indexmap() const
{
    return pmap_;
}

Permutation*
Permutation::product(Permutation* p) const
{
    uli idx = products_[nindex_][id_index_][p->id_index_];
    short sign = sign_ * p->sign_;
    return Permutation::get_permutation(sign, nindex_, idx);

    usi product_map[NINDEX];
    permute(p->pmap_, product_map);

    return get_permutation(sign_ * p->sign_, nindex_, product_map);

    Permutation* newp = new Permutation(nindex_, sign, product_map);
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

Permutation*
Permutation::compress(const usi* subset, usi nindex) const
{
    usi vals[NINDEX];
    const usi* subsetptr = subset;

    //map the subset indices
    usi subsetmap[NINDEX];
    for (usi i=0; i < nindex; ++i)
        subsetmap[subset[i]] = i;

    for (usi i=0; i < nindex; ++i, ++subsetptr)
        vals[i] = subsetmap[pmap_[subset[i]]];

    return Permutation::get_permutation(sign_, nindex, vals);
}

bool
Permutation::equals(Permutation* p) const
{
#if YETI_SANITY_CHECK
    if (nindex_ != p->nindex_)
    {
        cerr << "permutations with different numbers of indices compared" << endl;
        abort();
    }
#endif

    return p->id_index_ == id_index_;

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
    os << " " << "[" << id_index_ << "]";
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
Permutation::fixes_arr(const uli* vals) const
{
    return fixes_arr<uli>(vals);
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
Permutation::lt(Permutation* p, const uli *vals) const
{
    for (usi i=0; i < nindex_; ++i)
    {
        uli val_me = vals[pmap_[i]];
        uli val_p = vals[p->pmap_[i]];
        if (val_me != val_p)
            return val_me < val_p;
    }
    return false; //exactly equal
}

bool
Permutation::improves_sort(const uli* indices) const
{
    if (rank_ == 0) //identity
        return false;

    const uli* valptr = indices;
    const usi* mapptr = pmap_;
    uli check, val;
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
    : nindex_(n),
      closed_(true)
{
}

PermutationSet::PermutationSet(const PermutationSetPtr &set)
    : nindex_(set->nindex()),
      perms_(set->perms_),
      closed_(set->closed_)
{
}

PermutationSet::~PermutationSet()
{
}

bool
PermutationSet::is_product(Permutation* p) const
{
    PermutationSet::iterator itq(this->begin());
    PermutationSet::iterator stop(this->end());
    for ( ; itq != stop; ++itq)
    {
        Permutation* q = *itq;
        if (q->is_identity())
            continue;
            
        PermutationSet::iterator itr(this->begin());
        for ( ; itr != stop; ++itr)
        {
            Permutation* r = *itr;
            if (r->is_identity())
                continue;
            
            Permutation* qr = q->product(r);
            if (qr->equals(p))
                return true;
        }
    }
    return false;
}

Permutation*
PermutationSet::get(uli index) const
{
    return perms_[index];
}

PermutationSetPtr
PermutationSet::get_generator_set() const
{
    PermutationSetPtr genset = new PermutationSet(nindex_);
    PermutationGroupPtr prodgrp = new PermutationGroup(nindex_);

    const_cast<PermutationSet*>(this)->sort();

    PermutationSet::iterator it(begin());
    PermutationSet::iterator stop(end());
    

    for ( ; it != stop; ++it)
    {
        Permutation* p = *it;
        if (p->is_identity())
            continue;

        if (!prodgrp->is_product(p))
            genset->add(p);
        
        if (!prodgrp->contains(p))
        {
            prodgrp->add(p);
            prodgrp->close();
        }
       
    }
    
    return genset;
}

Permutation*
PermutationSet::get_lowest_permutation(const uli* vals) const
{
    //if no permutations, just return null
    if (perms_.size() == 0)
        return 0;

    Permutation* p(perms_[0]);
    for (iterator it(begin() + 1); it != end(); ++it)
    {
        Permutation* next(*it);
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

uli
PermutationSet::order() const
{
    return perms_.size();
}


void
PermutationSet::add(Permutation* p)
{
    if (p->nindex() != nindex_) //this permutation is not valid for this group
        throw ValuesNotEqual<usi>("permutation group order", nindex_, p->nindex(), __FILE__, __LINE__);

    if (contains(p))
        return;

    perms_.push_back(p);
    closed_ = false;
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
        Permutation* p(*it);
        add(Permutation::get_permutation(p->sign(), p->nindex(), nindex_, p->indexmap(), type));
    }
}

void
PermutationSet::add(const PermutationSetPtr& p)
{
    iterator it(p->begin());
    for ( ; it != p->end(); ++it)
        add(*it);
}

PermutationSetPtr
PermutationSet::complement_set(const PermutationSetPtr& set) const
{
    PermutationSetPtr compset = new PermutationSet(nindex_);
    PermutationSet::iterator it(set->begin());
    PermutationSet::iterator stop(set->end());
    for ( ; it != stop; ++it)
    {
        Permutation* p = *it;
        if (!this->contains(p))
            compset->add(p);
    }
    return compset;
}

bool
PermutationSet::contains(Permutation* p) const
{
    vector<Permutation*>::const_iterator it(perms_.begin());
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

    uli workspace[NINDEX];

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
    return lowest;
}

Permutation*
PermutationSet::lowest_image(Permutation **perm_list, usi nperms, const uli *src, uli *dst, usi nindex)
{
    Permutation* p =
        tmpl_lowest_image<Permutation**,Permutation*>(perm_list, perm_list + nperms, src, dst, nindex);
    return p;
}

Permutation*
PermutationSet::lowest_image(const uli* src, uli* dst) const
{
    Permutation* p =
        tmpl_lowest_image<PermutationSet::iterator,Permutation*>(begin(), end(), src, dst, nindex());
    return p;
}

void
PermutationSet::sanity_check(
    const PermutationSetPtr& grp
) const
{
    if (nindex_ != grp->nindex())
    {
        throw ValuesNotEqual<usi>("permutation sets act on a different indices", nindex_, grp->nindex(), __FILE__, __LINE__);
    }
}

PermutationSetPtr
PermutationSet::intersection_set(const PermutationSetPtr& grp) const
{
    sanity_check(grp);

    make(intergrp, PermutationSet, nindex_);
    vector<Permutation*>::const_iterator it(perms_.begin());
    for ( ; it != perms_.end(); ++it)
    {
        Permutation* p(*it);
        if (grp->contains(p))
            intergrp->add(p);
    }
    return intergrp;
}

bool
PermutationSet::improves_sort(const uli* indices) const
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
PermutationSet::union_set(const PermutationSetPtr& grp) const
{
    sanity_check(grp);

    PermutationSet* uniongrp = new PermutationSet(nindex_);

    vector<Permutation*>::const_iterator it(perms_.begin());
    for ( ; it != perms_.end(); ++it) //add all the elements
        uniongrp->add(*it);

    it = grp->perms_.begin();
    for ( ; it != grp->perms_.end(); ++it)
    {
        Permutation* p(*it);
        if (!contains(p)) //add any elements we missed
            uniongrp->add(p);
    }
    return uniongrp;
}

PermutationSetPtr
PermutationSet::compressed_subset(const usi* set, usi n) const
{
    PermutationSet* subset = new PermutationSet(n);
    vector<Permutation*>::const_iterator it(perms_.begin());
    for ( ; it != perms_.end(); ++it)
    {
        Permutation* p(*it);
        if (p->valid_subgroup(set, n))
            subset->add(p->compress(set, n));
    }
    return subset;
}


PermutationSetPtr
PermutationSet::subset(const usi* set, usi n) const
{
    make(subset, PermutationSet, nindex_);
    vector<Permutation*>::const_iterator it(perms_.begin());
    for ( ; it != perms_.end(); ++it)
    {
        Permutation* p(*it);
        if (p->valid_subgroup(set, n))
            subset->add(p);
    }

    return subset;
}

PermutationSetPtr
PermutationSet::isotropy_set(const uli* vals) const
{   
    PermutationSet* me = const_cast<PermutationSet*>(this);
    if (perms_.size() == 0)
        return me; //just send me back

    make(subset, PermutationSet, nindex_);
    vector<Permutation*>::const_iterator it(perms_.begin());
    for ( ; it != perms_.end(); ++it)
    {
        Permutation* p(*it);
        if (p->fixes_arr(vals))
            subset->add(p);
    }

    //to conserve memory, if the subgrp is the same, just return me
    if (subset->order() == (usi) perms_.size())
        return me;

    return subset;
}

PermutationSetPtr
PermutationSet::orbit(Permutation* p) const
{
    PermutationSet::iterator it(begin());
    PermutationSet::iterator stop(end());
    PermutationSetPtr newset(new PermutationSet(nindex_));
    for ( ; it != stop; ++it)
    {
        Permutation* next(*it);
        Permutation* product = p->product(next);
        newset->add(product);
    }
    return newset;
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
}

void
PermutationSet::close()
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
            Permutation* p(perms_[i]);
            for (int j=0; j < size; ++j)
            {
                Permutation* q(perms_[j]);
                Permutation* product(p->product(q));
                add(product);
            }
        }
    }
    sort();
    closed_ = true;
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
    os << Env::indent << "Permutation Set " << " of Order " << perms_.size();
    vector<Permutation*>::const_iterator it(perms_.begin());
    ++Env::indent;
    for (; it != perms_.end(); ++it)
    {
        os << endl;
        (*it)->print(os);
    }
    --Env::indent;
}

void
PermutationSet::sort()
{
    permutation_less p_lt_q;
    ::sort(perms_.begin(), perms_.end(), p_lt_q);
}

bool
PermutationSet::nontrivial_isotropy(const uli *indices) const
{
    PermutationSet::iterator it(begin());
    PermutationSet::iterator stop(end());
    for ( ; it != stop; ++it)
    {
        Permutation* p(*it);
        if (!p->is_identity() && p->fixes_arr(indices))
            return true;
    }
    return false;
}


PermutationGroup::PermutationGroup(usi n)
    : PermutationSet(n)
{
    //make the identity permutation
    Permutation* identity = Permutation::get_identity(n);
    perms_.push_back(identity);
    closed_ = true; //only identity, so closed
}

PermutationGroup::PermutationGroup(const PermutationSetPtr& set)
    : PermutationSet(set)
{
    closed_ = false;
}

PermutationGroup::PermutationGroup(
    usi nindex,
    const PermutationGroupPtr& grp,
    const usi* subset
) : PermutationSet(nindex)
{
    for ( iterator it(grp->begin()); it != grp->end(); ++it )
    {
        Permutation* psub = *it;
        Permutation* p = Permutation::get_permutation(
                                    psub->sign(),
                                    psub->nindex(),
                                    nindex,
                                    psub->indexmap(),
                                    subset
                             );
        add(p);
    }
    closed_ = true;
}

PermutationGroup::~PermutationGroup()
{
}

void
PermutationGroup::sanity_check(
    const PermutationGroupPtr& grp
) const
{
    PermutationSet::sanity_check(grp);
    if (!grp->closed())
        yeti_throw(SanityCheckError, "parameter group is not closed in method call");

    sanity_check();
}

void
PermutationGroup::sanity_check() const
{
    if (!closed_)
        yeti_throw(SanityCheckError, "cannot call method on non-closed group");
}

PermutationGroupPtr
PermutationGroup::conjugate(Permutation* p) const
{
    Permutation* pinv = p->inverse();

    PermutationGroupPtr newgroup(new PermutationGroup(this->nindex_));
    PermutationGroup::iterator it(begin());
    PermutationGroup::iterator stop(end());
    for ( ; it != stop; ++it)
    {
        Permutation* next(*it);
        Permutation* conj = (p->product(next))->product(pinv);
        newgroup->add(conj);
    }
    newgroup->closed_ = true;
    return newgroup;
}

PermutationGroupPtr
PermutationGroup::intersection_grp(const PermutationGroupPtr& grp) const
{
    PermutationSetPtr set(this->intersection_set(grp));
    make(pgrp, PermutationGroup, set);
    pgrp->closed_ = true;
    return pgrp;
}

PermutationGroupPtr
PermutationGroup::union_grp(const PermutationSetPtr& grp) const
{
    PermutationSetPtr set(this->union_set(grp));
    PermutationGroup* pgrp = new PermutationGroup(set);
    pgrp->close();
    return pgrp;
}

Permutation*
PermutationGroup::get_identity() const
{
    for ( iterator it(begin()); it != end(); ++it )
    {
        Permutation* p(*it);
        if (p->is_identity())
            return p;
    }

    yeti_throw(SanityCheckError, "how does the permutation group not contain the identity?");
    return 0;
}

PermutationGroupPtr
PermutationGroup::compressed_subgrp(const usi* subset, usi n) const
{
    PermutationSetPtr set(this->compressed_subset(subset, n));
    PermutationGroup* grp = new PermutationGroup(set);
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
PermutationGroup::improves_sort(const uli* indices) const
{
    if (nindex() == 1)
        return false; //only identity element

    return PermutationSet::improves_sort(indices);
}

PermutationGroupPtr
PermutationGroup::isotropy_grp(const uli* vals) const
{
    PermutationSetPtr set(this->isotropy_set(vals));
    make(grp, PermutationGroup, set);
    grp->closed_ = true;
    return grp;
}

bool
PermutationSet::closed() const
{
    return closed_;
}

void
PermutationSet::set_closed()
{
    closed_ = true;
}

PermutationSetPtr
PermutationGroup::quotient_set(const PermutationGroupPtr& pgrp) const
{
    sanity_check(pgrp);

    //validate that the passed in group is actually a subgrp
    if (!contains(pgrp))
    {
        cerr << "Current group: " << endl;
        this->print(cerr); cerr << endl;
        cerr << "Quotient group: " << endl;
        cerr << pgrp << endl;
        yeti_throw(SanityCheckError, "parameter subgroup for quotient is not a subgroup");
    }

    make(qgrp, PermutationSet, nindex());
    PermutationSet::iterator it_numer(begin());

    for ( ; it_numer != end(); ++it_numer)
    {
        Permutation* e(*it_numer);

        bool has_match = false;
        PermutationSet::iterator it_denom(pgrp->begin());
        for ( ; it_denom != pgrp->end(); ++it_denom)
        {
            Permutation* f(*it_denom);
            Permutation* product(f->product(e));
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

Permutation*
Permutation::get_identity(usi nindex)
{
    return identities_[nindex];
}

Permutation*
Permutation::get_permutation(short sign, usi nindex, uli idx)
{
    return sign == 1 ? plus_perms_[nindex][idx] : minus_perms_[nindex][idx];
}

Permutation*
Permutation::get_permutation(short sign, usi nindex, const usi* indexmap)
{
    uli idx = Permutation::id_index(nindex, indexmap);
    return get_permutation(sign, nindex, idx);
}

Permutation*
Permutation::get_permutation(
        short sign,
        usi nsubindex,
        usi nindex,
        const usi* indexmap,
        const usi* subset
)
{
    usi pmap[NINDEX];
    for (usi i=0; i < nindex; ++i)
        pmap[i] = i;

    for (usi i=0; i < nsubindex; ++i)
        pmap[subset[i]] = subset[indexmap[i]];

    return get_permutation(sign, nindex, pmap);
}

Permutation*
Permutation::get_permutation(
    short sign,
    usi nsubindex,
    usi nindex,
    const usi* indexmap,
    expand_t type
)
{
    usi pmap[NINDEX];
    for (usi i=0; i < nindex; ++i)
        pmap[i] = i;

    usi offset = 0;
    usi stop = nsubindex;
    if (type == ExpandBackward)
    {
        offset = nindex - nsubindex;
        stop = nindex;
    }

    const usi* pmapptr = indexmap;
    for (usi i=offset; i < stop; ++i, ++pmapptr)
        pmap[i] = *pmapptr + offset;

    return get_permutation(sign, nindex, pmap);
}

void
Permutation::build_permutations(usi nindex)
{
    ::memset(minus_perms_[nindex], 0, sizeof(Permutation*) * NPERMUTATIONS);
    ::memset(plus_perms_[nindex], 0, sizeof(Permutation*) * NPERMUTATIONS);

    usi available_indices[NINDEX];
    usi index_map[NINDEX];
    for (usi i=0; i < nindex; ++i)
        available_indices[i] = i;

    Permutation::next_permutation_index(nindex, 0, available_indices, index_map);

    //now compute all inverses
    usi new_index_map[NINDEX];
    uli nperms = factorial_[nindex];
    for (uli i=0; i < nperms; ++i)
    {
        Permutation* p = plus_perms_[nindex][i];
        const usi* pmap = p->indexmap();
        for (usi j=0; j < nindex; ++j)
        {
            new_index_map[pmap[j]] = j;
        }
        uli index = Permutation::id_index(nindex, new_index_map);
        inverses_[nindex][i] = index;
    }

    Permutation* identity = plus_perms_[nindex][0];
    identities_[nindex] = identity;


    usi product_map[NINDEX];
    for (uli i=0; i < nperms; ++i)
    {
        Permutation* pi = plus_perms_[nindex][i];
        for (uli j=0; j < nperms; ++j)
        {
            Permutation* pj = plus_perms_[nindex][j];
            pi->permute(pj->indexmap(), product_map);
            uli index = Permutation::id_index(nindex, product_map);
            products_[nindex][i][j] = index;
        }
    }
}

void
Permutation::next_permutation_index(
    usi nindex,
    usi index_number,
    const usi* available_indices,
    usi* index_map
)
{
    if (nindex == index_number)
        build_permutation(nindex, index_map);

    usi nremaining = nindex - index_number;
    usi next_available_indices[NINDEX];
    for (usi i=0; i < nremaining; ++i)
    {
        index_map[index_number] = available_indices[i];

        usi avail_idx = 0;
        for (usi j=0;  j < i; ++j, ++avail_idx)
        {
            next_available_indices[avail_idx] = available_indices[j];
        }
        for (usi j=i+1;  j < nindex; ++j, ++avail_idx)
        {
            next_available_indices[avail_idx] = available_indices[j];
        }

        next_permutation_index(nindex, index_number + 1, next_available_indices, index_map);
    }
}

void
Permutation::build_permutation(
    usi nindex,
    const usi* index_map
)
{
    {
        short plus = 1;
        Permutation* p = new (Permutation::ValidPermutationMalloc)
                        Permutation(nindex, plus, index_map);
        add_permutation(p);
    }

    {
        short minus = -1;
        Permutation* p = new (Permutation::ValidPermutationMalloc)
                        Permutation(nindex, minus, index_map);
        add_permutation(p);
    }
}

void
Permutation::add_permutation(Permutation* p)
{
    if (p->sign() == 1)
        plus_perms_[p->nindex()][p->id_index()] = p;
    else
        minus_perms_[p->nindex()][p->id_index()] = p;
    p->incref();
}

void*
Permutation::operator new(size_t size)
{
    cerr << "cannot construct a permutation this way!" << endl;
    abort();
}

void*
Permutation::operator new(size_t size, Permutation::permutation_malloc_flag_t flag)
{
    return ::malloc(size);
}

void
Permutation::operator delete(void* ptr, Permutation::permutation_malloc_flag_t flag)
{
    ::free(ptr);
}

void
Permutation::operator delete(void* ptr)
{
    ::free(ptr);
}
