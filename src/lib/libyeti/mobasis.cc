#include "index.h"
#include "mobasis.h"
#include "tensorparser.h"
#include "tensor.h"

#define N_META_INDICES 3

using namespace yeti;
using namespace std;

MOBasisRangeBuilder::MOBasisRangeBuilder(
    uli nocc_tile,
    uli nvir_tile,
    uli ncore,
    uli ndocc,
    uli nsocc,
    uli nvir,
    uli ncabs
)
    :
    ncore_pi_(new uli[1]),
    ndocc_pi_(new uli[1]),
    nsocc_pi_(new uli[1]),
    nvir_pi_(new uli[1]),
    ncabs_pi_(new uli[1]),
        nocc_tile_(nocc_tile),
        nvir_tile_(nvir_tile),
        core_(0),
        act_docc_(0),
        docc_(0),
        socc_(0),
        occ_(0),
        act_a_occ_(0),
        act_b_occ_(0),
        orb_(0),
        vir_(0),
        act_a_vir_(0),
        act_b_vir_(0),
        cabs_(0),
        ri_(0),
        nirrep_(1)
{
    ncore_pi_[0] = ncore;
    ndocc_pi_[0] = ndocc;
    nsocc_pi_[0] = nsocc;
    nvir_pi_[0] = nvir;
    ncabs_pi_[0] = ncabs;
    build();
}

MOBasisRangeBuilder::MOBasisRangeBuilder(
    uli nocc_tile,
    uli nvir_tile,
    usi nirrep,
    const uli* ncore_pi,
    const uli* ndocc_pi,
    const uli* nsocc_pi,
    const uli* nvir_pi,
    const uli* ncabs_pi

)
    :
    ncore_pi_(new uli[nirrep]),
    ndocc_pi_(new uli[nirrep]),
    nsocc_pi_(new uli[nirrep]),
    nvir_pi_(new uli[nirrep]),
    ncabs_pi_(new uli[nirrep]),
        nocc_tile_(nocc_tile),
        nvir_tile_(nvir_tile),
        core_(0),
        act_docc_(0),
        docc_(0),
        socc_(0),
        occ_(0),
        act_a_occ_(0),
        act_b_occ_(0),
        orb_(0),
        vir_(0),
        act_a_vir_(0),
        act_b_vir_(0),
        cabs_(0),
        ri_(0),
        nirrep_(nirrep)
{
    for (usi i=0; i < nirrep; ++i)
    {
        ncore_pi_[i] = ncore_pi[i];
        ndocc_pi_[i] = ndocc_pi[i];
        nsocc_pi_[i] = nsocc_pi[i];
        nvir_pi_[i] = nvir_pi[i];

        //optional array
        ncabs_pi_[i] = ncabs_pi ? ncabs_pi[i] : 0;
    }
    build();
}


MOBasisRangeBuilder::MOBasisRangeBuilder(
    uli nocc_tile,
    uli nvir_tile,
    uli ndocc,
    uli nvir
)
    :
    ncore_pi_(new uli[1]),
    ndocc_pi_(new uli[1]),
    nsocc_pi_(new uli[1]),
    nvir_pi_(new uli[1]),
    ncabs_pi_(new uli[1]),
        nocc_tile_(nocc_tile),
        nvir_tile_(nvir_tile),
        core_(0),
        act_docc_(0),
        docc_(0),
        socc_(0),
        occ_(0),
        orb_(0),
        vir_(0),
        cabs_(0),
        ri_(0),
        nirrep_(1)
{
    ncore_pi_[0] = 0;
    ndocc_pi_[0] = ndocc;
    nsocc_pi_[0] = 0;
    nvir_pi_[0] = nvir;
    ncabs_pi_[0] = 0;
    build();
}

MOBasisRangeBuilder::~MOBasisRangeBuilder()
{
    delete[] ncore_pi_;
    delete[] ndocc_pi_;
    delete[] nsocc_pi_;
    delete[] nvir_pi_;
    delete[] ncabs_pi_;
}

void
MOBasisRangeBuilder::align_spaces()
{
    //at this point we need to align all of the spaces
    //to be at the same depth due to symmetry
    usi maxdepth = 0;
    if (core_ && core_->depth() > maxdepth)
        maxdepth = core_->depth();
    if (act_docc_ && act_docc_->depth() > maxdepth)
        maxdepth = act_docc_->depth();
    if (socc_ && socc_->depth() > maxdepth)
        maxdepth = socc_->depth();
    if (vir_ && vir_->depth() > maxdepth)
        maxdepth = vir_->depth();
    if (cabs_ && cabs_->depth() > maxdepth)
        maxdepth = cabs_->depth();

    usi subdepth = maxdepth - 1;

    if (core_ && core_->depth() < maxdepth)
        core_->expand_subrange_depth(subdepth);
    if (act_docc_ && act_docc_->depth() < maxdepth)
        act_docc_->expand_subrange_depth(subdepth);
    if (socc_ && socc_->depth() < maxdepth)
        socc_->expand_subrange_depth(subdepth);
    if (vir_ && vir_->depth() < maxdepth)
        vir_->expand_subrange_depth(subdepth);
    if (cabs_ && cabs_->depth() < maxdepth)
        cabs_->expand_subrange_depth(subdepth);

    //store the symmetry depth
    symmdepth_ = subdepth;
}


void
MOBasisRangeBuilder::build()
{
    build(nocc_tile_, ncore_pi_, core_);
    build(nocc_tile_, ndocc_pi_, act_docc_);
    build(nocc_tile_, nsocc_pi_, socc_);
    build(nvir_tile_, nvir_pi_, vir_);
    build(nvir_tile_, ncabs_pi_, cabs_);

    align_spaces();



    uli ndocc_spaces = 0;
    if (core_) ++ndocc_spaces;
    if (act_docc_) ++ndocc_spaces;
    //attempt to merge docc spaces
    SubindexTuplePtr docc_tuple(new SubindexTuple(ndocc_spaces));

    uli idx = 0;
    if (core_)
    {
        docc_tuple->set(idx, core_.get());
        ++idx;
    }
    if (act_docc_)
    {
        docc_tuple->set(idx, act_docc_.get());
        ++idx;
    }
    if (docc_tuple->size() == 1)
    {
        docc_ = docc_tuple->get(0);
    }
    else if (docc_tuple->size() > 1)
    {
        docc_ = new IndexRange(docc_tuple);
    }

    //make the occ space
    uli nocc_spaces = 0;
    if (core_) ++nocc_spaces;
    if (act_docc_) ++nocc_spaces;
    if (socc_) ++nocc_spaces;

    SubindexTuplePtr occ_tuple(new SubindexTuple(nocc_spaces));
    idx = 0;
    if (core_)
    {
        occ_tuple->set(idx, core_.get());
        ++idx;
    }
    if (act_docc_)
    {
        occ_tuple->set(idx, act_docc_.get());
        ++idx;
    }
    if (socc_)
    {
        occ_tuple->set(idx, socc_.get());
        ++idx;
    }
    if (occ_tuple->size() == 1)
        occ_ = occ_tuple->get(0);
    else if (occ_tuple->size() > 1)
        occ_ = new IndexRange(occ_tuple);

    //make the act a occ space
    uli naocc_spaces = 0;
    if (act_docc_) ++naocc_spaces;
    if (socc_) ++naocc_spaces;
    idx = 0;
    SubindexTuplePtr act_a_occ_tuple(new SubindexTuple(naocc_spaces));
    if (act_docc_)
    {
        act_a_occ_tuple->set(idx, act_docc_.get());
        ++idx;
    }
    if (socc_)
    {
        act_a_occ_tuple->set(idx, socc_.get());
        ++idx;
    }
    if (act_a_occ_tuple->size() == 1)
        act_a_occ_ = act_a_occ_tuple->get(0);
    else if (act_a_occ_tuple->size() > 1)
        act_a_occ_ = new IndexRange(act_a_occ_tuple);

    //make the act b occ space
    uli nbocc_spaces = 0;
    if (act_docc_) ++nbocc_spaces;
    idx = 0;
    SubindexTuplePtr act_b_occ_tuple(new SubindexTuple(nbocc_spaces));
    if (act_docc_)
    {
        act_b_occ_tuple->set(idx, act_docc_.get());
        ++idx;
    }
    if (act_b_occ_tuple->size() == 1)
        act_b_occ_ = act_b_occ_tuple->get(0);
    else if (act_b_occ_tuple->size() > 1)
        act_b_occ_ = new IndexRange(act_b_occ_tuple);

    //make the act a vir space
    uli navir_spaces = 0;
    if (vir_) ++navir_spaces;
    idx = 0;
    SubindexTuplePtr act_a_vir_tuple(new SubindexTuple(navir_spaces));
    if (vir_){
        act_a_vir_tuple->set(idx, vir_.get());
        ++idx;
    }
    if (act_a_vir_tuple->size() == 1)
        act_a_vir_ = act_a_vir_tuple->get(0);
    else if (act_a_vir_tuple->size() > 1)
        act_a_vir_ = new IndexRange(act_a_vir_tuple);

    //make the act b vir space
    uli nbvir_spaces = 0;
    if (socc_) ++nbvir_spaces;
    if (vir_) ++nbvir_spaces;
    idx = 0;
    SubindexTuplePtr act_b_vir_tuple(new SubindexTuple(nbvir_spaces));
    if (socc_){
        act_b_vir_tuple->set(idx, socc_.get());
        idx++;
    }
    if (vir_)
    {
        act_b_vir_tuple->set(idx, vir_.get());
        idx++;
    }
    if (act_b_vir_tuple->size() == 1)
        act_b_vir_ = act_b_vir_tuple->get(0);
    else if (act_b_vir_tuple->size() > 1)
        act_b_vir_ = new IndexRange(act_b_vir_tuple);

    //make the orb space
    uli norb_index = 0;
    if (core_) ++norb_index;
    if (act_docc_) ++norb_index;
    if (socc_) ++norb_index;
    if (vir_) ++norb_index;


    SubindexTuplePtr orb_tuple(new SubindexTuple(norb_index));
    idx = 0;
    if (core_)
    {
        orb_tuple->set(idx, core_.get());
        ++idx;
    }
    if (act_docc_)
    {
        orb_tuple->set(idx, act_docc_.get());
        ++idx;
    }
    if (socc_)
    {
        orb_tuple->set(idx, socc_.get());
        ++idx;
    }
    if (vir_)
    {
        orb_tuple->set(idx, vir_.get());
        ++idx;
    }

    if (orb_tuple->size() == 1)
        orb_ = orb_tuple->get(0);
    else if (orb_tuple->size() > 1)
        orb_ = new IndexRange(orb_tuple);

    //make the orb space
    uli nri_index = 0;
    if (core_) ++nri_index;
    if (act_docc_) ++nri_index;
    if (socc_) ++nri_index;
    if (vir_) ++nri_index;
    if (cabs_) ++nri_index;
    SubindexTuplePtr ri_tuple(new SubindexTuple(nri_index));

    idx = 0;
    if (core_)
    {
        ri_tuple->set(idx, core_.get());
        ++idx;
    }
    if (act_docc_)
    {
        ri_tuple->set(idx, act_docc_.get());
        ++idx;
    }
    if (socc_)
    {
        ri_tuple->set(idx, socc_.get());
        ++idx;
    }
    if (vir_)
    {
        ri_tuple->set(idx, vir_.get());
        ++idx;
    }
    if (cabs_)
    {
        ri_tuple->set(idx, cabs_.get());
        ++idx;
    }
    if (ri_tuple->size() == 1)
        ri_ = orb_tuple->get(0);
    else if (ri_tuple->size() > 1)
        ri_ = new IndexRange(ri_tuple);

    ri_->increment_offsets();
}

void
MOBasisRangeBuilder::build(
    uli nper,
    uli ntot,
    IndexRangePtr& space
)
{
    if (ntot)
    {
        MORangeBuilder morange(nper, ntot);
        space = morange.get();
    }
}

void
MOBasisRangeBuilder::build(
    uli nper,
    uli* n_per_irrep,
    IndexRangePtr& space
)
{
    if (nirrep_ == 1)
    {
        build(nper, n_per_irrep[0],  space);
        return;
    }

    uli ntot = 0;
    for (usi h=0; h < nirrep_; ++h)
    {
        ntot += n_per_irrep[h];
    }
    if (ntot == 0)
    {
        space = 0; //set to null
        return;
    }

    SubindexTuplePtr subranges(new SubindexTuple(nirrep_,0));
    for (usi h=0; h < nirrep_; ++h)
    {
        uli nirrep = n_per_irrep[h];
        MORangeBuilder morange(nper, nirrep);
        IndexRange* irrep_range = morange.get();
        subranges->set(h, irrep_range);
    }
    space = new IndexRange(subranges);
    space->increment_offsets();
}

IndexRange*
MOBasisRangeBuilder::get_act_docc_range() const
{
    return act_docc_.get();
}

IndexRange*
MOBasisRangeBuilder::get_occ_range() const
{
    return occ_.get();
}

IndexRange*
MOBasisRangeBuilder::get_act_a_occ_range() const
{
    return act_a_occ_;
}

IndexRange*
MOBasisRangeBuilder::get_act_b_occ_range() const
{
    return act_b_occ_;
}

IndexRange*
MOBasisRangeBuilder::get_socc_range() const
{
    return socc_.get();
}

IndexRange*
MOBasisRangeBuilder::get_docc_range() const
{
    return docc_.get();
}

IndexRange*
MOBasisRangeBuilder::get_core_range() const
{
    return core_.get();
}

IndexRange*
MOBasisRangeBuilder::get_orb_range() const
{
    return orb_.get();
}

IndexRange*
MOBasisRangeBuilder::get_vir_range() const
{
    return vir_.get();
}

IndexRange*
MOBasisRangeBuilder::get_act_a_vir_range() const
{
    return act_a_vir_;
}

IndexRange*
MOBasisRangeBuilder::get_act_b_vir_range() const
{
    return act_b_vir_;
}

IndexRange*
MOBasisRangeBuilder::get_cabs_range() const
{
    return cabs_.get();
}

IndexRange*
MOBasisRangeBuilder::get_ri_range() const
{
    return ri_.get();
}

MOSymmetryMapPtr
MOBasisRangeBuilder::get_symmetry_map() const
{
    MOSymmetryMapPtr symm_map = new MOSymmetryMap(
                                    core_.get(),
                                    docc_.get(),
                                    socc_.get(),
                                    vir_.get(),
                                    nirrep_,
                                    symmdepth_
                                );
}

MORangeBuilder::MORangeBuilder(uli nper, uli ntot)
    : nper_(nper), ntot_(ntot)
{
}

IndexRange*
MORangeBuilder::get()
{
    //create the first index range
    uli start = 0;
    IndexRange* range = new IndexRange(start, ntot_, nper_);
    SubindexTuplePtr tuple = range->get_subranges();
    if (!tuple) //no subranges
        return range;

    while (tuple->size() > 2 * N_META_INDICES)
    {
        range = new IndexRange(start, tuple, N_META_INDICES);
        tuple = range->get_subranges();
    }
    return range;
}

MOSymmetryMap::MOSymmetryMap(
    IndexRange *core,
    IndexRange *docc,
    IndexRange *socc,
    IndexRange *vir,
    usi nirrep,
    usi symmdepth
) :
    nirrep_(nirrep),
    symmetry_depth_(symmdepth),
    nspaces_(0),
    irrep_map_(0)
{
    if (core)   ++nspaces_;
    if (docc)   ++nspaces_;
    if (socc)   ++nspaces_;
    if (vir)    ++nspaces_;

    irrep_map_ = new usi[nspaces_ * nirrep];

    usi idx = 0;
    for (usi i=0; i < nspaces_; ++i)
    {
        for (usi h=0; h < nirrep; ++h, ++idx)
        {
            irrep_map_[idx] = h;
        }
    }

}

