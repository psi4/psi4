#include "index.h"
#include "mobasis.h"
#include "tensorparser.h"
#include "tensor.h"
#include "exception.h"
#include "runtime.h"

#define N_META_INDICES 3

using namespace yeti;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

MOBasisRangeBuilder::MOBasisRangeBuilder(
    usi nlayers_extra,
    uli nidx_node_layer,
    uli nidx_thread_layer,
    uli nidx_occ_data_layer,
    uli nidx_vir_data_layer,
    uli ncore,
    uli ndocc,
    uli nsocc,
    uli nvir,
    uli ncabs
)
    :
    spin_orbital_debug_(false),
    nlayers_extra_(nlayers_extra),
    nidx_per_tile_node_layer_(nidx_node_layer),
    nidx_per_tile_thread_layer_(nidx_thread_layer),
    nidx_per_tile_occ_data_layer_(nidx_occ_data_layer),
    nidx_per_tile_vir_data_layer_(nidx_vir_data_layer),
    ncore_pi_(new uli[1]),
    ndocc_pi_(new uli[1]),
    nsocc_pi_(new uli[1]),
    nvir_pi_(new uli[1]),
    ncabs_pi_(new uli[1]),
        core_(0),
        act_docc_(0),
        docc_(0),
        socc_(0),
        occ_(0),
        core_a_(0),
        act_docc_a_(0),
        vir_a_(0),
        cabs_a_(0),
        orb_a_(0),
        ri_a_(0),
        core_b_(0),
        act_docc_b_(0),
        vir_b_(0),
        cabs_b_(0),
        orb_b_(0),
        ri_b_(0),
        orb_(0),
        vir_(0),
        cabs_(0),
        ri_(0),
        nirrep_(1)
{
    ncore_pi_[0] = ncore;
    ndocc_pi_[0] = ndocc;
    nsocc_pi_[0] = nsocc;
    nvir_pi_[0] = nvir;
    ncabs_pi_[0] = ncabs;
}

MOBasisRangeBuilder::MOBasisRangeBuilder(
    usi nlayers_extra,
    uli nidx_node_layer,
    uli nidx_thread_layer,
    uli nidx_occ_data_layer,
    uli nidx_vir_data_layer,
    usi nirrep,
    const uli* ncore_pi,
    const uli* ndocc_pi,
    const uli* nsocc_pi,
    const uli* nvir_pi,
    const uli* ncabs_pi
)
    :
      spin_orbital_debug_(false),
      nlayers_extra_(nlayers_extra),
      nidx_per_tile_node_layer_(nidx_node_layer),
      nidx_per_tile_thread_layer_(nidx_thread_layer),
      nidx_per_tile_occ_data_layer_(nidx_occ_data_layer),
      nidx_per_tile_vir_data_layer_(nidx_vir_data_layer),
    ncore_pi_(new uli[nirrep]),
    ndocc_pi_(new uli[nirrep]),
    nsocc_pi_(new uli[nirrep]),
    nvir_pi_(new uli[nirrep]),
    ncabs_pi_(new uli[nirrep]),
        core_(0),
        act_docc_(0),
        docc_(0),
        socc_(0),
        occ_(0),
        orb_(0),
        vir_(0),
        cabs_(0),
        ri_(0),
        core_a_(0),
        act_docc_a_(0),
        docc_a_(0),
        occ_a_(0),
        vir_a_(0),
        cabs_a_(0),
        orb_a_(0),
        ri_a_(0),
        core_b_(0),
        act_docc_b_(0),
        docc_b_(0),
        occ_b_(0),
        vir_b_(0),
        cabs_b_(0),
        orb_b_(0),
        ri_b_(0),
        nirrep_(nirrep)
{
    for (usi i=0; i < nirrep; ++i)
    {
        ncore_pi_[i] = ncore_pi ? ncore_pi[i] : 0;
        ndocc_pi_[i] = ndocc_pi ? ndocc_pi[i] : 0;
        nsocc_pi_[i] = nsocc_pi ? nsocc_pi[i] : 0;
        nvir_pi_[i] =  nvir_pi ? nvir_pi[i] : 0;

        //optional array
        ncabs_pi_[i] = ncabs_pi ? ncabs_pi[i] : 0;
    }
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
MOBasisRangeBuilder::append(
    uli& offset,
    const SubindexTuplePtr& tuple,
    const IndexRangePtr& range
)
{
    SubindexTuplePtr subtuple = range->get_subranges();
    for (uli idx=0; idx < range->nelements(); ++idx, ++offset)
    {
        tuple->set(offset, subtuple->get(idx));
    }
}

void
MOBasisRangeBuilder::build_spin_orbital()
{
    spin_orbital_debug_ = true;
    build();
}

void
MOBasisRangeBuilder::build_spin_restricted()
{
    spin_orbital_debug_ = false;
    build();
}

void
MOBasisRangeBuilder::build(
    IndexRangePtr& composite_range,
    const IndexRangePtr& sp1,
    const IndexRangePtr& sp2,
    const IndexRangePtr& sp3,
    const IndexRangePtr& sp4,
    const IndexRangePtr& sp5,
    const IndexRangePtr& sp6,
    const IndexRangePtr& sp7,
    const IndexRangePtr& sp8,
    const IndexRangePtr& sp9,
    const IndexRangePtr& sp10,
    const IndexRangePtr& sp11,
    const IndexRangePtr& sp12
)
{

    //make the occ space
    uli nsubranges = 0;
    if (sp1) nsubranges += sp1->nelements();
    if (sp2) nsubranges += sp2->nelements();
    if (sp3) nsubranges += sp3->nelements();
    if (sp4) nsubranges += sp4->nelements();
    if (sp5) nsubranges += sp5->nelements();
    if (sp6) nsubranges += sp6->nelements();
    if (sp7) nsubranges += sp7->nelements();
    if (sp8) nsubranges += sp8->nelements();
    if (sp9) nsubranges += sp9->nelements();
    if (sp10) nsubranges += sp10->nelements();
    if (sp11) nsubranges += sp11->nelements();
    if (sp12) nsubranges += sp12->nelements();

    SubindexTuplePtr tuple = new SubindexTuple(nsubranges);
    uli offset = 0;
    if (sp1) append(offset, tuple, sp1);
    if (sp2) append(offset, tuple, sp2);
    if (sp3) append(offset, tuple, sp3);
    if (sp4) append(offset, tuple, sp4);
    if (sp5) append(offset, tuple, sp5);
    if (sp6) append(offset, tuple, sp6);
    if (sp7) append(offset, tuple, sp7);
    if (sp8) append(offset, tuple, sp8);
    if (sp9) append(offset, tuple, sp9);
    if (sp10) append(offset, tuple, sp10);
    if (sp11) append(offset, tuple, sp11);
    if (sp12) append(offset, tuple, sp12);

    composite_range = new IndexRange(tuple);
    //set as totally symmetric
    composite_range->set_irrep(0);
}

void
MOBasisRangeBuilder::build()
{
    for (usi h=0; h < nirrep_; ++h)
    {
        ncore_ += ncore_pi_[h];
        nact_docc_ += ndocc_pi_[h];
        nsocc_ += nsocc_pi_[h];
        nvir_ += nvir_pi_[h];
        ncabs_ += ncabs_pi_[h];
    }
    ndocc_ = ncore_ + nact_docc_;
    norb_ = ncore_ + nact_docc_ + nsocc_ + nvir_;
    nri_ = norb_ + ncabs_;

    uli topoffset = 0;
    build(topoffset, nidx_per_tile_occ_data_layer_, ncore_pi_, core_);
    build(topoffset, nidx_per_tile_occ_data_layer_, ndocc_pi_, act_docc_);
    build(topoffset, nidx_per_tile_occ_data_layer_, nsocc_pi_, socc_);
    build(topoffset, nidx_per_tile_vir_data_layer_, nvir_pi_, vir_);
    build(topoffset, nidx_per_tile_vir_data_layer_, ncabs_pi_, cabs_);




    if (core_ || act_docc_)
        build(docc_, core_, act_docc_);
    if (core_ || act_docc_ || socc_)
        build(occ_, core_, act_docc_, socc_);
    build(orb_, core_, act_docc_, socc_, vir_);
    build(ri_, core_, act_docc_, socc_, vir_, cabs_);



    {
        IndexRangePtr nullrange = 0;
        std::vector<IndexRangePtr> nonnull_ranges(6, nullrange);
        usi nonnull_count = 0;
        if (core_) { nonnull_ranges[nonnull_count] = core_; ++nonnull_count;}
        if (act_docc_) { nonnull_ranges[nonnull_count] = act_docc_; ++nonnull_count;}
        if (socc_) { nonnull_ranges[nonnull_count] = socc_; ++nonnull_count;}
        if (vir_) { nonnull_ranges[nonnull_count] = vir_; ++nonnull_count;}
        if (cabs_) { nonnull_ranges[nonnull_count] = cabs_; ++nonnull_count;}
        for (usi i=1; i < nonnull_count; ++i)
        {
            nonnull_ranges[i]->offset(nonnull_ranges[i-1]);
        }
    }

    if (spin_orbital_debug_)
    {
        topoffset = 0;
        build(topoffset, nidx_per_tile_occ_data_layer_, ncore_pi_, core_a_);
        build(topoffset, nidx_per_tile_occ_data_layer_, ncore_pi_, core_b_);
        build(topoffset, nidx_per_tile_occ_data_layer_, ndocc_pi_, act_docc_a_);
        build(topoffset, nidx_per_tile_occ_data_layer_, ndocc_pi_, act_docc_b_);
        build(topoffset, nidx_per_tile_occ_data_layer_, nsocc_pi_, socc_a_);
        build(topoffset, nidx_per_tile_occ_data_layer_, nsocc_pi_, socc_b_);
        build(topoffset, nidx_per_tile_vir_data_layer_, nvir_pi_, vir_a_);
        build(topoffset, nidx_per_tile_vir_data_layer_, nvir_pi_, vir_b_);
        build(topoffset, nidx_per_tile_vir_data_layer_, ncabs_pi_, cabs_a_);
        build(topoffset, nidx_per_tile_vir_data_layer_, ncabs_pi_, cabs_b_);

        if (core_a_ || core_b_)
            build(spin_orb_core_, core_a_, core_b_);
        build(spin_orb_docc_, act_docc_a_, act_docc_b_);
        build(spin_orb_vir_, vir_a_, vir_b_);
        build(spin_orb_orb_, core_a_, core_b_, act_docc_a_, act_docc_b_, vir_a_, vir_b_);
        if (cabs_a_ || cabs_b_)
            build(spin_orb_cabs_, cabs_a_, cabs_b_);
        build(spin_orb_ri_, core_a_, core_b_, act_docc_a_, act_docc_b_, vir_a_, vir_b_, cabs_a_, cabs_b_);

        if (spin_orb_core_)
        {
            spin_orb_docc_->set_offset(spin_orb_core_->nelements());
            spin_orb_vir_->set_offset(spin_orb_docc_->nelements() + spin_orb_core_->nelements());
        }
        else
        {
            spin_orb_vir_->set_offset(spin_orb_docc_->nelements());
        }

        YetiRuntime::register_subranges(
            spin_orb_orb_.get(),
            spin_orb_docc_.get(),
            spin_orb_vir_.get(),
            act_docc_a_.get(),
            act_docc_b_.get(),
            vir_a_.get(),
            vir_b_.get()
        );

        YetiRuntime::register_subranges(
            spin_orb_docc_.get(),
            act_docc_a_.get(),
            act_docc_b_.get()
        );

        YetiRuntime::register_subranges(
            spin_orb_vir_.get(),
            vir_a_.get(),
            vir_b_.get()
        );

        IndexRangePtr nullrange = 0;
        std::vector<IndexRangePtr> nonnull_ranges(12, nullrange);
        usi nonnull_count = 0;
        if (core_a_) { nonnull_ranges[nonnull_count] = core_a_; ++nonnull_count;}
        if (core_b_) { nonnull_ranges[nonnull_count] = core_b_; ++nonnull_count;}
        if (act_docc_a_) { nonnull_ranges[nonnull_count] = act_docc_a_; ++nonnull_count;}
        if (act_docc_b_) { nonnull_ranges[nonnull_count] = act_docc_b_; ++nonnull_count;}
        if (socc_a_) { nonnull_ranges[nonnull_count] = socc_a_; ++nonnull_count;}
        if (socc_b_) { nonnull_ranges[nonnull_count] = socc_b_; ++nonnull_count;}
        if (vir_a_) { nonnull_ranges[nonnull_count] = vir_a_; ++nonnull_count;}
        if (vir_b_) { nonnull_ranges[nonnull_count] = vir_b_; ++nonnull_count;}
        if (cabs_a_) { nonnull_ranges[nonnull_count] = cabs_a_; ++nonnull_count;}
        if (cabs_b_) { nonnull_ranges[nonnull_count] = cabs_b_; ++nonnull_count;}
        for (usi i=1; i < nonnull_count; ++i)
        {
            nonnull_ranges[i]->offset(nonnull_ranges[i-1]);
        }

    }
    else
    {
        topoffset = 0;
        build(topoffset, nidx_per_tile_occ_data_layer_, ncore_pi_, core_a_);
        build(topoffset, nidx_per_tile_occ_data_layer_, ndocc_pi_, act_docc_a_);
        build(topoffset, nidx_per_tile_occ_data_layer_, nsocc_pi_, socc_a_);
        build(topoffset, nidx_per_tile_vir_data_layer_, nvir_pi_, vir_a_);
        build(topoffset, nidx_per_tile_vir_data_layer_, ncabs_pi_, cabs_a_);

        build(topoffset, nidx_per_tile_occ_data_layer_, ncore_pi_, core_b_);
        build(topoffset, nidx_per_tile_occ_data_layer_, ndocc_pi_, act_docc_b_);
        build(topoffset, nidx_per_tile_occ_data_layer_, nsocc_pi_, socc_b_);
        build(topoffset, nidx_per_tile_vir_data_layer_, nvir_pi_, vir_b_);
        build(topoffset, nidx_per_tile_vir_data_layer_, ncabs_pi_, cabs_b_);

        if (core_a_ || act_docc_a_)
            build(docc_a_, core_a_, act_docc_a_);
        if (core_a_ || act_docc_a_ || socc_a_)
            build(occ_a_, core_a_, act_docc_a_, socc_a_);
        build(orb_a_, core_a_, act_docc_a_, socc_a_, vir_a_);
        build(ri_a_, core_a_, act_docc_a_, socc_a_, vir_a_, cabs_a_);

        if (core_b_ || act_docc_b_)
            build(docc_b_, core_b_, act_docc_b_);
        if (core_b_ || act_docc_b_ || socc_b_)
            build(occ_b_, core_b_, act_docc_b_, socc_b_);
        build(orb_b_, core_b_, act_docc_b_, socc_b_, vir_b_);
        build(ri_b_, core_b_, act_docc_b_, socc_b_, vir_b_, cabs_b_);

        IndexRangePtr nullrange = 0;
        std::vector<IndexRangePtr> nonnull_ranges(12, nullrange);
        usi nonnull_count = 0;
        if (core_a_) { nonnull_ranges[nonnull_count] = core_a_; ++nonnull_count;}
        if (act_docc_a_) { nonnull_ranges[nonnull_count] = act_docc_a_; ++nonnull_count;}
        if (socc_a_) { nonnull_ranges[nonnull_count] = socc_a_; ++nonnull_count;}
        if (vir_a_) { nonnull_ranges[nonnull_count] = vir_a_; ++nonnull_count;}
        if (cabs_a_) { nonnull_ranges[nonnull_count] = cabs_a_; ++nonnull_count;}
        if (core_b_) { nonnull_ranges[nonnull_count] = core_b_; ++nonnull_count;}
        if (act_docc_b_) { nonnull_ranges[nonnull_count] = act_docc_b_; ++nonnull_count;}
        if (socc_b_) { nonnull_ranges[nonnull_count] = socc_b_; ++nonnull_count;}
        if (vir_b_) { nonnull_ranges[nonnull_count] = vir_b_; ++nonnull_count;}
        if (cabs_b_) { nonnull_ranges[nonnull_count] = cabs_b_; ++nonnull_count;}
        for (usi i=1; i < nonnull_count; ++i)
        {
            nonnull_ranges[i]->offset(nonnull_ranges[i-1]);
        }
    }


    if (core_)
    {
        YetiRuntime::register_subranges(docc_.get(), core_.get());
        YetiRuntime::register_subranges(orb_.get(), core_.get());
        YetiRuntime::register_subranges(ri_.get(), core_.get());
    }
    if (act_docc_)
    {
        YetiRuntime::register_subranges(docc_.get(), act_docc_.get());
        YetiRuntime::register_subranges(orb_.get(), act_docc_.get());
        YetiRuntime::register_subranges(ri_.get(), act_docc_.get());
    }
    if (docc_)
    {
        YetiRuntime::register_subranges(orb_.get(), docc_.get());
        YetiRuntime::register_subranges(ri_.get(), docc_.get());
    }
    if (socc_)
    {
        YetiRuntime::register_subranges(orb_.get(), socc_.get());
        YetiRuntime::register_subranges(ri_.get(), socc_.get());    
    }
    if (vir_)
    {
        YetiRuntime::register_subranges(orb_.get(), vir_.get());
        YetiRuntime::register_subranges(ri_.get(), vir_.get());    
    }
    if (cabs_)
    {
        YetiRuntime::register_subranges(ri_.get(), cabs_.get());    
    }
    YetiRuntime::register_subranges(ri_.get(), orb_.get());
}



void
MOBasisRangeBuilder::build(
    uli& topoffset,
    uli nidx_per_tile_data_layer,
    uli* norbs_per_irrep,
    IndexRangePtr& space
)
{
    uli ntot = 0;
    for (usi h=0; h < nirrep_; ++h)
    {
        ntot += norbs_per_irrep[h];
    }
    if (ntot == 0)
    {
        space = 0; //set to null
        return;
    }

    SubindexTuplePtr irrep_ranges(new SubindexTuple(nirrep_));
    uli nsubranges = 0;
    for (usi h=0; h < nirrep_; ++h)
    {
        uli norbs_in_irrep = norbs_per_irrep[h];
        MORangeBuilder morange(
            nlayers_extra_,
            nidx_per_tile_node_layer_,
            nidx_per_tile_thread_layer_,
            nidx_per_tile_data_layer,
            norbs_in_irrep,
            h
        );
        IndexRangePtr irrep_range = morange.get();
        if (irrep_range)
        {
            nsubranges += irrep_range->nelements();
            irrep_range->set_irrep(h);
            irrep_ranges->set(h, irrep_range.get());
        }
    }

    IndexRange* prev = irrep_ranges->get(0);
    for (usi h=1; h < nirrep_; ++h)
    {
        IndexRange* next = irrep_ranges->get(h);
        if (next && prev)
        {
            next->offset(prev);
            prev = next;
        }
    }

    SubindexTuplePtr subranges = 0;
    if (nirrep_ > 1)
    {
        subranges = new SubindexTuple(nsubranges);
        uli idx = 0;
        for (usi h=0; h < nirrep_; ++h)
        {
            IndexRange* irrep_range = irrep_ranges->get(h);
            if (!irrep_range)
                continue;

            uli nsubrange = irrep_range->nelements();
            for (uli subidx=0; subidx < nsubrange; ++subidx, ++idx)
            {
                IndexRange* subrange = irrep_range->get_subindex(subidx);
                subranges->set(idx, subrange);
                subrange->set_irrep(h);
            }
        }
    }
    else
    {
        subranges = irrep_ranges->get(0)->get_subranges();
    }

    space = new IndexRange(topoffset, subranges);
    //set as totally symmetric
    space->set_irrep(0);

    topoffset += space->nelements();


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
MOBasisRangeBuilder::get_cabs_range() const
{
    return cabs_.get();
}

IndexRange*
MOBasisRangeBuilder::get_ri_range() const
{
    return ri_.get();
}

IndexRange*
MOBasisRangeBuilder::get_spin_orbital_docc() const
{
    return spin_orb_docc_.get();
}

IndexRange*
MOBasisRangeBuilder::get_spin_orbital_vir() const
{
    return spin_orb_vir_.get();
}

IndexRange*
MOBasisRangeBuilder::get_spin_orbital_orb() const
{
    return spin_orb_orb_.get();
}

IndexRange*
MOBasisRangeBuilder::get_spin_orbital_cabs() const
{
    return spin_orb_cabs_.get();
}

IndexRange*
MOBasisRangeBuilder::get_spin_orbital_ri() const
{
    return spin_orb_ri_.get();
}

IndexRange*
MOBasisRangeBuilder::get_act_docc_alpha() const
{
    return act_docc_a_.get();
}

IndexRange*
MOBasisRangeBuilder::get_act_docc_beta() const
{
    return act_docc_b_.get();
}

IndexRange*
MOBasisRangeBuilder::get_vir_alpha() const
{
    return vir_a_.get();
}

IndexRange*
MOBasisRangeBuilder::get_vir_beta() const
{
    return vir_b_.get();
}

IndexRange*
MOBasisRangeBuilder::get_orb_alpha() const
{
    return orb_a_.get();
}

IndexRange*
MOBasisRangeBuilder::get_orb_beta() const
{
    return orb_b_.get();
}

IndexRange*
MOBasisRangeBuilder::get_cabs_alpha() const
{
    return cabs_a_.get();
}

IndexRange*
MOBasisRangeBuilder::get_cabs_beta() const
{
    return cabs_b_.get();
}

IndexRange*
MOBasisRangeBuilder::get_ri_alpha() const
{
    return ri_a_.get();
}

IndexRange*
MOBasisRangeBuilder::get_ri_beta() const
{
    return ri_b_.get();
}



void
MOBasisRangeBuilder::set_debug(bool flag)
{
    spin_orbital_debug_ = flag;
}

MORangeBuilder::MORangeBuilder(
    usi nlayers_extra,
    uli nidx_node_layer,
    uli nidx_thread_layer,
    uli nidx_data_layer,
    uli ntot,
    usi irrep
)
    :
   nlayers_extra_(nlayers_extra),
   nidx_per_tile_node_layer_(nidx_node_layer),
   nidx_per_tile_thread_layer_(nidx_thread_layer),
   nidx_per_tile_data_layer_(nidx_data_layer),
   ntot_(ntot),
   irrep_(irrep)
{
}

IndexRangePtr
MORangeBuilder::get()
{
    //create the first index range

    if (ntot_ == 0) //null range
        return 0;

    uli start = 0;
    uli nranges_data = ntot_ / nidx_per_tile_data_layer_;
    if (nranges_data == 0) ++nranges_data;

    IndexRange* range = new IndexRange(start, ntot_, nidx_per_tile_data_layer_);
    SubindexTuplePtr tuple = range->get_subranges();
    if (!tuple) //no subranges
    {
        raise(SanityCheckError, "mo range must have subranges!");
    }

    for (usi d=0; d < nlayers_extra_; ++d)
    {
        SubindexTuplePtr subranges = range->get_subranges();
        for (usi idx=0; idx < range->nelements(); ++idx)
        {
            SubindexTuplePtr metatuple = new SubindexTuple(1);
            IndexRange* metarange = new IndexRange(idx, metatuple);
            metatuple->set(0, subranges->get(idx));
            subranges->set(idx, metarange);
        }
        range = new IndexRange(subranges);
    }

    IndexRangePtr thread_range = new IndexRange(0, range->get_subranges(), nidx_per_tile_thread_layer_);

    IndexRangePtr node_range = new IndexRange(0, thread_range->get_subranges(), nidx_per_tile_node_layer_);
    node_range->set_irrep(irrep_);

    return node_range;
}
