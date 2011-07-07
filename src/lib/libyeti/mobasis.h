#ifndef yeti_MOBASIS_H
#define yeti_MOBASIS_H


#include "class.h"

#include "index.hpp"
#include "mobasis.hpp"
#include "tensorparser.hpp"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

class MOBasisRangeBuilder :
    public smartptr::Countable
{

    private:
        bool spin_orbital_debug_;

        usi build_depth_;

        uli nblocks_top_layer_;

        usi nirrep_;

        usi symmdepth_;

        uli ncore_;

        uli nact_docc_;

        uli nsocc_;

        uli ndocc_;

        uli norb_;

        uli nvir_;

        uli ncabs_;

        uli nri_;

        uli* ncore_pi_;

        uli* ndocc_pi_;

        uli* nsocc_pi_;

        uli* nvir_pi_;

        uli* ncabs_pi_;

        uli nocc_tile_;

        uli nvir_tile_;

        IndexRangePtr cabs_;

        IndexRangePtr ri_;

        IndexRangePtr orb_;

        IndexRangePtr core_;

        IndexRangePtr docc_;

        IndexRangePtr occ_;

        IndexRangePtr act_docc_a_;

        IndexRangePtr core_a_;

        IndexRangePtr socc_a_;

        IndexRangePtr vir_a_;

        IndexRangePtr orb_a_;

        IndexRangePtr cabs_a_;

        IndexRangePtr ri_a_;

        IndexRangePtr occ_a_;

        IndexRangePtr docc_a_;

        IndexRangePtr act_docc_b_;

        IndexRangePtr core_b_;

        IndexRangePtr socc_b_;

        IndexRangePtr vir_b_;

        IndexRangePtr cabs_b_;

        IndexRangePtr ri_b_;

        IndexRangePtr orb_b_;

        IndexRangePtr occ_b_;

        IndexRangePtr docc_b_;

        IndexRangePtr vir_;

        IndexRangePtr act_docc_;

        IndexRangePtr socc_;

        IndexRangePtr spin_orb_core_;

        IndexRangePtr spin_orb_docc_;

        IndexRangePtr spin_orb_vir_;

        IndexRangePtr spin_orb_orb_;

        IndexRangePtr spin_orb_cabs_;

        IndexRangePtr spin_orb_ri_;

        void build();

        void build(
            uli& topoffset,
            uli nper,
            uli* n_per_irrep,
            IndexRangePtr& space
        );

        void
        append(
            uli& offset,
            const SubindexTuplePtr& tuple,
            const IndexRangePtr& range
        );

        void
        build(
            IndexRangePtr& composite_range,
            const IndexRangePtr& sp1 = 0,
            const IndexRangePtr& sp2 = 0,
            const IndexRangePtr& sp3 = 0,
            const IndexRangePtr& sp4 = 0,
            const IndexRangePtr& sp5 = 0,
            const IndexRangePtr& sp6 = 0,
            const IndexRangePtr& sp7 = 0,
            const IndexRangePtr& sp8 = 0,
            const IndexRangePtr& sp9 = 0,
            const IndexRangePtr& sp10 = 0,
            const IndexRangePtr& sp11 = 0,
            const IndexRangePtr& sp12 = 0
        );

    public:
        /**
            @param build_depth See MORangeBuilder
            @param nblocks_top See MORangeBuilder
            @param nocc_tile  The number of occupied orbitals per tile
            @param nvir_tile  The number of virtual orbitals per tile
            @param ncore
            @param ndocc
            @param nsocc
            @param nvir
        */
        MOBasisRangeBuilder(
            usi build_depth,
            uli nblocks_top,
            uli nocc_tile,
            uli nvir_tile,
            uli ncore,
            uli ndocc,
            uli nsocc,
            uli nvir,
            uli ncabs = 0
        );

        /**
            @param build_depth See MORangeBuilder
            @param nblocks_top See MORangeBuilder
            @param nocc_tile  The number of occupied orbitals per tile
            @param nvir_tile  The number of virtual orbitals per tile
            @param ncore
            @param ndocc
            @param nsocc
            @param nvir
        */
        MOBasisRangeBuilder(
            usi build_depth,
            uli nblocks_top,
            uli nocc_tile,
            uli nvir_tile,
            uli ndocc,
            uli nvir
        );

        MOBasisRangeBuilder(
            usi build_depth,
            uli nblocks_top,
            uli nocc_tile,
            uli nvir_tile,
            usi nirrep,
            const uli* ncore_pi,
            const uli* ndocc_pi,
            const uli* nsocc_pi,
            const uli* nvir_pi,
            const uli* ncabs_pi
        );

        ~MOBasisRangeBuilder();

        IndexRange* get_core_range() const;

        IndexRange* get_act_docc_range() const;

        IndexRange* get_occ_range() const;

        IndexRange* get_socc_range() const;

        IndexRange* get_docc_range() const;

        IndexRange* get_orb_range() const;

        IndexRange* get_vir_range() const;

        IndexRange* get_ri_range() const;

        IndexRange* get_cabs_range() const;

        IndexRange* get_spin_orbital_docc() const;

        IndexRange* get_spin_orbital_vir() const;

        IndexRange* get_spin_orbital_orb() const;

        IndexRange* get_spin_orbital_cabs() const;

        IndexRange* get_spin_orbital_ri() const;

        IndexRange* get_act_docc_alpha() const;

        IndexRange* get_vir_alpha() const;

        IndexRange* get_orb_alpha() const;

        IndexRange* get_cabs_alpha() const;

        IndexRange* get_ri_alpha() const;

        IndexRange* get_act_docc_beta() const;

        IndexRange* get_vir_beta() const;

        IndexRange* get_orb_beta() const;

        IndexRange* get_cabs_beta() const;

        IndexRange* get_ri_beta() const;

        MOSymmetryMapPtr get_symmetry_map() const;

        void set_debug(bool flag);

        void build_spin_orbital();

        void build_spin_restricted();

};

class MORangeBuilder : public smartptr::Countable {

    private:
        usi build_depth_;

        uli nblocks_top_layer_;

        uli nper_;

        uli ntot_;

        usi irrep_;

    public:
        /**
            @param build_depth The number of metadata layers to build. At minimum
                    two metadata layers are required - a "top" and "bottom" layer...
                    bottom being the data layer.  Extra metadata layers can be
                    inserted between the top and bottom levels.  Build depth is
                    equal to 2 plus the number of extra metadata layers.
            @param The number of index subranges to have per tile in the top layer
            @param The number of indices to have tiled on the bottom layer
            @param The total number of indices
        */
        MORangeBuilder(
            usi build_depth,
            uli ntop,
            uli nbottom,
            uli ntot,
            usi irrep
        );

        IndexRangePtr get();

};

}

#ifdef redefine_size_t
#undef size_t
#endif

#endif // MOBASIS_H
