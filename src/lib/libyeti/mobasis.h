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

        usi nlayers_extra_;

        uli nidx_per_tile_node_layer_;

        uli nidx_per_tile_thread_layer_;

        uli nidx_per_tile_occ_data_layer_;

        uli nidx_per_tile_vir_data_layer_;

        usi nirrep_;

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

        /**
            @param topoffset Reference return of the index offset for the top (node layer) index offset
            @param nidx_per_tile_data_layer
            @param norbs_per_irrep
            @param space Reference return of the index range corresponding to the given subspace
        */
        void build(
            uli& topoffset,
            uli nidx_per_tile_data_layer,
            uli* norbs_per_irrep,
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
            @param nlayers_extra See MORangeBuilder
            @param nidx_node_layer  See MORangeBuilder
            @param nidx_thread_layer  See MORangeBuilder
            @param nidx_occ_data_layer  The number of occupied orbitals per tile
            @param nidx_vir_data_layer  The number of virtual orbitals per tile
            @param ncore
            @param ndocc
            @param nsocc
            @param nvir
        */
        MOBasisRangeBuilder(
            usi nlayers_extra,
            uli nidx_node_layer,
            uli nidx_thread_layer,
            uli nidx_occ_data_layer,
            uli nidx_vir_data_layer,
            uli ncore,
            uli ndocc,
            uli nsocc,
            uli nvir,
            uli ncabs = 0
        );


        MOBasisRangeBuilder(
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

class MORangeBuilder :
        public smartptr::Countable
{

    private:
        uli nlayers_extra_;

        uli nidx_per_tile_node_layer_;

        uli nidx_per_tile_thread_layer_;

        uli nidx_per_tile_data_layer_;

        uli ntot_;

        usi irrep_;

    public:
        /**
            @param  nlayers_extra
                    At minimum three metadata layers are required.  The ''node'' layer creates
                    a grid for distribution of tensor blocks across nodes.  The ''thread'' layer
                    distributes tasks for threads within a given node.  The final layer or
                    ''day'' layer generates a grid for each individual data block.  Extra
                    metadata layers can be inserted between the bottom ''data'' layer and
                    the ''thread'' layer.
            @param nidx_node_layer   The number of indices per tile in the node layer
            @param nidx_thread_layer The number of indices per tile in the thread layer
            @param nidx_data_layer   The number of indices per tile in the data layer
            @param ntot              The total number of indices
            @param irrep             The irrep number for the given range
        */
        MORangeBuilder(
            usi nlayers_extra,
            uli nidx_node_layer,
            uli nidx_thread_layer,
            uli nidx_data_layer,
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
