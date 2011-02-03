#ifndef yeti_MOBASIS_H
#define yeti_MOBASIS_H


#include "class.h"

#include "index.hpp"
#include "mobasis.hpp"
#include "tensorparser.hpp"


namespace yeti {

class MOBasisRangeBuilder : public smartptr::Countable {

    private:

        usi nirrep_;

        usi symmdepth_;

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

        IndexRange* act_a_occ_;

        IndexRange* act_b_occ_;

        IndexRange* act_a_vir_;

        IndexRange* act_b_vir_;

        IndexRangePtr vir_;

        IndexRangePtr act_docc_;

        IndexRangePtr socc_;

        void build();

        void align_spaces();

        void build(uli nper, uli ntot, IndexRangePtr& space);

        void build(uli nper, uli* n_per_irrep, IndexRangePtr& space);

    public:
        /**
            @param nocc_tile  The number of occupied orbitals per tile
            @param nvir_tile  The number of virtual orbitals per tile
            @param ncore
            @param ndocc
            @param nsocc
            @param nvir
        */
        MOBasisRangeBuilder(
            uli nocc_tile,
            uli nvir_tile,
            uli ncore,
            uli ndocc,
            uli nsocc,
            uli nvir,
            uli ncabs = 0
        );

        /**
            @param nocc_tile  The number of occupied orbitals per tile
            @param nvir_tile  The number of virtual orbitals per tile
            @param ndocc
            @param nvir
        */
        MOBasisRangeBuilder(
            uli nocc_tile,
            uli nvir_tile,
            uli ndocc,
            uli nvir
        );

        MOBasisRangeBuilder(
            uli nocc_tile,
            uli nvir_tile,
            usi nirrep,
            const uli* ncore_pi,
            const uli* ndocc_pi,
            const uli* nsocc_pi,
            const uli* nvir_pi,
            const uli* ncabs_pi = 0
        );

        ~MOBasisRangeBuilder();

        IndexRange* get_core_range() const;

        IndexRange* get_act_docc_range() const;

        IndexRange* get_occ_range() const;

        IndexRange* get_act_a_occ_range() const;

        IndexRange* get_act_b_occ_range() const;

        IndexRange* get_socc_range() const;

        IndexRange* get_docc_range() const;

        IndexRange* get_orb_range() const;

        IndexRange* get_vir_range() const;

        IndexRange* get_act_a_vir_range() const;

        IndexRange* get_act_b_vir_range() const;

        IndexRange* get_ri_range() const;

        IndexRange* get_cabs_range() const;

        MOSymmetryMapPtr get_symmetry_map() const;

};

class MORangeBuilder : public smartptr::Countable {

    private:
        uli nper_;

        uli ntot_;

    public:
        MORangeBuilder(uli nper, uli ntot);

        IndexRange* get();

};

class MOSymmetryMap :
    public smartptr::Countable
{

    private:
        usi nirrep_;

        usi nspaces_;

        usi* irrep_map_;

        usi symmetry_depth_;

    public:
        MOSymmetryMap(
            IndexRange* core,
            IndexRange* docc,
            IndexRange* socc,
            IndexRange* vir,
            usi nirrep,
            usi symmdepth
        );

};


}


#endif // MOBASIS_H
