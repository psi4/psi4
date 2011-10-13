#ifndef yeti_index_h
#define yeti_index_h

#include <vector>

#include "class.h"
#include "tuple.h"

#include "index.hpp"
#include "permutation.hpp"
#include "sort.hpp"


#include "mallocimpl.h"
#include "yetiobject.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

#define DEFAULT_INDEX_DESCR_ID 0

class IndexDescr :
    public smartptr::Countable
{
    
    private:

        /** A sentence-length description of the index*/
        std::string descr_;

        /** The complete list of tiles at a given metadata depth */
        SubindexTuplePtr range_list_;

        usi depth_;

        uli nelements_total_;

        uli* nelements_per_index_;

        uli* data_start_indices_;

        uli max_nelements_metadata_;

        uli average_nelements_metadata_;

        uli data_offset_;

        uli top_index_offset_;

        uli ntop_index_;

        uli max_nelements_data_;

        uli average_nelements_data_;

        uli total_data_size_;

        uli total_metadata_size_;

        uli nranges_data_;
        
        IndexRangePtr parent_range_;

        IndexDescrPtr subdescr_;

        usi* index_to_irrep_;

        size_t* total_data_sizes_;

        uli descr_id_;

    public:
        /**
         * @param id See #id_
         * @param descr See #descr_
         * @param
         */
        IndexDescr(
            const std::string& descr,
            const IndexRangePtr& range
        );

        virtual ~IndexDescr();

        IndexRange* get_top_range() const;

        IndexRange* get_range(usi depth, uli index) const;

        IndexRange* get_range(usi depth) const;
        
        IndexRange* get_parent_range() const;

        IndexDescr* get_subdescr() const;

        std::string descr() const;

        usi depth() const;

        uli data_index_start() const;

        bool is_equivalent(IndexDescr* descr) const;

        uli nelements_data() const;

        uli nranges_data() const;

        uli start(usi depth) const;

        void print(std::ostream& os = std::cout) const;

        uli nelements(usi depth, uli index) const;

        uli nelements(usi depth) const;

        uli nelements(uli index) const;

        uli max_nelements_metadata() const;

        uli average_nelements_metadata() const;

        uli index_start(uli index) const;

        uli index_increment(uli index) const;

        uli nelements_total() const;

        uli max_nelements_data() const;

        uli average_nelements_data() const;

        usi irrep(uli idx) const;

        size_t total_data_size(uli idx) const;

        bool has_symmetry() const;

        void set_descr_id(uli id);

        uli get_descr_id() const;
};

class TensorIndexDescr :
    public smartptr::Countable
{
    private:
        usi nindex_;

        CountableArray<IndexDescr> indices_;

        bool has_symmetry_;

    public:
        TensorIndexDescr(usi nindex);

        void set(usi index, const IndexDescrPtr& descr);

        usi depth() const;

        TensorIndexDescr* get_subdescr() const;

        IndexDescr* get(usi index) const;

        void permute(Permutation* p);

       void get_nelements(usi depth, const uli* indexset, uli* nelements) const;

       void get_nelements(const uli* indexset, uli* nelements) const;

        /**
            Get the number of data elements per index for a given set of unique data indices.
            @param indexset
            @param nelements
            @param p A permutation to apply to the descr before getting nelements
        */
        void get_nelements(const uli* indexset, uli* nelements, Permutation* p) const;

        void get_index_starts(const uli* indexset, uli* starts) const;

        void get_index_increments(const uli* indexset, uli* increments) const;

        void get_index_starts(uli* starts) const;

        void get_index_increments(uli* starts) const;

        uli nelements(const uli* indices);

        void get_nelements_data(uli* nelements) const;

        size_t get_nelements_data(usi depth, const uli* nelements) const;

        usi nindex() const;

        uli size(const uli* indices, usi depth);

        ulli totalsize() const;

        size_t total_data_size(const uli* indices) const;

        size_t max_nelements_data() const;

        size_t average_nelements_data() const;

        size_t max_nelements_metadata() const;

        size_t average_nelements_metadata() const;

        uli max_nblocks_data() const;

        uli nblocks_tot_at_depth(usi depth) const;

        void print(std::ostream& os = std::cout) const;

        bool has_symmetry() const;

};

/**
 @class IndexRange
 Encapsulates a specific range of tile indices
*/
class IndexRange :
    public smartptr::Countable
{
    
    private:
        
        /** The starting index for the tile*/
        uli start_;

        /** The number of indice in the tile.  This is the
         total number, not just the number of subindices. */
        uli nelements_;

        /**
            The index range number in the parent's list of subranges.
            Calling parent_->get(index_) should return this index range.
        */
        uli index_;

        usi irrep_;

        bool has_symmetry_;

        IndexRange* parent_;

        /** A set of subindices. This is null for bottom level */
        SubindexTuplePtr subranges_;

        void validate();

        /** Given a set of subindex sizes, initialize the
            subindices */
        void init(const std::vector<uli>& subsizes);

        /** Given a total number of indices and the desired number
            of indices per tile, construct the set of sizes.  This
            is really only relevant for routines with irregular n.
            For example, if n is 13 and nper is 5, this would create
            blocks of size 5, 4, 4.
            @param n The total number of indices
            @param nper The requested number of indices per tile
            @param sizes Reference return
        */
        static void
        form_sizes(
            uli n,
            uli nper,
            std::vector<uli>& sizes
        );

        uli nsubranges(usi depth) const;

        void increment_offset(usi depth, uli offset);

        bool _subrange_index_location(IndexRange* range, uli* indices) const;

        IndexRange* _squeeze_together_bottom_ranges(uli nper, uli* start_offsets);

    public:

        /**
         * Create a set of tile index as a composite of subindices
         * @param tuple
         * @parm idx The index number for the tile in its group of indices
         */
        IndexRange(
            uli start,
            const SubindexTuplePtr& tuple
        );

        IndexRange(
            const SubindexTuplePtr& tuple
        );

        /**
            Add a single extra metadata layer to the index range
        */
        IndexRange(
            IndexRange*,
            uli start = 0
        );

        /**
            Pad the index range with extra layers to match the given depth
        */
        IndexRange(
            IndexRange* range,
            usi depth,
            uli start = 0
        );

        /**
           Create an index range by picking out the subranges
           in the tuple.  For example, if start is 1 and nslice
           is 2, this would create an index range with subranges
           tuple[1] + tuple[2].
           @param start     The strarting index of the tuple slice
           @param nslice    The number of indices to include
           @param tuple     The tuple to base subranges on
           @param descr
         */
        IndexRange(
            uli start,
            uli nslice,
            const SubindexTuplePtr& tuple
        );

        /**
         Builds an index range with sliced subranges.  Based on nper,
         a set of subranges are sliced from tuple to create a subindex.
         For example, if tuple had 8 elements and nper was 3, this would
         create an index range with 3 subindices.  Subindex 0 would
         contain tuple[0]->tuple[2] as subindices, subindex 1 would contain
         tuple[3]->tuple[5], subindex 2 would have tuple[6]->tuple[7].
         @param start The starting index for the index range
         @param tuple The set of index ranges to build from
         @param descr
         @param nper  The desired number of indices per subindex
         */
        IndexRange(
            uli start,
            const SubindexTuplePtr& tuple,
            uli nper
        );

        /**
         * Create a new tile index with no subindices
         * @param start The index number start
         * @param n The number of indices in the range
         */
        IndexRange(
            uli start,
            uli n
         );

        /**
         * Create a new tile index with subindices of definite size
         * @param start The index number start
         * @param subsizes The sizes for each subindex
         */
        IndexRange(
            uli start,
            const std::vector<uli>& subsizes
        );

        /**
         * Create a new tile index
         * @param start The index number start
         * @param n The number of indices in the range
         */
        IndexRange(
            uli start,
            uli n,
            uli nper
        );

        IndexRange(
            usi depth,
            uli nper,
            uli start = 0
        );

        virtual ~IndexRange();

        void acquire_subranges(const SubindexTuplePtr& tuple);

        void increment_offsets();

        void increment_offsets(IndexRange* range);

        void set_offsets();

        void set_offsets(IndexRange* range);

        void set_parent(IndexRange* range);

        void get_subranges(std::list<IndexRange*>& sublist, usi depth) const;

        bool has_symmetry() const;

        bool has_subrange(IndexRange* range) const;

        bool contains(IndexRange* range) const;

        bool contains(uli index, usi depth) const;

        bool contains(uli index) const;

        IndexRange* squeeze_together_bottom_ranges(uli nper);

        /**
            @return The depth the given subrange exists at
        */
        usi get_subdepth_alignment(IndexRange* range);

        static IndexRange* get_merged_range(IndexRange* r1, IndexRange* r2);

        uli index() const;

        void set_index(uli index);

        void offset(const IndexRangePtr& index);

        /**
            @return The number of indices in the range
         */
        uli nelements() const;

        /**
            @param The depth to compute the total number of tiles at
            @return The total number of indices in the range at a given depth
        */
        uli ntot(usi depth = 0) const;

        /**
            @param The depth to compute the maximum size at
            @return The maximum range size at a given depth
        */
        uli nmax(usi depth = 0) const;

        /**
        */
        uli start(usi depth) const;

        /**
            @return The index (not inclusive) that defines the end of the index range.
                    Given three subranges (0-3), (6-8), (9-12) finish would be 13.
        */
        uli finish(usi depth) const;

        /**
            @return The index number the range begins on
         */
        uli start() const;

        uli stop() const;

        usi depth() const;

        void expand_subrange_depth(usi maxdepth);

        /**
            @param idx The subindex
            @return The index range defining the subindex idx
         */
        IndexRange* get_subindex(uli idx) const;

        IndexRange* get_parent() const;

        IndexRange* get_first_child() const;

        IndexRange* split_bottom_range() const;

        IndexRange* get_composite_range(usi depth) const;

        /**
            Create an index range in which
        */
        IndexRange* shift_bottom_range() const;

        void
        get_subindices(
            std::list<IndexRange*>& ranges,
            usi depth
        );

        /**
         * @return The set of index ranges for the subindices
         */
        SubindexTuple* get_subranges() const;

        /**
         * Create a subset index range for each index in the range
         * @param The size of each subrange
         */
        void split(uli range);

        /** */
        void print(std::ostream& os = std::cout) const;

        /**
            Whether the two indices are equal
            @param idx
            @return
        */
        bool equals(IndexRange* idx) const;

        /**
            Figure out the number of bottom level index ranges
            with a given size.  Map key is the size and map value
            is the size.
        */
        void sizes(std::map<uli, uli>& sizes) const;

        /**
          * @return Whether the index range has subindices
          */
        bool is_parent() const;

        void subrange_index_location(IndexRange* range, uli* indexset) const;

        bool is_contiguous() const;

        void set_offset(uli offset);

        usi get_irrep() const;

        void set_irrep(usi irrep);

};


class SubindexTuple :
    public CountableArray<IndexRange>
{

    public:
        SubindexTuple(uli n);

        SubindexTuplePtr slice(uli start, uli stop);

        /**
            @param The subindex to find the index number for
            @return The index of a given subindex
        */
        uli index(IndexRange* subidx) const;

};

std::ostream&
operator<<(std::ostream& os, IndexRange* range);

} //end namespace

#ifdef redefine_size_t
#undef size_t
#endif

#endif

