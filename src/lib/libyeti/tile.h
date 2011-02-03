#ifndef yeti_tile_h
#define yeti_tile_h

#include <list>

#include "class.h"
#include "index.h"

#include "tile.hpp"
#include "permutation.hpp"
#include "contraction.hpp"
#include "data.hpp"
#include "tensor.hpp"
#include "cache.hpp"
#include "sort.hpp"
#include "matrix.hpp"
#include "elementop.hpp"
#include "filler.hpp"

#include "mallocimpl.h"
#include "yetiobject.h"


#define DISTRIBUTION_DEPTH 1

namespace yeti {

/**
    @struct tile_estimate_t
    Struct for holding estimates of how large a tile's maximum element is
*/
struct tile_estimate_t {

    float maxlog;

    usi depth;

};

/**
    @class TileDistributer
    For a given set of tiles, distributes the storage of the tiles
    amongst all processors
*/
class TileDistributer :
    public smartptr::Serializable
{

    private:

        int* sorted_task_list_;

        usi* distr_indices_;

        usi ndistr_index_;

        uli ntasks_;

        uli nworkers_;

    public:
        /**
            @param tuple
            @param dist_indices The subset of indices to be used in distribution
            @param ndist_indices
        */
        TileDistributer(
            const usi* dist_indices,
            usi ndist_indices
        );

        ~TileDistributer();

        /**
            Register a tile within a tensor for parallel distribution
            @param tile
        */
        void tally(const TilePtr& tile);

        /**
            @param tile
            @return The process number in a parallel environment
                    that should compute tasks involving this tile
        */
        uli process_number(const TilePtr& tile) const;

        /**
            This clears the tally.  This should only be called after all
            tiles have been tallied.
        */
        void sort_tasks();

};

class TileFilter :
    public smartptr::Countable
{
    public:
        virtual bool is_zero(
            const uli* indices,
            usi depth
        ) const = 0;

        virtual TileFilter* copy() const = 0;

};

class PermutationalSymmetryFilter :
    public TileFilter
{
    private:
        PermutationGroupPtr pgrp_;

    public:
        PermutationalSymmetryFilter(const PermutationGroupPtr& grp);

        bool is_zero(const uli* indices, usi depth) const;

        TileFilter* copy() const;
};

/**
  @class TileMap
  Encapsulates a mapping from a composite index to a tile
*/
class TileMap :
    public smartptr::Serializable,
    public Malloc<TileMap>
{
    public:
        typedef Tile** iterator;

        /** Whether the tile map should be generated dynamically
            and cleared after every use or whether the the tile
            map should be permanently stored */
        typedef enum { permanent, temporary } storage_t;
        
    private:

        IndexerPtr indexer_;

        CountableArray<Tile>* tiles_;

        char* rigorously_zero_tiles_;

        IndexRangeTuplePtr tuple_;

        PermutationGroupPtr pgrp_;

        ThreadedTileElementComputerPtr filler_;

        /**
            Private variable used in recursive loop for
            generating all possible indices
        */
        uli* tmpindices_;

        /**
            Private variable used in recursive loop for
            generating all possible IndexRange subtuples
        */
        IndexRange** tmptuple_;

        /**
        */
        uli* image_;

        /**
            The index range offsets.  Not all indices are tiled since
            recursion depths may be mismatched.  For untiled
            indices, the offset is always zero.
        */
        uli* offsets_;

        /**
            The index range sizes.  Not all indices are tiled since
            recursion depths may be mismatched.  For untiled
            indices, the index range size is always 1.
        */
        uli* nindices_;

        uli* index_strides_;

        usi* is_tiled_index_;

        Tile** parents_;

        usi nparents_;

        /**
            The cursion depth of the highest level index range.
            Only index ranges at this level
            will be involved in the tiling scheme
        */
        usi maxdepth_;

        usi mindepth_;

        /**
        */
        bool computed_;

        bool allocated_;

        /**
            The number of tile indices
        */
        uli nindex_;

        /**
         * Recursion method
         * @param index The current index position in the recursion
         * @param idxmap The map being created
         */
        void iterate(usi index);

        void insert_new_tile();
        /**
            Constructor used for making copies of a TileMap
        */
        TileMap();

    public:

        /**
            Build an empty tile map with indexing based
            on the provided index range tuple
            @param tuple
            @param pgrp
            @param parent
        */
        TileMap(
            const IndexRangeTuplePtr& tuple,
            const PermutationGroupPtr& pgrp,
            const TilePtr& parent
         );

        /**
            Constructor for creating a tile map with a single element
            @param The single tile in the tile map
            @param parent
        */
        TileMap(
            Tile* tile,
            Tile* parent
         );


        ~TileMap();

        void allocate();

        bool is_allocated() const;

        /**
            @return The number of elements tiled per index. For indices
            that match the max metadata depth, this will be the number
            of indices in an index range.  For indices which are
            at a lower metadata depth, this will just be 1 since tiling
            only occurs on nindices with the correct metadata depth
        */
        const uli* nindices() const;

        /**
            @return The offsets at which each index range begins
        */
        const uli* offsets() const;

        const uli* index_strides() const;

        const usi* is_tiled_index() const;

        bool depths_aligned() const;

        /**
          Compute a composite index, map it, and return tile
          @param indexset
          @return Tile specified by indexset.  May be NULL.
         */
        TilePtr get(const constIndexSetPtr& indexset) const;

        /**
            Map the index and return associated tile
            @param index
            @return Tile specified by index. May be NULL.
        */
        Tile* get(uli index) const;

        /**
            Compute a composite index and return associate tile
            @param indices
            @return Tile specified by indices
        */
        Tile* get(const uli* indices) const;

        /**
            Check to see if the index number is at the correct depth to be tiled.
            If yes, return the appropriate subindex. If not, return the parent index range.
            @param index_number
            @param index
        */
        IndexRange* get_subrange(usi index_number, uli index) const;

        /**
            Given a composite index, compute the individual indices and place
            them in the parameter array.
        */
        void compute_indices(uli comp_index, uli* index_arr);

        IndexRangeTuplePtr get_index_ranges() const;

        /**
            Create a copy of the current tile map linked to a new parent tile
            @param parent The parent tile for the copied map
            @return Tile map copy
        */
        TileMap* copy(Tile* parent);

        Tile* get_parent_tile() const;

        TensorConfiguration* config() const;

        /**
            @param A set of indices specifying a particular tile
            @return Whether the tile exists
        */
        bool exists(const uli* indices) const;

        void retrieve();

        void release();

        void insert(uli index, Tile* tile);

        bool is_rigorously_zero(uli index) const;

        bool is_rigorously_zero(Tile** tile) const;

        uli index(const uli* indices) const;

        /**
            Create a composite index, map the index, and append the
            new tile to the map.
            @param tile
        */
        void insert(Tile* tile);

        void print(std::ostream& os = std::cout) const;

        /**
         * @return The number of non-null tiles
         */
        uli ntiles() const;

        /**
            @return The maximum possible number of tiles
        */
        uli size() const;

        /**
            @param nindex The index number
            @return The number of indices tiled for this index
        */
        uli nrange(usi nindex) const;

        /**
            @return The maximum recursion depth of all index ranges tiled here
        */
        usi maxdepth() const;

        iterator begin() const;

        iterator end() const;

        /**
            @param sort
            @param buffer Workspace buffer for temporarily holding sorted values
        */
        void sort(
            const SortPtr& sort,
            void* buffer
        );

        void set_max_log(float log, usi depth);

        void add_parent(Tile* tile);

        void remove_parent(Tile* tile);
        
};

class TileIterator : public smartptr::Countable {

    private:
        TileMap::iterator* iters_;

        TileMap::iterator* stops_;

        TileMap::iterator stop_;

        TileMap::iterator iter_;

        PermutationPtr perm_;

        void* buffer_;

        Tile* tile_;

        Tile** link_;

        TilePtr parent_;

        size_t me_;

        usi maxdepth_;

        int depth_;

        bool done_;

        bool local_;

    public:
        TileIterator(
            Tile** link,
            const TilePtr& parent,
            bool local,
            size_t me,
            const PermutationPtr& p = 0
        );

        ~TileIterator();

        void start();

        void next();

        bool done();

        Tile* tile() const;

};

class TileIteratorWorkspace {

    public:
        TileIteratorWorkspace();

        ~TileIteratorWorkspace();

        char** buffers;

};


/**
    @class Tile
    Parent class for tiles
*/
class Tile :
    public Malloc<Tile>,
    public YetiRuntimeSerializable
{

    public:
        typedef enum { meta_data = 0, data = 1 } tile_t;

        typedef enum { replicated = 0, distributed = 1 } distribution_t;

    protected:
        template <class T>
        void
        data_action(const T& obj);

        /** The set of index ranges for each index in the tile */
        IndexRangeTuplePtr ranges_;

        /**
            The index set defining which tile block this covers
        */
        uli* indices_;

        /** The parent tile.  Might be NULL. */
        TileMap* parent_;

        /**
            The top level tensor that owns this tile
        */
        TensorConfigurationPtr config_;

        TileMapPtr tilemap_;

        DataBlockPtr data_;

        /** Whether or not this tile belongs to a specific processor */
        bool has_owner_process_;

        bool computed_;

        /** If the tile has an owner process, this is the process number */
        uli owner_process_;

        /** An estimate of the maximum data value in the tensor */
        tile_estimate_t estimate_;

        bool _equals(const char** data);

        /**
            Used in making tile copies. Creates tile
            with empty values that will be filled in later.
        */
        Tile(const TensorConfigurationPtr& config);

        /**
            Called from tensor constructor
            @param ranges
            @param indices
        */
         Tile(
            const IndexRangeTuplePtr& ranges,
            size_t* indices
         );

        /**
            Transfer copied valus to empty copy tile
            @param tile
        */
        void copy_to(Tile* tile);

        /**
            Set up tile map or data block depending on whether
            this is a data or metadata tile
        */
        void init();

        void init_estimate();

        void _release(uli threadnum);

        void _retrieve(uli threadnum);


    public:

        /**
         * Constructor
         * @param ranges
         * @param indices
         * @param pgrp
         */
        Tile(
            const IndexRangeTuplePtr& ranges,
            size_t* indices,
            //const PermutationGroupPtr& pgrp,
            const TensorConfigurationPtr& config
         );

        /**
         * Constructor
         * @param ranges
         * @param indices
         * @param pgrp
         */
        Tile(
            const IndexRangeTuplePtr& ranges,
            size_t* indices,
            //const PermutationGroupPtr& pgrp,
            TileMap* parent
         );


        /**
            Debug constructor
            @param indices
            @param sizes
            @param nidx
        */
        Tile(
            TensorConfiguration* config,
            const size_t* indices,
            const size_t* sizes,
            usi nidx
        );

        virtual ~Tile();

        /**
            Accumulate tile values, potentially following a sort and scale
            @param tile
            @param sort
            @param scale
        */
        void accumulate(
            const TilePtr& tile,
            const SortPtr& sort,
            double scale = 1
        );

        /**
            Create all metadata, compute the total data size of the tensor,
            create a memory pool of the appropriate size, and distribute
            memory blocks to all data tiles
        */
        void allocate();

        /**
            @return The total size in bytes of all data tiles owned by
            the current processor
        */
        uli total_data_block_size();

        /**
            @return The top level tensor which contains all the configuration information
        */
        TensorConfiguration* config() const;

        /**
            Create a copy of the current tile, but linked to a new parent tile
            @param parent
            @return Copy tile
        */
        Tile* copy(TileMap* parent);

        /**
            @return The maximum recursion depth of all index ranges tiled here
        */
        usi depth() const;

        /**
            Distribute all subtiles on the given tile distribution object
            @param distr
        */
        void distribute(const TileDistributerPtr& distr);

        /**
            Loop all data tiles and perform an element operation in place
            on the data
            @param op
        */
        void element_op(const ElementOpPtr& op);

        /**
            @param data
            @return Whether
        */
        bool equals(const void* data);

        /**
            Figure out whether a tile object exists with the given indices
            @param indices
        */
        bool exists(const uli* indices) const;

        /**
            @param filler
            @return Whether the tile values match those generated by the filler
        */
        bool equals(const ThreadedTileElementComputerPtr& filler);

        /**
            Create all data blocks, but do not allocate or fill. Perform
            fill operation but with "null" filler
        */
        void fill();

        /**
            @param indices The indices specifying a given tile location
            @return The tile at the given index location. May be NULL.
        */
        Tile* get(const uli* indices) const;

        /**
            @return The data block for this tile. This will be NULL
                    for metadata tiles.
        */
        DataBlock* get_data() const;

        /**
         * @return The index ranges specifying the tile map for this tile
         */
        IndexRangeTuple* get_index_ranges() const;

        /** Get the set of indices defining the tile locatio
         * @return index set
         */
        const uli* indices() const;

        bool is_allocated() const;

        /**
         * Recursively walk through parents to get exact location
         * @return
         */
        TileLocationPtr get_location() const;

        /**
         * Return the map of subtiles. May be NULL
         * @return tile map
         */
        TileMap* get_map() const;

        const tile_estimate_t& get_max_log() const;

        /**
            @return Whether the given tile is stored on a given processor
        */
        bool has_owner_process() const;

        /**
            @return The number of indices tiled on each index.  See TileMap for details
                    as whether a given index range is tiled depends on recursion depth.
        */
        const uli* index_sizes() const;

        /**
            @return The offsets for the index ranges tiled on each index.  See TileMap for details
                    as whether a given index range is tiled depends on recursion depth.
        */
        const uli* index_offsets() const;

        bool is_aligned() const;

        /**
            Determines whether two tensors have the same structure.  This
            checks only the MetaDataTile structure, and does not compare
            values.  Tensors built from the same index ranges would be
            equivalent.
            @param tensor  The tensor to compare to
            @return Whether the tensors are equivalent
        */
        bool is_equivalent(const TilePtr& tile) const;

        /**
            @return Whether the tile has metadata subblocks
        */
        bool is_grand_parent() const;

        /**
            @return Whether the tile has subblocks, data or metadata
        */
        bool is_parent() const;

        /**
            @return Whether the tile has a parent
        */
        bool is_top_tile() const;

        /**
            Determine whether a given set of indices are permutationally unique
            based on metadata recursion depth and permutational symmetry
        */
        static bool is_unique(
            const size_t* indices,
            const PermutationGroupPtr& grp,
            usi depth
        );

        /**
            @return The log of max absolute value
        */
        float max_log() const;

        /**
            @return The maximum number of elements in either a data block or metadata map
        */
        uli max_nblock() const;

        /**
            @return The maximum number of bytes in either a data block or metadata map block
        */
        uli max_blocksize() const;

        /**
            @return The data size for this block
        */
        uli ndata() const;

        /**
            @return The number of indices being tiled
        */
        usi nindex() const;

        /**
            @param The index number
            @return The number of elements in a given index.  Depending on alignment
                     of recursion depths, this will be 1 if index is tiled or other value.
        */
        uli nrange(usi nindex) const;

        /** Determine the number of indices at a given index
            at a given recursion depth.
            @param depth Desired recursion depth
            @return Number of indices in range
            @throw SanityCheckError if depth is greater than maximum allowed
        */
        uli nrange(usi index, usi depth);

        /** Determine the number of tiles at a given level. The level
            given here correspond to TileMap::maxlevel_. Recursion depth
            starts with 0 being the deepest level.
            @param depth Desired recursion depth
            @return Number of tiles at a given recursion depth.
            @throw SanityCheckError if depth is greater than maximum allowed
        */
        uli ntiles(usi depth);

        /**
            @return The number of tiles owned on this processor by this tile
        */
        uli ntiles_owned();

        /** Determine the maximum number of tiles possible this might map to
        */
        uli ntiles_max() const;

        /** Determine the number of tiles owned by this process. See ntiles
            for more details on recursion depth.
            @param depth Desired recursion depth
            @return Number of tiles at a given recursion depth
            @throw SanityCheckError if depth is greater than maximum allowed
        */
        uli ntiles_owned(usi depth);

        virtual void print(std::ostream& os = std::cout) const;

        /**
            @return The process number that owns this tile
        */
        uli owner_process() const;

        void register_tiles(
            const TileRegistryPtr& reg,
            usi depth
        );

        /**
            @param sort
            @param buffer Buffer for holding temporary sorted values
        */
        void _sort(
            const SortPtr& sort,
            void* buffer
        );

        /**
            Permutate all metadata structures and data values
            @param p
        */
        void sort(const PermutationPtr& p);

        /**
            Scale all values by a given factor
            @param factor
        */
        void scale(int factor);

        /**
            Scale all values by a given factor
            @param factor
        */
        void scale(double factor);

        /**
            Update the tile to let it know that tile values have been scaled.
            This is called automatically by scale.
            @param log_factor  The log of the scale factor
        */
        void scale_max_log(float log_factor);

        void set_max_log(const tile_estimate_t& est);

        /**
            Set the maximum log of the tile, informing the tile
            at what recursion depth the estimate came from
        */
        void set_max_log(float max, usi depth);

        /**
            Set the process number in a parallel environment that will own this tile
        */
        void set_owner_process(uli owner);

        /**
            Register this tile with a tile distributer to configure a tile distribution
            @param distr
        */
        void tally(const TileDistributerPtr& distr);

        /**
            Perform a variety of cleanup/compression operations.
            Currently this just loops the tensor and redetermines the
            maximum element.
        */
        void update();

        tile_t type() const;

};

class TileMapBuilder :
    public smartptr::Countable
{

    public:
        virtual TileMapPtr build_map(
            const IndexRangeTuplePtr& tuple,
            const TilePtr& parent
        ) = 0;

};

class DefaultTileMapBuilder :
    public TileMapBuilder
{
    private:
        PermutationGroupPtr pgrp_;

    public:
       DefaultTileMapBuilder(const PermutationGroupPtr& grp);

       TileMapPtr build_map(
            const IndexRangeTuplePtr& tuple,
            const TilePtr& parent
        );
};

class TileRegistry :
    public smartptr::Countable
{
    private:
        typedef std::map<IndexRangeLocationPtr, Tile**, IndexRangeLocationCompare> node_map;

        node_map tilemap_;

        void delete_tile(Tile** tileptr);

    public:
        TileRegistry();

        void register_tile(Tile** tile);

        Tile** get_tile(const IndexRangeTuplePtr& tuple);

        bool has_tile(const IndexRangeTuplePtr& tuple);

        void remove_tile(const IndexRangeTuplePtr& tuple);

        void print(std::ostream& os = std::cout) const;

        void delete_tiles();

        uli size() const;

};


}

//YETI_RUNTIME_HEADER(Tile)

#endif
