#ifndef yeti_tensor_h
#define yeti_tensor_h

#include <libsmartptr/serialize.h>

#include "tile.h"

#include "data.hpp"
#include "index.hpp"
#include "tensor.hpp"
#include "matrix.hpp"
#include "permutation.hpp"


namespace yeti {

/**
  @class Tensor
  Encapsulates a top-level tile.  This essentially a #MetaDataTile
  with indices (0 0 0 0, etc) that starts the recursive structure of
  a tile.
*/
class Tensor :
    public Tile
{

    private:
        void init();

        void* mem_;

        uli memsize_;

        PermutationGroupPtr pgrp_;

        /**
            The maximum depth at which all index ranges will finally
            have their metadata recursions depths "aligned"
        */
        usi alignment_depth_;

        void finalize_mode();

        std::list<IndexRangeTuplePtr> nonnull_tuples_;

        TileRegistryPtr registry_;

    public:  //tensor configuration options
        typedef enum { alpha_tensor = 0, beta_tensor = 1, gamma_tensor = 2 } contraction_priority_t;

    protected:
         /**
          When computing location, push this tiles index set
          on to specify the tile location at this tile level
          @param loc
         */
        void push_location(const TileLocationPtr& loc) const;

        Tensor(
            const TensorConfigurationPtr& tensor
        );

    public:
        /**
            Create a tile map based on the index ranges provided.
            @param tuple  The index ranges definining the tensor's
                          underlying #MetaDataTile structure
            @param pgrp   The permutation group defining equivalent quantities
        */
        Tensor(
            const std::string& name,
            const IndexRangeTuplePtr& tuple,
            const PermutationGroupPtr& pgrp
        );

        /**
            Create a tile map based on the index ranges provided, but
            link to some parent tensor's configuration
            @param tuple  The index ranges definining the tensor's
                          underlying #MetaDataTile structure
            @param pgrp   The permutation group defining equivalent quantities
        */
        Tensor(
            const std::string& name,
            const IndexRangeTuplePtr& tuple,
            const PermutationGroupPtr& pgrp,
            const TensorConfigurationPtr& config
        );

        /**
            Create a tensor with an extra layer of metadata padded on top.
        */
        Tensor(
            const TensorPtr& tensor
        );

        virtual ~Tensor();

        void accumulate(
            const TensorPtr &tile,
            const PermutationPtr &p,
            double scale
        );

        usi alignment_depth() const;

        void allocate();

        static PermutationGroup*
        antisymmetric_permutation_group(
            const IndexRangeTuplePtr& bratuple,
            const IndexRangeTuplePtr& kettuple
        );

        static Tensor* build(
            const std::string& name,
            usi nidx,
            IndexRange* range
        );

        Tensor* copy(const std::string& name);

        void distribute();

        template <class data_t>
        data_t
        dot_product(const TensorPtr& tensor);

        double dot_product_double(const TensorPtr& tensor);

        float dot_product_float(const TensorPtr& tensor);

        int dot_product_int(const TensorPtr& tensor);

        quad dot_product_quad(const TensorPtr& tensor);

        /**
            Determines whether two tensors have the same structure
            and also that all the data values are equal.
            @param tensor  The tensor to compare to
            @return Whether the tensors are equivalent
        */
        bool equals(const TensorPtr& tensor);

        bool equals(const ThreadedTileElementComputerPtr& filler);

        bool equals(const void* data);

        /**
            Fill the tensor using the pre-defined tensor filler.
        */
        void fill();

        /**
            Assign a new tensor filler and fill.
        */
        void fill(const ThreadedTileElementComputerPtr &filler);

        void fill(const TileElementComputerPtr& val);

        /**
         * @return permutation group for the indices
         */
        PermutationGroup* get_permutation_grp() const;

        const std::string& name() const;

        void print(std::ostream& os = std::cout) const;

        void reconfigure(const TensorConfigurationPtr& config);

        /**
            Inform the t
        */
        void register_nonnull_tiles(const IndexRangeTuplePtr& tuple);

        /**
            Determine the number of blocks of a given size
            at the lowest level. Map key is
            the block size and map value is the count.
            @param sizes The map that will store block sizes
        */
        void sizes(std::map<size_t, size_t>& sizes) const;

        /**
            Sort all data/metadata by the given permutation
            @param p
        */
        void sort(const PermutationPtr& p);

        void set_read_mode();

        void set_write_mode();

        void set_accumulate_mode();

        /**
            @return The total number of elements in the tensor at the data level assuming
                    no non-null structure
        */
        uli total_size() const;

        void* operator new(uli n);

        void operator delete(void* ptr);

        Tensor::contraction_priority_t priority() const;

};

class TensorConfiguration :
    public smartptr::Countable
{

    private:
        friend class Tensor;

        DataMode* data_mode;

        TensorConfiguration();

    public:
        TensorConfiguration(
            const PermutationGroupPtr& grp,
            const std::string& name
        );

        ~TensorConfiguration();

        TensorConfigurationPtr copy(const std::string& name) const;

        DataBlockFactoryPtr data_factory;

        Tile::distribution_t tile_distribution_type;

        TileMap::storage_t map_storage_type;

        usi* distr_indices;

        usi nindex_distr;

        Tensor::contraction_priority_t priority;

        std::string name;

        LayeredDataCachePtr cache;

        ThreadedTileElementComputerPtr filler;

        TileMapBuilderPtr tile_map_builder;

        bool print_data;

        std::list<TileFilterPtr> filters;

        const DataMode* get_data_mode() const;

};

class SpecificTileMapBuilder :
    public TileMapBuilder
{
};


class SubtensorTileMapBuilder :
    public TileMapBuilder
{
    private:
        Tensor* tensor_;

        IndexRangeTuplePtr subtuple_;

        TileRegistryPtr tile_registry_;

        usi subdepth_;

    public:
        SubtensorTileMapBuilder(
            Tensor* tensor,
            const IndexRangeTuplePtr& subtuple
        );

        ~SubtensorTileMapBuilder();

        TileMapPtr build_map(
            const IndexRangeTuplePtr& tuple,
            const TilePtr& parent
        );

};


}

#endif
