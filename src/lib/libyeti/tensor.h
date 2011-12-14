#ifndef yeti_tensor_h
#define yeti_tensor_h

#include "class.h"
#include "mapimpl.h"
#include "mallocimpl.h"
#include "messenger.h"
#include "contraction.h"

#include "data.hpp"
#include "index.hpp"
#include "tensor.hpp"
#include "tensor.hpp"
#include "tensorparser.hpp"
#include "tensorblock.hpp"
#include "node.hpp"
#include "permutation.hpp"
#include "tensoraction.hpp"
#include "filler.hpp"
#include "elementop.hpp"
#include "matrix.hpp"
#include "contraction.hpp"


#include "cache.h"

#include "gigmatrix.h"

#include <list>

#ifdef redefine_size_t
#define size_t custom_size_t
#endif



namespace yeti {

class Tensor :
    public Malloc<Tensor>,
    public YetiRuntimeCountable
{

    public:
        typedef enum {
            in_core = 0,
            on_disk = 1,
            recomputed = 2,
            action = 3,
            default_storage = 4
        } tensor_storage_t;

        typedef enum {
            gamma_tensor = 0,
            beta_tensor = 1,
            alpha_tensor = 2,
            override_priority = 3
        } tensor_priority_t;

    private:
        uli sort_weight_;

        uli tensor_number_;

        std::string name_;

        bool is_open_subtensor_;

        TensorBlockMap* tensor_blocks_;

        Indexer* main_indexer_;

        std::list<TensorElementFilterPtr> filters_;

        PermutationGroupPtr tensor_grp_;

        PermutationGroupPtr original_grp_;

        DataCachePtr data_cache_;

        DataCachePtr metadata_cache_;

        PermutationSetPtr tensor_grp_generator_set_;

        TensorIndexDescrPtr block_descr_;

        TensorIndexDescrPtr parent_descr_;

        MultiThreadLock cache_lock_;

        usi depth_;
        
        uli index_end_[NINDEX];
        
        uli index_start_[NINDEX];

        bool is_distributed_;

        tensor_storage_t storage_type_;

        tensor_priority_t priority_;
        
        Tensor* parent_tensor_;

        size_t data_block_size_;

        size_t max_data_node_size_;

        size_t metadata_block_size_;

        uli nelements_metadata_av_;

        uli malloc_number_;

        uli nsubtensors_;

        void increment_subtensor_count();

        void decrement_subtensor_count();

        Contraction::tensor_position_t cxn_position_;

        void accumulate(
            TensorBlock* srcblock,
            Tensor* tensor,
            double scale,
            const SortPtr& sort,
            uli block_index
        );

        /**
            Assign a new tensor filler and fill.
        */
        void fill(const ThreadedTensorElementComputerPtr& filler);

        void init_symmetry_filter();

        void print_block(TensorBlock* block, std::ostream& os) const;
        
        ThreadedTensorElementComputerPtr filler_;

        void accumulate_post_process();

        void fill_blocks();

        TensorBlock* get_first_nonnull_block() const;

        uli* block_to_node_indexmap_;

    public:
        typedef TensorBlockMap::iterator iterator;

        Tensor(
            const std::string& name,
            const TensorIndexDescrPtr& descr,
            const PermutationGroupPtr& grp
        );

        Tensor(
            const TensorIndexDescrPtr& descr,
            Tensor* parent_tensor
        );

        ~Tensor();

        static void configure_cache(
            TensorIndexDescr* descr,
            uli& nblocks_metadata,
            size_t& metadata_block_size,
            uli& nblocks_data,
            size_t& data_block_size,
            size_t& max_data_node_size
        );

        void accumulate(
            uli threadnum,
            Tensor* tensor,
            double scale,
            const SortPtr& sort
        );

        void accumulate_no_barrier(
            Tensor* tensor,
            double scale,
            Permutation* p = 0
        );

        void accumulate(
            Tensor* tensor,
            double scale,
            Permutation* p = 0
        );

        /**
            A subtensor which has a subblock which is
            not permutationally unique but the unique parent
            is not contained in the subtensor
        */
        bool is_open_subtensor() const;

        uli get_subtensor_count() const;

        iterator begin() const;

        iterator end() const;

        void reset();

        /**
            Whether these tensors are the same or have the same parent
        */
        bool intersects(Tensor* tensor) const;

        bool contains(const uli* indices) const;

        void configure(BlockRetrieveAction* action);

        void configure(const TensorElementComputerPtr& filler);

        void configure(tensor_storage_t storage_type);

        void configure_degeneracy(const PermutationGroupPtr& pgrp);

        void distribute(const usi* indices, usi nindex);

        TensorConfiguration* config() const;

        void recompute_permutation();

        void get_matrix(RectMatrixPtr& matrix, MatrixConfiguration* config);

        void get_matrix(SymmMatrixPtr& matrix, MatrixConfiguration* config);

        void convert(const SymmMatrixPtr& matrix, MatrixConfiguration* config);

        void convert(const RectMatrixPtr& matrix, MatrixConfiguration* config);

        void convert(Matrix* matrix, MatrixConfiguration* config);

        void accumulate(const SymmMatrixPtr& matrix, MatrixConfiguration* config);

        void accumulate(const RectMatrixPtr& matrix, MatrixConfiguration* config);

        void accumulate(Matrix* matrix, MatrixConfiguration* config);

        void lowest_valid_indices(uli* indices);

        void assign(Tensor* tensor);

        void fill(const double* data);

        void unpack();

        /**
            This cannot be declared const since retrieve/release
            calls may modify the tensor.
        */
        bool equals(const void* data);

        bool equals(Tensor* tensor);

        void element_op(
            uli threadnum,
            ElementOp *op
        );

        void element_op(ElementOp* op);

        void fill(
            uli threadnum,
            const ThreadedTensorElementComputerPtr& filler
        );

        void fill(const MatrixPtr& matrix);

        void accumulate(const MatrixPtr& matrix, MatrixConfiguration* config);

        void fill(const TensorElementComputerPtr& val = 0);

        void global_sum();

        PermutationGroup* get_tensor_grp() const;

        PermutationGroup* get_original_grp() const;

        PermutationGroup* get_matrix_grp() const;
        
        const std::string& get_name() const;

        DataCache* get_data_cache() const;

        DataCache* get_metadata_cache() const;

        tensor_storage_t get_storage_type() const;

        bool is_distributed() const;

        bool is_in_core() const;

        bool is_recomputed() const;

        bool is_replicated() const;

        bool has_retrieved_block() const;

        TensorBlockMap* get_block_map() const;

        ThreadedTensorElementComputer* get_element_computer() const;

        TensorBlock* get_block(
            const uli* indices,
            Permutation* p = 0
        ) const;

        TensorBlock* get_block(
            uli index
        ) const;

        TensorIndexDescr* get_block_descr() const;

        TensorIndexDescr* get_parent_descr() const;

        usi get_depth() const;

        uli get_index(const uli* indices) const;

        uli get_unique_id(
            const uli* indices,
            Permutation* fetch_perm
        ) const;

        Contraction::tensor_position_t get_tensor_position() const;

        Tensor::tensor_priority_t get_priority() const;

        ulli get_totalsize() const;

        TensorBlock* get_make_block(uli idx);

        void insert_new_block(TensorBlock* block);

        bool is_closed_tensor() const;

        bool is_subtensor() const;

        void metadata_sort(const SortPtr& p);

        void internal_contraction(Tensor* dst_tensor, MatrixConfiguration* config);

        TensorBlock* make_block(const uli* indexset, Tensor::tensor_storage_t type = Tensor::default_storage);

        ulli nelements() const;

        ulli nelements_unique() const;

        uli nelements_metadata_av() const;

        uli nblocks_retrieved() const;

        void obsolete();

        bool nonzero() const;

        bool is_parent_block(TensorBlock* block) const;

        double norm();

        void print(std::ostream& os = std::cout);

        void print_blocks(std::ostream& os = std::cout);

        void scale(double scale);

        void set_tensor_position(Contraction::tensor_position_t position);

        void set_priority(Tensor::tensor_priority_t priority);

        void set_sort_size_weight(uli weight);

        uli get_sort_size_weight() const;

        void sort(uli threadnum, const SortPtr& sort);

        void sort(Permutation* p);

        void sync();

        void flush();
        
        /**
            If updates have been performed on a subtensor portion of the given
            subtensor, some blocks which have permutational symmetry may not be
            synced when the subtensor is synced.  This ensures that all permutationally
            related blocks are sync when the subtensor is synced.
        */
        void update_remote_blocks();

        void remote_wait();

        void update_from_subtensor();

        void update();

        void update(uli threadnum);

        void update_to_subtensor();

        void reset_degeneracy();

        bool is_cache_coherent() const;

        Permutation* get_fetch_permutation(const uli* indices) const;

        void zero();

        size_t data_block_size() const;

        size_t max_data_node_size() const;

        size_t metadata_block_size() const;

        Tensor* copy();

        Tensor* clone();

};


} //end namespace

#ifdef redefine_size_t
#undef size_t
#endif

#endif
