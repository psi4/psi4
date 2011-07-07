#ifndef yeti_tensor_h
#define yeti_tensor_h

#include "class.h"
#include "mapimpl.h"
#include "mallocimpl.h"

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
    public YetiRuntimeCountable,
    public Malloc<Tensor>
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
            replicated = 0,
            distributed = 1
        } tensor_distribution_t;

        typedef enum {
            gamma_tensor = 0,
            beta_tensor = 1,
            alpha_tensor = 2
        } tensor_priority_t;

        typedef enum {
            keep_zero_blocks = 0,
            delete_zero_blocks = 1
        } tensor_erase_policy_t;

    private:
        std::string name_;

        BlockMap<TensorBlock>* blocks_;

        Indexer* main_indexer_;

        std::list<TensorElementFilterPtr> filters_;

        DiskBufferPtr disk_buffer_;

        PermutationGroupPtr tensor_grp_;

        PermutationGroupPtr original_grp_;

        DataCachePtr data_cache_;

        DataCachePtr metadata_cache_;

        PermutationSetPtr tensor_grp_generator_set_;

        TensorIndexDescrPtr descr_;

        usi depth_;
        
        uli index_end_[NINDEX];
        
        uli index_start_[NINDEX];

        tensor_storage_t storage_type_;

        tensor_distribution_t distribution_type_;

        tensor_erase_policy_t erase_type_;

        tensor_priority_t priority_;
        
        Tensor* parent_tensor_;

        Permutation* sort_perm_;

        Permutation* unsort_perm_;

        size_t data_block_size_;

        size_t metadata_block_size_;

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

        void erase_zero_blocks();

        void erase_if_null(TensorBlock* block, uli idx);
        
        ThreadedTensorElementComputerPtr filler_;

        void accumulate_post_process();

        void fill_blocks(Tensor::tensor_storage_t);

        TensorBlock* get_first_nonnull_block() const;

    public:
        typedef NodeMap<TensorBlock>::iterator iterator;

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

        static bool in_destructor;
        
        static bool in_sync;

        void accumulate(
            uli threadnum,
            Tensor* tensor,
            double scale,
            const SortPtr& sort
        );

        void accumulate(
            Tensor* tensor,
            double scale,
            Permutation* p = 0
        );

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

        void configure(tensor_distribution_t dist_type);

        void configure(tensor_erase_policy_t erase_type);

        void configure_degeneracy(const PermutationGroupPtr& pgrp);

        TensorConfiguration* config() const;

        void get_matrix(RectMatrixPtr& matrix, MatrixConfiguration* config);

        void get_matrix(SymmMatrixPtr& matrix, MatrixConfiguration* config);

        void convert(const SymmMatrixPtr& matrix, MatrixConfiguration* config);

        void convert(const RectMatrixPtr& matrix, MatrixConfiguration* config);

        void convert(Matrix* matrix, MatrixConfiguration* config);

        void accumulate(const SymmMatrixPtr& matrix, MatrixConfiguration* config);

        void accumulate(const RectMatrixPtr& matrix, MatrixConfiguration* config);

        void accumulate(Matrix* matrix, MatrixConfiguration* config);

        void fill(const double* data);

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

        DiskBuffer* get_disk_buffer() const;

        Permutation* get_unsort_permutation() const;

        Permutation* get_sort_permutation() const;

        PermutationGroup* get_tensor_grp() const;

        PermutationGroup* get_original_grp() const;

        PermutationGroup* get_matrix_grp() const;
        
        const std::string& get_name() const;

        DataCache* get_data_cache() const;

        DataCache* get_metadata_cache() const;

        tensor_storage_t get_storage_type() const;

        tensor_distribution_t get_distribution_type() const;

        tensor_erase_policy_t get_erase_type() const;

        TensorBlockMap* get_block_map() const;

        ThreadedTensorElementComputer* get_element_computer() const;

        TensorBlock* get_block(
            const uli* indices,
            Permutation* p = 0
        ) const;

        TensorBlock* get_block(
            uli index
        ) const;

        TensorIndexDescr* get_descr() const;

        usi get_depth() const;

        uli get_index(const uli* indices) const;

        uli get_unique_id(
            const uli* indices,
            Permutation* fetch_perm
        ) const;

        Tensor::tensor_priority_t get_priority() const;

        ulli get_totalsize() const;

        TensorBlock* get_make_block(uli idx);

        void insert_new_block(TensorBlock* block);

        void insert_new_block(uli idx, TensorBlock* block);

        bool is_closed_tensor() const;

        bool is_subtensor() const;

        void metadata_sort(const SortPtr& p);

        void internal_contraction(Tensor* dst_tensor, MatrixConfiguration* config);

        TensorBlock* make_block(const uli* indexset, Tensor::tensor_storage_t type = Tensor::default_storage);

        ulli nelements() const;

        ulli nelements_unique() const;

        uli nblocks_retrieved() const;


        bool nonzero() const;

        bool unique_nonzero() const;

        bool is_parent_block(TensorBlock* block) const;

        double norm();

        void print(std::ostream& os = std::cout);

        void scale(double scale);

        void set_priority(Tensor::tensor_priority_t priority);

        void set_read_mode();

        void set_write_mode();

        void set_accumulate_mode();

        void sort(uli threadnum, const SortPtr& sort);

        void sort(Permutation* p);

        void sync();
        
        void sync_to_subtensor();

        void update();

        void update(uli threadnum);

        void update_to_subtensor();

        void reset_degeneracy();

        bool is_cache_coherent() const;

        Permutation* get_fetch_permutation(const uli* indices) const;

        bool is_unique(const uli* indices) const;

        void zero();

        size_t data_block_size() const;

        size_t metadata_block_size() const;

        Tensor* copy();

};


} //end namespace

#ifdef redefine_size_t
#undef size_t
#endif

#endif
