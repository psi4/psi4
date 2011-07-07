#ifndef yeti_tensorblock_h
#define yeti_tensorblock_h

#include "class.h"
#include "mapimpl.h"
#include "tensor.h"

#include "data.hpp"
#include "index.hpp"
#include "node.hpp"
#include "tensorbranch.hpp"
#include "tensorblock.hpp"
#include "permutation.hpp"
#include "tensorparser.hpp"
#include "filler.hpp"
#include "elementop.hpp"
#include "matrix.hpp"
#include "contraction.hpp"
#include "tensoraction.hpp"

#include "mapimpl.h"
#include "cache.h"
#include "taskqueue.h"

#include "gigmatrix.h"

#include <list>

#ifdef redefine_size_t
#define size_t custom_size_t
#endif


#define DATA_BLOCK_SIZE_IN_CORE_CUTOFF 50000

namespace yeti {


struct DataStorageNode :
    public Malloc<DataStorageNode>
{

    StorageBlock* block;

    DataStorageNode* next;

};

class DataStorageBlockAllocator :
    public Cachable
{
    protected:
        uli nstorage_blocks_;

        DataStorageNode* first_node_;

        DataStorageNode* last_node_;

    public:
        DataStorageBlockAllocator();

        DataStorageBlockAllocator(YetiRuntimeObject::thread_safety_flag_t flag);

        ~DataStorageBlockAllocator();

        DataStorageNode* first_data_node() const;

        CachedStorageBlock* allocate_cache_block(DataCache* cache);

        InCoreBlock* allocate_core_block(size_t size);

        LocalDiskBlock* allocate_disk_block(
            DiskBuffer* buffer,
            Cachable* item,
            DataCache* cache
        );

        void store(StorageBlock* block);

        uli get_num_data_storage_blocks() const;

        void clear_storage();
        
        virtual void flush_from_cache();
};

class TensorController :
    public smartptr::Countable,
    public Malloc<TensorController>
{
    protected:
        TensorBlock* parent_block_;

    public:
        TensorController(TensorBlock* block);

        virtual void retrieve(TensorBranch* branch) = 0;

        virtual void clear() = 0;

        virtual void flush(TensorBranch* branch) = 0;
        
        virtual void release(TensorBranch* branch) = 0;

        virtual void retrieve(MetaDataNode* mdnode) = 0;

        virtual bool need_recompute_data() = 0;

        virtual void retrieve(
            DataNode* dnode,
            const MetaDataNode* parent,
            const uli* indices
        ) = 0;

        virtual void flush_data() = 0;

        virtual void validate() = 0;

        virtual void obsolete() = 0;

        virtual void sync() = 0;

};

class TensorBlock :
    public DataStorageBlockAllocator,
    public Malloc<TensorBlock>,
    public TaskOwner
{
    private:
        uli block_number_;

        /** Multiple cache entries can be linked to the same tensor block.
            This ensures that these multiple entries do not both attempt
            to flush the same block. */
        bool is_flushed_;

        bool is_synced_;

        void init_mempool();

        void init_in_core();

        void init_on_disk();

        void init_recomputed();

        void init_remote();

        void init_action();

        void init_in_core_no_sort();

        void init_in_core_resort();

        void init_on_disk_no_sort();

        void init_on_disk_resort();

        void init_recomputed_no_sort();

        void init_recomputed_resort();

        void init_action_no_sort();

        void init_action_resort();

        void init_remote_no_sort();

        void init_remote_resort();

        void set_synced();

        void set_unsynced();

        friend class TensorBranch;

        /*********************************************
        This information is replicated across all nodes
        *********************************************/

        /**
            The permutation which maps the node indices
            onto the lowest possible "indexing" in the
            full permutational symmetry of the parent tensor
        */
        Permutation* resort_perm_;

        Permutation* fetch_perm_;

        /**
            The permutation of the block relative to its original
            unpermuted form
        */
        Permutation* block_perm_;

        TensorIndexDescr* descr_;

        TensorRetrieveActionPtr action_;

        Tensor* parent_tensor_;

        uli indices_[NINDEX];

        size_t total_nelements_data_;

        size_t data_block_size_;

        size_t metadata_block_size_;

        usi nindex_;

        usi depth_;

        usi degeneracy_;
        
        float maxlog_;

        bool is_subblock_;

        StorageBlockPtr metadata_block_;

        MemoryPoolPtr metadata_mempool_;

        TensorBranch* obsolete_branch_;

        TensorController* controller_;

        TensorController* accumulate_controller_;

        TensorController* read_controller_;

        TensorController* write_controller_;

        /**
            This is a mode for mixed read/write operations
            where BOTH the actual values need to be known
            and the ability to write new values after
            reading is necessary.
        */
        TensorController* verbatim_controller_;

        bool permutationally_unique_;

        /********************************************
        This information is loaded when fetching data
        remotely, and is therefore mallocd from the
        metadata pool
        *********************************************/
        TensorBranch* branch_;

        typedef enum { metadata_block, data_block } storage_block_type_t;

        StorageBlock*
        allocate_storage_block(size_t size, storage_block_type_t blocktype);

        void _obsolete();

        void _retrieve();

        void _release();


    public:
        TensorBlock(
            const uli* indexset,
            Tensor* parent,
            Tensor::tensor_storage_t type
        );

        TensorBlock(
            TensorBlock* parent,
            Permutation* p
        );

        TensorBlock(TensorBlock* parent);

        ~TensorBlock();

        void accumulate(
            TensorBlock* src,
            double scale,
            Sort* sort
        );

        StorageBlock*
        allocate_data_storage_block();

        StorageBlock*
        allocate_metadata_storage_block();

        void accumulate(
            TensorBlock* lblock,
            TensorBlock* rblock,
            Contraction* cxn
        );

        void convert(
            Matrix* matrix,
            MatrixConfiguration* config,
            TensorIndexDescr* descr
        );

        void accumulate(
            Matrix* matrix,
            MatrixConfiguration* config,
            TensorIndexDescr* descr
        );

        /**
        */
        void recompute_tensor_permutation();

        void recompute_uniqueness();

        void reinit(Tensor::tensor_storage_t type);

        void set_element_size(size_t size);

        void sort_tensor_permutation(Permutation* p);

        void sort_resort_permutation(Permutation* p);

        void internal_contraction(TensorBlock* dst_block, MatrixConfiguration* config);

        void element_op(ElementOp* op);

        void configure(const TensorRetrieveActionPtr& action);

        TensorRetrieveAction* get_retrieve_action() const;

        /**
            Passed in by pointer since the pointer gets
            updated to new position upon getting passed in
        */
        bool equals(const void** data);
        
        bool equals(TensorBlock* block);

        void fill(TensorElementComputer* filler);

        void free_data(char* data, size_t size);

        void free_metadata(char* data, size_t size);

        TensorBranch* get_branch() const;

        uli get_block_number() const;

        usi get_depth() const;

        Tensor* get_parent_tensor() const;

        TensorIndexDescr* get_descr() const;

        DataNode* get_first_data_node() const;

        usi get_degeneracy() const;

        const uli* get_indices() const;

        std::string get_index_string() const;

        TensorController* get_tensor_controller() const;

        TensorBlock* get_symmetry_unique_block() const;

        MemoryPool* get_metadata_mempool() const;

        StorageBlock* get_metadata_storage_block() const;
        
        void clear_metadata();

        void free_metadata();

        Permutation* get_block_permutation() const;

        TensorBranch* get_obsolete_branch() const;

        Permutation* get_resort_permutation() const;

        Permutation* get_fetch_permutation() const;

        float get_max_log() const;

        void init_branch();

        bool in_destructor() const;

        bool is_cache_coherent() const;

        bool is_permutationally_unique() const;

        bool is_parent_in_current_tensor() const;

        bool is_remote_block() const;

        bool in_cache() const;

        bool in_read_mode() const;

        bool in_write_mode() const;

        bool in_accumulate_mode() const;

        bool in_verbatim_mode() const;

        bool is_subblock() const;

        void set_subblock(bool flag);
        
        void set_degeneracy(usi degeneracy);

        bool set_accumulate_mode();

        void set_read_mode();

        void set_write_mode();

        void set_verbatim_mode();

        void set_obsolete_branch(TensorBranch* branch);
        
        void set_max_log(float maxlog);

        void data_sort(Sort* sort);

        void metadata_sort(Sort* sort);

        bool is_synced() const;

        void sync();
        
        void sync_max_log();

        void reset_degeneracy();
        
        void update();

        /*****
        Cachable interface
        *****/
        void flush_from_cache();

        void print(std::ostream& os = std::cout);
};

}

#endif
