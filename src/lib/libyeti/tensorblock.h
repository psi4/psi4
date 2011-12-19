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
#include "messenger.h"

#include "mapimpl.h"
#include "cache.h"
#include "taskqueue.h"

#include "gigmatrix.h"

#include <list>
#include <csignal>

#ifdef redefine_size_t
#define size_t custom_size_t
#endif


#define DATA_BLOCK_SIZE_IN_CORE_CUTOFF 50000

#define NODE_NUMBER_UNINITIALIZED   999999999
#define MALLOC_NUMBER_UNINITIALIZED 999999999

#define TRACK_TENSOR_BLOCK_HISTORY 1

namespace yeti {


struct DataStorageNode :
    public Malloc<DataStorageNode>
{

    StorageBlock* block;

    DataStorageNode* next;

};

class TensorBlockAccumulateTask :
    public Task
{
    public:
        TensorBlock* tmp_block;

        TensorBlock* parent_block;

        void run(uli threadnum);
        
        void print(std::ostream& os) const;

        void prefetch(uli threadnum);

        void finalize_task_subset(uli threadnum);

        uli append_info(uli* data) const;

        void* operator new(size_t size);

        void operator delete(void* ptr);
};

class DataStorageBlockAllocator :
    public Cachable
{
    protected:
        DataStorageNode* head_node_;
        
#if TRACK_TENSOR_BLOCK_HISTORY
        std::list<const char*> history_;

        void print_history() const;
#endif

    public:
        DataStorageBlockAllocator();

        DataStorageBlockAllocator(YetiRuntimeObject::thread_safety_flag_t flag);

        ~DataStorageBlockAllocator();

        DataStorageNode* head_data_node() const;

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

#if TRACK_TENSOR_BLOCK_HISTORY
        void add_event(const char* str);

        void add_event(TensorBlock* block);
#endif
        
        virtual void flush_from_cache();
};

class TensorController :
    public smartptr::Countable,
    public Malloc<TensorController>
{

    public:
        TensorController();

        virtual void retrieve(TensorBlock* block) = 0;

        virtual void clear(TensorBlock* block) = 0;

        virtual void flush(TensorBlock* block) = 0;
        
        virtual void release(TensorBlock* block) = 0;

        virtual void retrieve(MetaDataNode* mdnode) = 0;

        virtual bool need_recompute_data() = 0;

        virtual void retrieve(
            DataNode* dnode,
            const MetaDataNode* parent,
            const uli* indices
        ) = 0;

        virtual void renew(TensorBlock* block) = 0;

        virtual void flush_data(TensorBlock* block) = 0;

        virtual void validate(TensorBlock* block) = 0;

        virtual void obsolete(TensorBlock* block) = 0;

        virtual void sync(TensorBlock* block) = 0;

        virtual void finalize(TensorBlock* block) = 0;

        virtual void out_of_core_prefetch(TensorBlock* block) = 0;

        virtual void in_core_prefetch(TensorBlock* current, TensorBlock* prev) = 0;

};


class TensorBlock :
    public DataStorageBlockAllocator,
    public Malloc<TensorBlock>,
    public Sendable
{
    protected:
        SendStatus* _send_data(
            YetiMessenger* messenger,
            Message::message_data_type_t type,
            Message::message_action_t remote_action,
            uli nmessages
        );

    private:
        friend class ResortInCorePrefetch;

        static TensorController* in_core_controller_;

        static TensorController* in_core_sorted_write_controller_;

        static TensorController* in_core_sorted_read_controller_;

        static TensorController* in_core_sorted_verbatim_controller_;

        static TensorController* in_core_sorted_accumulate_controller_;

        static TensorController* disk_read_controller_;

        static TensorController* disk_write_controller_;

        static TensorController* disk_accumulate_controller_;

        static TensorController* disk_sorted_read_controller_;

        static TensorController* remote_read_controller_;

        static TensorController* remote_write_controller_;

        static TensorController* remote_verbatim_controller_;

        static TensorController* remote_accumulate_controller_;

        static TensorController* remote_sorted_read_controller_;

        static TensorController* action_read_controller_;

        static TensorController* action_sorted_read_controller_;

        static TensorController* abort_read_controller_;

        static TensorController* abort_write_controller_;

        static TensorController* abort_accumulate_controller_;

        static TensorController* abort_verbatim_controller_;

        static TensorController* recompute_read_controller_;

        static TensorController* tmp_thread_accumulate_controller_;

        static TensorController* tmp_remote_accumulate_controller_;

        static bool statics_done_;

        bool remote_wait_;

        bool tmp_block_;

        bool is_prefetched_;

        uli prefetch_count_;

        TensorBlock* unique_block_;

        uli node_number_;

        uli malloc_number_;

        bool is_synced_;

        const char* controller_fail(TensorController* controller) const;

        void init_mempool();

        void init_in_core();

        void init_on_disk();

        void init_recomputed();

        void init_remote();

        void init_local();

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

        void init_tmp_accumulate();

        void set_synced();

        void set_unsynced();

        bool need_parent_retrieve() const;

        void
        receive_branch(
            YetiMessenger* messenger,
            uli proc_number,
            size_t metadata_size,
            uli initial_tag,
            uli nmessages
        );

        friend class TensorBranch;
        friend class FlushOldBranchRenew;

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

        TensorRetrieveActionPtr action_;

        Tensor::tensor_storage_t storage_type_;

        Tensor* parent_tensor_;

        uli indices_[NINDEX];

        usi degeneracy_;

        bool is_subblock_;

        bool up_to_date_;

        bool task_owner_;

        StorageBlock* metadata_block_;

        MemoryPool* metadata_mempool_;

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

        bool fast_read_;

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

        void _renew();

        void _initialize();

        void _finalize();

        bool set_accumulate_mode();

        void set_read_mode();

        void set_write_mode();

        void set_verbatim_mode();

        void enter_recv(Message::message_action_t action);

        void exit_recv(Message::message_action_t action);

        TensorBlock(TensorBlock* block);

    public:
        TensorBlock(
            const uli* indexset,
            Tensor* parent,
            Tensor::tensor_storage_t type
        );

        /** Used for temporary accumulation */
        TensorBlock(Tensor* parent);


        ~TensorBlock();

        void accumulate(
            TensorBlock* src,
            double scale,
            Sort* sort
        );

        StorageBlock*
        allocate_data_storage_block(size_t size);

        StorageBlock*
        allocate_metadata_storage_block();

        void controller_fail() const;

        void recompute_permutation();

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


        void assign(TensorBlock* block);

        /**
        */
        void recompute_tensor_permutation();

        void init();

        void reinit();

        void uninit();

        static void init_statics();

        static void delete_statics();

        void set_as_unique();

        void set_element_size(size_t size);

        void sort_tensor_permutation(Permutation* p);

        void sort_resort_permutation(Permutation* p);

        void internal_contraction(TensorBlock* dst_block, MatrixConfiguration* config);

        void element_op(ElementOp* op);

        void configure(const TensorRetrieveActionPtr& action);

        void configure(Tensor::tensor_storage_t type);

        TensorRetrieveAction* get_retrieve_action() const;

        usi nindex() const;

        usi depth() const;

        size_t data_block_size() const;

        size_t metadata_block_size() const;

        TensorIndexDescr* descr() const;

        void finalize();

        void retrieve_accumulate();

        void retrieve_accumulate_no_lock();

        void retrieve_read();

        void retrieve_write();

        void retrieve_verbatim();

        void retrieve_verbatim_no_lock();

        bool is_nonnull() const;

        void release_accumulate();

        void release_read();

        void release_write();

        void release_verbatim();

        void wait_on_remote() const;

        void set_remote_wait(bool flag);

        bool is_task_owner() const;

        void set_task_owner(bool flag);

        bool is_flushable() const;

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

        uli get_malloc_number() const;

        usi get_depth() const;

        uli get_nelements_data() const;

        Tensor* get_parent_tensor() const;

        TensorIndexDescr* get_descr() const;

        DataNode* get_first_data_node() const;

        usi get_degeneracy() const;

        const uli* get_indices() const;

        std::string get_index_string() const;

        std::string get_block_name() const;

        TensorController* get_tensor_controller() const;

        TensorBlock* get_symmetry_unique_block() const;

        MemoryPool* get_metadata_mempool() const;

        StorageBlock* get_metadata_storage_block() const;
        
        void clear_metadata();

        void free_metadata();

        Permutation* get_block_permutation() const;

        Permutation* get_resort_permutation() const;

        Permutation* get_fetch_permutation() const;

        void init_branch();

        bool has_fast_read() const;

        bool is_prefetched() const;

        void set_prefetched(bool flag);

        bool is_up_to_date() const;

        bool in_destructor() const;

        bool is_cache_coherent() const;

        bool is_permutationally_unique() const;

        bool is_parent_in_current_tensor() const;

        bool is_remote_block() const;

        bool is_local_block() const;

        bool is_recomputed_block() const;

        bool has_controller() const;

        bool has_remote_controller() const;

        bool is_cached() const;

        bool in_read_mode() const;

        bool in_write_mode() const;

        bool in_accumulate_mode() const;

        bool in_verbatim_mode() const;

        bool is_subblock() const;

        void set_subblock(bool flag);
        
        void set_degeneracy(usi degeneracy);

        void data_sort(Sort* sort);

        void metadata_sort(Sort* sort);

        bool is_synced() const;

        void sync();

        void reset_degeneracy();
        
        void update();

        void prefetch_read();

        void prefetch_accumulate();

        void prefetch_verbatim();

        void prefetch_write();

        void complete_read();

        void prefetch_read(TensorBlock* prev_block);

        void prefetch_accumulate(TensorBlock* prev_block);

        void prefetch_verbatim(TensorBlock* prev_block);

        void prefetch_write(TensorBlock* prev_block);

        void configure_tmp_block(TensorBlock* block);

        /*****
        Cachable interface
        *****/
        void flush_from_cache();

        void print(std::ostream& os = std::cout);

        uli get_node_number() const;

        uli get_node_owner() const;

        void set_node_number(uli node);
        
        void wtf_retrieve();

        void wtf_release();

        void recv_data(
            YetiMessenger* messenger,
            Message::message_data_type_t type,
            Message::message_action_t action,
            uli proc_number,
            size_t metadata_size,
            uli initial_tag,
            uli nmessages
        );

        void redistribute(
            YetiMessenger* messenger,
            Message::message_data_type_t type,
            uli proc_number,
            size_t metadata_size,
            uli initial_tag,
            uli nmessages
        );


        bool is_waiting_on_redistribute() const;

        void accumulate_data(
            YetiMessenger* messenger,
            Message::message_data_type_t type,
            uli proc_number,
            size_t metadata_size,
            uli initial_tag,
            uli nmessages
        );

        


};

}

#endif
