#include "tensor.h"
#include "tensorblock.h"
#include "tensorbranch.h"
#include "tensoraction.h"
#include "node.h"
#include "index.h"
#include "class.h"
#include "exception.h"
#include "permutationimpl.h"
#include "env.h"
#include "data.h"
#include "runtime.h"
#include "dataimpl.h"
#include "threadimpl.h"
#include "filler.h"
#include "sortimpl.h"
#include "tensorimpl.h"
#include "datapolicies.h"
#include "metadatapolicies.h"
#include "branchpolicies.h"
#include "elementop.h"
#include "matrix.h"
#include "contraction.h"

#include <libsmartptr/strop.h>

using namespace yeti;
using namespace std;

#define DEBUG_TENSOR_BLOCK_RETRIEVE 0
#define PRINT_MALLOC_NUMBERS 0

#define block_fail(x,y) { controller_fail(); yeti_throw(x,y); }

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

TensorController* TensorBlock::in_core_controller_ = 0;
TensorController* TensorBlock::in_core_sorted_write_controller_ = 0;
TensorController* TensorBlock::in_core_sorted_read_controller_ = 0;
TensorController* TensorBlock::in_core_sorted_verbatim_controller_ = 0;
TensorController* TensorBlock::in_core_sorted_accumulate_controller_ = 0;
TensorController* TensorBlock::disk_read_controller_ = 0;
TensorController* TensorBlock::disk_write_controller_ = 0;
TensorController* TensorBlock::disk_accumulate_controller_ = 0;
TensorController* TensorBlock::disk_sorted_read_controller_ = 0;
TensorController* TensorBlock::remote_read_controller_ = 0;
TensorController* TensorBlock::remote_write_controller_ = 0;
TensorController* TensorBlock::remote_verbatim_controller_ = 0;
TensorController* TensorBlock::remote_accumulate_controller_ = 0;
TensorController* TensorBlock::remote_sorted_read_controller_ = 0;
TensorController* TensorBlock::action_read_controller_ = 0;
TensorController* TensorBlock::action_sorted_read_controller_ = 0;
TensorController* TensorBlock::abort_read_controller_ = 0;
TensorController* TensorBlock::abort_write_controller_ = 0;
TensorController* TensorBlock::abort_accumulate_controller_ = 0;
TensorController* TensorBlock::abort_verbatim_controller_ = 0;
TensorController* TensorBlock::recompute_read_controller_ = 0;
TensorController* TensorBlock::tmp_thread_accumulate_controller_ = 0;
TensorController* TensorBlock::tmp_remote_accumulate_controller_ = 0;
bool TensorBlock::statics_done_ = false;



DECLARE_MALLOC(TensorBlock);
DECLARE_MALLOC(DataStorageNode);


#define NBLOCKS_TMP_ACCUMULATE 1000
static char tmp_malloc_data[NBLOCKS_TMP_ACCUMULATE * sizeof(TensorBlock)];
static char tmp_malloc_mallocd[NBLOCKS_TMP_ACCUMULATE];
static FastMalloc tensor_block_tmp_malloc(
                    tmp_malloc_data, 
                    tmp_malloc_mallocd, 
                    sizeof(TensorBlock), 
                    NBLOCKS_TMP_ACCUMULATE, 
                    "tensor block tmp"
                  );
static TensorBlockAccumulateTask accumulate_tasks[NBLOCKS_TMP_ACCUMULATE];
static uli accumulate_task_num = 0;

typedef TensorControllerTemplate<
        ValidBranchController,
        DoNothingBranchRenew,
        DoNothingMempool,
        DoNothingBranchRetrieve,
        DoNothingDataControllerRetrieve,
        DoNothingDataControllerInit,
        DoNothingFinalize,
        DoNothingBranchFlush,
        DoNothingBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        DoNothingDataStorageFlush,
        AbortOnObsolete,
        DoNothingSync,
        DeleteAllDataClear,
        UnlockAfterPrefetch,
        DoNothingInCorePrefetch
    > InCoreController;


typedef TensorControllerTemplate<
        ValidBranchController,
        CacheBranchRenew,
        ResetMempool,
        SortedBranchRetrieve,
        ReallocateDataControllers,
        SortDataControllers,
        DoNothingFinalize,
        ClearBranchFlush,
        CacheBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        DeleteAllDataStorageFlush,
        SetFinalizedOnObsolete,
        DoNothingSync,
        ClearAllData,
        InitializePrefetch,
        DoNothingInCorePrefetch
    > InCoreSortedReadController;

typedef TensorControllerTemplate<
        //ValidBranchController,
        AbortAccumulateBranchValidation,
        CacheBranchRenew,
        ResetMempool,
        NewBranchRetrieve,
        DoNothingDataControllerRetrieve,
        DoNothingDataControllerInit,
        DoNothingFinalize,
        SortedAccumulateBranchFlush,
        CacheBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        DeleteAllDataStorageFlush,
        SetFinalizedOnObsolete,
        FlushOnSync,
        ClearAllData,
        InitializePrefetch,
        DoNothingInCorePrefetch
   > InCoreSortedAccumulateController;

typedef TensorControllerTemplate<
        AbortReadBranchValidation,
        DoNothingBranchRenew,
        DoNothingMempool,
        DoNothingBranchRetrieve,
        DoNothingDataControllerRetrieve,
        DoNothingDataControllerInit,
        DoNothingFinalize,
        DoNothingBranchFlush,
        DoNothingBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        DoNothingDataStorageFlush,
        AbortOnObsolete,
        DoNothingSync,
        DeleteAllDataClear,
        UnlockAfterPrefetch,
        DoNothingInCorePrefetch
    > AbortReadController;

typedef TensorControllerTemplate<
        AbortWriteBranchValidation,
        DoNothingBranchRenew,
        DoNothingMempool,
        DoNothingBranchRetrieve,
        DoNothingDataControllerRetrieve,
        DoNothingDataControllerInit,
        DoNothingFinalize,
        DoNothingBranchFlush,
        DoNothingBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        DoNothingDataStorageFlush,
        AbortOnObsolete,
        DoNothingSync,
        DeleteAllDataClear,
        UnlockAfterPrefetch,
        DoNothingInCorePrefetch
    > AbortWriteController;

typedef TensorControllerTemplate<
        AbortAccumulateBranchValidation,
        DoNothingBranchRenew,
        DoNothingMempool,
        DoNothingBranchRetrieve,
        DoNothingDataControllerRetrieve,
        DoNothingDataControllerInit,
        DoNothingFinalize,
        DoNothingBranchFlush,
        DoNothingBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        DoNothingDataStorageFlush,
        AbortOnObsolete,
        DoNothingSync,
        DeleteAllDataClear,
        UnlockAfterPrefetch,
        DoNothingInCorePrefetch
    > AbortAccumulateController;

typedef TensorControllerTemplate<
        AbortVerbatimBranchValidation,
        DoNothingBranchRenew,
        DoNothingMempool,
        DoNothingBranchRetrieve,
        DoNothingDataControllerRetrieve,
        DoNothingDataControllerInit,
        DoNothingFinalize,
        DoNothingBranchFlush,
        DoNothingBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        DoNothingDataStorageFlush,
        AbortOnObsolete,
        DoNothingSync,
        DeleteAllDataClear,
        UnlockAfterPrefetch,
        DoNothingInCorePrefetch
    > AbortVerbatimController;

typedef TensorControllerTemplate<
        ValidBranchController,
        CacheBranchRenew,
        ResetMempool,
        ActionBranchRetrieve,
        ReuseStorageBlocks,
        DoNothingDataControllerInit,
        DoNothingFinalize,
        ClearBranchFlush,
        CacheBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        DeleteAllDataStorageFlush,
        AbortOnObsolete,
        DoNothingSync,
        ClearAllData,
        UnlockAfterPrefetch,
        DoNothingInCorePrefetch
    > ActionReadController;


typedef TensorControllerTemplate<
        ValidBranchController,
        CacheBranchRenew,
        ResetMempool,
        SortedBranchRetrieve,
        ReallocateDataControllers,
        SortDataControllers,
        DoNothingFinalize,
        ClearBranchFlush,
        CacheBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        DeleteAllDataStorageFlush,
        SetFinalizedOnObsolete,
        DoNothingSync,
        DeleteAllDataClear,
        UnlockAfterPrefetch,
        DoNothingInCorePrefetch
    > ActionSortedReadController;

typedef TensorControllerTemplate<
        ValidBranchController,
        CacheBranchRenew,
        DoNothingMempool, //this is taken care of in the prefetch
        RemoteBlockBranchRetrieve,
        DoNothingDataControllerRetrieve,
        DoNothingDataControllerInit,
        DoNothingFinalize,
        ClearBranchFlush,
        CacheBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        DeleteAllDataStorageFlush,
        SetFinalizedOnObsolete,
        DoNothingSync,
        DeleteAllDataClear,
        RemoteBlockPrefetch,
        DoNothingInCorePrefetch
    > RemoteReadController;

typedef TensorControllerTemplate<
        ValidBranchController,
        CacheBranchRenew,
        ResetMempool,
        DoNothingBranchRetrieve,
        ReuseStorageBlocks,
        DoNothingDataControllerInit,
        RemoteAccumulateFinalize,
        ClearBranchFlush,
        FinalizeOnRelease,
        //CacheBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        DeleteAllDataStorageFlush,
        SetFinalizedOnObsolete,
        FlushOnSync,
        DeleteAllDataClear,
        UnlockAfterPrefetch,
        DoNothingInCorePrefetch
    > RemoteAccumulateController;

typedef TensorControllerTemplate<
        ValidBranchController,
        DoNothingBranchRenew,
        ResetMempool,
        DoNothingBranchRetrieve,
        ReuseStorageBlocks,
        DoNothingDataControllerInit,
        RemoteAccumulateFinalize,
        ClearBranchFlush,
        DoNothingBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        DeleteAllDataStorageFlush,
        SetFinalizedOnObsolete,
        AbortOnSync,
        DeleteAllDataClear,
        UnlockAfterPrefetch,
        DoNothingInCorePrefetch
    > TmpRemoteAccumulateController;

typedef TensorControllerTemplate<
        ValidBranchController,
        DoNothingBranchRenew,
        ResetMempool,
        DoNothingBranchRetrieve,
        ReuseStorageBlocks,
        DoNothingDataControllerInit,
        ThreadAccumulateFinalize,
        ClearBranchFlush,
        DoNothingBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        DeleteAllDataStorageFlush,
        SetFinalizedOnObsolete,
        AbortOnSync,
        DeleteAllDataClear,
        UnlockAfterPrefetch,
        DoNothingInCorePrefetch
    > TmpThreadAccumulateController;


typedef TensorControllerTemplate<
        ValidBranchController,
        CacheBranchRenew,
        ResetMempool,
        SortedBranchRetrieve,
        ReallocateDataControllers,
        SortDataControllers,
        DoNothingFinalize,
        ClearBranchFlush,
        CacheBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        DeleteAllDataStorageFlush,
        SetFinalizedOnObsolete,
        DoNothingSync,
        ClearAllData,
        InitializePrefetch,
        DoNothingInCorePrefetch
    > RemoteSortedReadController;

typedef TensorControllerTemplate<
        ValidBranchController,
        ConfigureElementComputerRenew,
        ResetMempool,
        ConfigureElementComputerRetrieve,
        ReuseStorageBlocks,
        DoNothingDataControllerInit,
        DoNothingFinalize,
        ClearBranchFlush,
        CacheBranchRelease,
        RecomputeMetaDataRetrieve,
        RecomputeDataRetrieve,
        DeleteAllDataStorageFlush,
        SetFinalizedOnObsolete,
        DoNothingSync,
        ClearAllData,
        InitializePrefetch,
        ResortInCorePrefetch
    > RecomputeReadController;

typedef TensorControllerTemplate<
        ValidBranchController,
        DoNothingBranchRenew,
        DoNothingMempool,
        RealignMemoryPoolBranchRetrieve,
        ReuseDataControllers,
        DoNothingDataControllerInit,
        DoNothingFinalize,
        ClearBranchFlush,
        CacheBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        ClearAllDataStorageFlush,
        AbortOnObsolete,
        DoNothingSync,
        ClearAllData,
        UnlockAfterPrefetch,
        DoNothingInCorePrefetch
    > DiskReadController;

typedef TensorControllerTemplate<
        ValidBranchController,
        DoNothingBranchRenew,
        DoNothingMempool,
        RealignMemoryPoolBranchRetrieve,
        ReuseDataControllers,
        DoNothingDataControllerInit,
        DoNothingFinalize,
        CommitBranchFlush,
        CacheBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        CommitAllDataStorageFlush,
        AbortOnObsolete,
        DoNothingSync,
        ClearAllData,
        UnlockAfterPrefetch,
        DoNothingInCorePrefetch
    > DiskWriteController;
    

typedef TensorControllerTemplate<
        ValidBranchController,
        DoNothingBranchRenew,
        DoNothingMempool,
        RealignMemoryPoolBranchRetrieve,
        ReuseDataControllers,
        DoNothingDataControllerInit,
        DoNothingFinalize,
        CommitBranchFlush,
        CacheBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        CommitAllDataStorageFlush,
        AbortOnObsolete,
        DoNothingSync,
        ClearAllData,
        UnlockAfterPrefetch,
        DoNothingInCorePrefetch
    > DiskAccumulateController;

typedef TensorControllerTemplate<
        ValidBranchController,
        DoNothingBranchRenew,
        ResetMempool,
        SortedBranchRetrieve,
        ReallocateDataControllers,
        SortDataControllers,
        DoNothingFinalize,
        ClearBranchFlush,
        CacheBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        DeleteAllDataStorageFlush,
        SetFinalizedOnObsolete,
        DoNothingSync,
        ClearAllData,
        UnlockAfterPrefetch,
        DoNothingInCorePrefetch
    > DiskSortedReadController;

    

DECLARE_PARENT_MALLOC(TensorController);
DECLARE_SUB_MALLOC(TensorController,AbortReadController);

DataStorageBlockAllocator::DataStorageBlockAllocator()
    :
  head_node_(0)
{
}

DataStorageBlockAllocator::DataStorageBlockAllocator(YetiRuntimeObject::thread_safety_flag_t flag)
    :
  Cachable(flag),
  head_node_(0)
{
}

DataStorageBlockAllocator::~DataStorageBlockAllocator()
{
    clear_storage();
}

uli
DataStorageBlockAllocator::get_num_data_storage_blocks() const
{
    uli n = 0;
    DataStorageNode* node = head_node_;
    while(node)
    {
        ++n;
        node = node->next;
    }
    return n;
}

DataStorageNode*
DataStorageBlockAllocator::head_data_node() const
{
    return head_node_;
}

void
DataStorageBlockAllocator::store(StorageBlock* block)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("store block");
#endif
    DataStorageNode* node = new DataStorageNode;
    node->block = block;
    node->next = 0;

    if (head_node_ == 0)
    {
        head_node_ = node;
    }
    else
    {
        node->next = head_node_;
        head_node_ = node;
    }
}

void
DataStorageBlockAllocator::flush_from_cache()
{
    //do nothing
}

CachedStorageBlock*
DataStorageBlockAllocator::allocate_cache_block(DataCache* cache)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("allocate cache block");
#endif
    CachedStorageBlock* block = new CachedStorageBlock(this, cache);
    block->retrieve();

    return block;
}

void
DataStorageBlockAllocator::clear_storage()
{
    DataStorageNode* node = head_node_;
    while (node)
    {
#if TRACK_TENSOR_BLOCK_HISTORY
        add_event("delete block");
#endif
        StorageBlock* block = node->block;
        delete block;
        DataStorageNode* next = node->next;
        DataStorageNode* prev = node;
        node->block = 0;
        node->next = 0;
        node = next;
        delete prev;
    }
    head_node_ = 0;
}

InCoreBlock*
DataStorageBlockAllocator::allocate_core_block(
    size_t size
)
{
    char* data = YetiRuntime::malloc(size);
    InCoreBlock* block = new InCoreBlock(data, size);
    block->retrieve();
    return block;
}

LocalDiskBlock*
DataStorageBlockAllocator::allocate_disk_block(
    DiskBuffer* buffer,
    Cachable* item,
    DataCache* cache
)
{
    LocalDiskBlock* block = new LocalDiskBlock(item,cache,buffer);
    block->retrieve();
    return block;
}

#define IN_CORE_CUTOFF 10000
TensorBlock::TensorBlock(
    const uli* indexset,
    Tensor* parent,
    Tensor::tensor_storage_t store_type
) :
    controller_(0),
    unique_block_(0),
    read_controller_(0),
    write_controller_(0),
    accumulate_controller_(0),
    verbatim_controller_(0),
    resort_perm_(0),
    fetch_perm_(0),
    branch_(0),
    metadata_mempool_(0),
    prefetch_count_(0),
    parent_tensor_(parent),
    degeneracy_(1),
    is_synced_(true), //no data yet, so true
    permutationally_unique_(false), //begin with false
    is_subblock_(false),
    block_perm_(0),
    action_(0),
    node_number_(NODE_NUMBER_UNINITIALIZED),
    malloc_number_(0),
    up_to_date_(true),
    storage_type_(Tensor::in_core),
    task_owner_(false),
    metadata_block_(0),
    tmp_block_(false),
    fast_read_(false),
    remote_wait_(false)
{
    malloc_number_ = Malloc<TensorBlock>::get_malloc_number(this);

    ::memcpy(indices_, indexset, descr()->nindex() * sizeof(uli));

    if (parent_tensor_->is_distributed())
    {
        block_fail(SanityCheckError, "cannot distribute tensor before creating all tensor blocks");
    }

    recompute_permutation();

    block_perm_ = parent_tensor_->get_tensor_grp()->get_identity();

    set_element_size(sizeof(double));

#if PRINT_MALLOC_NUMBERS
    cout << stream_printf("On node %d, malloc number %ld is %s\n", 
                          YetiRuntime::me(), malloc_number_,
                          get_block_name().c_str());
#endif
}

TensorBlock::TensorBlock(
    Tensor* parent
) :
    controller_(0),
    read_controller_(0),
    unique_block_(0),
    write_controller_(0),
    accumulate_controller_(0),
    verbatim_controller_(0),
    resort_perm_(0),
    fetch_perm_(0),
    branch_(0),
    metadata_mempool_(0),
    prefetch_count_(0),
    parent_tensor_(parent),
    degeneracy_(1),
    is_synced_(true), //no data yet, so true
    permutationally_unique_(true),
    is_subblock_(false),
    block_perm_(0),
    action_(0),
    node_number_(YetiRuntime::me()),
    malloc_number_(0),
    up_to_date_(true),
    storage_type_(Tensor::in_core),
    metadata_block_(0),
    task_owner_(false),
    tmp_block_(true),
    fast_read_(false),
    remote_wait_(false)
{
    malloc_number_ = MALLOC_NUMBER_UNINITIALIZED;
    for (uli i=0; i < descr()->nindex(); ++i)
        indices_[i] = 0;

    set_element_size(sizeof(double));

    init_tmp_accumulate();
}

TensorBlock::TensorBlock(
    TensorBlock* parent
) :
    //this object is created only within a thread
    DataStorageBlockAllocator(YetiRuntimeObject::thread_safe),
    controller_(0),
    read_controller_(0),
    unique_block_(0),
    write_controller_(0),
    accumulate_controller_(0),
    verbatim_controller_(0),
    resort_perm_(parent->get_resort_permutation()),
    fetch_perm_(parent->get_fetch_permutation()),
    branch_(0),
    metadata_mempool_(0),
    parent_tensor_(parent->get_parent_tensor()),
    degeneracy_(parent->get_degeneracy()),
    is_synced_(true), //no data yet, so true
    permutationally_unique_(true),
    is_subblock_(parent->is_subblock()),
    block_perm_(parent->get_block_permutation()),
    node_number_(YetiRuntime::me()),
    malloc_number_(0),
    prefetch_count_(0),
    up_to_date_(true),
    storage_type_(Tensor::in_core),
    task_owner_(false),
    tmp_block_(true),
    fast_read_(false),
    remote_wait_(false)
{
    malloc_number_ = parent->get_malloc_number();
    ::memcpy(indices_, parent->get_indices(), descr()->nindex() * sizeof(uli));

    //this is only ever used for temporary accumulation
    init_in_core_no_sort();
    controller_ = accumulate_controller_;

}

TensorBlock::~TensorBlock()
{
#if YETI_SANITY_CHECK
    if (!is_locked())
        block_fail(SanityCheckError, "Cannot delete an unlocked tensor block");
#endif
    /** The block might never have been initialized */
    if (read_controller_)
    {
        in_destructor_ = true;

        if (is_cached())
        {
            //always lock an object before flushing it
            //once the tensor is deleted, it goes out of cache
#if TRACK_TENSOR_BLOCK_HISTORY
            add_event("destructor flush");
#endif
            flush_from_cache();
            if (has_send_status())
                block_fail(SanityCheckError, "Cannot delete tensor block that has pending send");
        }

        wait_on_send();

        read_controller_->clear(this);

        clear_storage();
    }

    if (metadata_mempool_)
        delete metadata_mempool_;
    if (metadata_block_)
        delete metadata_block_;

}

usi
TensorBlock::nindex() const
{
    return descr()->nindex();
}

usi
TensorBlock::depth() const
{
    return descr()->depth() - 1;
}

size_t
TensorBlock::metadata_block_size() const
{
    return parent_tensor_->metadata_block_size();
}

size_t
TensorBlock::data_block_size() const
{
    return parent_tensor_->data_block_size();
}

TensorIndexDescr*
TensorBlock::descr() const
{
    return parent_tensor_->get_block_descr();
}

void
TensorBlock::set_as_unique()
{
    permutationally_unique_ = true;
    fetch_perm_ = parent_tensor_->get_tensor_grp()->get_identity();
    resort_perm_ = fetch_perm_;
}

void
TensorBlock::recompute_permutation()
{
    //erase any memory of permutational uniqueness
    //otherwise recompute tensor permutation just returns without doing
    //any work
    permutationally_unique_ = false;

    recompute_tensor_permutation();
    permutationally_unique_ = fetch_perm_->is_identity();
}

void
TensorBlock::reinit()
{
    uninit();
    init();
}

void
TensorBlock::uninit()
{
    if (!controller_) //already uninitialized
        return;

#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("uninit");
#endif

#if YETI_SANITY_CHECK
    if (!is_locked())
        block_fail(SanityCheckError, "Block must be locked before initializing");

    if (!metadata_block_)
        block_fail(SanityCheckError, "Block is already uninitialized");
#endif

    if (is_cached())
        flush_from_cache();

    read_controller_->clear(this);
    clear_storage();

    if (metadata_block_)
    {
        delete metadata_block_;
        metadata_block_ = 0;
    }
    if (metadata_mempool_)
    {
        delete metadata_mempool_;
        metadata_mempool_ = 0;
    }


    controller_ = 0;
    read_controller_ = 0;
    write_controller_ = 0;
    verbatim_controller_ = 0;
    accumulate_controller_ = 0;
    fast_read_ = false;
    set_initialized(false);
    set_prefetched(false);
    set_finalized(true);
}

void
TensorBlock::init()
{
    if (controller_) /** already initialized */
        return;

#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("init");
#endif

#if YETI_SANITY_CHECK
    if (!is_locked())
        block_fail(SanityCheckError, "Block must be locked before initializing");

    if (metadata_block_)
        block_fail(SanityCheckError, "Block is already initialized");
#endif

    if (is_remote_block() && storage_type_ != Tensor::recomputed)
    {
        init_remote();
    }
    else
    {
        init_local();
    }
    controller_ = read_controller_;
}

void
TensorBlock::init_local()
{
    switch(storage_type_)
    {
        case Tensor::in_core:
            init_in_core();
            break;
        case Tensor::on_disk:
            init_on_disk();
            break;
        case Tensor::recomputed:
            init_recomputed();
            break;
        case Tensor::action:
            init_action();
            break; //other initialization stuff happens elsewhere
    }
}

void
TensorBlock::init_statics()
{
    if (statics_done_)
        return;

    in_core_controller_ = new InCoreController;

    abort_read_controller_ = new AbortReadController;
    abort_write_controller_ = new AbortWriteController;
    abort_verbatim_controller_ = new AbortVerbatimController;
    abort_accumulate_controller_ = new AbortAccumulateController;

    in_core_sorted_read_controller_ = new InCoreSortedReadController;
    in_core_sorted_accumulate_controller_ = new InCoreSortedAccumulateController;

    action_read_controller_ = new ActionReadController;
    action_sorted_read_controller_ = new ActionSortedReadController;

    remote_read_controller_ = new RemoteReadController;
    remote_accumulate_controller_ = new RemoteAccumulateController;

    remote_sorted_read_controller_ = new RemoteSortedReadController;

    recompute_read_controller_ = new RecomputeReadController;

    disk_read_controller_ = new DiskReadController;
    disk_accumulate_controller_ = new DiskAccumulateController;
    disk_write_controller_ = new DiskWriteController;

    disk_sorted_read_controller_ = new DiskSortedReadController;

    tmp_thread_accumulate_controller_ = new TmpThreadAccumulateController;
    tmp_remote_accumulate_controller_ = new TmpRemoteAccumulateController;

    statics_done_ = true;
}

void
TensorBlock::delete_statics()
{
    if (!statics_done_)
        return;

    delete in_core_controller_;

    delete abort_read_controller_;
    delete abort_write_controller_;
    delete abort_verbatim_controller_;
    delete abort_accumulate_controller_;

    delete in_core_sorted_read_controller_;
    delete in_core_sorted_accumulate_controller_;

    delete action_read_controller_;
    delete action_sorted_read_controller_;

    delete remote_read_controller_;
    delete remote_accumulate_controller_;

    delete remote_sorted_read_controller_;

    delete recompute_read_controller_;

    delete disk_read_controller_;
    delete disk_accumulate_controller_;
    delete disk_write_controller_;

    delete disk_sorted_read_controller_;

    delete tmp_thread_accumulate_controller_;
    delete tmp_remote_accumulate_controller_;

    statics_done_ = false;
}

/*************
Begin public interface definition
*************/

void
TensorBlock::accumulate(
    TensorBlock* src,
    double scale,
    Sort* sort
)
{
#if 0
    for (uli i=0; i < descr()->nindex(); ++i)
    {
        if (src->get_indices()[i] != indices_[i])
        {
            cerr << "Accumulate error" << endl;
            cerr << src->get_block_name().c_str() << endl;
            cerr << get_block_name().c_str() << endl;
            abort();
        }
    }
#endif
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("accumulate block");
#endif

    up_to_date_ = false;

    branch_->set_element_type(src->get_branch()->element_type());
    if (src == this) //self-accumulation
    {
        if (!sort)
            block_fail(SanityCheckError, "you should never self-accumulate unless sorting!");

        CachedStorageBlock* block
            = new CachedStorageBlock(this, parent_tensor_->get_metadata_cache());
        block->retrieve();
        MemoryPool* mempool = new MemoryPool(block->size(), block->data());
        mempool->memcpy(metadata_mempool_);
        TensorBranch* sorted_branch = reinterpret_cast<TensorBranch*>(mempool->data());
        sorted_branch->set_metadata_mempool(mempool);
        branch_->sort_metadata_into(sort, sorted_branch);

        //now that the sort of metadata has taken place,
        //realign the memory pool
        sorted_branch->realign_memory_pool(
            branch_, // "old" branch
            sorted_branch // "new" branch
        );

        TensorDataController* controller = sorted_branch->first_data_controller();

        if (controller == 0)
        {
            delete mempool;
            delete block;
            return; //nothing doin
        }


        DataStorageBlockAllocator allocator;
        while (controller)
        {
            CachedStorageBlock* data_block
                   = allocator.allocate_cache_block(parent_tensor_->get_data_cache());
            allocator.store(data_block);
            //block is retrieved at this point
            controller->set_data(data_block->data());
            controller = controller->next;
        }

        branch_->get_node()->sort_data_into(
            sort,
            branch_,
            sorted_branch
        );

        branch_->get_node()->accumulate(
            sorted_branch->get_node(),
            1.0,
            0 //no sort
        );

        uli nblocks = allocator.get_num_data_storage_blocks();
        DataStorageNode* node = allocator.head_data_node();
        while (node)
        {
            CachedStorageBlock* block = static_cast<CachedStorageBlock*>(node->block);
        }
        delete mempool;
        block->release();
        delete block;
    }
    else
    {
#if YETI_SANITY_CHECK
        if (branch_->get_node() == 0)
        {
            cerr << "block " << get_index_string() << " has no metadata node" << endl;
            abort();
        }
        else if (src->get_branch()->get_node() == 0)
        {
            cerr << "block " << src->get_index_string() 
                << " on tensor " << src->get_parent_tensor()->get_name()
                << " has no metadata node" << endl;
            abort();
        }
#endif
        branch_->get_node()->accumulate(
            src->get_branch()->get_node(),
            scale,
            sort
        );
    }
    
    set_unsynced();
}

void
TensorBlock::accumulate(
    TensorBlock* lblock,
    TensorBlock* rblock,
    Contraction* cxn
)
{
#if YETI_SANITY_CHECK
    if (   lblock->get_degeneracy() == 0
        || rblock->get_degeneracy() == 0)
        block_fail(SanityCheckError, "block has zero degeneracy");
#endif
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("accumulate cxn");
#endif

    branch_->set_element_type(lblock->get_branch()->element_type());
    controller_->retrieve(branch_->get_node());

    branch_->get_node()->accumulate(
        lblock->get_branch()->get_node(),
        rblock->get_branch()->get_node(),
        lblock->get_degeneracy(),
        cxn
    );


    /** This is no longer compatible with its source destination */
    set_unsynced();
}


StorageBlock*
TensorBlock::allocate_data_storage_block(size_t size)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("allocate data storage block");
#endif
#if DEBUG_TENSOR_BLOCK_RETRIEVE 
    if (tmp_block_)
    {
        cout << stream_printf("Allocating data block of size %ld on %s\n",
                    size,
                    get_block_name().c_str());
        cout.flush();
    }
#endif
    if (parent_tensor_->nelements_metadata_av() < 20)//just assign the entire block a data segment
        return allocate_storage_block(size, TensorBlock::data_block);
    else //tile many blocks to minimize number of malloc calls and storage blocks
        return allocate_storage_block(parent_tensor_->data_block_size(), TensorBlock::data_block);
}

StorageBlock*
TensorBlock::allocate_storage_block(
    size_t size,
    TensorBlock::storage_block_type_t blocktype
)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("allocate storage block");
#endif
    Tensor::tensor_storage_t storage_type
        = parent_tensor_->get_storage_type();

    StorageBlock* block = 0;
    DataCache* cache = blocktype == TensorBlock::metadata_block ?
                parent_tensor_->get_metadata_cache() :
                parent_tensor_->get_data_cache();

    if (tmp_block_)
    {
        block = allocate_core_block(size);
    }
    else if (!permutationally_unique_ ||
        storage_type == Tensor::recomputed ||
        action_ ||
        is_remote_block())
    {
        if (cache->blocksize() < size)
        {
            cerr << stream_printf("Cache blocks are not large enough!\n"
                                "Cache block size is %ld but requested size is %ld\n",
                                cache->blocksize(), size);
            controller_fail();
            abort();
        }
        try{
            block = allocate_cache_block(cache);
        } catch (int e) {
            cerr << stream_printf("Error allocating new cache block on %s.  %d blocks allocated so far.\n",
                                    get_block_name().c_str(), get_num_data_storage_blocks()
                                );
            cerr.flush();
            throw e;
        }
    }
    else if (storage_type == Tensor::in_core)
    {
        block = this->allocate_core_block(size);
        YetiRuntime::register_allocation(parent_tensor_, size);
    }
    else if (storage_type == Tensor::on_disk)
    {

        if (cache->blocksize() < size)
        {
            cerr << "cache blocks are not large enough!"
                 << endl
                 << " cache block size is " << cache->blocksize()
                 << " but requested block size is " << size
                 << endl;
        }
        yeti_throw(SanityCheckError, "No disk usage yet");
    }
    else
    {
        cerr << "now you've screwed it up" << endl;
        abort();
    }

    if (blocktype == TensorBlock::data_block)
    {
        store(block);
    }
    else
    {
    }
    

    return block;
}

void
TensorBlock::convert(
    Matrix* matrix,
    MatrixConfiguration* config,
    TensorIndexDescr* descr
)
{
    branch_->get_node()->convert(matrix,config,descr);
}

void
TensorBlock::accumulate(
    Matrix* matrix,
    MatrixConfiguration* config,
    TensorIndexDescr* descr
)
{
    branch_->get_node()->accumulate(matrix,config,descr);
}

void
TensorBlock::element_op(ElementOp* op)
{
    branch_->get_node()->element_op(op);

    if (!in_read_mode())
        set_unsynced();
}

bool
TensorBlock::has_controller() const
{
    return controller_;
}

bool
TensorBlock::has_remote_controller() const
{
    bool check1 = read_controller_ == remote_read_controller_;
    bool check2 = read_controller_ == remote_sorted_read_controller_;
    return check1 && check2;
}

bool
TensorBlock::equals(const void** data)
{
    bool equals = branch_->get_node()->equals(data);   
    return equals;
}

bool
TensorBlock::equals(TensorBlock* block)
{
    return branch_->get_node()->equals(block->get_branch()->get_node());
}

void
TensorBlock::free_data(char* data, size_t size)
{
    YetiRuntime::free(data, size);
}

void
TensorBlock::free_metadata(char* data, size_t size)
{
    YetiRuntime::free(data, size);
}

StorageBlock*
TensorBlock::allocate_metadata_storage_block()
{
    return allocate_storage_block(parent_tensor_->metadata_block_size(), TensorBlock::metadata_block);

}

void
TensorBlock::set_element_size(size_t size)
{
}

void
TensorBlock::fill(TensorElementComputer* filler)
{
    if (!resort_perm_->is_identity()) //this is not a unique tile
        return;

    TemplateInfo::type_t element_type
        = filler->element_type(indices_, depth());

    branch_->set_element_type(element_type);
    branch_->get_node()->fill(filler);

    set_unsynced();

}

TensorIndexDescr*
TensorBlock::get_descr() const
{
    return parent_tensor_->get_block_descr();
}

TensorBranch*
TensorBlock::get_branch() const
{
    return branch_;
}

uli
TensorBlock::get_nelements_data() const
{
    return descr()->get_nelements_data(depth() + 1, get_indices());
}

uli
TensorBlock::get_block_number() const
{
    return parent_tensor_->get_index(indices_);
}

usi
TensorBlock::get_degeneracy() const
{
    return degeneracy_;
}

usi
TensorBlock::get_depth() const
{
    return depth();
}

const uli*
TensorBlock::get_indices() const
{
    return indices_;
}

std::string
TensorBlock::get_index_string() const
{
    return ClassOutput<const uli*>::str(nindex(), indices_);
}

std::string
TensorBlock::get_block_name() const
{
    std::string name = get_index_string();
    if (tmp_block_)
        name += " tmp block ";
    name += " on tensor ";
    std::string tensor_name = parent_tensor_->get_name();
    name += tensor_name;
    return name;
}

DataNode*
TensorBlock::get_first_data_node() const
{
    return branch_->get_node()->get_first_data_node();
}

MemoryPool*
TensorBlock::get_metadata_mempool() const
{
    return metadata_mempool_;
}

void
TensorBlock::free_metadata()
{
    if (metadata_block_ && metadata_block_->data())
        free_metadata(metadata_block_->data(), metadata_block_->size());
    clear_metadata();
}

void
TensorBlock::configure(Tensor::tensor_storage_t storage_type)
{
    storage_type_ = storage_type;
}

void
TensorBlock::configure(const TensorRetrieveActionPtr& action)
{
    configure(Tensor::action);
    action_ = action;
}

void
TensorBlock::clear_metadata()
{
    if (metadata_block_)
    {
        delete metadata_block_;
        metadata_block_ = 0;
    }
}

StorageBlock*
TensorBlock::get_metadata_storage_block() const
{
    return metadata_block_;
}

Tensor*
TensorBlock::get_parent_tensor() const
{
    return parent_tensor_;
}

Permutation*
TensorBlock::get_block_permutation() const
{
    return block_perm_;
}

Permutation*
TensorBlock::get_resort_permutation() const
{
    return resort_perm_;
}

Permutation*
TensorBlock::get_fetch_permutation() const
{
    return fetch_perm_;
}

TensorBlock*
TensorBlock::get_symmetry_unique_block() const
{
    if (unique_block_)
        return unique_block_;

    TensorBlock* block = parent_tensor_->get_block(indices_, fetch_perm_);
#if YETI_SANITY_CHECK
    if (permutationally_unique_ && block && block != this)
    {
        cerr << "permutationally unique block " << this->get_index_string()
             << " did not map to self -> "
             << block->get_index_string()
             << endl
             << fetch_perm_
             << endl
             << (void*) block << " " << (void*) this
             << endl;

        abort();
    }

    if (!block)
    {
        cerr << "block " << get_index_string()
             << " on tensor " << parent_tensor_->get_name()
             << " is not null, but parent block is!" << endl;
        cerr << fetch_perm_ << endl;
        abort();
    }
    else if (!block->is_permutationally_unique())
    {
        cerr << "fetched non-unique block" << endl;
        cerr << block->get_index_string() << endl;
        block->print(cerr); cerr << endl;
        cerr << "fetching from " << get_index_string() << endl;
        fetch_perm_->print(cerr); cerr << endl;
        abort();
    }
#endif
    return block;
}

TensorController*
TensorBlock::get_tensor_controller() const
{
    return controller_;
}

TensorRetrieveAction*
TensorBlock::get_retrieve_action() const
{
    return action_.get();
}

uli
TensorBlock::get_malloc_number() const
{
    return malloc_number_;
}

bool
TensorBlock::in_read_mode() const
{
    return controller_ == read_controller_;
}

bool
TensorBlock::in_accumulate_mode() const
{
    return controller_ == accumulate_controller_;
}

bool
TensorBlock::in_verbatim_mode() const
{
    return controller_ == verbatim_controller_;
}

bool
TensorBlock::in_write_mode() const
{
    return controller_ == write_controller_;
}

bool
TensorBlock::is_cached() const
{
    return metadata_block_->is_cached();
}

void
TensorBlock::init_mempool()
{
    metadata_mempool_ = new MemoryPool(metadata_block_->size(), metadata_block_->data());
}

void
TensorBlock::init_branch()
{
    /**
        Create the tensor branch in the memory pool
        owned by the tensor block
    */
    branch_ = new (metadata_mempool_)
        TensorBranch(
            indices_,
            depth(),
            this
        );
    branch_->set_descr(descr());
    branch_->set_malloc_number(malloc_number_);
}

void
TensorBlock::init_in_core()
{
    if (permutationally_unique_)
        init_in_core_no_sort();
    else
        init_in_core_resort();
}


void
TensorBlock::init_tmp_accumulate()
{
    char* metadata_data = YetiRuntime::malloc(metadata_block_size());
    metadata_block_ = new InCoreBlock(metadata_data, metadata_block_size());

    init_mempool();

    init_branch();

    read_controller_ = abort_read_controller_;

    accumulate_controller_ = tmp_thread_accumulate_controller_;

    write_controller_ = abort_write_controller_;

    verbatim_controller_ = abort_verbatim_controller_;

    controller_ = accumulate_controller_;

    set_initialized(true);

    finalized_ = true;
}

void
TensorBlock::init_in_core_no_sort()
{
    char* metadata_data = YetiRuntime::malloc(metadata_block_size());
    metadata_block_ = new InCoreBlock(metadata_data, metadata_block_size());
    init_mempool();
    init_branch();

    read_controller_ = in_core_controller_;

    write_controller_ = in_core_controller_;

    accumulate_controller_ = in_core_controller_;

    verbatim_controller_ = in_core_controller_;

    set_initialized(true);

    //set as retrieved
    fast_read_ = true;

    set_prefetched(true);
    set_finalized(false);
}

void
TensorBlock::init_in_core_resort()
{
    CachedStorageBlock* block = new CachedStorageBlock(this, parent_tensor_->get_metadata_cache());
    metadata_block_ = block;
    init_mempool();


    read_controller_ = in_core_sorted_read_controller_;

    write_controller_ = abort_write_controller_;

    accumulate_controller_ = in_core_sorted_accumulate_controller_;

    verbatim_controller_ = abort_verbatim_controller_;
}

void
TensorBlock::init_action()
{
    if (permutationally_unique_)
        init_action_no_sort();
    else
        init_in_core_resort();
}

void
TensorBlock::init_action_no_sort()
{
    CachedStorageBlock* block = new CachedStorageBlock(this, parent_tensor_->get_metadata_cache());
    metadata_block_ = block;
    init_mempool();

    read_controller_ = action_read_controller_;

    write_controller_ = abort_write_controller_;

    accumulate_controller_ = abort_accumulate_controller_;

    verbatim_controller_ = abort_verbatim_controller_;
}

void
TensorBlock::init_action_resort()
{
    CachedStorageBlock* block = new CachedStorageBlock(this, parent_tensor_->get_metadata_cache());
    metadata_block_ = block;
    init_mempool();

    read_controller_ = action_sorted_read_controller_;

    write_controller_ = abort_write_controller_;

    accumulate_controller_ = abort_accumulate_controller_;

    verbatim_controller_ = abort_verbatim_controller_;
}


void
TensorBlock::init_on_disk()
{
    if (permutationally_unique_)
        init_on_disk_no_sort();
    else
        init_on_disk_resort();
}

void
TensorBlock::init_on_disk_no_sort()
{
    yeti_throw(SanityCheckError, "no disk usage currently allowed");

    metadata_block_->retrieve();

    init_mempool();
    init_branch();

    read_controller_ = disk_read_controller_;
    
    write_controller_ = disk_write_controller_;

    accumulate_controller_ = disk_accumulate_controller_;

    verbatim_controller_ = abort_verbatim_controller_;

    metadata_block_->release();

}

void
TensorBlock::init_on_disk_resort()
{
   metadata_block_ = new CachedStorageBlock(this, parent_tensor_->get_metadata_cache());
   init_mempool();

   read_controller_ = disk_sorted_read_controller_;

    write_controller_ = abort_write_controller_;

    accumulate_controller_ = abort_accumulate_controller_;

    verbatim_controller_ = abort_verbatim_controller_;
}

void
TensorBlock::configure_tmp_block(TensorBlock* block)
{
    try{
        wait_on_send(); //very important...
    } catch (int e) {
        cerr << stream_printf("Wait on send failure for %s being reconfigured to %s on thread %d node %d\n",
                        get_block_name().c_str(), block->get_block_name().c_str(),
                        YetiRuntime::get_thread_number(), YetiRuntime::me());
        throw e;
    }
    for (uli i=0; i < nindex(); ++i)
        indices_[i] = block->indices_[i];

    malloc_number_ = block->malloc_number_;
    node_number_ = block->node_number_;
    branch_->set_malloc_number(malloc_number_);

    if (block->is_remote_block())
    {
        accumulate_controller_ = tmp_remote_accumulate_controller_;
    }
    else
    {
        accumulate_controller_ = tmp_thread_accumulate_controller_;
    }

    controller_ = accumulate_controller_;
}

bool
TensorBlock::is_up_to_date() const
{
    return up_to_date_;
}

bool
TensorBlock::is_subblock() const
{
    return is_subblock_;
}

bool
TensorBlock::is_flushable() const
{
    if (!metadata_block_->is_retrieved())
        return true;

    DataStorageNode* node = head_data_node();
    while(node)
    {
        StorageBlock* block = node->block;
        if (!block->is_retrieved())
            return true;
        node = node->next;
    }
    
    return false;
}

bool
TensorBlock::is_nonnull() const
{
    wait_on_remote(); //something else might be off working on initialize - which means this is technically nonnull
    return is_initialized();
}

void
TensorBlock::wait_on_remote() const
{
    while (remote_wait_)
        usleep(1);
}

void
TensorBlock::set_remote_wait(bool flag)
{
#if YETI_SANITY_CHECK
    if (!is_locked())
    {
        controller_fail();
        block_fail(SanityCheckError, "Block is not locked in set remote wait");
    }
#endif
    remote_wait_ = flag;
}

void
TensorBlock::set_task_owner(bool flag)
{
    task_owner_ = flag;
}


bool
TensorBlock::is_task_owner() const
{
    return task_owner_;
}

void
TensorBlock::set_subblock(bool flag)
{
    is_subblock_ = flag;
}

void
TensorBlock::init_remote()
{
    if (permutationally_unique_)
    {
        init_remote_no_sort();
    }
    else
    {
        init_remote_resort();
    }
}


void
TensorBlock::init_remote_no_sort()
{
    metadata_block_ = new CachedStorageBlock(this, parent_tensor_->get_metadata_cache());

    init_mempool();

    read_controller_ = remote_read_controller_;

    write_controller_ = abort_write_controller_;

    accumulate_controller_ = remote_accumulate_controller_;

    verbatim_controller_ = abort_verbatim_controller_;

}

void
TensorBlock::init_remote_resort()
{
    metadata_block_ = new CachedStorageBlock(this, parent_tensor_->get_metadata_cache());

    init_mempool();

    read_controller_ = remote_sorted_read_controller_;

    write_controller_ = abort_write_controller_;

    accumulate_controller_ = abort_accumulate_controller_;

    verbatim_controller_ = abort_verbatim_controller_;
}

void
TensorBlock::internal_contraction(
    TensorBlock* dst_block,
    MatrixConfiguration* config
)
{
    branch_->get_node()->internal_contraction(
        dst_block->get_branch()->get_node(),
        config
    );
}

void
TensorBlock::init_recomputed()
{
    if (permutationally_unique_)
        init_recomputed_no_sort();
    else
        init_recomputed_resort();
}

void
TensorBlock::init_recomputed_no_sort()
{
    metadata_block_ = new CachedStorageBlock(this, parent_tensor_->get_metadata_cache());
    init_mempool();

    read_controller_ = recompute_read_controller_;
    
    write_controller_ = abort_write_controller_;
    
    accumulate_controller_ = abort_accumulate_controller_;

    verbatim_controller_ = abort_verbatim_controller_;
}

void
TensorBlock::init_recomputed_resort()
{
    //no distinction is to be made between these
    init_recomputed_no_sort();
}

bool
TensorBlock::is_permutationally_unique() const
{
    return permutationally_unique_;
}

bool
TensorBlock::is_recomputed_block() const
{
    return storage_type_ == Tensor::recomputed;
}

bool
TensorBlock::is_local_block() const
{
    return YetiRuntime::me() == node_number_;
}

bool
TensorBlock::is_remote_block() const
{
    bool remote_block = YetiRuntime::me() != get_node_owner();
    return remote_block;
}

bool
TensorBlock::is_synced() const
{
    return is_synced_;
}

void
TensorBlock::set_synced()
{
    is_synced_ = true;
}

void
TensorBlock::set_unsynced()
{
    is_synced_ = false;
}

void
TensorBlock::_obsolete()
{
    if (!is_locked())
        block_fail(SanityCheckError, "Cannot obsolete unlocked tensor block");

#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("obsolete");
#endif

    wait_on_remote();

    //it's possible this is not symmetry unique and was never used
    if (controller_)
    {
#if DEBUG_TENSOR_BLOCK_RETRIEVE
        cout << stream_printf("Obsoleting %s on thread %ld node %ld\n", 
            get_block_name().c_str(), YetiRuntime::get_thread_number(), YetiRuntime::me());
        cout.flush();
#endif
        controller_->obsolete(this);
    }
}

void
TensorBlock::print(std::ostream& os)
{
    if (branch_ == 0)
        return;
    if (branch_->get_node() == 0)
        return;

    NormElementOp* op = new NormElementOp;
    element_op(op);
    double norm = op->norm();
    double maxlog = branch_->get_node()->get_max_log();

    set_synced();

    os << indexstr(descr()->nindex(), indices_);
    os << stream_printf(" Unique=%d    Malloc Numbers=(%d,%d)", 
                    permutationally_unique_, malloc_number_, branch_->get_malloc_number()) << endl;
    os   << "Branch: " << branch_ 
        << " Degeneracy: " << degeneracy_ 
        << stream_printf("norm: %8.4e  mempool: %p  size: %ld  remain: %ld", 
                         norm, metadata_mempool_,
                         metadata_mempool_->total_size(), metadata_mempool_->remaining())
        << endl;

    branch_->get_node()->print(os);
}


/*************
End public interface definition
*************/

/*************
Begin private interface definition
*************/

bool
TensorBlock::in_destructor() const
{
    return in_destructor_;
}

bool
TensorBlock::is_cache_coherent() const
{
    if (is_permutationally_unique())
        return true;

    bool coherent = metadata_block_->data() == 0;

    return coherent;
}


void
TensorBlock::complete_read()
{
    if (fast_read_)
        return;

    wait_on_remote(); //something else might be off working on initialize
    lock();

#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("complete read");
#endif

    set_read_mode();
    if (fast_read_) //this might have been set by initialization
    {
        unlock();
        return;
    }

#if YETI_SANITY_CHECK
    if (prefetch_count_ == 0 || !is_prefetched())
    {
        controller_fail();
        yeti_throw(SanityCheckError, "Tensor block is not prefetched. There is no read to complete");
    }
#endif
    --prefetch_count_; //decrement the prefetch count

    if (need_parent_retrieve())
        get_symmetry_unique_block()->complete_read();

    YetiRuntimeObject::retrieve();
    if (storage_type_ != Tensor::recomputed)
        unlock();
}

bool
TensorBlock::need_parent_retrieve() const
{
    return !permutationally_unique_ && storage_type_ != Tensor::recomputed;
}

void
TensorBlock::retrieve_read()
{
    if (fast_read_)
        return;

    //the wait on remote must happen outside the lock
    //thus we cannot lock or the usability restore
    //will deadlock
    wait_on_remote(); //something else might be off working on initialize
    lock();

#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("retrieve read");
#endif
    set_read_mode();
    if (fast_read_) //this might have been set by initialization
    {
        unlock();
        return;
    }

    if (need_parent_retrieve())
        get_symmetry_unique_block()->retrieve_read();

    YetiRuntimeObject::retrieve();
    if (storage_type_ != Tensor::recomputed)
        unlock();

    if (need_parent_retrieve())
        get_symmetry_unique_block()->release_read();
}

void
TensorBlock::release_read()
{
    if (fast_read_)
        return;


    if (storage_type_ == Tensor::recomputed)
    {
        #if TRACK_TENSOR_BLOCK_HISTORY
            add_event("release read");
        #endif
        YetiRuntimeObject::release();
        unlock();
    }
    else
    {
        lock();

        if (!permutationally_unique_)
            get_symmetry_unique_block()->release_read();

        #if TRACK_TENSOR_BLOCK_HISTORY
            add_event("release read");
        #endif
        YetiRuntimeObject::release();
        unlock();
    }

}

void
TensorBlock::retrieve_accumulate_no_lock()
{
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("retrieve accumulate");
#endif

    set_accumulate_mode();
    YetiRuntimeObject::retrieve();
}

void
TensorBlock::retrieve_accumulate()
{
    wait_on_remote(); //something else might be off working on initialize
    lock();
    retrieve_accumulate_no_lock();
}

void
TensorBlock::release_accumulate()
{
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("release accumulate");
#endif
    YetiRuntimeObject::release();
    unlock();
}

void
TensorBlock::wtf_retrieve()
{
    //wait_on_usable();
    //lock();
    //set_verbatim_mode();
    //YetiRuntimeObject::retrieve();
    //wait_on_send();
    //_retrieve();
}

void
TensorBlock::wtf_release()
{
    //_release();
    //YetiRuntimeObject::release();
    //unlock();
}

void
TensorBlock::retrieve_verbatim_no_lock()
{
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("retrieve verbatim");
#endif
    set_verbatim_mode();
    YetiRuntimeObject::retrieve();
}

void
TensorBlock::retrieve_verbatim()
{
    wait_on_remote(); //something else might be off working on initialize
    lock();
    retrieve_verbatim_no_lock();
}

void
TensorBlock::release_verbatim()
{
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("release verbatim");
#endif
    YetiRuntimeObject::release();
    unlock();
}

void
TensorBlock::retrieve_write()
{
    wait_on_remote(); //something else might be off working on initialize
    lock();
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("retrieve write");
#endif
    set_write_mode();
    YetiRuntimeObject::retrieve();
}

void
TensorBlock::release_write()
{
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("release write");
#endif
    YetiRuntimeObject::release();
    unlock();
}

void
TensorBlock::finalize()
{
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("finalize");
#endif
    if (finalized_)
    {
        set_prefetched(false);
        return;
    }
    
#if YETI_SANITY_CHECK
    if (!is_locked())
        block_fail(SanityCheckError, "cannot finalize unlocked object");
#endif
#if DEBUG_TENSOR_BLOCK_RETRIEVE
    cout << stream_printf("Finalizing %s on thread %ld node %ld\n", 
        get_block_name().c_str(), YetiRuntime::get_thread_number(), YetiRuntime::me());
    cout.flush();
#endif
    controller_->finalize(this);
    YetiRuntimeObject::set_finalized(true);
    set_prefetched(false);
}

void
TensorBlock::_initialize()
{
    if (!controller_)
        block_fail(SanityCheckError, "tensor block not yet initialized");

#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("initialize");
#endif

#if DEBUG_TENSOR_BLOCK_RETRIEVE
    cout << stream_printf("Initialize %s on thread %ld node %ld\n", 
        get_block_name().c_str(), YetiRuntime::get_thread_number(), YetiRuntime::me());
    cout.flush();
#endif

#if YETI_SANITY_CHECK
    //cache block should not already be retrieved
    if (metadata_block_->is_cached() && metadata_block_->is_retrieved())
        block_fail(SanityCheckError, "metadata block is retrieved already in initialize");
#endif

#if TRACK_TENSOR_BLOCK_HISTORY
    if (!permutationally_unique_)
    {
        get_symmetry_unique_block()->add_event(this);
        get_symmetry_unique_block()->add_event("non-unique child block initialize");
    }
#endif

	try{
		controller_->validate(this);
		bool in_core = metadata_block_->retrieve();
#if TRACK_TENSOR_BLOCK_HISTORY
        add_event("metadata retrieve");
#endif

		branch_ = reinterpret_cast<TensorBranch*>(metadata_block_->data());
		metadata_mempool_->set(metadata_block_->data());
		branch_->set_metadata_mempool(metadata_mempool_);
		branch_->set_descr(descr());
		branch_->set_parent(this);

		DataStorageNode* node = head_data_node(); 
		while(node)
		{
			StorageBlock* block = node->block;
			block->retrieve();
			node = node->next;
#if TRACK_TENSOR_BLOCK_HISTORY
            add_event("data retrieve");
#endif
		}
	} catch (int e) {
		if (e == TENSOR_BLOCK_POLICY_EXCEPTION)
			controller_fail();
		throw e;
	}
}

void
TensorBlock::_renew()
{
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("renew");
    if (!permutationally_unique_)
    {
        get_symmetry_unique_block()->add_event(this);
        get_symmetry_unique_block()->add_event("non-unique child block renew");
    }
#endif

    wait_on_send();
#if YETI_SANITY_CHECK
    if (metadata_block_->is_retrieved())
    {
        
        DataStorageNode* node = head_data_node();
        while (node)
        {
            if (!node->block->is_retrieved())
            {
                cerr << stream_printf("Tensor block %s retrieve failure\n", get_block_name().c_str());
                block_fail(SanityCheckError, "Data block is not retrieved at but metadata block is retrieved");
            }
            node = node->next;
        }
    }
#endif

  try{
#if DEBUG_TENSOR_BLOCK_RETRIEVE
    cout << stream_printf("Renew %s on thread %ld node %ld\n", 
        get_block_name().c_str(), YetiRuntime::get_thread_number(), YetiRuntime::me());
    cout.flush();
#endif


    controller_->renew(this);
    set_synced();

#if YETI_SANITY_CHECK
    if (!metadata_block_->is_retrieved())
    {
        controller_fail();
        block_fail(SanityCheckError, "Metadata block is not retrieved at end of renew");
    }

    DataStorageNode* node = head_data_node();
    while (node)
    {
        if (!node->block->is_retrieved())
            block_fail(SanityCheckError, "Data block is not retrieved at end of renew");
        node = node->next;
    }
#endif
  } catch (int e) {
    if (e == TENSOR_BLOCK_POLICY_EXCEPTION)
        controller_fail();
    throw e;
  }
}

void
TensorBlock::_retrieve()
{
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("retrieve");
    if (!permutationally_unique_)
    {
        get_symmetry_unique_block()->add_event(this);
        get_symmetry_unique_block()->add_event("non-unique child block retrieve");
    }
#endif

    wait_on_send();

#if DEBUG_TENSOR_BLOCK_RETRIEVE
    cout << stream_printf("Retrieve %s on thread %ld node %ld\n", 
        get_block_name().c_str(), YetiRuntime::get_thread_number(), YetiRuntime::me());
    cout.flush();
#endif


#if YETI_SANITY_CHECK
    if (metadata_block_->is_retrieved())
    {
        
        DataStorageNode* node = head_data_node();
        while (node)
        {
            if (!node->block->is_retrieved())
            {
                cerr << stream_printf("Tensor block %s retrieve failure\n", get_block_name().c_str());
                block_fail(SanityCheckError, "Data block is not retrieved at but metadata block is retrieved");
            }
            node = node->next;
        }
    }
#endif

  try{
    bool in_core = metadata_block_->retrieve();
#if TRACK_TENSOR_BLOCK_HISTORY
        add_event("metadata retrieve");
#endif
    if (!in_core)
        block_fail(SanityCheckError, "metadata block is not in core despite having entered _retrieve");

    if (!permutationally_unique_ && controller_ == read_controller_)
    {
        TensorBlock* block = get_symmetry_unique_block();
        block->retrieve_read(); //ensure retrieved throughout entire retrieve
        //reset the send status after the retrieve
        controller_->retrieve(this);
        block->release_read();
#if TRACK_TENSOR_BLOCK_HISTORY
        block->add_event("parent block retrieve");
#endif
    }
    else
    {
        //reset the send status after the retrieve
        controller_->retrieve(this);
    }

    //after retrieving, we are synced with the parent storage location
    set_synced();
    //we are also prefetched
    set_prefetched(true);

#if YETI_SANITY_CHECK
    if (!metadata_block_->is_retrieved())
    {
        block_fail(SanityCheckError, "Metadata block is not retrieved at end of retrieve");
    }

    DataStorageNode* node = head_data_node();
    while (node)
    {
        if (!node->block->is_retrieved())
        {
            cerr << stream_printf("Tensor block %s retrieve failure\n", get_block_name().c_str());
            block_fail(SanityCheckError, "Data block is not retrieved at end of retrieve");
        }
        node = node->next;
    }
#endif
  } catch (int e)
  {
    if (e == TENSOR_BLOCK_POLICY_EXCEPTION)
        controller_fail();
    throw e;
  }

}

const char*
TensorBlock::controller_fail(TensorController* controller) const
{
#define _controller_fail(x) if (controller == x) return #x;
    _controller_fail(in_core_controller_);
    _controller_fail(in_core_sorted_write_controller_);
    _controller_fail(in_core_sorted_read_controller_);
    _controller_fail(in_core_sorted_verbatim_controller_);
    _controller_fail(in_core_sorted_accumulate_controller_);
    _controller_fail(disk_read_controller_);
    _controller_fail(disk_write_controller_);
    _controller_fail(disk_accumulate_controller_);
    _controller_fail(disk_sorted_read_controller_);
    _controller_fail(remote_read_controller_);
    _controller_fail(remote_write_controller_);
    _controller_fail(remote_verbatim_controller_);
    _controller_fail(remote_accumulate_controller_);
    _controller_fail(remote_sorted_read_controller_);
    _controller_fail(action_read_controller_);
    _controller_fail(action_sorted_read_controller_);
    _controller_fail(abort_read_controller_);
    _controller_fail(abort_write_controller_);
    _controller_fail(abort_accumulate_controller_);
    _controller_fail(abort_verbatim_controller_);
    _controller_fail(recompute_read_controller_);
    _controller_fail(tmp_thread_accumulate_controller_);
    _controller_fail(tmp_remote_accumulate_controller_);
    return "";
}

void
TensorBlock::controller_fail() const
{
#if TRACK_TENSOR_BLOCK_HISTORY
    print_history();
#endif
    YetiRuntime::stack_print();

    if (tmp_block_)
        cout << stream_printf("Controller fail on tmp block %s on thread %d on node %d\n", 
                    get_block_name().c_str(), YetiRuntime::get_thread_number(), YetiRuntime::me());
    else
        cout << stream_printf("Controller fail on block %s on thread %d on node %d\n", 
                    get_block_name().c_str(), YetiRuntime::get_thread_number(), YetiRuntime::me());
    
    cerr << stream_printf("Read: %s\nWrite: %s\nVerbatim: %s\nAccumulate: %s\n",
                            controller_fail(read_controller_),
                            controller_fail(write_controller_),
                            controller_fail(verbatim_controller_),
                            controller_fail(accumulate_controller_)
                          );

}

#if TRACK_TENSOR_BLOCK_HISTORY
void
DataStorageBlockAllocator::print_history() const
{

    std::list<const char*>::const_iterator it = history_.begin();
    std::list<const char*>::const_iterator stop = history_.end();
    cerr << "History:\n";
    for ( ; it != stop; ++it)
    {
        const char* msg = *it;
        cerr << stream_printf("Event on node %d: %s\n", YetiRuntime::me(), msg);
    }
    cerr << endl;
}

void
DataStorageBlockAllocator::add_event(const char* event)
{
    char* str = new char[100];
    sprintf(str, "%s on thread %ld on node %ld", event, 
        YetiRuntime::get_thread_number(),
        YetiRuntime::me());
    history_.push_back(str);
}

void
DataStorageBlockAllocator::add_event(TensorBlock* block)
{
    char* str = new char[200];
    sprintf(str, "%s", block->get_block_name().c_str());
    history_.push_back(str);
}
#endif

void
TensorBlock::_release()
{
    if (prefetch_count_) //there are more blocks waiting
        return; //don't release

#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("release");
    if (!permutationally_unique_)
    {
        get_symmetry_unique_block()->add_event(this);
        get_symmetry_unique_block()->add_event("non-unique child block release");
    }
#endif

#if DEBUG_TENSOR_BLOCK_RETRIEVE
    cout << stream_printf("Release %s on thread %ld node %ld\n", 
        get_block_name().c_str(), YetiRuntime::get_thread_number(), YetiRuntime::me());
    cout.flush();
#endif

    controller_->release(this);

}
/*************
End private interface definition
*************/

/**********
Begin Cachable interface definition
**********/
void
TensorBlock::flush_from_cache()
{
    if (!controller_) //this was never actually used for anything
    {
        set_initialized(false);
        set_prefetched(false);
        return;
    }

#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("flush from cache");
#endif

#if DEBUG_TENSOR_BLOCK_RETRIEVE
    cout << stream_printf("Flushing %s on node %ld thread %ld\n", get_block_name().c_str(), YetiRuntime::me(), YetiRuntime::get_thread_number());
    cout.flush();
#endif

#if YETI_SANITY_CHECK
    if (metadata_block_->is_retrieved())
    {
        cerr << stream_printf("Block %s was not properly release before calling flush\n",
                        get_block_name().c_str());
        controller_fail();
        block_fail(SanityCheckError, "Cannot flush retrieved block");
    }
#endif

    wait_on_send(); 
    if (!finalized_)
    {
        finalize(); //might generate a new send
        wait_on_send();
    }

#if YETI_SANITY_CHECK
    if (!is_locked())
        block_fail(SanityCheckError, "Cannot flush unlocked tensor block");
#endif

    try{
        controller_->flush(this);
    } catch (int e) {
        cerr << stream_printf("Error in metadata flush of block %s", get_block_name().c_str());
        throw e;
    }
    try{
        controller_->flush_data(this);
    } catch (int e) {
        cerr << stream_printf("Error in data flush of block %s", get_block_name().c_str());
        throw e;
    }


#if YETI_SANITY_CHECK
    CachedStorageBlock* block = dynamic_cast<CachedStorageBlock*>(metadata_block_);
    if (block)
    {
        if (block->get_entry())
        {
            block_fail(SanityCheckError, "Cache block not cleared in flush");
        }
        if (block->data())
        {
            block_fail(SanityCheckError, "Block not cleared in flush");
        }
    }
#endif

    set_initialized(false);
    set_prefetched(false);

}


void
TensorBlock::recompute_tensor_permutation()
{
    if (permutationally_unique_)
    {
        resort_perm_ = parent_tensor_->get_tensor_grp()->get_identity();
        fetch_perm_ = resort_perm_;
        return;
    }

    fetch_perm_ = parent_tensor_->get_fetch_permutation(indices_);
    resort_perm_ = fetch_perm_->inverse();

}

void
TensorBlock::sort_tensor_permutation(Permutation* p)
{
    if (permutationally_unique_)
        return;

    fetch_perm_ = p->product(fetch_perm_->product(p->inverse()));
    resort_perm_ = p->product(resort_perm_->product(p->inverse()));
    block_perm_ = p->product(block_perm_);
}

void
TensorBlock::sort_resort_permutation(Permutation* p)
{
    if (permutationally_unique_)
    {
        cerr << "cannot map resort permutation of unique block!" << endl;
        abort();
    }
    resort_perm_ = resort_perm_->product(p);
}

void
TensorBlock::set_node_number(uli node)
{
    node_number_ = node;
}

uli
TensorBlock::get_node_number() const
{
    return node_number_;
}

uli
TensorBlock::get_node_owner() const
{
    return node_number_ == NODE_NUMBER_UNINITIALIZED ? YetiRuntime::me() : node_number_;
}

void
TensorBlock::set_degeneracy(usi degeneracy)
{
    degeneracy_ = degeneracy;
}

bool
TensorBlock::set_accumulate_mode()
{
    if (!controller_)
        init();

    //already in accumulate mode
    if (controller_ == accumulate_controller_)
        return true;

    set_prefetched(false);

#if YETI_SANITY_CHECK
    if (is_retrieved())
        block_fail(SanityCheckError, "cannot set accumulate mode while tensor block is retrieved");
#endif
    
    if (!is_synced())
    {
        cerr << stream_printf("%p %p %p %p %p",
                             controller_,
                             read_controller_,
                              write_controller_,
                              accumulate_controller_,
                              verbatim_controller_) << endl;

        cerr << "Cannot switch to accumulate mode.  Tensor block "
             << get_index_string() << " in tensor " << parent_tensor_->get_name()
             << " has not been synced" << endl;
        abort();
    }
    controller_ = accumulate_controller_;
    return true;
}

void
TensorBlock::set_verbatim_mode()
{    
    if (!controller_)
        init();

    //already in correct mode
    if (controller_ == verbatim_controller_)
        return;

    set_prefetched(false);
    if (!is_synced())
    {
        cerr << "Cannot switch to verbatim mode.  Tensor block has"
             << " not been synced from previous operation" << endl;
        cerr << ClassOutput<const uli*>
                ::str(descr()->nindex(), indices_) << endl;
        abort();
    }

    controller_ = verbatim_controller_;
}

void
TensorBlock::set_read_mode()
{
    if (!controller_)
    {
        init();
    }

    //already in correct mode
    if (controller_ == read_controller_)
        return;

    set_prefetched(false);
    
    if (!is_synced())
    {
        cerr << "Cannot switch to read mode.  Tensor block " << get_index_string() 
            << " on node " << YetiRuntime::me() << " has "
             << " not been synced from previous operation" << endl;
        abort();
    }
    controller_ = read_controller_;
}

void
TensorBlock::set_write_mode()
{
    if (!controller_)
        init();

    //already in correct mode
    if (controller_ == write_controller_)
        return;
    
    set_prefetched(false);

    if (!is_synced())
    {
        cerr << "Cannot switch to write mode.  Tensor block has"
             << " not been synced from previous operation" << endl;
        cerr << ClassOutput<const uli*>
                ::str(descr()->nindex(), indices_) << endl;
        abort();
    }
    controller_ = write_controller_;
}

void
TensorBlock::metadata_sort(Sort* sort)
{
    uli* tmp = reinterpret_cast<uli*>(sort->metadata_buffer(0));
    sort->get_permutation()->permute<uli>(indices_, tmp);
    ::memcpy(indices_, tmp, nindex() * sizeof(uli));

    sort_tensor_permutation(sort->get_permutation());
}

void
TensorBlock::data_sort(Sort* sort)
{
    block_perm_ = sort->get_permutation()->product(block_perm_);

    if (!permutationally_unique_)
    {
        cerr << "Cannot resort non-permutationally unique block" << endl;
        abort();
    }

    branch_->get_node()->sort(sort);
}

void
TensorBlock::sync()
{
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("sync");
#endif
    if (is_synced())
        return;

#if DEBUG_TENSOR_BLOCK_RETRIEVE
    cout << stream_printf("Syncing %s\n", get_block_name().c_str());
    cout.flush();
#endif
    controller_->sync(this);
    set_synced();
}

void
TensorBlock::update()
{
    reset_degeneracy();

    if (is_permutationally_unique() && is_local_block())
    {
        branch_->get_node()->update();
    }

    up_to_date_ = true;
}

void
TensorBlock::reset_degeneracy()
{
    degeneracy_ = 1;
}

void
TensorBlock::set_prefetched(bool flag)
{
    is_prefetched_ = flag;
}

bool
TensorBlock::is_prefetched() const
{
    return is_prefetched_;
}

bool
TensorBlock::has_fast_read() const
{
    return fast_read_;
}

void
TensorBlock::prefetch_read()
{
    if (fast_read_)
        return;

    lock();
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("prefetch read");
    if (!permutationally_unique_)
    {
        get_symmetry_unique_block()->add_event(this);
        get_symmetry_unique_block()->add_event("non-unique child block prefetch");
    }
#endif
    set_read_mode(); //always set the mode first as this may change prefetch state
    if (fast_read_)
    {
        unlock();
        return;
    }

    if (need_parent_retrieve())
        get_symmetry_unique_block()->prefetch_read();

    ++prefetch_count_;
    if (is_prefetched())
    {
#if TRACK_TENSOR_BLOCK_HISTORY
        add_event("skip prefetch");
#endif
        if (prefetch_count_ == 1 && !metadata_block_->is_retrieved())
        //there are no pending prefetchs - make sure all blocks are retrieved
        {
            metadata_block_->retrieve();
            #if TRACK_TENSOR_BLOCK_HISTORY
            add_event("metadata retrieve");
            #endif
            DataStorageNode* node = head_data_node(); 
            while(node)
            {
                StorageBlock* block = node->block;
                block->retrieve();
                node = node->next;
                #if TRACK_TENSOR_BLOCK_HISTORY
                add_event("data retrieve");
                #endif
            }
        }
        unlock();
        return;
    }

#if YETI_SANITY_CHECK
    if (!finalized_)
    {
        cerr << stream_printf("Invalid unfinalized prefetch state %s on thread %ld node %ld\n", 
            get_block_name().c_str(), YetiRuntime::get_thread_number(), YetiRuntime::me());
        cerr.flush();
        controller_fail();
        block_fail(SanityCheckError, "Invalid prefetch");
    }
    if (is_retrieved())
    {
        cerr << stream_printf("Invalid retrieved prefetch state %s on thread %ld node %ld\n", 
            get_block_name().c_str(), YetiRuntime::get_thread_number(), YetiRuntime::me());
        cerr.flush();
        controller_fail();
        block_fail(SanityCheckError, "Invalid prefetch");
    }
#endif

#if DEBUG_TENSOR_BLOCK_RETRIEVE
    cout << stream_printf("Prefetch read %s on thread %ld node %ld\n", 
            get_block_name().c_str(), YetiRuntime::get_thread_number(), YetiRuntime::me());
    cout.flush();
#endif

    //do NOT unlock - it us up to the controller to unlock
    controller_->out_of_core_prefetch(this);
}

void
TensorBlock::prefetch_accumulate()
{
    return; //right now... no prefetch accumulates
}

void
TensorBlock::prefetch_read(TensorBlock* prev_block)
{
    lock();
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("prefetch in core read");
#endif
#if DEBUG_TENSOR_BLOCK_RETRIEVE
    cout << stream_printf("Prefetch in core read %s on thread %ld node %ld\n", 
            get_block_name().c_str(), YetiRuntime::get_thread_number(), YetiRuntime::me());
    cout.flush();
#endif
    set_read_mode();
    
    /** Here it us up to the controller to ensure metadata and data are retrieved */
    controller_->in_core_prefetch(this, prev_block);

#if YETI_SANITY_CHECK
    CachedStorageBlock* meta_block = dynamic_cast<CachedStorageBlock*>(metadata_block_);
    if (meta_block)
    {
        bool meta_retrieved = meta_block->is_retrieved();
        DataStorageNode* node = head_data_node();
        while (node)
        {
            CachedStorageBlock* block = dynamic_cast<CachedStorageBlock*>(node->block);
            bool block_retrieved = block->is_retrieved();
            if (block_retrieved != meta_retrieved)
            {
                cerr << stream_printf("Tensor block %s in core prefetch failure\n", get_block_name().c_str());
                block_fail(SanityCheckError, "Metadata and data retrieve statuses differ");
            }
            node = node->next;
        }
    }
#endif

    unlock();
}

void
TensorBlock::prefetch_accumulate(TensorBlock* prev_block)
{
    lock();
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("prefetch in core accumulate");
#endif

#if DEBUG_TENSOR_BLOCK_RETRIEVE
    cout << stream_printf("Prefetch in core accumulate %s on thread %ld node %ld\n", 
            get_block_name().c_str(), YetiRuntime::get_thread_number(), YetiRuntime::me());
    cout.flush();
#endif
    set_accumulate_mode();
    YetiRuntimeObject::initialize();
    controller_->in_core_prefetch(this, prev_block);
    unlock();
}

void
TensorBlock::redistribute(
    YetiMessenger* messenger,
    Message::message_data_type_t type,
    uli proc_number,
    size_t metadata_size,
    uli initial_tag,
    uli nmessages
)
{
#if DEBUG_TENSOR_BLOCK_RETRIEVE 
    cout << stream_printf("Receiving redistribute data block %s on %ld\n",
                          get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
    lock();
    uninit();
    set_node_number(YetiRuntime::me());
    init_local();
    controller_ = read_controller_;
    //treat this as a retrieve - this means the block is locked and waiting when the message arrives
    recv_data(
        messenger, 
        type, 
        Message::ReadRequestedRetrieve,
        proc_number, 
        metadata_size, 
        initial_tag, 
        nmessages
    );
    unlock();
}

bool
TensorBlock::is_waiting_on_redistribute() const
{
    if (!controller_)
        return true;

    bool remote_controller = read_controller_ == remote_read_controller_ || read_controller_ == remote_sorted_read_controller_;

    return remote_controller;
}

void
TensorBlock::enter_recv(Message::message_action_t action)
{
    if (action == Message::GlobalSum)
        lock();
}

void
TensorBlock::exit_recv(Message::message_action_t action)
{
    set_remote_wait(false);
    if (action != Message::ReadRequestedRetrieve)
        unlock();
}


void
TensorBlock::recv_data(
    YetiMessenger* messenger,
    Message::message_data_type_t type,
    Message::message_action_t action,
    uli proc_number,
    size_t metadata_size,
    uli initial_tag,
    uli nmessages
)
{
    enter_recv(action);

#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("recv data");
#endif

#if YETI_SANITY_CHECK
    if (!is_locked())
        block_fail(SanityCheckError, "Block is not locked in recv_data");
    /** not really sure of the best way to fix this
    it can happen that a block sends itself to a remote node
    and before it has had a chance to mark itself as released
    the data comes back from the parent in a global sum.  this is technically
    an error as a retrieved object cannot receive. just put a
    small wait */
    uli ncheck = 1;
    while (is_retrieved())
    {
        usleep(1);
        ++ncheck;

        if (ncheck > MICROSECOND_LOCK_WAIT)
            block_fail(SanityCheckError, "block is retrieved and cannot receive data after 1 second");
    }
    if (type != Message::TensorBranch)
        block_fail(SanityCheckError, "Tensor block can only receive branch data type");
#endif

        
    if (initial_tag == 0)
    {
#if TRACK_TENSOR_BLOCK_HISTORY
        add_event("skip receive branch");
#endif
#if YETI_SANITY_CHECK
        if (nmessages != 0)
            block_fail(SanityCheckError, "initial tag is zero but nmessages > 1");
#endif
        if (is_initialized())
        {
            //zero the branch 
            ZeroBranchRenew zeroer;
            zeroer.renew(this);
        }
        else; //do nothing - leave as a zero branch
        
        exit_recv(action);
        return;
    }

    if (!is_initialized())
    {
        init(); 
        //set up read mode
        controller_ = read_controller_;
        YetiRuntimeObject::initialize();
    }

    try{
        receive_branch(messenger, proc_number, metadata_size, initial_tag, nmessages);
    } catch (int e) {
        cerr << "error in TensorBlock::recv_data\n";
        throw e;
    }

    exit_recv(action);
}

void
TensorBlock::accumulate_data(
    YetiMessenger* messenger,
    Message::message_data_type_t type,
    uli proc_number,
    size_t metadata_size,
    uli initial_tag,
    uli nmessages
)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("accumulate remote data");
#endif
    if (initial_tag == 0)
    {
#if YETI_SANITY_CHECK
        if (nmessages != 0)
            block_fail(SanityCheckError, "initial tag is zero but nmessages > 0");
#endif
        return; //nothing got sent
    }

    if (type != Message::TensorBranch)
        block_fail(SanityCheckError, "tensor block can only receive branch data type");

    //allocate a temporary block
    TensorBlock* tmp_block = new (tensor_block_tmp_malloc) TensorBlock(this);
    tmp_block->retrieve_read();
    try{
        tmp_block->receive_branch(messenger, proc_number, metadata_size, initial_tag, nmessages);
    }
    catch(int e) {
        cerr << "error in TensorBlock::accumulate_data\n";
        throw e;
    }

    if (!GlobalQueue::get_task_queue()->worker_threads_active())
    {
        Sort* sort = 0;
        retrieve_accumulate();
        try{
            accumulate(tmp_block, 1.0, 0);
        } catch (int e) {
            cerr << stream_printf("Accumulate tmp block exception in TensorBlock::accumulate_data"
                                  " on node %ld on thread %ld for malloc number %ld\n",
                                  YetiRuntime::me(), YetiRuntime::get_thread_number(),
                                  malloc_number_);
            cerr.flush();
            throw e;
        }
        release_accumulate();
        tmp_block->release_read();
        tmp_block->lock();
        tmp_block->TensorBlock::~TensorBlock();
        tensor_block_tmp_malloc.free(tmp_block);
    }
    else
    {
        TensorBlockAccumulateTask* task = new TensorBlockAccumulateTask;

        task->tmp_block = tmp_block;
        task->parent_block = this;
        GlobalQueue::get_task_queue()->add_unexpected_task(task);
    }
}

void
TensorBlock::receive_branch(
    YetiMessenger* messenger,
    uli proc_number,
    size_t metadata_size,
    uli initial_tag,
    uli nmessages
)
{        
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("receive remote branch");
#endif
    YetiRuntime::start_timer("receive_branch");
#if YETI_SANITY_CHECK
    if (!is_initialized())
    {
        cerr << stream_printf("Block number %d is not initialized in receive branch\n", malloc_number_);
        controller_fail();
        throw TENSOR_BLOCK_COMMUNICATION_EXCEPTION;
    }
#endif
    if (metadata_size > metadata_block_->size())
    {
        std::string msg = stream_printf("Metadata block %s has size %ld and is not"
                                        "large enough to receive %ld bytes",
                                        get_block_name().c_str(),
                                        metadata_block_->size(),
                                        metadata_size);
        cerr << msg << endl;
        abort();
    }

    messenger->lock();

    uli tag = initial_tag;
#if TRACK_TENSOR_BLOCK_HISTORY
    add_event("recv tensor branch");
#endif
    messenger->recv_no_lock(
        proc_number,
        tag, 
        branch_,
        metadata_size
    );
    //parent block gets overwritten

    TensorBranch* obsolete_branch = branch_->get_last_branch_location();
    if (branch_ != obsolete_branch)
        branch_->realign_memory_pool(obsolete_branch, branch_);

    branch_->set_parent(this);
    branch_->set_descr(descr());
    branch_->set_metadata_mempool(metadata_mempool_);
    metadata_mempool_->set_data_size(metadata_size);

    TensorDataController* controller = branch_->first_data_controller();
    DataStorageNode* dnode = head_node_;
    while(controller)
    {
        size_t data_size = controller->get_data_size();

        StorageBlock* data_block = 0;
        if (dnode)
        {
            data_block = dnode->block;
            dnode = dnode->next;
        }
        else
        {
#if YETI_SANITY_CHECK
            if (data_size > parent_tensor_->get_data_cache()->blocksize())
            {
                cerr << stream_printf("Remote receive got invalid block size %ld "
                                     "for block %s.  Max block size is %ld.\n"
                                     "Trying to receive initial tag %ld for %ld messages.\n"
                                     "Metadata recv size is %ld. Actual metadata "
                                     "size is %ld\n",
                                    data_size, 
                                    get_block_name().c_str(),
                                    parent_tensor_->get_data_cache()->blocksize(),
                                    initial_tag, nmessages,
                                    metadata_size,
                                    parent_tensor_->get_metadata_cache()->blocksize()
                                );
                controller_fail();
                abort();
            }
#endif
            try {
                data_block = allocate_data_storage_block(data_size);
            } catch (int e) {
                cerr << "Error allocating cache block in receive branch\n";
                cerr.flush();
                throw e;
            }
        }

        controller->set_data(data_block->data());

        if (data_size > 0) //there is a message to receive
        {
            ++tag;
#if TRACK_TENSOR_BLOCK_HISTORY
            add_event("recv data block");
#endif
            messenger->recv_no_lock(
                proc_number,
                tag,
                data_block->data(),
                data_size
            );
        }
        controller = controller->next;
    }
    messenger->unlock();

    /** if the number of data blocks exceeds the number of controllers */
    while(dnode)
    {
        StorageBlock* data_block = dnode->block;
        branch_->allocate_data_controller(data_block);
        dnode = dnode->next;
    }
    


#if YETI_SANITY_CHECK
    uli ncheck = branch_->ncontrollers_nonzero() + 1;
    if (ncheck != nmessages)
    {
        std::string str = stream_printf("Wrong number of messages (%ld sent, should be %ld) on"
                                 " tensor block %s",
                                 ncheck,
                                 nmessages,
                                 get_block_name().c_str());
        block_fail(SanityCheckError, str);
    }

    if (branch_->get_malloc_number() != malloc_number_)
    {
            cerr << stream_printf("Received the wrong tensor branch on node %d thread %d, you twat.  Asked for %d and got %d\n"
                                  "My malloc number is %ld and the branch malloc number is %ld\n",
                                  YetiRuntime::me(), YetiRuntime::get_thread_number(),
                                  malloc_number_, branch_->get_malloc_number(),
                                  malloc_number_, branch_->get_malloc_number());
            cerr.flush();
            YetiRuntime::stack_print();
            YetiRuntime::get_messenger()->lock();
            YetiRuntime::get_messenger()->send_data_header(
                proc_number,
                GenericDataHeader,
                Message::TensorBranch,
                Message::Print,
                malloc_number_,
                0,
                0,
                0
            );
            if (branch_->get_malloc_number() != MALLOC_NUMBER_UNINITIALIZED)
            {
                YetiRuntime::get_messenger()->send_data_header(
                    proc_number,
                    GenericDataHeader,
                    Message::TensorBranch,
                    Message::Print,
                    branch_->get_malloc_number(),
                    0,
                    0,
                    0
                );
            }
            sleep(3);
            abort();
    }

#if 0
    const uli* mdindices = branch_->get_node()->get_indices();
    usi nindex = descr()->nindex();
    for (uli i=0; i < nindex; ++i)
    {
        if (indices_[i] != mdindices[i])
        {
            std::string str1 = indexstr(nindex, indices_);
            std::string str2 = indexstr(nindex, mdindices);
            cerr << stream_printf("Received the wrong tensor branch on node %d thread %d, you twat.  Asked for %s and got %s\n"
                                  "My malloc number is %ld and the branch malloc number is %ld\n",
                                  YetiRuntime::me(), YetiRuntime::get_thread_number(),
                                  str1.c_str(), str2.c_str(),
                                  malloc_number_, branch_->get_malloc_number());
            cerr.flush();
            YetiRuntime::stack_print();
            YetiRuntime::get_messenger()->lock();
            YetiRuntime::get_messenger()->send_data_header(
                proc_number,
                GenericDataHeader,
                Message::TensorBranch,
                Message::Print,
                malloc_number_,
                0,
                0,
                0
            );
            sleep(3);
            abort();
        }
    }
#endif

#endif

    YetiRuntime::stop_timer("receive_branch");
}


SendStatus*
TensorBlock::_send_data(
    YetiMessenger* messenger,
    Message::message_data_type_t type,
    Message::message_action_t action,
    uli proc_number
)
{
    uli header_tag = GenericDataHeader;
    if (action == Message::ReadRequestedRetrieve || action == Message::ReadRequestedPrefetch)
        header_tag = RequestedDataHeader;

    YetiRuntime::start_timer("send data");
    if (type != Message::TensorBranch)
    {
        Env::errn() << "Tensor block cannot send data of type " << type << endl;
        abort();
    }

    if (!metadata_block_) //this was never even initialized ... nothing to send
    {
        if (action == Message::AccumulateUnexpected) //nothing to do
        {
            //the node is not expecting anything... send nothing
            return 0;
        }
        else if (action == Message::GlobalSum)
        {
            //the node expects something... so explicitly send it the zero
            messenger->lock();
            uli initial_tag = 0;
            uli nmessages = 0;
            uli block_number = get_malloc_number();
            size_t data_size = 0;
            SendStatus* status = messenger->send_data_header( //do not track status on global sums
                proc_number,
                header_tag,
                Message::TensorBranch,
                action,
                block_number,
                data_size,
                initial_tag,
                nmessages
            );
            messenger->unlock();
            return status;
        }
    }

#if YETI_SANITY_CHECK
    if (in_read_mode() && fast_read_) 
    {
        //no need to wait on send status
    }
    else if (has_send_status())
    {
        cerr << stream_printf("Block %s already has send status in send_data\n",
                              get_block_name().c_str());
        cerr.flush();
        throw TENSOR_BLOCK_COMMUNICATION_EXCEPTION;
    }

    if (!Sendable::is_waiting() && !is_locked() && !metadata_block_->is_retrieved())
    {
        cerr << stream_printf("Block %s is not retrieved in send\n",
                                get_block_name().c_str());
        cerr.flush();
        abort();
    }
#endif

    StorageBlock* storage_block = get_metadata_storage_block();
    TensorBranch* current_branch = get_branch();
    size_t data_size = current_branch->get_metadata_mempool()->data_size();

    uli block_number = get_malloc_number();
    uli nmessages = current_branch->ncontrollers_nonzero() + 1; //+1 for the metadata
    uli initial_tag = 0;
    if (nmessages == 1) //nothing to send... this has no data
    {
        nmessages = 0;
        if (action == Message::AccumulateUnexpected) //nothing to do
        {
            return 0;
        }
        //even if the branch is zero you still need to send the data header
        //for a global sum because the remote node will expect some
        //signal for the global sum on that block
    }
    else
    {
        initial_tag = messenger->allocate_initial_tag(nmessages, proc_number);
    }

    messenger->lock();

    SendStatus* next_status = messenger->send_data_header(
        proc_number,
        header_tag,
        Message::TensorBranch,
        action,
        block_number,
        data_size,
        initial_tag,
        nmessages
    );

    if (initial_tag == 0) //no data to send
    {
        messenger->unlock();
        return next_status; //nothing to send
    }

    uli tag = initial_tag;
    next_status = messenger->send_no_lock(
        proc_number,
        tag,
        storage_block->data(),
        data_size
    );

    TensorDataController* controller = current_branch->first_data_controller();
    while (controller)
    {
        size_t data_size = controller->get_data_size();
        if (data_size == 0) //nothing to send
        {
            controller = controller->next;
            continue;
        }

#if YETI_SANITY_CHECK
        if (data_size > parent_tensor_->get_data_cache()->blocksize())
        {
            cerr << stream_printf("Send got invalid block size %ld "
                                 "for block %s.  Max block size is %ld.\n",
                                data_size, 
                                get_block_name().c_str(),
                                parent_tensor_->get_data_cache()->blocksize());
            controller_fail();
            abort();
        }
#endif

        ++tag;
        char* data = controller->get_data();
        next_status = messenger->send_no_lock(
            proc_number,
            tag,
            data,
            data_size
        );
        controller = controller->next;
    }
    messenger->unlock();

#if YETI_SANITY_CHECK
    uli ncheck = tag - initial_tag + 1;
    if (ncheck != nmessages)
    {
        std::string str = stream_printf("Wrong number of messages (%ld sent, should be %ld) on"
                                 " tensor block %s",
                                 ncheck,
                                 nmessages,
                                 get_block_name().c_str());
        block_fail(SanityCheckError, str);
    }
#endif
    YetiRuntime::stop_timer("send data");

    return next_status;
}

TensorDataController::TensorDataController(
    TensorBranch* branch,
    StorageBlock* storage_block
) :
    total_size_(0),
    remaining_(0),
    branch_(branch),
    first_data_node_(0),
    last_data_node_(0),
    next(0)
{
    MemoryPool* mempool = branch->get_metadata_mempool();

    total_size_ = storage_block->size();
    remaining_ = total_size_;

    data_ = storage_block->data();
}

TensorDataController::~TensorDataController()
{
}

DataNode*
TensorDataController::first() const
{
    return first_data_node_;
}

bool
TensorDataController::register_node(
    DataNode* node,
    size_t blocksize
)
{
    if (blocksize > remaining_)
        return false;

    size_t current_pos = total_size_ - remaining_;
    node->set_offset(current_pos);
    node->set_data_block(data_);

    node->next_node = 0;

    if (first_data_node_ == 0)
    {
        first_data_node_ = node;
        last_data_node_ = node;
    }
    else
    {
        last_data_node_->next_node = node;
        last_data_node_ = node;
    }

    ::memset(node->data(), 0, blocksize);

    remaining_ -= blocksize;
    return true; //this was valid
}

void
TensorDataController::realign_memory_pool(
    TensorBranch* old_branch,
    TensorBranch* new_branch
)
{
    branch_ = new_branch;
    if (first_data_node_)
    {
        realign_pointer(first_data_node_, old_branch, new_branch);
        realign_pointer(last_data_node_, old_branch, new_branch);

        DataNode* node = first_data_node_;
        while (node->next_node)
        {
            realign_pointer(node->next_node, old_branch, new_branch);
            node = node->next_node;
        }
    }

    if (next)
        realign_pointer(next, old_branch, new_branch);
}

char*
TensorDataController::get_data() const
{
    return data_;
}

size_t
TensorDataController::get_total_size() const
{
    return total_size_;
}

size_t
TensorDataController::get_remaining() const
{
    return remaining_;
}

size_t
TensorDataController::get_data_size() const
{
    return total_size_ - remaining_;
}

void
TensorDataController::set_data(char *data)
{
    data_ = data;
    DataNode* node = first_data_node_;
    while (node)
    {
        node->set_data_block(data_);
        node = node->next_node;
    }
}

TensorController::TensorController()
{
}

void
TensorBlockAccumulateTask::run(uli threadnum)
{
    Sort* no_sort = 0;
    parent_block->retrieve_accumulate();
    parent_block->accumulate(tmp_block, 1.0, no_sort);
    parent_block->release_accumulate();
    tmp_block->release_read();
    tmp_block->lock();
    tmp_block->TensorBlock::~TensorBlock();
    tensor_block_tmp_malloc.free(tmp_block);
}

void
TensorBlockAccumulateTask::prefetch(uli threadnum)
{
}

void
TensorBlockAccumulateTask::finalize_task_subset(uli threadnum)
{
}

uli
TensorBlockAccumulateTask::append_info(uli* data) const
{
    yeti_throw(SanityCheckError, "Tensor block accumulate task cannot append info for dynamic load balance");
    return 0;
}

void*
TensorBlockAccumulateTask::operator new(size_t size)
{
    TensorBlockAccumulateTask* task = &accumulate_tasks[accumulate_task_num];
    ++accumulate_task_num;
    if (accumulate_task_num == NBLOCKS_TMP_ACCUMULATE)
        accumulate_task_num = 0;
    return task;
}

void
TensorBlockAccumulateTask::operator delete(void* ptr)
{
}

void
TensorBlockAccumulateTask::print(std::ostream& os) const
{
    os << "TensorBlockAccumulateTask";
}

