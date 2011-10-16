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

typedef TensorControllerTemplate<
        ValidBranchController,
        DoNothingBranchRenew,
        DoNothingMempool,
        DoNothingBranchRetrieve,
        DoNothingDataControllerRetrieve,
        DoNothingDataControllerInit,
        DoNothingBranchFlush,
        DoNothingBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        DoNothingDataStorageFlush,
        AbortOnObsolete,
        DoNothingSync,
        DeleteAllDataClear,
        DoNothingOutOfCorePrefetch,
        DoNothingInCorePrefetch
    > InCoreController;


typedef TensorControllerTemplate<
                    ValidBranchController,
                    DoNothingBranchRenew,
                    ResetMempool,
                    SortedBranchRetrieve,
                    ReallocateDataControllers,
                    SortDataControllers,
                    ClearBranchFlush,
                    CacheBranchRelease,
                    DoNothingMetaDataRetrieve,
                    DoNothingDataRetrieve,
                    DeleteAllDataStorageFlush,
                    ClearMetaDataOnObsolete,
                    DoNothingSync,
                    ObsoleteAllDataClear,
                    DoNothingOutOfCorePrefetch,
                    DoNothingInCorePrefetch
                > InCoreSortedReadController;

typedef TensorControllerTemplate<
                        ValidBranchController,
                        DoNothingBranchRenew,
                        ResetMempool,
                        NewBranchRetrieve,
                        DoNothingDataControllerRetrieve,
                        DoNothingDataControllerInit,
                        SortedAccumulateBranchFlush,
                        CacheBranchRelease,
                        DoNothingMetaDataRetrieve,
                        DoNothingDataRetrieve,
                        DeleteAllDataStorageFlush,
                        ClearMetaDataOnObsolete,
                        FlushOnSync,
                        ObsoleteAllDataClear,
                        DoNothingOutOfCorePrefetch,
                        DoNothingInCorePrefetch
                   > InCoreSortedAccumulateController;

typedef TensorControllerTemplate<
                    AbortReadBranchValidation,
                    DoNothingBranchRenew,
                    DoNothingMempool,
                    DoNothingBranchRetrieve,
                    DoNothingDataControllerRetrieve,
                    DoNothingDataControllerInit,
                    DoNothingBranchFlush,
                    DoNothingBranchRelease,
                    DoNothingMetaDataRetrieve,
                    DoNothingDataRetrieve,
                    DoNothingDataStorageFlush,
                    AbortOnObsolete,
                    DoNothingSync,
                    DeleteAllDataClear,
                    DoNothingOutOfCorePrefetch,
                    DoNothingInCorePrefetch
                > AbortReadController;

typedef TensorControllerTemplate<
                    AbortWriteBranchValidation,
                    DoNothingBranchRenew,
                    DoNothingMempool,
                    DoNothingBranchRetrieve,
                    DoNothingDataControllerRetrieve,
                    DoNothingDataControllerInit,
                    DoNothingBranchFlush,
                    DoNothingBranchRelease,
                    DoNothingMetaDataRetrieve,
                    DoNothingDataRetrieve,
                    DoNothingDataStorageFlush,
                    AbortOnObsolete,
                    DoNothingSync,
                    DeleteAllDataClear,
                    DoNothingOutOfCorePrefetch,
                    DoNothingInCorePrefetch
                > AbortWriteController;

typedef TensorControllerTemplate<
                    AbortAccumulateBranchValidation,
                    DoNothingBranchRenew,
                    DoNothingMempool,
                    DoNothingBranchRetrieve,
                    DoNothingDataControllerRetrieve,
                    DoNothingDataControllerInit,
                    DoNothingBranchFlush,
                    DoNothingBranchRelease,
                    DoNothingMetaDataRetrieve,
                    DoNothingDataRetrieve,
                    DoNothingDataStorageFlush,
                    AbortOnObsolete,
                    DoNothingSync,
                    DeleteAllDataClear,
                    DoNothingOutOfCorePrefetch,
                    DoNothingInCorePrefetch
                > AbortAccumulateController;

typedef TensorControllerTemplate<
                    AbortVerbatimBranchValidation,
                    DoNothingBranchRenew,
                    DoNothingMempool,
                    DoNothingBranchRetrieve,
                    DoNothingDataControllerRetrieve,
                    DoNothingDataControllerInit,
                    DoNothingBranchFlush,
                    DoNothingBranchRelease,
                    DoNothingMetaDataRetrieve,
                    DoNothingDataRetrieve,
                    DoNothingDataStorageFlush,
                    AbortOnObsolete,
                    DoNothingSync,
                    DeleteAllDataClear,
                    DoNothingOutOfCorePrefetch,
                    DoNothingInCorePrefetch
                > AbortVerbatimController;

typedef TensorControllerTemplate<
                    ValidBranchController,
                    DoNothingBranchRenew,
                    ResetMempool,
                    ActionBranchRetrieve,
                    DoNothingDataControllerRetrieve,
                    DoNothingDataControllerInit,
                    ClearBranchFlush,
                    CacheBranchRelease,
                    DoNothingMetaDataRetrieve,
                    DoNothingDataRetrieve,
                    DeleteAllDataStorageFlush,
                    AbortOnObsolete,
                    DoNothingSync,
                    ObsoleteAllDataClear,
                    DoNothingOutOfCorePrefetch,
                    DoNothingInCorePrefetch
                > ActionReadController;


typedef TensorControllerTemplate<
                    ValidBranchController,
                    DoNothingBranchRenew,
                    ResetMempool,
                    SortedBranchRetrieve,
                    ReallocateDataControllers,
                    SortDataControllers,
                    ClearBranchFlush,
                    CacheBranchRelease,
                    DoNothingMetaDataRetrieve,
                    DoNothingDataRetrieve,
                    DeleteAllDataStorageFlush,
                    ClearMetaDataOnObsolete,
                    DoNothingSync,
                    ObsoleteAllDataClear,
                    DoNothingOutOfCorePrefetch,
                    DoNothingInCorePrefetch
                > ActionSortedReadController;

typedef TensorControllerTemplate<
        ValidBranchController,
        DoNothingBranchRenew,
        DoNothingMempool, //this is taken care of in the prefetch
        RemoteBlockBranchRetrieve,
        DoNothingDataControllerRetrieve,
        DoNothingDataControllerInit,
        ClearBranchFlush,
        CacheBranchRelease,
        DoNothingMetaDataRetrieve,
        DoNothingDataRetrieve,
        DeleteAllDataStorageFlush,
        AbortOnObsolete,
        DoNothingSync,
        DeleteAllDataClear,
        RemoteBlockPrefetch,
        DoNothingInCorePrefetch
    > RemoteReadController;

typedef TensorControllerTemplate<
                    ValidBranchController,
                    RemoteAccumulateRenew,
                    //ZeroBranchRenew,
                    //FlushOldBranchRenew,
                    //DoNothingBranchRenew,
                    ResetMempool,
                    DoNothingBranchRetrieve,
                    DoNothingDataControllerRetrieve,
                    DoNothingDataControllerInit,
                    //RemoteAccumulateFlush,
                    //DoNothingBranchFlush,
                    ClearBranchFlush,
                    RemoteAccumulateBranchRelease,
                    //CacheBranchRelease,
                    DoNothingMetaDataRetrieve,
                    DoNothingDataRetrieve,
                    DeleteAllDataStorageFlush,
                    ClearMetaDataOnObsolete,
                    //FlushOnSync,
                    RemoteAccumulateSync,
                    DeleteAllDataClear,
                    DoNothingOutOfCorePrefetch,
                    DoNothingInCorePrefetch
                > TmpRemoteAccumulateController;

typedef TensorControllerTemplate<
                    ValidBranchController,
                    DoNothingBranchRenew,
                    ResetMempool,
                    DoNothingBranchRetrieve,
                    DoNothingDataControllerRetrieve,
                    DoNothingDataControllerInit,
                    RemoteAccumulateFlush,
                    CacheBranchRelease,
                    DoNothingMetaDataRetrieve,
                    DoNothingDataRetrieve,
                    DeleteAllDataStorageFlush,
                    ClearMetaDataOnObsolete,
                    RemoteAccumulateSync,
                    DeleteAllDataClear,
                    DoNothingOutOfCorePrefetch,
                    DoNothingInCorePrefetch
                > RemoteAccumulateController;

typedef TensorControllerTemplate<
                    ValidBranchController,
                    ZeroBranchRenew,
                    ResetMempool,
                    DoNothingBranchRetrieve,
                    DoNothingDataControllerRetrieve,
                    DoNothingDataControllerInit,
                    ClearBranchFlush,
                    ThreadAccumulateBranchRelease,
                    DoNothingMetaDataRetrieve,
                    DoNothingDataRetrieve,
                    DeleteAllDataStorageFlush,
                    ClearMetaDataOnObsolete,
                    AbortOnSync,
                    DeleteAllDataClear,
                    DoNothingOutOfCorePrefetch,
                    DoNothingInCorePrefetch
                > TmpThreadAccumulateController;

typedef TensorControllerTemplate<
                    ValidBranchController,
                    DoNothingBranchRenew,
                    ResetMempool,
                    SortedBranchRetrieve,
                    ReallocateDataControllers,
                    SortDataControllers,
                    ClearBranchFlush,
                    CacheBranchRelease,
                    DoNothingMetaDataRetrieve,
                    DoNothingDataRetrieve,
                    DeleteAllDataStorageFlush,
                    ClearMetaDataOnObsolete,
                    DoNothingSync,
                    ObsoleteAllDataClear,
                    ParentBlockReadPrefetch,
                    DoNothingInCorePrefetch
                > RemoteSortedReadController;

typedef TensorControllerTemplate<
                    ValidBranchController,
                    DoNothingBranchRenew,
                    ResetMempool,
                    ConfigureElementComputerBranchRetrieve,
                    DoNothingDataControllerRetrieve,
                    DoNothingDataControllerInit,
                    ClearBranchFlush,
                    CacheBranchRelease,
                    RecomputeMetaDataRetrieve,
                    RecomputeDataRetrieve,
                    DeleteAllDataStorageFlush,
                    ClearMetaDataOnObsolete,
                    DoNothingSync,
                    ObsoleteAllDataClear,
                    DoNothingOutOfCorePrefetch,
                    ResortInCorePrefetch
                > RecomputeReadController;

typedef TensorControllerTemplate<
                    ValidBranchController,
                    DoNothingBranchRenew,
                    DoNothingMempool,
                    RealignMemoryPoolBranchRetrieve,
                    ReuseDataControllers,
                    DoNothingDataControllerInit,
                    ClearBranchFlush,
                    CacheBranchRelease,
                    DoNothingMetaDataRetrieve,
                    DoNothingDataRetrieve,
                    ClearAllDataStorageFlush,
                    AbortOnObsolete,
                    DoNothingSync,
                    ObsoleteAllDataClear,
                    DoNothingOutOfCorePrefetch,
                    DoNothingInCorePrefetch
                > DiskReadController;

typedef TensorControllerTemplate<
                    ValidBranchController,
                    DoNothingBranchRenew,
                    DoNothingMempool,
                    RealignMemoryPoolBranchRetrieve,
                    ReuseDataControllers,
                    DoNothingDataControllerInit,
                    CommitBranchFlush,
                    CacheBranchRelease,
                    DoNothingMetaDataRetrieve,
                    DoNothingDataRetrieve,
                    CommitAllDataStorageFlush,
                    AbortOnObsolete,
                    DoNothingSync,
                    ObsoleteAllDataClear,
                    DoNothingOutOfCorePrefetch,
                    DoNothingInCorePrefetch
                > DiskWriteController;
                

typedef TensorControllerTemplate<
                    ValidBranchController,
                    DoNothingBranchRenew,
                    DoNothingMempool,
                    RealignMemoryPoolBranchRetrieve,
                    ReuseDataControllers,
                    DoNothingDataControllerInit,
                    CommitBranchFlush,
                    CacheBranchRelease,
                    DoNothingMetaDataRetrieve,
                    DoNothingDataRetrieve,
                    CommitAllDataStorageFlush,
                    AbortOnObsolete,
                    DoNothingSync,
                    ObsoleteAllDataClear,
                    DoNothingOutOfCorePrefetch,
                    DoNothingInCorePrefetch
                > DiskAccumulateController;

typedef TensorControllerTemplate<
                    ValidBranchController,
                    DoNothingBranchRenew,
                    ResetMempool,
                    SortedBranchRetrieve,
                    ReallocateDataControllers,
                    SortDataControllers,
                    ClearBranchFlush,
                    CacheBranchRelease,
                    DoNothingMetaDataRetrieve,
                    DoNothingDataRetrieve,
                    DeleteAllDataStorageFlush,
                    ClearMetaDataOnObsolete,
                    DoNothingSync,
                    ObsoleteAllDataClear,
                    DoNothingOutOfCorePrefetch,
                    DoNothingInCorePrefetch
                > DiskSortedReadController;

    

DECLARE_PARENT_MALLOC(TensorController);
DECLARE_SUB_MALLOC(TensorController,AbortReadController);

DataStorageBlockAllocator::DataStorageBlockAllocator()
    :
  nstorage_blocks_(0),
  first_node_(0),
  last_node_(0)
{
}

DataStorageBlockAllocator::DataStorageBlockAllocator(YetiRuntimeObject::thread_safety_flag_t flag)
    :
  Cachable(flag),
  nstorage_blocks_(0),
  first_node_(0),
  last_node_(0)
{
}

DataStorageBlockAllocator::~DataStorageBlockAllocator()
{
    clear_storage();
}

uli
DataStorageBlockAllocator::get_num_data_storage_blocks() const
{
    return nstorage_blocks_;
}

DataStorageNode*
DataStorageBlockAllocator::first_data_node() const
{
    return first_node_;
}

void
DataStorageBlockAllocator::store(StorageBlock* block)
{
    DataStorageNode* node = new DataStorageNode;
    node->block = block;
    node->next = 0;

    if (first_node_ == 0)
    {
        first_node_ = node;
        last_node_ = node;
    }
    else
    {
        last_node_->next = node;
        last_node_ = node;
    }
    ++nstorage_blocks_;
}

void
DataStorageBlockAllocator::flush_from_cache()
{
    //do nothing
}

CachedStorageBlock*
DataStorageBlockAllocator::allocate_cache_block(DataCache* cache)
{
    CachedStorageBlock* block = new CachedStorageBlock(this, cache);
    block->retrieve();

    return block;
}

void
DataStorageBlockAllocator::clear_storage()
{
    DataStorageNode* node = first_node_;
    while (node)
    {
        StorageBlock* block = node->block;
        delete block;
        DataStorageNode* next = node->next;
        DataStorageNode* prev = node;
        node->block = 0;
        node->next = 0;
        node = next;
        delete prev;
    }
    nstorage_blocks_ = 0;
    first_node_ = 0;
    last_node_ = 0;
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
    read_controller_(0),
    write_controller_(0),
    accumulate_controller_(0),
    verbatim_controller_(0),
    resort_perm_(0),
    fetch_perm_(0),
    branch_(0),
    block_number_(parent->get_index(indexset)),
    depth_(parent->get_block_descr()->depth() - 1),
    nindex_(parent->get_block_descr()->nindex()),
    metadata_mempool_(0),
    data_block_size_(parent->data_block_size()),
    metadata_block_size_(parent->metadata_block_size()),
    total_nelements_data_(0),
    parent_tensor_(parent),
    descr_(parent->get_block_descr()),
    degeneracy_(1),
    is_synced_(true), //no data yet, so true
    permutationally_unique_(false), //begin with false
    is_subblock_(false),
    block_perm_(0),
    action_(0),
    node_number_(YetiRuntime::me()),
    malloc_number_(0),
    up_to_date_(true),
    disk_buffer_(0),
    storage_type_(Tensor::in_core),
    metadata_block_(0),
    tmp_block_(false)
{
    malloc_number_ = Malloc<TensorBlock>::get_malloc_number(this);

    ::memcpy(indices_, indexset, descr_->nindex() * sizeof(uli));

    if (parent_tensor_->is_distributed())
        node_number_ = parent_tensor_->get_node_number(block_number_);

    recompute_permutation();

    block_perm_ = parent_tensor_->get_tensor_grp()->get_identity();

    total_nelements_data_ = descr_->get_nelements_data(depth_ + 1, indices_);
    set_element_size(sizeof(double));
}

TensorBlock::TensorBlock(
    Tensor* parent
) :
    controller_(0),
    read_controller_(0),
    write_controller_(0),
    accumulate_controller_(0),
    verbatim_controller_(0),
    resort_perm_(0),
    fetch_perm_(0),
    branch_(0),
    block_number_(0),
    depth_(parent->get_block_descr()->depth() - 1),
    nindex_(parent->get_block_descr()->nindex()),
    metadata_mempool_(0),
    data_block_size_(parent->data_block_size()),
    metadata_block_size_(parent->metadata_block_size()),
    total_nelements_data_(5e5), //set to a ridiculous size
    parent_tensor_(parent),
    descr_(parent->get_block_descr()),
    degeneracy_(1),
    is_synced_(true), //no data yet, so true
    permutationally_unique_(true),
    is_subblock_(false),
    block_perm_(0),
    action_(0),
    node_number_(YetiRuntime::me()),
    malloc_number_(0),
    up_to_date_(true),
    disk_buffer_(0),
    storage_type_(Tensor::in_core),
    metadata_block_(0),
    tmp_block_(true)
{
    malloc_number_ = Malloc<TensorBlock>::get_malloc_number(this);
    for (uli i=0; i < descr_->nindex(); ++i)
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
    write_controller_(0),
    accumulate_controller_(0),
    verbatim_controller_(0),
    resort_perm_(parent->get_resort_permutation()),
    fetch_perm_(parent->get_fetch_permutation()),
    branch_(0),
    block_number_(0),
    depth_(parent->get_descr()->depth() - 1),
    nindex_(parent->get_descr()->nindex()),
    metadata_mempool_(0),
    parent_tensor_(parent->get_parent_tensor()),
    data_block_size_(parent->data_block_size_),
    metadata_block_size_(parent->metadata_block_size_),
    descr_(parent->get_descr()),
    degeneracy_(parent->get_degeneracy()),
    is_synced_(true), //no data yet, so true
    permutationally_unique_(true),
    is_subblock_(parent->is_subblock()),
    block_perm_(parent->get_block_permutation()),
    node_number_(YetiRuntime::me()),
    malloc_number_(0),
    up_to_date_(true),
    storage_type_(Tensor::in_core),
    tmp_block_(true)
{
    malloc_number_ = parent->get_malloc_number();
    ::memcpy(indices_, parent->get_indices(), descr_->nindex() * sizeof(uli));
    task_owner_number_ = parent->get_task_owner_number();

    //this is only ever used for temporary accumulation
    init_in_core_no_sort();
    controller_ = accumulate_controller_;

}

TensorBlock::~TensorBlock()
{
#if YETI_SANITY_CHECK
    if (!is_locked())
        raise(SanityCheckError, "cannot delete an unlocked tensor block");
#endif
    /** The block might never have been initialized */
    if (read_controller_)
    {
        in_destructor_ = true;

        if (is_cached())
        {
            //always lock an object before flushing it
            //once the tensor is deleted, it goes out of cache
            flush_from_cache();
        }

        //DeleteAllDataClear clear;
        //clear.clear(this, metadata_block_, storage_blocks_, nstorage_blocks_);
        read_controller_->clear(this);

        clear_storage();
    }

}

void
TensorBlock::recompute_permutation()
{
    recompute_tensor_permutation();
    permutationally_unique_ = fetch_perm_->is_identity();
    task_owner_number_ = parent_tensor_->get_unique_id(indices_, fetch_perm_);
}

void
TensorBlock::init()
{
    if (controller_) /** already initialized */
        return;

#if YETI_SANITY_CHECK
    if (!is_locked())
        raise(SanityCheckError, "block must be locked before initializing");

    if (metadata_block_)
        raise(SanityCheckError, "block is already initialized");
#endif


    nstorage_blocks_ = 0;

    if (is_remote_block())
    {
        init_remote();
    }
    else
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

    controller_ = read_controller_;
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
    up_to_date_ = false;

    branch_->set_element_type(src->get_branch()->element_type());
    if (src == this) //self-accumulation
    {
        if (!sort)
            raise(SanityCheckError, "you should never self-accumulate unless sorting!");

        CachedStorageBlock* block
            = new CachedStorageBlock(this, parent_tensor_->get_metadata_cache());
        block->retrieve();
        MemoryPool* mempool = new MemoryPool(block->size(), block->data());
        mempool->memcpy(metadata_mempool_.get());
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
        DataStorageNode* node = allocator.first_data_node();
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
        raise(SanityCheckError, "block has zero degeneracy");
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
TensorBlock::allocate_data_storage_block()
{
    return allocate_storage_block(data_block_size_, TensorBlock::data_block);
}

StorageBlock*
TensorBlock::allocate_storage_block(
    size_t size,
    TensorBlock::storage_block_type_t blocktype
)
{
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
            cerr << "cache blocks are not large enough!"
                 << endl
                 << " cache block size is " << cache->blocksize()
                 << " but requested block size is " << size
                 << endl;
            abort();
        }
        block = allocate_cache_block(cache);
    }
    else if (storage_type == Tensor::in_core)
    {
        block = this->allocate_core_block(size);
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
        block = this->allocate_disk_block(
                    disk_buffer_.get(),
                    this,
                    cache
                   );
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
    return allocate_storage_block(metadata_block_size_, TensorBlock::metadata_block);

}

void
TensorBlock::set_element_size(size_t size)
{
    size_t max_datasize = total_nelements_data_ * size;
#if 0
    if (max_datasize < data_block_size_)
    {
        data_block_size_ = max_datasize;
    }
#endif
    if (permutationally_unique_ && parent_tensor_->get_storage_type() == Tensor::in_core)
    {
        if (max_datasize < data_block_size_)
            data_block_size_ = max_datasize;
    }
}

void
TensorBlock::fill(TensorElementComputer* filler)
{
    if (!resort_perm_->is_identity()) //this is not a unique tile
        return;

    TemplateInfo::type_t element_type
        = filler->element_type(indices_, depth_);

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
TensorBlock::get_block_number() const
{
    return block_number_;
}

usi
TensorBlock::get_degeneracy() const
{
    return degeneracy_;
}

usi
TensorBlock::get_depth() const
{
    return depth_;
}

const uli*
TensorBlock::get_indices() const
{
    return indices_;
}

std::string
TensorBlock::get_index_string() const
{
    return ClassOutput<const uli*>::str(nindex_, indices_);
}

std::string
TensorBlock::get_block_name() const
{
    std::string name = get_index_string();
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
    return metadata_mempool_.get();
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
TensorBlock::configure(const DiskBufferPtr& disk_buffer)
{
    configure(Tensor::on_disk);
    disk_buffer_ = disk_buffer;
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

    metadata_block_ = 0;
}

StorageBlock*
TensorBlock::get_metadata_storage_block() const
{
    return metadata_block_.get();
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

DiskBuffer*
TensorBlock::get_disk_buffer() const
{
    return disk_buffer_.get();
}

TensorBlock*
TensorBlock::get_symmetry_unique_block() const
{
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
    branch_ = new (metadata_mempool_.get())
        TensorBranch(
            indices_,
            depth_,
            this
        );
    branch_->set_descr(descr_);
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
    char* metadata_data = YetiRuntime::malloc(metadata_block_size_);
    metadata_block_ = new InCoreBlock(metadata_data, metadata_block_size_);

    init_mempool();

    init_branch();

    read_controller_ = abort_read_controller_;

    accumulate_controller_ = tmp_thread_accumulate_controller_;

    write_controller_ = abort_write_controller_;

    verbatim_controller_ = abort_verbatim_controller_;

    controller_ = accumulate_controller_;

    set_initialized(true);
}

void
TensorBlock::init_in_core_no_sort()
{
    char* metadata_data = YetiRuntime::malloc(metadata_block_size_);
    metadata_block_ = new InCoreBlock(metadata_data, metadata_block_size_);
    init_mempool();
    init_branch();

    read_controller_ = in_core_controller_;

    write_controller_ = in_core_controller_;

    accumulate_controller_ = in_core_controller_;

    verbatim_controller_ = in_core_controller_;

    set_initialized(true);
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
    metadata_block_ = new LocalDiskBlock(
                        this,
                        parent_tensor_->get_metadata_cache(),
                        disk_buffer_.get()
                      );
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
    wait_on_send(); //very important...
    if (block->is_remote_block())
    {
        node_number_ = block->node_number_;
        malloc_number_ = block->malloc_number_;
        accumulate_controller_ = tmp_remote_accumulate_controller_;
    }
    else
    {
        accumulate_controller_ = tmp_thread_accumulate_controller_;
    }

    for (uli i=0; i < nindex_; ++i)
        indices_[i] = block->indices_[i];

    block_number_ = block->block_number_;
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

    accumulate_controller_ = in_core_sorted_accumulate_controller_;

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
TensorBlock::is_remote_block() const
{
    bool remote_block = YetiRuntime::me() != node_number_;
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
        raise(SanityCheckError, "cannot obsolete unlocked tensor block");

    //it's possible this is not symmetry unique and was never used
    if (controller_)
    {
        controller_->obsolete(this);
    }
}

void
TensorBlock::print(std::ostream& os)
{
    NormElementOp* op = new NormElementOp;
    element_op(op);
    double norm = op->norm();
    double maxlog = branch_->get_node()->get_max_log();

    set_synced();

    os    << "Degeneracy: " << degeneracy_ << " NBlocks: " << nstorage_blocks_ << endl
        << stream_printf("max log: %8.4f    norm: %8.4e", maxlog, norm)
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
TensorBlock::retrieve_read()
{
    lock();
    set_read_mode();
    YetiRuntimeObject::retrieve();
    if (storage_type_ != Tensor::recomputed)
        unlock();
}

void
TensorBlock::release_read()
{
    if (storage_type_ == Tensor::recomputed)
    {
        YetiRuntimeObject::release();
        unlock();
    }
    else
    {
        lock();
        YetiRuntimeObject::release();
        unlock();
    }

}

void
TensorBlock::retrieve_accumulate_no_lock()
{
    set_accumulate_mode();
    YetiRuntimeObject::retrieve();
}

void
TensorBlock::retrieve_accumulate()
{
    lock();
    retrieve_accumulate_no_lock();
}

void
TensorBlock::release_accumulate()
{
    YetiRuntimeObject::release();
    unlock();
}

void
TensorBlock::retrieve_verbatim()
{
    lock();
    set_verbatim_mode();
    YetiRuntimeObject::retrieve();
}

void
TensorBlock::release_verbatim()
{
    YetiRuntimeObject::release();
    unlock();
}

void
TensorBlock::retrieve_write()
{
    lock();
    set_write_mode();
    YetiRuntimeObject::retrieve();
}

void
TensorBlock::release_write()
{
    YetiRuntimeObject::release();
    unlock();
}

void
TensorBlock::finalize()
{
    YetiRuntimeObject::finalize();
}

void
TensorBlock::_initialize()
{
    if (!controller_)
        raise(SanityCheckError, "tensor block not yet initialized");

#if YETI_SANITY_CHECK
    //cache block should not already be retrieved
    if (metadata_block_->is_cached() && metadata_block_->is_retrieved())
        raise(SanityCheckError, "metadata block is retrieved already in initialize");
#endif

    controller_->validate(this);
    bool in_core = metadata_block_->retrieve();

    branch_ = reinterpret_cast<TensorBranch*>(metadata_block_->data());
    metadata_mempool_->set(metadata_block_->data());
    branch_->set_metadata_mempool(metadata_mempool_.get());
    branch_->set_descr(descr_);
    branch_->set_parent(this);

    DataStorageNode* node = first_data_node(); 
    while(node)
    {
        StorageBlock* block = node->block;
        block->retrieve();
        node = node->next;
    }
}

void
TensorBlock::_renew()
{
    wait_on_send();
    reset_send_status();

    controller_->renew(this);
    set_synced();
}

void
TensorBlock::_retrieve()
{
    wait_on_send();

    //reset the send status after the retrieve
    controller_->retrieve(this);

    reset_send_status();

    //after retrieving, we are synced with the parent storage location
    set_synced();
}

void
TensorBlock::_release()
{
    controller_->release(this);
}
/*************
End private interface definition
*************/

void
TensorBlock::release_callback()
{
    bool locked = trylock();
    if (!locked) //another thread came along to renew this
        return; //just leave everything alone

    metadata_block_->release();
    DataStorageNode* node = first_data_node();
    while(node)
    {
        StorageBlock* data_block = node->block;
        data_block->release();
        node = node->next;
    }

    reset_send_status();
    finalize();

    unlock();
}

void
TensorBlock::preflush()
{
#if YETI_SANITY_CHECK
    if (!is_locked())
        raise(SanityCheckError, "cannot flush unlocked tensor block");
#endif
    controller_->preflush(this);
}

/**********
Begin Cachable interface definition
**********/
void
TensorBlock::flush_from_cache()
{
    if (!controller_) //this was never actually used for anything
        return;

#if YETI_SANITY_CHECK
    if (metadata_block_->is_retrieved())
        raise(SanityCheckError, "cannot flush retrieved block");
#endif

    wait_on_send();

#if YETI_SANITY_CHECK
    if (!is_locked())
        raise(SanityCheckError, "cannot flush unlocked tensor block");
#endif

    controller_->flush(this);
    controller_->flush_data(this);

    reset_send_status();

#if YETI_SANITY_CHECK
    CachedStorageBlock* block = dynamic_cast<CachedStorageBlock*>(metadata_block_.get());
    if (block)
    {
        if (block->get_entry())
        {
            raise(SanityCheckError, "cache block not cleared in flush");
        }
        if (block->data())
        {
            raise(SanityCheckError, "block not cleared in flush");
        }
    }
#endif

    finalize();
}

/********
End Cachable interface definition
*********/

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
        return false;
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

    if (!is_synced())
    {
        cerr << "Cannot switch to verbatim mode.  Tensor block has"
             << " not been synced from previous operation" << endl;
        cerr << ClassOutput<const uli*>
                ::str(descr_->nindex(), indices_) << endl;
        abort();
    }

    controller_ = verbatim_controller_;
}

void
TensorBlock::set_read_mode()
{

    if (!controller_)
    {
        if (!controller_) //another thread may have come in and done this
            init();
    }

    //already in correct mode
    if (controller_ == read_controller_)
        return;
    
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
    
    if (!is_synced())
    {
        cerr << "Cannot switch to write mode.  Tensor block has"
             << " not been synced from previous operation" << endl;
        cerr << ClassOutput<const uli*>
                ::str(descr_->nindex(), indices_) << endl;
        abort();
    }
    controller_ = write_controller_;
}

void
TensorBlock::metadata_sort(Sort* sort)
{
    uli* tmp = reinterpret_cast<uli*>(sort->metadata_buffer(0));
    sort->get_permutation()->permute<uli>(indices_, tmp);
    ::memcpy(indices_, tmp, nindex_ * sizeof(uli));

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
    if (is_synced())
        return;

    controller_->sync(this);
    set_synced();
}

void
TensorBlock::update()
{
    reset_degeneracy();

    if (is_permutationally_unique())
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
TensorBlock::prefetch_read()
{
    lock();
    set_read_mode();
    controller_->out_of_core_prefetch(this);
    unlock();
}

void
TensorBlock::prefetch_accumulate()
{
    lock();
    set_accumulate_mode();
    controller_->out_of_core_prefetch(this);
    unlock();
}

void
TensorBlock::prefetch_read(TensorBlock* prev_block)
{
    lock();
    set_read_mode();
    controller_->in_core_prefetch(this, prev_block);
    unlock();
}

void
TensorBlock::prefetch_accumulate(TensorBlock* prev_block)
{
    lock();
    set_accumulate_mode();
    controller_->in_core_prefetch(this, prev_block);
    unlock();
}


void
TensorBlock::recv_data(
    YetiMessenger* messenger,
    Message::message_data_type_t type,
    uli proc_number,
    size_t metadata_size,
    uli initial_tag,
    uli nmessages
)
{

#if YETI_SANITY_CHECK
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

        if (ncheck > 1e7)
            raise(SanityCheckError, "block is retrieved and cannot receive data after 10 seconds");
    }
#endif

    if (type != Message::TensorBranch)
        raise(SanityCheckError, "tensor block can only receive branch data type");


    receive_branch(messenger, proc_number, metadata_size, initial_tag, nmessages);
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
    heisenfxn(TensorBlock::accumulate_data);
    if (type != Message::TensorBranch)
        raise(SanityCheckError, "tensor block can only receive branch data type");

    MallocOverride m;

    //allocate a temporary block
    TensorBlock tmp_block(this);

    tmp_block.retrieve_read();
    tmp_block.receive_branch(messenger, proc_number, metadata_size, initial_tag, nmessages);

    Sort* sort = 0;
    retrieve_accumulate();
    accumulate(&tmp_block, 1.0, 0);
    release_accumulate();

    tmp_block.release_read();
    tmp_block.lock();
    //currently allocated in core
    //tmp_block.flush_from_cache();
    heisenfxn(TensorBlock::accumulate_data);
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
    heisenfxn(TensorBlock::receive_branch);
#if YETI_SANITY_CHECK
    if (!is_initialized())
        raise(SanityCheckError, "block is not initialized");
#endif
    if (metadata_size > metadata_block_->size())
    {
        Env::errn() << "Metadata block is not allocated large enough for recv" << endl;
        abort();
    }

    messenger->lock();

    heisenfxn(TensorBlock::receive_branch);
    uli tag = initial_tag;
    messenger->recv_no_lock(
        proc_number,
        tag, 
        branch_,
        metadata_size
    );
    //parent block gets overwritten

    heisenfxn(TensorBlock::receive_branch);
    TensorBranch* obsolete_branch = branch_->get_last_branch_location();
    if (branch_ != obsolete_branch)
        branch_->realign_memory_pool(obsolete_branch, branch_);

    branch_->set_parent(this);
    branch_->set_descr(descr_);
    branch_->set_metadata_mempool(metadata_mempool_.get());
    metadata_mempool_->set_data_size(metadata_size);

    TensorDataController* controller = branch_->first_data_controller();
    DataStorageNode* dnode = first_node_;
    heisenfxn(TensorBlock::receive_branch);
    while(controller)
    {
        ++tag;
        size_t data_size = controller->get_data_size();
        StorageBlock* data_block = 0;
        if (dnode)
        {
            data_block = dnode->block;
            dnode = dnode->next;
        }
        else
        {
            data_block = allocate_data_storage_block();
        }

        messenger->recv_no_lock(
            proc_number,
            tag,
            data_block->data(),
            data_size
        );
        controller->set_data(data_block->data());
        controller = controller->next;
        heisenfxn(TensorBlock::receive_branch);
    }
    messenger->unlock();


#if YETI_SANITY_CHECK
    uli ncheck = branch_->ncontrollers() + 1;
    if (ncheck != nmessages)
        raise(SanityCheckError, "wrong number of messages sent to tensor block");
#endif


    set_waiting(false);

    heisenfxn(TensorBlock::receive_branch);
}


SendStatus*
TensorBlock::send_data(
    YetiMessenger* messenger,
    Message::message_data_type_t type,
    Message::message_action_t action,
    uli proc_number
)
{
    heisenfxn(TensorBlock::send_data);
    if (type != Message::TensorBranch)
    {
        Env::errn() << "Tensor block cannot send data of type " << type << endl;
        abort();
    }

#if YETI_SANITY_CHECK
    if (!Sendable::is_waiting() && !is_locked() && !metadata_block_->is_retrieved())
    {
        cerr << "block " << get_block_name() << " is not retrieved in send " << endl;
        abort();
    }
#endif

    StorageBlock* storage_block = get_metadata_storage_block();
    TensorBranch* current_branch = get_branch();
    size_t data_size = current_branch->get_metadata_mempool()->data_size();

    uli block_number = get_malloc_number();
    uli nmessages = current_branch->ncontrollers() + 1;
    uli initial_tag = messenger->allocate_initial_tag(nmessages);

    cout << stream_printf("allocate on node %d initial tag %d for %d messages",
			YetiRuntime::me(), initial_tag, nmessages) << endl;

    messenger->lock();

    heisenfxn(TensorBlock::send_data);
    SendDataHeader header;
    SendStatus* next_status = messenger->send_data_header(
        &header,
        proc_number,
        Message::TensorBranch,
        action,
        block_number,
        data_size,
        initial_tag,
        nmessages
    );

    heisenfxn(TensorBlock::send_data);
    uli tag = initial_tag;
    SendStatus* prev_status = next_status;
    next_status = messenger->send_no_lock(
        proc_number,
        tag,
        storage_block->data(),
        data_size
    );
    next_status->set_dependent(prev_status);

    TensorDataController* controller = current_branch->first_data_controller();
    while (controller)
    {
        heisenfxn(TensorBlock::send_data);

        ++tag;

        char* data = controller->get_data();
        size_t data_size = controller->get_data_size();

        next_status = messenger->send_no_lock(
            proc_number,
            tag,
            data,
            data_size
        );
        next_status->set_dependent(prev_status);
        controller = controller->next;
    }
    messenger->unlock();

#if YETI_SANITY_CHECK
    uli ncheck = tag - initial_tag + 1;
    if (ncheck != nmessages)
        raise(SanityCheckError, "wrong number of messages sent to tensor block");
#endif
    
    heisenfxn(TensorBlock::send_data);
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
