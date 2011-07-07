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

#define LEAK_CHECK 0

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

//static std::map<void*,int> blocks;


DECLARE_MALLOC(TensorBlock);
DECLARE_MALLOC(DataStorageNode);

typedef TensorControllerTemplate<
                    AbortReadBranchValidation,
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
                    DeleteAllDataClear  
                > AbortReadTensorController;

typedef TensorControllerTemplate<
                    AbortWriteBranchValidation,
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
                    DeleteAllDataClear
                > AbortWriteTensorController;

typedef TensorControllerTemplate<
                    AbortAccumulateBranchValidation,
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
                    DeleteAllDataClear
                > AbortAccumulateTensorController;

typedef TensorControllerTemplate<
                    AbortVerbatimBranchValidation,
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
                    DeleteAllDataClear
                > AbortVerbatimTensorController;

DECLARE_PARENT_MALLOC(TensorController);
DECLARE_SUB_MALLOC(TensorController,AbortReadTensorController);

typedef TensorControllerTemplate<
                    ValidBranchController,
                    ResetMempool,
                    SortedBranchRetrieve,
                    ReallocateDataControllers,
                    DoNothingDataControllerInit,
                    ClearBranchFlush,
                    DoNothingBranchRelease,
                    DoNothingMetaDataRetrieve,
                    DoNothingDataRetrieve,
                    DeleteAllDataStorageFlush,
                    ClearMetaDataOnObsolete,
                    DoNothingSync,
                    DeleteAllDataClear  
                > SortedCopyController;

typedef   TensorControllerTemplate<
                        ValidBranchController,
                        DoNothingMempool,
                        RealignMemoryPoolBranchRetrieve,
                        ReallocateDataControllers,
                        DoNothingDataControllerInit,
                        ClearBranchFlush,
                        DoNothingBranchRelease,
                        DoNothingMetaDataRetrieve,
                        DoNothingDataRetrieve,
                        ClearAllDataStorageFlush,
                        AbortOnObsolete,
                        DoNothingSync,
                        DeleteAllDataClear  
                    > CopyController;
                    
                    
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

#if DEBUG_CACHE_USAGE
    TensorBlock* tblock = static_cast<TensorBlock*>(this);
    dout << "allocated data cache entry " << block->get_entry()->offset
         << " to block " << tblock->get_index_string()
         << " on tensor " << tblock->get_parent_tensor()->get_name()
         << " at location " << (void*) block
         << endl;
#endif

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
    obsolete_branch_(0),
    block_number_(parent->get_index(indexset)),
    depth_(parent->get_descr()->depth() - 1),
    nindex_(parent->get_descr()->nindex()),
    metadata_mempool_(0),
    data_block_size_(parent->data_block_size()),
    metadata_block_size_(parent->metadata_block_size()),
    total_nelements_data_(0),
    parent_tensor_(parent),
    descr_(parent->get_descr()),
    is_flushed_(false),
    degeneracy_(1),
    maxlog_(LOG_UNITIALIZED),
    is_synced_(true), //no data yet, so true
    permutationally_unique_(parent->is_unique(indexset)),
    is_subblock_(false),
    block_perm_(0),
    action_(0)
{
    ::memcpy(indices_, indexset, descr_->nindex() * sizeof(uli));

#if DEBUG_CACHE_USAGE
    dout << "constructed block " << get_index_string() << " on tensor "
         << parent_tensor_->get_name()
         << " at location " << (void*) this
         << endl;
#endif

    recompute_uniqueness();
    recompute_tensor_permutation();
    reinit(store_type);

    block_perm_ = parent_tensor_->get_sort_permutation();

    task_owner_number_ = parent_tensor_->get_unique_id(indices_, fetch_perm_);

    total_nelements_data_ = descr_->get_nelements_data(depth_ + 1, indices_);
    set_element_size(sizeof(double));

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
    obsolete_branch_(0),
    block_number_(0),
    depth_(parent->get_descr()->depth() - 1),
    nindex_(parent->get_descr()->nindex()),
    metadata_mempool_(0),
    parent_tensor_(parent->get_parent_tensor()),
    data_block_size_(parent->data_block_size_),
    metadata_block_size_(parent->metadata_block_size_),
    descr_(parent->get_descr()),
    is_flushed_(false),
    degeneracy_(parent->get_degeneracy()),
    maxlog_(LOG_UNITIALIZED),
    is_synced_(true), //no data yet, so true
    permutationally_unique_(true),
    is_subblock_(parent->is_subblock()),
    block_perm_(parent->get_block_permutation())
{
    ::memcpy(indices_, parent->get_indices(), descr_->nindex() * sizeof(uli));
    task_owner_number_ = parent->get_task_owner_number();
    init_in_core_no_sort();

}

TensorBlock::~TensorBlock()
{
    if (!read_controller_)
    {
        cerr << "block " << get_index_string() << " on tensor " << parent_tensor_->get_name()
                << " never initialized" << endl;
        abort();
    }
    in_destructor_ = true;

    if (in_cache())
        //once the tensor is deleted, it goes out of cache
        flush_from_cache();

    //DeleteAllDataClear clear;
    //clear.clear(this, metadata_block_, storage_blocks_, nstorage_blocks_);
    read_controller_->clear();

    clear_storage();

    delete read_controller_;
    delete write_controller_;
    delete accumulate_controller_;
    delete verbatim_controller_;
}

void
TensorBlock::reinit(Tensor::tensor_storage_t store_type)
{
    if (metadata_block_)
    {
#if DEBUG_CACHE_USAGE
        dout << "reinit block " << get_index_string()
            << " on tensor " << parent_tensor_->get_name() 
            << " at location " << (void*) this
            << endl;
#endif
        controller_->clear();
        nstorage_blocks_ = 0;
    }
        
    controller_ = 0;
    
    if (read_controller_)
    {
        delete read_controller_;
        read_controller_ = 0;
    }

    if (write_controller_)
    {
        delete write_controller_;
        write_controller_ = 0;
    }

    if (accumulate_controller_)
    {
        delete accumulate_controller_;
        accumulate_controller_ = 0;
    }

    if (is_remote_block())
    {
        init_remote();
    }
    else
    {
        if (store_type == Tensor::default_storage)
            store_type = parent_tensor_->get_storage_type();

        switch(store_type)
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
            //nothing to accumulate
            block->free_cache_entry();
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

        Permutation* descr_perm = 0;
        branch_->get_node()->sort_data_into(
            sort,
            branch_,
            sorted_branch,
            descr_perm
        );
;

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
            block->free_cache_entry();
        }
        delete mempool;
        block->release();
        delete block;
    }
    else
    {
        branch_->get_node()->accumulate(
            src->get_branch()->get_node(),
            scale,
            sort
        );
    }
    
    maxlog_ = branch_->get_node()->get_max_log();
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

    if (!permutationally_unique_ ||
        storage_type == Tensor::recomputed ||
        action_ ||
        is_remote_block())
    {
#if MALLOC_CACHE_BLOCK
        block = this->allocate_core_block(size);
#else
        if (cache->blocksize() < size)
        {
            cerr << "cache blocks are not large enough!"
                 << endl
                 << " cache block size is " << cache->blocksize()
                 << " but requested block size is " << size
                 << endl;
        }
        block = allocate_cache_block(cache);
#endif
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
                    parent_tensor_->get_disk_buffer(),
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
        store(block);

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

    maxlog_ = branch_->get_node()->get_max_log();
}

TensorIndexDescr*
TensorBlock::get_descr() const
{
    return parent_tensor_->get_descr();
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
TensorBlock::configure(const TensorRetrieveActionPtr& action)
{
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

TensorBranch*
TensorBlock::get_obsolete_branch() const
{
    return obsolete_branch_;
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

float
TensorBlock::get_max_log() const
{
    return maxlog_;
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
TensorBlock::in_cache() const
{
    return metadata_block_->data();
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
TensorBlock::init_in_core_no_sort()
{
#if LEAK_CHECK
    dout << "allocating metadata block on tensor "
        << parent_tensor_->get_name() << " in init" << endl;
#endif

    char* metadata_data = YetiRuntime::malloc(metadata_block_size_);
    metadata_block_ = new InCoreBlock(metadata_data, metadata_block_size_);
    init_mempool();
    init_branch();


    typedef TensorControllerTemplate<
            ValidBranchController,
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
            DeleteAllDataClear
        > mem_controller_t;

    read_controller_ = new mem_controller_t(this);

    write_controller_ = new mem_controller_t(this);

    accumulate_controller_ = new mem_controller_t(this);

    verbatim_controller_ = new mem_controller_t(this);
}

void
TensorBlock::init_in_core_resort()
{
#if MALLOC_CACHE_BLOCK
    char* metadata_data = 0;
    metadata_block_ = new InCoreBlock(metadata_data, metadata_block_size_);
    metadata_mempool_ = new MemoryPool(metadata_block_->size(), metadata_data);
#else
    CachedStorageBlock* block = new CachedStorageBlock(this, parent_tensor_->get_metadata_cache());
    metadata_block_ = block;
    init_mempool();
#endif

#if MALLOC_CACHE_BLOCK
    read_controller_ = new TensorControllerTemplate<
                        ValidBranchController,
                        ResetMempool,
                        SortedBranchRetrieve,
                        ReallocateDataControllers,
                        SortDataControllers,
                        ClearBranchFlush,
                        DoNothingBranchRelease,
                        DoNothingMetaDataRetrieve,
                        DoNothingDataRetrieve,
                        DeleteAllDataStorageFlush,
                        ClearMetaDataOnObsolete,
                        DoNothingSync,
                        ObsoleteAllDataClear
                    >(this);
#else
    read_controller_ = new TensorControllerTemplate<
                        ValidBranchController,
                        ResetMempool,
                        SortedBranchRetrieve,
                        ReallocateDataControllers,
                        SortDataControllers,
                        ClearBranchFlush,
                        DoNothingBranchRelease,
                        DoNothingMetaDataRetrieve,
                        DoNothingDataRetrieve,
                        DeleteAllDataStorageFlush,
                        ClearMetaDataOnObsolete,
                        DoNothingSync,
                        ObsoleteAllDataClear
                    >(this);
#endif

    write_controller_ = new AbortWriteTensorController(this);



#if 0
    accumulate_controller_ = new AbortAccumulateTensorController(this);
#else
    accumulate_controller_ = new TensorControllerTemplate<
                            ValidBranchController,
                            ResetMempool,
                            NewBranchRetrieve,
                            DoNothingDataControllerRetrieve,
                            DoNothingDataControllerInit,
                            SortedAccumulateBranchFlush,
                            DoNothingBranchRelease,
                            DoNothingMetaDataRetrieve,
                            DoNothingDataRetrieve,
                            DeleteAllDataStorageFlush,
                            ClearMetaDataOnObsolete,
                            FlushOnSync,
                            ObsoleteAllDataClear
                        >(this);
#endif

    verbatim_controller_ = new AbortVerbatimTensorController(this);
}

void
TensorBlock::init_action()
{
    if (permutationally_unique_)
        init_action_no_sort();
    else
        reinit(Tensor::default_storage);
}

void
TensorBlock::init_action_no_sort()
{
    CachedStorageBlock* block = new CachedStorageBlock(this, parent_tensor_->get_metadata_cache());
    metadata_block_ = block;
    init_mempool();

    read_controller_ = new TensorControllerTemplate<
                        ValidBranchController,
                        ResetMempool,
                        ActionBranchRetrieve,
                        DoNothingDataControllerRetrieve,
                        DoNothingDataControllerInit,
                        ClearBranchFlush,
                        DoNothingBranchRelease,
                        DoNothingMetaDataRetrieve,
                        DoNothingDataRetrieve,
                        DeleteAllDataStorageFlush,
                        AbortOnObsolete,
                        DoNothingSync,
                        ObsoleteAllDataClear
                    >(this);

    write_controller_ = new AbortWriteTensorController(this);

    accumulate_controller_ = new AbortAccumulateTensorController(this);

    verbatim_controller_ = new AbortVerbatimTensorController(this);
}

void
TensorBlock::init_action_resort()
{
    CachedStorageBlock* block = new CachedStorageBlock(this, parent_tensor_->get_metadata_cache());
    metadata_block_ = block;
    init_mempool();

    read_controller_ = new TensorControllerTemplate<
                        ValidBranchController,
                        ResetMempool,
                        SortedBranchRetrieve,
                        ReallocateDataControllers,
                        SortDataControllers,
                        ClearBranchFlush,
                        DoNothingBranchRelease,
                        DoNothingMetaDataRetrieve,
                        DoNothingDataRetrieve,
                        DeleteAllDataStorageFlush,
                        ClearMetaDataOnObsolete,
                        DoNothingSync,
                        ObsoleteAllDataClear
                    >(this);

    write_controller_ = new AbortWriteTensorController(this);

    accumulate_controller_ = new AbortAccumulateTensorController(this);

    verbatim_controller_ = new AbortVerbatimTensorController(this);
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
    DiskBuffer* buffer = parent_tensor_->get_disk_buffer();
    metadata_block_ = new LocalDiskBlock(
                        this,
                        parent_tensor_->get_metadata_cache(),
                        buffer
                      );
    metadata_block_->retrieve();

    init_mempool();
    init_branch();

    read_controller_ = new TensorControllerTemplate<
                        ValidBranchController,
                        DoNothingMempool,
                        RealignMemoryPoolBranchRetrieve,
                        ReuseDataControllers,
                        DoNothingDataControllerInit,
                        ClearBranchFlush,
                        DoNothingBranchRelease,
                        DoNothingMetaDataRetrieve,
                        DoNothingDataRetrieve,
                        ClearAllDataStorageFlush,
                        AbortOnObsolete,
                        DoNothingSync,
                        ObsoleteAllDataClear
                    >(this);
                    
    write_controller_ = new TensorControllerTemplate<
                        ValidBranchController,
                        DoNothingMempool,
                        RealignMemoryPoolBranchRetrieve,
                        ReuseDataControllers,
                        DoNothingDataControllerInit,
                        CommitBranchFlush,
                        DoNothingBranchRelease,
                        DoNothingMetaDataRetrieve,
                        DoNothingDataRetrieve,
                        CommitAllDataStorageFlush,
                        AbortOnObsolete,
                        DoNothingSync,
                        ObsoleteAllDataClear
                    >(this);
                    
    accumulate_controller_ = new TensorControllerTemplate<
                        ValidBranchController,
                        DoNothingMempool,
                        RealignMemoryPoolBranchRetrieve,
                        ReuseDataControllers,
                        DoNothingDataControllerInit,
                        CommitBranchFlush,
                        DoNothingBranchRelease,
                        DoNothingMetaDataRetrieve,
                        DoNothingDataRetrieve,
                        CommitAllDataStorageFlush,
                        AbortOnObsolete,
                        DoNothingSync,
                        ObsoleteAllDataClear
                    >(this);

    verbatim_controller_ = new AbortVerbatimTensorController(this);

    metadata_block_->release();
}

void
TensorBlock::init_on_disk_resort()
{
   metadata_block_ = new CachedStorageBlock(this, parent_tensor_->get_metadata_cache());
   init_mempool();
   read_controller_ = new TensorControllerTemplate<
                        ValidBranchController,
                        ResetMempool,
                        SortedBranchRetrieve,
                        ReallocateDataControllers,
                        SortDataControllers,
                        ClearBranchFlush,
                        DoNothingBranchRelease,
                        DoNothingMetaDataRetrieve,
                        DoNothingDataRetrieve,
                        DeleteAllDataStorageFlush,
                        ClearMetaDataOnObsolete,
                        DoNothingSync,
                        ObsoleteAllDataClear
                    >(this);

    write_controller_ = new AbortWriteTensorController(this);

    accumulate_controller_ = new AbortAccumulateTensorController(this);

    verbatim_controller_ = new AbortVerbatimTensorController(this);
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

    read_controller_ = new TensorControllerTemplate<
                        ValidBranchController,
                        ResetMempool,
                        ConfigureElementComputerBranchRetrieve,
                        DoNothingDataControllerRetrieve,
                        DoNothingDataControllerInit,
                        ClearBranchFlush,
                        DoNothingBranchRelease,
                        RecomputeMetaDataRetrieve,
                        RecomputeDataRetrieve,
                        DeleteAllDataStorageFlush,
                        ClearMetaDataOnObsolete,
                        DoNothingSync,
                        ObsoleteAllDataClear  
                    >(this);
    
    write_controller_ = new AbortWriteTensorController(this);
    
    accumulate_controller_ = new AbortAccumulateTensorController(this);

    verbatim_controller_ = new AbortVerbatimTensorController(this);
}

void
TensorBlock::init_recomputed_resort()
{
    metadata_block_ = new CachedStorageBlock(this, parent_tensor_->get_metadata_cache());
    init_mempool();

    read_controller_ = new TensorControllerTemplate<
                        ValidBranchController,
                        ResetMempool,
                        ConfigureElementComputerBranchRetrieve,
                        DoNothingDataControllerRetrieve,
                        DoNothingDataControllerInit,
                        ClearBranchFlush,
                        DoNothingBranchRelease,
                        RecomputeMetaDataRetrieve,
                        RecomputeDataRetrieve,
                        DeleteAllDataStorageFlush,
                        ClearMetaDataOnObsolete,
                        DoNothingSync,
                        ObsoleteAllDataClear  
                    >(this);

    write_controller_ = new AbortWriteTensorController(this);
    
    accumulate_controller_ = new AbortAccumulateTensorController(this);

    verbatim_controller_ = new AbortVerbatimTensorController(this);

}

bool
TensorBlock::is_permutationally_unique() const
{
    return permutationally_unique_;
}

bool
TensorBlock::is_remote_block() const
{
    bool remote_block = parent_tensor_->get_distribution_type() == Tensor::distributed
              && YetiRuntime::me() != block_number_ & YetiRuntime::nproc();
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
TensorBlock::set_max_log(float maxlog)
{
    maxlog_ = maxlog;
}

void
TensorBlock::_obsolete()
{
    //it's possible this is not symmetry unique and was never used
    if (controller_)
    {
#if DEBUG_CACHE_USAGE
        dout << "obsoleting block " << get_index_string()
            << " on tensor " << parent_tensor_->get_name() 
            << " at location " << (void*) this
            << endl;
#endif
        controller_->obsolete();
    }
}

void
TensorBlock::print(std::ostream& os)
{
    NormElementOp* op = new NormElementOp;
    element_op(op);
    double norm = op->norm();


    set_synced();

    os    << "Degeneracy: " << degeneracy_ << endl
        << stream_printf("max log: %8.4f    norm: %8.4e", maxlog_, norm)
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
TensorBlock::_retrieve()
{
#if DEBUG_CACHE_USAGE
    dout << nstorage_blocks_ << " storage blocks on beginning retrieve of block " << get_index_string()
          << " on tensor " << parent_tensor_->get_name() 
          << " at location " << (void*) this
          << endl;
#endif

    controller_->validate();
    bool in_core = metadata_block_->retrieve();
    if (!in_core)
    {
#if DEBUG_CACHE_USAGE
        CachedStorageBlock* block = static_cast<CachedStorageBlock*>(metadata_block_.get());
        dout << "allocated metadata entry " << block->get_entry()->offset
             << " to block " << get_index_string()
            << " on tensor " << parent_tensor_->get_name() 
            << " at location " << (void*) this
            << endl;

        dout << "sorted branch retrieve for " << nstorage_blocks_ << " pre-existing blocks" 
            << " on tensor " << parent_tensor_->get_name() 
            << " at location " << (void*) this
            << endl;
#endif

        branch_ = reinterpret_cast<TensorBranch*>(metadata_block_->data());
        metadata_mempool_->set(metadata_block_->data());
        branch_->set_metadata_mempool(metadata_mempool_.get());
        controller_->retrieve(branch_);
    }
    else
    {
#if DEBUG_CACHE_USAGE
        dout << "regular branch retrieve on block "
             << ClassOutput<const uli*>::str(nindex_, indices_) 
            << " on tensor " << parent_tensor_->get_name() 
             << " at location " << (void*) this
             << endl;


        CachedStorageBlock* cache_block = dynamic_cast<CachedStorageBlock*>(metadata_block_.get());
        if (cache_block)
        {
            dout << "retrieve metadata entry " << cache_block->get_entry()->offset
            << " on block " << get_index_string() 
            << " on tensor " << parent_tensor_->get_name()
            << " at location " << (void*) this
            << endl;
        }
#endif

        //just retrieve the data blocks
        DataStorageNode* node = first_node_;
        while (node)
        {
            node->block->retrieve();
            node = node->next;
        }

    }

#if DEBUG_CACHE_USAGE
    dout << nstorage_blocks_ << " storage blocks on ending retrieve"
        << " of block " << get_index_string()
        << " on tensor " << parent_tensor_->get_name() 
        << " at location " << (void*) this
        << endl;
#endif

    //after retrieving, we are synced with the parent storage location
    set_synced();
}

void
TensorBlock::_release()
{
    controller_->release(branch_);

#if DEBUG_CACHE_USAGE
    dout << "there are " << nstorage_blocks_ << " storage blocks on release"
          << " of block " << get_index_string()
          << " on tensor " << parent_tensor_->get_name() 
          << " at location " << (void*) this
          << endl;
#endif

    metadata_block_->release();
    DataStorageNode* node = first_data_node();
    while (node)
    {
        node->block->release();

        CachedStorageBlock* data_block = dynamic_cast<CachedStorageBlock*>(node->block);
        if (data_block && data_block->get_entry() && data_block->get_entry()->is_pulled())
        {
            raise(
               SanityCheckError,
               "released cached storage block, but cache entry is still pulled"
            );
        }

        node = node->next;
    }
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
        return;

    if (is_flushed_)
    {
        cerr << "somehow flushed from elsewhere" << endl;
        abort();
    }

#if DEBUG_CACHE_USAGE
    dout << "start flushing block " << get_index_string()
         << " on tensor " << parent_tensor_->get_name()
         << " at location " << (void*) this
         << endl;

    dout << nstorage_blocks_ << " storage block on begin flush of block " << get_index_string()
          << " on tensor " << parent_tensor_->get_name()
          << " at location " << (void*) this
          << endl;
#endif

    is_flushed_ = true;
    controller_->flush(branch_);
    controller_->flush_data();
    is_flushed_ = false;

#if DEBUG_CACHE_USAGE
    dout << nstorage_blocks_ << " storage block on end flush of block " << get_index_string()
          << " on tensor " << parent_tensor_->get_name() 
          << " at location " << (void*) this
          << endl;

    dout << "done flushing block " << get_index_string()
        << " on tensor " << parent_tensor_->get_name()
        << " at location " << (void*) this
        << endl;
#endif
}

/********
End Cachable interface definition
*********/

void
TensorBlock::recompute_uniqueness()
{
    permutationally_unique_ = parent_tensor_->is_unique(indices_);
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
TensorBlock::set_degeneracy(usi degeneracy)
{
    degeneracy_ = degeneracy;
}

bool
TensorBlock::set_accumulate_mode()
{
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
    //already in correct mode
    if (controller_ == read_controller_)
        return;
    
    if (!is_synced())
    {
        cerr << "Cannot switch to read mode.  Tensor block has"
             << " not been synced from previous operation" << endl;
        abort();
    }
    controller_ = read_controller_;
}

void
TensorBlock::set_write_mode()
{
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
TensorBlock::set_obsolete_branch(TensorBranch* branch)
{
    obsolete_branch_ = branch;
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

    controller_->sync();
    set_synced();
}

void
TensorBlock::sync_max_log()
{
    TensorBlock* block = get_symmetry_unique_block();
    maxlog_ = block ? block->maxlog_ : LOG_ZERO;
}

void
TensorBlock::update()
{
    if (!is_permutationally_unique())
    {
        cerr << "Cannot update a block which is not" 
             << " permutationally unique" << endl;
        abort();
    }
    branch_->get_node()->update();

    maxlog_ = branch_->get_node()->get_max_log();

    reset_degeneracy();
}

void
TensorBlock::reset_degeneracy()
{
    degeneracy_ = 1;
}

TensorDataController::TensorDataController(
    TensorBranch* branch,
    StorageBlock* storage_block
) :
    size_(0),
    remaining_(0),
    branch_(branch),
    first_data_node_(0),
    last_data_node_(0),
    next(0)
{
    MemoryPool* mempool = branch->get_metadata_mempool();

    size_ = storage_block->size();
    remaining_ = size_;

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

    size_t current_pos = size_ - remaining_;
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
TensorDataController::get_size() const
{
    return size_;
}

size_t
TensorDataController::get_remaining() const
{
    return remaining_;
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

TensorController::TensorController(TensorBlock* block)
    : parent_block_(block)
{
}
