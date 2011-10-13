#include "tensor.h"
#include "tensorbranch.h"
#include "tensoraction.h"
#include "tensorblock.h"
#include "node.h"
#include "index.h"
#include "class.h"
#include "exception.h"
#include "permutation.h"
#include "env.h"
#include "data.h"
#include "runtime.h"
#include "dataimpl.h"
#include "threadimpl.h"
#include "filler.h"
#include "sortimpl.h"
#include "tensorimpl.h"
#include "branchpolicies.h"
#include "datapolicies.h"
#include "metadatapolicies.h"

#include <libsmartptr/strop.h>

using namespace yeti;
using namespace std;

void
DoNothingMempool::retrieve(TensorBlock* block)
{
}

void
ResetMempool::retrieve(TensorBlock* block)
{
    TensorBranch* branch = block->get_branch();
    MemoryPool* mempool = branch->get_metadata_mempool();
    mempool->reset();
    TensorBranch* newbranch = new (mempool) TensorBranch(
                            block->get_indices(),
                            block->get_depth(),
                            block
                           );

    if (newbranch != branch)
    {
        raise(SanityCheckError, "tensor branch not reset properly");
    }
}

void
FlushOldBranchRenew::renew(TensorBlock* block)
{
    block->flush_from_cache();
    block->_initialize();
    block->set_initialized(true);
    block->set_finalized(false);
    block->_retrieve();
}

void
RemoteAccumulateRenew::renew(TensorBlock* block)
{
    block->wait_on_send(); //wait on any previous sends
    ZeroBranchRenew::renew(block);
}

void
ZeroBranchRenew::renew(TensorBlock* block)
{
    TensorBranch* branch = block->get_branch();

    MemoryPool* mempool = branch->get_metadata_mempool();
    mempool->reset();

    /** reinitialize an empty branch */
    TensorBranch* newbranch = new (mempool) TensorBranch(
                            block->get_indices(),
                            block->get_depth(),
                            block
                           );

    /** Allocate new controllers for each of the already existing storage blocks */
    DataStorageNode* node = block->first_data_node();
    while(node)
    {
        StorageBlock* data_block = node->block;
        branch->allocate_data_controller(data_block);
        node = node->next;
    }
}

void
DoNothingBranchRenew::renew(TensorBlock* block)
{
}

void
ActionBranchRetrieve::retrieve(TensorBlock* block)
{
    TensorRetrieveAction* action = block->get_retrieve_action();
    if (!action)
    {
        cerr << "Action tensor has null action" << endl;
        abort();
    }

    action->accumulate_to(block);
}

void
ConfigureElementComputerBranchRetrieve::retrieve(TensorBlock* block)
{
    ThreadedTensorElementComputer* threaded_filler =
        block->get_parent_tensor()->get_element_computer();

    if (!threaded_filler)
    {
        cerr << "Recomputed tensor has null element computer" << endl;
        abort();
    }
    
    uli threadnum = YetiRuntime::get_thread_number();
    TensorElementComputer* filler = threaded_filler->get_computer(threadnum);
    block->get_branch()->set_element_computer(filler);
}

void
ConfigureElementComputerAndSortBranchRetrieve::retrieve(TensorBlock* block)
{
    ConfigureElementComputerBranchRetrieve::retrieve(block);
    SortedBranchRetrieve::retrieve(block);
}

void
DoNothingBranchRetrieve::retrieve(TensorBlock* block)
{
}

void
ValidBranchController::validate(TensorBlock* block)
{
}

void
AbortReadBranchValidation::validate(TensorBlock* block)
{
    cerr << "Block " << block->get_index_string() << " on tensor "
         << block->get_parent_tensor()->get_name()
         << " cannot be configured for read" << endl;
    abort();
}

void
AbortWriteBranchValidation::validate(TensorBlock* block)
{
    cerr << "Block " << block->get_index_string() << " on tensor "
         << block->get_parent_tensor()->get_name()
         << " cannot be configured for write" << endl;
    abort();
}

void
AbortAccumulateBranchValidation::validate(TensorBlock* block)
{
    cerr << "Block " << block->get_index_string() << " on tensor "
         << block->get_parent_tensor()->get_name()
         << " cannot be configured for accumulate" << endl;
    abort();
}

void
AbortVerbatimBranchValidation::validate(TensorBlock* block)
{
    cerr << "Block " << block->get_index_string() << " on tensor "
         << block->get_parent_tensor()->get_name()
         << " cannot be configured for verbatim" << endl;
    abort();
}

void
RealignMemoryPoolBranchRetrieve::retrieve(
    TensorBlock* block
)
{
    TensorBranch* branch = block->get_branch();
    TensorBranch* old_branch = branch->get_last_branch_location();
    TensorBranch* new_branch = branch;

    if (!old_branch)
    {
        raise(SanityCheckError, "no old branch to realign");
    }
    new_branch->realign_memory_pool(old_branch, new_branch);
}

void
NewBranchRetrieve::retrieve(TensorBlock* block)
{
    /** This should have been reinitialized by memory pool reset */

    //TensorBlock* block = branch->get_parent_block();
    //block->init_branch();
}

void
RemoteBlockBranchRetrieve::retrieve(TensorBlock* block)
{
    block->set_waiting(true);
    YetiRuntime::get_messenger()->send_data_request(
        Message::TensorBranch,
        block->get_node_number(),
        block->get_malloc_number()
    );
    block->wait_on_send();
}

void
SortedBranchRetrieve::retrieve(
    TensorBlock* block
)
{
    TensorBlock* unique_block = block->get_symmetry_unique_block();
    retrieve(block, unique_block);
}

void
SortedBranchRetrieve::retrieve(
    TensorBlock* block,
    TensorBlock* unique_block,
    Permutation* perm
)
{
    if (!unique_block)
    {
        ResetMempool reset;
        reset.retrieve(block);
        return;
    }

    cout << "sorting " << unique_block->get_index_string()
        << " into " << block->get_index_string() << endl;

    TensorBranch* branch = block->get_branch();
    MemoryPool* mempool = branch->get_metadata_mempool();

    unique_block->retrieve_read();
    mempool->memcpy(unique_block->get_metadata_mempool());
    //reset the parent - gets messed up by memcpy
    branch->set_parent(block);
    branch->set_metadata_mempool(mempool);

    Sort* sort = new Sort(perm);
    unique_block->get_branch()
        ->sort_metadata_into(sort, branch);

    //now that the sort of metadata has taken place,
    //realign the memory pool
    branch->realign_memory_pool(
        unique_block->get_branch(),
        branch
    );

    unique_block->release_read();
    delete sort;
}

void
SortedBranchRetrieve::retrieve(
    TensorBlock* block,
    TensorBlock* unique_block
)
{
    retrieve(block, unique_block, block->get_resort_permutation());
}


void
DoNothingDataControllerRetrieve::retrieve(
    TensorBlock* block
)
{
}

void
ReuseDataControllers::retrieve(TensorDataController* controller)
{
    DataNode* node = controller->first();
    char* data = controller->get_data();
    while(node)
    {
         node->set_data_block(data);
         node = node->next_node;
    }
}


void
ReuseDataControllers::retrieve(
    TensorBlock* block
)
{
    TensorDataController* controller = block->get_branch()->first_data_controller();
    DataStorageNode* node = block->first_data_node();
    while(controller)
    {
        StorageBlock* block = node->block;
        block->retrieve();
        controller->set_data(block->data());
        retrieve(controller);
        node = node->next;
        controller = controller->next;
    }
}

void
MemsetDataControllers::retrieve(
    TensorBlock* block
)
{
    TensorDataController* controller = block->get_branch()->first_data_controller();
    DataStorageNode* node = block->first_data_node();
    while(controller)
    {
        StorageBlock* block = node->block;
        block->memset();
        node = node->next;
        controller = controller->next;
    }
}

void
SortDataControllers::retrieve(
    TensorBlock* block,
    TensorBlock* unique_block,
    Permutation* perm
)
{
    unique_block->retrieve_read();
    Sort* sort = new Sort(perm);
    unique_block->get_branch()->sort_data_into(sort, block->get_branch());
    unique_block->release_read();

    delete sort;
}

void
SortDataControllers::retrieve(
    TensorBlock* block
)
{
    TensorBlock* unique_block = block->get_symmetry_unique_block();
    if (!unique_block)
        return;

    retrieve(block, unique_block, block->get_resort_permutation());
}

void
ReallocateDataControllers::retrieve(
    TensorBlock* block
)
{
    TensorDataController* controller = block->get_branch()->first_data_controller();
    block->clear_storage();
    while(controller)
    {
        StorageBlock* storage_block = block->allocate_data_storage_block();
        controller->set_data(storage_block->data());
        controller = controller->next;
    }
}

void
DoNothingDataControllerInit::retrieve(TensorBlock* block)
{
}

void
DoNothingBranchFlush::flush(TensorBlock* block)
{
}

void
DoNothingPreflush::preflush(TensorBlock* block)
{
    //do nothing
}

void
ClearBranchFlush::flush(TensorBlock* block)
{
    StorageBlock* storage_block = block->get_metadata_storage_block();
    storage_block->clear();
}

void
RemoteAccumulateFlush::preflush(TensorBlock* block)
{
    YetiRuntime::get_messenger()->queue_send_request(
        Message::TensorBranch,
        Message::Accumulate,
        block->get_node_number(),
        block
    );
}

void
RemoteAccumulateFlush::flush(TensorBlock* block)
{
    if (!block->has_send_status()) //I need to do the flushing
    {
        SendStatus* status = block->send_data(
            YetiRuntime::get_messenger(),
            Message::TensorBranch,
            Message::Accumulate,
            block->get_node_number()
        );

        status->wait();

        block->reset_send_status();
    }
    else //if there was a previously exisiting send, we have already waited on it
    {
    }
    ClearBranchFlush::flush(block);
}

void
SortedAccumulateBranchFlush::flush(TensorBlock* block)
{   
    if (block->in_destructor())
        return;

    TensorBlock* unique_block = block->get_symmetry_unique_block();
    Permutation* unsort_perm = block->get_resort_permutation()->inverse();
    Sort* sort = new Sort(unsort_perm);
    unique_block->retrieve_accumulate();
    unique_block->accumulate(block, 1.0, sort);
    unique_block->release_accumulate();
    delete sort;
}

void
CommitBranchFlush::flush(TensorBlock* block)
{
    StorageBlock* storage_block = block->get_metadata_storage_block();
    storage_block->commit();
    storage_block->clear();
}

void
DoNothingBranchRelease::release(TensorBlock* block)
{
}

void
SetFinalizedBranchRelease::release(TensorBlock* block)
{
    block->finalize();
}

void
CacheBranchRelease::release(TensorBlock* block)
{
    block->get_metadata_storage_block()->release();
    DataStorageNode* node = block->first_data_node();
    while(node)
    {
        StorageBlock* data_block = node->block;
        data_block->release();
        node = node->next;
    }
}

void
ThreadAccumulateBranchRelease::release(TensorBlock* block)
{
    TensorBlock* parent = block->get_parent_tensor()
                               ->get_block(block->get_block_number());
    Sort* no_sort = 0;
    parent->retrieve_accumulate();
    parent->accumulate(block, 1.0, no_sort);
    parent->release_accumulate();
}

void
RemoteAccumulateBranchRelease::release(TensorBlock* block)
{
    YetiRuntime::get_messenger()
            ->queue_send_request(
                Message::TensorBranch,
                Message::Accumulate,
                block->get_node_number(),
                block,
                &TensorBlock::release_callback
            );

#if 0
    SendStatus* status = block->send_data(
        YetiRuntime::get_messenger(),
        Message::TensorBranch,
        Message::Accumulate,
        block->get_node_number()
    );
    block->set_send_status(status);
    status->wait(); //leave the status so the flush knows we sent our data already

    block->get_metadata_storage_block()->release();
    DataStorageNode* node = block->first_data_node();
    while(node)
    {
        StorageBlock* data_block = node->block;
        data_block->release();
        node = node->next;
    }
#endif
}

void
DoNothingMetaDataRetrieve::retrieve(MetaDataNode* mdnode)
{
}

void
RecomputeMetaDataRetrieve::retrieve(MetaDataNode* mdnode)
{
    mdnode->init_subnodes();
    
    TileMap::iterator it(mdnode->begin());
    TileMap::iterator stop(mdnode->end());
    for ( ; it != stop; ++it)
    {
        TileNode* node = *it;
        node->set_max_log(LOG_NONZERO);
    }
    
}

void
DoNothingMetaDataRelease::release(MetaDataNode* mdnode)
{
}


void
DoNothingDataRetrieve::retrieve(
    DataNode* node,
    const MetaDataNode* parent,
    const uli* indexset
)
{
}

bool
DoNothingDataRetrieve::need_recompute() const
{
    return false;
}

template <typename data_t>
void
RecomputeDataRetrieve::compute(
    TensorElementComputer* filler,
    DataNode* node,
    const uli* indexset
)
{
    data_t* data = reinterpret_cast<data_t*>(node->data());
    filler->compute(indexset, data, node->nelements());
}

void
RecomputeDataRetrieve::retrieve(
    DataNode* node,
    const MetaDataNode* parent,
    const uli* indexset
)
{
    TensorBranch* branch = parent->get_parent_branch();
    TensorElementComputer* filler = branch->get_element_computer();
    TemplateInfo::type_t elem_type = branch->element_type();

    branch->allocate(node);

    data_type_switch(
        elem_type,
        this->compute,
        filler, 
        node,
        indexset
    );
}

bool
RecomputeDataRetrieve::need_recompute() const
{
    return true;
}

void
DoNothingDataRelease::release(DataNode* node, const MetaDataNode* parent)
{
}

void
DoNothingDataStorageFlush::flush(TensorBlock* block)
{
}

void
CommitAllDataStorageFlush::flush(TensorBlock* block)
{
    DataStorageBlockAllocator* allocator = block;
    DataStorageNode* node = allocator->first_data_node();
    while(node)
    {
        StorageBlock* data_block = node->block;
        data_block->commit();
        data_block->clear();
        node = node->next;
    }
}

void
ClearAllDataStorageFlush::flush(TensorBlock* block)
{
    uli nblocks = block->get_num_data_storage_blocks();
    DataStorageBlockAllocator* allocator = block;
    DataStorageNode* node = allocator->first_data_node();
    while(node)
    {
        StorageBlock* data_block = node->block;
        data_block->clear();
        node = node->next;
    }
}

void
DeleteAllDataStorageFlush::flush(TensorBlock* block)
{
    block->clear_storage();
}

void
AbortOnObsolete::obsolete(TensorBlock* block)
{
    cerr << "invalid obsolete declaration on tensor block" << endl;
    abort();
}

void
ClearMetaDataOnObsolete::obsolete(TensorBlock* block)
{
    block->get_metadata_storage_block()->obsolete(block);
    DataStorageNode* node = block->first_data_node();
    while (node)
    {
        node->block->obsolete(block);
        node = node->next;
    }
}

void
DoNothingSync::sync(TensorBlock* block)
{
}

void
RemoteAccumulateSync::sync(TensorBlock* block)
{
    /** just wait until the block has finished its send */
    block->wait_on_send();
}

void
AbortOnSync::sync(TensorBlock* block)
{
    cerr << "block should never be synced" << endl;
    abort();
}

void
FlushOnSync::sync(TensorBlock* block)
{
    if (block->is_cached())
    {
        block->lock();
        if (!block->is_cached()) //race condition!
        {
            //block got flushed in the meanwhile
            block->unlock();
            return;
        }
        block->flush_from_cache();
        block->unlock();
    }

    block->obsolete();
    //not in cache.. our work is done
}

void
DeleteAllDataClear::clear(
    TensorBlock* tensor_block
)
{
    tensor_block->free_metadata();
    DataStorageNode* node = tensor_block->first_data_node();
    while (node)
    {
        StorageBlock* storage_block = node->block;
        tensor_block->free_data(storage_block->data(), storage_block->size());
        node = node->next;
    }
    tensor_block->clear_storage();
}

void
ObsoleteAllDataClear::clear(
    TensorBlock* tensor_block
)
{
    tensor_block->get_metadata_storage_block()->obsolete(tensor_block);
    tensor_block->clear_metadata();
    DataStorageNode* node = tensor_block->first_data_node();
    while (node)
    {
        StorageBlock* storage_block = node->block;
        storage_block->obsolete(tensor_block);
        node = node->next;
    }
}

void
DoNothingOutOfCorePrefetch::prefetch(TensorBlock* block)
{
}

void
ParentBlockReadPrefetch::prefetch(TensorBlock* block)
{
    cout << "parent block prefetch" << endl;
    block->get_symmetry_unique_block()->prefetch_read();
}

void
DoNothingInCorePrefetch::prefetch(
    TensorBlock* current_block,
    TensorBlock* prev_block
)
{

}

void
ResortInCorePrefetch::prefetch(
    TensorBlock* current_block,
    TensorBlock* prev_block
)
{
    //check to see if blocks are related by symmetry
    if      (current_block->get_task_owner_number() != prev_block->get_task_owner_number())
        return;
    else if (current_block->get_block_number() == prev_block->get_block_number())
        return; //previous block is the same!

    //go ahead and do all the sorting work
    current_block->retrieve();
    prev_block->retrieve();


    if (current_block->get_num_data_storage_blocks() > 0)
    { //this is still in corejust use what we have
        prev_block->release();
        current_block->release();
        return;
    }

    Permutation* sort_perm = current_block->get_resort_permutation()
                                ->product(prev_block->get_resort_permutation()->inverse());

    SortedBranchRetrieve::retrieve(current_block, prev_block, sort_perm);
    ReallocateDataControllers::retrieve(current_block);
    SortDataControllers::retrieve(current_block, prev_block, sort_perm);

    prev_block->release();
    current_block->release();
}

void
RemoteBlockPrefetch::prefetch(TensorBlock* block)
{
    block->set_waiting(true);
    block->YetiRuntimeObject::initialize();
    YetiRuntime::get_messenger()->send_data_request(
        Message::TensorBranch,
        block->get_node_number(),
        block->get_malloc_number()
    );
}


#ifdef redefine_size_t
#define size_t custom_size_t
#endif