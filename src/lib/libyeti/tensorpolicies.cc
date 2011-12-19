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

#define DEBUG_TENSOR_BLOCK_POLICIES 0

#define DEBUG_PARALLEL 0

void
DoNothingMempool::retrieve(TensorBlock* block)
{
}

void
ResetMempool::retrieve(TensorBlock* block)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("reset mempool");
#endif

#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Reset mempool retrieve on block %s on node %d\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
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
        yeti_throw(SanityCheckError, "tensor branch not reset properly");
    }
}

void
FlushOldBranchRenew::renew(TensorBlock* block)
{
#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Flush old branch renew on block %s on node %d\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
    block->flush_from_cache();
    block->_initialize();
    block->set_initialized(true);
    block->set_finalized(false);
    block->_retrieve();
}

void
RemoteAccumulateRenew::renew(TensorBlock* block)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("remote accumulate renew");
#endif

#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Remote accumulate renew on block %s on node %d\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
    try{
        block->wait_on_send();
    } catch (int e) {
        if (e == TENSOR_WAIT_ERROR)
            cerr << stream_printf("Remote accumulate renew is hanging on %s\n",
                        block->get_block_name().c_str());
        throw e;
    }
    ZeroBranchRenew::renew(block);
}

void
CacheBranchRenew::renew(TensorBlock* block)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("cache branch renew");
#endif

#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Cache branch renew on block %s on node %d\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
    block->get_metadata_storage_block()->retrieve();
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("metadata retrieve");
#endif

    /** Allocate new controllers for each of the already existing storage blocks */
    DataStorageNode* node = block->head_data_node();
    while(node)
    {
        StorageBlock* data_block = node->block;
        data_block->retrieve();
#if TRACK_TENSOR_BLOCK_HISTORY
        block->add_event("data retrieve");
#endif
        node = node->next;
    }
}

void
ZeroBranchRenew::renew(TensorBlock* block)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("zero branch renew");
#endif

#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Zero branch renew on block %s on node %d\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
    TensorBranch* branch = block->get_branch();

    block->get_metadata_storage_block()->retrieve();
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("metadata retrieve");
#endif

    MemoryPool* mempool = branch->get_metadata_mempool();
    mempool->reset();

    /** reinitialize an empty branch */
    TensorBranch* newbranch = new (mempool) TensorBranch(
                            block->get_indices(),
                            block->get_depth(),
                            block
                           );

    /** Allocate new controllers for each of the already existing storage blocks */
    DataStorageNode* node = block->head_data_node();
    while(node)
    {
        StorageBlock* data_block = node->block;
        data_block->retrieve();
#if TRACK_TENSOR_BLOCK_HISTORY
        block->add_event("data retrieve");
#endif
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
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("action branch retrieve");
#endif
    TensorRetrieveAction* action = block->get_retrieve_action();
    if (!action)
    {
        cerr << "Action tensor has null action" << endl;
        abort();
    }

    action->accumulate_to(block);
}

void
ConfigureElementComputerRenew::renew(TensorBlock* block)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("configure element computer renew");
#endif
    CacheBranchRenew::renew(block);
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
ConfigureElementComputerRetrieve::retrieve(TensorBlock* block)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("configure element computer retrieve");
#endif
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
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("configure element computer and sort retrieve");
#endif
    ConfigureElementComputerRetrieve::retrieve(block);
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
    throw TENSOR_BLOCK_POLICY_EXCEPTION;
}

void
AbortWriteBranchValidation::validate(TensorBlock* block)
{
    cerr << "Block " << block->get_index_string() << " on tensor "
         << block->get_parent_tensor()->get_name()
         << " cannot be configured for write" << endl;
    throw TENSOR_BLOCK_POLICY_EXCEPTION;
}

void
AbortAccumulateBranchValidation::validate(TensorBlock* block)
{
    cerr << "Block " << block->get_index_string() << " on tensor "
         << block->get_parent_tensor()->get_name()
         << " cannot be configured for accumulate" << endl;
    throw TENSOR_BLOCK_POLICY_EXCEPTION;
}

void
AbortVerbatimBranchValidation::validate(TensorBlock* block)
{
    cerr << "Block " << block->get_index_string() << " on tensor "
         << block->get_parent_tensor()->get_name()
         << " cannot be configured for verbatim" << endl;
    throw TENSOR_BLOCK_POLICY_EXCEPTION;
}

void
RealignMemoryPoolBranchRetrieve::retrieve(
    TensorBlock* block
)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("realign memory pool");
#endif

#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Realign branch retrieve on block %s on node %d\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
    TensorBranch* branch = block->get_branch();
    TensorBranch* old_branch = branch->get_last_branch_location();
    TensorBranch* new_branch = branch;

    if (!old_branch)
    {
        yeti_throw(SanityCheckError, "no old branch to realign");
    }
    new_branch->realign_memory_pool(old_branch, new_branch);
}

void
NewBranchRetrieve::retrieve(TensorBlock* block)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("new branch retrieve");
#endif

#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("New branch retrieve on block %s on node %d\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
    /** This should have been reinitialized by memory pool reset */

    //TensorBlock* block = branch->get_parent_block();
    //block->init_branch();
}

void
RemoteBlockBranchRetrieve::retrieve(TensorBlock* block)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("remote block retrieve");
#endif

#if DEBUG_PARALLEL  || DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Remote block branch retrieve on block %s on node %d\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
    ReuseStorageBlocks initer; initer.retrieve(block);
    if (!block->is_prefetched())
    {
        block->set_remote_wait(true);
        YetiRuntime::get_messenger()->send_data_request(
            Message::TensorBranch,
            block->get_node_number(),
            block->get_malloc_number(),
            Message::ReadRequestedRetrieve
        );
        block->set_prefetched(true);
        try{
            block->wait_on_remote();
        } catch (int e) {
            if (e == TENSOR_WAIT_ERROR)
                cerr << stream_printf("Remote block branch retrieve is hanging on %s\n",
                            block->get_block_name().c_str());
            throw e;
        }
    }
}

void
SortedBranchRetrieve::retrieve(
    TensorBlock* block
)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("sorted branch retrieve");
#endif

#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Sorted branch retrieve on block %s on node %d\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
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

#if YETI_SANITY_CHECK
    if (!unique_block->has_fast_read() && !unique_block->is_retrieved())
    {
        cerr << stream_printf("Unique block is not retrieved in sorted branch retrieve\n");
        unique_block->controller_fail();
        abort();
    }
#endif
    sort(block, unique_block, perm);
}

void
SortedBranchRetrieve::sort(
    TensorBlock* block,
    TensorBlock* unique_block,
    Permutation* perm
)
{
    TensorBranch* branch = block->get_branch();
    MemoryPool* mempool = block->get_metadata_mempool();
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

#if TRACK_TENSOR_BLOCK_HISTORY
    unique_block->add_event("sort from parent block");
#endif

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
ReuseStorageBlocks::retrieve(
    TensorBlock* block
)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("reuse storage blocks");
#endif

#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Reuse storage blocks on block %s on node %d\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
    /** Allocate new controllers for each of the already existing storage blocks */
    DataStorageNode* node = block->head_data_node();
    TensorBranch* branch = block->get_branch();
    while(node)
    {
        StorageBlock* data_block = node->block;
        data_block->retrieve();
#if TRACK_TENSOR_BLOCK_HISTORY
        block->add_event("data retrieve");
#endif
        branch->allocate_data_controller(data_block);
        node = node->next;
    }
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
    TensorBlock* tensor_block
)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    tensor_block->add_event("reuse data controllers");
#endif

#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Reuse data controllers on block %s on node %d\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
    TensorDataController* controller = tensor_block->get_branch()->first_data_controller();
    DataStorageNode* node = tensor_block->head_data_node();
    while(controller)
    {
        StorageBlock* data_block = node->block;
        data_block->retrieve();
#if TRACK_TENSOR_BLOCK_HISTORY
        tensor_block->add_event("data retrieve");
#endif
        controller->set_data(data_block->data());
        retrieve(controller);
        node = node->next;
        controller = controller->next;
    }
}

void
MemsetDataControllers::retrieve(
    TensorBlock* tensor_block
)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    tensor_block->add_event("memset data controllers");
#endif

#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Memset data controllers on block %s on node %d\n",
                          tensor_block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
    TensorDataController* controller = tensor_block->get_branch()->first_data_controller();
    DataStorageNode* node = tensor_block->head_data_node();
    while(controller)
    {
        StorageBlock* data_block = node->block;
        data_block->retrieve();
#if TRACK_TENSOR_BLOCK_HISTORY
        tensor_block->add_event("data retrieve");
#endif
        data_block->memset();
        node = node->next;
        controller = controller->next;
    }
}

void
SortDataControllers::sort(
    TensorBlock* block,
    TensorBlock* unique_block,
    Permutation* perm
)
{
#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Sort data controllers on block %s on node %ld\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
    TensorIndexDescr* parent_descr = 0;
    if (block->get_block_permutation() != unique_block->get_block_permutation())
    {
        /** the descr for this tensor block will be wrong */
        parent_descr = unique_block->get_parent_tensor()->get_parent_descr();
        unique_block->get_branch()->set_descr(parent_descr);
    }

    Sort* sort = new Sort(perm);
    unique_block->get_branch()->sort_data_into(sort, block->get_branch());

    if (parent_descr) //resort the old descr
    {
        unique_block->get_branch()->set_descr(unique_block->get_descr());
    }

    delete sort;
}

void
SortDataControllers::retrieve(
    TensorBlock* block,
    TensorBlock* unique_block,
    Permutation* perm
)
{
#if YETI_SANITY_CHECK
    if (!unique_block->has_fast_read() && !unique_block->is_retrieved())
    {
        block->controller_fail();
        abort();
    }
#endif
    sort(block, unique_block, perm);
}

void
SortDataControllers::retrieve(
    TensorBlock* block
)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("sort data controllers");
#endif

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
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("reallocate data controllers");
#endif

#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Reallocate data controllers on block %s on node %d\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
    TensorDataController* controller = block->get_branch()->first_data_controller();
    DataStorageNode* node = block->head_data_node();
    TensorBranch* branch = block->get_branch();
    while(node)
    {
        StorageBlock* sblock = node->block;
        sblock->retrieve();
#if TRACK_TENSOR_BLOCK_HISTORY
        block->add_event("data retrieve");
#endif
        if (controller)
        {
            controller->set_data(sblock->data());
            controller = controller->next;
        }
        else //we have an extra data block
        {
            branch->allocate_data_controller(sblock);
        }
        node = node->next;
    }

    /** if we need more storage blocks */
    while(controller)
    {
        StorageBlock* storage_block = block->allocate_data_storage_block(controller->get_data_size());
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
ClearBranchFlush::flush(TensorBlock* block)
{
    StorageBlock* storage_block = block->get_metadata_storage_block();
    storage_block->clear();
}

void
FinalizeOnRelease::release(TensorBlock* block)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("finalize on release");
#endif
    block->finalize();
    CacheBranchRelease::release(block);
}

void
DoNothingFinalize::finalize(TensorBlock* block)
{
    //do nothing
}

void
RemoteAccumulateFinalize::finalize(TensorBlock* block)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("remote accumulate finalize");
#endif

#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Remote accumulate finalize on block %s on node %d\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
    block->send_data(
        YetiRuntime::get_messenger(),
        Message::TensorBranch,
        Message::AccumulateUnexpected,
        block->get_node_number()
    );
}

void
SortedAccumulateBranchFlush::flush(TensorBlock* block)
{   
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("sorted accumulate flush");
#endif

#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Sorted accumulate flush on block %s on node %d\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
    if (block->in_destructor())
        return;

    TensorBlock* unique_block = block->get_symmetry_unique_block();
    Permutation* unsort_perm = block->get_resort_permutation()->inverse();
    Sort* sort = new Sort(unsort_perm);
    unique_block->retrieve_accumulate();
    try{
        unique_block->accumulate(block, 1.0, sort);
    } catch (int e) {
        cerr << "Exception in sorted accumlate branch flush\n";
        cerr.flush();
        throw e;
    }
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
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("cache branch release");
#endif

#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Cache branch release on block %s on node %d\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif

#if YETI_SANITY_CHECK
    if (!block->is_locked())
    {
        cerr << stream_printf("Block %s not locked in cache branch release on thread %d on node %d\n",
            block->get_block_name().c_str(),
            YetiRuntime::get_thread_number(),
            YetiRuntime::me());
        throw TENSOR_BLOCK_POLICY_EXCEPTION;
    }
#endif

    block->get_metadata_storage_block()->release();
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("metadata release");
#endif
    DataStorageNode* node = block->head_data_node();
    while(node)
    {
        StorageBlock* data_block = node->block;
        data_block->release();
#if TRACK_TENSOR_BLOCK_HISTORY
        block->add_event("data release");
#endif
        node = node->next;
    }
}

void
ThreadAccumulateFinalize::finalize(TensorBlock* block)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("thread accumulate finalize");
#endif

#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Thread accumulate finalize on block %s on node %d\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
	YetiRuntime::start_timer("tmp block accumulate");
    TensorBlock* parent = Malloc<TensorBlock>::get_object(block->get_malloc_number());
    Sort* no_sort = 0;
    parent->retrieve_accumulate();
    try{
        parent->accumulate(block, 1.0, no_sort);
    } catch (int e) {
        cerr << "Thread accumulate finalize accumulate exception\n";
        cerr.flush();
        throw e;
    }
    parent->release_accumulate();
	YetiRuntime::stop_timer("tmp block accumulate");
}


void
DoNothingMetaDataRetrieve::retrieve(MetaDataNode* mdnode)
{
}

void
RecomputeMetaDataRetrieve::retrieve(MetaDataNode* mdnode)
{
    usi depth = mdnode->get_depth();
    TensorValueEstimater* estimator = mdnode->get_parent_branch()->get_element_computer()->get_estimater(depth);

    mdnode->init_subnodes();
    if (depth == 1)
    {
        uli indexset[NINDEX];
        TileMap::iterator it(mdnode->begin());
        TileMap::iterator stop(mdnode->end());
        TileMap* map = mdnode->get_node_map();
        uli idx = 0;
        for ( ; it != stop; ++it, ++idx)
        {
            DataNode* node = static_cast<DataNode*>(*it);
            map->indices(idx, indexset);
            float maxlog = estimator ? estimator->max_log(indexset) : LOG_NONZERO;
            node->set_max_log(maxlog);
        }
    }
    else
    {
        TileMap::iterator it(mdnode->begin());
        TileMap::iterator stop(mdnode->end());
        for ( ; it != stop; ++it)
        {
            MetaDataNode* node = static_cast<MetaDataNode*>(*it);
            const uli* indexset = node->get_indices();
            float maxlog = estimator ? estimator->max_log(indexset) : LOG_NONZERO;
            node->set_max_log(maxlog);
        }
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
    TensorValueEstimater* estimator = filler->get_estimater(1);
    float maxlog = estimator ? estimator->max_log(indexset) : LOG_NONZERO;
    if(maxlog < YetiRuntime::nonnull_cutoff) { // This whole block can be memset
#if COUNT_SCREENING_SKIPS
        YetiRuntime::increment_screening_skips(0);
#endif
#if PRINT_SCREENING_SKIPS
        Env::out0() << "Skipped computation of block by memsetting since maxlog was " << maxlog << endl;
#endif
        ::memset(data, 0, node->nelements() * sizeof(data_t));
    }
    else
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
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("commit all data");
#endif
    DataStorageBlockAllocator* allocator = block;
    DataStorageNode* node = allocator->head_data_node();
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
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("clear all data");
#endif
    uli nblocks = block->get_num_data_storage_blocks();
    DataStorageBlockAllocator* allocator = block;
    DataStorageNode* node = allocator->head_data_node();
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
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("delete all data flush");
#endif
    block->clear_storage();
}

void
AbortOnObsolete::obsolete(TensorBlock* block)
{
    cerr << stream_printf("Invalid obsolete declaration on tensor block %p %s\n",
                          block, block->get_block_name().c_str());
    cerr.flush();
    block->controller_fail();
    abort();
}

void
SetFinalizedOnObsolete::obsolete(TensorBlock* block)
{
#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << "Finalize obsoleting " << block->get_block_name() << endl;
#endif
    block->finalize();
}

void
DoNothingSync::sync(TensorBlock* block)
{
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
#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Flush on sync on block %s on node %d\n",
                          block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
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
#if TRACK_TENSOR_BLOCK_HISTORY
    tensor_block->add_event("delete all data clear");
#endif
    tensor_block->free_metadata();
    DataStorageNode* node = tensor_block->head_data_node();
    while (node)
    {
        StorageBlock* storage_block = node->block;
        tensor_block->free_data(storage_block->data(), storage_block->size());
        node = node->next;
    }
    tensor_block->clear_storage();
}

void
ClearAllData::clear(
    TensorBlock* tensor_block
)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    tensor_block->add_event("clear all data clear");
#endif
    tensor_block->get_metadata_storage_block()->clear();
    tensor_block->clear_metadata();
    DataStorageNode* node = tensor_block->head_data_node();
    while (node)
    {
        StorageBlock* storage_block = node->block;
        storage_block->clear();
        node = node->next;
    }
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
#if TRACK_TENSOR_BLOCK_HISTORY
    current_block->add_event("resort in core prefetch");
#endif

#if YETI_SANITY_CHECK
    if (current_block->get_symmetry_unique_block() != prev_block->get_symmetry_unique_block())
        yeti_throw(SanityCheckError, "Unrelated blocks in core prefetched together");
    if (current_block == prev_block)
        yeti_throw(SanityCheckError, "Identical blocks prefetched together");
#endif

#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Resort in core from %s to %s on node\n",
                          prev_block->get_block_name().c_str(),
                          current_block->get_block_name().c_str(), YetiRuntime::me());
    cout.flush();
#endif
    current_block->YetiRuntimeObject::initialize();

    //go ahead and do all the sorting work
#if YETI_SANITY_CHECK
    if (!prev_block->is_retrieved())
    {
        if (!prev_block->has_fast_read()) //a fast read block doesn't bother with retrieves
        {
            cerr << stream_printf("On me=%ld,\nCurrent block %s on node %ld\nPrev block %s on node %ld", 
                YetiRuntime::me(),
                current_block->get_block_name().c_str(), current_block->get_node_number(),
                prev_block->get_block_name().c_str(), prev_block->get_node_number()) << endl;
            yeti_throw(SanityCheckError, "Parent block is not retrieved in in core prefetch");
        }
    }
#endif

    Permutation* sort_perm = current_block->get_resort_permutation()
                                ->product(prev_block->get_resort_permutation()->inverse());

    //ensure the metadata is retrieved
    current_block->get_metadata_storage_block()->retrieve();
#if TRACK_TENSOR_BLOCK_HISTORY
    current_block->add_event("metadata retrieve");
#endif

    SortedBranchRetrieve::sort(current_block, prev_block, sort_perm);
    ReallocateDataControllers::retrieve(current_block);
    SortDataControllers::sort(current_block, prev_block, sort_perm);
    
    current_block->set_finalized(false); //renew on next retrieve
    current_block->set_prefetched(true);
}

void
RemoteBlockPrefetch::prefetch(TensorBlock* block)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("remote block prefetch");
#endif

#if DEBUG_TENSOR_BLOCK_POLICIES
    cout << stream_printf("Remote block prefetch %s\n", block->get_block_name().c_str());
    cout.flush();
#endif
    block->YetiRuntimeObject::initialize();
    CacheBranchRenew renewer; renewer.renew(block);

    block->set_remote_wait(true);
    YetiRuntime::get_messenger()->send_data_request(
        Message::TensorBranch,
        block->get_node_number(),
        block->get_malloc_number(),
        Message::ReadRequestedPrefetch
    );
    block->set_prefetched(true);
}

void
InitializePrefetch::prefetch(TensorBlock* tensor_block)
{
    if (!tensor_block->is_initialized())
    {
        tensor_block->YetiRuntimeObject::initialize();
    }
    else if (!tensor_block->get_metadata_storage_block()->is_retrieved()) //make sure all data is retrieved
    {
        tensor_block->get_metadata_storage_block()->retrieve();
        #if TRACK_TENSOR_BLOCK_HISTORY
        tensor_block->add_event("metadata retrieve");
        #endif
        DataStorageNode* node = tensor_block->head_data_node(); 
        while(node)
        {
            StorageBlock* data_block = node->block;
            data_block->retrieve();
            node = node->next;
            #if TRACK_TENSOR_BLOCK_HISTORY
            tensor_block->add_event("data retrieve");
            #endif
        }
    }
    tensor_block->set_prefetched(true);
    tensor_block->unlock();
}

void
UnlockAfterPrefetch::prefetch(TensorBlock* block)
{
    block->set_prefetched(true);
    block->unlock();
}

void
ParentBlockReadPrefetch::prefetch(TensorBlock* block)
{
#if TRACK_TENSOR_BLOCK_HISTORY
    block->add_event("remote block prefetch");
    block->get_symmetry_unique_block()->add_event("parent block prefetch");
#endif
    block->get_symmetry_unique_block()->prefetch_read();

    /** Make sure block is initialized and waiting */
    block->YetiRuntimeObject::initialize();
    CacheBranchRenew renewer; renewer.renew(block);
    block->set_prefetched(true);
    block->unlock();
}

#ifdef redefine_size_t
#define size_t custom_size_t
#endif
