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
DoNothingMempool::retrieve(TensorBranch* branch)
{
}

void
ResetMempool::retrieve(TensorBranch* branch)
{
    TensorBlock* block = branch->get_parent_block();
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
ActionBranchRetrieve::retrieve(TensorBranch* branch)
{
    TensorBlock* parent_block = branch->get_parent_block();
    TensorRetrieveAction* action = parent_block->get_retrieve_action();
    if (!action)
    {
        cerr << "Action tensor has null action" << endl;
        abort();
    }

    action->accumulate_to(parent_block);
}

void
ConfigureElementComputerBranchRetrieve::retrieve(TensorBranch* branch)
{
    ThreadedTensorElementComputer* threaded_filler =
        branch->get_parent_block()->get_parent_tensor()->get_element_computer();

    if (!threaded_filler)
    {
        cerr << "Recomputed tensor has null element computer" << endl;
        abort();
    }
    
    uli threadnum = YetiRuntime::get_thread_number();
    TensorElementComputer* filler = threaded_filler->get_computer(threadnum);
    branch->set_element_computer(filler);
}

void
ResortAndConfigureElementComputerBranchRetrieve::retrieve(TensorBranch* branch)
{
    ConfigureElementComputerBranchRetrieve::retrieve(branch);

    //loop through the permutation group and try to find
    //an already existing tensor block
    TensorBlock* block = branch->get_parent_block();
    Tensor* tensor = block->get_parent_tensor();
    PermutationGroupPtr tensor_grp = tensor->get_tensor_grp();
    PermutationGroup::iterator it(tensor_grp->begin());
    PermutationGroup::iterator stop(tensor_grp->end());
    uli indices[NINDEX];
    
    TensorBlock* unique_block = 0;
    for ( ; it != stop; ++it)
    {
        Permutation* p(*it);
        if (p->is_identity())
            continue;
        
        p->permute(block->get_indices(), indices);
        TensorBlock* test_block = tensor->get_block(indices);
        if (test_block->is_retrieved())
        {
            test_block->retrieve(); //increment retrieve count
            unique_block = test_block;
            break;
        }
    }
    
    //sort what you can into the already existing block
    if (unique_block)
        SortedBranchRetrieve::retrieve(block, unique_block);
}

void
DoNothingBranchRetrieve::retrieve(TensorBranch* branch)
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
    TensorBranch* branch
)
{
    TensorBlock* block = branch->get_parent_block();
    //this is the location of the old branch in memory
    TensorBranch* old_branch = block->get_obsolete_branch();
    TensorBranch* new_branch = block->get_branch();

    if (!old_branch)
    {
        raise(SanityCheckError, "no old branch to realign");
    }
    new_branch->realign_memory_pool(old_branch, new_branch);
}

void
NewBranchRetrieve::retrieve(TensorBranch* branch)
{
    /** This should have been reinitialized by memory pool reset */

    //TensorBlock* block = branch->get_parent_block();
    //block->init_branch();
}

void
SortedBranchRetrieve::retrieve(
    TensorBranch* branch
)
{
    TensorBlock* block = branch->get_parent_block();
    TensorBlock* unique_block = block->get_symmetry_unique_block();

    retrieve(block, unique_block);
}

void
SortedBranchRetrieve::retrieve(
    TensorBlock* block,
    TensorBlock* unique_block
)
{
    TensorBranch* branch = block->get_branch();
    if (!unique_block)
    {
        ResetMempool reset;
        reset.retrieve(branch);
        return;
    }

    MemoryPool* mempool = branch->get_metadata_mempool();



    unique_block->set_read_mode();
    unique_block->retrieve();
    mempool->memcpy(unique_block->get_metadata_mempool());
    //reset the parent - gets messed up by memcpy
    branch->set_parent(block);
    branch->set_metadata_mempool(mempool);

    Sort* sort = new Sort(block->get_resort_permutation());
    unique_block->get_branch()
        ->sort_metadata_into(sort, branch);

    //now that the sort of metadata has taken place,
    //realign the memory pool
    branch->realign_memory_pool(
        unique_block->get_branch(),
        branch
    );

    unique_block->release();
    delete sort;
}


void
DoNothingDataControllerRetrieve::retrieve(
    TensorBranch* branch
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
    TensorBranch* branch
)
{
    TensorDataController* controller = branch->first_data_controller();
    DataStorageNode* node = branch->get_parent_block()->first_data_node(); 
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
    TensorBranch* branch
)
{
    TensorDataController* controller = branch->first_data_controller();
    DataStorageNode* node = branch->get_parent_block()->first_data_node(); 
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
    TensorBranch* branch
)
{
    TensorBlock* block = branch->get_parent_block();
    TensorBlock* unique_block = block->get_symmetry_unique_block();
    if (!unique_block)
        return;

    unique_block->set_read_mode();
    unique_block->retrieve();
    Sort* sort = new Sort(block->get_resort_permutation());
    unique_block->get_branch()->sort_data_into(sort, branch);
    unique_block->release();

    delete sort;
}

void
ReallocateDataControllers::retrieve(
    TensorBranch* branch
)
{
    TensorBlock* block = branch->get_parent_block();
    TensorDataController* controller = branch->first_data_controller();

    block->clear_storage();
    while(controller)
    {
        StorageBlock* storage_block = branch->get_parent_block()
                                ->allocate_data_storage_block();
        controller->set_data(storage_block->data());
        controller = controller->next;
    }
}

void
DoNothingDataControllerInit::retrieve(TensorBranch* branch)
{
}

void
DoNothingBranchFlush::flush(TensorBranch* branch)
{
}

void
ClearBranchFlush::flush(TensorBranch* branch)
{
    TensorBlock* tensor_block = branch->get_parent_block();
    StorageBlock* block = tensor_block->get_metadata_storage_block();
    block->clear(tensor_block);
}

void
SortedAccumulateBranchFlush::flush(TensorBranch* branch)
{   
    TensorBlock* block = branch->get_parent_block();
    if (block->in_destructor())
        return;

    TensorBlock* unique_block = block->get_symmetry_unique_block();
    Permutation* unsort_perm = block->get_resort_permutation()->inverse();
    Sort* sort = new Sort(unsort_perm);
    unique_block->set_accumulate_mode();
    unique_block->retrieve();
    unique_block->accumulate(block, 1.0, sort);
    unique_block->release();
    delete sort;
}

void
CommitBranchFlush::flush(TensorBranch* branch)
{
    TensorBlock* tensor_block = branch->get_parent_block();
    tensor_block->set_obsolete_branch(branch);
    StorageBlock* block = tensor_block->get_metadata_storage_block();
    block->commit();
    block->clear(tensor_block);
}

void
DoNothingBranchRelease::release(TensorBranch* branch)
{
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
        data_block->clear(block);
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
        data_block->clear(block);
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
FlushOnSync::sync(TensorBlock* block)
{
    if (block->in_cache())
    {
        block->lock();
        if (!block->in_cache()) //race condition!
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

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

