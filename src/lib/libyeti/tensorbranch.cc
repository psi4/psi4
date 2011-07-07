#include "tensor.h"
#include "tensorblock.h"
#include "tensorbranch.h"
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
#include "datapolicies.h"
#include "metadatapolicies.h"
#include "branchpolicies.h"
#include "elementop.h"
#include "matrix.h"

#include <libsmartptr/strop.h>

using namespace yeti;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif



TensorBranch::TensorBranch(
    const uli* indexset,
    usi depth,
    TensorBlock* parent
)
    :
    element_type_(TemplateInfo::double_type), //default to double for now
    parent_block_(parent),
    mdnode_(0),
    filler_(0),
    element_size_(sizeof(double)),
    first_controller_(0),
    last_controller_(0),
    mempool_(parent->get_metadata_mempool())
{
    mdnode_ = new (mempool_)
        MetaDataNode(indexset, depth, this);

}

TensorBranch::~TensorBranch()
{
    delete mdnode_;
}

void
TensorBranch::allocate_data_controller()
{
    StorageBlock* block = parent_block_->allocate_data_storage_block();
    TensorDataController* controller
        = new (mempool_)
        TensorDataController(
            this, block
        );
    controller->next = 0;
    //add the newly allocated controller to the list of controllers
    if (first_controller_ == 0)
    {
        first_controller_ = controller;
        last_controller_ = controller;
    }
    else
    {
        last_controller_->next = controller;
        last_controller_ = controller;
    }
}

void
TensorBranch::allocate(DataNode* node)
{
    size_t blocksize = node->nelements() * element_size();
    if (first_controller_ == 0)
        allocate_data_controller();

    bool enough_space = last_controller_->register_node(node, blocksize);

    if (!enough_space)
    {
        TensorDataController* controller = first_controller_;
        while (controller && !enough_space)
        {
            enough_space = controller->register_node(node, blocksize);
            controller = controller->next;
        }
    }

    if (!enough_space) //allocate a new memory block
    {
        allocate_data_controller();
        enough_space = last_controller_->register_node(node, blocksize);
        if (!enough_space)
        {
            std::string str = stream_printf(
                "Not enough space for data block of size %ld\n"
                "nelements=%ld size=%ld\n"
                "A larger data block size is needed!",
                blocksize, node->nelements(), element_size()
            );
            raise(SanityCheckError, str);
        }
    }

    ::memset(node->data(), 0, blocksize);
}

TensorDataController*
TensorBranch::first_data_controller() const
{
    return first_controller_;
}


size_t
TensorBranch::element_size() const
{
    return element_size_;
}

TemplateInfo::type_t
TensorBranch::element_type() const
{
    return element_type_;
}

size_t
TensorBranch::size_data_wasted() 
{
    size_t wasted = 0;
    TensorDataController* controller = first_controller_;
    while (controller)
    {
        wasted += controller->get_remaining();
        controller = controller->next;
    }
    return wasted;
}

TensorIndexDescr*
TensorBranch::get_descr() const
{
    return parent_block_->descr_;
}

TensorElementComputer*
TensorBranch::get_element_computer() const
{
    return filler_;
}

MemoryPool*
TensorBranch::get_metadata_mempool() const
{
    return mempool_;
}

MetaDataNode*
TensorBranch::get_node() const
{
    return mdnode_;
}

TensorBlock*
TensorBranch::get_parent_block() const
{
    return parent_block_;
}

TensorController*
TensorBranch::get_tensor_controller() const
{
    return parent_block_->get_tensor_controller();
}

void
TensorBranch::realign_memory_pool(
    TensorBranch* old_branch,
    TensorBranch* new_branch
)
{
    realign_pointer(mdnode_, old_branch, new_branch);

    if (first_controller_)
    {
        realign_pointer(first_controller_, old_branch, new_branch);
    }
    if (last_controller_)
    {
        realign_pointer(last_controller_, old_branch, new_branch);
    }

    TensorDataController* controller = first_controller_;
    while(controller)
    {
        controller->realign_memory_pool(old_branch, new_branch);
        controller = controller->next;
    }

    mdnode_->realign_memory_pool(old_branch, new_branch);

}

void
TensorBranch::set_element_computer(TensorElementComputer* filler)
{
    filler_ = filler;
    
    element_type_ = filler->element_type(
                        parent_block_->get_indices(), 
                        mdnode_->get_depth()
                    );

}

void
TensorBranch::set_parent(TensorBlock* parent)
{
    parent_block_ = parent;
}

void
TensorBranch::set_metadata_mempool(MemoryPool* mempool)
{
    mempool_ = mempool;
}

void
TensorBranch::set_element_type(TemplateInfo::type_t type)
{
    element_type_ = type;
    data_type_switch(element_type_, size_of, element_size_);
    parent_block_->set_element_size(element_size_);
}

void
TensorBranch::sort_data_into(
    Sort* sort,
    TensorBranch* new_branch
)
{
    Permutation* descr_perm = 0;
    Permutation* block_perm = this->get_parent_block()->get_block_permutation();
    Permutation* unsort_perm = new_branch->get_parent_block()->get_parent_tensor()->get_unsort_permutation();
    Permutation* fetch_perm = new_branch->get_parent_block()->get_fetch_permutation();
    Permutation* resort_perm = new_branch->get_parent_block()->get_resort_permutation();
#if 1
    if (!block_perm->is_inverse(unsort_perm))
    {
        descr_perm = block_perm->product(unsort_perm);
    }
#else
    if (!fetch_perm->is_inverse(resort_perm))
    {
        descr_perm = resort_perm->inverse();
    }
#endif
    mdnode_->sort_data_into(sort, this, new_branch, descr_perm);
}

void
TensorBranch::sort_metadata_into(
    Sort* sort,
    TensorBranch* new_branch
)
{
   mdnode_->sort_metadata_into(sort, this, new_branch);
}
