#include "tensor.h"
#include "node.h"
#include "tensorblock.h"
#include "tensorbranch.h"
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
#include "gigmatrix.h"
#include "contraction.h"
#include "tensoraction.h"

#include <libsmartptr/strop.h>

using namespace yeti;
using namespace std;


BlockRetrieveAction::~BlockRetrieveAction()
{
}

AccumulateBlockRetrieveAction::~AccumulateBlockRetrieveAction()
{
}

ContractionBlockRetrieveAction::~ContractionBlockRetrieveAction()
{
}

TensorRetrieveAction::TensorRetrieveAction(Tensor* tensor)
    : tensor_(tensor)
{
}

TensorRetrieveAction::~TensorRetrieveAction()
{
    
    std::list<BlockRetrieveAction*>::const_iterator it = actions_.begin();
    std::list<BlockRetrieveAction*>::const_iterator stop = actions_.end();
    for ( ; it != stop; ++it)
    {
        delete *it;
    }
}

void
TensorRetrieveAction::add(BlockRetrieveAction* action)
{
    actions_.push_back(action);
}

void
TensorRetrieveAction::accumulate_to(TensorBlock* block)
{
    std::list<BlockRetrieveAction*>::const_iterator it = actions_.begin();
    std::list<BlockRetrieveAction*>::const_iterator stop = actions_.end();
    for ( ; it != stop; ++it)
    {
        (*it)->accumulate_to(block);
    }
}

ContractionBlockRetrieveAction::ContractionBlockRetrieveAction(
    const YetiContractionPtr& cxn,
    Tensor* ptensor
) : 
    yeticxn_(cxn),
    cxn_(0),
    block_indexer_(0)
{
    cxn_ = new Contraction(
        yeticxn_->contraction_scale,
        yeticxn_->ltensor->get_tensor(),
        yeticxn_->rtensor->get_tensor(),
        ptensor,
        yeticxn_->lindex,
        yeticxn_->rindex,
        yeticxn_->ltensor->get_permutation(),
        yeticxn_->rtensor->get_permutation()
    );

    block_indexer_ = new Indexer(ptensor->get_block_descr());

    for (uli i=0; i < YetiRuntime::nthread(); ++i)
    {
        cxn_->get_configuration(i)->configure_left_block(yeticxn_->ltensor->get_tensor());
        cxn_->get_configuration(i)->configure_right_block(yeticxn_->rtensor->get_tensor());
    }
    ContractionConfiguration* config = cxn_->get_configuration(0);
}

void
ContractionBlockRetrieveAction::accumulate_to(TensorBlock* block)
{
    Permutation* l_presort = yeticxn_->left_presort; 
    Tensor* ltensor = yeticxn_->ltensor->get_tensor();
    Permutation* r_presort = yeticxn_->right_presort; 
    Tensor* rtensor = yeticxn_->rtensor->get_tensor();

    if (l_presort || r_presort)
    {
        yeti_throw(SanityCheckError, "contraction block retrieve not valid with sorts");
    }
    
    uli threadnum = YetiRuntime::get_thread_number();
    ContractionConfiguration* cxn_config = cxn_->get_configuration(threadnum);
    uli nrows = cxn_config->ncxn_rows_left();
    uli ncols = cxn_config->ncxn_cols_right();
    uli nlink = cxn_config->ncxn_cols_left();
#if YETI_SANITY_CHECK
    uli check = cxn_config->ncxn_rows_right();
    if (check != nlink)
        yeti_throw(SanityCheckError, "matrices not aligned for multiplication");
#endif
    uli index = block_indexer_->index(block->get_indices());
    uli row = index / ncols;
    uli col = index % ncols;

    for (uli link=0; link < nlink; ++link)
    {
        TensorBlock* lblock = cxn_config->get_left_block(ltensor, row, link);
        if (!lblock || lblock->get_degeneracy() == 0)
        {
            continue;
        }

        TensorBlock* rblock = cxn_config->get_right_block(rtensor, link, col);
        if (!rblock || rblock->get_degeneracy() == 0)
        {
            continue;
        }

        if (!block->is_permutationally_unique())
        {
            yeti_throw(SanityCheckError, "contraction task accumulating to non-unique block");
        }

        lblock->retrieve_read();
        rblock->retrieve_read();
        block->accumulate(lblock, rblock, cxn_.get());
        lblock->release_read();
        rblock->release_read();
    }
}

AccumulateBlockRetrieveAction::AccumulateBlockRetrieveAction(
    const YetiTensorPtr& tensor
) : tensor_(tensor)
{
}

void
AccumulateBlockRetrieveAction::accumulate_to(TensorBlock* block)
{
    Sort* sort = 0;
    double coef = tensor_->get_scale();
    Permutation* p = tensor_->get_permutation();
    Tensor* src_tensor = tensor_->get_tensor();
    TensorBlock* src_block = 0;
    if (p) 
    {
        Permutation* inv = p->inverse();
        src_block = src_tensor->get_block(block->get_indices(), inv);
        if (!src_block)
            return;
        sort = new Sort(p);
    }
    else
    {
        src_block = src_tensor->get_block(block->get_indices());
        if (!src_block)
            return;
    }

    src_block->retrieve_read();
    block->accumulate(src_block, coef, sort);

    src_block->release_read();

    if (p)
    {   
        delete sort;
    }
}


ElementOpRetrieveAction::ElementOpRetrieveAction(
    ElementOp* op
) : element_op_(op)
{
}

void
ElementOpRetrieveAction::accumulate_to(TensorBlock* block)
{
    block->element_op(element_op_);
}

SortBlockRetrieveAction::SortBlockRetrieveAction(
    Permutation* p
) : perm_(p)
{
}

void
SortBlockRetrieveAction::accumulate_to(TensorBlock* block)
{
    Sort* sort = new Sort(perm_);
    block->metadata_sort(sort);
    block->data_sort(sort);
    delete sort;
}
