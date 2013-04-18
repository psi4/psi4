#ifndef yeti_tensoraction_h
#define yeti_tensoraction_h

#include "class.h"
#include "mapimpl.h"
#include "mallocimpl.h"
#include "tensorparser.h"

#include "data.hpp"
#include "index.hpp"
#include "tensor.hpp"
#include "tensorbranch.hpp"
#include "tensoraction.hpp"
#include "tensorblock.hpp"
#include "node.hpp"
#include "permutation.hpp"
#include "filler.hpp"
#include "elementop.hpp"
#include "matrix.hpp"
#include "contraction.hpp"


#include "cache.h"

#include "gigmatrix.h"

#include <list>

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

class BlockRetrieveAction
{
    public:
        virtual void accumulate_to(TensorBlock* block) = 0;

        virtual ~BlockRetrieveAction();

};

class SortBlockRetrieveAction  :
    public BlockRetrieveAction
{
    private:
        Permutation* perm_;

    public:
        SortBlockRetrieveAction(
            Permutation* p
        );

        void accumulate_to(TensorBlock* block);

};

class ElementOpRetrieveAction  :
    public BlockRetrieveAction
{
    private:
        ElementOp* element_op_;

    public:
        ElementOpRetrieveAction(
            ElementOp* op
        );

        void accumulate_to(TensorBlock* block);

};
    
class AccumulateBlockRetrieveAction  :
    public BlockRetrieveAction
{
    private:
        YetiTensorPtr tensor_;

    public:
        AccumulateBlockRetrieveAction(
            const YetiTensorPtr& tensor
        );

        ~AccumulateBlockRetrieveAction();

        void accumulate_to(TensorBlock* block);

};

class ContractionBlockRetrieveAction :
    public BlockRetrieveAction
{
    private:
        YetiContractionPtr yeticxn_;

        Indexer* block_indexer_;

        ContractionPtr cxn_;

    public:
        ContractionBlockRetrieveAction(
            const YetiContractionPtr& cxn,
            Tensor* ptensor
        );

        ~ContractionBlockRetrieveAction();

        void accumulate_to(TensorBlock* block);

};

class TensorRetrieveAction :
    public smartptr::Countable
{
    private:
        Tensor* tensor_;

        std::list<BlockRetrieveAction*> actions_;

    public:
        TensorRetrieveAction(Tensor* tensor);

        ~TensorRetrieveAction();

        void accumulate_to(TensorBlock* block);

        void add(BlockRetrieveAction* action);

};

}

#endif
