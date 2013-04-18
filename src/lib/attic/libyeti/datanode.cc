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

DataNode::DataNode(
    const uli* indexset,
    TensorBranch* parent
) : nelements_(parent->get_descr()->nelements(indexset)),
    data_(0),
    next_node(0)
{
}


DataNode::~DataNode()
{
}

void
DataNode::set_data_block(char* data_block_start)
{
    cout.flush();
    data_ = data_block_start + offset_;
}

void
DataNode::set_offset(data_block_vir_addr_unsigned_offset_t offset)
{
    offset_ = offset;
}


uli
DataNode::nelements() const
{
    return nelements_;
}

char*
DataNode::data() const
{
    return data_;
}

TileNode::TileNode()
    : maxlog_(LOG_UNITIALIZED)
{
}

float
TileNode::get_max_log() const
{
    return maxlog_;
}

void
TileNode::set_max_log(float maxlog)
{
    maxlog_ = maxlog;
}

void
TileNode::update_max_log(float maxlog)
{
    if (maxlog_ == LOG_UNITIALIZED || maxlog > maxlog_)
        maxlog_ = maxlog;
}
