#include "tensor.h"
#include "tensorblock.h"
#include "tensorbranch.h"
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

#define PRINT_VALUES 1

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

MetaDataNode::MetaDataNode(
    const uli* indexset,
    usi depth,
    TensorBranch* parent
) :
    depth_(depth),
    parent_branch_(parent),
    nindex_(parent->get_descr()->nindex()),
    mempool_(parent->get_metadata_mempool()),
    nodes_(0)
{
    ::memcpy(indices_, indexset, nindex_ * sizeof(uli));
}

MetaDataNode::~MetaDataNode()
{
    if (nodes_)
        delete nodes_;
}

MetaDataNode::iterator
MetaDataNode::begin() const
{
    return nodes_->begin();
}

MetaDataNode::iterator
MetaDataNode::end() const
{
    return nodes_->end();
}

void
MetaDataNode::realign_memory_pool(
    TensorBranch* old_branch,
    TensorBranch* new_branch
)
{
    realign_pointer(parent_branch_, old_branch, new_branch);
    if (nodes_)
    {
        realign_pointer(nodes_, old_branch, new_branch);
        nodes_->realign_memory_pool(old_branch, new_branch);
    }

    if (depth_ == 1)
        return;

    if (nodes_ == 0)
        return;

    TileMap::iterator it(nodes_->begin());
    TileMap::iterator stop(nodes_->end());
    for ( ; it != stop; ++it)
    {
        MetaDataNode* node = static_cast<MetaDataNode*>(*it);
        if (node)
            node->realign_memory_pool(old_branch, new_branch);
    }
}

template <typename data_t>
void
MetaDataNode::_accumulate(
    DataNode* src,
    DataNode* dst,
    double scale,
    Sort* sort
)
{
    uli nelements = src->nelements();

    data_t* srcdata = reinterpret_cast<data_t*>(src->data());
    data_t* dstdata = reinterpret_cast<data_t*>(dst->data());
    if (sort && src == dst)
    {
        data_t* buffer = reinterpret_cast<data_t*>(sort->data_buffer());
        sort->sort<data_t>(srcdata, buffer);
        //now accumulate values from the buffer
        for (uli i=0; i < nelements; ++i, ++srcdata, ++dstdata)
        {
            (*dstdata) += (*srcdata); //these are already scaled
        }
    }
    else if (sort)
    {
        //directly accumulate
        sort->accumulate<data_t, data_t>(srcdata, dstdata, scale);
    }
    else
    {
        for (uli i=0; i < nelements; ++i, ++srcdata, ++dstdata)
        {
            (*dstdata) += scale * (*srcdata);
        }
    }
}

template <typename data_t>
void
MetaDataNode::_convert(
    data_t*& dataptr,
    data_t* matrixptr,
    usi index,
    const uli* nelements,
    const uli* strides
)
{
    uli n = nelements[index];
    uli stride = strides[index];
    if (index == nindex_ - 1)
    {
        for (uli i=0; i < n; ++i, ++matrixptr, ++dataptr)
        {
            *matrixptr = *dataptr;
        }
        return;
    }
    else
    {
        for (uli i=0; i < n; ++i, matrixptr += stride)
            _convert<data_t>(dataptr, matrixptr, index + 1, nelements, strides);
    }    
}

template <typename data_t>
void
MetaDataNode::_convert(
    DataNode* node,
    Matrix* matrix,
    uli matrix_offset,
    const uli* nelements,
    const uli* strides
)
{
    data_t* dataptr = reinterpret_cast<data_t*>(node->data());
    data_t* matrixptr = const_cast<data_t*>(reinterpret_cast<const data_t*>(matrix->data()));
    matrixptr += matrix_offset;
    _convert<data_t>(dataptr, matrixptr, 0, nelements, strides);
}

template <typename data_t>
void
MetaDataNode::_accumulate(
    data_t*& dataptr,
    data_t* matrixptr,
    usi index,
    const uli* nelements,
    const uli* strides
)
{
    uli n = nelements[index];
    uli stride = strides[index];
    if (index == nindex_ - 1)
    {
        for (uli i=0; i < n; ++i, ++matrixptr, ++dataptr)
        {
            *dataptr += *matrixptr;
        }
        return;
    }
    else
    {
        for (uli i=0; i < n; ++i, matrixptr += stride)
            _accumulate<data_t>(dataptr, matrixptr, index + 1, nelements, strides);
    }    
}

template <typename data_t>
void
MetaDataNode::_accumulate(
    DataNode* node,
    Matrix* matrix,
    uli matrix_offset,
    const uli* nelements,
    const uli* strides
)
{
    data_t* dataptr = reinterpret_cast<data_t*>(node->data());
    data_t* matrixptr = const_cast<data_t*>(reinterpret_cast<const data_t*>(matrix->data()));
    matrixptr += matrix_offset;
    _accumulate<data_t>(dataptr, matrixptr, 0, nelements, strides);
}

template <typename data_t>
void
MetaDataNode::_element_op(
    DataNode* node,
    const uli* index_starts,
    const uli* sizes,
    ElementOp* op
)
{
    op->element_op(
        node->nelements(),
        index_starts,
        sizes,
        reinterpret_cast<data_t*>(node->data())
    );
}

template <typename data_t>
void
MetaDataNode::_fill(
    DataNode* node,
    const uli* indices,
    TensorElementComputer* filler
)
{

    data_t* dataptr = reinterpret_cast<data_t*>(node->data());
    filler->compute(indices, dataptr, node->nelements());

    node->compute_max_log<data_t>();
}


template <typename data_t>
void
MetaDataNode::_print_data_node(
    DataNode* node,
    TensorBranch* branch,
    std::ostream& os
) const
{
    os << Env::indent;
    uli n = node->nelements();
    char* data = node->data();
    os << stream_printf("node: %p data: %p", node, (void*) data);
    const data_t* dataptr = reinterpret_cast<const data_t*>(data);
    for (uli i=0; i < n; ++i, ++dataptr)
    {
        if (i % 6 == 0)
            os << endl << Env::indent;
        os << stream_printf(TypeInfo<data_t>::printf_str, *dataptr)
            << " ";
     }
}

template <typename data_t>
void
MetaDataNode::_equals(
    DataNode* lnode,
    DataNode* rnode,
    bool& equals
)
{
    data_t* lptr = reinterpret_cast<data_t*>(lnode->data());
    data_t* rptr = reinterpret_cast<data_t*>(rnode->data());

    uli n = lnode->nelements();
    for (uli i=0; i < n; ++i, ++lptr, ++rptr)
    {
        data_t l = *lptr;
        data_t r = *rptr;
        if (!TestEquals<data_t>::equals(l, r))
        {
            equals = false;
            return;
        }
    }

    equals = true;
}

template <typename data_t>
void
MetaDataNode::_internal_contraction(
    DataNode* src_node,
    DataNode* dst_node,
    uli nr, uli nc
)
{
    data_t* dst_ptr = reinterpret_cast<data_t*>(dst_node->data());
    data_t* src_ptr = reinterpret_cast<data_t*>(src_node->data());
    data_t* data = dst_ptr;

    for (uli r=0; r < nr; ++r, ++dst_ptr)
    {
        for (uli c=0; c < nc; ++c, ++src_ptr)
        {
            (*dst_ptr) += (*src_ptr);
        }
    }
}

template <typename data_t>
void
MetaDataNode::_sort(
    Sort* sort,
    char* src,
    char* dst
)
{
    sort->sort<data_t>(
        reinterpret_cast<data_t*>(src),
        reinterpret_cast<data_t*>(dst)
    );
}

template <typename data_t>
void
MetaDataNode::_sort(
    DataNode* node,
    Sort* sort
)
{
    data_t* data = reinterpret_cast<data_t*>(node->data());
    data_t* buffer = reinterpret_cast<data_t*>(sort->data_buffer());
    sort->sort<data_t>(data, buffer);

    ::memcpy(data, buffer, node->nelements() * sizeof(data_t));
}

void
MetaDataNode::init_tile_map()
{
    if (!nodes_)
    {
        nodes_ = new (mempool_) NodeMap<TileNode>(
                        mempool_,
                        parent_branch_->get_descr(),
                        indices_,
                        depth_
                     );
    }
}

void
MetaDataNode::internal_contraction(
    MetaDataNode* dst_mdnode,
    MatrixConfiguration* config
)
{
    config->configure_block(
        nodes_->sizes(),
        depth_
    );

    uli nr = config->nrows();
    uli nc = config->ncols();
    uli dst_indices[NINDEX];
    uli src_indices[NINDEX];
    uli nelements[NINDEX];

    TileNode** it_src = nodes_->begin();
    TileMap* dst_map = dst_mdnode->get_node_map();
    TileMap* src_map = this->get_node_map();

    TemplateInfo::type_t dst_type = dst_mdnode->get_parent_branch()->element_type();
    TemplateInfo::type_t src_type = parent_branch_->element_type();

    if (dst_type != src_type)
    {
        cerr << "Internal contraction has different data types" << endl;
        abort();
    }

    TensorIndexDescr* descr = parent_branch_->get_descr();

    uli idx = 0;
    for (uli r=0; r < nr; ++r)
    {
        TileNode* dst_node = 0;
        if (dst_map)
            dst_node = dst_map->get(r);

        for (uli c=0; c < nc; ++c, ++it_src, ++idx)
        {
            TileNode* src_node = *it_src;
            if (!src_node)
                continue;

            if (!dst_node)
            {
                if (!dst_map)
                {
                    dst_mdnode->init_tile_map();
                    dst_map = dst_mdnode->get_node_map();
                }

                dst_map->indices(r, dst_indices);

                if (depth_ == 1)
                {
                    DataNode* dnode = new (mempool_)
                        DataNode(dst_indices, dst_mdnode->get_parent_branch());
                    dst_mdnode->get_parent_branch()->allocate(dnode);
                    dst_node = dnode;
                }
                else
                {
                    dst_node = new (mempool_)
                        MetaDataNode(dst_indices, depth_ - 1,
                            dst_mdnode->get_parent_branch());
                }

                dst_map->set(r, dst_node);
            }

            if (depth_ == 1)
            {
                DataNode* node_dst = static_cast<DataNode*>(dst_node);
                DataNode* node_src = static_cast<DataNode*>(src_node);
                src_map->indices(idx, src_indices);
                descr->get_nelements(src_indices, nelements);
                config->configure_block(nelements, 0);
                uli ndata_rows = config->nrows();
                uli ndata_cols = config->ncols();
                data_type_switch(
                    dst_type,
                    this->_internal_contraction,
                    node_src,
                    node_dst,
                    ndata_rows,
                    ndata_cols
                );

            }
            else
            {
                MetaDataNode* node_dst = static_cast<MetaDataNode*>(dst_node);
                MetaDataNode* node_src = static_cast<MetaDataNode*>(src_node);
                node_src->internal_contraction(node_src, config);
            }
        }
    }

    config->reset_contraction_depth(depth_ + 1);
}

void
MetaDataNode::init_subnodes()
{
    if (nodes_) //nothing to do
        return;

    init_tile_map();

    NodeMap<TileNode>::iterator it(nodes_->begin());
    NodeMap<TileNode>::iterator stop(nodes_->end());
    uli idx = 0;
    uli indexset[NINDEX];
    TileNode* node = 0;
    usi next_depth = depth_ - 1;
    for ( ; it != stop; ++it, ++idx)
    {
        nodes_->indices(idx, indexset);
        if      (depth_ == 1)
        {
            node = new (mempool_) DataNode(indexset, parent_branch_);
        }
        else
        {
            node = new (mempool_) MetaDataNode(indexset, next_depth, parent_branch_);
        }
        nodes_->set(idx, node);
    }
}

void
MetaDataNode::accumulate(
    MetaDataNode* src_parent,
    double scale,
    Sort* sort
)
{
    //nothing here to accumulate
    if (src_parent->is_empty())
        return;

    TensorController* src_controller = src_parent->get_parent_branch()
                                        ->get_tensor_controller();
    TensorController* dst_controller = parent_branch_->get_tensor_controller();

    src_controller->retrieve(src_parent);
    dst_controller->retrieve(this);

    if (nodes_ == 0) //needs to be allocated
    {
        init_tile_map();
    }

    NodeMap<TileNode>::iterator itdst(nodes_->begin());
    NodeMap<TileNode>::iterator itsrc(src_parent->begin());
    NodeMap<TileNode>::iterator stop(src_parent->end());

    uli indexset[NINDEX];
    uli nelements[NINDEX];
    uli src_indices[NINDEX];

    TileMap* src_map = src_parent->get_node_map();

    if (sort)
    {
        src_parent->get_nelements(nelements);
        sort->configure(nelements);
        TileNode** buffer
            = reinterpret_cast<TileNode**>(sort->metadata_buffer(depth_ - 1));
        sort->sort_noscale<TileNode*>(itsrc, buffer);
        stop = buffer + (stop - itsrc);
        itsrc = buffer;
    }

    bool need_recompute_dst = dst_controller->need_recompute_data();
    if (need_recompute_dst)
    {
        cerr << "Cannot accumulate to a recomputed tensor!" << endl;
        abort();
    }
    bool need_recompute_src = src_controller->need_recompute_data();


    TemplateInfo::type_t elem_type = parent_branch_->element_type();
    usi next_depth = depth_ - 1;
    uli idx = 0;
    Permutation* perm = 0;
    if (sort)
        perm = sort->get_permutation();
    for ( ; itsrc != stop; ++itsrc, ++itdst, ++idx)
    {
        TileNode* srcnode = *itsrc;
        if (!srcnode)
            continue;

        TileNode* dstnode = *itdst;
        nodes_->indices(idx, indexset);

        if (depth_ == 1)
        {
            DataNode* dst_dnode = static_cast<DataNode*>(dstnode);
            DataNode* src_dnode = static_cast<DataNode*>(srcnode);
            if (!dst_dnode)
            {
                dst_dnode = new (mempool_)
                    DataNode(indexset, parent_branch_);
                parent_branch_->allocate(dst_dnode);
                nodes_->set(idx, dst_dnode);
            }

            if (perm)
            {
                IndexableMap<TileNode>::permuted_indices(
                    idx, nindex_, src_indices,
                    nodes_->cumulative_sizes(),
                    src_map->offsets(),
                    perm
                );
            }
            else
            {
                src_map->indices(idx, src_indices);
            }

            if (need_recompute_src && src_dnode->data() == 0)
                src_controller->retrieve(src_dnode, src_parent, src_indices);

            if (sort)
            {
                src_parent->get_descr()->get_nelements(src_indices, nelements);
                sort->configure(nelements);
            }

#if YETI_SANITY_CHECK
            if (src_dnode->nelements() != dst_dnode->nelements())
            {
                string srcname = src_parent->get_parent_branch()->get_parent_block()
                                    ->get_parent_tensor()->get_name();
                string dstname = parent_branch_->get_parent_block()
                                    ->get_parent_tensor()->get_name();
                cerr << "misaligned data blocks in accumulate on tensor "
                     << srcname
                     << " from tensor "
                     << dstname
                     << endl;
                cerr << indexstr(nindex_, indexset) 
                     << " on tensor " << dstname
                     << " has "
                     << dst_dnode->nelements() << " elements" << endl;
                cerr << indexstr(nindex_, src_indices) 
                     << " on tensor " << srcname
                     << " has "
                     << src_dnode->nelements() << " elements" << endl;
                abort();
            }

            if (sort)
            {
                if (src_dnode->nelements() != sort->nelements())
                {
                    cerr << idx << " indexed as " << indexstr(nindex_, src_indices) << endl;
                    cerr << "sort not configured properly for source node" << endl;
                    cerr << "sort has " << sort->nelements() << " elements" << endl;
                    cerr << "node has " << src_dnode->nelements() << " elements" << endl;
                    cerr << "lengths: " << indexstr(sort->nindex(), sort->lengths()) << endl;
                    cerr << "nstrides: " << indexstr(sort->nindex(), sort->nstrides()) << endl;
                    cerr << "src sizes: " << indexstr(nindex_, src_map->sizes()) << endl;
                    cerr << "src offsets: " << indexstr(nindex_, src_map->offsets()) << endl;
                    cerr << "dst sizes: " << indexstr(nindex_, nodes_->sizes()) << endl;
                    cerr << "dst offsets: " << indexstr(nindex_, nodes_->offsets()) << endl;
                    sort->get_permutation()->print(cerr); cerr << endl;
                    cerr << indexstr(nindex_, indexset) << " is accumulating " << indexstr(nindex_, src_indices) << endl;
                    abort();
                }
                else if ((src_dnode->nelements() - 1) != sort->total_stride())
                {
                    cerr << "sort not configured properly for source node" << endl;
                    cerr << "sort has total stride " << sort->total_stride() << endl;
                    cerr << "node has " << src_dnode->nelements() << " elements" << endl;
                    cerr << "lengths: " << indexstr(sort->nindex(), sort->lengths()) << endl;
                    cerr << "nstrides: " << indexstr(sort->nindex(), sort->nstrides()) << endl;
                    sort->get_permutation()->print(cerr); cerr << endl;
                    abort();
                }
            }
#endif


            data_type_switch(
                elem_type,
                this->_accumulate,
                src_dnode,
                dst_dnode,
                scale,
                sort
            );

             dstnode = dst_dnode;
        }
        else
        {
            MetaDataNode* dst_mdnode = static_cast<MetaDataNode*>(dstnode);
            MetaDataNode* src_mdnode = static_cast<MetaDataNode*>(srcnode);
            if (!dst_mdnode)
            {
                nodes_->indices(idx, indexset);
                dst_mdnode = new (mempool_)
                    MetaDataNode(indexset, next_depth, parent_branch_);
                nodes_->set(idx, dst_mdnode);
            }

            dst_mdnode->accumulate(src_mdnode, scale, sort);
            dstnode = dst_mdnode;
        }
    }
}

bool
MetaDataNode::equals(const void** data)
{
    parent_branch_->get_tensor_controller()->retrieve(this);
    if (depth_ == 1)
        return data_equals(data);
    else
        return metadata_equals(data);
}

bool
MetaDataNode::equals(MetaDataNode* mdnode)
{
    TensorController* controller_me = parent_branch_->get_tensor_controller();
    controller_me->retrieve(this);

    TensorController* controller_other = mdnode->get_parent_branch()->get_tensor_controller();
    controller_other->retrieve(mdnode);

    TileMap* nodes_other = mdnode->get_node_map();

    bool i_have_null_nodes = (nodes_ == 0) || maxlog_ < YetiRuntime::nonnull_cutoff;
    bool other_has_null_nodes = (nodes_ == 0)
                    || mdnode->get_max_log() < YetiRuntime::nonnull_cutoff;

    if (i_have_null_nodes != other_has_null_nodes)
        return false;
    else if (i_have_null_nodes)
        return true;


    NodeMap<TileNode>::iterator it_me(nodes_->begin());
    NodeMap<TileNode>::iterator stop(nodes_->end());
    NodeMap<TileNode>::iterator it_other(nodes_other->begin());

    size_t element_size = parent_branch_->element_size();

    bool need_recompute_me = controller_me->need_recompute_data();
    bool need_recompute_other = controller_other->need_recompute_data();
    uli indexset[NINDEX];


    TensorIndexDescr* descr = parent_branch_->get_descr();

    uli idx = 0;
    for ( ; it_me != stop; ++it_me, ++it_other, ++idx)
    {
        TileNode* node_me = *it_me;
        TileNode* node_other = *it_other;

        bool null_me = node_me ? node_me->get_max_log() < YetiRuntime::nonnull_cutoff : true;
        bool null_other = node_other ?  node_other->get_max_log() < YetiRuntime::nonnull_cutoff : true;

        if (null_me && !null_other)
            return false;
        else if (null_other && !null_me)
            return false;
        else if (null_me && null_other)
            continue;

        if (depth_ == 1)
        {
            DataNode* dnode_me = static_cast<DataNode*>(node_me);
            if (need_recompute_me && dnode_me->data() == 0)
            {
                nodes_->indices(idx, indexset);
                controller_me->retrieve(dnode_me, this, indexset);
            }

            DataNode* dnode_other = static_cast<DataNode*>(node_other);
            if (need_recompute_other && dnode_other->data() == 0)
            {
                nodes_other->indices(idx, indexset);
                controller_other->retrieve(dnode_other, mdnode, indexset);
            }

            size_t blocksize = dnode_me->nelements() * element_size;

            bool equals = false;
            data_type_switch(
                get_parent_branch()->element_type(),
                this->_equals,
                dnode_me,
                dnode_other,
                equals
            );

            //int check = ::memcmp(dnode_me->data(), dnode_other->data(), blocksize);
            if (!equals) //not equal
            {
                //return false;

                uli indexset[NINDEX];
                nodes_->indices(idx, indexset);
                cout << ClassOutput<const uli*>::str(
                    descr->nindex(),
                    indexset
                ) << endl;
                data_type_switch(
                    get_parent_branch()->element_type(),
                    this->_print_data_node,
                    dnode_me,
                    parent_branch_,
                    cout
                );
                cout << endl;

                mdnode->get_node_map()->indices(idx, indexset);
                cout << ClassOutput<const uli*>::str(
                    descr->nindex(),
                    indexset
                ) << endl;
                data_type_switch(
                    get_parent_branch()->element_type(),
                    this->_print_data_node,
                    dnode_other,
                    mdnode->get_parent_branch(),
                    cout
                );


                return false;
            }
            else
            {
            }
        }
        else
        {
            MetaDataNode* mdnode_me = static_cast<MetaDataNode*>(node_me);
            MetaDataNode* mdnode_other = static_cast<MetaDataNode*>(node_other);
            bool check = mdnode_me->equals(mdnode_other);
            if (!check)
                return false;
        }
    }
    return true;
}


bool
MetaDataNode::data_equals(const void** data)
{
    TensorController* controller = parent_branch_->get_tensor_controller();
    size_t element_size = parent_branch_->element_size();
    NodeMap<TileNode>::iterator it(nodes_->begin());
    NodeMap<TileNode>::iterator stop(nodes_->end());
    uli idx = 0;
    const char* dataptr = reinterpret_cast<const char*>(*data);
    bool need_recompute = controller->need_recompute_data();
    uli indexset[NINDEX];
    for ( ; it != stop; ++it, ++idx)
    {
        DataNode* node = static_cast<DataNode*>(*it);
        if (!node || node->get_max_log() < YetiRuntime::nonnull_cutoff)
            continue;


        if (need_recompute && node->data() == 0)
        {
            nodes_->indices(idx, indexset);
            controller->retrieve(node, this, indexset);
        }

        size_t blocksize = node->nelements() * element_size;
        int check = ::memcmp(dataptr, node->data(), blocksize);
        if (check) //not equal
            return false;

        dataptr += node->nelements() * element_size;
    }

    //update the data ptr
    *data = dataptr;

    return true;
}

void
MetaDataNode::convert(
    Matrix* matrix,
    MatrixConfiguration* config,
    TensorIndexDescr* descr
)
{
    if (depth_ == 1)
        convert_data(matrix, config, descr);
    else
        convert_metadata(matrix, config, descr);
}

void
MetaDataNode::accumulate(
    Matrix* matrix,
    MatrixConfiguration* config,
    TensorIndexDescr* descr
)
{
    if (depth_ == 1)
        accumulate_data(matrix, config, descr);
    else
        accumulate_metadata(matrix, config, descr);
}

void
MetaDataNode::convert_data(
    Matrix* matrix,
    MatrixConfiguration* config,
    TensorIndexDescr* descr
)
{
    uli total_sizes[NINDEX];
    uli strides[NINDEX];
    uli nelements[NINDEX];
    uli index_starts[NINDEX];
    uli indices[NINDEX];

    for (usi i=0; i < nindex_; ++i)
    {
        total_sizes[i] = descr->get(i)->nelements_data();
    }

    uli stride = 1;
    for (int i=nindex_-1; i >=0; --i)
    {
        strides[i] = stride;
        stride *= total_sizes[i]; 
    }

    NodeMap<TileNode>::iterator it(nodes_->begin());
    NodeMap<TileNode>::iterator stop(nodes_->end());
    uli idx = 0;
    TemplateInfo::type_t elem_type = parent_branch_->element_type();
    for ( ; it != stop; ++it, ++idx)
    {
        DataNode* node = static_cast<DataNode*>(*it);
        if (!node)
            continue;

        nodes_->indices(idx, indices);
        descr->get_index_increments(indices, index_starts);
        uli matrix_offset = 0;
        for (usi i=0; i < nindex_; ++i)
            matrix_offset += index_starts[i] * strides[i];
        descr->get_nelements(indices, nelements);


        data_type_switch(
            elem_type,
           this->_convert,
           node, matrix, matrix_offset,
           nelements, strides
        );
    }
}

void
MetaDataNode::accumulate_data(
    Matrix* matrix,
    MatrixConfiguration* config,
    TensorIndexDescr* descr
)
{
    if (nodes_ == 0) //needs to be allocated
        init_tile_map();

    uli total_sizes[NINDEX];
    uli strides[NINDEX];
    uli nelements[NINDEX];
    uli index_starts[NINDEX];
    uli indices[NINDEX];

    for (usi i=0; i < nindex_; ++i)
    {
        total_sizes[i] = descr->get(i)->nelements_data();
    }

    uli stride = 1;
    for (int i=nindex_-1; i >=0; --i)
    {
        strides[i] = stride;
        stride *= total_sizes[i]; 
    }

    NodeMap<TileNode>::iterator it(nodes_->begin());
    NodeMap<TileNode>::iterator stop(nodes_->end());
    uli idx = 0;
    TemplateInfo::type_t elem_type = parent_branch_->element_type();
    for ( ; it != stop; ++it, ++idx)
    {
        DataNode* node = static_cast<DataNode*>(*it);
        nodes_->indices(idx, indices);
        if (!node)
        {
            node = new (mempool_) DataNode(indices, this->parent_branch_);
            parent_branch_->allocate(node);
            nodes_->set(idx, node);
        }

        descr->get_index_increments(indices, index_starts);
        descr->get_nelements(indices, nelements);
        uli matrix_offset = 0;
        for (usi i=0; i < nindex_; ++i)
            matrix_offset += index_starts[i] * strides[i];

        data_type_switch(
            elem_type,
            this->_accumulate,
            node, matrix, matrix_offset, 
            nelements, strides 
        );
    }
}

void
MetaDataNode::convert_metadata(
    Matrix* matrix,
    MatrixConfiguration* config,
    TensorIndexDescr* descr
)
{
    NodeMap<TileNode>::iterator it(nodes_->begin());
    NodeMap<TileNode>::iterator stop(nodes_->end());
    uli idx = 0;
    TemplateInfo::type_t elem_type = parent_branch_->element_type();
    for ( ; it != stop; ++it, ++idx)
    {
        MetaDataNode* node = static_cast<MetaDataNode*>(*it);
        if (!node)
            continue;

        node->convert(matrix, config, descr);
    }
}

void
MetaDataNode::accumulate_metadata(
    Matrix* matrix,
    MatrixConfiguration* config,
    TensorIndexDescr* descr
)
{
    if (nodes_ == 0) //needs to be allocated
        init_tile_map();

    NodeMap<TileNode>::iterator it(nodes_->begin());
    NodeMap<TileNode>::iterator stop(nodes_->end());
    uli idx = 0;
    uli indexset[NINDEX];
    usi next_depth = depth_ - 1;
    for ( ; it != stop; ++it, ++idx)
    {
        MetaDataNode* node = static_cast<MetaDataNode*>(*it);
        nodes_->indices(idx, indexset);
        if (!node)
        {
            node = new (mempool_) MetaDataNode(indexset, next_depth, parent_branch_);
            nodes_->set(idx, node);
        }
        node->accumulate(matrix, config, descr);
    }
}

bool
MetaDataNode::metadata_equals(const void** data)
{
    NodeMap<TileNode>::iterator it(nodes_->begin());
    NodeMap<TileNode>::iterator stop(nodes_->end());
    uli idx = 0;
    for ( ; it != stop; ++it, ++idx)
    {
        MetaDataNode* node = static_cast<MetaDataNode*>(*it);
        if (!node || node->get_max_log() < YetiRuntime::nonnull_cutoff)
            continue;

        bool check = node->equals(data);
        if (!check)
            return false;
    }

    return true;
}

void
MetaDataNode::element_op(ElementOp* op)
{
    TensorController* controller = parent_branch_->get_tensor_controller();
    controller->retrieve(this);
    if (nodes_ == 0)
        return;

    NodeMap<TileNode>::iterator it(nodes_->begin());
    NodeMap<TileNode>::iterator stop(nodes_->end());

    uli idx = 0;
    uli sizes[NINDEX];
    uli index_starts[NINDEX];
    uli indexset[NINDEX];
    TensorIndexDescr* descr = parent_branch_->get_descr();
    TemplateInfo::type_t elem_type = parent_branch_->element_type();
    for ( ; it != stop; ++it, ++idx)
    {
        TileNode* node = *it;
        if (!node)
        {
            nodes_->indices(idx, indexset);
            continue;
        }

        if (depth_ == 1)
        {
            DataNode* dnode = static_cast<DataNode*>(node);
            nodes_->indices(idx, indexset);
            descr->get_nelements(indexset, sizes);
            descr->get_index_starts(indexset, index_starts);
            data_type_switch(
                elem_type,
                this->_element_op,
                dnode,
                index_starts,
                sizes,
                op
            );
        }
        else
        {
            MetaDataNode* mdnode = static_cast<MetaDataNode*>(node);
            mdnode->element_op(op);
        }
    }
}

void
MetaDataNode::sort(
    Sort* sort
)
{
    if (nodes_ == 0)
    {
        uli* tmp = reinterpret_cast<uli*>(sort->metadata_buffer(0));
        sort->get_permutation()->permute<uli>(indices_, tmp);
        ::memcpy(indices_, tmp, nindex_ * sizeof(uli));
        return;
    }

    TensorController* controller = parent_branch_->get_tensor_controller();
    controller->retrieve(this);


    NodeMap<TileNode>::iterator it(nodes_->begin());
    NodeMap<TileNode>::iterator stop(nodes_->end());

    uli indexset[NINDEX];
    uli nelements[NINDEX];

    TensorIndexDescr* descr = parent_branch_->get_descr();

    bool need_recompute = controller->need_recompute_data();
    if (need_recompute)
    {
        cerr << "Cannot sort a recomputed tensor block!" << endl;
        abort();
    }

    TemplateInfo::type_t elem_type = parent_branch_->element_type();
    uli idx = 0;
    for ( ; it != stop; ++it, ++idx)
    {
        TileNode* node = *it;
        if (!node)
            continue;

        nodes_->indices(idx, indexset);

        if (depth_ == 1)
        {
            DataNode* dnode = static_cast<DataNode*>(node);

            descr->get_nelements(indexset, nelements);
            sort->configure(nelements);

            data_type_switch(
                elem_type,
                this->_sort,
                dnode,
                sort
            );
        }
        else
        {
            MetaDataNode* mdnode = static_cast<MetaDataNode*>(node);
            mdnode->sort(sort);

        }
    }

    nodes_->sort(sort);
    uli* tmp = reinterpret_cast<uli*>(sort->metadata_buffer(0));
    sort->get_permutation()->permute<uli>(indices_, tmp);
    ::memcpy(indices_, tmp, nindex_ * sizeof(uli));
}

void
MetaDataNode::sort_data_into(
    Sort* sort,
    TensorBranch* old_branch,
    TensorBranch* new_branch
)
{
    if (nodes_ == 0)
        return;

    uli* dstindices = get_realigned_pointer(indices_, old_branch, new_branch);
    sort->get_permutation()->permute(indices_, dstindices);

    TileMap::iterator it(nodes_->begin());
    TileMap::iterator stop(nodes_->end());

    TensorIndexDescr* descr = parent_branch_->get_descr();

    size_t element_size = old_branch->element_size();
    uli src_indices[NINDEX];
    TileMap* src_map = this->get_node_map();
    if (depth_ == 1)
    {
        uli nelements[NINDEX];
        TemplateInfo::type_t elemtype = parent_branch_->element_type();
        uli idx = 0;
        for ( ; it != stop; ++it, ++idx)
        {
            DataNode* srcnode = static_cast<DataNode*>(*it);
            if (!srcnode)
                continue;

            src_map->indices(idx, src_indices);
            descr->get_nelements(src_indices, nelements);
            sort->configure(nelements);
            DataNode* dstnode = get_realigned_pointer(
                                  srcnode, old_branch, new_branch
                                );

#if YETI_SANITY_CHECK
            if (srcnode->nelements() != sort->nelements())
            {
                cerr << "sort not configured properly for source node "
                     << indexstr(nindex_, src_indices)
                     << " on tensor "
                  << parent_branch_->get_parent_block()->get_parent_tensor()->get_name()
                     << endl;
                cerr << "node " << indexstr(nindex_, src_indices) << endl;
                cerr << "sort permutation: ";
                sort->get_permutation()->print(cerr); cerr << endl;
                cerr << "fetch permutation: ";
                new_branch->get_parent_block()->get_fetch_permutation()->print(cerr);
                cerr << endl;

                cerr << "block permutation: ";
                //this->get_parent_branch()
                new_branch->get_parent_block()->get_block_permutation()->print(cerr);
                cerr << endl;

                cerr << "sort has " << sort->nelements() << " elements" << endl;
                cerr << "node has " << srcnode->nelements() << " elements" << endl;
                cerr << "parent node " << indexstr(nindex_, indices_) << endl;
                cerr << "map offsets: " << indexstr(nindex_, src_map->offsets()) << endl;
                cerr << "map sizes: " << indexstr(nindex_, src_map->sizes()) << endl;
                cerr << "sizes: " << indexstr(nindex_, nelements) << endl;
                cerr << "lengths: " << indexstr(sort->nindex(), sort->lengths()) << endl;
                cerr << "nstrides: " << indexstr(sort->nindex(), sort->nstrides()) << endl;
                descr->print(cerr); cerr << endl;

                get_parent_branch()->get_parent_block()->print(cerr);
                cerr << endl;
                abort();
            }
            else if ((srcnode->nelements() - 1) != sort->total_stride())
            {
                cerr << "sort not configured properly for source node" << endl;
                cerr << indexstr(nindex_, src_indices) << endl;
                cerr << "sort has total stride " << sort->total_stride() << endl;
                cerr << "node has " << srcnode->nelements() << " elements" << endl;
                cerr << "lengths: " << indexstr(sort->nindex(), sort->lengths()) << endl;
                cerr << "nstrides: " << indexstr(sort->nindex(), sort->nstrides()) << endl;
                abort();
            }
#endif

            data_type_switch(
                elemtype,
                this->_sort,
                sort,
                srcnode->data(),
                dstnode->data()
            );


        }
    }
    else
    {
        for ( ; it != stop; ++it)
        {
            MetaDataNode* node = static_cast<MetaDataNode*>(*it);
            if (!node)
                continue;
            node->sort_data_into(sort, old_branch, new_branch);
        }
    }
}

void
MetaDataNode::sort_metadata_into(
    Sort* sort,
    TensorBranch* old_branch,
    TensorBranch* new_branch
)
{
    Permutation* p = sort->get_permutation();
    uli* sort_array
        = get_realigned_pointer<uli>(indices_, old_branch, new_branch);
    p->permute(indices_, sort_array);
    if (nodes_ == 0)
        return;

    //sort the nodes themselves
    TileNode** dst = get_realigned_pointer<TileNode*>(
        nodes_->begin(), old_branch, new_branch
    );
    sort->configure(nodes_->sizes());
    sort->sort_noscale<TileNode*>(nodes_->begin(), dst);

    //now sort the array data
    sort_array = get_realigned_pointer<uli>(
        nodes_->sizes(), old_branch, new_branch
    );
    p->permute(nodes_->sizes(), sort_array);

    sort_array = get_realigned_pointer<uli>(
        nodes_->offsets(), old_branch, new_branch
    );
    p->permute(nodes_->offsets(), sort_array);

    TileMap* newmap = get_realigned_pointer<TileMap>(
        nodes_, old_branch, new_branch
    );

    newmap->init_sizes();

    if (depth_ == 1)
        return;

    TileMap::iterator it(nodes_->begin());
    TileMap::iterator stop(nodes_->end());
    for ( ; it != stop; ++it)
    {
        MetaDataNode* node = static_cast<MetaDataNode*>(*it);
        if (!node)
            continue;
        node->sort_metadata_into(sort, old_branch, new_branch);
    }
}

void
MetaDataNode::fill(TensorElementComputer* filler)
{
    usi next_depth = depth_ - 1;

    uli indexset[NINDEX];

    if (nodes_ == 0) //needs to be allocated
        init_tile_map();

    TensorValueEstimater* estimater = filler->get_estimater(depth_);

    NodeMap<TileNode>::iterator it(nodes_->begin());
    NodeMap<TileNode>::iterator stop(nodes_->end());


    uli idx = 0;
    maxlog_ = LOG_ZERO;
    for ( ; it != stop; ++it, ++idx)
    {
        TileNode* node = *it;
        if (node)
        {
            raise(SanityCheckError, "cannot refill tensor node!");
        }

        nodes_->indices(idx, indexset);

        float node_max_log = estimater ? estimater->max_log(indexset) : LOG_NONZERO;
        if (node_max_log < YetiRuntime::nonnull_cutoff)
            continue;

        if (depth_ == 1) //create data nodes
        {

            DataNode* dnode = new (mempool_) DataNode(indexset, parent_branch_);

            parent_branch_->allocate(dnode);

            data_type_switch(
                parent_branch_->element_type(),
                this->_fill,
                dnode,
                indexset,
                filler
            );

            node = dnode;
        }
        else
        {
            MetaDataNode* mdnode = new (mempool_)
                MetaDataNode(indexset, next_depth, parent_branch_);
            mdnode->fill(filler);
            node =  mdnode;
        }

        update_max_log(node->get_max_log());
        nodes_->set(idx, node);
    }
}


void
MetaDataNode::update()
{
    if (!nodes_)
    {
        set_max_log(LOG_ZERO);
        return;
    }

    set_max_log(LOG_UNITIALIZED);

    NodeMap<TileNode>::iterator it(nodes_->begin());
    NodeMap<TileNode>::iterator stop(nodes_->end());

    uli idx = 0;
    for ( ; it != stop; ++it, ++idx)
    {
        TileNode* node = *it;
        if (!node)
            continue;

        if (depth_ == 1) //create data nodes
        {
            DataNode* dnode = static_cast<DataNode*>(node);
            data_type_switch(
                parent_branch_->element_type(),
                dnode->compute_max_log
            );
        }
        else
        {
            MetaDataNode* mdnode = static_cast<MetaDataNode*>(node);
            mdnode->update();
        }
        this->update_max_log(node->get_max_log());
    }
}

DataNode*
MetaDataNode::get_first_data_node() const
{
    if (!nodes_)
        return 0;

    uli idx = 0;
    if (depth_ == 1)
        return static_cast<DataNode*>(nodes_->get(idx));
    else
        return static_cast<MetaDataNode*>(nodes_->get(idx))->get_first_data_node();

}

void
MetaDataNode::get_nelements(uli* nelements) const
{
    //currently an off-by-1 error I simply can't be bothered to fix
    parent_branch_->get_descr()->get_nelements(depth_ + 1, indices_, nelements);
}

TensorIndexDescr*
MetaDataNode::get_descr() const
{
    return parent_branch_->get_descr();
}

usi
MetaDataNode::get_depth() const
{
    return depth_;
}

TileNode*
MetaDataNode::get_node(uli index) const
{
    return nodes_->get(index);
}

NodeMap<TileNode>*
MetaDataNode::get_node_map() const
{
    return nodes_;
}

TensorBranch*
MetaDataNode::get_parent_branch() const
{
    return parent_branch_;
}

TensorController*
MetaDataNode::get_tensor_controller() const
{
    return parent_branch_->get_tensor_controller();
}

bool
MetaDataNode::is_empty() const
{
    return nodes_ == 0;
}

void
MetaDataNode::print_data(std::ostream& os)
{
    TensorController* controller = parent_branch_->get_tensor_controller();
    controller->retrieve(this);

    if (nodes_ == 0)
        return;

    bool need_recompute = controller->need_recompute_data();

    uli idx = 0;
    uli indexset[NINDEX];
    ++Env::indent;
    NodeMap<TileNode>::iterator it(nodes_->begin());
    NodeMap<TileNode>::iterator stop(nodes_->end());
    for ( ; it != stop; ++it, ++idx)
    {
        DataNode* node = static_cast<DataNode*>(*it);
        //if (!node || node->get_max_log() < YetiRuntime::print_cutoff)
        //    continue;

        nodes_->indices(idx, indexset);

        if (need_recompute && node->data() == 0)
            controller->retrieve(node, this, indexset);

        os << ClassOutput<const uli*>::str(nindex_, indexset)
            << stream_printf(" max=10^(%7.4f)", node->get_max_log())
            << endl;

#if PRINT_VALUES
        data_type_switch(
            parent_branch_->element_type(),
            this->_print_data_node,
            node, parent_branch_, os
        )
        os << endl;
#endif
    }
    --Env::indent;
}

void
MetaDataNode::print_metadata(std::ostream &os)
{
    TensorController* controller = parent_branch_->get_tensor_controller();
    controller->retrieve(this);

    if (nodes_ == 0)
        return;

    uli idx = 0;
    uli indexset[NINDEX];
    ++Env::indent;
    NodeMap<TileNode>::iterator it(nodes_->begin());
    NodeMap<TileNode>::iterator stop(nodes_->end());
    for ( ; it != stop; ++it, ++idx)
    {
        MetaDataNode* node = static_cast<MetaDataNode*>(*it);
        //if (!node || node->get_max_log() < YetiRuntime::print_cutoff)
        //    continue;

        nodes_->indices(idx, indexset);
        os << Env::indent << "MetaDataNode " << (void*) this << " "
            << ClassOutput<const uli*>::str(nindex_, indexset)
             << stream_printf(" max=10^(%7.4f)", node->get_max_log())
            << endl;

        node->print(os);
    }
    --Env::indent;
}

void
MetaDataNode::print(std::ostream& os)
{
    if (depth_ == 1)
        this->print_data(os);
    else
        this->print_metadata(os);
}

void
MetaDataNode::accumulate_subnode(
    uli index,
    DataNode* lnode,
    DataNode* rnode,
    uli nrows,
    uli ncols,
    uli nlink,
    double scale,
    Contraction* cxn
)
{
    uli indexset[NINDEX];
    DataNode* pnode = static_cast<DataNode*>(nodes_->get(index));
    if (!pnode)
    {
        nodes_->indices(index, indexset);
        pnode = new (mempool_) DataNode(indexset, parent_branch_);
        nodes_->set(index, pnode);
        parent_branch_->allocate(pnode);
    }

    cxn->get_engine()->contract(
        lnode, rnode, pnode,
        nrows, ncols, nlink,
        scale * cxn->scale_factor(),
        parent_branch_->element_type()
    );

    pnode->set_max_log(1);
}

void
MetaDataNode::accumulate_subnode(
    uli index,
    MetaDataNode* lnode,
    MetaDataNode* rnode,
    double scale,
    Contraction* cxn
)
{
    uli indexset[NINDEX];
    MetaDataNode* pnode = static_cast<MetaDataNode*>(nodes_->get(index));
    if (!pnode)
    {
        nodes_->indices(index, indexset);
        pnode = new (mempool_) MetaDataNode(indexset, depth_ - 1, parent_branch_);
        nodes_->set(index, pnode);
    }
    pnode->accumulate(
        lnode,
        rnode,
        scale,
        cxn
    );
    pnode->set_max_log(1);
}

void
MetaDataNode::accumulate(
    MetaDataNode* l_mdnode,
    MetaDataNode* r_mdnode,
    double scale,
    Contraction* cxn
)
{
    uli threadnum = YetiRuntime::get_thread_number();
    ContractionConfiguration* cxn_config = cxn->get_configuration(threadnum);

    uli indices[NINDEX];
    uli nelements[NINDEX];

    TensorIndexDescr* ldescr= l_mdnode->get_descr();
    TensorController* lcontroller = l_mdnode->get_tensor_controller();
    lcontroller->retrieve(l_mdnode);
    if (!l_mdnode->get_node_map())
        return;

    TensorIndexDescr* rdescr= r_mdnode->get_descr();
    TensorController* rcontroller = r_mdnode->get_tensor_controller();
    rcontroller->retrieve(r_mdnode);
    if (!r_mdnode->get_node_map())
        return;

    cxn_config->configure_left_block(l_mdnode);
    cxn_config->configure_right_block(r_mdnode);
    //cxn->configure_product_block(this);

    uli nrows = cxn_config->ncxn_rows_left();
    uli nlink = cxn_config->ncxn_cols_left();
    uli ncols = cxn_config->ncxn_cols_right();
#if YETI_SANITY_CHECK
    uli check = cxn_config->ncxn_rows_right();
    if (nlink != check)
    {
        cerr << "Depth = " << depth_ << endl;
        cerr << "Left tensor "
            << l_mdnode->get_parent_branch()->get_parent_block()->get_parent_tensor()->get_name() << endl;
        cerr << "Left sizes = "
             << indexstr(l_mdnode->get_descr()->nindex(), l_mdnode->get_node_map()->sizes()) << endl;
        l_mdnode->get_descr()->print(cerr); cerr << endl;
        cxn_config->get_left_config()->get_index()->print(cerr); cerr << endl;
        for (usi i=0; i < l_mdnode->get_descr()->nindex(); ++i)
        {
            uli idx = l_mdnode->indices_[i];
            cerr << stream_printf("size[%d] = %d",
                                  idx,
                                  l_mdnode->get_descr()->get(i)->nelements(depth_, idx)
                                 ) << endl;
        }
        Tensor* rtensor = r_mdnode->get_parent_branch()->get_parent_block()->get_parent_tensor();
        cerr << "Right tensor "
            << rtensor->get_name() << endl;
        rtensor->get_tensor_grp()->print(cerr); cerr << endl;

        cerr << "Right sizes = "
             << indexstr(r_mdnode->get_descr()->nindex(), r_mdnode->get_node_map()->sizes()) << endl;
        r_mdnode->get_descr()->print(cerr); cerr << endl;
        cxn_config->get_right_config()->get_index()->print(cerr); cerr << endl;
        for (usi i=0; i < r_mdnode->get_descr()->nindex(); ++i)
        {
            uli idx = r_mdnode->indices_[i];
            cerr << stream_printf("size[%d] = %d",
                                  idx,
                                  r_mdnode->get_descr()->get(i)->nelements(depth_, idx)
                                 ) << endl;
        }
        std::string str = stream_printf(
            "nlink is %d on left and %d on right.\n"
            "Metadata nodes not configured properly for multiplication!\n"
            "Left node %s\n"
            "Right node %s\n"
            , nlink, check,
            indexstr(l_mdnode->nindex_, l_mdnode->indices_),
            indexstr(r_mdnode->nindex_, r_mdnode->indices_)
        );
        raise(SanityCheckError, str);
    }
#endif

    bool need_recompute_left = lcontroller->need_recompute_data();
    bool need_recompute_right = rcontroller->need_recompute_data();


    uli lidx = 0; uli ridx = 0;
    for (uli row=0; row < nrows; ++row)
    {
        for (uli col=0; col < ncols; ++col)
        {
            for (uli link=0; link < nlink; ++link)
            {
                TileNode* lnode = cxn_config->get_left_node(l_mdnode, row, link, lidx);
                if (!lnode)
                    continue;

                TileNode* rnode = cxn_config->get_right_node(r_mdnode, link, col, ridx);
                if (!rnode)
                    continue;

#if USE_SCREENING
                double maxlog_l = lnode->get_max_log();
                double maxlog_r = rnode->get_max_log();
                if (maxlog_l + maxlog_r < YetiRuntime::matrix_multiply_cutoff)
                   continue; //nothing to accumulate
#endif

                if (!nodes_)
                    init_tile_map();


                uli pidx = row * ncols + col;


                if (depth_ == 1)
                {
                    DataNode* l_subnode = static_cast<DataNode*>(lnode);
                    DataNode* r_subnode = static_cast<DataNode*>(rnode);

                    l_mdnode->get_node_map()->indices(lidx, indices);
                    if (need_recompute_left && l_subnode->data() == 0)
                        lcontroller->retrieve(l_subnode, l_mdnode, indices);

                    ldescr->get_nelements(indices, nelements);
                    uli nrows_data = cxn_config->nrows_left(nelements);
                    uli nlink_data = cxn_config->ncols_left(nelements);

                    r_mdnode->get_node_map()->indices(ridx, indices);
                    if (need_recompute_right && r_subnode->data() == 0)
                        rcontroller->retrieve(r_subnode, r_mdnode, indices);

                    rdescr->get_nelements(indices, nelements);
                    uli ncols_data = cxn_config->ncols_right(nelements);


                    accumulate_subnode(
                        pidx, l_subnode, r_subnode,
                        nrows_data, ncols_data, nlink_data,
                        scale,
                        cxn
                    );
                }
                else
                {

                    MetaDataNode* l_subnode
                        = static_cast<MetaDataNode*>(lnode);

                    MetaDataNode* r_subnode
                        = static_cast<MetaDataNode*>(rnode);

                    accumulate_subnode(pidx, l_subnode, r_subnode, scale, cxn);
                }
            }
        }
    }

    //send the contraction back to the previous depth
    cxn_config->reset_contraction_depth(depth_ + 1);
}

template <typename data_t>
void
DataNode::compute_max_log()
{
    data_t* dataptr = reinterpret_cast<data_t*>(data_);
    float max = 0;
    for (uli i=0; i < nelements_; ++i, ++dataptr)
    {
        float absmax = fabs(*dataptr);
        if (absmax > max)
            max = absmax;
    }
    if (max == 0)
        maxlog_ = LOG_ZERO;
    else
        maxlog_ = log10(max);
}