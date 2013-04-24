#ifndef yeti_node_h
#define yeti_node_h

#include <libsmartptr/serialize.h>
#include "mapimpl.h"

#include "data.hpp"
#include "index.hpp"
#include "tensor.hpp"
#include "tensorbranch.hpp"
#include "tensorblock.hpp"
#include "node.hpp"
#include "permutation.hpp"
#include "tensorparser.hpp"
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


class TileNode
{
    protected:
        float maxlog_;
        
    public:
        TileNode();
        
        float get_max_log() const;
        
        void set_max_log(float maxlog);
        
        void update_max_log(float maxlog);
};

class DataNode :
    public TileNode,
    public MempoolVirtualAddressMalloc
{

    private:
        data_block_vir_addr_unsigned_offset_t offset_;

        uli nelements_;

        char* data_;

    public:
        DataNode(
            const uli* indexset,
            TensorBranch* parent
        );

        ~DataNode();

        template <typename data_t>
        void compute_max_log();

        void set_offset(data_block_vir_addr_unsigned_offset_t offset);

        char* data() const;

        uli nelements() const;

        /**
            Assign a data location to the data node. The
            actual data pointer for the node will be placed
            at block_start + #DataNode::offset_;
        */
        void set_data_block(char* data_block_start);

        DataNode* next_node;
};

class MetaDataNode :
    public TileNode,
    public MempoolVirtualAddressMalloc
{
    private:
        NodeMap<TileNode>* nodes_;

        TensorBranch* parent_branch_;

        usi depth_;

        uli indices_[NINDEX];

        usi nindex_;

        MemoryPool* mempool() const;

        void
        accumulate_subnode(
            uli index,
            DataNode* lnode,
            DataNode* rnode,
            uli nrows,
            uli ncols,
            uli nlink,
            double scale,
            Contraction* cxn
        );

        void
        accumulate_subnode(
            uli index,
            MetaDataNode* lnode,
            MetaDataNode* rnode,
            double scale,
            Contraction* cxn
        );

        template <typename data_t>
        void
        _equals(
            DataNode* lnode,
            DataNode* rnode,
            bool& equals
        );

        template <typename data_t>
        void
        _element_op(
            DataNode* node,
            const uli* index_starts,
            const uli* sizes,
            ElementOp* op
        );

        template <typename data_t>
        void
        _accumulate(
            DataNode* src,
            DataNode* dst,
            double scale,
            Sort* sort
        );

        template <typename data_t>
        void
        _assign(
            DataNode* src,
            DataNode* dst
        );

        template <typename data_t>
        void
        _fill(
            DataNode* node,
            const uli* indices,
            TensorElementComputer* filler
        );
        
        template <typename data_t>
        void
        _update(
            DataNode* node
        );

        template <typename data_t>
        void
        _sort(
            DataNode* node,
            Sort* sort
        );

        template <typename data_t>
        void
        _print_data_node(
            DataNode* node,
            TensorBranch* branch,
            std::ostream& os
        ) const;

        template <typename data_t>
        void
        _sort(
            Sort* sort,
            char* src,
            char* dst
        );

        template <typename data_t>
        void _convert(
            data_t*& dataptr,
            data_t* matrixptr,
            usi index,
            const uli* nelements,
            const uli* strides
        );

        template <typename data_t>
        void _convert(
            DataNode* node,
            Matrix* matrix,
            uli matrix_offset,
            const uli* nelements,
            const uli* strides
        );

        template <typename data_t>
        void _accumulate(
            data_t*& dataptr,
            data_t* matrixptr,
            usi index,
            const uli* nelements,
            const uli* strides
        );

        template <typename data_t>
        void _accumulate(
            DataNode* node,
            Matrix* matrix,
            uli matrix_offset,
            const uli* nelements,
            const uli* strides
        );

        template <typename data_t>
        void
        _internal_contraction(
            DataNode* src_node,
            DataNode* dst_node,
            uli nrows, uli ncols
        );

        bool data_equals(const void** data);
 
        void init_tile_map();
        
        bool metadata_equals(const void** data);

        /**
            The template parameter determines either
            data or metadata nodes
        */
        void print_data(std::ostream& os = std::cout);

        void print_metadata(std::ostream& os = std::cout);

        void convert_data(
            Matrix* matrix,
            MatrixConfiguration* config,
            TensorIndexDescr* descr
        );

        void convert_metadata(
            Matrix* matrix,
            MatrixConfiguration* config,
            TensorIndexDescr* descr
        );

        void accumulate_data(
            Matrix* matrix,
            MatrixConfiguration* config,
            TensorIndexDescr* descr
        );

        void accumulate_metadata(
            Matrix* matrix,
            MatrixConfiguration* config,
            TensorIndexDescr* descr
        );

    public:
        typedef NodeMap<TileNode>::iterator iterator;

        MetaDataNode(
            const uli* indexset,
            usi depth,
            TensorBranch* parent
        );

        ~MetaDataNode();

        void accumulate(
            MetaDataNode* src,
            double scale,
            Sort* sort
        );

        void accumulate(
            MetaDataNode* l_mdnode,
            MetaDataNode* r_mdnode,
            double scale,
            Contraction* cxn
        );

        void convert(
            Matrix* matrix,
            MatrixConfiguration* config,
            TensorIndexDescr* descr
        );

        void accumulate(
            Matrix* matrix,
            MatrixConfiguration* config,
            TensorIndexDescr* descr
        );

        void assign(MetaDataNode* node);

        iterator begin() const;

        iterator end() const;

        void element_op(ElementOp* op);

        bool equals(const void** data);
        
        bool equals(MetaDataNode* mdnode);

        void fill(TensorElementComputer* filler);

        usi get_depth() const;
        
        DataNode* get_first_data_node() const;

        TensorIndexDescr* get_descr() const;

        TileNode* get_node(uli index) const;

        NodeMap<TileNode>* get_node_map() const;

        TensorBranch* get_parent_branch() const;

        const uli* get_indices() const { return indices_; }

        void get_nelements(uli* nelements) const;

        TensorController* get_tensor_controller() const;
        
        void init_subnodes();

        void internal_contraction(
            MetaDataNode* dst_node,
            MatrixConfiguration* config
        );

        bool is_empty() const;

        void load_metadata(long offset);

        void print(std::ostream& os = std::cout);

        void realign_memory_pool(
            TensorBranch* old_branch,
            TensorBranch* new_branch
        );

        void sort(Sort* sort);

        void sort_metadata_into(
            Sort* sort,
            TensorBranch* old_branch,
            TensorBranch* new_branch
        );

        void sort_data_into(
            Sort* sort,
            TensorBranch* old_branch,
            TensorBranch* new_branch
        );

        void update();

};

}

#endif
