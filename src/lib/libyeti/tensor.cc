#include "tensor.h"
#include "tensoraction.h"
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

#include <libsmartptr/strop.h>

using namespace yeti;
using namespace std;

DECLARE_MALLOC(Tensor);

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

bool
yeti::tensor_less(Tensor* l, Tensor* r)
{
    uli lpriority = l->get_priority();
    uli rpriority = r->get_priority();
    if (lpriority != rpriority)
        return lpriority < rpriority;


    bool l_dist_type = l->is_distributed();
    bool r_dist_type = r->is_distributed();
    if (l_dist_type != r_dist_type)
        return r_dist_type;
    //if true, r is distributed and l is not so l < r
    //if false, l is distributed and r is not so r > l

    //in core is always faster than recomputing which is probably faster than going to disk
    Tensor::tensor_storage_t l_store_type = l->get_storage_type();
    Tensor::tensor_storage_t r_store_type = r->get_storage_type();
    if (l_store_type != r_store_type)
        return l_store_type < r_store_type;

    //everything is the same now, select based on size
    ulli lsize = l->get_totalsize() * l->get_sort_size_weight();
    uli rsize = r->get_totalsize() * r->get_sort_size_weight();
    if (lsize != rsize)
        return lsize < rsize;

    return false; //exactly equal
}

/*********
End thread interface definition
*********/

class ElementOpThread :
    public Thread
{
    private:
        Tensor* tensor_;

        ElementOp* op_;

    public:
        ElementOpThread(
            uli threadnum,
            Tensor* tensor,
            ElementOp* op
        ) :
            Thread(threadnum),
            tensor_(tensor),
            op_(op)
        {
        }

        void run()
        {
            tensor_->element_op(threadnum_, op_);
        }

};

class TensorUpdateThread :
    public Thread
{
    private:
        Tensor* tensor_;

    public:
        TensorUpdateThread(
            uli threadnum,
            Tensor* tensor
        ) :
            Thread(threadnum),
            tensor_(tensor)
        {
        }

        void run()
        {
            tensor_->update(threadnum_);
        }

};

class TensorFillThread :
    public Thread
{
    private:
        Tensor* tensor_;

        ThreadedTensorElementComputerPtr filler_;

    public:
        TensorFillThread(
            uli threadnum,
            Tensor* tensor,
            const ThreadedTensorElementComputerPtr& filler
        ) :
            Thread(threadnum),
            tensor_(tensor),
            filler_(filler)
        {
        }

        void run()
        {
            tensor_->fill(threadnum_, filler_);
        }

};

class AccumulateThread :
    public Thread
{
    private:
        Permutation* perm_;
        
        SortPtr sort_;

        Tensor* src_tensor_;

        Tensor* dst_tensor_;

        double scale_;

    public:
        AccumulateThread(
            uli threadnum,
            Tensor* src,
            Tensor* dst,
            double scale,
            Permutation* p
        ) : Thread(threadnum),
            sort_(0),
            scale_(scale),
            src_tensor_(src),
            dst_tensor_(dst),
            perm_(p)
        {
        }

        ~AccumulateThread()
        {
            sort_ = 0;
        }

        void run()
        {
            sort_ = new Sort(perm_);
            dst_tensor_->accumulate(
                            threadnum_,
                            src_tensor_,
                            scale_,
                            sort_
                         );
        }

};

class SortThread :
    public Thread
{
    private:
        Tensor* tensor_;

        SortPtr sort_;

        Permutation* perm_;

    public:
        SortThread(
            uli threadnum,
            Tensor* tensor,
            Permutation* p
        ) :
            Thread(threadnum),
            tensor_(tensor),
            sort_(0),
            perm_(p)
        {
        }

        ~SortThread(){
            sort_ = 0;
        }

        void run()
        {
            sort_ = new Sort(perm_);
            tensor_->sort(threadnum_, sort_);
        }

        const SortPtr& get_sort() const
        {
            return sort_;
        }

};

#define DEFAULT_TENSOR_CACHE_SIZE 10000000 //10 MB
//#define MAX_TENSOR_METADATA_BLOCK_SIZE 200000
#define MAX_TENSOR_CACHE_BLOCKS 50
//#define MAX_TENSOR_CACHE_BLOCKS NCACHE_ENTRIES_MAX
Tensor::Tensor(
    const std::string& name,
    const TensorIndexDescrPtr& descr,
    const PermutationGroupPtr& pgrp
) :
    name_(name),
    sort_weight_(1),
    tensor_grp_(pgrp),
    original_grp_(pgrp),
    tensor_blocks_(0),
    block_descr_(descr),
    depth_(descr->depth()),
    storage_type_(Tensor::in_core),
    distribution_indices_(0),
    priority_(Tensor::gamma_tensor),
    filler_(0),
    parent_tensor_(0),
    data_cache_(0),
    metadata_cache_(0),
    metadata_block_size_(0),
    data_block_size_(0),
    main_indexer_(new Indexer(descr.get())),
    block_to_node_indexmap_(0),
    distr_indexer_(0),
    full_index_to_distribution_index_(0),
    malloc_number_(0),
    cxn_position_(Contraction::product_tensor)
{
    tensor_blocks_ = new TensorBlockMap(block_descr_.get());


    malloc_number_ = Malloc<Tensor>::get_malloc_number(this);

    for (usi i=0; i < block_descr_->nindex(); ++i)
    {
        IndexRange* range = block_descr_->get(i)->get_top_range();
        index_start_[i] = range->start();
        index_end_[i] = range->start() + range->nelements();
    }

    if (block_descr_->has_symmetry())
        init_symmetry_filter();

    usi overflow = 2;
    size_t max_nblocks_data = block_descr_->max_nblocks_data();
    if (max_nblocks_data < 4)
        max_nblocks_data = 4;


    usi md_depth = 2;
    size_t nblocks_tot_metadata = block_descr_->nblocks_tot_at_depth(md_depth);
    usi d_depth = 1;
    size_t nblocks_tot_data = block_descr_->nblocks_tot_at_depth(d_depth);
    size_t max_nblocks_metadata = max_nblocks_data * nblocks_tot_metadata / nblocks_tot_data;

    //bad things can happen with really small tile sizes
    size_t overflow_buffer = 10000;
    metadata_block_size_ =
              overflow * max_nblocks_data * sizeof(DataNode)
            + sizeof(TensorBranch)
            + max_nblocks_metadata * (sizeof(MetaDataNode) + sizeof(NodeMap<TileNode>))
            + overflow_buffer;

#ifdef MAX_TENSOR_METADATA_BLOCK_SIZE
    if (metadata_block_size_ > MAX_TENSOR_METADATA_BLOCK_SIZE)
        metadata_block_size_ = MAX_TENSOR_METADATA_BLOCK_SIZE;
#endif
    size_t nblocks_metadata = DEFAULT_TENSOR_CACHE_SIZE / metadata_block_size_;
    size_t total_metadata_cache_size = DEFAULT_TENSOR_CACHE_SIZE;
    if (nblocks_metadata > MAX_TENSOR_CACHE_BLOCKS)
        total_metadata_cache_size = metadata_block_size_ * MAX_TENSOR_CACHE_BLOCKS;
    else if (nblocks_metadata < 5)
        total_metadata_cache_size = 10 * metadata_block_size_;
    metadata_cache_ = new DataCache(total_metadata_cache_size, metadata_block_size_);

    size_t average_ndata = block_descr_->average_nelements_data();
    size_t max_ndata = block_descr_->max_nelements_data();
    size_t average_nmetadata = block_descr_->average_nelements_metadata();
    size_t max_nmetadata = block_descr_->max_nelements_metadata();
    size_t data_block_nelements = average_nmetadata * average_ndata < max_ndata ? 
                                    max_ndata : average_nmetadata * average_ndata;

    //we want five of the largest data blocks on one controller
    data_block_size_ = data_block_nelements * sizeof(double);

    size_t ndata_blocks = DEFAULT_TENSOR_CACHE_SIZE / data_block_size_;
    //scale cache size down to get to 10000 entries
    size_t total_data_cache_size = DEFAULT_TENSOR_CACHE_SIZE;
    if (ndata_blocks > MAX_TENSOR_CACHE_BLOCKS)
        total_data_cache_size = data_block_size_ * MAX_TENSOR_CACHE_BLOCKS;
    else if (ndata_blocks < 10)
        total_data_cache_size = 25 * data_block_size_;
    data_cache_ = new DataCache(total_data_cache_size, data_block_size_);

    fill_blocks();

#if 0
    dout << "Tensor " << name_ << " metadata block size  " << metadata_block_size_ << endl;
    dout << "Tensor " << name_ << " metadata cache size  " << total_metadata_cache_size << endl;
    dout << "Tensor " << name_ << " data block size      " << data_block_size_ << endl;
    dout << "Tensor " << name_ << " data cache size      " << total_data_cache_size << endl;
    dout << "Tensor " << name_ << " max num data         " << max_ndata << endl;
    dout << "Tensor " << name_ << " av num data          " << average_ndata << endl;
    dout << "Tensor " << name_ << " max num metadata     " << max_nmetadata << endl;
    dout << "Tensor " << name_ << " av num metadata      " << average_nmetadata << endl;
#endif
}

Tensor::Tensor(
    const TensorIndexDescrPtr& descr,
    Tensor* parent_tensor
) :
    name_(parent_tensor->get_name()),
    tensor_grp_(parent_tensor->get_tensor_grp()),
    original_grp_(parent_tensor->get_original_grp()),
    tensor_blocks_(0),
    block_descr_(descr),
    depth_(descr->depth()),
    sort_weight_(1),
    storage_type_(parent_tensor->get_storage_type()),
    distribution_indices_(parent_tensor->distribution_indices_),
    priority_(parent_tensor->get_priority()),
    filler_(parent_tensor->get_element_computer()),
    parent_tensor_(parent_tensor),
    tensor_grp_generator_set_(0),
    metadata_block_size_(parent_tensor->metadata_block_size()),
    data_block_size_(parent_tensor->data_block_size()),
    data_cache_(parent_tensor->get_data_cache()),
    metadata_cache_(parent_tensor->get_metadata_cache()),
    main_indexer_(0),
    block_to_node_indexmap_(0),
    full_index_to_distribution_index_(0),
    malloc_number_(0),
    cxn_position_(parent_tensor->cxn_position_)
{
    malloc_number_ = Malloc<Tensor>::get_malloc_number(this);

    parent_tensor_->incref();

    tensor_blocks_ = new BlockMap<TensorBlock>(descr.get());

    for (usi i=0; i < block_descr_->nindex(); ++i)
    {
        IndexRange* range = block_descr_->get(i)->get_top_range();
        index_start_[i] = range->start();
        index_end_[i] = range->start() + range->nelements();
    }


    TensorBlockMap::iterator it(tensor_blocks_->begin());
    TensorBlockMap::iterator stop(tensor_blocks_->end());
    //loop all of the blocks and attempt to fetch a block from the parent tensor
    uli idx = 0;
    uli indices[NINDEX];
    for ( ; it != stop; ++it, ++idx)
    {
        tensor_blocks_->indices(idx, indices);
        TensorBlock* block = parent_tensor->get_block(indices);
        if (block)
        {
            tensor_blocks_->set(idx, block);
            block->set_subblock(true);
        }
    }

    tensor_grp_generator_set_ = tensor_grp_->get_generator_set();

    if (block_descr_->has_symmetry())
        init_symmetry_filter();
}

Tensor::~Tensor()
{
    remote_wait();

    if (!parent_tensor_)
    {
        TensorBlock** it = tensor_blocks_->begin();
        TensorBlock** stop = tensor_blocks_->end();
        for ( ; it != stop; ++it)
        {
            TensorBlock* block = *it;
            if (!block)
                continue;

            if (!block->is_permutationally_unique())
                block->obsolete();

            block->set_subblock(false);
        }

        if (distribution_indices_)
            delete[] distribution_indices_;

        if (full_index_to_distribution_index_)
            delete[] full_index_to_distribution_index_;

    }

    //make sure all tensor blocks are released
    TensorBlock** it = tensor_blocks_->begin();
    TensorBlock** stop = tensor_blocks_->end();
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;

        //sleep until the block is retrieved
        while (block->is_retrieved())
            usleep(1);
    }

    uli idx = 0;
    delete tensor_blocks_;
    tensor_blocks_ = 0;
    data_cache_ = 0;
    metadata_cache_ = 0;

    if (parent_tensor_)
        boost::intrusive_ptr_release(parent_tensor_);
}

void
Tensor::init_symmetry_filter()
{
    return;
    TensorElementFilterPtr filter = new TensorSymmetryFilter;
    filter->set_index_descr(block_descr_.get());
    filters_.push_back(filter);
}

void
Tensor::reset()
{
    if (parent_tensor_) //I am the parent!
    {
        cerr << "cannot reset child tensor" << endl;
        abort();
    }

    TensorBlock** it = tensor_blocks_->begin();
    TensorBlock** stop = tensor_blocks_->end();
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;

        if (block && !block->is_permutationally_unique())
            block->obsolete();

        boost::intrusive_ptr_release(block);
        *it = 0;
    }
}

void
Tensor::insert_new_block(
    TensorBlock* block
)
{
    const uli* block_indices = block->get_indices();
    uli idx = tensor_blocks_->index(block_indices);
#if YETI_SANITY_CHECK
    TensorBlock* old_block = tensor_blocks_->get(idx);
    if (old_block)
    {
        cerr << "old block " << old_block->get_index_string()
             << " exists at index " << idx
             << " where block " << block->get_index_string()
             << " is attempting to be inserted"
             << endl;
        abort();
    }
#endif
    tensor_blocks_->set(idx, block);
}




TensorBlock*
Tensor::make_block(
    const uli* indexset,
    Tensor::tensor_storage_t store_type
)
{
    size_t size = block_descr_->total_data_size(indexset);
    if (size == 0)
    {
        //nothing here!
        return 0;
    }


    {
        //evaluate whether or not the block is zero
        std::list<TensorElementFilterPtr>::const_iterator it = filters_.begin();
        std::list<TensorElementFilterPtr>::const_iterator stop = filters_.end();
        for ( ; it != stop; ++it)
        {
            bool zero = (*it)->is_zero(indexset);
            if (zero)
                return 0;
        }
    }

    TensorBlock* block = 0;

    if (parent_tensor_)
    {
        block = new TensorBlock(indexset, parent_tensor_, store_type);
        parent_tensor_->insert_new_block(block);
        this->insert_new_block(block);
        block->set_subblock(true);
    }
    else
    {
        block = new TensorBlock(indexset, this, store_type);
        this->insert_new_block(block);
    }

#if 0
    PermutationGroup::iterator it(tensor_grp_->begin());
    PermutationGroup::iterator stop(tensor_grp_->end());
    uli block_indices[NINDEX];
    for ( ; it != stop; ++it)
    {
       Permutation* p = *it;
        if (p->is_identity() || p->fixes_arr(indexset))
            continue;
        
        p->permute(indexset, block_indices);

        if (parent_tensor_)
        {
            //a permutational symmetry might be valid, but the
            //parent tensor might not contain the corresponding block
            if (!parent_tensor_->contains(block_indices))
                continue;

            TensorBlock* nonunique_block = parent_tensor_->get_block_map()->get(block_indices);
            if (nonunique_block) //already exists
                continue;

            nonunique_block = new TensorBlock(block_indices, parent_tensor_, store_type);

            parent_tensor_->insert_new_block(nonunique_block);
            if (this->contains(block_indices))
                this->insert_new_block(nonunique_block);
        }
        else
        {
            if (!this->contains(block_indices))
                continue;

            TensorBlock* nonunique_block = tensor_blocks_->get(block_indices);
            if (nonunique_block) //already exists
                continue;

            nonunique_block = new TensorBlock(block_indices, this, store_type);

            insert_new_block(nonunique_block);
        }

    }
#endif

    return block;
}

void
Tensor::accumulate(
    uli threadnum,
    Tensor* tensor,
    double scale,
    const SortPtr& sort
)
{
    uli nthread = YetiRuntime::nthread_compute();
    NodeMap<TensorBlock>::iterator stop(tensor->end());
    TensorBlock** itsrc = 0;
    if (sort)
    {
        sort->configure(tensor->get_block_map()->sizes());
        TensorBlock** src_blocks = reinterpret_cast<TensorBlock**>(sort->metadata_buffer(depth_));
        sort->sort_noscale<TensorBlock*>(tensor->begin(), src_blocks);
        itsrc = src_blocks;
        stop = itsrc + tensor->get_block_map()->size();
    }
    else
    {
        itsrc = tensor->begin();
    }

    TensorBlockMap::iterator itdst(tensor_blocks_->begin());
    uli indices[NINDEX];
    uli permuted_indices[NINDEX];
    uli tasknum = 0;
    for ( ; itsrc != stop; ++itsrc, ++itdst)
    {

        TensorBlock* srcblock = *itsrc;
        if (!srcblock) //we don't need to accumulate anything
        {
            continue;
        }

        TensorBlock* dstblock = *itdst;
        if (!dstblock)
        {
            raise(SanityCheckError, "all destination blocks should already exist");
        }

        /** Only accumulate to permutationally unique blocks */
        if (dstblock->is_permutationally_unique())
        {
            if (tasknum % nthread != threadnum)
            {
                ++tasknum;
                continue;
            }

            //retrieve the block in case it needs to be fetched from
            //disk or recomputed
            srcblock->retrieve_read();
            dstblock->retrieve_accumulate();
            //fill the block with values
            dstblock->accumulate(srcblock, scale, sort.get());
            dstblock->update();
            //release the block back to disk or wherever it needs to go
            dstblock->release_accumulate();
            srcblock->release_read();

            ++tasknum;
        }
        else if (is_parent_block(dstblock))
        {
            if (tasknum % nthread != threadnum)
            {
                ++tasknum;
                continue;
            }

            //not permutationally unique, but should still be accumualted
            TensorBlock* accblock = dstblock->get_symmetry_unique_block();
            if (!accblock->is_up_to_date())
            {
                continue; //parent block stuff is tricky... hacky fix
            }

            Permutation* sort_perm = sort->get_permutation();
            Permutation* composite_perm = dstblock->get_fetch_permutation()->product(sort_perm);
            sort->configure(composite_perm);
            accblock->retrieve_accumulate();
            srcblock->retrieve_read();
            accblock->accumulate(srcblock, scale, sort.get());
            accblock->update();
            accblock->release_accumulate();
            srcblock->release_read();
            sort->configure(sort_perm);
        }

    }

}

void
Tensor::accumulate_post_process()
{
    TensorBlock** it = tensor_blocks_->begin();
    TensorBlock** stop = tensor_blocks_->end();
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (block && block->is_permutationally_unique())
            block->update();
    }

    it = tensor_blocks_->begin();
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (block && !block->is_permutationally_unique())
        {
            block->lock();
            block->update();
            block->obsolete();
            block->unlock();
        }
    }

    remote_wait();
}


void
Tensor::accumulate(
    Tensor* tensor,
    double scale,
    Permutation* p
)
{
    heisenfxn(Tensor::accumulate);

    if (tensor->get_block_map()->size() != tensor_blocks_->size())
    {
        cerr << "Cannot accumulate tensor.  Block map sizes do not match." << endl;
        abort();
    }

   Permutation* pacc = 0;
    if (p) pacc = p;
    else pacc = tensor->get_tensor_grp()->get_identity();

    ThreadGroup* grp = YetiRuntime::get_thread_grp();
    uli nthread = YetiRuntime::nthread_compute();
    for (uli i=0; i < nthread; ++i)
    {
        Thread* thr = new AccumulateThread(i, tensor, this, scale, pacc);
        grp->add(thr);
    }
    grp->run();
    grp->wait();
    grp->clear();

    if (parent_tensor_)
        parent_tensor_->accumulate_post_process();
    else
        accumulate_post_process();

    heisenfxn(Tensor::accumulate);
}

bool
Tensor::is_distributed() const
{
    return distribution_indices_;
}

bool
Tensor::is_replicated() const
{
    return !is_distributed();
}

bool
Tensor::is_subtensor() const
{
    return parent_tensor_;
}

uli
Tensor::get_unique_id(const uli* indices, Permutation* fetch_perm) const
{
    if (parent_tensor_)
        return parent_tensor_->get_unique_id(indices, fetch_perm);


    if (!fetch_perm || fetch_perm->is_identity())
    {
        return main_indexer_->index(indices);
    }
    else
    {
        uli unique_indices[NINDEX];
        fetch_perm->permute(indices, unique_indices);
        return main_indexer_->index(unique_indices);
    }
}

Permutation*
Tensor::get_fetch_permutation(const uli* indices) const
{

    if (parent_tensor_)
        return parent_tensor_->get_fetch_permutation(indices);

    Permutation* p = original_grp_->get_lowest_permutation(indices);
    return p;
}

uli
Tensor::get_sort_size_weight() const
{
    return sort_weight_;
}

void
Tensor::set_sort_size_weight(uli weight)
{
    sort_weight_ = weight;
}

bool
Tensor::is_closed_tensor() const
{
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;

        TensorBlock* unique_block = block->get_symmetry_unique_block();
        if (!unique_block)
        {
            return false;
        }

        if (!this->contains(unique_block->get_indices()))
        {
            return false;
        }
    }

    return true;
}

void
Tensor::sort(
    Permutation* p
)
{
    heisenfxn(Tensor::sort);
    ThreadGroup* grp = YetiRuntime::get_thread_grp();
    uli nthread = YetiRuntime::nthread_compute();
    SortThread* main_thread = new SortThread(0, this, p);
    grp->add(main_thread);
    for (uli i=1; i < nthread; ++i)
    {
        Thread* thr = new SortThread(i, this, p);
        grp->add(thr);
    }

    grp->run();
    grp->wait();

    metadata_sort(main_thread->get_sort());

    TensorBlockMap::iterator it = tensor_blocks_->begin();
    TensorBlockMap::iterator stop = tensor_blocks_->end();
    for (  ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;

        if (!block->is_permutationally_unique() || block->is_remote_block())
        {
            //we don't need to sort... but we do need to declare it obsolete
            block->obsolete();
            continue;
        }
    }

    grp->clear();

    remote_wait();
    heisenfxn(Tensor::sort);
}

void
Tensor::metadata_sort(const SortPtr& sort)
{
    //sort the tensor blocks
    tensor_blocks_->sort(sort);

    Permutation* p = sort->get_permutation();
    block_descr_->permute(p);

    if (parent_tensor_)
    {
        parent_tensor_->metadata_sort(sort);
        tensor_grp_ = parent_tensor_->get_tensor_grp();
    }
    else
    {
        tensor_grp_ = tensor_grp_->conjugate(p);
        iterator it = begin();
        iterator stop = end();
        for ( ; it != stop; ++it)
        {
            TensorBlock* block = *it;
            if (!block)
                continue;

            block->metadata_sort(sort.get());
        }
    }

    uli tmp[NINDEX];

    usi nindex = block_descr_->nindex();
    size_t cpysize = nindex * sizeof(uli);
    p->permute(index_start_, tmp);
    ::memcpy(index_start_, tmp, cpysize);
    p->permute(index_end_, tmp);
    ::memcpy(index_end_, tmp, cpysize);


    //validate subtensor
    if (parent_tensor_)
    {
        uli indexset[NINDEX];
        iterator it = begin();
        iterator stop = end();
        for ( ; it != stop; ++it)
        {
            TensorBlock* block = *it;
            if (!block || block->is_permutationally_unique())
                continue;

            block->get_fetch_permutation()->permute(block->get_indices(), indexset);
            if (!this->contains(indexset)) //doesn't have parent block!
            {
                block->sort_resort_permutation(sort->get_permutation());
            }


        }
    }
}

void
Tensor::sort(
    uli threadnum,
    const SortPtr& sort
)
{
    uli idx = 0;
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    uli nthread = YetiRuntime::nthread_compute();
    uli tasknum = 0;
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (!block) //we don't need to create another node
            continue;

        if (!block->is_permutationally_unique() || block->is_remote_block())
            continue;

        if (tasknum % nthread == threadnum)
        {
            //retrieve the block in case it needs to be fetched from
            //disk or recomputed
            block->retrieve_write();
            block->data_sort(sort.get());
            //release the block back to disk or wherever it needs to go
            block->release_write();
        }
        ++tasknum;
    }
}

Tensor::iterator
Tensor::begin() const
{
    return tensor_blocks_->begin();
}

Tensor::iterator
Tensor::end() const
{
    return tensor_blocks_->end();
}

bool
Tensor::contains(const uli* indices) const
{
    usi nindex = block_descr_->nindex();
    for (usi i=0; i < nindex; ++i)
    {
        uli idx = indices[i];
        if (idx < index_start_[i] || idx >= index_end_[i])
            return false;            
    }
    return true;
}

void
Tensor::configure(Tensor::tensor_storage_t storage_type)
{
    storage_type_ = storage_type;
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    if (storage_type == Tensor::on_disk)
    {
        std::string name = this->name_ + ".tensor";
        DiskBufferPtr disk_buffer = new DiskBuffer(name);
        for ( ; it != stop; ++it)
        {
            TensorBlock* block = *it;
            if (block)
                block->configure(disk_buffer);
        }

    }
    else
    {
        for ( ; it != stop; ++it)
        {
            TensorBlock* block = *it;
            if (block)
                block->configure(storage_type);
        }
    }
}

void
Tensor::configure(BlockRetrieveAction* block_action)
{
    TensorBlock* block = get_first_nonnull_block();
    if (!block)
    {
        raise(SanityCheckError, "how did you create a tensor with no blocks?");
    }

    configure(Tensor::action);

    if (!block->get_retrieve_action())
    {
        TensorRetrieveActionPtr action = new TensorRetrieveAction(this);
        action->add(block_action);
        NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
        NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
        for ( ; it != stop; ++it)
        {
            TensorBlock* block = *it;
            if (!block)
                continue;
            block->configure(action);
        }
    }
    else
    {
        TensorRetrieveAction* action = block->get_retrieve_action();
        action->add(block_action);
    }
}

void
Tensor::configure(const TensorElementComputerPtr& filler)
{
    filler->set_index_descr(block_descr_.get());
    filler_ = new ThreadedTensorElementComputer(filler);
    storage_type_ = Tensor::recomputed;
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;
        block->configure(Tensor::recomputed);
    }
}

void
Tensor::configure_degeneracy(const PermutationGroupPtr& grp)
{
    PermutationSetPtr pset = grp->get_generator_set();


    uli idx = 0;
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    uli tmpspace[NINDEX];
    uli indices_used[100];
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (!block) //we don't need to create another node
            continue;
            
        const uli* block_indices = block->get_indices();

        uli idx = tensor_blocks_->index(block_indices);
        indices_used[0] = idx;
        usi nperms = 1;
        PermutationSet::iterator it = grp->begin();
        PermutationSet::iterator stop = grp->end();
        for ( ; it != stop; ++it)
        {
            Permutation* p(*it);
            if (p->is_identity())
            {
                continue;
            }
            else if (p->fixes_arr(block_indices))
            {
                //this does not contribute to the degeneracy
                continue;
            }

            p->permute(block_indices, tmpspace);
            if (!this->contains(tmpspace)) //this block is not relevant
                continue;

            if (p->improves_sort(block_indices))
            {
                nperms = 0;
                break;
            }


            //ensure not a repeat
            idx = tensor_blocks_->index(tmpspace);

            bool repeat = false;
            for (usi i=0; i < nperms; ++i)
            {
                if (indices_used[i] == idx)
                {
                    repeat = true;
                    break;
                }
            }

            if (!repeat)
            {
                indices_used[nperms] = idx;
                ++nperms;
            }
        }

        block->set_degeneracy(nperms);
    }
}

void
Tensor::convert(
    const RectMatrixPtr& matrix,
    MatrixConfiguration* config
)
{
    this->convert(matrix.get(), config);
}

void
Tensor::get_matrix(
    RectMatrixPtr& matrix,
    MatrixConfiguration* config
)
{
    if (!matrix)
    {
        uli nelements[NINDEX];
        block_descr_->get_nelements_data(nelements);
        uli nrows = config->get_index()->nrows(nelements);
        uli ncols = config->get_index()->ncols(nelements);
        matrix = new Matrix(nrows,ncols);
    }
    matrix.zero();

    convert(matrix.get(), config);
}

void
Tensor::get_matrix(
    SymmMatrixPtr& matrix,
    MatrixConfiguration* config
)
{
    if (!matrix)
    {
        uli nelements[NINDEX];
        block_descr_->get_nelements_data(nelements);
        uli nrows = config->get_index()->nrows(nelements);
        uli ncols = config->get_index()->ncols(nelements);
        if (nrows != ncols)
        {
            cerr << "matrix is not symmetric" << endl;
            abort();
        }

        matrix = new Matrix(nrows,ncols);
    }

    matrix.zero();

    convert(matrix.get(), config);
}

Tensor*
Tensor::copy()
{

    Tensor* newtensor = clone();
    newtensor->accumulate(this, 1.0);
    return newtensor;
}

Tensor*
Tensor::clone()
{
    PermutationGroup* grp = new PermutationGroup(block_descr_->nindex());
    PermutationGroup::iterator it = grp->begin();
    PermutationGroup::iterator stop = grp->end();
    for ( ; it != stop; ++it)
    {
        grp->add(*it);
    }
    grp->set_closed();

    usi nindex = block_descr_->nindex();
    TensorIndexDescr* new_descr = new TensorIndexDescr(nindex);
    for (usi i=0; i < nindex; ++i)
    {
        new_descr->set(i, block_descr_->get(i));
    }

    Tensor* newtensor = new Tensor(name_, new_descr, grp);
    newtensor->configure(storage_type_);

    return newtensor;
}

void
Tensor::convert(
    Matrix* matrix,
    MatrixConfiguration* config
)
{
    uli idx = 0;
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (!block) //we don't need to create another node
            continue;

        block->retrieve_read();
        block->convert(matrix, config, block_descr_.get());
        block->release_read();
    }
}

void
Tensor::accumulate(
    const SymmMatrixPtr& matrix,
    MatrixConfiguration* config
)
{
    accumulate(matrix.get(), config);
}

void
Tensor::accumulate(
    const RectMatrixPtr& matrix,
    MatrixConfiguration* config
)
{
    accumulate(matrix.get(), config);
}

void
Tensor::accumulate(
    Matrix* matrix,
    MatrixConfiguration* config
)
{
    uli idx = 0;
    uli indexset[NINDEX];
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (!block)
            continue; //null

        if (!block->is_permutationally_unique())
            continue;

        block->retrieve_write();
        block->accumulate(matrix, config, block_descr_.get());
        block->release_write();
    }

    update();
}

void
Tensor::fill(const double* data)
{
    TensorElementComputerPtr comp = new DoubleArrayElementComputer(data);
    fill(comp);
}


void
Tensor::element_op(ElementOp* op)
{
    op->configure(this);

    ThreadGroup* grp = YetiRuntime::get_thread_grp();
    uli nthread = YetiRuntime::nthread_compute();
    for (uli i=0; i < nthread; ++i)
    {
        Thread* thr = new ElementOpThread(i, this, op);
        grp->add(thr);
    }
    grp->run();
    grp->wait();
    grp->clear();

    if (op->do_update_after())
    {
        update();
        if (parent_tensor_)
        {
            parent_tensor_->update_to_subtensor();
        }
    }
}

void
Tensor::element_op(
    uli threadnum,
    ElementOp* op
)
{
    uli idx = 0;
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    uli nthread = YetiRuntime::nthread_compute();
    bool read_mode = false;
    uli indexset[NINDEX];
    uli tasknum = 0;
    for ( ; it != stop; ++it, ++idx)
    {
        if (tasknum % nthread != threadnum)
            continue;

        TensorBlock* block = *it;
        if (!block) //we don't need to create another node
        {
            continue;
        }

        if (!block->is_permutationally_unique())
        {
            continue;
        }

        /**
            Depending on the type of element operation,
            you might need different modes - e.g.
        */
        op->retrieve(block);

        //retrieve the block in case it needs to be fetched from
        //disk or recomputed
        block->retrieve();
        block->element_op(op);

        read_mode = block->in_read_mode();
        if (!read_mode)
        {
            block->update();
            block->sync(); //make sure we are synced with parent storage
        }
        op->release(block);

        ++tasknum;
    }

}


void
Tensor::fill(
    const ThreadedTensorElementComputerPtr& filler
)
{
    timer::Timer::start("fill");

    ThreadGroup* grp = YetiRuntime::get_thread_grp();
    uli nthread = YetiRuntime::nthread_compute();
    for (uli i=0; i < nthread; ++i)
    {
        Thread* thr = new TensorFillThread(i, this, filler);
        grp->add(thr);
    }
    grp->run();
    grp->wait();
    grp->clear();


    remote_wait();

    timer::Timer::stop("fill");

}

bool
Tensor::equals(
    const void* data
)
{
    uli idx = 0;
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    
    const void* dataptr = data;
    
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (!block || !block->is_permutationally_unique())
            continue;

        block->retrieve_read();

        //dataptr gets incremented
        bool check = block->equals(&dataptr);
        
        if (!check)
        {
            block->release_read();
            return false;
        }
        
        block->release_read();
    }
    
    return true;
}

bool
Tensor::equals(Tensor* tensor)
{
    if (tensor->get_block_map()->size() != tensor_blocks_->size())
        return false;

    uli idx = 0;
    NodeMap<TensorBlock>::iterator it_me(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    NodeMap<TensorBlock>::iterator 
        it_other(tensor->get_block_map()->begin());
    
    for ( ; it_me != stop; ++it_me, ++it_other, ++idx)
    {
        TensorBlock* block_me = *it_me;
        TensorBlock* block_other = *it_other;
        
        bool null_me = block_me;
        bool null_other = block_other;

        if (null_me && !null_other)
            return false;
        else if (null_other && !null_me)
            return false;
        else if (null_other && null_me) //equal! both zero!
            continue;

        if (!block_me->is_permutationally_unique() 
            && !block_other->is_permutationally_unique())
            continue;
        
        block_me->retrieve_read();
        block_other->retrieve_read();
        bool check = block_me->equals(block_other);
        block_me->release_read();
        block_other->release_read();
        
        if (!check)
            return false;
    }
    
    return true;
}

void
Tensor::fill_blocks()
{
    uli idx = 0;
    uli indexset[NINDEX];
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (block)
        {
            if (block->is_permutationally_unique())
            {
                raise(SanityCheckError, "cannot re-fill tensor");
            }
        } //we don't need to create another node

        tensor_blocks_->indices(idx, indexset);
        /** Only do permutationally unique blocks! */
        if (tensor_grp_->improves_sort(indexset))
            continue;

        //do not initialize the blocks on a fill
        block = make_block(indexset, storage_type_);
    }
    remote_wait();
}

void
Tensor::global_sum()
{
    heisenfxn(Tensor::global_sum);
    if (is_distributed())
        raise(SanityCheckError, "cannot global sum a distributed tensor");

    if (YetiRuntime::nproc() == 1)
        return;

    YetiMessenger* messenger = YetiRuntime::get_messenger();

    TensorBlockMap::iterator it = tensor_blocks_->begin();
    TensorBlockMap::iterator stop = tensor_blocks_->end();
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;

        if (messenger->nchildren() == 0) //I am the bottom node in the binary tree
        {
            block->retrieve_read();
            /** It is very important first to retrieve the block
                and then reset the retrieve flag.  The retrieve
                will block if the received flag is set to false
            */
            block->reset_parent_signal();

            SendStatus* status = block->send_data(
                    messenger,
                    Message::TensorBranch,
                    Message::GlobalSum,
                    messenger->parent()
            );
            status->wait();

            block->release_read();
        }
        else
        {
            while (block->nchild_signals() < messenger->nchildren())
            {
                usleep(1);
            }
            block->reset_child_signal();
        }
    }

    it = tensor_blocks_->begin();
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;

        while (!block->received_parent_signal())
            usleep(1);

    }

    remote_wait();
    heisenfxn(Tensor::global_sum);
}

void
Tensor::fill(
    uli threadnum,
    const ThreadedTensorElementComputerPtr& threaded_filler
)
{
    uli tasknum = 0;
    uli nthread = YetiRuntime::nthread_compute();
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    TensorElementComputer* filler = threaded_filler->get_computer(threadnum);
    TensorValueEstimater* estimater = filler->get_estimater(depth_); 
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;  //some blocks can be determined rigorously zero ahead of time
        if (!block)
            continue;
        else if (!block->is_permutationally_unique())
            continue;
        else if (block->is_remote_block())
            continue;

        float maxlog = estimater ? estimater->max_log(block->get_indices()) : LOG_NONZERO;
        if (maxlog < YetiRuntime::nonnull_cutoff) //this contribution can be skipped
            continue;

        if (tasknum % nthread == threadnum)
        {
            //retrieve the block in case it needs to be fetched from
            //disk or recomputed
            block->retrieve_write();
            //fill the block with values
            block->fill(filler);

            //after the fill, we are done with the block
            //so we need to ensure that it is synced with the
            //parent storage area
            block->sync();

            //release the block back to disk or wherever it needs to go
            block->release_write();
        }

        ++tasknum;
   }


}

void
Tensor::fill(const MatrixPtr& matrix)
{
    TensorElementComputerPtr filler
            = new YetiMatrixElementComputer(matrix);
    fill(filler);
}

void
Tensor::fill(const TensorElementComputerPtr& val)
{
    if (val)
    {
        val->set_index_descr(block_descr_.get());
        fill(new ThreadedTensorElementComputer(val));
    }
    else
    {
        TensorElementComputerPtr filler
            = new MemsetElementComputer(TemplateInfo::double_type);
        filler->set_index_descr(block_descr_.get());
        fill(new ThreadedTensorElementComputer(filler));
    }
}

TensorBlock*
Tensor::get_block(uli index) const
{
    return tensor_blocks_->get(index);
}

TensorBlock*
Tensor::get_block(const uli* indices, Permutation* p) const
{
    if (parent_tensor_)
        return parent_tensor_->get_block(indices, p);

    if (p)
    {
        uli permuted_indices[NINDEX];
        p->permute(indices, permuted_indices);
        return tensor_blocks_->get(permuted_indices);
    }
    else
    {
        return tensor_blocks_->get(indices);
    }
}

TensorBlockMap*
Tensor::get_block_map() const
{
    return tensor_blocks_;
}

DataCache*
Tensor::get_data_cache() const
{
    return data_cache_.get();
}

DataCache*
Tensor::get_metadata_cache() const
{
    return metadata_cache_.get();
}

TensorIndexDescr*
Tensor::get_block_descr() const
{
    return block_descr_.get();
}

usi
Tensor::get_depth() const
{
    return depth_;
}

ThreadedTensorElementComputer*
Tensor::get_element_computer() const
{
    return filler_.get();
}

TensorBlock*
Tensor::get_first_nonnull_block() const
{
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (block)
            return block;
    }

    return 0;
}

uli
Tensor::get_node_number(uli block_number) const
{
    if (is_distributed())
    {
        uli subidx = full_index_to_distribution_index_[block_number];
        return subidx % YetiRuntime::nproc();
    }
    else
        return 0;
}

uli
Tensor::get_index(const uli *indices) const
{
    return tensor_blocks_->index(indices);
}

const std::string&
Tensor::get_name() const
{
    return name_;
}

Tensor::tensor_priority_t
Tensor::get_priority() const
{
    return priority_;
}

Tensor::tensor_storage_t
Tensor::get_storage_type() const
{
    return storage_type_;
}

PermutationGroup*
Tensor::get_tensor_grp() const
{
    return tensor_grp_.get();
}
PermutationGroup*
Tensor::get_original_grp() const
{
    return original_grp_.get();
}

ulli
Tensor::get_totalsize() const
{
    return block_descr_->totalsize();
}

size_t
Tensor::metadata_block_size() const
{
    return metadata_block_size_;
}

size_t
Tensor::data_block_size() const
{
    return data_block_size_;
}

void
Tensor::distribute(const usi* indices, usi nindex_distr)
{
    if (YetiRuntime::nproc() == 1)
        return;

    if (distribution_indices_)
    {
        raise(SanityCheckError, "cannot redestribute tensor");
    }

    distribution_indices_ = new usi[nindex_distr];
    nindex_distr_ = nindex_distr;
    for (usi i=0; i < nindex_distr; ++i)
        distribution_indices_[i] = indices[i];

    uli sizes[NINDEX];
    uli offsets[NINDEX];
    for (usi i=0; i < nindex_distr; ++i)
    {
        sizes[i] = tensor_blocks_->sizes()[distribution_indices_[i]];
        offsets[i] = tensor_blocks_->offsets()[distribution_indices_[i]];
    }

    distr_indexer_ = new Indexer(sizes, offsets, nindex_distr);
    uli nblocks = tensor_blocks_->size();
    full_index_to_distribution_index_ = new uli[nblocks];

    uli distr_indexset[NINDEX];
    uli full_indexset[NINDEX];
    uli permuted_indexset[NINDEX];
    for (uli idx=0; idx < nblocks; ++idx)
    {
        tensor_blocks_->indices(idx, full_indexset);
        Permutation* p = tensor_grp_->get_lowest_permutation(full_indexset);
        p->permute(full_indexset, permuted_indexset);

        for (usi i=0; i < nindex_distr; ++i)
        {
            distr_indexset[i] = permuted_indexset[distribution_indices_[i]];
        }

        uli subset_idx = distr_indexer_->index(distr_indexset);

        full_index_to_distribution_index_[idx] = subset_idx;

        TensorBlock* block = tensor_blocks_->get(idx);
        if (block)
        {
            uli node_assignment = subset_idx % YetiRuntime::nproc();
            block->set_node_number(node_assignment);
        }

    }

    remote_wait();
}

TensorBlock*
Tensor::get_make_block(uli index)
{
    if (index >= tensor_blocks_->size())
    {
        std::string err = stream_printf(
            "invalid block index %d passed to tensor %s",
            index, name_.c_str()
        );
        raise(SanityCheckError, err);
    }

    TensorBlock* block = tensor_blocks_->get(index);
    if (block)
        return block;

    uli indices[NINDEX];
    //doesn't exist yet... build it        
    tensor_blocks_->indices(index, indices);
    block = make_block(indices);

    return block;
}

bool
Tensor::intersects(Tensor* tensor) const
{
    if (this == tensor)
        return true;
    else if (this->parent_tensor_ == tensor)
        return true;
    else if (this == tensor->parent_tensor_)
        return true;
    else if (this->parent_tensor_ != 0 && this->parent_tensor_ == tensor->parent_tensor_)
        return true;
    else
        return false;
}

void
Tensor::internal_contraction(
    Tensor* dst_tensor,
    MatrixConfiguration* config
)
{
    config->configure_block(
        tensor_blocks_->sizes(),
        get_depth()
    );
    uli nr = config->nrows();
    uli nc = config->ncols();

    TensorBlockMap* dst_map = dst_tensor->get_block_map();
    TensorBlockMap* src_map = this->get_block_map();

    TensorBlock** it_dst = dst_map->begin();
    TensorBlock** it_src = src_map->begin();

    uli dst_indices[NINDEX];
    for (uli r=0; r < nr; ++r, ++it_dst)
    {
        TensorBlock* dst_block = *it_dst;
        dst_map->indices(r, dst_indices);
        for (uli c=0; c < nc; ++c, ++it_src)
        {
            TensorBlock* src_block = *it_src;
            if (!src_block)
                continue;

            if (!dst_block)
                dst_block = dst_tensor->make_block(dst_indices);
            if (!dst_block)
                continue; //filtered out

            src_block->retrieve_read();
            dst_block->retrieve_accumulate();
            src_block->internal_contraction(
                dst_block,
                config
            );
            src_block->release_read();
            dst_block->release_accumulate();
        }
    }

    dst_tensor->update();
}

bool
Tensor::nonzero() const
{
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    uli idx = 0;
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;

        block->retrieve();
        float maxlog = block->get_branch()->get_node()->get_max_log();
        block->release();

        if (maxlog > YetiRuntime::nonnull_cutoff)
        {
            return true;
        }
    }
    return false;
}

bool
Tensor::is_parent_block(TensorBlock* block) const
{
    if (block->is_permutationally_unique())
    {
        return true;
    }

    uli indices[NINDEX];
    block->get_fetch_permutation()->permute(block->get_indices(), indices);
    bool contains_parent = this->contains(indices);

    return !contains_parent;
}

double
Tensor::norm()
{
    NormElementOp* op = new NormElementOp;
    element_op(op);
    return op->norm();
}

ulli
Tensor::nelements() const
{
    return block_descr_->totalsize();
}

uli
Tensor::nblocks_retrieved() const
{
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    uli nretrieved = 0;
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (block && block->is_retrieved())
            ++nretrieved;
    }
    return nretrieved;
}

ulli
Tensor::nelements_unique() const
{
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    uli idx = 0;
    uli indexset[NINDEX];
    ulli n = 0;
    for ( ; it != stop; ++it, ++idx)
    {
        tensor_blocks_->indices(idx, indexset);
        if (tensor_grp_->improves_sort(indexset))
            continue;

        n += block_descr_->get_nelements_data(depth_, indexset);
    }
    return n;
}

void
Tensor::print_block(TensorBlock* block, std::ostream& os) const
{
    const uli* indexset = block->get_indices();
    os << endl << "Tensor Block " << block->get_index_string()
       << " on node " << block->get_node_number() << endl;
    os << "irreps:";
    for (usi i=0; i < block_descr_->nindex(); ++i)
        os << " " << block_descr_->get(i)->irrep(indexset[i]);
    os << endl;
    os << "Fetch permutation:  ";
    block->get_fetch_permutation()->print(os); os << endl;
    os << "Resort permutation: ";
    block->get_resort_permutation()->print(os); os << endl;
    ++Env::indent;
    block->retrieve_read();
    block->print(os);
    block->release_read();
    --Env::indent;
}

void
Tensor::print(std::ostream& os)
{
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());

    if (YetiRuntime::me() == 0)
    {
        os << "Tensor " << name_ << endl;
        os << "Permutation Group" << endl;
        os << tensor_grp_ << endl;
        os << "Original Group" << endl;
        os << original_grp_ << endl;
    }

    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;

        if (!is_distributed() || YetiRuntime::me() == block->get_node_number())
            print_block(block, os);
    }

}

void
Tensor::scale(double scale)
{
    ElementOp* op = new ScaleOp(scale);
    element_op(op);
    delete op;
}


void
Tensor::zero()
{
    ElementOp* op = new ZeroOp;
    element_op(op);
    delete op;
}

void
Tensor::set_tensor_position(Contraction::tensor_position_t position)
{
    cxn_position_ = position;
}

Contraction::tensor_position_t 
Tensor::get_tensor_position() const
{
    return cxn_position_;
}

void
Tensor::set_priority(Tensor::tensor_priority_t priority)
{
    priority_ = priority;
    if (parent_tensor_)
        parent_tensor_->set_priority(priority);
}

void
Tensor::remote_wait()
{
    heisenfxn(Tensor::remote_wait);
    if (YetiRuntime::nproc() > 1)
    {
        heisenfxn(Tensor::remote_wait);
        timer::Timer::start("tensor barrier");
        YetiRuntime::get_messenger()->wait_barrier();
        timer::Timer::stop("tensor barrier");
    }
    heisenfxn(Tensor::remote_wait);
}

void
Tensor::update_remote_blocks()
{
}

void
Tensor::update_from_subtensor()
{
    if (parent_tensor_)
    {
        cerr << "should not be called on subtensor" << endl;
        abort();
    }

    iterator it = tensor_blocks_->begin();
    iterator stop = tensor_blocks_->end();
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;

        //this is parent of parent tensor but maps onto a block
        //that was changed in the course of some operation
        if (!block->is_subblock() && !block->is_permutationally_unique())
        {
            block->obsolete();
        }
    }

}

void
Tensor::update_to_subtensor()
{
    if (parent_tensor_)
    {
        cerr << "should not be called on subtensor" << endl;
        abort();
    }

}

void
Tensor::flush()
{
    foreach_nonnull(block, tensor_blocks_, TensorBlock,
        block->flush_from_cache();
    )
}

void
Tensor::sync()
{
    foreach_nonnull(block, tensor_blocks_, TensorBlock,
        block->sync();
    )

    remote_wait();
}


void
Tensor::reset_degeneracy()
{
    foreach_nonnull(block, tensor_blocks_, TensorBlock,
        block->reset_degeneracy();
    )
}

void
Tensor::update(uli threadnum)
{
    uli idx = 0;
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    uli nthread = YetiRuntime::nthread_compute();
    uli tasknum = 0;
    for ( ; it != stop; ++it, ++idx)
    {
        if (tasknum % nthread != threadnum)
            continue;

        TensorBlock* block = *it;
        if (!block)
            continue;
        
        if (is_distributed() && block->get_node_number() != YetiRuntime::me())
            continue;

        if (block->is_permutationally_unique())
        {
            block->retrieve_read();
            block->update();
            block->release_read();
        }

        ++tasknum;
    }
}

void
Tensor::recompute_permutation()
{
    
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;
        block->recompute_permutation();
    }
}

void
Tensor::update()
{    
    ThreadGroup* grp = YetiRuntime::get_thread_grp();
    uli nthread = YetiRuntime::nthread_compute();
    for (uli i=0; i < nthread; ++i)
    {
        /** These threads will update the symmetry unique blocks */
        Thread* thr = new TensorUpdateThread(i, this);
        grp->add(thr);
    }
    grp->run();
    grp->wait();
    grp->clear();


    /** Now go through and work on the symmetry non-unique blocks */
    uli idx = 0;
    NodeMap<TensorBlock>::iterator it(tensor_blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(tensor_blocks_->end());
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;

        if (!block->is_permutationally_unique())
        {
            block->reset_degeneracy();
            block->obsolete();
        }
    }

    if (parent_tensor_)
        parent_tensor_->update_from_subtensor();

    remote_wait();
}


bool
Tensor::is_cache_coherent() const
{
    foreach_nonnull(block, tensor_blocks_, TensorBlock,
        if (!block->is_cache_coherent())
            return false;
    )
    return true;
}
