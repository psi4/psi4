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

bool Tensor::in_destructor = false;
bool Tensor::in_sync = false;

bool
yeti::tensor_less(Tensor* l, Tensor* r)
{
    uli lpriority = l->get_priority();
    uli rpriority = r->get_priority();
    if (lpriority != rpriority)
        return lpriority < rpriority;


    Tensor::tensor_distribution_t l_dist_type = l->get_distribution_type();
    Tensor::tensor_distribution_t r_dist_type = r->get_distribution_type();
    if (l_dist_type != r_dist_type)
        return l_dist_type < r_dist_type;

    //in core is always faster than recomputing which is probably faster than going to disk
    Tensor::tensor_storage_t l_store_type = l->get_storage_type();
    Tensor::tensor_storage_t r_store_type = r->get_storage_type();
    if (l_store_type != r_store_type)
        return l_store_type < r_store_type;

    //everything is the same now, select based on size
    ulli lsize = l->get_totalsize();
    uli rsize = r->get_totalsize();
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
Tensor::Tensor(
    const std::string& name,
    const TensorIndexDescrPtr& descr,
    const PermutationGroupPtr& pgrp
) :
    name_(name),
    tensor_grp_(pgrp),
    original_grp_(pgrp),
    blocks_(0),
    descr_(descr),
    depth_(descr->depth()),
    storage_type_(Tensor::in_core),
    distribution_type_(Tensor::replicated),
    erase_type_(Tensor::delete_zero_blocks),
    disk_buffer_(0),
    priority_(Tensor::gamma_tensor),
    filler_(0),
    parent_tensor_(0),
    sort_perm_(0),
    data_cache_(0),
    metadata_cache_(0),
    metadata_block_size_(0),
    data_block_size_(0),
    main_indexer_(new Indexer(descr.get()))
{
    blocks_ = new BlockMap<TensorBlock>(descr.get());
    
    for (usi i=0; i < descr_->nindex(); ++i)
    {
        IndexRange* range = descr_->get(i)->get_top_range();
        index_start_[i] = range->start();
        index_end_[i] = range->start() + range->nelements();
    }

    if (descr_->has_symmetry())
        init_symmetry_filter();

    sort_perm_ = pgrp->get_identity();
    unsort_perm_ = sort_perm_;

    usi overflow = 2;
    size_t max_nblocks_data = descr_->max_nblocks_data();
    if (max_nblocks_data < 4)
        max_nblocks_data = 4;


    usi md_depth = 2;
    size_t nblocks_tot_metadata = descr_->nblocks_tot_at_depth(md_depth);
    usi d_depth = 1;
    size_t nblocks_tot_data = descr_->nblocks_tot_at_depth(d_depth);
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
    if (nblocks_metadata > NCACHE_ENTRIES_MAX)
        total_metadata_cache_size = metadata_block_size_ * NCACHE_ENTRIES_MAX;
    else if (nblocks_metadata < 5)
        total_metadata_cache_size = 10 * metadata_block_size_;
    metadata_cache_ = new DataCache(total_metadata_cache_size, metadata_block_size_);

    size_t average_ndata = descr_->average_nelements_data();
    size_t max_ndata = descr_->max_nelements_data();
    size_t average_nmetadata = descr_->average_nelements_metadata();
    size_t max_nmetadata = descr_->max_nelements_metadata();
    size_t data_block_nelements = average_nmetadata * average_ndata < max_ndata ? 
                                    max_ndata : average_nmetadata * average_ndata;

    //we want five of the largest data blocks on one controller
    data_block_size_ = data_block_nelements * sizeof(double);

    size_t ndata_blocks = DEFAULT_TENSOR_CACHE_SIZE / data_block_size_;
    //scale cache size down to get to 10000 entries
    size_t total_data_cache_size = DEFAULT_TENSOR_CACHE_SIZE;
    if (ndata_blocks > NCACHE_ENTRIES_MAX)
        total_data_cache_size = data_block_size_ * NCACHE_ENTRIES_MAX;
    else if (ndata_blocks < 10)
        total_data_cache_size = 50 * data_block_size_;
    data_cache_ = new DataCache(total_data_cache_size, data_block_size_);

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
    blocks_(0),
    descr_(descr),
    depth_(descr->depth()),
    storage_type_(parent_tensor->get_storage_type()),
    distribution_type_(parent_tensor->get_distribution_type()),
    erase_type_(parent_tensor->get_erase_type()),
    disk_buffer_(parent_tensor->get_disk_buffer()),
    priority_(parent_tensor->get_priority()),
    filler_(parent_tensor->get_element_computer()),
    parent_tensor_(parent_tensor),
    tensor_grp_generator_set_(0),
    sort_perm_(parent_tensor->get_sort_permutation()),
    unsort_perm_(parent_tensor_->get_unsort_permutation()),
    metadata_block_size_(parent_tensor->metadata_block_size()),
    data_block_size_(parent_tensor->data_block_size()),
    data_cache_(parent_tensor->get_data_cache()),
    metadata_cache_(parent_tensor->get_metadata_cache()),
    main_indexer_(0)
{
    parent_tensor_->incref();
    blocks_ = new BlockMap<TensorBlock>(descr.get());
    
    for (usi i=0; i < descr_->nindex(); ++i)
    {
        IndexRange* range = descr_->get(i)->get_top_range();
        index_start_[i] = range->start();
        index_end_[i] = range->start() + range->nelements();
    }

    TensorBlockMap::iterator it(blocks_->begin());
    TensorBlockMap::iterator stop(blocks_->end());
    //loop all of the blocks and attempt to fetch a block from the parent tensor
    uli idx = 0;
    uli indices[NINDEX];
    for ( ; it != stop; ++it, ++idx)
    {
        blocks_->indices(idx, indices);
        TensorBlock* block = parent_tensor->get_block(indices);
        if (block)
        {
            blocks_->set(idx, block);
            block->set_subblock(true);
        }
    }

    tensor_grp_generator_set_ = tensor_grp_->get_generator_set();

    if (descr_->has_symmetry())
        init_symmetry_filter();
}

Tensor::~Tensor()
{
    if (YetiRuntime::is_threaded_runtime())
    {
        cerr << "threaded runtime in tensor destructor!" << endl;
        abort();
    }

    in_destructor = true;
    size_t wasted = 0;
    if (!parent_tensor_)
    {
        TensorBlock** it = blocks_->begin();
        TensorBlock** stop = blocks_->end();
        for ( ; it != stop; ++it)
        {
            TensorBlock* block = *it;
            if (!block)
                continue;

            if (!block->is_permutationally_unique())
                block->obsolete();

            //if (block->is_permutationally_unique())
            //    wasted += block->get_branch()->size_data_wasted();

            block->set_subblock(false);
        }
        size_t total = sizeof(double) * get_totalsize();
        //dout << stream_printf("wasted %ld bytes on tensor %s out of %ld", 
        //    wasted, name_.c_str(), total) << endl;
    }
    uli idx = 0;
    delete blocks_;
    blocks_ = 0;
    in_destructor = false;
    data_cache_ = 0;
    metadata_cache_ = 0;

    if (parent_tensor_)
        boost::intrusive_ptr_release(parent_tensor_);
}

void
Tensor::init_symmetry_filter()
{
    TensorElementFilterPtr filter = new TensorSymmetryFilter;
    filter->set_index_descr(descr_.get());
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

    TensorBlock** it = blocks_->begin();
    TensorBlock** stop = blocks_->end();
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
    uli idx = blocks_->index(block->get_indices());
    insert_new_block(idx, block);
}

void
Tensor::insert_new_block(uli idx, TensorBlock* block)
{
#if YETI_SANITY_CHECK
    TensorBlock* old_block = blocks_->get(idx);
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

    blocks_->set(idx, block);
}



TensorBlock*
Tensor::make_block(
    const uli* indexset,
    Tensor::tensor_storage_t store_type
)
{
    size_t size = descr_->total_data_size(indexset);
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
            TensorBlock* nonunique_block = blocks_->get(block_indices);
            if (nonunique_block) //already exists
                continue;

            nonunique_block = new TensorBlock(block_indices, this, store_type);

            insert_new_block(nonunique_block);
        }

    }

#if YETI_SANITY_CHECK
    //now do a sanity check
    if (!block->is_permutationally_unique())
    {
        block->get_fetch_permutation()->permute(block->get_indices(), block_indices);
        if (!this->contains(block_indices) &&
             (parent_tensor_ && !parent_tensor_->contains(block_indices)))
        {
            cerr << "created block " << block->get_index_string()
                 << " " << block->get_fetch_permutation()
                 << " but unique parent block does not exist!"
                 << endl;
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
    uli idx = 0;
    uli nthread = YetiRuntime::nthread();
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

    TensorBlockMap::iterator itdst(blocks_->begin());
    uli indices[NINDEX];
    for ( ; itsrc != stop; ++itsrc, ++itdst, ++idx)
    {
        if (INVALID_THREAD_TASK(idx, threadnum, nthread))
            continue;

        TensorBlock* srcblock = *itsrc;
        if (!srcblock) //we don't need to accumulate anything
            continue;

        TensorBlock* dstblock = *itdst;
        if (!dstblock)
        {
            blocks_->indices(idx, indices);
            if (tensor_grp_->improves_sort(indices)) //only permutationally unique!
                continue;

            dstblock = make_block(indices);
        }

        /** Only accumulate to permutationally unique blocks */
        if (dstblock->is_permutationally_unique())
        {
            //set the owner thread for the given task
            bool test = dstblock->set_accumulate_mode();
            if (!test)
            {
                cerr << stream_printf("set accumulate failed %d %s", __LINE__, __FILE__) << endl;
                abort();
            }
            srcblock->set_read_mode();
            //retrieve the block in case it needs to be fetched from
            //disk or recomputed
            srcblock->retrieve();
            dstblock->retrieve();
            //fill the block with values
            dstblock->accumulate(srcblock, scale, sort.get());
            dstblock->update();
            //release the block back to disk or wherever it needs to go
            dstblock->release();
            srcblock->release();
        }
    }

}

void
Tensor::accumulate_post_process()
{
    TensorBlock** it = blocks_->begin();
    TensorBlock** stop = blocks_->end();
    uli idx = 0;
    uli nthread = YetiRuntime::nthread();
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;

        if (block->is_permutationally_unique())
            continue;

        block->sync_max_log();
        block->obsolete();
    }
}


void
Tensor::accumulate(
    Tensor* tensor,
    double scale,
    Permutation* p
)
{


    if (tensor->get_block_map()->size() != blocks_->size())
    {
        cerr << "Cannot accumulate tensor.  Block map sizes do not match." << endl;
        abort();
    }

   Permutation* pacc = 0;
    if (p) pacc = p;
    else pacc = tensor->get_tensor_grp()->get_identity();

    ThreadGroup* grp = YetiRuntime::get_thread_grp();
    uli nthread = YetiRuntime::nthread();
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

#if 0//YETI_SANITY_CHECK
    if (!is_closed_tensor())
    {
        cerr << "Cannot accumulate to a tensor with a symmetry non-unique block "
             << " without a matching parent tensor block " << endl;
        abort();
    }
#endif
}

bool
Tensor::is_subtensor() const
{
    return parent_tensor_;
}

bool
Tensor::is_unique(const uli* indices) const
{
    if (parent_tensor_)
        return parent_tensor_->is_unique(indices);

    uli tmpspace[NINDEX];
    unsort_perm_->permute(indices, tmpspace);

    Permutation* p = original_grp_->get_lowest_permutation(tmpspace);
    if (p->is_identity())
        return true;

    uli tmp2[NINDEX];
    p->permute(tmpspace, tmp2);
    sort_perm_->permute(tmp2, tmpspace);
    bool unique_parent = this->contains(tmpspace);

    return !unique_parent;
}

uli
Tensor::get_unique_id(const uli* indices, Permutation* fetch_perm) const
{
    if (parent_tensor_)
        return parent_tensor_->get_unique_id(indices, fetch_perm);

    if (unsort_perm_->is_identity())
    {
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
    else
    {
        uli unsorted_unique_indices[NINDEX];
        if (!fetch_perm || fetch_perm->is_identity())
        {
            unsort_perm_->permute(indices, unsorted_unique_indices);
        }
        else
        {
            uli unique_indices[NINDEX];
            fetch_perm->permute(indices, unique_indices);
            unsort_perm_->permute(unique_indices, unsorted_unique_indices);
        }
        return main_indexer_->index(unsorted_unique_indices);
    }
}

Permutation*
Tensor::get_fetch_permutation(const uli* indices) const
{
   Permutation* p = 0;
    if (parent_tensor_)
    {
        p = parent_tensor_->get_fetch_permutation(indices);
    }
    else if (unsort_perm_->is_identity())
    {
        p = original_grp_->get_lowest_permutation(indices);
    }
    else
    {
        uli tmpspace[NINDEX];
        unsort_perm_->permute(indices, const_cast<uli*>(tmpspace));
       Permutation* q = original_grp_->get_lowest_permutation(tmpspace);
        p = sort_perm_->product(q->product(unsort_perm_));
    }
    return p;
}

bool
Tensor::is_closed_tensor() const
{
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;

        if (block->get_max_log() < YetiRuntime::nonnull_cutoff)
            continue;

        if (block->is_permutationally_unique())
        {
            //continue;
        }

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
    ThreadGroup* grp = YetiRuntime::get_thread_grp();
    uli nthread = YetiRuntime::nthread();
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

    grp->clear();


}

void
Tensor::metadata_sort(const SortPtr& sort)
{
    //sort the tensor blocks
    blocks_->sort(sort);

    Permutation* p = sort->get_permutation();
    descr_->permute(p);
    if (parent_tensor_)
    {
        /** Subtensors if resorted must be immediately resorted back
            to original form. Any other permutation is invalid. */
        bool subtensor_resort = !sort_perm_->is_identity();

        parent_tensor_->metadata_sort(sort);
        tensor_grp_ = parent_tensor_->get_tensor_grp();
        sort_perm_ = parent_tensor_->get_sort_permutation();
        unsort_perm_ = parent_tensor_->get_unsort_permutation();

        if (subtensor_resort && !sort_perm_->is_identity())
        {
        }
    }
    else
    {
        tensor_grp_ = tensor_grp_->conjugate(p);
        sort_perm_ = p->product(sort_perm_);
        unsort_perm_ = sort_perm_->inverse();
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

    usi nindex = descr_->nindex();
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

#if YETI_SANITY_CHECK
            Permutation* ptest = block->get_fetch_permutation()->product(block->get_resort_permutation());
            if (sort_perm_->is_identity() && !ptest->is_identity())
            {
                cerr << "nonunique subblock does not have sort realigned to parent block" << endl;
                abort();
            }
#endif

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
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
    uli nthread = YetiRuntime::nthread();
    for ( ; it != stop; ++it, ++idx)
    {
        if (INVALID_THREAD_TASK(idx, threadnum, nthread))
            continue;

        TensorBlock* block = *it;
        if (!block) //we don't need to create another node
            continue;


        if (!block->is_permutationally_unique())
        {
            //we don't need to sort... but we do need to declare it obsolete
            block->obsolete();
            continue;
        }

        block->set_write_mode();
        //retrieve the block in case it needs to be fetched from
        //disk or recomputed
        block->retrieve();
        block->data_sort(sort.get());
        //release the block back to disk or wherever it needs to go
        block->release();
    }
}

Tensor::iterator
Tensor::begin() const
{
    return blocks_->begin();
}

Tensor::iterator
Tensor::end() const
{
    return blocks_->end();
}

bool
Tensor::contains(const uli* indices) const
{
    usi nindex = descr_->nindex();
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
    if (storage_type == Tensor::on_disk)
    {
        std::string name = this->name_ + ".tensor";
        disk_buffer_ = new DiskBuffer(name);
    }
    
}

void
Tensor::configure(tensor_erase_policy_t erase_type)
{
    erase_type_ = erase_type;
}

void
Tensor::configure(tensor_distribution_t dist_type)
{
    distribution_type_ = dist_type;
}

void
Tensor::configure(BlockRetrieveAction* block_action)
{
    TensorBlock* block = get_first_nonnull_block();
    if (!block)
    {
        fill_blocks(Tensor::action);
        block = get_first_nonnull_block();
    }

    if (!block->get_retrieve_action())
    {
        TensorRetrieveActionPtr action = new TensorRetrieveAction(this);
        action->add(block_action);
        NodeMap<TensorBlock>::iterator it(blocks_->begin());
        NodeMap<TensorBlock>::iterator stop(blocks_->end());
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
    filler->set_index_descr(descr_.get());
    filler_ = new ThreadedTensorElementComputer(filler);
    storage_type_ = Tensor::recomputed;

    uli idx = 0;
    uli indexset[NINDEX];
    TensorValueEstimater* estimater = filler->get_estimater(depth_);
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (block)
        {
            if (block->is_permutationally_unique())
            {
                raise(SanityCheckError, "cannot re-fill tensor in configure filler");
            }
            else
            {
                block->sync_max_log();
                continue;
            }
        }
        
        blocks_->indices(idx, indexset);

        float maxlog = LOG_NONZERO;
        if (estimater)
            maxlog = estimater->max_log(indexset);
        
        if (maxlog < YetiRuntime::nonnull_cutoff)
            continue;
 
        //if we have gotten here we actually need to create a new node
        block = make_block(indexset);       
        if (!block) //rigorously zero
            continue;

        block->set_max_log(maxlog);  
    }
}

void
Tensor::configure_degeneracy(const PermutationGroupPtr& grp)
{
    PermutationSetPtr pset = grp->get_generator_set();


    uli idx = 0;
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
    uli tmpspace[NINDEX];
    uli indices_used[100];
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (!block) //we don't need to create another node
            continue;
            
        const uli* block_indices = block->get_indices();

        uli idx = blocks_->index(block_indices);
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
            idx = blocks_->index(tmpspace);

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
        descr_->get_nelements_data(nelements);
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
        descr_->get_nelements_data(nelements);
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
    PermutationGroup* grp = new PermutationGroup(descr_->nindex());
    PermutationGroup::iterator it = grp->begin();
    PermutationGroup::iterator stop = grp->end();
    for ( ; it != stop; ++it)
    {
        grp->add(*it);
    }
    grp->set_closed();

    usi nindex = descr_->nindex();
    TensorIndexDescr* newdescr = new TensorIndexDescr(nindex);
    for (usi i=0; i < nindex; ++i)
        newdescr->set(i, descr_->get(i));

    Tensor* newtensor = new Tensor(name_, newdescr, grp);
    newtensor->accumulate(this, 1.0);
    return newtensor;
}

void
Tensor::convert(
    Matrix* matrix,
    MatrixConfiguration* config
)
{
    uli idx = 0;
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (!block) //we don't need to create another node
            continue;

        block->set_read_mode();
        block->retrieve();
        block->convert(matrix, config, descr_.get());
        block->release();
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
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (!block)
        {
            blocks_->indices(idx, indexset);
            block = make_block(indexset);
        }

        if (!block)
            continue; //null

        if (!block->is_permutationally_unique())
            continue;

        block->set_write_mode();
        block->retrieve();
        block->accumulate(matrix, config, descr_.get());
        block->release();
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
    uli nthread = YetiRuntime::nthread();
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
        erase_zero_blocks();
        if (parent_tensor_)
        {
            parent_tensor_->update_to_subtensor();
            parent_tensor_->erase_zero_blocks();
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
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
    uli nthread = YetiRuntime::nthread();
    bool read_mode = false;
    uli indexset[NINDEX];
    for ( ; it != stop; ++it, ++idx)
    {
        //if (INVALID_THREAD_TASK(idx, threadnum, nthread))
        if (threadnum != 0)
            continue;

        TensorBlock* block = *it;
        if (!block) //we don't need to create another node
        {
            blocks_->indices(idx, indexset);
            continue;
        }

        if (!block->is_permutationally_unique())
        {
            block->sync_max_log();
            block->obsolete();
            continue;
        }

        /**
            Depending on the type of element operation,
            you might need different modes - e.g.
        */
        op->set_mode(block);

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
        block->release();
    }

}


void
Tensor::fill(
    const ThreadedTensorElementComputerPtr& filler
)
{
    if (!unsort_perm_->is_identity())
    {
        raise(SanityCheckError, "Cannot fill a tensor not in its original permutation state");
    }

    ThreadGroup* grp = YetiRuntime::get_thread_grp();
    uli nthread = YetiRuntime::nthread();
    for (uli i=0; i < nthread; ++i)
    {
        Thread* thr = new TensorFillThread(i, this, filler);
        grp->add(thr);
    }
    grp->run();
    grp->wait();
    grp->clear();

    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (block && !block->is_permutationally_unique())
        {
            block->sync_max_log();
        }
    }

    erase_zero_blocks();
}

bool
Tensor::equals(
    const void* data
)
{
    uli idx = 0;
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
    
    const void* dataptr = data;
    
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (!block || !block->is_permutationally_unique())
            continue;

        block->set_read_mode();

        block->retrieve();

        //dataptr gets incremented
        bool check = block->equals(&dataptr);
        
        if (!check)
        {
            block->release();
            return false;
        }
        
        block->release();
    }
    
    return true;
}

bool
Tensor::equals(Tensor* tensor)
{
    if (tensor->get_block_map()->size() != blocks_->size())
        return false;

    uli idx = 0;
    NodeMap<TensorBlock>::iterator it_me(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
    NodeMap<TensorBlock>::iterator 
        it_other(tensor->get_block_map()->begin());
    
    for ( ; it_me != stop; ++it_me, ++it_other, ++idx)
    {
        TensorBlock* block_me = *it_me;
        TensorBlock* block_other = *it_other;
        
        bool null_me = block_me ? block_me->get_max_log() < YetiRuntime::nonnull_cutoff : true;
        bool null_other = block_other ?  block_other->get_max_log() < YetiRuntime::nonnull_cutoff : true;

        if (null_me && !null_other)
            return false;
        else if (null_other && !null_me)
            return false;
        else if (null_other && null_me) //equal! both zero!
            continue;

        if (!block_me->is_permutationally_unique() 
            && !block_other->is_permutationally_unique())
            continue;
        
        block_me->set_read_mode();
        block_other->set_read_mode();
        
        block_me->retrieve();
        block_other->retrieve();
        bool check = block_me->equals(block_other);
        block_me->release();
        block_other->release();
        
        if (!check)
            return false;
    }
    
    return true;
}

void
Tensor::fill_blocks(Tensor::tensor_storage_t store_type)
{
    uli idx = 0;
    uli indexset[NINDEX];
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
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

        blocks_->indices(idx, indexset);
        /** Only do permutationally unique blocks! This is very important for threads! */
        if (tensor_grp_->improves_sort(indexset))
            continue;

        //do not initialize the blocks on a fill
        block = make_block(indexset, store_type);
        if (block)
            block->set_max_log(1);
    }
}

void
Tensor::fill(
    uli threadnum,
    const ThreadedTensorElementComputerPtr& new_filler
)
{
    uli idx = 0;
    uli nthread = YetiRuntime::nthread();
    uli indexset[NINDEX];
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
    TensorElementComputer* filler = new_filler->get_computer(threadnum);
    TensorValueEstimater* estimater = filler->get_estimater(depth_);
    for ( ; it != stop; ++it, ++idx)
    {
        if (INVALID_THREAD_TASK(idx, threadnum, nthread))
            continue;

        TensorBlock* block = *it;
        if (block)
        {
            if (block->is_permutationally_unique())
            {
                raise(SanityCheckError, "cannot re-fill tensor");
            }
            else //this block was created by a non-permutationally unique block
            {
                block->sync_max_log();
                continue;    
            }
        } //we don't need to create another node

        blocks_->indices(idx, indexset);
        /** Only do permutationally unique blocks! This is very important for threads! */
        if (tensor_grp_->improves_sort(indexset))
            continue;

        float maxlog = estimater ? estimater->max_log(indexset) : LOG_NONZERO;
        if (maxlog < YetiRuntime::nonnull_cutoff) //this contribution can be skipped
            continue;



        //if we have gotten here we actually need to create
        //a new node
        block = make_block(indexset);
        if (!block) //some blocks can be determined rigorously zero ahead of time
            continue;

        if (!block->is_permutationally_unique())
        {
            raise(SanityCheckError, "should not have reached here in fill");
        }
        block->set_write_mode();
        //retrieve the block in case it needs to be fetched from
        //disk or recomputed
        block->retrieve();
        //fill the block with values
        block->fill(filler);

        //after the fill, we are done with the block
        //so we need to ensure that it is synced with the
        //parent storage area
        block->sync();

        //release the block back to disk or wherever it needs to go
        block->release();
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
        val->set_index_descr(descr_.get());
        fill(new ThreadedTensorElementComputer(val));
    }
    else
    {
        TensorElementComputerPtr filler
            = new MemsetElementComputer(TemplateInfo::double_type);
        filler->set_index_descr(descr_.get());
        fill(new ThreadedTensorElementComputer(filler));
    }
}

TensorBlock*
Tensor::get_block(uli index) const
{
    return blocks_->get(index);
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
        return blocks_->get(permuted_indices);
    }
    else
    {
        return blocks_->get(indices);
    }
}

TensorBlockMap*
Tensor::get_block_map() const
{
    return blocks_;
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
Tensor::get_descr() const
{
    return descr_.get();
}

usi
Tensor::get_depth() const
{
    return depth_;
}

DiskBuffer*
Tensor::get_disk_buffer() const
{
    return disk_buffer_.get();
}

ThreadedTensorElementComputer*
Tensor::get_element_computer() const
{
    return filler_.get();
}

TensorBlock*
Tensor::get_first_nonnull_block() const
{
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (block)
            return block;
    }

    return 0;
}

uli
Tensor::get_index(const uli *indices) const
{
    return blocks_->index(indices);
}

const std::string&
Tensor::get_name() const
{
    return name_;
}

Tensor::tensor_distribution_t
Tensor::get_distribution_type() const
{
    return distribution_type_;
}

Tensor::tensor_erase_policy_t
Tensor::get_erase_type() const
{
    return erase_type_;
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

Permutation*
Tensor::get_sort_permutation() const
{
    return sort_perm_;
}

Permutation*
Tensor::get_unsort_permutation() const
{
    return unsort_perm_;
}

ulli
Tensor::get_totalsize() const
{
    return descr_->totalsize();
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

TensorBlock*
Tensor::get_make_block(uli index)
{
    if (index >= blocks_->size())
    {
        std::string err = stream_printf(
            "invalid block index %d passed to tensor %s",
            index, name_.c_str()
        );
        raise(SanityCheckError, err);
    }

    TensorBlock* block = blocks_->get(index);
    if (block)
        return block;

    uli indices[NINDEX];
    //doesn't exist yet... build it        
    blocks_->indices(index, indices);
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
        blocks_->sizes(),
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

            bool test = dst_block->set_accumulate_mode();
            if (!test)
            {
                cerr << stream_printf("set accumulate failed %d %s", __LINE__, __FILE__) << endl;
                abort();
            }
            src_block->set_read_mode();

            src_block->retrieve();
            dst_block->retrieve();
            src_block->internal_contraction(
                dst_block,
                config
            );
            src_block->release();
            dst_block->release();
        }
    }

    dst_tensor->update();
}

bool
Tensor::nonzero() const
{
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
    uli idx = 0;
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (block && block->get_max_log() > YetiRuntime::nonnull_cutoff)
        {
            return true;
        }
    }
    return false;
}

bool
Tensor::unique_nonzero() const
{
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
    uli idx = 0;
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (block
           && block->get_max_log() > YetiRuntime::nonnull_cutoff
           &&  block->is_permutationally_unique()
                //this->is_parent_block(block)
        )
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
        return true;

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
    return descr_->totalsize();
}

uli
Tensor::nblocks_retrieved() const
{
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
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
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
    uli idx = 0;
    uli indexset[NINDEX];
    ulli n = 0;
    for ( ; it != stop; ++it, ++idx)
    {
        blocks_->indices(idx, indexset);
        if (tensor_grp_->improves_sort(indexset))
            continue;

        n += descr_->get_nelements_data(depth_, indexset);
    }
    return n;
}

void
Tensor::print(std::ostream& os)
{
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());

    os << "Tensor " << name_ << endl;
    os << tensor_grp_ << endl;
    uli idx = 0;
    uli indexset[NINDEX];
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;

        if (block->get_max_log() < YetiRuntime::print_cutoff)
            continue;

        blocks_->indices(idx, indexset);
        os << endl << "Tensor Block "
            << ClassOutput<const uli*>::str(descr_->nindex(), indexset)
            << endl;
        os << "irreps:";
        for (usi i=0; i < descr_->nindex(); ++i)
            os << " " << descr_->get(i)->irrep(indexset[i]);
        os << endl;
        os << "Fetch permutation:  ";
        block->get_fetch_permutation()->print(os); os << endl;
        os << "Resort permutation: ";
        block->get_resort_permutation()->print(os); os << endl;
        ++Env::indent;
        block->set_read_mode();
        block->retrieve();
        block->print(os);
        block->release();
        --Env::indent;
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
Tensor::set_priority(Tensor::tensor_priority_t priority)
{
    priority_ = priority;
}

void
Tensor::set_read_mode()
{
    foreach_nonnull(block, blocks_, TensorBlock,
        block->set_read_mode();
    )
}

void
Tensor::set_write_mode()
{
    foreach_nonnull(block, blocks_, TensorBlock,
        block->set_write_mode();
    )
}

void
Tensor::set_accumulate_mode()
{
    foreach_nonnull(block, blocks_, TensorBlock,
        block->set_accumulate_mode();
    )
}

void
Tensor::sync_to_subtensor()
{
    if (parent_tensor_)
    {
        cerr << "should not be called on subtensor" << endl;
        abort();
    }

    iterator it = blocks_->begin();
    iterator stop = blocks_->end();
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;

        //this is parent of parent tensor but maps onto a block
        //that was changed in the course of some operation
        if (!block->is_subblock() && !block->is_permutationally_unique())
        {
            block->sync_max_log();
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

    iterator it = blocks_->begin();
    iterator stop = blocks_->end();
    for ( ; it != stop; ++it)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;

        //this is parent of parent tensor but maps onto a block
        //that was changed in the course of some operation
        if (!block->is_subblock() && !block->is_permutationally_unique())
        {
            block->sync_max_log();
        }
    }

    erase_zero_blocks();
}


void
Tensor::sync()
{
    if (parent_tensor_)
        parent_tensor_->sync_to_subtensor();

    in_sync = true;

    foreach_nonnull(block, blocks_, TensorBlock,
        block->sync();
    )

    in_sync = false;
}


void
Tensor::reset_degeneracy()
{
    foreach_nonnull(block, blocks_, TensorBlock,
        block->reset_degeneracy();
    )
}

void
Tensor::update(uli threadnum)
{
    uli idx = 0;
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
    uli nthread = YetiRuntime::nthread();

    for ( ; it != stop; ++it, ++idx)
    {
        if (INVALID_THREAD_TASK(idx, threadnum, nthread))
            continue;

        TensorBlock* block = *it;
        if (!block)
            continue;

        if (block->is_permutationally_unique())
        {
            block->set_read_mode();
            block->retrieve();
            block->update();
            block->release();
        }
    }
}

void
Tensor::update()
{    
    ThreadGroup* grp = YetiRuntime::get_thread_grp();
    uli nthread = YetiRuntime::nthread();
    for (uli i=0; i < nthread; ++i)
    {
        Thread* thr = new TensorUpdateThread(i, this);
        grp->add(thr);
    }
    grp->run();
    grp->wait();
    grp->clear();

    uli idx = 0;
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());
    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;

        if (!block->is_permutationally_unique())
        {
            block->sync_max_log();
            block->reset_degeneracy();
        }

#if YETI_SANITY_CHECK
            Permutation* ptest = block->get_fetch_permutation()->product(block->get_resort_permutation());
            if (sort_perm_->is_identity() && !ptest->is_identity())
            {
                cerr << "nonunique subblock does not have sort realigned to parent block" << endl;
                abort();
            }
#endif
    }


    erase_zero_blocks();
    if (parent_tensor_)
        parent_tensor_->update_to_subtensor();
}


void
Tensor::erase_zero_blocks()
{
    return;
    if (erase_type_ == Tensor::keep_zero_blocks)
        return;

    uli idx = 0;
    NodeMap<TensorBlock>::iterator it(blocks_->begin());
    NodeMap<TensorBlock>::iterator stop(blocks_->end());

    for ( ; it != stop; ++it, ++idx)
    {
        TensorBlock* block = *it;
        if (!block)
            continue;

        try{
            erase_if_null(block, idx);
        }
        catch(int e)
        {
            cerr << stream_printf("Caught at %d %s", __LINE__, __FILE__) << endl;
            abort();
        }

    }
}

void
Tensor::erase_if_null(
    TensorBlock* block,
    uli idx
)
{
    if (block->get_max_log() < YetiRuntime::nonnull_cutoff)
    {
        blocks_->erase(idx);
    }
}


bool
Tensor::is_cache_coherent() const
{
    foreach_nonnull(block, blocks_, TensorBlock,
        if (!block->is_cache_coherent())
            return false;
    )
    return true;
}


