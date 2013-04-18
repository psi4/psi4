#include "contraction.h"
#include "matrix.h"
#include "env.h"
#include "permutation.h"
#include "class.h"
#include "tensor.h"
#include "tensorblock.h"
#include "node.h"
#include "index.h"
#include "exception.h"
#include "contractionimpl.h"
#include "sortimpl.h"
#include "dataimpl.h"
#include "malloc.h"
#include "runtime.h"
#include "threadimpl.h"
#include <libsmartptr/printstream.h>
#include <algorithm>

using namespace yeti;
using namespace std;

#define DEBUG_CXN 0

#define PRINT_CXN_TYPE 1

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

DECLARE_MALLOC(Contraction);
DECLARE_MALLOC(ContractionTask);


ContractionTask::ContractionTask(
    uli row,
    uli col,
    uli link,
    Contraction* cxn
) :
    main_block_(0),
    row_(row),
    col_(col),
    link_(link),
    cxn_(cxn),
    follow_permutations_(true)
{
}

ContractionTask::ContractionTask(
    uli row,
    uli col,
    uli link,
    TensorBlock* main_block,
    Contraction* cxn
) :
    main_block_(main_block),
    row_(row),
    col_(col),
    link_(link),
    cxn_(cxn),
    follow_permutations_(false)
{
    main_block_->set_task_owner(true);
}

ContractionTask::~ContractionTask()
{
}

void
ContractionTask::print(ostream& os) const
{
    os << stream_printf("ContractionTask\n");
}

void
ContractionTask::prefetch(uli threadnum)
{
    switch(cxn_->order())
    {
        case Contraction::p_l_r_cxn: prefetch_p_l_r(threadnum); break;
        case Contraction::p_r_l_cxn: prefetch_p_r_l(threadnum); break;
        case Contraction::r_l_p_cxn: prefetch_r_l_p(threadnum); break;
        case Contraction::r_p_l_cxn: prefetch_r_p_l(threadnum); break;
        case Contraction::l_p_r_cxn: prefetch_l_p_r(threadnum); break;
        case Contraction::l_r_p_cxn: prefetch_l_r_p(threadnum); break;
    }
}

void
ContractionTask::prefetch_left_driven(uli threadnum)
{
    uli col = 0;
    if (!main_block_) /** Check to see if we need to be done */
    {
        uli indices[NINDEX];
        uli lidx = cxn_->get_left_index(row_, link_);
        cxn_->ltensor()->get_block_map()->indices(lidx, indices);
        if (cxn_->ltensor()->get_tensor_grp()->improves_sort(indices)) //this block is not to be done
            return;
        
        main_block_ = get_left_block(threadnum, row_, link_);
    }

    TensorBlock* lblock = main_block_;
    if (!lblock) //no work to be done here
        return;

    TensorBlock* rblock = get_right_block(threadnum, link_, col);
    TensorBlock* pblock = get_product_block(threadnum, row_, col);
    if (rblock && pblock && cxn_->do_task(lblock, rblock, pblock))
    {
        rblock->prefetch_read();
    }
    lblock->set_task_owner(true);
    lblock->prefetch_read();
}

void
ContractionTask::prefetch_l_r_p(uli threadnum)
{
    prefetch_left_driven(threadnum);
}

void
ContractionTask::prefetch_l_p_r(uli threadnum)
{
    prefetch_left_driven(threadnum);
}

void
ContractionTask::prefetch_right_driven(uli threadnum)
{
    uli row = 0;
    if (!main_block_) /** Check to see if we need to be done */
    {
        uli indices[NINDEX];
        uli ridx = cxn_->get_right_index(link_, col_);
        cxn_->rtensor()->get_block_map()->indices(ridx, indices);
        if (cxn_->rtensor()->get_tensor_grp()->improves_sort(indices)) //this block is not to be done
            return;
        
        main_block_ = get_right_block(threadnum, link_, col_);
    }

    TensorBlock* rblock = main_block_;
    if (!rblock) //no work to be done here
        return;

    TensorBlock* lblock = get_left_block(threadnum, row_, link_);
    TensorBlock* pblock = get_product_block(threadnum, row_, col_);
    if (pblock && lblock && cxn_->do_task(lblock, rblock, pblock))
    {
        lblock->prefetch_read();
    }
    rblock->set_task_owner(true);
    rblock->prefetch_read();
}

void
ContractionTask::prefetch_r_l_p(uli threadnum)
{
    prefetch_right_driven(threadnum);
}

void
ContractionTask::prefetch_r_p_l(uli threadnum)
{
    prefetch_right_driven(threadnum);
}

void
ContractionTask::prefetch_product_driven(uli threadnum)
{
    uli link = 0;
    TensorBlock* pblock = main_block_ = get_product_block(threadnum, row_, col_);
    if (!pblock) //no work to be done here
        return;

    TensorBlock* lblock = get_left_block(threadnum, row_, link);
    TensorBlock* rblock = get_right_block(threadnum, link, col_);
    if (rblock && lblock && cxn_->do_task(lblock, rblock, pblock))
    {
        rblock->prefetch_read();
        lblock->prefetch_read();
    }
}

void
ContractionTask::prefetch_p_l_r(uli threadnum)
{
    prefetch_product_driven(threadnum);
}

void
ContractionTask::prefetch_p_r_l(uli threadnum)
{
    prefetch_product_driven(threadnum);
}

void
ContractionTask::run(uli threadnum)
{
    if (!main_block_) //nothing to do
        return;

    switch(cxn_->order())
    {
        case Contraction::p_l_r_cxn: run_p_l_r(threadnum); break;
        case Contraction::p_r_l_cxn: run_p_r_l(threadnum); break;
        case Contraction::r_l_p_cxn: run_r_l_p(threadnum); break;
        case Contraction::r_p_l_cxn: run_r_p_l(threadnum); break;
        case Contraction::l_p_r_cxn: run_l_p_r(threadnum); break;
        case Contraction::l_r_p_cxn: run_l_r_p(threadnum); break;
    }
}


uli
ContractionTask::append_info(uli* data) const
{
    yeti_throw(SanityCheckError, "Dynamic load balancing not yet available");
    return 0;
}

void
ContractionTask::run(
    uli threadnum,
    TensorBlock* lblock,
    TensorBlock* rblock,
    TensorBlock* pblock
)
{
#if DEBUG_CXN
    cout << stream_printf("Accumulating contraction task %s : %s = %s * %s on node %d on thread %d [%d]=[%d][%d]\n",
    pblock->get_parent_tensor()->get_name().c_str(),
    pblock->get_index_string().c_str(),
    lblock->get_index_string().c_str(),
    rblock->get_index_string().c_str(),
    YetiRuntime::me(),
    threadnum,
    pblock->get_node_owner(),
    lblock->get_node_owner(),
    rblock->get_node_owner());
    Env::outn().flush();
#endif


    uli threadtest = YetiRuntime::get_thread_number();
    if (threadtest != threadnum)
    {
        YetiRuntime::lock_print();
        cerr << "thread number " << threadtest << " returned is wrong for thread " << threadnum << endl;
        YetiRuntime::unlock_print();
        abort();
    }

#if YETI_SANITY_CHECK
    if (!lblock->has_fast_read() && !lblock->is_retrieved())
        yeti_throw(SanityCheckError, "Left block is not retrieved in cxn");
    if (!rblock->has_fast_read() && !rblock->is_retrieved())
        yeti_throw(SanityCheckError, "Right block is not retrieved in cxn");
    if (!pblock->is_retrieved())
        yeti_throw(SanityCheckError, "Product block is not retrieved in cxn");
#endif

    pblock->accumulate(lblock, rblock, cxn_);
}

Contraction::tensor_position_t
get_tensor_position(
    Tensor* ltensor,
    Tensor* rtensor,
    Tensor* ptensor,
    Tensor* tensor
)
{
    if (tensor == ltensor)
        return Contraction::left_tensor;
    else if (tensor == rtensor)
        return Contraction::right_tensor;
    else
        return Contraction::product_tensor;
}

Contraction::Contraction(
    double scale,
    Tensor* ltensor,
    Tensor* rtensor,
    Tensor* ptensor,
    const MatrixIndexPtr& lindex,
    const MatrixIndexPtr& rindex,
    Permutation* lperm,
    Permutation* rperm
) :
    contraction_grp_(0),
    default_grp_(0),
    final_grp_(0),
    lconfig_(new MatrixConfiguration(lindex, ltensor->get_tensor_grp())),
    rconfig_(new MatrixConfiguration(rindex, rtensor->get_tensor_grp())),
    pconfig_(0),
    lindex_(lindex),
    rindex_(rindex),
    pindex_(0),
    nrows_(0),
    ncols_(0),
    nlink_(0),
    row_(0),
    col_(0),
    link_(0),
    ltensor_(ltensor),
    rtensor_(rtensor),
    ptensor_(ptensor),
    alpha_type_(Contraction::left_tensor),
    engine_(0),
    scale_(scale),
    cxn_configs_(new ContractionConfiguration*[YetiRuntime::nthread()]),
    use_replicated_product_(false),
    product_thread_clash_(false),
    product_node_clash_(false),
    all_tensors_replicated_(false)
{

    if (ltensor_->get_priority() == Tensor::override_priority)
        yeti_throw(SanityCheckError,
            "left tensor in contraction cannot be given override priority");
    if (rtensor_->get_priority() == Tensor::override_priority)
        yeti_throw(SanityCheckError,
            "right tensor in contraction cannot be given override priority");
    if (ptensor_->get_priority() == Tensor::override_priority)
        yeti_throw(SanityCheckError,
            "product tensor in contraction cannot be given override priority");

    if (ptensor_->get_subtensor_count() > 1)
        yeti_throw(SanityCheckError,
            "product tensor in contraction cannot have more than one subtensor");


    YetiRuntime::start_timer("cxn constructor");
    ltensor->set_tensor_position(Contraction::left_tensor);
    rtensor->set_tensor_position(Contraction::right_tensor);
    ptensor->set_tensor_position(Contraction::product_tensor);

    //if rows are the back indices, matrix is "transposed"
    transpose_left_ = lindex_->is_transpose();
    transpose_right_ = rindex_->is_transpose();
    
    usi nrow_idx = transpose_left_ ? lindex->ncolindex() : lindex->nrowindex();
    usi ncol_idx = transpose_right_ ? rindex->nrowindex() : rindex->ncolindex();
    usi ncxn_idx = transpose_left_ ? lindex->nrowindex() : lindex->ncolindex();
    usi ntarget_idx = nrow_idx + ncol_idx;
    
    if (ntarget_idx == 0) //dot product
        init_dot_product();
    else if (ncxn_idx == 0) //direct product
        init_direct_product();
    else
        init_cxn();

    //configure the priorities for the contraction
    Tensor* null_tensor = 0;
    vector<Tensor*> tensors(3, null_tensor);
    tensors[0] = ltensor_;
    tensors[1] = rtensor_;
    tensors[2] = ptensor_;

    /** the tensors can have a size weighting in the sort */
    ltensor_->set_sort_size_weight(1);
    rtensor_->set_sort_size_weight(1);
    ptensor_->set_sort_size_weight(3);

    if (ptensor_->is_open_subtensor()) //this MUST be the alpha tensor
        ptensor_->set_priority(Tensor::override_priority);

    std::sort(tensors.begin(), tensors.end(), tensor_less);

    alpha_tensor_ = tensors[2];
    beta_tensor_ = tensors[1];
    gamma_tensor_ = tensors[0];

    if      (alpha_tensor_->is_replicated() && YetiRuntime::nproc() > 1)
    {
        all_tensors_replicated_ = true;
        //if the product tensor has enough blocks, distribute based on that one
        uli nblocks = ptensor_->get_block_map()->size();
        if (ptensor_->get_block_map()->size() > YetiRuntime::nproc())
        {
            distribution_type_ = product_tensor;
        }
        else if (ltensor_->get_block_map()->size() > YetiRuntime::nproc())
        {
            distribution_type_ = right_tensor;
        }
        else if (rtensor_->get_block_map()->size() > YetiRuntime::nproc())
        {
            distribution_type_ = left_tensor;
        }
        else //no tensor is big enough to distribute tasks on...
            //as long as the product tensor isn't a dot product... use it
        {
            if (ntarget_idx > 0)
                distribution_type_ = product_tensor;
            else
                distribution_type_ = left_tensor;
        }
    }
    else if (alpha_tensor_ == ptensor_)
    {
        distribution_type_ = product_tensor;
    }
    else if (alpha_tensor_ == ltensor_)
    {
        distribution_type_ = left_tensor;
    }
    else if (alpha_tensor_ == rtensor_)
    {
        distribution_type_ = right_tensor;
    }

    for (uli i=0; i < YetiRuntime::nthread(); ++i)
    {
        cxn_configs_[i] = new ContractionConfiguration(lconfig_, rconfig_, pconfig_, this);
    }

    bool l_ext_local, l_cxn_local, r_ext_local, r_cxn_local, p_rows_local, p_cols_local;
    cxn_configs_[0]->check_distribution(
        ltensor, rtensor, ptensor,
        l_ext_local, l_cxn_local, 
        r_cxn_local, r_ext_local,
        p_rows_local, p_cols_local
    );

    bool col_external_index_distribution = 
         l_cxn_local && ltensor->is_distributed() 
      && r_cxn_local && rtensor->is_distributed()
      && ptensor->is_distributed() 
      && p_rows_local;

    bool replicated_col_external_index_distribution =
         ltensor->is_replicated()
      && r_cxn_local && rtensor->is_distributed()
      && ptensor->is_distributed() 
      && p_rows_local;

    bool row_external_index_distribution = 
         l_cxn_local && ltensor->is_distributed() 
      && r_cxn_local && rtensor->is_distributed()
      && ptensor->is_distributed() 
      && p_cols_local;

    bool replicated_row_external_index_distribution =
          l_cxn_local && ltensor->is_distributed()
       && rtensor->is_replicated()
       && ptensor->is_distributed()
       && p_cols_local;
#if 0
    Env::out0() << stream_printf(
        "left locality:  ext=%d cxn=%d\n"
        "right locality: ext=%d cxn=%d\n"
        "prod locality:  row=%d col=%d",
        l_ext_local, l_cxn_local, 
        r_ext_local, r_cxn_local,
        p_rows_local, p_cols_local
    ) << endl;
#endif

    //no dot products!
    bool cxn_index_distribution = 
        YetiRuntime::nproc() > 1
        && l_ext_local && !ltensor_->is_recomputed()
        && r_ext_local && !rtensor_->is_recomputed()
        && ntarget_idx > 0; //no dot product

    bool minimal_recompute_col_distribution =
           rtensor->is_recomputed()
        && p_rows_local && ptensor->is_distributed();

    bool minimal_recompute_row_distribution =
           ltensor->is_recomputed()
        && p_cols_local && ptensor->is_distributed();
        


    bool alpha_tensor_downgradable = false;
    bool beta_tensor_upgradable = false;
    bool gamma_tensor_upgradable = false;
    if (row_external_index_distribution)
    {
        alpha_tensor_ = rtensor_;
        beta_tensor_ = ptensor_;
        gamma_tensor_ = ltensor_;
        distribution_type_ = Contraction::product_tensor;
    }
    else if (replicated_row_external_index_distribution)
    {
        alpha_tensor_ = ptensor_;
        beta_tensor_ = ltensor_;
        gamma_tensor_ = rtensor_;
        distribution_type_ = Contraction::product_tensor;
    }
    else if (col_external_index_distribution)
    {
        alpha_tensor_ = ltensor_;
        beta_tensor_ = ptensor_;
        gamma_tensor_ = rtensor_;
        distribution_type_ = Contraction::product_tensor;
    }
    else if (replicated_col_external_index_distribution)
    {
        alpha_tensor_ = ptensor_;
        beta_tensor_ = rtensor_;
        gamma_tensor_ = ltensor_;
        distribution_type_ = Contraction::product_tensor;
    }
    else if (minimal_recompute_col_distribution)
    {
        alpha_tensor_ = rtensor_;
        beta_tensor_ = ptensor_;
        gamma_tensor_ = ltensor_;
        distribution_type_ = Contraction::product_tensor;
    }
    else if (minimal_recompute_row_distribution)
    {
        alpha_tensor_ = ltensor_;
        beta_tensor_ = ptensor_;
        gamma_tensor_ = rtensor_;
        distribution_type_ = Contraction::product_tensor;
    }
    else if (cxn_index_distribution)
    {
        alpha_tensor_ = ptensor_;
        beta_tensor_ = rtensor_;
        gamma_tensor_ = ltensor_;
        distribution_type_ = Contraction::product_tensor;
    }
    else
    {
        alpha_tensor_downgradable
            = alpha_tensor_->is_distributed() && alpha_tensor_->get_storage_type() == Tensor::in_core;
        beta_tensor_upgradable
            = beta_tensor_->get_storage_type() != Tensor::in_core || beta_tensor_->is_distributed();
        gamma_tensor_upgradable
            = gamma_tensor_->get_storage_type() != Tensor::in_core || gamma_tensor_->is_distributed();

        if (alpha_tensor_downgradable && beta_tensor_upgradable)
        {
            Tensor* tmp = alpha_tensor_;
            alpha_tensor_ = beta_tensor_;
            beta_tensor_ = tmp;
            if (gamma_tensor_upgradable)
            {
                tmp = gamma_tensor_;
                gamma_tensor_ = beta_tensor_;
                beta_tensor_ = tmp;
            }
        }
    }




    if (alpha_tensor_ == ltensor_)
        alpha_type_ = left_tensor;
    else if (alpha_tensor_ == rtensor_)
        alpha_type_ = right_tensor;
    else if (alpha_tensor_ == ptensor_)
        alpha_type_ = product_tensor;
    else
        yeti_throw(SanityCheckError, "no tensor is the alpha tensor");

    if (beta_tensor_ == ltensor_)
        this->beta_type_ = left_tensor;
    else if (beta_tensor_ == rtensor_)
        this->beta_type_ = right_tensor;
    else if (beta_tensor_ == ptensor_)
        this->beta_type_ = product_tensor;
    else
        yeti_throw(SanityCheckError, "no tensor is the beta tensor");

    if (gamma_tensor_ == ltensor_)
        this->gamma_type_ = left_tensor;
    else if (gamma_tensor_ == rtensor_)
        this->gamma_type_ = right_tensor;
    else if (gamma_tensor_ == ptensor_)
        this->gamma_type_ = product_tensor;
    else
        yeti_throw(SanityCheckError, "no tensor is the gamma tensor");

    alpha_tensor_->set_priority(Tensor::alpha_tensor);
    beta_tensor_->set_priority(Tensor::beta_tensor);
    gamma_tensor_->set_priority(Tensor::gamma_tensor);


    if (transpose_left_ && transpose_right_)
    {
#if PRINT_CXN_TYPE
        dout << "Contraction type T-T" << endl;
#endif
        engine_ = new ContractionTemplate<Contraction_tt>;
    }
    else if (transpose_left_)
    {
#if PRINT_CXN_TYPE
        dout << "Contraction type T-N" << endl;
#endif
        engine_ = new ContractionTemplate<Contraction_tn>;
    }
    else if (transpose_right_)
    {
#if PRINT_CXN_TYPE
        dout << "Contraction type N-T" << endl;
#endif
        engine_ = new ContractionTemplate<Contraction_nt>;
    }
    else
    {
#if PRINT_CXN_TYPE
        dout << "Contraction type N-N" << endl;
#endif
        engine_ = new ContractionTemplate<Contraction_nn>;
    }

    //we have to make sure all of the matrices are at the same depth
    //determine the maximum recursion depth
    usi ldepth = ltensor_->get_depth();
    usi rdepth = rtensor_->get_depth();
    usi pdepth = ptensor_->get_depth();

    if (ldepth != rdepth || rdepth != pdepth)
        yeti_throw(SanityCheckError, "tensors not aligned for multiplication");

    ltensor_->configure_degeneracy(lconfig_->get_cxn_permutation_grp());
    if (ltensor_ != rtensor_)
        rtensor_->configure_degeneracy(rconfig_->get_cxn_permutation_grp());

    uli nthread = YetiRuntime::nthread();
    uli nproc = YetiRuntime::nproc();
    bool thread_replicate = nthread > 1 && ptensor_->get_block_map()->size() < 1000;
    bool node_replicate = nproc > 1 && ptensor_->is_replicated();
    use_replicated_product_ = thread_replicate || node_replicate;
    /** we might have a thread clash, but if the tensor is replicated we will build an entire tmp tensor */
    product_thread_clash_ = nthread > 1 && alpha_type_ != product_tensor && ptensor_->is_distributed();
    
    product_node_clash_ = ptensor->is_distributed() && YetiRuntime::nproc() > 1 && distribution_type_ != product_tensor;
    bool need_tmp_product_block = product_thread_clash_ || product_node_clash_;
    for (uli i=0; i < nthread; ++i)
    {
        if (use_replicated_product_)
            cxn_configs_[i]->clone_product_tensor_for_threads(ptensor_);

        if (need_tmp_product_block && !use_replicated_product_)
        {
            TensorBlock* tmpblock = new TensorBlock(ptensor_);
            cxn_configs_[i]->set_tmp_accumulate_block(tmpblock);
        }
    }


#if 1//PRINT_CXN_DETAILS
    Env::out0() << endl
                << "Product:        " << ptensor_->get_name() << endl
                << "Left:           " << ltensor_->get_name() << endl
                << "Right:          " << rtensor_->get_name() << endl
                << "Alpha:          " << alpha_tensor_->get_name() << endl
                << "Beta:           " << beta_tensor_->get_name() << endl
                << "Gamma:          " << gamma_tensor_->get_name() << endl
                << "Distr:          " << distribution_type_ << endl
                << "Repl Prod:      " << use_replicated_product_ << endl
                << "Tmp Prod:       " << need_tmp_product_block << endl
                << "Thread Clash:   " << product_thread_clash_ << endl
                << "Node Clash:     " << product_node_clash_ << endl
                << "Replicated:     " << all_tensors_replicated_ << endl
                << "Open Subtensor: " << ptensor_->is_open_subtensor() << endl
                << "Row Ext Distr:  " << row_external_index_distribution << endl
                << "Re Row Ext:     " << minimal_recompute_row_distribution << endl
                << "Col Ext Distr:  " << col_external_index_distribution << endl
                << "Re Col Ext:     " << minimal_recompute_col_distribution << endl
                << "Cxn Distr:      " << cxn_index_distribution << endl
                << "Alpha Down:     " << alpha_tensor_downgradable << endl
                << "Beta Up:        " << beta_tensor_upgradable << endl
                << "Gamma Up:       " << gamma_tensor_upgradable << endl
                << endl;
#endif

    GlobalTileQueue::add(this);

    uli no_offset = 0;
    uli no_max = 0;

    for (uli i=0; i < YetiRuntime::nthread(); ++i)
    {
        cxn_configs_[i]->configure_left_block(ltensor_);
        cxn_configs_[i]->configure_right_block(rtensor_);
    }

    ContractionConfiguration* cxn_config = cxn_configs_[0];
    nrows_ = cxn_config->ncxn_rows_left();
    ncols_ = cxn_config->ncxn_cols_right();
    nlink_ = cxn_config->ncxn_cols_left();
#if YETI_SANITY_CHECK
    uli check = cxn_config->ncxn_rows_right();
    if (check != nlink_)
        yeti_throw(SanityCheckError, "matrices not aligned for multiplication");
#endif

    if  (alpha_type_ == left_tensor)
    {
        if (beta_type_ == right_tensor)
            order_ = Contraction::l_r_p_cxn; 
        else
            order_ = Contraction::l_p_r_cxn;
    }
    else if (alpha_type_ == right_tensor)
    {
        if (beta_type_ == left_tensor)
            order_ = Contraction::r_l_p_cxn;
        else
            order_ = Contraction::r_p_l_cxn;
    }
    else if (alpha_type_ == product_tensor)
    {
        if (beta_type_ == left_tensor)
            order_ = Contraction::p_l_r_cxn;
        else
            order_ = Contraction::p_r_l_cxn;
    }
    else
        yeti_throw(SanityCheckError, "Contraction is clearly misconfigured");

#if PRINT_CXN_TYPE
    switch(order_)
    {
        case p_l_r_cxn: cout << "p_l_r" << endl; break;
        case p_r_l_cxn: cout << "p_r_l" << endl; break;
        case r_l_p_cxn: cout << "r_l_p" << endl; break;
        case r_p_l_cxn: cout << "r_p_l" << endl; break;
        case l_p_r_cxn: cout << "l_p_r" << endl; break;
        case l_r_p_cxn: cout << "l_r_p" << endl; break;
    }
#endif

    YetiRuntime::stop_timer("cxn constructor");
}

Contraction::~Contraction()
{
    delete engine_;
    uli nthread = YetiRuntime::nthread();
    for (uli i=0; i < nthread; ++i)
    {
        delete cxn_configs_[i];
    }
    delete[] cxn_configs_;
}

void
Contraction::add_dynamic_task(
    ContractionTask* task
)
{
    yeti_throw(SanityCheckError, "No dynamic load balancing currently allowed");
}

Task*
Contraction::get_next_task()
{
    switch(order_)
    {
        case p_l_r_cxn: return get_next_p_l_r(); break;
        case p_r_l_cxn: return get_next_p_r_l(); break;
        case r_l_p_cxn: return get_next_r_l_p(); break;
        case r_p_l_cxn: return get_next_r_p_l(); break;
        case l_p_r_cxn: return get_next_l_p_r(); break;
        case l_r_p_cxn: return get_next_l_r_p(); break;
    }
    return 0;
}

bool
Contraction::do_task(
    TensorBlock* lblock,
    TensorBlock* rblock,
    TensorBlock* pblock
)
{
    uli proc_number = YetiRuntime::me();
    if      (YetiRuntime::nproc() == 1)
    {
        return true;
    }
    else if (all_tensors_replicated_)
    {
        if (distribution_type_ == product_tensor)
        {
            uli tasknum = pblock->get_block_number() % YetiRuntime::nproc();
            return tasknum == proc_number;
        }
        else if (distribution_type_ == right_tensor)
        {
            uli tasknum = rblock->get_block_number() % YetiRuntime::nproc();
            return tasknum == proc_number;
        }
        else if (distribution_type_ == left_tensor)
        {
            uli tasknum = lblock->get_block_number() % YetiRuntime::nproc();
            return tasknum == proc_number;
        }
        else
        {
            yeti_throw(SanityCheckError, "All replicated tensor cxn has no configuration. WTF mate?");
        }
    }
    else if (distribution_type_ == product_tensor)
    {
        return pblock->get_node_number() == proc_number;
    }
    else if (distribution_type_ == left_tensor)
    {
        return lblock->get_node_number() == proc_number;
    }
    else if (distribution_type_ == right_tensor)
    {
        return rblock->get_node_number() == proc_number;
    }
    else
    {
        yeti_throw(SanityCheckError, "invalid contraction distribution type");
    }
    return false;
}

Contraction::tensor_position_t
Contraction::get_alpha_type() const
{
    return alpha_type_;
}

void
Contraction::init_cxn_grp()
{
    lcxn_grp_ = transpose_left_ ?
        lconfig_->get_row_permutation_grp() :
        lconfig_->get_col_permutation_grp();
    rcxn_grp_ = transpose_right_ ?
        rconfig_->get_col_permutation_grp() :
        rconfig_->get_row_permutation_grp();
        
    contraction_grp_ = lcxn_grp_->intersection_grp(rcxn_grp_);


    if (transpose_left_)
    {
        PermutationGroupPtr lcxn_grp_full 
            = lindex_->expand_on_rows(contraction_grp_);
            
        lconfig_->configure_symmetrization_set(lcxn_grp_full);
       
        bool cxn_on_cols = false;
        lconfig_->configure_isotropy(lcxn_grp_full, cxn_on_cols);
    }
    else
    {
        PermutationGroupPtr lcxn_grp_full 
            = lindex_->expand_on_cols(contraction_grp_);
            
        lconfig_->configure_symmetrization_set(lcxn_grp_full);
        
        bool cxn_on_cols = true;
        lconfig_->configure_isotropy(lcxn_grp_full, cxn_on_cols);
    }

    PermutationGroupPtr rcxn_grp_full;
    if (transpose_right_)
    {
        PermutationGroupPtr rcxn_grp_full 
            = rindex_->expand_on_cols(contraction_grp_);
        
        rconfig_->configure_symmetrization_set(rcxn_grp_full);
        
        bool cxn_on_cols = true;
        rconfig_->configure_isotropy(rcxn_grp_full, cxn_on_cols);
    }
    else
    {
        PermutationGroupPtr rcxn_grp_full 
            = rindex_->expand_on_rows(contraction_grp_);

        rconfig_->configure_symmetrization_set(rcxn_grp_full);

        bool cxn_on_cols = false;
        rconfig_->configure_isotropy(rcxn_grp_full, cxn_on_cols);
    }
}

void
Contraction::init_cxn()
{
    init_cxn_grp();

    usi nrows_product = transpose_left_ ? lconfig_->get_index()->ncolindex() :
                                    lconfig_->get_index()->nrowindex();
    usi ncols_product = transpose_right_ ? rconfig_->get_index()->nrowindex() :
                                    rconfig_->get_index()->ncolindex();
    usi nprod_idx = nrows_product + ncols_product;


    default_grp_ = new PermutationGroup(nprod_idx);

    default_grp_->add_and_expand(
        transpose_left_ ?
            lconfig_->get_col_permutation_grp() :
            lconfig_->get_row_permutation_grp(),
        Permutation::ExpandForward
    );

    default_grp_->add_and_expand(
        transpose_right_ ?
            rconfig_->get_row_permutation_grp() :
            rconfig_->get_col_permutation_grp(),
        Permutation::ExpandBackward
    );
    default_grp_->close();

    final_grp_ = ptensor_->get_tensor_grp();

    pindex_ = new MatrixIndex(nrows_product, ncols_product);
    pconfig_ = new MatrixConfiguration(pindex_, final_grp_);
    pconfig_->configure_symmetrization_set(default_grp_);
    pconfig_->configure_matrix_grp(final_grp_);
}

void
Contraction::init_dot_product()
{
    init_cxn_grp();

    usi nprod_idx = lconfig_->get_index()->ncolindex()
                    + rconfig_->get_index()->ncolindex();

    //no permutational symmetry for a single element
    final_grp_ = ptensor_->get_tensor_grp();
    default_grp_ = final_grp_;

    //create the matrix index for the product configuration
    pindex_ = new MatrixIndex(0,1);
    pconfig_ = new MatrixConfiguration(pindex_, final_grp_);
    pconfig_->configure_matrix_grp(final_grp_);
}

void
Contraction::init_direct_product()
{
    lcxn_grp_ = ltensor_->get_tensor_grp();
    rcxn_grp_ = rtensor_->get_tensor_grp();

    usi nrow_idx = lconfig_->get_index()->nrowindex();
    usi ncol_idx = rconfig_->get_index()->nrowindex();
    pindex_ = new MatrixIndex(nrow_idx,ncol_idx);

    PermutationGroupPtr lcxn_grp_full
        = pindex_->expand_on_rows(lcxn_grp_);
    PermutationGroupPtr rcxn_grp_full
        = pindex_->expand_on_cols(rcxn_grp_);

    default_grp_ = lcxn_grp_full->union_grp(rcxn_grp_full);
    final_grp_ = default_grp_;

    //create the matrix index for the product configuration
    pconfig_ = new MatrixConfiguration(pindex_, final_grp_);

    bool cxn_on_cols = true;
    lcxn_grp_full = new PermutationGroup(nrow_idx);
    lconfig_->configure_isotropy(lcxn_grp_full, cxn_on_cols);

    rcxn_grp_full = new PermutationGroup(ncol_idx);
    rconfig_->configure_isotropy(rcxn_grp_full, cxn_on_cols);

}

Tensor*
Contraction::get_product_tensor() const
{
    return ptensor_;
}

Tensor*
Contraction::get_right_tensor() const
{
    return rtensor_;
}

Tensor*
Contraction::get_left_tensor() const
{
    return ltensor_;
}

ContractionEngine*
Contraction::get_engine() const
{
    return engine_;
}

void
Contraction::finalize()
{
    for (uli t=0; t < YetiRuntime::nthread(); ++t)
    {
        TensorBlock* tmp_block = cxn_configs_[t]->get_tmp_block();
        if (tmp_block)
        {
            tmp_block->lock();
            tmp_block->finalize();
            tmp_block->unlock();
        }
    }

    YetiRuntime::start_timer("cxn finalize");
    if (use_replicated_product_) //no syncing needs to be done
    {
        Tensor* main_tensor = cxn_configs_[0]->get_product_tensor();
        for (uli t=1; t < YetiRuntime::nthread(); ++t)
        {
            main_tensor->accumulate_no_barrier(cxn_configs_[t]->get_product_tensor(), 1.0);
        }
        main_tensor->global_sum();
        ptensor_->accumulate(main_tensor, 1.0);
    }
    else
    {
        ltensor_->sync();
        rtensor_->sync();
        ptensor_->sync();
    }

    ltensor_->reset_degeneracy();
    rtensor_->reset_degeneracy();

    ptensor_->set_priority(Tensor::gamma_tensor);
    ltensor_->set_priority(Tensor::gamma_tensor);
    rtensor_->set_priority(Tensor::gamma_tensor);
    YetiRuntime::stop_timer("cxn finalize");

    YetiRuntime::get_messenger()->wait_barrier();
}

void
Contraction::print(std::ostream& os) const
{
    os << Env::indent << "Contraction" << endl;
    ++Env::indent;

    os << Env::indent << "Contraction Group" << endl;
    os << contraction_grp_ << endl;

    os << Env::indent << "Default Group" << endl;
    os << default_grp_ << endl;

    --Env::indent;
}

MatrixConfiguration*
Contraction::get_product_matrix_configuration() const
{
    return pconfig_.get();
}

MatrixConfiguration*
Contraction::get_left_matrix_configuration() const
{
    return lconfig_.get();
}

MatrixConfiguration*
Contraction::get_right_matrix_configuration() const
{
    return rconfig_.get();
}

ContractionConfiguration*
Contraction::get_configuration(uli threadnum) const
{
    return cxn_configs_[threadnum];
}

PermutationGroup*
Contraction::get_contraction_grp() const
{
    return contraction_grp_.get();
}

PermutationGroup*
Contraction::get_default_grp() const
{
    return default_grp_.get();
}

PermutationGroup*
Contraction::get_final_grp() const
{
    return final_grp_.get();
}

double
Contraction::scale_factor() const
{
    return scale_;
}

void
Contraction::run()
{
    /** blocks will start communicating from here */
    YetiRuntime::get_messenger()->wait_barrier();

    YetiRuntime::start_timer("cxn configure");
    GlobalTileQueue::configure();
    YetiRuntime::stop_timer("cxn configure");
    YetiRuntime::start_timer("cxn run");
    GlobalTileQueue::run();
    YetiRuntime::stop_timer("cxn run");
}

Contraction::cxn_order_t
Contraction::order() const
{
    return order_;
}

uli
ContractionConfiguration::nrows_left(const uli *sizes) const
{
    return transpose_left_ ?
               lindex_->ncols(sizes) :
               lindex_->nrows(sizes);
}

uli
ContractionConfiguration::ncols_left(const uli *sizes) const
{
    return transpose_left_ ?
               lindex_->nrows(sizes) :
               lindex_->ncols(sizes);
}

uli
ContractionConfiguration::nrows_right(const uli *sizes) const
{
    return transpose_right_ ?
               rindex_->ncols(sizes) :
               rindex_->nrows(sizes);
}

uli
ContractionConfiguration::ncols_right(const uli *sizes) const
{
    return transpose_right_ ?
               rindex_->nrows(sizes) :
               rindex_->ncols(sizes);
}

uli
ContractionConfiguration::ncxn_rows_left() const
{
    return transpose_left_ ? lconfig_->ncols() : lconfig_->nrows();
}

uli
ContractionConfiguration::ncxn_cols_left() const
{
    return transpose_left_ ? lconfig_->nrows() : lconfig_->ncols();
}

uli
ContractionConfiguration::ncxn_rows_right() const
{
    return transpose_right_ ? rconfig_->ncols() : rconfig_->nrows();
}

uli
ContractionConfiguration::ncxn_cols_right() const
{
    return transpose_right_ ? rconfig_->nrows() : rconfig_->ncols();
}

void
ContractionConfiguration::configure_left_block(Tensor* tensor)
{
    configure_left_block(tensor->get_block_map()->sizes(), tensor->get_depth());
}

void
ContractionConfiguration::configure_left_block(MetaDataNode* node)
{
    configure_left_block(node->get_node_map()->sizes(), node->get_depth());
}

void
ContractionConfiguration::configure_left_block(const uli* sizes, usi depth)
{
    lconfig_->configure_block(sizes, depth);
}

void
ContractionConfiguration::configure_right_block(Tensor* tensor)
{
    configure_right_block(tensor->get_block_map()->sizes(), tensor->get_depth());
}

void
ContractionConfiguration::configure_right_block(MetaDataNode* node)
{
    configure_right_block(node->get_node_map()->sizes(), node->get_depth());
}

void
ContractionConfiguration::configure_right_block(const uli* sizes, usi depth)
{
    rconfig_->configure_block(sizes, depth);
}

MatrixConfiguration*
ContractionConfiguration::get_left_config() const
{
    return lconfig_.get();
}

MatrixConfiguration*
ContractionConfiguration::get_right_config() const
{
    return rconfig_.get();
}

MatrixConfiguration*
ContractionConfiguration::get_product_config() const
{
    return pconfig_.get();
}

uli
ContractionConfiguration::get_left_index(uli r, uli c) const
{

    uli idx = transpose_left_ ?
                c * lconfig_->ncols() + r :
                r * lconfig_->ncols() + c;

    return idx;
}

uli
ContractionConfiguration::get_right_index(uli r, uli c) const
{
    uli idx = transpose_right_ ?
                c * rconfig_->ncols() + r :
                r * rconfig_->ncols() + c;

    return idx;
}

void
Contraction::get_row_col_product(uli idx, uli& r, uli& c) const
{
    r = idx / ncols_;
    c = idx % ncols_;
}

uli
Contraction::get_product_index(uli r, uli c) const
{
    uli idx = r * ncols_ + c;
    return idx;
}

TensorBlock*
Contraction::get_product_block(uli r, uli c) 
{
    uli idx = get_product_index(r,c);
    return ptensor_->get_block(idx);
}

void
Contraction::get_row_col_left(uli idx, uli& r, uli& c) const
{
    if (transpose_left_)
    {
        c = idx / nrows_;
        r = idx % nrows_;
    }
    else
    {
        r = idx / nlink_;
        c = idx % nlink_;
    }
}

uli
Contraction::get_left_index(uli r, uli c) const
{
    uli nrows = nrows_;
    uli ncols = nlink_;
    uli idx = transpose_left_ ?
                c * nrows + r :
                r * ncols + c;
    return idx;
}

TensorBlock*
Contraction::get_left_block(uli r, uli c)
{
    uli idx = get_left_index(r, c);
    return ltensor_->get_block(idx);
}

void
Contraction::get_row_col_right(uli idx, uli& r, uli& c) const
{
    if (transpose_right_)
    {
        c = idx / nlink_;
        r = idx % nlink_;
    }
    else
    {
        r = idx / ncols_;
        c = idx % ncols_;
    }
}

uli
Contraction::get_right_index(uli r, uli c) const
{
    uli nrows = nlink_;
    uli ncols = ncols_;
    uli idx = transpose_right_ ?
                c * nrows + r :
                r * ncols + c;
    return idx;
}

TensorBlock*
Contraction::get_right_block(uli r, uli c)
{
    uli idx = get_right_index(r,c);
    return rtensor_->get_block(idx);
}

uli
ContractionConfiguration::get_product_index(uli r, uli c) const
{
    uli idx = r * pconfig_->ncols() + c;
    return idx;
}

void
ContractionConfiguration::reset_contraction_depth(usi depth)
{
    lconfig_->reset_contraction_depth(depth);
    rconfig_->reset_contraction_depth(depth);
}

TensorBlock*
ContractionConfiguration::get_left_block(
    Tensor* tensor,
    uli r, uli c
) const
{
    uli idx = get_left_index(r,c);
    TensorBlock* block = tensor->get_block(idx);
    return block;
}

TileNode*
ContractionConfiguration::get_left_node(
    MetaDataNode* node,
    uli r, uli c, uli& idx
) const
{
    idx = get_left_index(r,c);
    return node->get_node(idx);
}

TensorBlock*
ContractionConfiguration::get_right_block(
    Tensor* tensor,
    uli r, uli c
) const
{
    uli idx = get_right_index(r,c);
    TensorBlock* block = tensor->get_block(idx);
    return block;
}

TileNode*
ContractionConfiguration::get_right_node(
    MetaDataNode* node,
    uli r, uli c, uli& idx
) const
{
    idx = get_right_index(r,c);
    return node->get_node(idx);
}

TensorBlock*
ContractionConfiguration::get_product_block(
    Tensor* tensor,
    uli r, uli c
) const
{
    uli idx = r * ncxn_cols_right() + c;
    return tensor->get_make_block(idx);
}

TileNode*
ContractionConfiguration::get_product_node(
    MetaDataNode* node,
    uli r, uli c
) const
{
    //uli idx = get_product_index(r,c);
    uli idx = r * ncxn_cols_right() + c;
    return node->get_node(idx);
}

TensorBlock*
ContractionConfiguration::get_tmp_block() const
{
    return tmp_block_;
}

bool
Contraction::product_thread_clash() const
{
    return product_thread_clash_;
}

bool
Contraction::product_node_clash() const
{
    return product_node_clash_;
}

TensorBlock*
ContractionTask::get_product_block(
    uli threadnum,
    uli row, 
    uli col
)
{
    TensorBlock* pblock = cxn_->get_configuration(threadnum)
                              ->get_product_block(cxn_->ptensor(), row, col);
    if (pblock && cxn_->ptensor()->is_parent_block(pblock))
        return pblock;
    else
        return 0;
}

TensorBlock*
ContractionTask::get_left_block(
    uli threadnum,
    uli row, 
    uli link
)
{
    TensorBlock* lblock = cxn_->get_configuration(threadnum)
                              ->get_left_block(cxn_->ltensor(), row, link);
    if (lblock && lblock->get_degeneracy() != 0)
        return lblock;
    else
        return 0;
}

TensorBlock*
ContractionTask::get_right_block(
    uli threadnum,
    uli link,
    uli col
)
{
    TensorBlock* rblock = cxn_->get_configuration(threadnum)
                              ->get_right_block(cxn_->rtensor(), link, col);
    if (rblock && rblock->get_degeneracy() != 0)
        return rblock;
    else
        return 0;
}

void
ContractionTask::run_left_driven(uli threadnum)
{
    if (!main_block_)
        return;

    uli ncols = cxn_->ncols();
    TensorBlock* lblock = main_block_;
    lblock->complete_read();

    TensorBlock* rblock = get_right_block(threadnum, link_, 0);
    TensorBlock* pblock = get_product_block(threadnum, row_, 0);
    bool do_current = rblock && pblock && cxn_->do_task(lblock, rblock, pblock);
    for (uli col=1; col < ncols; ++col)
    {
        TensorBlock* next_rblock = get_right_block(threadnum, link_, col);
        TensorBlock* next_pblock = get_product_block(threadnum, row_, col);
        bool do_next = next_rblock && next_pblock && cxn_->do_task(lblock, next_rblock, next_pblock);
        if (do_next)
            next_rblock->prefetch_read();


        if (do_current)
        {
            rblock->complete_read();
            TensorBlock*  my_pblock = cxn_->get_configuration(threadnum)->get_product_block(pblock);
            my_pblock->retrieve_accumulate_no_lock();
            run(threadnum, lblock, rblock, my_pblock);
            my_pblock->release_accumulate();
            rblock->release_read();
        }

        rblock = next_rblock;
        pblock = next_pblock;
        do_current = do_next;
    }


    ContractionTask* task = 0;
    if (follow_permutations_)
    {
        PermutationGroup::iterator it = cxn_->ltensor()->get_tensor_grp()->begin();
        PermutationGroup::iterator stop = cxn_->ltensor()->get_tensor_grp()->end();
        uli permuted_indices[NINDEX];
        uli link, row;
        for ( ; it != stop; ++it)
        {
            Permutation* p = *it;
            if (p->is_identity())
                continue;
            
            p->permute(lblock->get_indices(), permuted_indices);
            uli idx = cxn_->ltensor()->get_index(permuted_indices);
            TensorBlock* new_block = 0;
            if (cxn_->ltensor()->contains(permuted_indices))
                new_block = cxn_->ltensor()->get_block(idx);

            if (new_block == 0 || new_block->is_task_owner()) // Already did this
                continue;

            new_block->prefetch_read(lblock);
            cxn_->get_row_col_left(idx, row, link);
            ContractionTask* next = new ContractionTask(row, 0, link, new_block, cxn_);
            next->next = task;
            task = next;
        }
    }

    /** The last task in the group is not yet done */
    if (task)   
        task->prefetch(threadnum);

    if (do_current)
    {
        rblock->complete_read();
        TensorBlock*  my_pblock = cxn_->get_configuration(threadnum)->get_product_block(pblock);
        my_pblock->retrieve_accumulate_no_lock();
        run(threadnum, lblock, rblock, my_pblock);
        my_pblock->release_accumulate();
        rblock->release_read();
    }

    /** Clear out the task owners */
    while(task)
    {
        Task* old = task;
        if (task->next)
            task->next->prefetch(threadnum);
        task->run(threadnum);
        task = static_cast<ContractionTask*>(task->next);
        delete old;
    }

    lblock->release_read();
    lblock->set_task_owner(false);
}

void
ContractionTask::run_l_r_p(uli threadnum)
{
    run_left_driven(threadnum);
}

Task*
Contraction::get_next_l_r_p()
{
    //link, row, col to ensure maximum reuse of left tensor
    //and reasonable reuse of right tensor

    if (link_ == nlink_)
        return 0;

    Task* task = new ContractionTask(row_, 0, link_, this);

    ++row_;
    if (row_ == nrows_)
    {
        row_ = 0;
        ++link_;
    }

    return task;
}

void
ContractionTask::run_l_p_r(uli threadnum)
{
    run_left_driven(threadnum);
}

Task*
Contraction::get_next_l_p_r()
{
    if (row_ == nrows_)
        return 0;

    Task* task = new ContractionTask(row_, 0, link_, this);
    ++link_;
    if (link_ == nlink_)
    {
        link_ = 0;
        ++row_;
    }

    return task;
}

void
ContractionTask::run_right_driven(uli threadnum)
{
    if (!main_block_)
        return;

    ContractionConfiguration* cxn_config = cxn_->get_configuration(threadnum);
    uli nrows = cxn_->nrows();
    TensorBlock* rblock = main_block_;
    rblock->complete_read();

    TensorBlock* pblock = get_product_block(threadnum, 0, col_);
    TensorBlock* lblock = get_left_block(threadnum, 0, link_);
    bool do_current = lblock && pblock && cxn_->do_task(lblock, rblock, pblock);
    for (uli row=1; row < nrows; ++row)
    {
        TensorBlock* next_lblock = get_left_block(threadnum, row, link_);
        TensorBlock* next_pblock = get_product_block(threadnum, row, col_);
        bool do_next = next_lblock && next_pblock && cxn_->do_task(next_lblock, rblock, next_pblock);
        if (do_next)
            next_lblock->prefetch_read();

        if (do_current)
        {
            lblock->complete_read();
            TensorBlock*  my_pblock = cxn_->get_configuration(threadnum)->get_product_block(pblock);
            my_pblock->retrieve_accumulate_no_lock();
            run(threadnum, lblock, rblock, my_pblock);
            my_pblock->release_accumulate();
            lblock->release_read();
        }

        lblock = next_lblock;
        pblock = next_pblock;
        do_current = do_next;
    }


    ContractionTask* task = 0;
    if (follow_permutations_)
    {
        PermutationGroup::iterator it = cxn_->rtensor()->get_tensor_grp()->begin();
        PermutationGroup::iterator stop = cxn_->rtensor()->get_tensor_grp()->end();
        uli permuted_indices[NINDEX];
        uli link, col;
        for ( ; it != stop; ++it)
        {
            Permutation* p = *it;
            if (p->is_identity())
                continue;
            
            p->permute(rblock->get_indices(), permuted_indices);
            uli idx = cxn_->rtensor()->get_index(permuted_indices);
            TensorBlock* new_block = 0;
            if (cxn_->rtensor()->contains(permuted_indices))
                new_block = cxn_->rtensor()->get_block(idx);

            if (new_block == 0 || new_block->is_task_owner()) // Already did this - this gets set when we prefetch previous tasks
                continue;

            new_block->prefetch_read(rblock);

            cxn_->get_row_col_right(idx, link, col);
            ContractionTask* next = new ContractionTask(0, col, link, new_block, cxn_);

            next->next = task;
            task = next;
        }
    }

    if (task)
        task->prefetch(threadnum);

    if (do_current)
    {
        lblock->complete_read();
        TensorBlock*  my_pblock = cxn_->get_configuration(threadnum)->get_product_block(pblock);
        my_pblock->retrieve_accumulate_no_lock();
        run(threadnum, lblock, rblock, my_pblock);
        my_pblock->release_accumulate();
        lblock->release_read();
    }

    while (task)
    {
        if (task->next)
            task->next->prefetch(threadnum);

        task->run(threadnum);
        Task* old = task;
        task = static_cast<ContractionTask*>(task->next);
        delete old;
    }

    rblock->set_task_owner(false);
    rblock->release_read();
}

void
ContractionTask::run_r_l_p(uli threadnum)
{
    run_right_driven(threadnum);
}

Task*
Contraction::get_next_r_l_p()
{
    //link, col, row to ensure maximum reuse of left tensor
    //and reasonable reuse of right tensor

    if (link_ == nlink_)
        return 0;

    Task* task = new ContractionTask(0, col_, link_, this);
    ++col_;
    if (col_ == ncols_)
    {
        ++link_;
        col_ = 0;
    }
    return task;
}


void
ContractionTask::run_r_p_l(uli threadnum)
{
    run_right_driven(threadnum);
}

Task*
Contraction::get_next_r_p_l()
{
    //col, link, row to ensure maximum reuse of left tensor
    //and reasonable reuse of right tensor

    if (col_ == ncols_)
        return 0;

    Task* task = new ContractionTask(0, col_, link_, this);

    ++link_;
    if (link_ == nlink_)
    {
        ++col_;
        link_ = 0;
    }

    return task;
}

void
ContractionTask::run_product_driven(uli threadnum)
{
    if (!main_block_)
        return;

    ContractionConfiguration* cxn_config = cxn_->get_configuration(threadnum);
    uli nlink = cxn_->nlink();
    TensorBlock* pblock = main_block_;
    TensorBlock* my_pblock = cxn_->get_configuration(threadnum)->get_product_block(pblock);
    my_pblock->retrieve_accumulate_no_lock();

    TensorBlock* rblock = get_right_block(threadnum, 0, col_);
    TensorBlock* lblock = get_left_block(threadnum, row_, 0);
    bool do_current = lblock && rblock && cxn_->do_task(lblock, rblock, pblock);
    for (uli link=1; link < nlink; ++link)
    {
        TensorBlock* next_lblock = get_left_block(threadnum, row_, link);
        TensorBlock* next_rblock = get_right_block(threadnum, link, col_);
        bool do_next = next_lblock && next_rblock && cxn_->do_task(next_lblock, next_rblock, pblock);
        if (do_next)
        {
            next_lblock->prefetch_read();
            next_rblock->prefetch_read();
        }

        if (do_current)
        {
            lblock->complete_read();
            rblock->complete_read();
            run(threadnum, lblock, rblock, my_pblock);
            lblock->release_read();
            rblock->release_read();
        }

        lblock = next_lblock;
        rblock = next_rblock;
        do_current = do_next;
    }

    if (do_current)
    {
        lblock->complete_read();
        rblock->complete_read();
        run(threadnum, lblock, rblock, my_pblock);
        lblock->release_read();
        rblock->release_read();
    }

    my_pblock->release_accumulate();
}

void
ContractionTask::run_p_r_l(uli threadnum)
{
    run_product_driven(threadnum);
}

Task*
Contraction::get_next_p_r_l()
{
    if (col_ == ncols_)
        return 0;

    Task* task = new ContractionTask(row_, col_, 0, this);

    ++row_;
    if (row_ == nrows_)
    {
        ++col_;
        row_ = 0;
    }

    return task;
}


void
ContractionTask::run_p_l_r(uli threadnum)
{
    run_product_driven(threadnum);
}

Task*
Contraction::get_next_p_l_r()
{
    if (row_ == nrows_)
        return 0;

    Task* task = new ContractionTask(row_, col_, 0, this);

    ++col_;
    if (col_ == ncols_)
    { 
        ++row_;
        col_ = 0;
    }

    return task;
}

uli
Contraction::ncols() const
{
    return ncols_;
}

uli
Contraction::nrows() const
{
    return nrows_;
}

uli
Contraction::nlink() const
{
    return nlink_;
}

Tensor*
Contraction::ptensor() const
{
    return ptensor_;
}

Tensor*
Contraction::ltensor() const
{
    return ltensor_;
}

Tensor*
Contraction::rtensor() const
{
    return rtensor_;
}

ContractionConfiguration::ContractionConfiguration(
    const MatrixConfigurationPtr& lconfig,
    const MatrixConfigurationPtr& rconfig,
    const MatrixConfigurationPtr& pconfig,
    Contraction* cxn
) :
    lconfig_(lconfig->copy()),
    rconfig_(rconfig->copy()),
    pconfig_(pconfig->copy()),
    lindex_(lconfig->get_index()),
    rindex_(rconfig->get_index()),
    pindex_(pconfig->get_index()),
    transpose_left_(lconfig->get_index()->is_transpose()),
    transpose_right_(rconfig->get_index()->is_transpose()),
    product_tensor_(0),
    tmp_block_(0),
    next_tmp_block_(0),
    cxn_(cxn)
{
}

ContractionConfiguration::~ContractionConfiguration()
{
    if (product_tensor_)
        delete product_tensor_;

    if (tmp_block_)
    {
        tmp_block_->lock(); //destructor expects lock
        delete tmp_block_;
    }

    if (next_tmp_block_)
    {
        next_tmp_block_->lock(); //destructor expects lock
        delete next_tmp_block_;
    }
}

void
ContractionConfiguration::clone_product_tensor_for_threads(Tensor* tensor)
{
    product_tensor_ = tensor->clone();
}

Tensor*
ContractionConfiguration::get_product_tensor()
{
    return product_tensor_;
}

void
ContractionConfiguration::set_tmp_accumulate_block(TensorBlock* pblock)
{
#if YETI_SANITY_CHECK
    if (tmp_block_)
        yeti_throw(SanityCheckError, "cannot reset tmp accumulate block");
#endif
    tmp_block_ = pblock;
    next_tmp_block_ = new TensorBlock(pblock->get_parent_tensor());
}

TensorBlock*
ContractionConfiguration::get_product_block(TensorBlock* block)
{
    if (product_tensor_)
    {
        TensorBlock* new_block = product_tensor_->get_block(block->get_indices());
        if (!new_block)
        {
            cerr << "Product block " << block->get_block_name() << " is null!" << endl;
            abort();
        }
        new_block->lock();
        return new_block;
    }

    if (cxn_->product_thread_clash() || cxn_->product_node_clash()) 
    {
        if (block->get_malloc_number() == tmp_block_->get_malloc_number())
        {
            //reuse the old tmp block
            tmp_block_->lock();
            return tmp_block_;
        }
        tmp_block_->lock();
        tmp_block_->finalize(); //finalize the old one
        tmp_block_->unlock();
        //switch to the new one
        next_tmp_block_->configure_tmp_block(block);
        TensorBlock* tmp = tmp_block_;
        tmp_block_ = next_tmp_block_;
        next_tmp_block_ = tmp;
        tmp_block_->lock();
        return tmp_block_;
    }
    else
    {
        block->wait_on_remote();
        block->lock();
        return block;
    }
}

void
ContractionConfiguration::check_distribution(
    Tensor* ltensor,
    Tensor* rtensor,
    Tensor* ptensor,
    bool& ltensor_all_rows_local,
    bool& ltensor_all_cols_local,
    bool& rtensor_all_rows_local,
    bool& rtensor_all_cols_local,
    bool& ptensor_all_rows_local,
    bool& ptensor_all_cols_local
)
{
    configure_left_block(ltensor);
    configure_right_block(rtensor);
    uli nrows = ncxn_rows_left();
    uli ncols = ncxn_cols_right();
    uli nlink = ncxn_cols_left();
    TensorBlock* rblock = get_right_block(rtensor, 0, 0);
    rtensor_all_rows_local = 1;
    for (uli r=0; r < nlink; ++r)
    {
        TensorBlock* block = get_right_block(rtensor, r, 0);
        if (block->get_node_owner() != rblock->get_node_owner())
        {
            rtensor_all_rows_local = 0;
            break;
        }
    }

    rtensor_all_cols_local = 1;
    for (uli c=0; c < ncols; ++c)
    {
        TensorBlock* block = get_right_block(rtensor, 0, c);
        if (block->get_node_owner() != rblock->get_node_owner())
        {
            rtensor_all_cols_local = 0;
            break;
        }
    }

    TensorBlock* lblock = get_left_block(ltensor, 0, 0);
    ltensor_all_rows_local = 1;
    for (uli r=0; r < nrows; ++r)
    {
        TensorBlock* block = get_left_block(ltensor, r, 0);
        if (block->get_node_owner() != lblock->get_node_owner())
        {
            ltensor_all_rows_local = 0;
            break;
        }
    }

    ltensor_all_cols_local = 1;
    for (uli c=0; c < nlink; ++c)
    {
        TensorBlock* block = get_left_block(ltensor, 0, c);
        if (block->get_node_owner() != lblock->get_node_owner())
        {
            ltensor_all_cols_local = 0;
            break;
        }
    }
    
    TensorBlock* pblock = get_product_block(ptensor, 0, 0);
    ptensor_all_rows_local = 1;
    for (uli r=0; r < nrows; ++r)
    {
        TensorBlock* block = get_product_block(ptensor, r, 0);
        if (block->get_node_owner() != pblock->get_node_owner())
        {
            ptensor_all_rows_local = 0;
            break;
        }
    }

    ptensor_all_cols_local = 1;
    for (uli c=0; c < ncols; ++c)
    {
        TensorBlock* block = get_product_block(ptensor, 0, c);
        if (block->get_node_owner() != pblock->get_node_owner())
        {
            ptensor_all_cols_local = 0;
            break;
        }
    }


}
