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
#include <algorithm>

using namespace yeti;
using namespace std;

#define DEBUG_CXN 0

#define PRINT_CXN_TYPE 0

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

ContractionTask::ContractionTask(
    TensorBlock* lblock,
    TensorBlock* rblock,
    TensorBlock* pblock,
    Contraction* cxn
) :
    lblock_(lblock),
    rblock_(rblock),
    pblock_(pblock),
    cxn_(cxn)
{
    if (lblock->get_degeneracy() != rblock->get_degeneracy())
    {
        cerr << "Degeneracy of contraction blocks do not match!" << endl
            << "Permutational symmetry has not been properly configured." << endl;
        cerr << lblock->get_parent_tensor()->get_name() 
             << lblock->get_index_string() << " " << lblock->get_degeneracy() << endl;
        cerr << rblock->get_parent_tensor()->get_name() 
             << rblock->get_index_string() << " " << rblock->get_degeneracy() << endl;
        abort();
    }
}

ContractionTask::~ContractionTask()
{
}

void
ContractionTask::print(ostream& os) const
{
    os << Env::indent << "Contraction Task" << endl;
    ++Env::indent;
    --Env::indent;
}

void
ContractionTask::run(uli threadnum)
{
#if DEBUG_CXN
                dout << "Accumulating contraction task "
            << ClassOutput<const uli*>::str(pblock_->get_descr()->nindex(), pblock_->get_indices()) << " = "
            << ClassOutput<const uli*>::str(lblock_->get_descr()->nindex(), lblock_->get_indices()) << " * "
            << ClassOutput<const uli*>::str(rblock_->get_descr()->nindex(), rblock_->get_indices())
            << endl;
#endif

    uli threadtest = YetiRuntime::get_thread_number();
    if (threadtest != threadnum)
    {
        cerr << "thread number wrong" << endl;
    }



    lblock_->set_read_mode();
    rblock_->set_read_mode();

    lblock_->retrieve();
    rblock_->retrieve();


    bool has_lock = pblock_->trylock();
    if (!has_lock)
    {
        TensorBlock* tmpblock = new TensorBlock(pblock_);
        tmpblock->set_accumulate_mode();
        tmpblock->retrieve();
        tmpblock->accumulate(lblock_, rblock_, cxn_);
        pblock_->lock();
        pblock_->set_accumulate_mode();
        pblock_->retrieve_nolock();
        pblock_->accumulate(tmpblock, 1.0, 0);
        pblock_->release_lock();
        tmpblock->release();
        delete tmpblock;

    }
    else
    {
        pblock_->set_accumulate_mode();
        pblock_->retrieve_nolock();
        pblock_->accumulate(lblock_, rblock_, cxn_);
        pblock_->release_nolock();
        pblock_->unlock();
    }

    lblock_->release();
    rblock_->release();



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
    ltensor_(ltensor),
    rtensor_(rtensor),
    ptensor_(ptensor),
    alphatype_(Contraction::left_tensor),
    engine_(0),
    scale_(scale),
    cxn_configs_(new ContractionConfiguration*[YetiRuntime::nthread()])
{

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
    std::sort(tensors.begin(), tensors.end(), tensor_less);

    Tensor*  atensor(tensors[2]);
    if (atensor == ltensor_)
        this->alphatype_ = left_tensor;
    else if (atensor == rtensor_)
        this->alphatype_ = right_tensor;
    else if (atensor == ptensor_)
        this->alphatype_ = product_tensor;
    else
        raise(SanityCheckError, "no tensor is the alpha tensor");

    Tensor* btensor(tensors[1]);
    if (btensor == ltensor_)
        this->betatype_ = left_tensor;
    else if (btensor == rtensor_)
        this->betatype_ = right_tensor;
    else if (btensor == ptensor_)
        this->betatype_ = product_tensor;
    else
        raise(SanityCheckError, "no tensor is the beta tensor");

    Tensor* gtensor(tensors[0]);
    if (gtensor == ltensor_)
        this->gammatype_ = left_tensor;
    else if (gtensor == rtensor_)
        this->gammatype_ = right_tensor;
    else if (gtensor == ptensor_)
        this->gammatype_ = product_tensor;
    else
        raise(SanityCheckError, "no tensor is the gamma tensor");



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
        raise(SanityCheckError, "tensors not aligned for multiplication");

    ltensor_->configure_degeneracy(lconfig_->get_cxn_permutation_grp());
    if (ltensor_ != rtensor_)
        rtensor_->configure_degeneracy(rconfig_->get_cxn_permutation_grp());

    uli nthread = YetiRuntime::nthread();
    for (uli i=0; i < nthread; ++i)
    {
        cxn_configs_[i] = new ContractionConfiguration(lconfig_, rconfig_, pconfig_);
    }
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
Contraction::build()
{
    //register this task with the task queue
    incref(); //keep from getting deleted
    GlobalTileQueue::add(this);
    decref();

    ContractionConfiguration* cxn_config = cxn_configs_[0];

    //configure everything
    cxn_config->configure_left_block(ltensor_);
    cxn_config->configure_right_block(rtensor_);
    //configure_product_block(ptensor_);

    if  (alphatype_ == left_tensor)
    {
        if (betatype_ == right_tensor)
            build_l_r_p();
        else
            build_l_p_r();
    }
    else if (alphatype_ == right_tensor)
    {
        if (betatype_ == left_tensor)
            build_r_l_p();
        else
            build_r_p_l();
    }
    else if (alphatype_ == product_tensor)
    {
        if (betatype_ == left_tensor)
            build_p_l_r();
        else
            build_p_r_l();
    }
    else
        raise(SanityCheckError, "contraction is clearly misconfigured");
}

void
Contraction::finalize()
{
    ptensor_->sync();
    ptensor_->update();

    ltensor_->reset_degeneracy();
    rtensor_->reset_degeneracy();

    //cout << ltensor_->nblocks_retrieved() << " blocks retrieved on left tensor" << endl;
    //cout << rtensor_->nblocks_retrieved() << " blocks retrieved on right tensor" << endl;
    //cout << ptensor_->nblocks_retrieved() << " blocks retrieved on product tensor" << endl;
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
    GlobalTileQueue::configure();
    GlobalTileQueue::run();
}

void
Contraction::clear()
{
    GlobalTileQueue::clear();
}

uli
Contraction::ntasks()
{
    return GlobalTileQueue::ntasks();
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

void
Contraction::build_l_r_p()
{
    ContractionConfiguration* cxn_config = cxn_configs_[0];
#if PRINT_CXN_TYPE
    dout << "lrp" << endl;
#endif
    uli nrows = cxn_config->ncxn_rows_left();
    uli ncols = cxn_config->ncxn_cols_right();
    uli nlink = cxn_config->ncxn_cols_left();
#if YETI_SANITY_CHECK
    uli check = cxn_config->ncxn_rows_right();
    if (check != nlink)
        raise(SanityCheckError, "matrices not aligned for multiplication");
#endif
    //link, row, col to ensure maximum reuse of left tensor
    //and reasonable reuse of right tensor
    for (uli link=0; link < nlink; ++link)
    {
        for (uli row=0; row < nrows; ++row)
        {
            TensorBlock* lblock = cxn_config->get_left_block(ltensor_, row, link);
            if (!lblock || lblock->get_degeneracy() == 0)
            {
                continue;
            }

            for (uli col=0; col < ncols; ++col)
            {
                TensorBlock* rblock = cxn_config->get_right_block(rtensor_, link, col);
                if (!rblock || rblock->get_degeneracy() == 0)
                {
                    continue;
                }

#if USE_SCREENING
                if (lblock->get_max_log() + rblock->get_max_log() 
                    < YetiRuntime::matrix_multiply_cutoff)
                    continue;
#endif

                TensorBlock* pblock = cxn_config->get_product_block(ptensor_, row, col);
                if (!pblock) //might be rigorously zero
                    continue;

                if (!ptensor_->is_parent_block(pblock))
                {
                    continue;
                }

                ContractionTask* task =
                    new ContractionTask(lblock, rblock, pblock, this);

                //the cxn priority we get here may not be the same as
                //the cxn priority set above as previous contractions may
                //have already configured contraction priorities
                GlobalTileQueue::add(lblock, task);
            }
        }
    }
}

void
Contraction::build_l_p_r()
{
    ContractionConfiguration* cxn_config = cxn_configs_[0];
#if PRINT_CXN_TYPE
    dout << "lpr" << endl;
#endif
    uli nrows = cxn_config->ncxn_rows_left();
    uli ncols = cxn_config->ncxn_cols_right();
    uli nlink = cxn_config->ncxn_cols_left();
#if YETI_SANITY_CHECK
    uli check = cxn_config->ncxn_rows_right();
    if (check != nlink)
        raise(SanityCheckError, "matrices not aligned for multiplication");
#endif
    //row, link, col to ensure maximum reuse of left tensor
    //and reasonable reuse of right tensor
    for (uli row=0; row < nrows; ++row)
    {
        for (uli link=0; link < nlink; ++link)
        {
            TensorBlock* lblock = cxn_config->get_left_block(ltensor_, row, link);
            if (!lblock || lblock->get_degeneracy() == 0)
                continue;
                
            for (uli col=0; col < ncols; ++col)
            { 
                TensorBlock* rblock = cxn_config->get_right_block(rtensor_, link, col);
                if (!rblock || rblock->get_degeneracy() == 0)
                    continue;

#if USE_SCREENING
                if (lblock->get_max_log() + rblock->get_max_log() 
                    < YetiRuntime::matrix_multiply_cutoff)
                {
                    continue;
                }
#endif

                TensorBlock* pblock = cxn_config->get_product_block(ptensor_, row, col);
                if (!pblock) //might be rigorously zero
                    continue;

                if (!ptensor_->is_parent_block(pblock))
                {
                    continue;
                }

                ContractionTask* task =
                    new ContractionTask(lblock, rblock, pblock, this);

                //the cxn priority we get here may not be the same as
                //the cxn priority set above as previous contractions may
                //have already configured contraction priorities
                GlobalTileQueue::add(lblock, task);
            }
        }
    }
}

void
Contraction::build_r_l_p()
{
    ContractionConfiguration* cxn_config = cxn_configs_[0];
#if PRINT_CXN_TYPE
    dout << "rlp" << endl;
#endif
    uli nrows = cxn_config->ncxn_rows_left();
    uli ncols = cxn_config->ncxn_cols_right();
    uli nlink = cxn_config->ncxn_cols_left();
#if YETI_SANITY_CHECK
    uli check = cxn_config->ncxn_rows_right();
    if (check != nlink)
        raise(SanityCheckError, "matrices not aligned for multiplication");
#endif
    //link, col, row to ensure maximum reuse of left tensor
    //and reasonable reuse of right tensor

    for (uli link=0; link < nlink; ++link)
    {
        for (uli col=0; col < ncols; ++col)
        {
            TensorBlock* rblock = cxn_config->get_right_block(rtensor_, link, col);
            if (!rblock || rblock->get_degeneracy() == 0)
            {
                continue;
            }
                
            for (uli row=0; row < nrows; ++row)
            {
                TensorBlock* lblock = cxn_config->get_left_block(ltensor_, row, link);
                if (!lblock || lblock->get_degeneracy() == 0)
                {
                    continue;
                }

#if USE_SCREENING
                if (lblock->get_max_log() + rblock->get_max_log() 
                    < YetiRuntime::matrix_multiply_cutoff)
                    continue;
#endif

                TensorBlock* pblock = cxn_config->get_product_block(ptensor_, row, col);
                if (!pblock) //might be rigorously zero
                    continue;


                if (!ptensor_->is_parent_block(pblock))
                    continue;

                ContractionTask* task =
                    new ContractionTask(lblock, rblock, pblock, this);

                //the cxn priority we get here may not be the same as
                //the cxn priority set above as previous contractions may
                //have already configured contraction priorities
                GlobalTileQueue::add(rblock, task);




            }
        }
    }
}

void
Contraction::build_r_p_l()
{
    ContractionConfiguration* cxn_config = cxn_configs_[0];
#if PRINT_CXN_TYPE
    dout << "rpl" << endl;
#endif
    uli nrows = cxn_config->ncxn_rows_left();
    uli ncols = cxn_config->ncxn_cols_right();
    uli nlink = cxn_config->ncxn_cols_left();

#if YETI_SANITY_CHECK
    uli check = cxn_config->ncxn_rows_right();
    if (check != nlink)
        raise(SanityCheckError, "matrices not aligned for multiplication");
#endif
    //col, link, row to ensure maximum reuse of left tensor
    //and reasonable reuse of right tensor
    for (uli col=0; col < ncols; ++col)
    {
        for (uli link=0; link < nlink; ++link)
        {
            
            TensorBlock* rblock = cxn_config->get_right_block(rtensor_, link, col);
            if (!rblock || rblock->get_degeneracy() == 0)
                continue;

            for (uli row=0; row < nrows; ++row)
            {
                TensorBlock* lblock = cxn_config->get_left_block(ltensor_, row, link);
                if (!lblock || lblock->get_degeneracy() == 0)
                    continue;
#if USE_SCREENING
                if (lblock->get_max_log() + rblock->get_max_log() 
                    < YetiRuntime::matrix_multiply_cutoff)
                    continue;
#endif

                TensorBlock* pblock = cxn_config->get_product_block(ptensor_, row, col);
                if (!pblock) //might be rigorously zero
                    continue;

                if (!ptensor_->is_parent_block(pblock))
                    continue;

                ContractionTask* task =
                    new ContractionTask(lblock, rblock, pblock, this);

                //the cxn priority we get here may not be the same as
                //the cxn priority set above as previous contractions may
                //have already configured contraction priorities
                GlobalTileQueue::add(rblock, task);
            }
        }
    }
}

void
Contraction::build_p_r_l()
{
    ContractionConfiguration* cxn_config = cxn_configs_[0];
#if PRINT_CXN_TYPE
    dout << "prl" << endl;
#endif
    uli nrows = cxn_config->ncxn_rows_left();
    uli ncols = cxn_config->ncxn_cols_right();
    uli nlink = cxn_config->ncxn_cols_left();
#if YETI_SANITY_CHECK
    uli check = cxn_config->ncxn_rows_right();
    if (check != nlink)
        raise(SanityCheckError, "matrices not aligned for multiplication");
#endif
    //col, row, link to ensure maximum reuse of left tensor
    //and reasonable reuse of right tensor
    for (uli col=0; col < ncols; ++col)
    {
        for (uli row=0; row < nrows; ++row)
        {
            for (uli link=0; link < nlink; ++link)
            {
                TensorBlock* rblock = cxn_config->get_right_block(rtensor_, link, col);
                if (!rblock || rblock->get_degeneracy() == 0)
                    continue;

                TensorBlock* lblock = cxn_config->get_left_block(ltensor_, row, link);
                if (!lblock || lblock->get_degeneracy() == 0)
                    continue;

#if USE_SCREENING
                if (lblock->get_max_log() + rblock->get_max_log() 
                    < YetiRuntime::matrix_multiply_cutoff)
                {
                    continue;
                }
#endif

                TensorBlock* pblock = cxn_config->get_product_block(ptensor_, row, col);
                if (!pblock) //might be rigorously zero
                {
                    continue;
                }

                if (!ptensor_->is_parent_block(pblock))
                {
                    continue;
                }

                ContractionTask* task =
                        new ContractionTask(lblock, rblock, pblock, this);

                //the cxn priority we get here may not be the same as
                //the cxn priority set above as previous contractions may
                //have already configured contraction priorities
                GlobalTileQueue::add(pblock, task);
            }
        }
    }
}

void
Contraction:: build_p_l_r()
{
    ContractionConfiguration* cxn_config = cxn_configs_[0];
#if PRINT_CXN_TYPE
    dout << "plr" << endl;
#endif
    uli nrows = cxn_config->ncxn_rows_left();
    uli ncols = cxn_config->ncxn_cols_right();
    uli nlink = cxn_config->ncxn_cols_left();
#if YETI_SANITY_CHECK
    uli check = cxn_config->ncxn_rows_right();
    if (check != nlink)
        raise(SanityCheckError, "matrices not aligned for multiplication");
#endif
    //row, col, link to ensure maximum reuse of left tensor
    //and reasonable reuse of right tensor
    for (uli row=0; row < nrows; ++row)
    {
        for (uli col=0; col < ncols; ++col)
        {
            for (uli link=0; link < nlink; ++link)
            {
                uli prod_idx = row * ncols + col;

                TensorBlock* lblock = cxn_config->get_left_block(ltensor_, row, link);
                if (!lblock || lblock->get_degeneracy() == 0)
                {
                    continue;
                }

                TensorBlock* rblock = cxn_config->get_right_block(rtensor_, link, col);
                if (!rblock || rblock->get_degeneracy() == 0)
                {
                    continue;
                }

#if USE_SCREENING
                if (lblock->get_max_log() + rblock->get_max_log() 
                    < YetiRuntime::matrix_multiply_cutoff)
                {

                    continue;
                }
#endif

                TensorBlock* pblock = cxn_config->get_product_block(ptensor_, row, col);
                if (!pblock) //might be rigorously zero
                {
                    continue;
                }

                if (!ptensor_->is_parent_block(pblock))
                {
                     continue;
                }

                ContractionTask* task =
                        new ContractionTask(lblock, rblock, pblock, this);

                //the cxn priority we get here may not be the same as
                //the cxn priority set above as previous contractions may
                //have already configured contraction priorities
                GlobalTileQueue::add(pblock, task);

            }
        }
    }
}

ContractionConfiguration::ContractionConfiguration(
    const MatrixConfigurationPtr& lconfig,
    const MatrixConfigurationPtr& rconfig,
    const MatrixConfigurationPtr& pconfig
) :
    lconfig_(lconfig->copy()),
    rconfig_(rconfig->copy()),
    pconfig_(pconfig->copy()),
    lindex_(lconfig->get_index()),
    rindex_(rconfig->get_index()),
    pindex_(pconfig->get_index()),
    transpose_left_(lconfig->get_index()->is_transpose()),
    transpose_right_(rconfig->get_index()->is_transpose())
{
}

