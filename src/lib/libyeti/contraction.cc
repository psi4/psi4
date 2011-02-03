#include "contraction.h"
#include "matrix.h"
#include "env.h"
#include "permutation.h"
#include "tile.h"
#include "class.h"
#include "tensor.h"
#include "index.h"
#include "exception.h"
#include "contractionimpl.h"
#include "sortimpl.h"
#include "dataimpl.h"
#include "malloc.h"
#include "runtime.h"
#include "threadimpl.h"

using namespace yeti;
using namespace std;


ContractionTask::ContractionTask(
        Matrix* lmatrix,
        Matrix* rmatrix,
        Matrix* product_matrix,
        Contraction* cxn
) :
    lmatrix_(lmatrix),
    rmatrix_(rmatrix),
    product_matrix_(product_matrix),
    cxn_(cxn)
{
    lmatrix->set_as_multiplicand();
    rmatrix->set_as_multiplicand();
    product_matrix->set_as_product();
}

ContractionTask::~ContractionTask()
{
}

void
ContractionTask::print(ostream& os) const
{
    os << Env::indent << "Contraction Task" << endl;
    ++Env::indent;
    os << Env::indent << "Left Matrix    =   " << lmatrix_ << endl;
    os << Env::indent << "Right Matrix   =   " << rmatrix_ << endl;
    os << Env::indent << "Product Matrix =   " << product_matrix_;
    --Env::indent;
}

void
ContractionTask::run(uli threadnum)
{
    lmatrix_->retrieve(threadnum);
    rmatrix_->retrieve(threadnum);
    product_matrix_->retrieve(threadnum);

    product_matrix_->accumulate_product(lmatrix_, rmatrix_, cxn_, threadnum);

    lmatrix_->release(threadnum);
    rmatrix_->release(threadnum);
    product_matrix_->release(threadnum);
}


Contraction::Contraction(
    double scale,
    const TensorPtr& ltensor,
    const TensorPtr& rtensor,
    const MatrixIndexPtr& lindex,
    const MatrixIndexPtr& rindex,
    const PermutationGroupPtr& required_grp,
    const PermutationPtr& lperm,
    const PermutationPtr& rperm,
    const TensorPtr& ptensor
) :
    contraction_grp_(0),
    required_grp_(required_grp),
    product_set_(0),
    special_set_(0),
    default_grp_(0),
    final_grp_(0),
    lconfig_(new MatrixConfiguration(lindex, ltensor->get_permutation_grp(), lperm)),
    rconfig_(new MatrixConfiguration(rindex, rtensor->get_permutation_grp(), rperm)),
    ltensor_(ltensor),
    rtensor_(rtensor),
    product_tensor_(ptensor),
    alphatype_(Contraction::left_tensor),
    alpha_tensor_(0),
    ltmp_(0),
    rtmp_(0),
    ptmp_(0),
    engine_(0),
    lmatrix_(0),
    rmatrix_(0),
    pmatrix_(0),
    scale_(scale)
{
    contraction_grp_ = lconfig_->get_row_permutation_grp()
                        ->intersection_grp(rconfig_->get_row_permutation_grp());

    PermutationGroupPtr lcxngrp(lindex->expand_on_rows(contraction_grp_));
    lconfig_->configure_quotient_set(lcxngrp);
    lconfig_->configure_isotropy(lcxngrp);

    PermutationGroupPtr rcxngrp(rindex->expand_on_rows(contraction_grp_));
    rconfig_->configure_quotient_set(rcxngrp);
    rconfig_->configure_isotropy(rcxngrp);

    usi nprod_idx = lconfig_->get_index()->ncolindex()
                    + rconfig_->get_index()->ncolindex();

    default_grp_ = new PermutationGroup(nprod_idx);

    default_grp_->add_and_expand(
        lconfig_->get_col_permutation_grp(),
        Permutation::ExpandForward
    );
    default_grp_->add_and_expand(
        rconfig_->get_col_permutation_grp(),
        Permutation::ExpandBackward
    );


    /** Configure the special set of permutations that the product tensor
        happens to have.  These are symmetries that arise "accidentally,"
        and cannot be directly determined from the symmetries of the l,r tensors */
    if (product_tensor_) //we have pre-programmed special symmetries
        special_set_ = product_tensor_->get_permutation_grp();
    else //no special symmetries given
        special_set_ = new PermutationGroup(default_grp_->nindex());

    //now add the special set to the default group since these symmetries
    //are automatically obtained
    default_grp_->add(special_set_);
    default_grp_->close();

    final_grp_ = default_grp_->union_grp(required_grp_);

    product_set_ = final_grp_->quotient_set(default_grp_);

    //create the matrix index for the product configuration
    make(product_index, MatrixIndex,
           MatrixIndex::front, lindex->ncolindex(),
           MatrixIndex::back, rindex->ncolindex());
    product_config_ = new MatrixConfiguration(product_index, final_grp_, final_grp_->get_identity());
    product_config_->configure_quotient_set(default_grp_);

    PermutationGroupPtr emptygrp(new PermutationGroup(default_grp_->nindex()));
    product_config_->configure_isotropy(emptygrp);

    if (!product_tensor_) //create a new one
    {
        IndexRangeTuplePtr tuple;
        size_t* indexset = yeti_malloc_indexset();  //not actually relevant
        MatrixConfiguration::product(
                lindex.get(),
                rindex.get(),
                ltensor->get_index_ranges(),
                rtensor->get_index_ranges(),
                ltensor->indices(),
                rtensor->indices(),
                ltensor->get_permutation_grp()->get_identity().get(),
                rtensor->get_permutation_grp()->get_identity().get(),
                indexset,
                tuple
                );
        product_tensor_ = new Tensor("product", tuple, final_grp_);
        //distribute, but do not allocate
        product_tensor_->distribute();
        yeti_free_indexset(indexset);
    }
    else if (!product_tensor_->is_allocated())
    {
        //tensor already exists... just make sure the permutation group is correct
        product_tensor_->get_permutation_grp()->add(final_grp_);
        product_tensor_->get_permutation_grp()->close();

        product_tensor_->retrieve(NOT_THREADED);

        //distribute, but do not allocate
        //allocation comes later
        product_tensor_->distribute();
        product_tensor_->release(NOT_THREADED);
    }
    else
    {
        //for now these permutation groups must be equal
        if (!product_tensor_->get_permutation_grp()->contains(final_grp_) ||
            !final_grp_->contains(product_tensor_->get_permutation_grp()) )
            raise(SanityCheckError, "Contraction cannot accumulate to tensor. Tensor permutation group is not correct.");
    }

    //configure the priorities for the contraction
    vector<TensorPtr> tensors;
    tensors.push_back(ltensor);
    tensors.push_back(rtensor);
    tensors.push_back(product_tensor_);
    sort(tensors.begin(), tensors.end(), tensor_less);

    TensorPtr atensor(tensors[2]);
    atensor->config()->priority = Tensor::alpha_tensor;
    if (atensor == ltensor)
    {
        this->alphatype_ = left_tensor;
        alpha_tensor_ = ltensor_;
    }
    else if (atensor == rtensor)
    {
        this->alphatype_ = right_tensor;
        alpha_tensor_ = rtensor_;
    }
    else if (atensor == product_tensor_)
    {
        this->alphatype_ = product_tensor;
        alpha_tensor_ = product_tensor_;
    }
    else
    {
        raise(SanityCheckError, "no tensor is the alpha tensor");
    }

    tensors[1]->config()->priority = Tensor::beta_tensor;
    tensors[0]->config()->priority = Tensor::gamma_tensor;

    //depending on the matrix types involved
    MatrixIndex::matrix_index_t ltype = lindex->col_index_type();
    MatrixIndex::matrix_index_t rtype = rindex->col_index_type();

    //if rows are the contraction indices, matrix is "transposed"
    //for the contraction
    bool transpose_left = ltype == MatrixIndex::back;
    bool transpose_right = rtype == MatrixIndex::back;

    if (transpose_left && transpose_right)
        engine_ = new ContractionTemplate<Contraction_tt>;
    else if (transpose_left)
        engine_ = new ContractionTemplate<Contraction_tn>;
    else if (transpose_right)
        engine_ = new ContractionTemplate<Contraction_nt>;
    else
        engine_ = new ContractionTemplate<Contraction_nn>;

    //we have to make sure all of the matrices are at the same depth
    //determine the maximum recursion depth
    usi ldepth = ltensor->depth();
    usi rdepth = rtensor->depth();
    usi pdepth = product_tensor_->depth();

    //determine the maximum depth for the 3 matrices
    usi maxdepth = ldepth > rdepth ? ldepth : rdepth;
    maxdepth = pdepth > maxdepth ? pdepth : maxdepth;

   if (   lindex->ncolindex() == 0
       && rindex->ncolindex() == 0
       && product_tensor_->depth() == 0
   ) //dot product
   {
        product_tensor_ = new Tensor(product_tensor_);
   }

    //configure the data mode for each tensor
    ltensor_->set_read_mode();
    rtensor_->set_read_mode();
    product_tensor_->set_accumulate_mode();

    //and now create the matrices to be used
    lmatrix_ = new Matrix(ltensor_, lconfig_);

    rmatrix_ = new Matrix(rtensor_, rconfig_);

    incref();
    pmatrix_ = new Matrix(lmatrix_, rmatrix_, this);
    decref(); //ensure non-deletion

    //now that the product tensor has been built up,
    //allocate the actual memory for it
    product_tensor_->allocate();

    ltmp_ = new char[ltensor_->max_blocksize()];
    rtmp_ = new char[rtensor_->max_blocksize()];
    ptmp_ = new char[product_tensor_->max_blocksize()];

    //register this task with the task queue
    incref(); //keep from getting deleted
    GlobalTileQueue::add(this);
    decref();

    ltensor_->set_read_mode();
    rtensor_->set_read_mode();
    product_tensor_->set_accumulate_mode();
}

Contraction::~Contraction()
{
    delete[] rtmp_;
    delete[] ptmp_;
    delete[] ltmp_;
    delete engine_;
}

Tensor*
Contraction::get_product_tensor() const
{
    return product_tensor_.get();
}

Tensor*
Contraction::get_right_tensor() const
{
    return rtensor_.get();
}

Tensor*
Contraction::get_left_tensor() const
{
    return ltensor_.get();
}

Matrix*
Contraction::get_right_matrix() const
{
    return rmatrix_.get();
}

Matrix*
Contraction::get_left_matrix() const
{
    return lmatrix_.get();
}

Matrix*
Contraction::get_product_matrix() const
{
    return pmatrix_.get();
}

Tensor*
Contraction::get_alpha_tensor() const
{
    return alpha_tensor_.get();
}

void
Contraction::add_task(
    Matrix* lmatrix,
    Matrix* rmatrix,
    Matrix* pmatrix
)
{
    Task* task = new ContractionTask(lmatrix, rmatrix, pmatrix, this);
    if  (alphatype_ == Contraction::left_tensor)
        GlobalTileQueue::add(lmatrix->get_unique_tile(), task);
    else if  (alphatype_ == Contraction::right_tensor)
        GlobalTileQueue::add(rmatrix->get_unique_tile(), task);
    else
        GlobalTileQueue::add(pmatrix->get_unique_tile(), task);
}

ContractionEngine*
Contraction::get_engine() const
{
    return engine_;
}

void
Contraction::finalize()
{
    //all tensors should be in read mode
    ltensor_->set_read_mode();
    rtensor_->set_read_mode();
    product_tensor_->set_read_mode();
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

    if (product_set_->order() > 1)
    {
        os << Env::indent << "Product Set" << endl;
        os << product_set_ << endl;
    }

    --Env::indent;
}

void
Contraction::print_cxn(std::ostream& os)
{

}

MatrixConfiguration*
Contraction::get_product_matrix_configuration() const
{
    return product_config_.get();
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

PermutationSet*
Contraction::get_product_set() const
{
    return product_set_.get();
}

void*
Contraction::left_tmp_space() const
{
    return ltmp_;
}

void*
Contraction::right_tmp_space() const
{
    return rtmp_;
}

void*
Contraction::product_tmp_space() const
{
    return ptmp_;
}

double
Contraction::scale_factor() const
{
    return scale_;
}

void
Contraction::configure()
{
    GlobalTileQueue::configure();
}

void
Contraction::run()
{
    GlobalTileQueue::run();
}

void
Contraction::clear()
{
    GlobalTileQueue::clear();
}

uli
Contraction::nowners()
{
    return GlobalTileQueue::nowners();
}

uli
Contraction::ntasks()
{
    return GlobalTileQueue::ntasks();
}

