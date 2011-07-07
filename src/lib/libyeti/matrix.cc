#include "matrix.h"
#include "index.h"
#include "data.h"
#include "class.h"
#include "permutation.h"
#include "exception.h"
#include "env.h"
#include "tensor.h"
#include "sort.h"
#include "malloc.h"
#include "runtime.h"

using namespace yeti;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

MatrixConfiguration::MatrixConfiguration(
    const MatrixIndexPtr& index,
    const PermutationGroupPtr& pgrp
) : index_(index),
    pgrp_(pgrp),
    rowgrp_(0),
    colgrp_(0),
    fullrowgrp_(0),
    fullcolgrp_(0),
    nrows_(0),
    ncols_(0),
    cxn_grp_(0),
    matrix_grp_(0),
    tilegrp_(pgrp)
{
    fullrowgrp_ = pgrp_->subgrp(index->rows(), index->nrowindex());
    fullcolgrp_ = pgrp_->subgrp(index->cols(), index->ncolindex());

    rowgrp_ = pgrp_->compressed_subgrp(index->rows(), index->nrowindex());
    colgrp_ = pgrp_->compressed_subgrp(index->cols(), index->ncolindex());
}

MatrixConfiguration::MatrixConfiguration()
    : index_(0),
    pgrp_(0),
    rowgrp_(0),
    colgrp_(0),
    fullrowgrp_(0),
    fullcolgrp_(0),
    nrows_(0),
    ncols_(0),
    cxn_grp_(0),
    matrix_grp_(0),
    tilegrp_(0)
{
}

MatrixConfiguration*
MatrixConfiguration::copy() const
{
    MatrixConfiguration* cpy = new MatrixConfiguration;
    cpy->index_ = index_;
    cpy->pgrp_ = pgrp_;
    cpy->rowgrp_ = rowgrp_;
    cpy->colgrp_ = colgrp_;
    cpy->fullrowgrp_ = fullrowgrp_;
    cpy->fullcolgrp_ = fullcolgrp_;
    cpy->matrix_grp_ = matrix_grp_;
    cpy->tilegrp_ = tilegrp_;
    return cpy;
}

MatrixConfiguration::~MatrixConfiguration()
{
}

void
MatrixConfiguration::configure_block(
    const uli *sizes,
    usi depth
)
{
    nrows_ = index_->nrows(sizes);
    ncols_ = index_->ncols(sizes);
    nrows_at_depth_[depth] = nrows_;
    ncols_at_depth_[depth] = ncols_;
}

void
MatrixConfiguration::reset_contraction_depth(usi depth)
{
    nrows_ = nrows_at_depth_[depth];
    ncols_ = ncols_at_depth_[depth];
};

uli
MatrixConfiguration::nrows() const
{
    return nrows_;
}

uli
MatrixConfiguration::ncols() const
{
    return ncols_;
}

void
MatrixConfiguration::configure_isotropy(
    const PermutationGroupPtr& grp,
    bool cxn_on_cols
)
{
    cxn_grp_ = grp;
    
    if (cxn_on_cols)
        matrix_grp_ = cxn_grp_->union_grp(this->fullrowgrp_);
    else
        matrix_grp_ = cxn_grp_->union_grp(this->fullcolgrp_);
}

void
MatrixConfiguration::configure_matrix_grp(const PermutationGroupPtr& grp)
{
    cxn_grp_ = grp;
    matrix_grp_ = grp;
}

void
MatrixConfiguration::configure_symmetrization_set(
    const PermutationGroupPtr& cxn_grp
)
{
    PermutationGroupPtr target_grp(cxn_grp->union_grp(fullcolgrp_));
    PermutationGroupPtr quotient_grp = target_grp->intersection_grp(tilegrp_);
    if (tilegrp_->order() > quotient_grp->order())
    {
        //symmetrization is required
        PermutationSetPtr quotient_set = tilegrp_->quotient_set(quotient_grp);
        symmetrization_set_ = quotient_set->get_generator_set();
    }
    else
    {
        symmetrization_set_ = new PermutationSet(quotient_grp->nindex());
    }
}

PermutationGroup*
MatrixConfiguration::get_col_permutation_grp() const
{
    return colgrp_.get();
}

PermutationGroup*
MatrixConfiguration::get_cxn_permutation_grp() const
{
    return cxn_grp_.get();
}

PermutationGroup*
MatrixConfiguration::get_full_col_permutation_grp() const
{
    return fullcolgrp_.get();
}

PermutationGroup*
MatrixConfiguration::get_full_permutation_grp() const
{
    return pgrp_.get();
}

PermutationGroup*
MatrixConfiguration::get_full_row_permutation_grp() const
{
    return fullrowgrp_.get();
}

MatrixIndex*
MatrixConfiguration::get_index() const
{
    return index_.get();
}

PermutationGroup*
MatrixConfiguration::get_matrix_permutation_grp() const
{
    return matrix_grp_.get();
}

PermutationSet*
MatrixConfiguration::get_symmetrization_set() const
{
    return symmetrization_set_.get();
}

PermutationGroup*
MatrixConfiguration::get_row_permutation_grp() const
{
    return rowgrp_.get();
}

void
MatrixConfiguration::print(std::ostream &os) const
{
    os << "Matrix Index Configuration" << endl;
    os << this->index_ << endl;
}

MatrixIndex::MatrixIndex(
    usi nrows,
    usi ncols,
    bool transpose
) :
    nrowindex_(nrows),
    ncolindex_(ncols),
    transpose_(transpose)
{
    init();
}


MatrixIndex::~MatrixIndex()
{
}

void
MatrixIndex::init()
{

    usi pmap[NINDEX];
    for (usi i=0; i < nrowindex_; ++i)
    {
        rowindices_[i] = i; //rowindices are offset by the number of params
        pmap[i] = i; //parameters are only a single index and do not affect permutations
    }
    for (usi i=0; i < ncolindex_; ++i)
    {
        colindices_[i] = i + nrowindex_;
        pmap[i + nrowindex_] = i + nrowindex_;
    }

}

usi
MatrixIndex::colindex(usi i) const
{
    return colindices_[i];
}

bool
MatrixIndex::is_transpose() const
{
    return transpose_;
}

const usi*
MatrixIndex::cols() const
{
    return colindices_;
}

uli
MatrixIndex::ncols(const uli* sizes) const
{
    uli nc = 1;
    for (usi i=0; i < ncolindex_; ++i)
        nc *= sizes[colindices_[i]];
    return nc;
}

uli
MatrixIndex::nrows(const uli* sizes) const
{
    uli nr = 1;
    for (usi i=0; i < nrowindex_; ++i)
    {
        nr *= sizes[rowindices_[i]];
    }

    return nr;
}

PermutationGroupPtr
MatrixIndex::expand_on_cols(const PermutationGroupPtr& pgrp) const
{
    PermutationGroup* newgrp = new PermutationGroup(nindex(), pgrp, colindices_);
    return newgrp;
}

PermutationGroupPtr
MatrixIndex::expand_on_rows(const PermutationGroupPtr& pgrp) const
{
    PermutationGroup* newgrp = new PermutationGroup(nindex(), pgrp, rowindices_);
    return newgrp;
}

usi
MatrixIndex::ncolindex() const
{
    return ncolindex_;
}

usi
MatrixIndex::nindex() const
{
    return ncolindex_ + nrowindex_;
}

usi
MatrixIndex::nrowindex() const
{
    return nrowindex_;
}

void
MatrixIndex::print(ostream& os) const
{
    os << Env::indent << "Matrix Index" << endl;

    ++Env::indent;

    os << Env::indent << "Row Indices: ";
    for (usi i=0; i < nrowindex_; ++i)
        os << " " << rowindices_[i];

    os << Env::indent << "Col Indices: ";
    for (usi i=0; i < ncolindex_; ++i)
        os << " " << colindices_[i];

    --Env::indent;
}

usi
MatrixIndex::rowindex(usi i) const
{
    return rowindices_[i];
}

const usi*
MatrixIndex::rows() const
{
    return rowindices_;
}

uli
MatrixIndex::rowstart(
    const uli* indexstarts,
    const uli* totalsizes
) const
{
    uli idx = 0;
    uli stride = 1;
    for (usi i=1; i <= nrowindex_; ++i)
    {
        usi idxnumber = nrowindex_ - i;
        idx += indexstarts[idxnumber] * stride;
        stride *= totalsizes[idxnumber];
    }
    return idx;
}


uli
MatrixIndex::colstart(
    const uli* indexstarts,
    const uli* totalsizes
) const
{
    usi nindex = ncolindex_ + nrowindex_;
    uli idx = 0;
    uli stride = 1;
    for (usi i=1; i <= ncolindex_; ++i)
    {
        usi idxnumber = nindex - i;
        idx += indexstarts[idxnumber] * stride;
        stride *= totalsizes[idxnumber];
    }
    return idx;
}
