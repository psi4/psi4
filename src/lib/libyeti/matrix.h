#ifndef yeti_matrix_h
#define yeti_matrix_h

#include <list>

#include "class.h"

#include "index.hpp"
#include "matrix.hpp"
#include "tensor.hpp"
#include "permutation.hpp"
#include "data.hpp"
#include "sort.hpp"

#include "mallocimpl.h"
#include "yetiobject.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

/**
    @class MatrixIndex
    Class which defines which tensor indices are to be considered as row,col
    indices in a tensor contraction
*/
class MatrixIndex :
    public smartptr::Serializable
{


    private:
        /** The tensor indices to be considered as the composite row index */
        usi rowindices_[NINDEX];

        /** The tensor indices to be considered as the composite col index */
        usi colindices_[NINDEX];

        /** The number of tensor indices in the composite row index */
        usi nrowindex_;

        /** The number of tensor indices in the composite col index */
        usi ncolindex_;

        bool transpose_;

        void init();

    public:
        /**
            @param row_type
            @param nrows
            @param col_type
            @param ncols
        */
        MatrixIndex(
            usi nrows,
            usi ncols,
            bool transpose = false
        );

        virtual ~MatrixIndex();

        /**
            @return The number of tensor indices to be used in computing composite row index
        */
        usi nrowindex() const;

        /**
            @return The number of tensor indices to be used in computing composite col index
        */
        usi ncolindex() const;

        /**
            @return The total number of row,col indices
        */
        usi nindex() const;

        /**
            @return The ith tensor index to be considered in the composite col index.
                    For tensor T(i,j,a,b) with ij rows, ab cols colindex(0) = a, colindex(1) = b
        */
        usi colindex(usi i) const;

        /**
            @return The ith tensor index to be considered in the composite row index.
                    For tensor T(i,j,a,b) with ij rows, ab cols rowindex(0) = i, rowindex(1) = j
        */
        usi rowindex(usi i) const;

        bool is_transpose() const;

        /**
            Takes a permutation group for the row indices and expands into a
            permutation group for all the tensor indices
            @param pgrp The group to expand
            @return The "expanded" permutation group
        */
        PermutationGroupPtr expand_on_rows(const PermutationGroupPtr& pgrp) const;

        /**
            Takes a permutation group for the column indices and expands into a
            permutation group for all the tensor indices
            @param pgrp The group to expand
            @return The "expanded" permutation group
        */
        PermutationGroupPtr expand_on_cols(const PermutationGroupPtr& pgrp) const;

        /**
            @return The array defining col indices. See #colindex.
        */
        const usi* cols() const;

        /**
            @return The array defining row indices. See #rowindex.
        */
        const usi* rows() const;

        uli nrows(const uli* sizes) const;

        uli ncols(const uli* sizes) const;

        void print(std::ostream& os = std::cout) const;

        uli rowstart(const uli* indexstarts, const uli* totalsizes) const;

        uli colstart(const uli* indexstarts, const uli* totalsizes) const;
};

/**
    @class MatrixConfiguration
    Class defining which tensor indices compose rows,cols of matrix
    and also permutational symmetry of matrix.  Additionally, the matrix
    might be generated as a sort of the underlying data.
*/
class MatrixConfiguration :
        public smartptr::Serializable
{

    private:
        uli nrows_;

        uli ncols_;

        uli nrows_at_depth_[MAX_DEPTH];

        uli ncols_at_depth_[MAX_DEPTH];

        MatrixIndexPtr index_;

        /** The full permutation group for the quantity in the matrix */
        PermutationGroupPtr pgrp_;

        /** The permutation group for the underlying tile */
        PermutationGroupPtr tilegrp_;

        /** The permutation group for the row indices,
            considering only the row indices */
        PermutationGroupPtr rowgrp_;

        /** The permutation group for the col indices */
        PermutationGroupPtr colgrp_;

        /** The permutation group for the row indices,
            considering all tensor indices */
        PermutationGroupPtr fullrowgrp_;

        /** The permutation group for the col indices,
            considering all tensor indices */
        PermutationGroupPtr fullcolgrp_;
        
        PermutationSetPtr symmetrization_set_;

        PermutationGroupPtr cxn_grp_;

        PermutationGroupPtr matrix_grp_;

        MatrixConfiguration();

    public:
        /**
            @param index
            @param pgrp The permutational symmetry of the matrix
        */
        MatrixConfiguration(
           const MatrixIndexPtr& index,
           const PermutationGroupPtr& pgrp
        );

        virtual ~MatrixConfiguration();

        PermutationGroup* get_cxn_permutation_grp() const;

        /**
            @return The permtuation group for row indices
        */
        PermutationGroup* get_row_permutation_grp() const;

        /**
            @return The permtuation group for col indices
        */
        PermutationGroup* get_col_permutation_grp() const;

        /**
            @return The permtuation group for row indices considered
                    as a subset of all indices. For T(i,j,a,b) with
                    ij rows, ab cols the permuations here are of size 4
                    but act on only 2 elements
        */
        PermutationGroup* get_full_row_permutation_grp() const;

        /**
            @return The permtuation group for col indices considered
                    as a subset of all indices. For T(i,j,a,b) with
                    ij rows, ab cols the permuations here are of size 4
                    but act on only 2 elements
        */
        PermutationGroup* get_full_col_permutation_grp() const;

        /**
            @return The union of the full col/row permutation groups
        */
        PermutationGroup* get_full_permutation_grp() const;

        /**
            @return The full permutation group, minus the elements
                    which are not valid based on a given contraction.
                    Only permutations compatible with a given contraction
                    are kept.  For example, T(i,j,a,b) could have i->j,a->b
                    symmetry but the i->j symmetry is lost in the contraction
                    T(i,j,a,b) * F(b,c)
        */
        PermutationGroup* get_matrix_permutation_grp() const;

        PermutationSet* get_symmetrization_set() const;

        /**
            Configure the quotient set used in unpacking the matrix from
            the underlying tensor based on known permutational symmetry
            @param cxn_grp
        */
        void configure_symmetrization_set(const PermutationGroupPtr& cxn_grp);

        /**
            Configure the matrix permutation group and contraction groups
        */
        void configure_isotropy(const PermutationGroupPtr& grp, bool cxn_on_cols);
        
        void configure_matrix_grp(const PermutationGroupPtr& grp);

        void print(std::ostream& os = std::cout) const;

        MatrixIndex* get_index() const;

        void configure_block(
            const uli* sizes,
            usi depth
        );

        void reset_contraction_depth(usi depth);

        uli nrows() const;

        uli ncols() const;

        MatrixConfiguration* copy() const;
};



} //end namespace

#ifdef redefine_size_t
#undef size_t
#endif

#endif
