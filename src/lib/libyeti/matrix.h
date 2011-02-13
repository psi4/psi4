#ifndef yeti_matrix_h
#define yeti_matrix_h

#include <list>

#include "class.h"

#include "index.hpp"
#include "matrix.hpp"
#include "tensor.hpp"
#include "tile.hpp"
#include "permutation.hpp"
#include "contraction.hpp"
#include "data.hpp"

#include "mallocimpl.h"
#include "yetiobject.h"


namespace yeti {

/**
    @class MatrixMap
    Encapsulates an array of matrices used for metadata layers of the matrix
*/
class MatrixMap :
    public smartptr::Countable,
    public Malloc<MatrixMap>
{

    private:
        uli nrows_;

        uli ncols_;

        CountableArray<Matrix>* blocks_;

        /**
            Insert a matrix at a given row and column
            @param row
            @param col
            @param m
        */
        void insert(uli row, uli col, Matrix* m);

    public:
        /**
            @param nrows
            @param ncols
        */
        MatrixMap(uli nrows, uli ncols);

        ~MatrixMap();

        /**
            Insert a matrix into the map.  The matrix
            knows its own row,col position
            @param m
        */
        void append(Matrix* m);

        /**
            @param row
            @param col
            @return The matrix at a given row,col position
        */
        Matrix* get(uli row, uli col) const;

        Matrix** begin() const;

        Matrix** end() const;

};

/**
  @class Matrix
        Encapsulates a tensor that has been cast into a matrix.  Essentially, a tensor
        is cast into a matrix by picking which tensor indices should be considered row and col
        indices in a matrix multiplication.
  */
class Matrix :
    public Malloc<Matrix>,
    public YetiRuntimeSerializable
{

    public:
        typedef enum {none, multiplicand, product} matrix_type_t;

    private:
        /**
            The configuration of the matrix defining row, col indices of the underlying tile
            as well as the permutational symmetry of the indices.
        */
        MatrixConfigurationPtr config_;

        IndexRangeTuplePtr tuple_;

        /**
            The set of tiles that generate a matrix with a particular set of
            row,col indices.  Different tiles can be mapped onto a given set of
            row,col indices by different permutations. No tile should be
            repeated here.
        */
        Tile** tiles_;

        /**
            The set of permutations that map #tiles_ onto the given row,col indices.
            No permutation should be repeated here.
        */
        Permutation** perms_;

        /**
            The number of generating tiles and permutations
        */
        usi ngenerator_;

        /**
            Allocate the arrays for the tile/perm generators
        */
        void allocate_generators(usi ngen);

        /**
            Whether this matrix is generated from a uniqe tile. This is always
            true for data matrices.
        */
        bool unique_tile_;

        /**
        */
        void verify_unique_generator();

        /**
            The map of subblocks
        */
        MatrixMapPtr Imap_;

        /**
        */
        uli rowindex_;

        /**
        */
        uli colindex_;

        /**
            The parent matrix. NULL for top matrix.
        */
        Matrix* parent_;

        /**
            Workspace array used for determining permutational uniqueness of certain blocks
        */
        uli* tmpindices_;

        /**

        */
        DataBlockPtr data_;

        /**
            Whether this matrix is considered as a product
            or this matrix is considered as a multiplicand
        */
        matrix_type_t mtype_;

        uli ncols_;

        uli nrows_;

        /**
            The permutational symmetry prefactor to associate with this matrix block.
            For example, in the contraction D(p,q) * G(p,q,r,s) where both tensors
            have permutation symmetry p->q, all blocks with p != q will have
            nisotropy = 2.  All blocks with p=q will only have 1 since they
            are already symmetry unique.
        */
        int nisotropy_;

        /**
            The log of the maximum absolute value
        */
        float maxlog_;

        /**
            Whether the matrix has already been constructed. This skips
            work on subsequent retrievals.
        */
        bool constructed_;

        /**
            Build a subproduct matrix based on a given contraction.  If the
            multiplicand matrices are sufficiently small, this may return
            a NULL pointer.
            @param rowindex
            @param colindex
            @param lmatrix
            @param rmatrix
            @param cxn
        */
        Matrix*
        build_subproduct(
            uli composite_index,
            uli rowindex,
            uli colindex,
            Matrix* lmatrix,
            Matrix* rmatrix,
            Contraction* cxn
        );

        /**
            After generators have been allocated, this computes
            the matrix dimensions from the underlying tile generator
        */
        void compute_dimensions();

        /**
        */
        bool data_equals(const MatrixPtr& matrix);

        /**
            Based on the permutational symmetry of the
            matrix within a given contraction, determine
            whether the next set of matrix indices are necessary
            based on symmetry.
            @param indexset
            @param p
            @param grp
            @param depth
            @param matrix
        */
        bool is_unique_block(
            const size_t* indexset,
            Permutation* p,
            PermutationGroup* grp,
            usi depth,
            Matrix* matrix
        );

        /**
            When retrieving a multiplicand matrix, build
            all of the matrix subblocks and insert them into the matrix map
        */
        void make_blocks();

        bool metadata_equals(const MatrixPtr& matrix);

        void metadata_retrieve(usi depth);

        void retrieve_as_data_multiplicand(DataBlock* dblock, uli threadnum);

        void retrieve_as_metadata_multiplicand(uli threadnum);

        void retrieve_as_data_product(DataBlock* dblock, uli threadnum);

        void retrieve_as_metadata_product(uli threadnum);

        void release_as_multiplicand(uli threadnum);

        void release_as_product(uli threadnum);

        void retrieve_as_multiplicand(uli threadnum);

        void retrieve_as_product(uli threadnum);

        void set_tuple(Tile* tile, Permutation* p);

    protected:
        void _retrieve(uli threadnum);

        void _release(uli threadnum);

    public:
        /**
            Subblock constructor
            @param config
            @param parent
            @param rowindex
            @param colindex
            @param nisotropy  See #nisotropy_.
        */
        Matrix(
            const MatrixConfigurationPtr& config,
            const MatrixPtr& parent,
            size_t rowindex,
            size_t colindex,
            usi nisotropy
        );

        /**
            Contraction constructor for accumulating matrix multiplication
            @param lmatrix
            @param rmatrix
            @param cxn
        */
        Matrix(
            const MatrixPtr& lmatrix,
            const MatrixPtr& rmatrix,
            const ContractionPtr &cxn
        );

        /**
            Top level matrix constructor.
            @param tile The tile to base the matrix on
            @param config The configuration specifying permutational symmetry
                           and which tensor indices are row,col
        */
        Matrix(
            const TensorPtr& tile,
            const MatrixConfigurationPtr& config
        );

        ~Matrix();

        /**
            This can be called from within a thread
            @param lmatrix
            @param rmatrix
            @param cxn
            @param threadnum
        */
        void accumulate_product(
            Matrix* lmatrix,
            Matrix* rmatrix,
            Contraction* cxn,
            uli threadnum
        );

        /**
            This can be called from within a thread
            @param lmatrix
            @param rmatrix
            @param cxn
            @param threadnum
        */
        void accumulate_data_product(
            Matrix* lmatrix,
            Matrix* rmatrix,
            Contraction* cxn,
            uli threadnum
        );

        /**
            This can be called from within a thread
            @param lmatrix
            @param rmatrix
            @param cxn
            @param threadnum
        */
        void accumulate_metadata_product(
            Matrix* lmatrix,
            Matrix* rmatrix,
            Contraction* cxn,
            uli threadnum
        );

        /**
            Add a new tile, permutation pair for generating matrix blocks
            @param tile
            @param p
        */
        void add_generator(
            Tile* tile,
            Permutation* p
        );

        /**
            @return The column index
        */
        size_t colindex() const;

        /**
            @return The metadata depth of the given matrix.  Depth of 0 means data matrix.
        */
        usi depth() const;

        /**
        */
        bool equals(const MatrixPtr& matrix);

        /**
            @param row
            @param col
            @return Matrix subblock at given position
        */
        Matrix* get(size_t row, size_t col) const;

        /**
            @return The underyling data block
        */
        DataBlock* get_data() const;

        IndexRangeTuple* get_index_ranges() const;

        /**
            @return The matrix configuration specifying permutational symmetry and which
                    tensor indices correspond to row,col
        */
        MatrixConfiguration* get_matrix_config() const;

        /**
        */
        MatrixMap* get_map() const;

        /**
            @return The permutation corresponding the first generating tile
        */
        Permutation* get_unique_permutation() const;

        /**
            @return The first generating tile
        */
        Tile* get_unique_tile() const;        

        /**
            The log of the max absolute value of all underlying blocks
        */
        float max_log() const;

        /**
            Make all blocks for a given tile, permutation pair
            @param tile
            @param p
        */
        void make_blocks(
            Tile* tile,
            Permutation* p
        );

        /**
            @return The number of data blocks in the cols of the matrix
        */
        uli ncol_data_blocks();

        /**
            @return The number of data blocks across rows of this matrix
        */
        uli nrow_data_blocks();

        /**
            In certain cases, the number of rows will be "zero" if no row indices.
            This returns 1 instead of zero for contraction purposes.
            @return The number of rows to consider for a contraction.
        */
        uli ncxn_rows() const;

        /**
            In certain cases, the number of cols will be "zero" if no col indices.
            This returns 1 instead of zero for contraction purposes.
            @return The number of cols to consider for a contraction.
        */
        uli ncxn_cols() const;

        /**
            @return The number of cols
        */
        uli ncols() const;

        /**
            @return The number of symmetry unique tiles that generate the given block
        */
        usi ngenerator() const;

        /**
            @return See #nisotropy_
        */
        usi nisotropy() const;

        /**
            @return The number of rows
        */
        uli nrows() const;

        void print(std::ostream& os = std::cout) const;

        void retrieve_to_depth(usi depth);

        /**
            @return The column index for the given matrix
        */
        size_t rowindex() const;

        /**
            Set a unique tile/permutation generator for this matrix
            @param tile
            @param p
        */
        void set_generator(
            Tile* tile,
            Permutation* p
        );

        /**
            Set the matrix to be considered as a product matrix in a contraction
        */
        void set_as_product();

        /**
            Set the matrix to be considered as a multiplicand matrix in a contraction
        */
        void set_as_multiplicand();

};

/**
    @class MatrixIndex
    Class which defines which tensor indices are to be considered as row,col
    indices in a tensor contraction
*/
class MatrixIndex : public smartptr::Serializable {

    public:
        /**
            This determines how the underlying tensor should be considered.
            The row indices are always the contraction indices.  Therefore,
            if the tensor contraction is T(ijab) * T(abkl), the tensors must
            be interpreted as T(ab,ij) * (Tab,kl).  Therefore, the left tensor
            would be specified by rows = {2,3} and cols = {0,1} while
            the right tensor would be specified by rows = {0,1} and cols = {2,3}.
            In this case, the left tensor rows are of type "back" while the left
            tensor cols are of type "front"
        */
        typedef enum { front, back } matrix_index_t;

    private:
        /** The tensor indices to be considered as the composite row index */
        usi* rowindices_;

        /** The tensor indices to be considered as the composite col index */
        usi* colindices_;

        /** The number of tensor indices in the composite row index */
        usi nrowindex_;

        /** The number of tensor indices in the composite col index */
        usi ncolindex_;

        PermutationPtr row_to_col_perm_;

        matrix_index_t rowtype_;

        matrix_index_t coltype_;

        void init();

    public:
        /**
            @param row_type
            @param nrows
            @param col_type
            @param ncols
        */
        MatrixIndex(
            matrix_index_t row_type,
            usi nrows,
            matrix_index_t col_type,
            usi ncols
        );

        /**
            Simplified constructor defining rows as "front" and cols as "back" indices
            @param nrows
            @param ncols
        */
        MatrixIndex(
            usi nrows,
            usi ncols
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

        /**
        */
        matrix_index_t row_index_type() const;

        /**
        */
        matrix_index_t col_index_type() const;

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

        /**
            @return The permutation defining the transpose nature of row,col indices.
                    For tensor T(i,j,a,b) with ij cols and ab rows, this would be
                    P(2,3,0,1) denoting a transpose.
        */
        Permutation* get_row_column_permutation() const;

        void print(std::ostream& os = std::cout) const;

};

/**
    @class MatrixConfiguration
    Class defining which tensor indices compose rows,cols of matrix
    and also permutational symmetry of matrix.  Additionally, the matrix
    might be generated as a sort of the underlying data.
*/
class MatrixConfiguration : public smartptr::Serializable {

    private:
        MatrixIndexPtr index_;

        PermutationPtr indexperm_;

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

        /** The permutation group containing those permutations
            that are neither row nor col permutations */
        PermutationSetPtr quotientset_;

        PermutationGroupPtr cxngrp_;

        PermutationGroupPtr matrix_grp_;

        size_t* workspace_;

    public:
        /**
            @param index
            @param pgrp The permutational symmetry of the matrix
            @param indexperm The permutation defining how the tensor values
                              are to be sorted to form the matrix. This is usually
                              just the identity permutation
        */
        MatrixConfiguration(
           const MatrixIndexPtr& index,
           const PermutationGroupPtr& pgrp,
           const PermutationPtr& indexperm
        );

        virtual ~MatrixConfiguration();

        /**
            @return The permtuation group for row indices
        */
        PermutationGroupPtr get_row_permutation_grp() const;

        /**
            @return The permtuation group for col indices
        */
        PermutationGroupPtr get_col_permutation_grp() const;

        /**
            @return The permtuation group for row indices considered
                    as a subset of all indices. For T(i,j,a,b) with
                    ij rows, ab cols the permuations here are of size 4
                    but act on only 2 elements
        */
        PermutationGroupPtr get_full_row_permutation_grp() const;

        /**
            @return The permtuation group for col indices considered
                    as a subset of all indices. For T(i,j,a,b) with
                    ij rows, ab cols the permuations here are of size 4
                    but act on only 2 elements
        */
        PermutationGroupPtr get_full_col_permutation_grp() const;

        /**
            @return The union of the full col/row permutation groups
        */
        PermutationGroupPtr get_full_permutation_grp() const;

        /**
            @return The full permutation group, minus the elements
                    which are not valid based on a given contraction.
                    Only permutations compatible with a given contraction
                    are kept.  For example, T(i,j,a,b) could have i->j,a->b
                    symmetry but the i->j symmetry is lost in the contraction
                    T(i,j,a,b) * F(b,c)
        */
        PermutationGroupPtr get_matrix_permutation_grp() const;

        /**
            @return The set of permutations needed to generate all unique
                    matrix blocks from the underlying tensor
        */
        PermutationSetPtr get_quotient_set() const;

        /**
            @param tile
            @param p
            @return The number of tile indices fixed by the permutation
        */
        usi nisotropy(
            const TilePtr& tile,
            const PermutationPtr& p
        );

        Matrix* merge_columns();

        /**
            Configure the quotient set used in unpacking the matrix from
            the underlying tensor based on known permutational symmetry
            @param cxn_grp
        */
        void configure_quotient_set(const PermutationGroupPtr& cxn_grp);

        /**
            Configure the matrix permutation group and contraction groups
        */
        void configure_isotropy(const PermutationGroupPtr& grp);

        void print(std::ostream& os = std::cout) const;

        /**
        */
        static void
        product(
            const MatrixIndex* lindex,
            const MatrixIndex* rindex,
            const IndexRangeTuple* ltuple,
            const IndexRangeTuple* rtuple,
            const size_t* lindices,
            const size_t* rindices,
            const Permutation* lperm,
            const Permutation* rperm,
            size_t* indices,
            IndexRangeTuplePtr& tuple
        );

        /**
            Indexer which extracts subset row indices from complete set of indices
            @param tuple
            @param p
            @return An indexer object based upon a given tuple and permutation
        */
        Indexer* rowindex(
            const IndexRangeTuplePtr& tuple,
            const PermutationPtr& p
        ) const;

        /**
            Indexer which extracts subset row indices from complete set of indices
            @param sizes
            @param tuple
            @param p
            @return An indexer object based upon a given tuple and permutation
        */
        Indexer* rowindex(
            const size_t* sizes,
            const size_t* offsets,
            const PermutationPtr& p = 0
        ) const;

        /**
            Indexer which extracts subset col indices from complete set of indices
            @param tuple
            @param p
            @return An indexer object based upon a given tuple and permutation
        */
        Indexer* colindex(
            const IndexRangeTuplePtr& tuple,
            const PermutationPtr& p
        ) const;

        /**
            Indexer which extracts subset col indices from complete set of indices
            @param sizes
            @param tuple
            @param p
            @return An indexer object based upon a given tuple and permutation
        */
        Indexer* colindex(
            const size_t* sizes,
            const size_t* offsets,
            const PermutationPtr& p = 0
        ) const;

        MatrixIndex* get_index() const;


};

class MatrixMergerVector :
    public smartptr::Countable
{

    private:
        CountableArray<DataBlock> blocks_;

        uli nblocks_;

    public:
        MatrixMergerVector(uli nblocks);

};

class MatrixMerger :
    public smartptr::Countable
{

    private:

};

}

#endif
