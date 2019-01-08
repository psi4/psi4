/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _psi_src_lib_libmints_matrix_h_
#define _psi_src_lib_libmints_matrix_h_

#include "psi4/libmints/dimension.h"
#include "psi4/libmints/typedefs.h"
#include "psi4/libpsi4util/exception.h"

#include <cstdio>
#include <string>
#include <vector>
#include <memory>

namespace psi {

struct dpdfile2;

class PSIO;
class Vector;
class Dimension;
class Molecule;
class Vector3;

enum diagonalize_order { evals_only_ascending = 0, ascending = 1, evals_only_descending = 2, descending = 3 };

/*! \ingroup MINTS
 *  \class Matrix
 *  \brief Makes using matrices just a little easier.
 *
 * Using a matrix factory makes creating these a breeze.
 */
class PSI_API Matrix : public std::enable_shared_from_this<Matrix> {
   protected:
    /// Matrix data
    double*** matrix_;
    /// Number of irreps
    int nirrep_;
    /// Rows per irrep array
    Dimension rowspi_;
    /// Columns per irrep array
    Dimension colspi_;
    /// Name of the matrix
    std::string name_;
    /// Symmetry of this matrix (in most cases this will be 0 [totally symmetric])
    int symmetry_;

    /// Allocates matrix_
    void alloc();
    /// Release matrix_
    void release();

    /// Copies data from the passed matrix to this matrix_
    void copy_from(double***);

    /// allocate a block matrix -- analogous to libciomr's block_matrix
    static double** matrix(int nrow, int ncol);
    /// free a (block) matrix -- analogous to libciomr's free_block
    static void free(double** Block);

    void print_mat(const double* const* const a, int m, int n, std::string out) const;

    /// Numpy Shape
    std::vector<int> numpy_shape_;

   public:
    /// Default constructor, zeros everything out
    Matrix();
    /**
     * Constructor, zeros everything out, sets name_
     *
     * @param name Name of the matrix, used in saving and printing.
     */
    Matrix(const std::string& name, int symmetry = 0);
    /// copy reference constructor
    Matrix(const Matrix& copy);
    /// Explicit shared point copy constructor
    explicit Matrix(const SharedMatrix& copy);
    /// copy pointer constructor
    explicit Matrix(const Matrix* copy);
    /**
     * Constructor, sets up the matrix
     *
     * @param nirrep Number of blocks.
     * @param rowspi Array of length nirreps giving row dimensionality.
     * @param colspi Array of length nirreps giving column dimensionality.
     */
    Matrix(int nirrep, const int* rowspi, const int* colspi, int symmetry = 0);
    /**
     * Constructor, sets name_, and sets up the matrix
     *
     * @param name Name of the matrix.
     * @param nirrep Number of blocks.
     * @param rowspi Array of length nirreps giving row dimensionality.
     * @param colspi Array of length nirreps giving column dimensionality.
     */
    Matrix(const std::string& name, int nirrep, const int* rowspi, const int* colspi, int symmetry = 0);
    /**
     * Constructor, forms non-standard matrix.
     * @param nirrep Number of blocks.
     * @param rows Singular value. All blocks have same number of rows.
     * @param colspi Array of length nirreps. Defines blocking scheme for columns.
     */
    Matrix(int nirrep, int rows, const int* colspi);

    /**
     * Constructor, forms non-standard matrix.
     * @param nirrep Number of blocks.
     * @param rowspi Array of length nirreps. Defines blocking scheme for rows.
     * @param cols Singular value. All blocks have same number of columns.
     */
    Matrix(int nirrep, const int* rowspi, int cols);

    /**
     * Constructor, sets up the matrix
     * Convenience case for 1 irrep
     *
     * @param rows Row dimensionality.
     * @param cols Column dimensionality.
     */
    Matrix(int rows, int cols);
    /**
     * Constructor, sets up the matrix
     * Convenience case for 1 irrep
     *
     * @param name Name of the matrix.
     * @param rows Row dimensionality.
     * @param cols Column dimensionality.
     */
    Matrix(const std::string& name, int rows, int cols);

    /**
     * Contructs a Matrix from a dpdfile2
     *
     * @param inFile dpdfile2 object to replicate (must already be initialized).
     */
    Matrix(dpdfile2* inFile);

    /**
     * Constructor using Dimension objects to define order and dimensionality.
     *
     * @param name Name of the matrix.
     * @param rows Dimension object providing row information.
     * @param cols Dimension object providing column information.
     * @param symmetry overall symmetry of the data.
     */
    Matrix(const std::string& name, const Dimension& rows, const Dimension& cols, int symmetry = 0);

    /**
     * Constructor using Dimension objects to define order and dimensionality.
     *
     * @param rows Dimension object providing row information.
     * @param cols Dimension object providing column information.
     * @param symmetry overall symmetry of the data.
     */
    Matrix(const Dimension& rows, const Dimension& cols, int symmetry = 0);

    /// Destructor, frees memory
    virtual ~Matrix();

    /**
     * Initializes a matrix
     *
     * @param nirrep Number of blocks in this matrix.
     * @param rowspi Array of length nirreps giving row dimensionality.
     * @param colspi Array of length nirreps giving column dimensionality.
     * @param name Name of the matrix.
     * @param symmetry Overall symmetry of the data.
     */
    void init(int nirrep, const int* rowspi, const int* colspi, const std::string& name = "", int symmetry = 0);

    void init(const Dimension& rowspi, const Dimension& colspi, const std::string& name = "", int symmetry = 0);

    /// Creates an exact copy of the matrix and returns it.
    SharedMatrix clone() const;

    /**
     * Convenient creation function return SharedMatrix
     */
    static SharedMatrix create(const std::string& name, const Dimension& rows, const Dimension& cols);

    /**
     * @{
     * Copies data onto this
     * @param cp Object to copy from.
     */
    void copy(const SharedMatrix& cp);
    void copy(const Matrix& cp);
    void copy(const Matrix* cp);
    /** @} */

    /**
     * Horizontally concatenate matrices
     * @param mats std::vector of Matrix objects to concatenate
     */
    static SharedMatrix horzcat(const std::vector<SharedMatrix>& mats);

    /**
     * Vertically concatenate matrices
     * @param mats std::vector of Matrix objects to concatenate
     */
    static SharedMatrix vertcat(const std::vector<SharedMatrix>& mats);

    /**
    ** For a matrix of 3D vectors (ncol==3), rotate a set of points around an
    ** arbitrary axis.  Vectors are the rows of the matrix.
    **
    ** @param  axis  Vector3   : axis around which to rotate (need not be normalized)
    ** @param  phi   double    : magnitude of rotation in rad
    ** @param  Sn    bool      : if true, then also reflect in plane through origin and
    **                           perpendicular to rotation
    ** @returns SharedMatrix with rotated points (rows)
    */
    SharedMatrix matrix_3d_rotation(Vector3 axis, double phi, bool Sn);

    /// Copies data to the row specified. Assumes data is of correct length.
    void copy_to_row(int h, int row, double const* const data);

    enum SaveType { Full, SubBlocks, LowerTriangle };

    /**
     * @{
     * Load a matrix from a PSIO object from fileno with tocentry of size nso
     *
     * @param psio PSIO object to read with.
     * @param fileno File to read from.
     * @param tocentry Table of contents entry to use.
     * @param nso Number of orbitals to use to read in.
     * @returns true if loaded, false otherwise.
     */
    bool load(psi::PSIO* psio, size_t fileno, const std::string& tocentry, int nso);
    bool load(std::shared_ptr<psi::PSIO>& psio, size_t fileno, const std::string& tocentry, int nso);
    /** @} */

    /**
     * @{
     * Loads the block matrix from PSIO object with fileno and with the toc position of the name of the matrix
     *  The matrix must be correctly sized and named for this to work
     *
     * @param psio PSIO object to read with.
     * @param fileno File to read from.
     * @param savetype Save information suffixing point group label.
     */
    void load(psi::PSIO* const psio, size_t fileno, SaveType savetype = LowerTriangle);
    void load(std::shared_ptr<psi::PSIO>& psio, size_t fileno, SaveType savetype = LowerTriangle);
    /** @} */

    /**
     * Loads a block matrix from an ASCII file (see tests/mints3 for file format).
     *
     * @param filename Name of the file to read in.
     */
    void load(const std::string& filename);

    /**
     * Loads a matrix from an ASCII file. The matrix data resembles MPQC matrix printing
     * with additional size data.
     *
     * @param filename Name of the file to read in.
     */
    void load_mpqc(const std::string& filename);

    /**
     * @{
     * Saves the matrix in ASCII format to filename
     *
     * @param filename Name of the file to write to.
     * @param append Append to the file?
     * @param saveLowerTriangle Save only the lower triangle?
     * @param saveSubBlocks Save three index quantities denoting symmetry block (true), or convert to a full matrix and
     * save that (false)?
     */
    void save(const std::string& filename, bool append = true, bool saveLowerTriangle = true,
              bool saveSubBlocks = false);
    /** @} */

    /**
     * @{
     * Saves the block matrix to PSIO object with fileno and with the toc position of the name of the matrix
     *
     * @param psio PSIO object to write with.
     * @param fileno File to write to.
     * @param savetype Save information suffixing point group label.
     */
    void save(psi::PSIO* const psio, size_t fileno, SaveType savetype = LowerTriangle);
    void save(std::shared_ptr<psi::PSIO>& psio, size_t fileno, SaveType savetype = LowerTriangle);
    /** @} */

    /**
     * Set every element of matrix_ to val
     *
     * @param val Value to apply to entire matrix.
     */
    void set(double val);

    /**
     * Copies lower triangle tri to matrix_, calls tri_to_sq
     *
     * @param tri Lower triangle matrix to set to.
     */
    void set(const double* const tri);

    /**
     * @{
     * Copies sq to matrix_
     *
     * @param sq Double matrix to copy over.
     */
    void set(const double* const* const sq);
    /** @} */

    /**
     * @{
     * Copies sq to a specific irrep block of matrix_
     *
     * @param sq Double matrix to copy
     * @param irrep irrep block into which we copy
     */
    void set(const double* const* const sq, int irrep);
    /** @} */

    /**
     * Set a single element of matrix_
     *
     * @param h Subblock to address
     * @param m Row
     * @param n Column
     * @param val Value
     */
    void set(int h, int m, int n, double val) { matrix_[h][m][n] = val; }

    /**
     * Set a single element of matrix_
     *
     * @param m Row
     * @param n Column
     * @param val Value
     */
    void set(int m, int n, double val) { matrix_[0][m][n] = val; }

    /**
     * @{
     * Set the diagonal of matrix_ to vec
     *
     * @param vec Vector to apply to the diagonal.
     */
    void set_diagonal(const Vector* const vec);
    void set_diagonal(const Vector& vec);
    void set_diagonal(const std::shared_ptr<Vector>& vec);
    /** @} */

    /**
     * Returns a single element of matrix_
     *
     * @param h Subblock
     * @param m Row
     * @param n Column
     * @return value at position (h, m, n)
     */
    double get(const int& h, const int& m, const int& n) const { return matrix_[h][m][n]; }

    /**
     * Returns a single element of matrix_
     *
     * @param m Row
     * @param n Column
     * @return value at position (m, n)
     */
    double get(const int& m, const int& n) const { return matrix_[0][m][n]; }

    /**
     * Returns a single row of matrix_
     *
     * @param h Subblock
     * @param m Row
     * @return SharedVector object
     */
    SharedVector get_row(int h, int m);

    /**
     * Returns a single column of matrix_
     *
     * @param h Subblock
     * @param m Column
     * @return SharedVector object
     */
    SharedVector get_column(int h, int m);

    /**
     * Set a single row of matrix_
     *
     * @param h Subblock
     * @param m Row
     * @param vec SharedVector object to set the row to
     */
    void set_row(int h, int m, SharedVector vec);

    /**
     * Set a single column of matrix_
     *
     * @param h Subblock
     * @param m Column
     * @param vec SharedVector object to set the column to
     */
    void set_column(int h, int m, SharedVector vec);

    /**
     * Get a matrix block
     *
     * @param rows Rows slice
     * @param cols Columns slice
     * @return SharedMatrix object
     */
    SharedMatrix get_block(const Slice& rows, const Slice& cols);

    /**
     * Set a matrix block
     *
     * @param rows Rows slice
     * @param cols Columns slice
     * @param block the SharedMatrix object block to set
     */
    void set_block(const Slice& rows, const Slice& cols, SharedMatrix block);

    /**
     * Returns the double** pointer to the h-th irrep block matrix
     * NOTE: This method is provided for convenience in advanced
     * BLAS/LAPACK calls, and should be used with caution. In particular,
     * operations performed with these pointers should be scoped to avoid
     * erroneous alteration of the objects primitive data. Moreover,
     * the memory location/size of the double** obtained with this method
     * should NEVER be resized, moved, or freed.
     *
     * @param h Subblock
     * @return pointer to h-th subblock in block-matrix form
     */
    double** pointer(const int& h = 0) const { return matrix_[h]; }
    const double** const_pointer(const int& h = 0) const { return const_cast<const double**>(matrix_[h]); }

    /**
     * Returns the double* pointer to the h-th irrep block matrix
     * NOTE: This method is provided for convenience in advanced
     * BLAS/LAPACK calls, and should be used with caution. In particular,
     * operations performed with these pointers should be scoped to avoid
     * erroneous alteration of the objects primitive data. Moreover,
     * the memory location/size of the double* obtained with this method
     * should NEVER be resized, moved, or freed.
     *
     * @param h Subblock
     * @return pointer to h-th subblock in block-matrix form
     */
    double* get_pointer(const int& h = 0) const {
        if (rowspi_[h] * (size_t)colspi_[h] > 0)
            return &(matrix_[h][0][0]);
        else
            return nullptr;
    }
    const double* get_const_pointer(const int& h = 0) const {
        if (rowspi_[h] * (size_t)colspi_[h] > 0)
            return const_cast<const double*>(&(matrix_[h][0][0]));
        else
            return nullptr;
    }

    size_t size(const int& h = 0) const { return colspi_[h] * (size_t)rowspi_[h]; }

    /// apply_denominators a matrix to this
    void apply_denominator(const Matrix* const);
    /// apply_denominators a matrix to this
    void apply_denominator(const Matrix&);
    /// apply_denominators a matrix to this
    void apply_denominator(const SharedMatrix&);

    /**
     * Returns a copy of the current matrix.
     *
     * @returns the matrix
     */
    double** to_block_matrix() const;
    /**
     * Returns a copy of the current matrix.
     *
     * @returns the SharedMatrix
     */
    SharedMatrix to_block_sharedmatrix() const;
    /**
     * Returns a copy of the current matrix in lower triangle form.
     *
     * @returns the matrix
     */
    double* to_lower_triangle() const;

    /**
     * Sets the name of the matrix, used in print(...) and save(...)
     *
     * @param name New name to use.
     */
    void set_name(const std::string& name) { name_ = name; }

    /**
     * Gets the name of the matrix.
     */
    const std::string& name() const { return name_; }

    /// Python compatible printer
    void print_out() const { print("outfile"); }

    /**
     * Print the matrix using print_mat
     *
     * @param outfile File point to use, defaults to Psi4's outfile.
     * @param extra When printing the name of the 'extra' will be printing after the name.
     */
    void print(std::string outfile = "outfile", const char* extra = nullptr) const;

    /// Prints the matrix with atom and xyz styling.
    void print_atom_vector(std::string out_fname = "outfile");

    /**
     * Prints the matrix so that it can be copied and pasted into Mathematica easily.
     */
    void print_to_mathematica();

    /**
     * Print the matrix with corresponding eigenvalues below each column
     *
     * @param values Eigenvalues to print associated with eigenvectors.
     * @param out Where to print to, defaults to Psi4's outfile.
     */
    void eivprint(const Vector* const values, std::string out = "outfile");
    /// Print the matrix with corresponding eigenvalues below each column
    void eivprint(const Vector& values, std::string out = "outfile");
    /// Print the matrix with corresponding eigenvalues below each column
    void eivprint(const std::shared_ptr<Vector>& values, std::string out = "outfile");

    /// Returns the rows in irrep h
    int rowdim(const int& h = 0) const { return rowspi_[h]; }
    /// Returns the cols in irrep h
    int coldim(const int& h = 0) const { return colspi_[h]; }

    /// Returns the rows per irrep array
    const Dimension& rowspi() const { return rowspi_; }
    /// Returns the rows per irrep array
    int rowspi(const int& h) const { return rowdim(h); }
    /// Returns the columns per irrep array
    const Dimension& colspi() const { return colspi_; }
    /// Returns the columns per irrep array
    int colspi(const int& h) const { return coldim(h); }
    /// Returns the number of irreps
    const int& nirrep() const { return nirrep_; }

    /// Returns the total number of rows.
    int nrow() const {
        int rows = 0;
        for (int h = 0; h < nirrep(); ++h) rows += rowdim(h);
        return rows;
    }

    /// Returns the total number of columns.
    int ncol() const {
        int cols = 0;
        for (int h = 0; h < nirrep(); ++h) cols += coldim(h);
        return cols;
    }

    /// Returns the row size of the largest block.
    int max_nrow() const {
        int row = 0;
        for (int h = 0; h < nirrep(); ++h)
            if (row < rowdim(h)) row = rowdim(h);
        return row;
    }

    /// Returns the column size of the largest block.
    int max_ncol() const {
        int col = 0;
        for (int h = 0; h < nirrep(); ++h)
            if (col < coldim(h)) col = coldim(h);
        return col;
    }

    /**
     * Returns the overall symmetry of the matrix.
     * For a totally-symmetric matrix this will be 0.
     * The value returned is compatible with bitwise XOR (^) math.
     */
    const int& symmetry() const { return symmetry_; }

    /**
     * Symmetrizes the a gradient like matrix (N, 3) using information
     * from the given Molecule.
     */
    void symmetrize_gradient(std::shared_ptr<Molecule> mol);

    /**
     * Symmetrizes the a Hessian like matrix (3 * N, 3 * N) using information
     * from the given Molecule.
     */
    void symmetrize_hessian(std::shared_ptr<Molecule> mol);

    /// Set this to identity
    void identity();
    /// Zeros this out
    void zero();
    /// Zeros the diagonal
    void zero_diagonal();

    // Math routines
    /// Returns the trace of this
    double trace();
    /// Creates a new matrix which is the transpose of this
    SharedMatrix transpose();

    /// In place transposition
    void transpose_this();

    /// Adds a matrix to this
    void add(const Matrix* const);
    /// Adds a matrix to this
    void add(const Matrix&);
    /// Adds a matrix to this
    void add(const SharedMatrix&);

    /// Subtracts a matrix from this
    void subtract(const Matrix* const);
    /// Subtracts a matrix from this
    void subtract(const SharedMatrix&);
    /// Multiplies the two arguments and adds their result to this
    void accumulate_product(const Matrix* const, const Matrix* const);
    void accumulate_product(const SharedMatrix&, const SharedMatrix&);
    /// Scales this matrix
    void scale(double);
    /// Returns the sum of the squares of this
    double sum_of_squares();
    /// Returns the rms of this
    double rms();
    /// Returns the absolute maximum value
    double absmax();
    /// Add val to an element of this
    void add(int h, int m, int n, double val) {
#ifdef PSIDEBUG
        if (m > rowspi_[h] || n > colspi_[h ^ symmetry_]) {
            outfile->Printf("out of bounds: symmetry_ = %d, h = %d, m = %d, n = %d\n", symmetry_, h, m, n);

            throw PSIEXCEPTION("What are you doing, Rob?");
        }
#endif
        matrix_[h][m][n] += val;
    }
    /// Add val to an element of this
    void add(int m, int n, double val) {
#ifdef PSIDEBUG
        if (m > rowspi_[0] || n > colspi_[0 ^ symmetry_]) {
            outfile->Printf("out of bounds: symmetry_ = %d, h = %d, m = %d, n = %d\n", symmetry_, 0, m, n);

            return;
        }
#endif
        matrix_[0][m][n] += val;
    }

    void element_add_mirror() {
        for (int h = 0; h < nirrep_; ++h) {
            for (int i = 0; i < rowspi_[h]; ++i) {
                for (int j = 0; j < i; ++j) {
                    matrix_[h][i][j] = matrix_[h][j][i] = (matrix_[h][i][j] + matrix_[h][j][i]);
                }
            }
        }
    }

    /// Scale row m of irrep h by a
    void scale_row(int h, int m, double a);
    /// Scale column n of irrep h by a
    void scale_column(int h, int n, double a);

    /** Special function to add symmetry to a Matrix .
     *
     *  \param a Matrix to transform
     *  \param transformer The matrix returned by PetiteList::aotoso() that acts as the transformer
     */
    void apply_symmetry(const SharedMatrix& a, const SharedMatrix& transformer);

    /** Special function to remove symmetry from a matrix.
     *
     *  \param a symmetry matrix to transform
     *  \param transformer The matrix returned by PetiteList::sotoao() that acts as the transformer
     */
    void remove_symmetry(const SharedMatrix& a, const SharedMatrix& transformer);
    /** Performs a the transformation L^ F R. Result goes to this.
     *
     * \param L left transformation matrix (will be transposed)
     * \param F matrix to apply transformation to
     * \param R right transformation matrix (will not be transposed)
     */
    void transform(const SharedMatrix& L, const SharedMatrix& F, const SharedMatrix& R);

    /// @{
    /// Transform a by transformer save result to this
    void transform(const Matrix* const a, const Matrix* const transformer);
    void transform(const SharedMatrix& a, const SharedMatrix& transformer);
    /// @}

    /// @{
    /// Transform this by transformer
    void transform(const Matrix* const transformer);
    void transform(const SharedMatrix& transformer);
    /// @}

    /// @{
    /// Back transform a by transformer save result to this
    void back_transform(const Matrix* const a, const Matrix* const transformer);
    void back_transform(const SharedMatrix& a, const SharedMatrix& transformer);
    /// @}

    /// @{
    /// Back transform this by transformer
    void back_transform(const Matrix* const transformer);
    void back_transform(const SharedMatrix& transformer);
    /// @}

    /// Returns the vector dot product of this by rhs
    double vector_dot(const Matrix* const rhs);
    double vector_dot(const SharedMatrix& rhs);
    double vector_dot(const Matrix& rhs);

    /// @{
    /** General matrix multiply, saves result to this
     * \param transa Transpose the left matrix
     * \param transb Transpose the right matrix
     * \param alpha Prefactor for the matrix multiplication
     * \param a Left matrix
     * \param b Right matrix
     * \param beta Prefactor for the resulting matrix
     */
    void gemm(bool transa, bool transb, double alpha, const Matrix* const a, const Matrix* const b, double beta);
    void gemm(bool transa, bool transb, double alpha, const SharedMatrix& a, const SharedMatrix& b, double beta);
    void gemm(bool transa, bool transb, double alpha, const SharedMatrix& a, const Matrix& b, double beta);
    void gemm(bool transa, bool transb, double alpha, const Matrix& a, const SharedMatrix& b, double beta);
    /// @}

    /// @{
    /** Raw access to the underlying dgemm call. Saves result to this.
     */
    void gemm(const char& transa, const char& transb, const std::vector<int>& m, const std::vector<int>& n,
              const std::vector<int>& k, const double& alpha, const SharedMatrix& a, const std::vector<int>& lda,
              const SharedMatrix& b, const std::vector<int>& ldb, const double& beta, const std::vector<int>& ldc,
              const std::vector<unsigned long>& offset_a = std::vector<unsigned long>(),
              const std::vector<unsigned long>& offset_b = std::vector<unsigned long>(),
              const std::vector<unsigned long>& offset_c = std::vector<unsigned long>());
    void gemm(const char& transa, const char& transb, const int& m, const int& n, const int& k, const double& alpha,
              const SharedMatrix& a, const int& lda, const SharedMatrix& b, const int& ldb, const double& beta,
              const int& ldc, const unsigned long& offset_a = 0, const unsigned long& offset_b = 0,
              const unsigned long& offset_c = 0);
    /// @}

    /** Simple doublet GEMM with on-the-fly allocation
     * \param A The first matrix
     * \param B The second matrix
     * \param transA Transpose the first matrix
     * \param transB Transpose the second matrix
     */
    static SharedMatrix doublet(const SharedMatrix& A, const SharedMatrix& B, bool transA = false, bool transB = false);

    /** Simple triplet GEMM with on-the-fly allocation
     * \param A The first matrix
     * \param B The second matrix
     * \param C The third matrix
     * \param transA Transpose the first matrix
     * \param transB Transpose the second matrix
     * \param transC Transpose the third matrix
     */
    static SharedMatrix triplet(const SharedMatrix& A, const SharedMatrix& B, const SharedMatrix& C,
                                bool transA = false, bool transB = false, bool transC = false);

    /**
     * Simple AXPY call with support for irreps Y = a * X + Y
     * @param a Scaling parameter
     * @param X Matrix to be be added
     */
    void axpy(double a, SharedMatrix X);

    /** Summation collapse along either rows (0) or columns (1), always producing a column matrix
     * \param dim 0 (row sum) or 1 (col sum)
     * \return \sum_{i} M_{ij} => T_j if dim = 0 or
     *         \sum_{j} M_{ij} => T_i if dim = 1
     */
    SharedMatrix collapse(int dim = 0);

    /// @{
    /// Diagonalizes this, eigvectors and eigvalues must be created by caller.  Only for symmetric matrices.
    void diagonalize(Matrix* eigvectors, Vector* eigvalues, diagonalize_order nMatz = ascending);
    void diagonalize(SharedMatrix& eigvectors, std::shared_ptr<Vector>& eigvalues, diagonalize_order nMatz = ascending);
    void diagonalize(SharedMatrix& eigvectors, Vector& eigvalues, diagonalize_order nMatz = ascending);
    /// @}

    /// @{
    /// Diagonalizes this, applying supplied metric, eigvectors and eigvalues must be created by caller.  Only for
    /// symmetric matrices.
    void diagonalize(SharedMatrix& metric, SharedMatrix& eigvectors, std::shared_ptr<Vector>& eigvalues,
                     diagonalize_order nMatz = ascending);
    /// @}

    /// @{
    /// General SVD, such that A = USV. U, S, and V must be allocated by the caller.
    void svd(SharedMatrix& U, SharedVector& S, SharedMatrix& V);
    /// @}

    /// @{
    /// General SVD, such that A = USV. U, S, and V must be allocated by the caller.
    /// all M columns of U and all N rows of V**T are returned in the arrays U and VT;
    /// This assumes totally symmetric quantities.
    void svd_a(SharedMatrix& U, SharedVector& S, SharedMatrix& V);
    /// @}

    ///@{
    /// Matrices/Vectors U (m x k), S (k), V (k x n) to feed to Matrix::svd
    std::tuple<SharedMatrix, SharedVector, SharedMatrix> svd_temps();
    ///@}

    ///@{
    /// Matrices/Vectors U (m x m), S (k), V (n x n) to feed to Matrix::svd_a
    std::tuple<SharedMatrix, SharedVector, SharedMatrix> svd_a_temps();
    ///@}

    ///@{
    /// Matrix of size (m x n) which is the conditioned pseudoinverse of this (m x n)
    SharedMatrix pseudoinverse(double condition, int& nremoved);
    ///@}

    /*! Extract a conditioned orthonormal basis from this SPD matrix
     *  via canonical orthogonalization.
     *  @param delta the relative condition to maintain
     *  @return X, a SharedMatrix with m x m' dimension (m' < m if conditioning occurred)
     */
    SharedMatrix canonical_orthogonalization(double delta = 0.0, SharedMatrix eigvec = SharedMatrix());

    /*! Computes the fully pivoted partial Cholesky factorization of a real symmetric
     * positive semidefinite matrix A, to numerical precision \delta.
     *
     * The results L is dimpi x sigpi, where dimpi is the dimension of the original
     * matrix, and sigpi is the number of columns required to reach accuracy delta.
     * Sigpi is strictly less than or equal to dimpi
     *
     * The decomposition is of the form:
     *   A' \approx LL^T,
     *  where A is the original matrix, and L is the resultant
     *  Cholesky factor. The error matrix is:
     *   D = A - A'
     *
     * This function recursively computes the Schur complement to determine the
     * optimal ordering of columns.
     *
     * This algorithm requires up to 3 total core matrices of the size of the original
     * These are 1) the original, 2) The resultant, and 3) a temporary matrix
     *
     * \param delta maximum allowed error in the error matrix D, which always
     * occurs on the diagonal. Defaults to 0.0, in which case the numerically
     * exact factor is returned.
     * \param throw_if_negative bool, throw if pivot <= 0.0 is detected?
     * \return L, SharedMatrix, with rows of dimension dimpi and columns of
     * dimension sigpi
     */
    SharedMatrix partial_cholesky_factorize(double delta = 0.0, bool throw_if_negative = false);

    /*! Computes a low-rank factorization <P,N> such that PP'-NN' \approx A in an optimal sense in the 2-norm.
     * Columns of P,N are truncated after the singular values fall below delta
     * P contains columns corresponding to positive eigenvalues, N to columns corresponding to negative eigenvalues/
     * This is the real Hermitian-indefinite analog of partial Cholesky factorization.
     *
     * This algorithm requires memory equivalent to this matrix plus the equivalent eigendecompositon
     * call via DSYEV
     *
     * \param delta maximum allowed 2-norm of the error matrix D,
     * Defaults to 0.0, in which case the numerically
     * exact square root is returned.
     * \return P positive part of square root, with only significant columns included
     * \return N negative part of square root, with only significant columns included
     */
    std::pair<SharedMatrix, SharedMatrix> partial_square_root(double delta = 0.0);

    /*! Computes the Cholesky factorization of a real symmetric
     *  positive definite matrix A.
     *
     *  This is the block version of the algorithm, calling Level 3 BLAS.
     */
    void cholesky_factorize();

    /*! Computes the inverse of a real symmetric positive definite
     *  matrix A using the Cholesky factorization A = L*L**T
     *  computed by cholesky_factorize().
     */
    void invert();

    /*! Computes the inverse of a matrix using the LU factorization.
     *  This method inverts U and then computes inv(A) by solving the system
     *  inv(A)*L = inv(U) for inv(A).
     */
    void general_invert();

    /*! Computes the pseudo power of a real symmetric matrix
     *   A using eigendecomposition. This operation is uniquely defined
     *   for all symmetric matrices for integral alpha, and for
     *   all symmetric positive definite matrices for all alpha.
     *
     *   A fractional power of a Hermitian non-SPD matrix is not uniquely
     *   defined due to the ambiguity of the complex roots of unity, and
     *   will often be returned as NaN due to the formation of an imaginary
     *   root of an eigenvalue. Fractional powers should only be called for SPD
     *   matrices, and integral powers should always be specified with literals.
     *
     *   For negative powers, this operation is very sensitive to condition,
     *   and will discard eigenvectors corresponding to small eigenvalues which
     *   contribute to a condition number smaller than cutoff.
     *   The resultant power is actually a pseudo-power
     *
     *   \param alpha  The power to raise the matrix to
     *   \param cutoff The smallest absolute value of a condition number to allow to
     *   contribute in the formation of a negative power of A
     *   \returns a Dimension object with the remaining sizes. Can be used in a View.
     */
    Dimension power(double alpha, double cutoff = 1.0E-12);

    /*!
     * Computes the approximate
     * exponential of a general real square matrix via Pade
     * symmetric Pade approximation (orthonormality guaranteed)
     * (defaults to a 2 x 2 Pade table, with no
     * scaling or balancing)
     */
    void expm(int n = 2, bool scale = false);

    /// Swap rows i and j
    void swap_rows(int h, int i, int j);
    /// Swap cols i and j
    void swap_columns(int h, int i, int j);

    /*! Average off-diagonal elements */
    void hermitivitize();
    /*! Copy lower triangle to upper triangle */
    void copy_lower_to_upper();
    /*! Copy upper triangle to lower triangle */
    void copy_upper_to_lower();
    /*! Zero lower triangle */
    void zero_lower();
    /*! Zero upper triangle */
    void zero_upper();
    /*! Zero row */
    void zero_row(int h, int i);
    /*! Zero column */
    void zero_column(int h, int i);

    // Reference versions of the above functions
    /// Transform a by transformer save result to this
    void transform(const Matrix& a, const Matrix& transformer);
    /// Transform this by transformer
    void transform(const Matrix& transformer);
    /// Back transform a by transformer save result to this
    void back_transform(const Matrix& a, const Matrix& transformer);
    /// Back transform this by transformer
    void back_transform(const Matrix& transformer);

    /**
     * Expands the row dimension by one, and then orthogonalizes vector v against
     * the current rows, before setting the new row to the orthogonalized copy of v
     */
    bool add_and_orthogonalize_row(const SharedVector v);

    /*! @{
     * Assume this is a orthogonal matrix.  This function Gram-Schmidt
     * orthogonalizes a new vector v and adds it to matrix A. This must contain
     * a free row pointer for a new row.  Don't add orthogonalized v' if
     * norm(v') < NORM_TOL.
     *
     * Adapted from libqt's version by David Sherrill, Feb 1994
     *
     * \param rows current number of valid rows in this
     *             (this must have space for 'rows+1' row.)
     * \param v vector to add to A after it has been made orthogonal
     *             to rest of A
     *
     * \returns true if a vector is added, false otherwise
     */
    bool schmidt_add_row(int h, int rows, Vector& v) noexcept;
    bool schmidt_add_row(int h, int rows, double* v) noexcept;
    /// @}

    /*! Calls libqt schmidt function */
    void schmidt();

    /*! Schmidt orthogonalize this. S is the overlap matrix.
     *  n is the number of columns to orthogonalize. */
    void schmidt_orthog(SharedMatrix S, int n);

    /*! Schmidt orthogonalize this. You'll likely want to View this matrix afterwards
     *  using the result to obtain a properly sized Matrix.
     *  \param S overlap matrix.
     *  \param tol is the tolerance.
     *  \returns A Dimension object tell you how many were removed in each irrep.
     */
    Dimension schmidt_orthog_columns(SharedMatrix S, double tol, double* res = nullptr);

    /*!
     * Project out the row vectors in the matrix provided out of this matrix.
     * Assumes all matrices are C1 in nature. Future version will handle irreps.
     * Note: this is destroyed.
     *
     * \param v Matrix to project out
     */
    void project_out(Matrix& v);

    /// General matrix multiply, saves result to this
    void gemm(bool transa, bool transb, double alpha, const Matrix& a, const Matrix& b, double beta);
    /// Diagonalize a symmetric matrix. Eigvectors and eigvalues must be created by caller.
    void diagonalize(Matrix& eigvectors, Vector& eigvalues, int nMatz = 1);

    /// @{
    /// Retrieves the i'th irrep
    double** operator[](int i) { return matrix_[i]; }
    double& operator()(int i, int j) { return matrix_[0][i][j]; }
    const double& operator()(int i, int j) const { return matrix_[0][i][j]; }
    double& operator()(int h, int i, int j) { return matrix_[h][i][j]; }
    const double& operator()(int h, int i, int j) const { return matrix_[h][i][j]; }
    /// @}

    /// Writes this to the dpdfile2 given
    void write_to_dpdfile2(dpdfile2* outFile);

    /// @{
    /// Checks matrix equality.
    /// @param rhs Matrix to compare to.
    /// @returns true if equal, otherwise false.
    bool equal(const Matrix& rhs, double TOL = 1.0e-10);
    bool equal(const SharedMatrix& rhs, double TOL = 1.0e-10);
    bool equal(const Matrix* rhs, double TOL = 1.0e-10);
    /// @}

    /// @{
    /// Checks matrix equality, but allows rows to be in a different order.
    /// @param rhs Matrix to compare to.
    /// @returns true if equal, otherwise false.
    bool equal_but_for_row_order(const Matrix& rhs, double TOL = 1.0e-10);
    bool equal_but_for_row_order(const SharedMatrix& rhs, double TOL = 1.0e-10);
    bool equal_but_for_row_order(const Matrix* rhs, double TOL = 1.0e-10);
    /// @}

    /**
     * Adds accessability to the matrix shape for numpy
     */
    void set_numpy_shape(std::vector<int> shape) { numpy_shape_ = shape; }
    std::vector<int> numpy_shape() { return numpy_shape_; }

    /**
     * Rotates columns i and j in irrep h, by an angle theta
     * @param h - the irrep in which the rotation will be applied
     * @param i - the zero-based (within irrep) column number for i
     * @param j - the zero-based (within irrep) column number for j
     * @param theta - the angle (in radians) about which to rotate
     */
    void rotate_columns(int h, int i, int j, double theta);
    friend class Vector;
};

}  // namespace psi

#endif  // MATRIX_H
