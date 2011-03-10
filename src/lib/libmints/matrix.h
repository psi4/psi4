#ifndef _psi_src_lib_libmints_matrix_h_
#define _psi_src_lib_libmints_matrix_h_

#include <cstdio>
#include <string>
#include <cstring>

#include <libparallel/serialize.h>
#include <libdpd/dpd.h>

namespace boost {
template<class T> class shared_ptr;
// Forward declarations for boost.python used in the extract_subsets

namespace python{
       class tuple;
}}

namespace psi {

class PSIO;
class Matrix;
class Vector;
class SimpleVector;
class MatrixFactory;
class SimpleMatrix;
class Dimension;

extern FILE *outfile;

/*! \ingroup MINTS
 *  \class Matrix
 *  \brief Makes using matrices just a little earlier.
 *
 * Using a matrix factory makes creating these a breeze.
 */
class Matrix : public Serializable {
protected:
    /// Matrix data
    double ***matrix_;
    /// Number of irreps
    int nirrep_;
    /// Rows per irrep array
    int *rowspi_;
    /// Columns per irrep array
    int *colspi_;
    /// Name of the matrix
    std::string name_;
    /// Symmetry of this matrix (in most cases this will be 0 [totally symmetric])
    int symmetry_;

    /// Allocates matrix_
    void alloc();
    /// Release matrix_
    void release();

    /// Copies data from the passed matrix to this matrix_
    void copy_from(double ***);

    /// allocate a block matrix -- analogous to libciomr's block_matrix
    static double** matrix(int nrow, int ncol) {
        double** mat = (double**) malloc(sizeof(double*)*nrow);
        const size_t size = sizeof(double)*nrow*ncol;
        mat[0] = (double*) malloc(size);
        memset((void *)mat[0], 0, size);
        for(int r=1; r<nrow; ++r) mat[r] = mat[r-1] + ncol;
        return mat;
    }
    /// free a (block) matrix -- analogous to libciomr's free_block
    static void free(double** Block) {
        ::free(Block[0]);  ::free(Block);
    }

    void print_mat(const double *const *const a, int m, int n, FILE *out) const;

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
    explicit Matrix(const boost::shared_ptr<Matrix>& copy);
    /// copy pointer constructor
    explicit Matrix(const Matrix* copy);
    /**
     * Constructor, sets up the matrix
     *
     * @param nirreps Number of blocks.
     * @param rowspi Array of length nirreps giving row dimensionality.
     * @param colspi Array of length nirreps giving column dimensionality.
     */
    Matrix(int nirrep, const int *rowspi, const int *colspi, int symmetry = 0);
    /**
     * Constructor, sets name_, and sets up the matrix
     *
     * @param name Name of the matrix.
     * @param nirreps Number of blocks.
     * @param rowspi Array of length nirreps giving row dimensionality.
     * @param colspi Array of length nirreps giving column dimensionality.
     */
    Matrix(const std::string& name, int nirrep, const int *rowspi, const int *colspi, int symmetry = 0);
    /**
     * Constructor, forms non-standard matrix.
     * @param nirreps Number of blocks.
     * @param rows Singular value. All blocks have same number of rows.
     * @param colspi Array of length nirreps. Defines blocking scheme for columns.
     */
    Matrix(int nirrep, int rows, const int *colspi);

    /**
     * Constructor, forms non-standard matrix.
     * @param nirreps Number of blocks.
     * @param rowspi Array of length nirreps. Defines blocking scheme for rows.
     * @param cols Singular value. All blocks have same number of columns.
     */
    Matrix(int nirrep, const int* rowspi, int cols);

    /**
     * Constructor, sets up the matrix
     * Convenience case for 1 irrep
     * Note: You should be using SimpleMatrix
     *
     * @param rows Row dimensionality.
     * @param cols Column dimensionality.
     */
    Matrix(int rows, int cols);
    /**
     * Constructor, sets up the matrix
     * Convenience case for 1 irrep
     * Note: You should be using SimpleMatrix
     *
     * @param name Name of the matrix.
     * @param rows Row dimensionality.
     * @param cols Column dimensionality.
     */
    Matrix(const std::string&, int rows, int cols);

    /**
     * Contructs a Matrix from a dpdfile2
     *
     * @param inFile dpdfile2 object to replicate (must already be initialized).
     */
    Matrix(dpdfile2 *inFile);

    /**
     * Constructor using Dimension objects to define order and dimensionality.
     *
     * @param name Name of the matrix.
     * @param rows Dimension object providing row information.
     * @param cols Dimension object providing column information.
     */
    Matrix(const std::string& name, const Dimension& rows, const Dimension& cols, int symmetry = 0);

    /// Destructor, frees memory
    ~Matrix();

    /**
     * Initializes a matrix
     *
     * @param nirreps Number of blocks in this matrix.
     * @param rowspi Array of length nirreps giving row dimensionality.
     * @param colspi Array of length nirreps giving column dimensionality.
     */
    void init(int nirrep, const int *rowspi, const int *colspi, const std::string& name = "", int symmetry = 0);

    /// Creates an exact copy of the matrix and returns it.
    Matrix* clone() const;

    /**
     * @{
     * Copies data onto this
     * @param cp Object to copy from.
     */
    void copy(const boost::shared_ptr<Matrix>& cp);
    void copy(const Matrix& cp);
    void copy(const Matrix* cp);
    /** @} */

    enum SaveType {
        Full,
        SubBlocks,
        LowerTriangle
    };

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
    bool load(psi::PSIO* psio, unsigned int fileno, const std::string& tocentry, int nso);
    bool load(boost::shared_ptr<psi::PSIO>& psio, unsigned int fileno, const std::string& tocentry, int nso) {
        return load(psio.get(), fileno, tocentry, nso);
    }
    /** @} */

    /**
     * @{
     * Loads the block matrix from PSIO object with fileno and with the toc position of the name of the matrix
     *  The matrix must be correctly sized and named for this to work
     *
     * @param psio PSIO object to read with.
     * @param fileno File to read from.
     * @param saveSubBlocks Save information suffixing point group label.
     */
    void load(psi::PSIO* const psio, unsigned int fileno, SaveType savetype=LowerTriangle);
    void load(boost::shared_ptr<psi::PSIO>& psio, unsigned int fileno, SaveType savetype=LowerTriangle) {
        load(psio.get(), fileno, savetype);
    }
    /** @} */

    /**
     * @{
     * Saves the matrix in ASCII format to filename
     *
     * @param filename Name of the file to write to.
     * @param append Append to the file?
     * @param saveLowerTriangle Save only the lower triangle?
     * @param saveSubBlocks Save three index quantities denoting symmetry block (true), or convert to a full matrix and save that (false)?
     */
    void save(const std::string& filename, bool append=true, bool saveLowerTriangle = true, bool saveSubBlocks=false);
    /** @} */

    /**
     * @{
     * Saves the block matrix to PSIO object with fileno and with the toc position of the name of the matrix
     *
     * @param psio PSIO object to write with.
     * @param fileno File to write to.
     * @param saveSubBlocks Save information suffixing point group label.
     */
    void save(psi::PSIO* const psio, unsigned int fileno, SaveType savetype=LowerTriangle);
    void save(boost::shared_ptr<psi::PSIO>& psio, unsigned int fileno, SaveType savetype=LowerTriangle) {
        save(psio.get(), fileno, savetype);
    }
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
    void set(const double * const tri);

    /**
     * @{
     * Copies sq to matrix_
     *
     * @param sq Double matrix to copy over.
     */
    void set(const double * const * const sq);
    /** @} */

    /**
     * @{
     * Copies sq to matrix_
     *
     * @param sq SimpleMatrix object to set this matrix to.
     */
    void set(const SimpleMatrix * const sq);
    void set(const boost::shared_ptr<SimpleMatrix>& sq);
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
     * @{
     * Set the diagonal of matrix_ to vec
     *
     * @param vec Vector to apply to the diagonal.
     */
    void set(const Vector * const vec);
    void set(const Vector& vec);
    void set(const boost::shared_ptr<Vector>& vec);
    /** @} */

    /**
     * Returns a single element of matrix_
     *
     * @param h Subblock
     * @param m Row
     * @param n Column
     * @returns value at position (h, m, n)
     */
    double get(int h, int m, int n) const { return matrix_[h][m][n]; }

    /**
     * Python wrapper for get
     */
    double pyget(const boost::python::tuple& key);
    /**
     * Python wrapper for set
     */
    void pyset(const boost::python::tuple& key, double value);

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
     * @returns pointer to h-th subblock in block-matrix form
     */
    double** pointer(int h = 0) const { return matrix_[h]; }
    const double** const_pointer(int h=0) const { return const_cast<const double**>(matrix_[h]); }

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
     * @returns pointer to h-th subblock in block-matrix form
     */
    double* get_pointer(int h = 0) const {
        if(rowspi_[h]*colspi_[h] > 0)
           return &(matrix_[h][0][0]);
        else
           return 0;}
    const double* get_const_pointer(int h=0) const {
        if(rowspi_[h]*colspi_[h] > 0)
           return const_cast<const double*>(&(matrix_[h][0][0]));
        else
           return 0;}

    /**
     * Returns a copy of the current matrix.
     *
     * @returns the matrix
     */
    double **to_block_matrix() const;
    /**
     * Returns a copy of the current matrix in lower triangle form.
     *
     * @returns the matrix
     */
    double *to_lower_triangle() const;

    /**
     * Converts this to a full non-symmetry-block matrix
     *
     * @returns The SimpleMatrix copy of the current matrix.
     */
    SimpleMatrix *to_simple_matrix() const;

    /**
     * Sets the name of the matrix, used in print(...) and save(...)
     *
     * @param name New name to use.
     */
    void set_name(const std::string& name) {
        name_ = name;
    }

    /**
     * Gets the name of the matrix.
     */
    std::string name() const {
        return name_;
    }

    /// Python compatible printer
    void print_out() const { print(outfile); }

    /**
     * Print the matrix using print_mat
     *
     * @param outfile File point to use, defaults to Psi4's outfile.
     * @param extra When printing the name of the 'extra' will be printing after the name.
     */
    void print(FILE *out = outfile, const char *extra=NULL) const;

    /**
     * Print the matrix with corresponding eigenvalues below each column
     *
     * @param values Eigenvalues to print associated with eigenvectors.
     * @param out Where to print to, defaults to Psi4's outfile.
     */
    void eivprint(const Vector * const values, FILE *out = outfile);
    /// Print the matrix with corresponding eigenvalues below each column
    void eivprint(const Vector& values, FILE *out = outfile);
    /// Print the matrix with corresponding eigenvalues below each column
    void eivprint(const boost::shared_ptr<Vector>& values, FILE *out = outfile);

    /// Returns the rows in irrep h
    int rowdim(int h = 0) const {
        return rowspi_[h];
    }
    /// Returns the cols in irrep h
    int coldim(int h = 0) const {
        return colspi_[h];
    }

    /// Returns the rows per irrep array
    int *rowspi() const {
        return rowspi_;
    }
    /// Returns the columns per irrep array
    int *colspi() const {
        return colspi_;
    }
    /// Returns the number of irreps
    int nirrep() const {
        return nirrep_;
    }

    /// Returns the total number of rows.
    int nrow() const {
        int rows = 0;
        for (int h=0; h<nirrep(); ++h)
            rows += rowdim(h);
        return rows;
    }

    /// Returns the total number of columns.
    int ncol() const {
        int cols = 0;
        for (int h=0; h<nirrep(); ++h)
            cols += coldim(h);
        return cols;
    }

    /**
     * Returns the overall symmetry of the matrix.
     * For a totally-symmetric matrix this will be 0.
     * The value returned is compatible with bitwise XOR (^) math.
     */
    int symmetry() const {
        return symmetry_;
    }

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
    Matrix *transpose();

    /// Adds a matrix to this
    void add(const Matrix* const);
    /// Adds a matrix to this
    void add(const Matrix&);
    /// Adds a matrix to this
    void add(const boost::shared_ptr<Matrix>&);

    /// Subtracts a matrix from this
    void subtract(const Matrix* const);
    /// Subtracts a matrix from this
    void subtract(const boost::shared_ptr<Matrix>&);
    /// Multiplies the two arguments and adds their result to this
    void accumulate_product(const Matrix* const, const Matrix* const);
    void accumulate_product(const boost::shared_ptr<Matrix>&, const boost::shared_ptr<Matrix>&);
    /// Scales this matrix
    void scale(double);
    /// Returns the sum of the squares of this
    double sum_of_squares();
    /// Returns the rms of this
    double rms();
    /// Add val to an element of this
    void add(int h, int m, int n, double val) {
        if (m > rowspi_[h] || n > colspi_[h^symmetry_]) {
            fprintf(outfile, "out of bounds: symmetry_ = %d, h = %d, m = %d, n = %d\n",
                    symmetry_, h, m, n);
            fflush(outfile);
            return;
        }
        matrix_[h][m][n] += val;
    }
    /// Scale row m of irrep h by a
    void scale_row(int h, int m, double a);
    /// Scale column n of irrep h by a
    void scale_column(int h, int n, double a);

    /** Special function to transform a SimpleMatrix (no symmetry) into
     *  a symmetry matrix.
     *
     *  \param a SimpleMatrix to transform
     *  \param transformer Matrix from PetiteList to act as transformer
     */
    void transform(const boost::shared_ptr<SimpleMatrix>& a, const boost::shared_ptr<Matrix>& transformer);

    /** Performs a the transformation L^ F R. Result goes to this.
     *
     * \param L left transformation matrix (will be transposed)
     * \param F matrix to apply transformation to
     * \param R right transformation matrix (will not be transposed)
     */
    void transform(const boost::shared_ptr<Matrix>& L,
                   const boost::shared_ptr<Matrix>& F,
                   const boost::shared_ptr<Matrix>& R);

    /// @{
    /// Transform a by transformer save result to this
    void transform(const Matrix* const a, const Matrix* const transformer);
    void transform(const boost::shared_ptr<Matrix>& a, const boost::shared_ptr<Matrix>& transformer);
    /// @}

    /// @{
    /// Transform this by transformer
    void transform(const Matrix* const transformer);
    void transform(const boost::shared_ptr<Matrix>& transformer);
    /// @}

    /// @{
    /// Back transform a by transformer save result to this
    void back_transform(const Matrix* const a, const Matrix* const transformer);
    void back_transform(const boost::shared_ptr<Matrix>& a, const boost::shared_ptr<Matrix>& transformer);
    /// @}

    /// @{
    /// Back transform this by transformer
    void back_transform(const Matrix* const transformer);
    void back_transform(const boost::shared_ptr<Matrix>& transformer);
    /// @}

    /// Returns the vector dot product of this by rhs
    double vector_dot(const Matrix* const rhs);
    double vector_dot(const boost::shared_ptr<Matrix>& rhs);
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
    void gemm(bool transa, bool transb, double alpha, const boost::shared_ptr<Matrix>& a, const boost::shared_ptr<Matrix>& b, double beta);
    void gemm(bool transa, bool transb, double alpha, const boost::shared_ptr<Matrix>& a, const Matrix& b, double beta);
    void gemm(bool transa, bool transb, double alpha, const Matrix& a, const boost::shared_ptr<Matrix>& b, double beta);
    /// @}

    /// @{
    /// Diagonalize this places eigvectors and eigvalues must be created by caller.
    void diagonalize(Matrix* eigvectors, Vector* eigvalues, int nMatz = 1);
    void diagonalize(boost::shared_ptr<Matrix>& eigvectors, boost::shared_ptr<Vector>& eigvalues, int nMatz = 1);
    void diagonalize(boost::shared_ptr<Matrix>& eigvectors, Vector& eigvalues, int nMatz = 1);
    /// @}

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

    /*! Copy lower triangle to upper triangle */
    void copy_lower_to_upper();
    /*! Copy upper triangle to lower triangle */
    void copy_upper_to_lower();

    // Reference versions of the above functions
    /// Transform a by transformer save result to this
    void transform(const Matrix& a, const Matrix& transformer);
    /// Transform this by transformer
    void transform(const Matrix& transformer);
    /// Back transform a by transformer save result to this
    void back_transform(const Matrix& a, const Matrix& transformer);
    /// Back transform this by transformer
    void back_transform(const Matrix& transformer);

    /// General matrix multiply, saves result to this
    void gemm(bool transa, bool transb, double alpha, const Matrix& a, const Matrix& b, double beta);
    /// Diagonalize. Eigvectors and eigvalues must be created by caller.
    void diagonalize(Matrix& eigvectors, Vector& eigvalues, int nMatz = 1);

    /// @{
    /// Retrieves the i'th irrep
    double** operator[](int i) { return matrix_[i]; }
    double& operator()(int h, int i, int j) { return matrix_[h][i][j]; }
    const double& operator()(int h, int i, int j) const { return matrix_[h][i][j]; }
    /// @}

    // Serializable pure virtual functions:
    void send(Communicator* comm);
    void recv(Communicator* comm);
    void bcast(Communicator* comm, int broadcaster);

    /// Writes this to the dpdfile2 given
    void write_to_dpdfile2(dpdfile2 *outFile);

    /// @{
    /// Checks matrix equality.
    /// @param rhs Matrix to compare to.
    /// @returns true if equal, otherwise false.
    bool equal(const Matrix& rhs);
    bool equal(const boost::shared_ptr<Matrix>& rhs);
    bool equal(const Matrix* rhs);
    /// @}

    friend class Vector;
};

/*! \ingroup MINTS
 *  \class SimpleMatrix
 *  \brief Simple matrix class. Not symmetry blocked.
 */
class SimpleMatrix
{
protected:
    /// Matrix data
    double **matrix_;
    /// Number of rows and columns
    int rows_, cols_;
    /// Nae of the matrix
    std::string name_;

    /// Allocates matrix_
    void alloc();
    /// Releases matrix_
    void release();

    /// Copies data from the passed matrix to this matrix_
    void copy_from(double **);

    /// allocate a block matrix -- analogous to libciomr's block_matrix
    static double** matrix(int nrow, int ncol) {
        double** mat = (double**) malloc(sizeof(double*)*nrow);
        const size_t size = sizeof(double)*nrow*ncol;
        mat[0] = (double*) malloc(size);
        //bzero((void*)mat[0],size);
        memset((void *)mat[0], '\0', size);
        for(int r=1; r<nrow; ++r) mat[r] = mat[r-1] + ncol;
        return mat;
    }
    /// free a (block) matrix -- analogous to libciomr's free_block
    static void free(double** Block) {
        ::free(Block[0]);  ::free(Block);
    }

public:
    /// Default constructor, zeros everything out
    SimpleMatrix();
    /// Constructor, zeros everything out, sets name_
    SimpleMatrix(std::string name);
    /// Explicit copy reference constructor
    SimpleMatrix(const SimpleMatrix& copy);
    /// Explicit copy pointer constructor
    explicit SimpleMatrix(const SimpleMatrix* copy);
    explicit SimpleMatrix(boost::shared_ptr<SimpleMatrix> copy);
    /// Constructor, sets up the matrix
    SimpleMatrix(int nrow, int ncol);
    /// Constructor, sets name_, and sets up the matrix
    SimpleMatrix(std::string name, int nrow, int ncol);
    /// Converts Matrix reference to SimpleMatrix
    SimpleMatrix(const Matrix& copy);
    /// Converts Matrix pointer to SimpleMatrix
    SimpleMatrix(const Matrix* copy);

    /// Constructor using Dimension's
    SimpleMatrix(std::string& name, const Dimension& nrow, const Dimension& ncol);

    /// Destructor, frees memory
    ~SimpleMatrix();

    /// Initializes a matrix
    void init(int rowspi, int colspi, std::string name = "");

    /// Creates an exact copy of the matrix and returns it.
    SimpleMatrix* clone() const;
    /// Copies cp's data onto this
    void copy(SimpleMatrix* cp);
    void copy(boost::shared_ptr<SimpleMatrix> cp);

    /// Set every element of this to val
    void set(double val);
    /// Copies lower triangle tri to matrix_
    void set(const double *tri);
    /// Set a single element of matrix_
    void set(int m, int n, double val) { matrix_[m][n] = val; }
    void set(SimpleVector *vec);
    void set(boost::shared_ptr<SimpleVector> vec);
    void set(double **mat);

    double** pointer() const { return matrix_; }
    const double** const_pointer() const { return const_cast<const double**>(matrix_); }

    /// Sets the diagonal of matrix_ to vec
    double get(int m, int n) { return matrix_[m][n]; }
    /// Returns matrix_
    double **to_block_matrix() const;
    /// Sets the name of the matrix
    void set_name(std::string name) { name_ = name; }
    /// Returns the pointer to the matrix
    double* ptr() { return &matrix_[0][0]; }

    /// Python compatible printer
    void print_out() { print(outfile); }

    /// Prints the matrix with print_mat
    void print(FILE *out = outfile);

    /// Prints the matrix with atom and xyz styling.
    void print_atom_vector(FILE *out = outfile);

    /// Print the matrix with corresponding eigenvalues below each column
    void eivprint(SimpleVector *values, FILE *out = outfile);
    void eivprint(boost::shared_ptr<SimpleVector> values, FILE *out = outfile);
    /// The number of rows
    int nrow() const { return rows_; }
    /// The number of columns
    int ncol() const { return cols_; }
    /// Set matrix to identity
    void identity();
    /// Zero out the matrix
    void zero();
    /// Zero out the diagonal
    void zero_diagonal();

    /// Returns the trace of this
    double trace() const;
    /*! Computes the inverse of a real symmetric positive definite
     *  matrix A using the Cholesky factorization A = L*L**T
     *  computed by cholesky_factorize().
     */
    void invert();
    /// Create a new SimpleMatrix which is the transpose of this
    SimpleMatrix *transpose();
    /// Transpose of this
    void transpose_this();

    void copy_lower_to_upper();
    void copy_upper_to_lower();

    void lu_factorize();

    /// Add a matrix to this
    void add(const SimpleMatrix*);
    void add(boost::shared_ptr<SimpleMatrix>);
    /// Subtracts a matrix from this
    void subtract(const SimpleMatrix*);
    void subtract(boost::shared_ptr<SimpleMatrix>);
    /// Multiples the two arguments and adds their result to this
    void accumulate_product(const SimpleMatrix*, const SimpleMatrix*);
    void accumulate_product(boost::shared_ptr<SimpleMatrix>, boost::shared_ptr<SimpleMatrix>);
    /// Scales this matrix
    void scale(double);
    /// Returns the sum of the squares of this
    double sum_of_squares();
    /// Add val to an element of this
    void add(int m, int n, double val) { matrix_[m][n] += val; }
    /// Scale row m by a
    void scale_row(int m, double a);
    /// Scale column n by a
    void scale_column(int n, double a);
    /// Transform a by transformer save result to this
    void transform(SimpleMatrix* a, SimpleMatrix* transformer);
    void transform(boost::shared_ptr<SimpleMatrix> a, boost::shared_ptr<SimpleMatrix> transformer);
    /// Transform this by transformer
    void transform(SimpleMatrix* transformer);
    void transform(boost::shared_ptr<SimpleMatrix> transformer);
    /// Back transform a by transformer save result to this
    void back_transform(SimpleMatrix* a, SimpleMatrix* transformer);
    void back_transform(boost::shared_ptr<SimpleMatrix> a, boost::shared_ptr<SimpleMatrix> transformer);
    /// Back transform this by transformer
    void back_transform(SimpleMatrix* transformer);
    void back_transform(boost::shared_ptr<SimpleMatrix> transformer);

    /// Return the vector dot product of rhs by this
    double vector_dot(SimpleMatrix* rhs);
    double vector_dot(boost::shared_ptr<SimpleMatrix> rhs);

    /// Element add mirror.
    /// Performs (i, j) = (j, i) = (i, j) + (j, i)
    void element_add_mirror();

    /// General matrix multiply, saves result to this
    void gemm(bool transa, bool transb, double alpha, const SimpleMatrix* a, const SimpleMatrix* b, double beta);
    void gemm(bool transa, bool transb, double alpha, boost::shared_ptr<SimpleMatrix> a, boost::shared_ptr<SimpleMatrix> b, double beta);
    /// Diagonalize this, eigvector and eigvalues must be created by caller.
    void diagonalize(SimpleMatrix* eigvectors, SimpleVector* eigvalues, int sort=1);
    void diagonalize(boost::shared_ptr<SimpleMatrix> eigvectors, boost::shared_ptr<SimpleVector> eigvalues, int sort=1);

    /// Saves the block matrix to PSIO object with fileno and with the toc position of the name of the matrix
    void save(psi::PSIO* psio, unsigned int fileno);
    void save(psi::PSIO& psio, unsigned int fileno);
    void save(boost::shared_ptr<psi::PSIO> psio, unsigned int fileno);
    void load(boost::shared_ptr<psi::PSIO> psio, unsigned int fileno);

    /// Saves the matrix in ASCII format to filename
    void save(const char *filename, bool append=true, bool saveLowerTriangle = true);
    /// Saves the matrix in ASCII format to filename
    void save(std::string filename, bool append=true, bool saveLowerTriangle = true);

    /// Retrieves the i'th row
    double* operator[](int i) { return matrix_[i]; }
    double& operator()(int i, int j) { return matrix_[i][j]; }
    const double& operator()(int i, int j) const { return matrix_[i][j]; }

    /// Swap rows i and j
    void swap_rows(int i, int j);
    /// Swap cols i and j
    void swap_columns(int i, int j);

    friend class Matrix;
};

/*! \ingroup MINTS
 *  \class View
 *  \brief Provides a view to a region of the matrix.
 */
class View
{
protected:
    /// Matrix we are viewing.
    boost::shared_ptr<Matrix> matrix_;
    /// Number of irreps
    int nirrep_;
    /// Starting offsets in matrix_;
    int *row_offset_per_irrep_;
    int *col_offset_per_irrep_;
    /// Number of rows we are viewing
    int *rows_per_irrep_;
    int *cols_per_irrep_;

private:
    View();  // No default constructor
    View(const View& );  // No copy constructor

public:
    virtual ~View();

    /** Constructor, assumes offsets for each irrep is 0
     *  @param nirrep Number of irreps
     *  @param rows How many rows per irrep are we interested in
     *  @param cols How many cols per irrep are we interested in
     */
    View(int nirrep, int *rows, int *cols);
    /** Constructor, user provides offsets and dimensions.
     *  @param nirrep Number of irreps
     *  @param rows How many rows per irrep
     *  @param cols How many cols per irrep
     *  @param row_offsets Row offset per irrep
     *  @param col_offsets Column offset per irrep
     */
    View(int nirrep, int *rows, int *cols, int *row_offsets, int *col_offsets);
    /** Constructor, user provides a Matrix to view and desired row count
     *  @param matrix Matrix we want to view, View obtains nirrep from it
     *  @param rows How many rows per irrep
     *  @param cols How many cols per irrep
     */
    View(boost::shared_ptr<Matrix> matrix, int *rows, int *cols);
    /** Constructor, user provides a Matrix to view and desired row count
     *  @param matrix Matrix we want to view, View obtains nirrep from it
     *  @param rows How many rows per irrep
     *  @param cols How many cols per irrep
     *  @param row_offsets Row offset per irrep
     *  @param col_offsets Column offset per irrep
     */
    View(boost::shared_ptr<Matrix> matrix, int *rows, int *cols, int *row_offsets, int *col_offsets);

    /** Operator () overload. Creates a new Matrix that only contains the view.
     *  @return New Matrix containing the view.
     */
    boost::shared_ptr<Matrix> operator()();

    /** Set the Matrix that we should be viewing.
     *  @param matrix Matrix to view.
     *  @return The old Matrix we were viewing.
     */
    boost::shared_ptr<Matrix> view(boost::shared_ptr<Matrix> matrix);
};


typedef boost::shared_ptr<Matrix>       SharedMatrix;
typedef boost::shared_ptr<SimpleMatrix> SharedSimpleMatrix;

}

#endif // MATRIX_H
