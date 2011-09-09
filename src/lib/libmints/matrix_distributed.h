#ifndef _psi_src_lib_libmints_matrix_distributed_h_
#define _psi_src_lib_libmints_matrix_distributed_h_

#include <cstdio>
#include <string>
#include <cstring>

#include <libparallel/parallel.h>
#include <libdpd/dpd.h>
#include "matrix.h"
#include "matrix_block.h"

namespace boost {
template<class T> class shared_ptr;
// Forward declarations for boost.python used in the extract_subsets

namespace python {
       class tuple;
}}

namespace psi {

class PSIO;
class Matrix;
//class Vector;
//class SimpleVector;
//class MatrixFactory;
//class SimpleMatrix;
class Dimension;

extern FILE *outfile;

typedef boost::shared_ptr<Matrix> SharedMatrix;

/*! \ingroup MINTS
 *  \class Distributed_Matrix
 *  \brief Makes using distributed matrices just a little earlier.
 *
 * Using a matrix factory makes creating these a breeze.
 * This has NOT been interfaced with the matrix factory yet.
 */
class Distributed_Matrix : public Block
#if HAVE_MADNESS == 1
        , public madness::WorldObject<Distributed_Matrix>
#endif
{
protected:

    /// The name of the distributed matrix
    std::string name_;
    /// The total number of blocks
    int nblocks_;
    /// The number of rows
    int nrows_;
    /// The number of columns
    int ncols_;
    /// The total number of sub_blocks
    int tsb_;
    /// The total number of elements
    int telements_;
    /// The offset to the divided row block;
    std::vector<int> block_roffset_;
    /// The offset to the divided column block;
    std::vector<int> block_coffset_;
    /// The size of each divided row block;
    std::vector<int> block_csize_;
    /// The size of each divided column block;
    std::vector<int> block_rsize_;

    /// Process ID
    int me_;
    /// Number of MPI Processes
    int nprocs_;
    /// Number of threads
    int nthreads_;
    /// Communicator Type
    std::string comm_;

#if HAVE_MADNESS == 1
    /// Madness world object
    SharedMadWorld madworld_;
    /// Mutex for printing
    SharedMutex print_mutex_;
    madness::Void clear();

public:
    /// Default constructor: clears everything out
    Distributed_Matrix();
    /// Constructs a distributed (rows x cols) matrix using the default distribution
    Distributed_Matrix(const int &rows, const int &cols);
    /// Constructs a distributed (rows x cols) matrix using the default distribution
    Distributed_Matrix(const std::string &name, const int &rows, const int &cols);
    /// copy reference constructor
    Distributed_Matrix(const Distributed_Matrix &copy);
    //    /// Explicit shared point copy constructor
    //    explicit Distributed_Matrix(const boost::shared_ptr<Distributed_Matrix> copy);
    //    /// Explicit copy pointer constructor
    //    explicit Distributed_Matrix(const Distributed_Matrix* copy);

    /// Return the number of rows in the distributed matrix
    int nrows() const {return nrows_;}
    /// Return the number of columns in the distributed matrix
    int ncols() const {return ncols_;}
    /// Return the name of the matrix
    std::string name() const {return name_;}

    ~Distributed_Matrix() { Communicator::world->sync(); }

    /// Initializes the distributed matrix
    madness::Void common_init(const int &rows, const int &cols);

    /// Print a give block
    madness::Void print_block(const int &i) const;
    /// Prints all of the blocks (proc 0 does the printing)
    void print_all_blocks() const;
    /// Prints all of the sub_blocks
    void print_all_sblocks() const;
    /// Prints a given sub_block
    madness::Void print_sblock(const int &i) const;

    /// Set the distributed matrix to the identity
    madness::Void identity();

    /**
     * Return the distributed matrix value
     *
     * @param row The row where the value is.
     * @param col The column where the value is.
     */
    madness::Future<double> get_val(const int &row, const int &col);

    /**
     * Set the matrix element (i,j) equal to val
     *
     * @param row The row where the value is to be put.
     * @param col The column where the value is to be put.
     * @param val The value to set the matrix element.
     */
    madness::Void set_val(const int &row, const int &col, const double &val);

    /// Overloaded equal operators
    Distributed_Matrix& operator= (const Distributed_Matrix &mat);
    Distributed_Matrix& operator= (const Distributed_Matrix *mat);
    Distributed_Matrix& operator= (const boost::shared_ptr<Distributed_Matrix> mat);
    madness::Void operator= (const double &val);
    /// Overloaded += operator
    madness::Void operator+= (const Distributed_Matrix &mat);
    /// Overloaded + operator
    Distributed_Matrix operator+ (const Distributed_Matrix &rhs);
    /// Overloaded == operator
    bool operator ==(const Distributed_Matrix &mat);
    /// Overloaded != operator
    bool operator !=(const Distributed_Matrix &mat);    

    /// Return the owner of the block
    int owner(const int &blk) const { return blk%nprocs_; }

    /// Prints a vector in matrix format
    madness::Void print_mat(const std::vector<double> &a, const int &m,
                            const int &n, const std::string &name) const;

    /// Set the name of the matrix
    madness::Void set_name(const std::string &nm) { name_ = nm; }

    /// Fill the entire matrix with the given value
    madness::Void fill(const double &val);

    /// Multiply two matrices
    Distributed_Matrix operator* (const Distributed_Matrix &rhs);

    /// This zeros out a distributed matrix
    madness::Void zero() { *this = 0.0; }

#else

public:
    /// Default constructor: clears everything out
    Distributed_Matrix()
    { throw PSIEXCEPTION("Distributed matrix only works with MADNESS.\n"); }
    /// Constructs a distributed (rows x cols) matrix using the default distribution
    Distributed_Matrix(const int &rows, const int &cols)
    { throw PSIEXCEPTION("Distributed matrix only works with MADNESS.\n"); }
    /// Constructs a distributed (rows x cols) matrix using the default distribution
    Distributed_Matrix(const std::string &name, const int &rows, const int &cols)
    { throw PSIEXCEPTION("Distributed matrix only works with MADNESS.\n"); }
    /// copy reference constructor
    Distributed_Matrix(const Distributed_Matrix &copy)
    { throw PSIEXCEPTION("Distributed matrix only works with MADNESS.\n"); }

    ~Distributed_Matrix() { Communicator::world->sync(); }

    int nrows() const {return nrows_;}
    int ncols() const {return ncols_;}
    std::string name() const {return name_;}
    void common_init(const int &rows, const int &cols);
    void print_block(const int &i);
    void print_all_blocks();
    void print_all_sblocks();
    void print_sblock(const int &i);
    void identity();
    double get_val(const int &row, const int &col);
    void set_val(const int &row, const int &col, const double &val);
    Distributed_Matrix& operator= (const Distributed_Matrix &mat);
    Distributed_Matrix& operator= (const Distributed_Matrix *mat);
    Distributed_Matrix& operator= (const boost::shared_ptr<Distributed_Matrix> mat);
    void operator= (const double &val);
    void operator+= (const Distributed_Matrix &mat);
    Distributed_Matrix operator+ (const Distributed_Matrix &rhs);
    bool operator ==(const Distributed_Matrix &mat);
    bool operator !=(const Distributed_Matrix &mat);
    int owner(const int &blk) { return blk%nprocs_; }
    void print_mat(const std::vector<double> &a, const int &m,
                   const int &n, const std::string &name);
    void set_name(const std::string &nm) { name_ = nm; }
    void fill(const double &val);
    Distributed_Matrix operator* (const Distributed_Matrix &rhs);
    void zero() { *this = 0.0; }
#endif

};

} // End of psi namespace



#endif // MATRIX_DISTRIBUTED_H

