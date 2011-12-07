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

typedef SharedMatrix SharedMatrix;

/*! \ingroup MINTS
 *  \class Distributed_Matrix
 *  \brief Makes using distributed matrices just a little earlier.
 *
 * Using a matrix factory makes creating these a breeze.
 * This has NOT been interfaced with the matrix factory yet.
 */
class Distributed_Matrix
#ifdef HAVE_MADNESS
        : public madness::WorldObject<Distributed_Matrix>
#endif
{
protected:

    /// A vector containing the tiles
    std::vector<std::vector<double> > tiles_;
    /// Name of the distributed matrix
    std::string name_;
    /// The tile size
    int tile_sz_;
    /// The total number of tiles
    int ntiles_;
    /// The total number of rows in the distributed matrix
    int nrows_;
    /// The total number of columns in the distributed matrix
    int ncols_;
    /// The total number of elements in the distributed matrix
    int nelements_;
    /// The size of the rows after they have been split
    std::vector<int> row_split_;
    /// The size of the columns after they have been split
    std::vector<int> col_split_;
    /// The number of rows in each tile: Each proc only keeps track of what it owns
    std::vector<int> tile_nrows_;
    /// The number of columns in each tile: Each proc only keeps track of what it owns
    std::vector<int> tile_ncols_;

    /// Process ID
    int me_;
    /// Number of MPI Processes
    int nprocs_;
    /// Number of threads
    int nthreads_;
    /// Communicator Type
    std::string comm_;

    // A global to local tile map
    std::map<int,int> global_local_tile_;
    // A local to global tile map
    std::vector<int> local_global_tile_;


#ifdef HAVE_MADNESS
    /// Madness world object
    SharedMadWorld madworld_;
    /// Mutex for printing
    SharedMutex print_mutex_;
    /// mutex for addition
    SharedMutex add_mutex_;

    /// Sets a specific tile to the identity
    madness::Void set_tile_to_identity(const int &t);
    /// Zeros out a specific tile
    madness::Void zero_tile(const int &t);
    /// Return the value from a specific tile
    madness::Future<double> return_tile_val(const int &t, const int &row, const int &col);
    /// Set the row,col value in the distributred matrix
    madness::Void set_tile_value(const int &t, const int &row,
                                 const int &col, const double &val);
    /// Copies a give tile
    madness::Void copy_tile(const int &t, const std::vector<double> &tile);
    /// Sums the rhs tile and this tile
    madness::Void sum_tile(const int &t, const std::vector<double> &tile);
    /// Fills a given tile with the value
    madness::Void fill_tile(const int &t, const double &val);
    /// Scale the tile by the given value
    madness::Void scale_tile(const int &t, const double &val);
    /// Zero the diagonal of a tile
    madness::Void zero_tile_diagonal(const int &t);
    /// Return the trace of a tile
    double trace_tile(const int &t);
    /// copies the inverse of the given tile
    madness::Void copy_invert_tile(const int &t, const std::vector<double> &tile,
                                   const int &stride);
    /// Returns the dot product of a given tile
    double vector_dot_tile(const int &t, const Distributed_Matrix &rmat);



    /// Do the matrix-matrix multiplication for the give tiles
    madness::Void mxm(const int &t,
                      const std::vector<double> &a,
                      const std::vector<double> &b,
                      const int &a_row,
                      const int &a_col,
                      const int &b_col,
                      const double &c_scale = 0.0);

    std::vector<double> mxm_dns(const int &t,
                                const std::vector<double> &a,
                                const std::vector<double> &b,
                                const int &a_row,
                                const int &a_col,
                                const int &b_col,
                                const double &c_scale);

    template<typename T>
    void free_vector(std::vector<T> &vec)
    {
        std::vector<T> tmp;
        vec.clear();
        vec.swap(tmp);
    }
    template<typename T1, typename T2>
    void free_map(std::map<T1,T2> &map)
    {
        std::map<T1,T2> tmp;
        map.clear();
        map.swap(tmp);
    }

    void clear_matrix();

public:
    /// Default constructor: clears everything out
    Distributed_Matrix();
    /// Constructs a distributed (rows x cols) matrix using the default distribution
    Distributed_Matrix(const int &nrows, const int &ncols, const int &tile_sz = 64,
                       const std::string &name = "");

    /// Initialize the distributed matrix
    void common_init(const int &nrows, const int &ncols, const int &tile_sz = 64,
                     const std::string &name = "");

    /// Return the owner of the block
    int owner(const int &tile) const { return tile%nprocs_; }

    /// Return a tile
    std::vector<double> get_tile(const int &t) { return tiles_[global_local_tile_[t]]; }
    /// Return the number of rows in a tile
    int t_nrow(const int &t) { return tile_nrows_[global_local_tile_[t]]; }
    /// Return the number of cols in a tile
    int t_ncol(const int &t) { return tile_ncols_[global_local_tile_[t]]; }

    /// Print all of the tiles
    madness::Void print_all_tiles() const;
    /// Print a given tile
    madness::Void print_tile(const int &t) const;

    /// Print a given tile (only process 0 should call this)
    madness::Void print_mat(const int &tile, const std::vector<double> &a,
                            const int &m, const int &n) const;

    /// Set the distributed matrix to the identity
    madness::Void identity();
    /// Zero the entire distributed matrix
    madness::Void zero();

    /// Returns the i,j value from the distributed matrix
    madness::Future<double> get_val(const int &row, const int &col) const;
    /// Set the i,j value in the distributed matrix
    madness::Void set_val(const int &row, const int &col, const double &val);

    /// Set the name of the distributed matrix
    void set_name(const std::string &name) { name_ = name; }

    /// Copies the rhs distributed matrix
    Distributed_Matrix& operator =(const Distributed_Matrix &rhs);
    Distributed_Matrix& operator =(const Distributed_Matrix *rhs);
    Distributed_Matrix& operator =(const SharedMatrix mat);

    /// Adds the rhs matrix to this distributed matrix
    madness::Void operator +=(const Distributed_Matrix &rhs);
    madness::Void operator +=(const Distributed_Matrix *rhs);

    /// Adds two matrices and returns the result
    Distributed_Matrix operator +(const Distributed_Matrix &rhs);
    Distributed_Matrix operator +(const Distributed_Matrix *rhs);


    /// Fill the distributed matris with a value
    madness::Void fill(const double &val);

    /// Check to see if the distributed matrices are the same size and tiled the same
    bool operator ==(const Distributed_Matrix &rhs);
    bool operator ==(const Distributed_Matrix *rhs);
    bool operator ==(const SharedMatrix mat) const;

    /// Check to see if the distributed matrices are not the same size and tiled the same
    bool operator !=(const Distributed_Matrix &rhs);
    bool operator !=(const Distributed_Matrix *rhs);

    /// Perform a matrix-matrix multiplication of the distributed matrices and returns the result
    Distributed_Matrix operator *(const Distributed_Matrix &rhs);

    /// Scale the matrix by the given value
    madness::Void scale(const double &val);
    /// Set the diagonal of the distributed matrix to zero
    madness::Void zero_diagonal();

    /// Compute the trace of the distributed matrix
    double trace();

    /// Return the transpose of the distributed matrix
    Distributed_Matrix transpose();

    /// Compute the product of lmat and rmat and add it to this by the c_scale
    madness::Void product(const Distributed_Matrix &lmat,
                          const Distributed_Matrix &rmat,
                          double c_scale = 1.0);

    /// Return the dot product of this and rmat
    double vector_dot(const Distributed_Matrix &rmat);

    /// Transform a distributed matrix with the transformer
    Distributed_Matrix transform(Distributed_Matrix &transformer);

    /// Perform a matrix multiplication using the DNS Algorithm
    madness::Void DNS_MXM(const Distributed_Matrix &lhs, const Distributed_Matrix &rhs);

};
#else
public:
    /// Default constructor: clears everything out
    Distributed_Matrix()
    { throw PSIEXCEPTION("Distributed matrix only works with MADNESS.\n"); }
    /// Constructs a distributed (rows x cols) matrix using the default distribution
    Distributed_Matrix(const int &nrows, const int &ncols, const int &tile_sz = 64,
                       const std::string &name = "")
    { throw PSIEXCEPTION("Distributed matrix only works with MADNESS.\n"); }
    /// Initialize the distributed matrix
    void common_init(const int &nrows, const int &ncols, const int &tile_sz = 64,
                     const std::string &name = "");

    /// Return the owner of the block
    int owner(const int &tile) const { return tile%nprocs_; }

    /// Return a tile
    std::vector<double> get_tile(const int &t) { return tiles_[global_local_tile_[t]]; }
    /// Return the number of rows in a tile
    int t_nrow(const int &t) { return tile_nrows_[global_local_tile_[t]]; }
    /// Return the number of cols in a tile
    int t_ncol(const int &t) { return tile_ncols_[global_local_tile_[t]]; }

    /// Print all of the tiles
    void print_all_tiles() const {}
    /// Print a given tile
    void print_tile(const int &t) const {}

    /// Print a given tile (only process 0 should call this)
    void print_mat(const int &tile, const std::vector<double> &a,
                            const int &m, const int &n) const {}

    /// Set the distributed matrix to the identity
    void identity() {}
    /// Zero the entire distributed matrix
    void zero() {}

    /// Returns the i,j value from the distributed matrix
    double get_val(const int &row, const int &col) const {}
    /// Set the i,j value in the distributed matrix
    void set_val(const int &row, const int &col, const double &val) {}

    /// Set the name of the distributed matrix
    void set_name(const std::string &name) { name_ = name; }

    /// Copies the rhs distributed matrix
    Distributed_Matrix& operator =(const Distributed_Matrix &rhs) {}
    Distributed_Matrix& operator =(const Distributed_Matrix *rhs) {}
    Distributed_Matrix& operator =(const SharedMatrix mat) {}

    /// Adds the rhs matrix to this distributed matrix
    void operator +=(const Distributed_Matrix &rhs) {}
    void operator +=(const Distributed_Matrix *rhs) {}

    /// Adds two matrices and returns the result
    Distributed_Matrix operator +(const Distributed_Matrix &rhs) {}
    Distributed_Matrix operator +(const Distributed_Matrix *rhs) {}


    /// Fill the distributed matris with a value
    void fill(const double &val) {}

    /// Check to see if the distributed matrices are the same size and tiled the same
    bool operator ==(const Distributed_Matrix &rhs) {}
    bool operator ==(const Distributed_Matrix *rhs) {}
    bool operator ==(const SharedMatrix mat) const {}

    /// Check to see if the distributed matrices are not the same size and tiled the same
    bool operator !=(const Distributed_Matrix &rhs) {}
    bool operator !=(const Distributed_Matrix *rhs) {}

    /// Perform a matrix-matrix multiplication of the distributed matrices and returns the result
    Distributed_Matrix operator *(const Distributed_Matrix &rhs) {}

    /// Scale the matrix by the given value
    void scale(const double &val) {}
    /// Set the diagonal of the distributed matrix to zero
    void zero_diagonal() {}

    /// Compute the trace of the distributed matrix
    double trace() {}

    /// Return the transpose of the distributed matrix
    Distributed_Matrix transpose() {}

    /// Compute the product of lmat and rmat and add it to this by the c_scale
    void product(const Distributed_Matrix &lmat,
                          const Distributed_Matrix &rmat,
                          double c_scale = 1.0) {}

    /// Return the dot product of this and rmat
    double vector_dot(const Distributed_Matrix &rmat) {}

    /// Transform a distributed matrix with the transformer
    Distributed_Matrix transform(Distributed_Matrix &transformer) {}

    /// Perform a matrix multiplication using the DNS Algorithm
    void DNS_MXM(const Distributed_Matrix &lhs, const Distributed_Matrix &rhs) {}
};
#endif

} // End of PSI namespace

#ifdef HAVE_MADNESS

namespace madness {  namespace archive {

    /// Serialize a psi Matrix
    template <class Archive>
    struct ArchiveStoreImpl< Archive, psi::Distributed_Matrix> {
        static void store(const Archive &ar, const psi::Distributed_Matrix &t) {
        };
    };

    /// Deserialize a psi Matrix ... existing psi Matrix is replaced
    template <class Archive>
    struct ArchiveLoadImpl< Archive, psi::Distributed_Matrix > {
        static void load(const Archive& ar, psi::Distributed_Matrix &t) {
        };
    };

}}
#endif // End of HAVE_MADNESS


#endif // MATRIX_DISTRIBUTED_H

