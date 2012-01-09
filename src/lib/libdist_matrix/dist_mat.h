#ifndef _psi_src_lib_libdist_array_dist_mat_h
#define _psi_src_lib_libdist_array_dist_mat_h

#include <stdio.h>
#include <cstdio>
#include <string>

#include <libdpd/dpd.h>

#include "tensor/tensor.h"
#include "process_grid.h"

#include <libparallel/parallel.h>


namespace boost {
template<class T> class shared_ptr;
// Forward declarations for boost.python used in the extract_subsets

namespace python {
       class tuple;
}}

namespace psi {

extern FILE *outfile;

//template <typename T>
//class Tile;
//template <typename T>
//class SubTile;

class Distributed_Matrix;
template <class DM>
class SubMatrix;
template <class DM>
class SubValue;


/**
 *  \ingroup Distributed Matrix
 *  \class SubValue
 *  \brief An object that gives the user the
 * ability to access elements using the double
 * bracket notation (i.e. A[i][j]).
 */
template <class DM>
class SubValue
{
private:
    /// This allows the SubValue class to call SubMatrix member functions. (May not need this friend.)
    friend class SubMatrix<DM>; // may not need this friend
    /// This allows the SubValue class to call Distributed_Matrix member functions.
    friend class Distributed_Matrix;

    typedef unsigned int size_type;

    /// The specific row the user requires access to.
    const size_type row_;
    /// The specific column the user requires access to.
    const size_type col_;
    /// A bool that specifies whether or not the user would like access to the entire column.
    const bool all_column_;
    /// A bool that specifies whether or not the user would like access to the entire row.
    const bool all_row_;
    /// The pointer to the distributed matrix object.
    DM *dm_;

    /**
     * The user does not need to use this class so hide default constructor. This
     * should never be used.
     */
    SubValue()
        : dm_(), row_(0), col_(0), all_row_(false), all_column_(false)
    { }

    /**
     * This constructor creates a SubValue object that can be used to access specific elements
     * using the double bracket notation (i.e. A[i][j]).
     *
     * @param dm This is the pointer to the distributed matrix object.
     * @param row The row index that the user requests access to.
     * @param col The column index that the user requests access to.
     */
    SubValue(DM *dm, const int &row, const int &col)
        : dm_(dm), row_(row), col_(col), all_row_(false), all_column_(false)
    { }

    /**
     * This constructor creates a SubValue object that can be used to access an entire row
     * using the double bracket notation (i.e. A[i]["*"]).
     *
     * @param dm This is the pointer to the distributed matrix object.
     * @param row The row index that the user requests access to.
     * @param col This must be the string "*" to signify that the entire row should be accessed.
     */
    SubValue(DM *dm, const int &row, const std::string &col)
        : dm_(dm), row_(row), col_(0), all_row_(true), all_column_(false)
    { }

    /**
     * This constructor creates a SubValue object that can be used to access an entire column
     * using the double bracket notation (i.e. A["*"][j]).
     *
     * @param dm This is the pointer to the distributed matrix object.
     * @param row This must be the string "*" to signify that the entire row should be accessed.
     * @param col The column index that the user requests access to.
     */
    SubValue(DM *dm, const std::string &row, const int &col)
        : dm_(dm), row_(0), col_(col), all_row_(false), all_column_(true)
    { }

public:

    /**
     * This returns a specific value that was specified by A[i][j].
     *
     * If the value is local, get() will return immediately.
     * If the value is not local, get() will wait until the value is
     * fetched from the remote node.
     */
    inline double get()
    {
        return dm_->task(Communicator::world->me(), &DM::get_val, row_, col_).get();
    }

    /**
     * This returns a specific value that was specified by A[i][j].
     *
     * This actually returns a future, which will return immediately.
     * The future can be utlized in another task, which is the ideal use for a future,
     * or it can be accessed by using .get(), which is used to access the future and
     * should not be confused with SubValue.get().
     */
    inline madness::Future<double> get_future() const
    {
        return dm_->task(Communicator::world->me(), &DM::get_val, row_, col_);
    }

    /**
     * This allows a user to set an element, row, or column equal to a value.
     *
     * Examples \n
     * \code
     * A[i][j] = 1.0;       // Set the element A[i][j] equal to 1.0
     * A[i]["*"] = 1.0;     // Set the row i equal to 1.0
     * A["*"][j] = 1.0;     // Set the column j equal to 1.0
     * \endcode
     *
     * @param val The value to set the element, row, or column.
     */
    inline void operator= (const double &val)
    {
        if (all_row_)
            dm_->set_row(row_, val);
        else if (all_column_)
            dm_->set_col(col_, val);
        else
            dm_->task(Communicator::world->me(), &DM::set_val, row_, col_, val);
    }

    /**
     * This allows a user to set an element, row, or column equal to a value
     * that is wrapped by a future.
     *
     * If the user is setting an element to a value,
     * this will submit a task to the taskq and return immediately.  The task may
     * or may not run immediately, so the user should be careful when using this.
     * The user may need to use a sync in order to be sure that the value has be
     * retieved and set.
     *
     * If the user is setting an entire row or column, this may or may not return immediately
     * because setting a row or column is a global operation, which requires the actual value.  The
     * user does not need to use a sync in this case.
     *
     * Examples \n
     * \code
     * A[i][j] = B[k][l].get_future();       // Set the element A[i][j] equal to B[k][l] (MAY NEED A SYNC AFTER)
     * A[i]["*"] = B[k][l].get_future();     // Set the row i equal to B[k][l]
     * A["*"][j] = B[k][l].get_future();     // Set the column j equal to B[k][l]
     * \endcode
     *
     * @param val A future that wraps the value to set the element, row, or column.
     */
    inline void operator= (const madness::Future<double> &val)
    {
        if (all_row_)
            dm_->set_row(row_, val.get());
        else if (all_column_)
            dm_->set_col(col_, val.get());
        else
            dm_->task(Communicator::world->me(), &DM::set_val, row_, col_, val);
    }

    /**
     * This allows a user to set an element, row, or column equal to an element, row, or column value from another distributed matrix.
     *
     * Examples \n
     * \code
     * A[i][j] = B[k][l];       // Set A[i][j] equal to B[k][l]
     * A[i]["*"] = B[k]["*"];   // Set row i of A equal to row k of B
     * A["*"][j] = B["*"][l];   // Set column j of A equal to column l of B
     * \endcode
     *
     * @param sval A SubValue object that points to a distributed matrix that will be copied (i.e. B[i][j]).
     */
    inline void operator= (const SubValue<DM> &sval)
    {
        if (all_row_) {
            if (sval.all_row_)
                dm_->copy_row(this->row_, (*sval.dm_), sval.row_);
            else
                dm_->set_row(this->row_, sval.get_future());
        }
        else if (all_column_) {
            if (sval.all_column_)
                dm_->copy_col(this->col_, (*sval.dm_), sval.col_);
            else
                dm_->set_col(this->col_, sval.get_future());
        }
        else
            this->dm_->task(Communicator::world->me(), &DM::set_val,
                            sval.row_, sval.col_, sval.get_future());
    }

};

/**
 *  \ingroup Distributed Matrix
 *  \class SubMatrix
 *  \brief An object that is returned after the first
 * bracket when using the double bracket notation (i.e. A[i]).
 */
template <class DM>
class SubMatrix
{
private:
    typedef unsigned int size_type;

    /// This allows the SubMatrix class to call Distributed_Matrix member functions.
    friend class Distributed_Matrix;

    /// The specific row the user requires access to.
    const size_type row_;
    /// A bool that specifies whether or not the user would like access to the entire column.
    const bool all_col_;
    /// The pointer to the distributed matrix object.
    DM *dm_;

    /** The user does not need to use this class so hide default constructor. This
     * should never be used.
     */
    SubMatrix()
        : dm_(), row_(0), all_col_(false)
    { }

    /** This constructor creates a SubMatrix object that is returned after the first bracket when
     * using the double bracket notation (i.e. A[i]).
     *
     * @param dm This is the pointer to the distributed matrix object.
     * @param row The row index that the user requests access to.
     *
     */
    SubMatrix(DM *dm, const int &row)
        : dm_(dm), row_(row), all_col_(false)
    { }

    /**
     * This constructor creates a SubMatrix object that is returned after the first bracket when
     * using the double bracket notation (i.e. A[i]).
     *
     * @param dm This is the pointer to the distributed matrix object.
     * @param row This must be the string "*" to signify that the entire row should be accessed.
     *
     */
    SubMatrix(DM *dm, const std::string &row)
        : dm_(dm), row_(0), all_col_(true)
    { }


public:

    /**
     * This overloads the bracket operator to return a SubValue object.  The first bracket
     * returns a SubMatrix object, and the second bracket (defined here) returns the SubValue
     * that can access an element or an entire column.
     *
     * @param col The column that the user needs access to.
    */
    SubValue<DM> operator[] (size_type col)
    {
        if (all_col_) return SubValue<DM>(dm_, "*", col);
        else {
            assert( col < dm_->ncols() );
            return SubValue<DM>(dm_, row_, col);
        }
    }

    /**
     * This overloads the bracket operator to return a SubValue object.  The first bracket
     * returns a SubMatrix object, and the second bracket (defined here) returns the SubValue
     * that can access an entire row.
     *
     * @param col The column that the user needs access to.
    */
    SubValue<DM> operator[] (const std::string &col)
    {
        if (col == "*" && all_col_ == false)
            return SubValue<DM>(dm_, row_, "*");
        else
            throw PSIEXCEPTION("To specify an entire row use Distributed_Matrix[row]['*'].\n");
        if (col == "*" && all_col_ == true)
            throw PSIEXCEPTION("If you want to copy a dist_matrix simply do A = B.\n");
    }

//    void operator= (const double &val)
//    {
//        dm_->set_row(row_, val);
//    }

//    void operator= (const SubMatrix<DM> &smat)
//    {
//        dm_->copy_row(row_, (*smat.dm_), smat.row_);
//    }

};


/**
 *  \ingroup Distributed Matrix
 *  \class Distributed_Matrix
 *  \brief Makes using distributed matrices just a little earlier.
 *
 * This is a disributed matrix class that provides functionality to the user and hides
 * a lot of parallelism.  It is built upon the parallel runtime library from the
 * Multiresolution Adaptive Numerical Environment for Scientific Simulation (MADNESS)
 * software package. The MADNESS parallel runtime was designed around a distributed/threaded
 * task-based parallel scheme, and has been shown to scale to over 100k cores. Specifically,
 * the Distributed_Matrix class heavily utilizes Futures and distributed objects from the
 * MADNESS parallel runtime. Futures are data that will be available in the future, which
 * allows overlap of communication and computation, and distributed objects are objects that
 * are aware of their distribution allowing tasks and data to be easily passed between nodes.
 *
 * \section Distribution
 * The Distributed_Matrix class divides a matrix into tiles, which then get distributed
 * to the processes using a process_grid.\n\n
 * Example:
 * A 5x5 matrix with a tile size of 2x2 is divided as shown below, (Note that there is no padding)
 * \code
 *  ----------------------
 * |  0   1 |  2   3 |  4 |
 * |  5   6 |  7   8 |  9 |
 * |----------------------|
 * | 10  11 | 12  13 | 14 |
 * | 15  16 | 17  18 | 19 |
 * |----------------------|
 * | 20  21 | 22  23 | 24 |
 *  ----------------------
 * \endcode
 * and if we have a 2x2 process_grid of,
 * \code
 *  -------
 * | 0 | 1 |
 * |-------|
 * | 2 | 3 |
 *  -------
 * \endcode
 * then the data distribution of our 5x5 matrix looks like the following
 * \code
 *  ----------------------
 * |  0     |  1     |  0 |
 * |        |        |    |
 * |----------------------|
 * |  2     |  3     |  2 |
 * |        |        |    |
 * |----------------------|
 * |  0     |  1     |  0 |
 *  ----------------------
 * \endcode
 * If we have a 3x2 process_grid as shown below,
 * \code
 *  -------
 * | 0 | 1 |
 * |-------|
 * | 2 | 3 |
 * |-------|
 * | 4 | 5 |
 *  -------
 * \endcode
 * then the data distribution of our 5x5 matrix looks like the following
 * \code
 *  ----------------------
 * |  0     |  1     |  0 |
 * |        |        |    |
 * |----------------------|
 * |  2     |  3     |  2 |
 * |        |        |    |
 * |----------------------|
 * |  4     |  5     |  4 |
 *  ----------------------
 * \endcode
 *
 *  \section Construction
 * The Distributed_Matrix class was designed to allow the user to easily create and manipulate a
 * Distributed_Matrix. The creation of a Distributed_Matrix is done first by creating a process_grid,
 * which is then passed to the constructor of the Distributed_Matrix.\n\n
 * Example:\n
 * (Note. This example uses liboption from PSI4 to read in the processor grid and tile size
 * from a PSI4 input file. However, PSI4 liboption is not required to create either a process
 * grid nor a Distributed_Matrix.
 * \code
 * #include <libdist_matrix/dist_mat.h>
 *
 * // Get the number of dimensions from the PSI4 input file
 * int pgrid_ndims = options["PROCESS_GRID"].size();
 * // Create a vector to contain the process_grid dimensions
 * std::vector<int> pgrid_dimension_sizes;
 * // Get and insert the size of each dimension from the PSI4 input file
 * for (int i=0; i < pgrid_ndims; i++) {
 *   pgrid_dimension_sizes.push_back( options["PROCESS_GRID"][i].to_integer() );
 * }
 *
 * // Create a process_grid object
 * process_grid pgrid(pgrid_dimension_sizes);
 *
 * // Get the tile size of the Distributed_Matrix from the PSI4 input file
 * int tile_sz = options.get_int("TILE_SZ");
 *
 * // Create a 10,000 x 10,000 Distributed_Matrix named "A"
 * Distributed_Matrix A(pgrid, 10000, 10000, tile_sz, "A");
 * \endcode
 * In the above example, the two constructors that are important are the creation of the process_grid
 * and the creation of a 10,000 x 10,000 Distributed_Matrix.
 *
 *  \section Accessing_Elements
 * Once a Distributed_Matrix object has been created, there are a couple of was to access elements.\n \n
 * \par Using the bracket notation
 * The easiest way to access elements is to use the double bracket notation. There is currently support
 * for accessing an individual element, entire row, or entire column. Setting an entire row or column
 * is a global routine that must be called by all processes. However, accessing individual elements
 * can be performed by any process regardless of locality. This is done by using the MADNESS taskq
 * and the main thread will return immediately, so that processing can continue.
 * This is both a blessing and a curse because, while it provides a powerful means to overlap communication with
 * computation, it has the downfall of giving some resposibility to the user to make sure that data
 * has been set before trying to access it. Lunch is never free! We try to indicate where a fence/sync
 * is needed in all of the examples, but be aware of this pitfall when accessing individual elements.
 * A few examples of setting and getting elements, row, or column to a specified value are shown below.\n\n
 * NOTE: The user does not need to manually access the MADNESS taskq. This is handled by the distributed
 * matrix class internally.\n\n
 * Example:
 * \code
 * A[i][j] = 1.0;                       // Set the element A[i][j] equal to 1.0 (MAY NEED A SYNC AFTER)
 * A[i]["*"] = 1.0;                     // Set the row i equal to 1.0 (ALL PROCESSES MUST CALL THIS)
 * A["*"][j] = 1.0;                     // Set the column j equal to 1.0 (ALL PROCESSES MUST CALL THIS)
 *
 * A[i][j] = B[k][l].get();             // Set the element A[i][j] equal to B[k][l] (MAY NEED A SYNC AFTER)
 * A[i]["*"] = B[k][l].get();           // Set the row i equal to B[k][l] (ALL PROCESSES MUST CALL THIS)
 * A["*"][j] = B[k][l].get();           // Set the column j equal to B[k][l] (ALL PROCESSES MUST CALL THIS)
 *
 * A[i][j] = B[k][l].get_future();      // Set the element A[i][j] equal to B[k][l] (MAY NEED A SYNC AFTER)
 * A[i]["*"] = B[k][l].get_future();    // Set the row i equal to B[k][l] (ALL PROCESSES MUST CALL THIS)
 * A["*"][j] = B[k][l].get_future();    // Set the column j equal to B[k][l] (ALL PROCESSES MUST CALL THIS)
 *
 * A[i][j] = B[k][l];                   // Set A[i][j] equal to B[k][l] (MAY NEED A SYNC AFTER)
 * A[i]["*"] = B[k]["*"];               // Set row i of A equal to row k of B (ALL PROCESSES MUST CALL THIS)
 * A["*"][j] = B["*"][l];               // Set column j of A equal to column l of B (ALL PROCESSES MUST CALL THIS)
 *
 * // Set all element in a 5x5 Distributed_Matrix
 * for (int i=0, ij=0; i < 5; j++) {
 *   for (int j=0; j < 5; j++, ij++) {
 *     A[i][j] = ij;
 *   }
 * }
 * // MAY NEED A SYNC here to make sure all of the data has been set
 * \endcode \n
 * The above examples shows two ways to get values from a Distributed_Matrix.  The first example uses
 * .get(), and the second example uses .get_future(), which have the behaviour as described below.\n\n
 * .get()\n
 * If the data is local the value is returned immediately.\n
 * If the data is not local, the main thread will pause until the data from the remote node is retrieved.\n\n
 * .get_future()\n
 * The main thread posts a future and returns immediately, <em>regardless of locality</em>.\n\n
 * The use of .get_future() gives the user the ability to take advantage of futures.  BY FAR, the best
 * way to use a future is to pass them into another task, which is handled by the Distributed_Matrix class internally
 * using the notation shown in the above examples. If the user needs access to the data, the user should use .get().
 *
 * \par Using member functions
 * The use of the double bracket notation has some overhead due to the C++ standard of not allowing brackets
 * to be easily strung together ( i.e. [][] ).  Therefore, to accomidate for two brackets ( i.e A[i][j] ), the first
 * bracket returns a SubMatrix object, which can subsequently take the second bracket, which returns a SubValue object
 * to give access to the actuall element.
 * (NOTE: There is always a trade-off between ease of use and efficiency) Therefore,
 * for those who can not tolerate this overhead, we have also provide member functions to access elements in a
 * Distributed_Matrix.  The same examples provided in the "Using the bracket notation" section are shown below.\n\n
 * NOTE: A fence/sync may be needed due to the reasons discussed in the "Using the bracket notation" section.\n\n
 * Example:
 * \code
 * A.set(i, j, 1.0);                    // Set the element A[i][j] equal to 1.0 (MAY NEED A SYNC AFTER)
 * A.set(i, "*", 1.0);                  // Set the row i equal to 1.0 (ALL PROCESSES MUST CALL THIS)
 * A.set("*", j, 1.0);                  // Set the column j equal to 1.0 (ALL PROCESSES MUST CALL THIS)
 *
 * A.set(i, j, B.get_val(k, l));        // Set the element A[i][j] equal to B[k][l] (MAY NEED A SYNC AFTER)
 * A.set(i, "*", B.get_val(k, l));      // Set the row i equal to B[k][l] (ALL PROCESSES MUST CALL THIS)
 * A.set("*", j, B.get_val(k, l));      // Set the column j equal to B[k][l] (ALL PROCESSES MUST CALL THIS)
 *
 * A.copy_row(i, B, k);                 // Set row i of A equal to row k of B (ALL PROCESSES MUST CALL THIS)
 * A.copy_col(j, B, l);                 // Set column j of A equal to column l of B (ALL PROCESSES MUST CALL THIS)
 *
 * // Set all element in a 5x5 Distributed_Matrix
 * for (int i=0, ij=0; i < 5; j++) {
 *   for (int j=0; j < 5; j++, ij++) {
 *     A.set(i, j, ij);
 *   }
 * }
 * // MAY NEED A SYNC here to make sure all of the data has been set
 * \endcode
 * The above example only has one way to get values from a Distributed_Matrix, which has the following behaviour.\n\n
 * .get_val()\n
 * The main thread posts a future and returns immediately \emph{regardless of locality}.\n\n
 * If the user needs access to the data, the user should use .get_future().get().
 *
 *  \section Programmers_Notes
 * \par Notation
 * There are several indices that need to be utilized in a Distributed_Matrix, so we introduce strict
 * guide lines to follow throughout the Distributed_Matrix class. All examples in this section use
 * the 5x5 Distributed_Matrix that was introduced in the "Distribution" section using a 2x2 process_grid.\n\n
 * Examples:\n\n
 * 1) The index values of the Distributed_Matrix are denoted using "i" for the row index, "j" for the column index,
 * and "ij" for a specific element.\n\n
 * \code
 * ij values are shown inside the matrix
 *
 *    j  0   1    2   3    4
 * i   ----------------------
 * 0  |  0   1 |  2   3 |  4 |
 * 1  |  5   6 |  7   8 |  9 |
 *    |----------------------|
 * 2  | 10  11 | 12  13 | 14 |
 * 3  | 15  16 | 17  18 | 19 |
 *    |----------------------|
 * 4  | 20  21 | 22  23 | 24 |
 *     ----------------------
 * \endcode \n
 * 2) The local index values for the local matrices are denoted using "a" for the local row index,
 * "b" for the local column index, and "ab" for a specific local element.\n\n
 * Examples:
 * \code
 * ab values are shown inside the matrix
 *
 *    b  0   1    0   1    0
 * a   ----------------------
 * 0  |  0   1 |  0   1 |  0 |
 * 1  |  2   3 |  2   3 |  1 |
 *    |----------------------|
 * 0  |  0   1 |  0   1 |  0 |
 * 1  |  2   3 |  2   3 |  1 |
 *    |----------------------|
 * 0  |  0   1 |  0   1 |  0 |
 *     ----------------------
 * \endcode \n
 * 3) The tile index values of the Distributed_Matrix are denoted using "ti" for the tile row index,
 * "tj" for the tile column index, and "tij" for a specific tile.\n\n
 * Examples:
 * \code
 * tij values are shown inside the matrix
 *
 *   tj  0        1        2
 * ti  ----------------------
 * 0  |  0     |  1     |  2 |
 *    |        |        |    |
 *    |----------------------|
 * 1  |  3     |  4     |  5 |
 *    |        |        |    |
 *    |----------------------|
 * 2  |  6     |  7     |  8 |
 *     ----------------------
 * \endcode \n
 * 4) The local tile index values are denoted using "tab" . The notation shown in
 * the example below is (p,tab), which refers to the process that owns the tile and the local tile index,
 * respectively.\n\n
 * Examples:
 * \code
 *   -------------------------
 *  | (0,0)  |  (1,0) | (0,1) |
 *  |        |        |       |
 *  |-------------------------|
 *  |  (2,0) |  (3,0) | (2,1) |
 *  |        |        |       |
 *  |-------------------------|
 *  | (0,2)  | (1,1)  | (0,3) |
 *   -------------------------
 * \endcode \n\n
 * \par Index_Conversion
 * To accomidate the conversion to and from different notations, we include several private conversion member functions.
 * The member functions are shown in the example code below.\n\n
 * Examples:\n
 * \code
 * convert_i_to_ti(i);          // convert index i to tile row ti
 * convert_j_to_tj(j);          // convert index j to tile column tj
 *
 * convert_i_to_a(i);           // convert index i to local matrix index a
 * convert_j_to_b(j);           // convert index j to local matrix index b
 *
 * convert_tij_to_ti(tij);      // convert tile tij to tile row ti
 * convert_tij_to_tj(tij);      // convert tile tij to tile column tj
 *
 * int tab = local(tij);        // convert tile tij to the locally owned tile index tab
 * \endcode
 *
 * \par Using MADNESS
 * The Distributed_Matrix class is actually a madness::WorldObject, which is an object that is aware of its distribution.
 * We say aware because the madness::WorldObject provides the ability to send and receive tasks that are member
 * functions of the object. This functionality is automatically available by simply inheriting from
 * madness::WorldObject<classname>. With this inheritance, the madness::WorldObject functionality can be used as shown
 * in the example code below.\n\n
 * \code
 * // Function:     get_tile_tij(tij) - Returns the tile tij from the owner of tij wrapped in a future
 * // Main Thread:  The main thread posts a future named tile_ij and keeps going.
 * // Thread Pool:  If tile tij is local, the local thread pool runs get_tile_tij(tij).
 * //               If tile tij is not local, the member function get_tile_tij(tij) is
 * //               sent to the owner of tij, and the remote thread runs get_tile_tij(tij).
 * madness::Future<madness::Tensor<double> > tile_ij = this->task(this->owner(tij), &Distributed_Matrix::get_tile_tij, tij);
 *
 * // Function:     copy_tile_tij(tkl, tile_ij) - Copies tile_tij into tile tkl.
 * // Main Thread:  The main thread submits a task (copy_tile_tij) to the taskq.
 * //               The main thread will continue regardless if the all of the
 * //               dependencies are fullfilled. If the main thread hits a fence/sync
 * //               The main thread will become part of the thread pool and run tasks.
 * // Thread Pool:  The following task will not be run by the thread pool until all of
 * //               the dependencies are handled (i.e. tile_ij is available). Once all
 * //               of the dependencies are handled, the task pool of the owner of tile tkl
 * //               will run the member function copy_tile_tij(tkl, tile_ij).
 * this->task(this->owner(tkl), &Distributed_Matrix::copy_tile_tij, tkl, tile_ij);
 * \endcode \n
 * The code above shows how the madness runtime works on a few different levels.  The first describes
 * the function that is being submitted to the task queue, the sencond describes how the main thread behaves, and the third
 * describes how the thread pool handles and runs the task.  To develop code that utilizes the madness runtime it is a good
 * idea to think about how these three things work and interact.  It is for this reason the Distributed_Matrix class has the
 * functionality implemented the way that it is.  Specifically, every function that is provided to the user actually has two
 * pieces.  The first is the function that the user sees, and the second is a private member function that gets submitted to
 * the task queue.  For example, the two parts of code to invert a matrix is shown in the example below.\n\n
 * \code
 *   Distributed_Matrix& Distributed_Matrix::transpose()
 *   {
 *       // This is a global operation so need a sync to make sure there are not any stray tasks
 *       Communicator::world->sync();
 *       // Set up a Distributed_Matrix to temporarily hold the transposed matrix
 *       Distributed_Matrix result(this->pgrid_, this->ncols_, this->nrows_,
 *                                 this->tile_sz_, this->name_);
 *       // Loop through the tiles
 *       for (int ti=0; ti < result.tile_nrows_; ti++) {
 *           for (int tj=0; tj < result.tile_ncols_; tj++) {
 *               // If I own tij
 *               if (me_ == result.owner(ti,tj)) {
 *                   // Get tile tji
 *                   madness::Future<madness::Tensor<double> > tile_ji = this->task(this->owner(tj,ti), &Distributed_Matrix::get_tile, tj, ti);
 *                   // Invert and copy tile tji into tile tij.  Calls the function copy_invert_tile(ti, tj, tile_ji).
 *                   result.task(me_, &Distributed_Matrix::copy_invert_tile, ti, tj, tile_ji);
 *               }
 *           }
 *       }
 *       // This sync makes sure all of the inverting and copying finish
 *       Communicator::world->sync();
 *       // swap *this and result
 *       std::swap(*this, result);
 *       // return *this, result is destroyed
 *       return *this;
 *   }
 *
 *   madness::Void Distributed_Matrix::copy_invert_tile(const int &ti, const int &tj, const madness::Tensor<double> &tile)
 *   {
 *       int loc = local(ti,tj);
 *       int stride = tile.dim(1);
 *       int nrows = data_[loc].dim(0);
 *       double *ptr, *end, *tptr;
 *
 *       // Loop through the rows and copy the inverted tile into this
 *       for (int t=0; t < nrows; t++) {
 *           ptr = &data_[loc](t,0);
 *           end = ptr + data_[loc].dim(1);
 *           tptr = const_cast<double*>(&tile(0,t));
 *           while (ptr < end) {
 *               *ptr++ = *tptr;
 *               tptr+=stride;
 *           }
 *       }
 *   }
 * \endcode \n
 * This example shows an important use of the madness parallel runtime, the use of the taskq. As multi-core processors matures,
 * there will be an increasing number of cores per node, which in terms of MADNESS means that there is an ever increasing number of threads
 * that are a part of the thread pool.  In order to take full advantage of these cores, the bulk of the work need to be running through the
 * taskq, leaving the main thread free to continue work (i.e. setting up more and more tasks). The example above shows two functions.
 * The first function, which is public to the user, is the transpose function. This function is run by the main thread resulting in a
 * number of tasks that are generated.  The first set of tasks that are generated are the get_tile_tij(tj,ti) functions.  This task is the same that
 * was presented in the previous example, and it simply returns the requested tile from whomever owns it.  The second set of tasks that are
 * generated are the cop_invert_tile(ti, tj, tile_ji) functions.  These tasks have one dependency (the future tile_ji).  Once this tile is
 * available (either it is local or it has been recieved from a remote node), the task will be queued up to run in the taskq.  It is the need
 * to utilize the taskq that results in the need for two functions.  The transpose function serves as an API to the user, while the copy_invert_tile is the
 * function needed to utilize the madness taskq.  While this adds more lines to the code, it serves to hide a lot of the parallelism from the
 * user.
 *
 */
class Distributed_Matrix
#ifdef HAVE_MADNESS
        : public madness::WorldObject<Distributed_Matrix>
#endif
{

private:

    friend class SubMatrix<Distributed_Matrix>;

    typedef unsigned int size_type;

    // Process grid stuff
    /// The process grid
    process_grid<2> pgrid_;
    /// The dimension sizes of the process grid
    std::vector<int> pgrid_dims_;

    /// An array that contains the data
    std::vector<madness::Tensor<double> > data_;

    /// Name of the distributed matrix
    std::string name_;

    // Global matrix information
    /// The global number of rows in the distributed matrix
    int nrows_;
    /// The global number of columns in the distributed matrix
    int ncols_;
    /// The total number of elements in the distributed matrix
    int nelements_;

    // Tile information
    /// The tile size
    int tile_sz_;
    /// The number of tile rows
    int tile_nrows_;
    /// The number of tile columns
    int tile_ncols_;
    /// The total number of tiles
    int ntiles_;
    /// The number of local tiles for each process
    int local_ntiles_;
    /// global to local tile map
    std::map<int,int> global_to_local_map_;
    /// local to global tile map
    std::map<int,int> local_to_global_map_;

    // Parallel information
    /// Process ID
    int me_;
    /// Number of MPI Processes
    int nprocs_;
    /// Number of threads
    int nthreads_;
    /// Communicator Type
    std::string comm_;
    /// Madness world object
    SharedMadWorld madworld_;
    /// Mutex for printing
    SharedMutex print_mutex_;
    /// mutex for addition
    SharedMutex add_mutex_;
    /// Mutex for multiplication
    SharedMutex mult_mutex_;
    /// Mutex for setting values
    SharedMutex set_mutex_;

    /**
     * Clears a vector and free's the memory.  This is needed because
     * clear does not actually free the memory of a std::vector.
     *
     * @param vec The std::vector that needs to be cleared.
     */
    template<typename T1>
    void free_vector(std::vector<T1> &vec)
    {
        std::vector<T1> tmp;
        vec.clear();
        vec.swap(tmp);
    }
    /**
     * Clears a std::map and free's the memory.  This is needed because
     * clear does not actually free the memory of a std::map.
     *
     * @param map The std::map that needs to be cleared.
     */
    template<typename T1, typename T2>
    void free_map(std::map<T1,T2> &map)
    {
        std::map<T1,T2> tmp;
        map.clear();
        map.swap(tmp);
    }

    /** Initialize the parallel stuff */
    void parallel_init();

    /**
     * Prints the given tile
     *
     * @param ti The tile row index.
     * @param tj The tile column index.
     * @param tile The tile to be printed.
     */
    madness::Void print_tile(const int &ti, const int &tj, const madness::Tensor<double> &tile) const;
    /**
     * Prints the given tile
     *
     * @param tij The tile index.
     * @param tile The tile to be printed.
     */
    madness::Void print_tile_tij(const int &tij, const madness::Tensor<double> &tile) const;

    /**
     * Zero a specific tile
     *
     * @param ti The tile row index.
     * @param tj The tile column index.
     */
    madness::Void zero_tile(const int &ti, const int &tj);
    /**
     * Zero a specific tile
     *
     * @param tij The tile index.
     */
    madness::Void zero_tile_tij(const int &tij);

    /**
     * Set a tile to identity
     *
     * @param ti The tile row index.
     * @param tj The tile column index.
     */
    madness::Void set_tile_to_identity(const int &ti, const int &tj);
    /**
     * Set a tile to identity
     *
     * @param tij The tile index.
     */
    madness::Void set_tile_to_identity_tij(const int &tij);

    /**
     * Set a tile diagonal to zero
     *
     * @param ti The tile row index.
     * @param tj The tile column index.
     */
    madness::Void set_tile_diagonal_zero(const int &ti, const int &tj);
    /**
     * Set a tile diagonal to zero
     *
     * @param tij The tile index.
     */
    madness::Void set_tile_diagonal_zero_tij(const int &tij);

    /**
     * Return a specified value
     *
     * @param tij The tile index.
     * @param a The local row index.
     * @param b The local column index.
     */
    madness::Future<double> return_tile_val(const int &tij, const int &a,
                                            const int &b);
    /**
     * Set a specific value with val
     *
     * @param tij The tile index.
     * @param a The local row index.
     * @param b The local column index.
     * @param val The value.
     */
    madness::Void set_tile_val(const int &tij, const int &a,
                               const int &b, const double &val);


    /**
     * Fill a tile with a given value
     *
     * @param tij The tile index.
     * @param val The value to fill the tile.
     */
    madness::Void fill_tile(const int &tij, const double &val);

    /**
     * Add tiles from two distributed matrices.
     *
     * @param ti The tile row index.
     * @param tj The tile column index.
     * @param tile The rhs tile that will be added to this(ti,tj) (lhs) tile.
     */
    madness::Void sum_tile(const int &ti, const int &tj, const madness::Tensor<double> &tile);
    /**
     * Add tiles from two distributed matrices.
     *
     * @param tij The tile index.
     * @param tile The rhs tile that will be added to this(ti,tj) (lhs) tile.
     */
    madness::Void sum_tile_tij(const int &tij, const madness::Tensor<double> &tile);

    /**
     * Scale a specific tile by a given value.
     *
     * @param ti The tile row index.
     * @param tj The tile column index.
     * @param val The value to scale the tile.
     */
    madness::Void scale_tile(const int &ti, const int &tj, const double &val);
    /**
     * Scale a specific tile by a given value.
     *
     * @param tij The tile index.
     * @param val The value to scale the tile.
     */
    madness::Void scale_tile_tij(const int &tij, const double &val);

    /**
     * Returns the trace of a specific tile
     *
     * @param ti The tile row index.
     * @param tj The tile column index.
     */
    double trace_tile(const int &ti, const int &tj);
    /**
     * Returns the trace of a specific tile
     *
     * @param tij The tile index.
     */
    double trace_tile_tij(const int &tij);

    /**
     * Invert and copy the given tile. If the tile to be inverted is local
     * to where the tile will be copied to, use this copy invert because a pointer
     * is passed instead of the Tensor.  Need to test this to see if it actually
     * saves any time as compared to passing a future<Tensor<double> >.
     *
     * @param ti The tile row index to copy the inverted tile.
     * @param tj The tile column index to copy the inverted tile.
     * @param tptr The pointer to the local tile to be inverted and copied.
     * @param stride The stride of the local tile to be inverted.
     */
    madness::Void local_copy_invert_tile(const int &ti, const int &tj,
                                         const double *tptr,
                                         const int &stride);
    /**
     * Invert and copy the given tile. This will work for local and non-local tiles.
     *
     * @param ti The tile row index to copy the inverted tile.
     * @param tj The tile column index to copy the inverted tile.
     * @param tile The tile to be inverted.
     */
    madness::Void non_local_copy_invert_tile(const int &ti, const int &tj,
                                             const madness::Tensor<double> &tile);

    /**
     * Returns the dot product of two tiles. This currently is only set up for
     * Distributed_Matrix that are the same size and tile size.
     * TODO: Need to implement a general dot function.
     * TODO: Need to implement one of these with a future<Tensor<double> > to see if the pointer saves any time.
     *
     * @param lptr The pointer to the lhs tile.
     * @param rptr The pointer to the rhs tile.
     * @param size The size of the tiles.
     */
    double vector_dot_tile(const double *lptr, const double *rptr, const int &size);

    /**
     * Set a row in the tile to a given value.
     *
     * @param tij The tile index.
     * @param a The local row index.
     * @param val The value to set the row to.
     */
    madness::Void set_tile_row(const int &tij, const int &a, const double &val);
    /**
     * Set a column in the tile to the a given value.
     *
     * @param tij The tile index.
     * @param b The local column index.
     * @param val The value to set the column to.
     */
    madness::Void set_tile_col(const int &tij, const int &b, const double &val);

    /**
     * Copy the tile into *this->data_[tij].
     *
     * @param tij The tile index where the tile will be copied to.
     * @param tile The tile that will be copied.
     */
    madness::Void copy_tile_tij(const int &tij, const madness::Tensor<double> &tile);
    /**
     * Copy the tile into *this->data_[ti, tj].
     *
     * @param ti The tile row index where the tile will be copied to.
     * @param tj The tile column index where the tile will be copied to.
     * @param tile The tile to be copied.
     */
    madness::Void copy_tile(const int &ti, const int &tj, const madness::Tensor<double> &tile);

    /**
     * Copy a tile row from X into Y (this).  This function is only used
     * if the tile row from X and Y are both local.
     * A pointer to the tile row X is passed instead of
     * a vector containing the row.
     *
     * @param ya the row
     * @param y_tij The tile index that the tile will be copied to.
     * @param *x A pointer to the tile row in X that is being copied.
     */
    madness::Void local_copy_tile_row(const int &ya, const int &y_tij,
                                      const double *x);
    /**
     * Copy a tile row from X into Y (this).  This function works if the tile
     * from X is either local or remote.
     *
     * @param X A vector containing the row that needs to be copied. Use get_tile_row to get this.
     * @param ya The local row index of Y.
     * @param y_tij The tile index of Y.
     */
    madness::Void non_local_copy_tile_row(const std::vector<double> &X,
                                          const int &ya,
                                          const int &y_tij);

    /**
     * Copy a tile column from X into Y (this).  This function is only used
     * if the tile column from X and Y are both local.
     * A pointer to the tile column X is passed instead of
     * a vector containing the column.
     *
     * @param yb the column
     * @param y_tij The tile index that the tile will be copied to.
     * @param *x A pointer to the tile column in X that is being copied.
     */
    madness::Void local_copy_tile_col(const int &b, const int &tij,
                                      const double *copy, const int &stride);
    /**
     * Copy a tile column from X into Y (this).  This function works if the tile
     * from X is either local or remote.
     *
     * @param X A vector containing the column that needs to be copied. Use get_tile_column to get this.
     * @param yb The local column index of Y.
     * @param y_tij The tile index of Y.
     */
    madness::Void non_local_copy_tile_col(const std::vector<double> &X,
                                          const int &y_tij,
                                          const int &yb);

//    /// Copy a tile into a distributed matrix
//    madness::Void copy_from_tile(const int &row, const int &col,
//                                 const int &inc, const madness::Tensor<double> tile);


//    /**
//     * Converts a row/column to a tile row/column
//     *
//     * @param rc The row or column that needs to be converted
//     * @param ya The local row index of Y.
//     * @param y_tij The tile index of Y.
//     */
//    inline int rc_to_tile(const int &rc) { return rc/tile_sz_; }

    /**
     * Return a vector containing a tile row.
     *
     * @param a The local row index.
     * @param tij The tile index.
     */
    std::vector<double> get_tile_row(const int &a, const int &tij);
    /**
     * Return a vector containing a tile column.
     *
     * @param b The local column index.
     * @param tij The tile index.
     */
    std::vector<double> get_tile_col(const int &b, const int &tij);

    /**
     * Do a matrix-matrix muliplication of the two tiles, and return
     * the result in a Tensor
     *
     * @param A The lhs of the mxm.
     * @param B The rhs of the mxm.
     */
    madness::Void mxm(const madness::Tensor<double> &a,
                      const madness::Tensor<double> &b,
                      const int &c,
                      const double &c_scale);

    /**
     * Convert an "ij" matrix index to an "i" matrix row index
     *
     * @param ij The "ij" matrix index to be converted.
     */
    inline int convert_ij_to_i(const int &ij) { return ij/ncols_; }
    /**
     * Convert an "ij" matrix index to a "j" matrix column index
     *
     * @param ij The "ij" matrix index to be converted.
     */
    inline int convert_ij_to_j(const int &ij) { return ij%ncols_; }

    /**
     * Convert an "i" matrix row index to a "ti" tile row index
     *
     * @param i The "i" matrix row index to be converted.
     */
    inline int convert_i_to_ti(const int &i) const { return i/tile_sz_; }
    /**
     * Convert a "j" matrix column index to a "tj" tile column index
     *
     * @param j The "j" matrix column index to be converted.
     */
    inline int convert_j_to_tj(const int &j) const { return j/tile_sz_; }

    /**
     * Convert an "i" matrix row index to an "a" local matrix row index
     *
     * @param i The "i" matrix row index to be converted.
     */
    inline int convert_i_to_a(const int &i) const { return i%tile_sz_; }
    /**
     * Convert a "j" matrix column index to a "b" local matrix column index
     *
     * @param j The "j" matrix column index to be converted.
     */
    inline int convert_j_to_b(const int &j) const { return j%tile_sz_; }

    /**
     * Convert a "tij" tile index to a "ti" tile row index
     *
     * @param tij The "tij" tile index to be converted.
     */
    inline int convert_tij_to_ti(const int &tij) const {return tij/tile_ncols_;}
    /**
     * Convert a "tij" tile index to a "tj" tile column index
     *
     * @param tij The "tij" tile index to be converted.
     */
    inline int convert_tij_to_tj(const int &tij) const {return tij%tile_ncols_;}

    /**
     * Gets and returns the owner of a given tile using the process_grid.
     *
     * @param tij The tile index.
     */
    int pg_map(const int &tij)
    {
        int ti = convert_tij_to_ti(tij);
        int tj = convert_tij_to_tj(tij);

        return pgrid_( (ti)%pgrid_.dim_size(0),
                       (tj)%pgrid_.dim_size(1) ); }

    /**
     * Returns a specific tile.
     *
     * @param ti The tile row index.
     * @param tj The tile column index.
     */
    madness::Future<madness::Tensor<double> > get_tile(const int &ti, const int &tj)
    { return madness::Future<madness::Tensor<double> >( data_[local(ti,tj)] ); }
    /**
     * Returns a specific tile.
     *
     * @param tij The tile index.
     */
    madness::Future<madness::Tensor<double> > get_tile_tij(const int &tij)
    { return get_tile(convert_tij_to_ti(tij), convert_tij_to_tj(tij)); }
    madness::Future<madness::Tensor<double> > get_remote_tile_tij(const int &tij, const int &ti, const int &tj)
    { return get_tile(convert_tij_to_ti(tij), convert_tij_to_tj(tij)); }

    /**
     * Converts the global tile index to a local tile index.
     *
     * @param ti The tile row index.
     * @param tj The tile column index.
     */
    inline int local(const int &ti, const int &tj)
    { return local(ti*tile_ncols_ + tj); }
    /**
     * Converts the global tile index to a local tile index.
     *
     * @param tij The tile index.
     */
    inline int local(const int &tij)
    { return global_to_local_map_[tij]; }

    /**
     * Returns a process grid row.
     *
     * @param row The process grid row index.
     */
    inline std::vector<int> pgrid_row(const int &row)
    { return pgrid_.process_row(row%pgrid_dims_[0]); }
    /**
     * Returns a process grid column.
     *
     * @param row The process grid column index.
     */
    inline std::vector<int> pgrid_col(const int &col)
    { return pgrid_.process_column(col%pgrid_dims_[1]); }

    /**
     * Process 0 prints some matrix information. Useful for debugging.
     *
     */
    void print_matrix_info();

public:

    ~Distributed_Matrix();

    /**
     * Default constructor. It initializes the parallel stuff
     * and zeros everything out.
     */
    Distributed_Matrix();

    /**
     * Constructor, sets up the matrix
     * @param pgrid A process grid.
     * @param nrows Number of rows in the distributed matrix.
     * @param ncols Number of columns in the distributed matrix.
     * @param tile_sz The size of the tiles. Defaults to 64.
     * @param name The name of the distributed matrix.
     */
    Distributed_Matrix(const process_grid<2> &pgrid,
                       const int &nrows, const int &ncols,
                       const int &tile_sz = 64,
                       const std::string &name = "");

    /**
     * Returns the owner of a specific tile
     *
     * @param ti The tile row index.
     * @param tj the tile column index.
     */
    inline int owner(const int &ti, const int &tj) const
    { return pgrid_(ti%pgrid_dims_[0],tj%pgrid_dims_[1]); }
    /**
     * Returns the owner of a specific tile
     *
     * @param tij The tile index.
     */
    inline int owner(const int &tij) const
    { return owner(convert_tij_to_ti(tij), convert_tij_to_tj(tij)); }

    /**
     * Initialize the distributed matrix and allocate the tiles.
     *
     * @param pgrid A process grid.
     * @param nrows Number of rows in the distributed matrix.
     * @param ncols Number of columns in the distributed matrix.
     * @param tile_sz The size of the tiles. Defaults to 64.
     * @param name The name of the distributed matrix.
     */
    void initialize(const process_grid<2> &pgrid,
                     const int &nrows, const int &ncols,
                     const int &tile_sz = 64,
                     const std::string &name = "");

    /**
     * Sets the Distributed_Matrix equal to copy and returns the
     * the new Distributed_Matrix. This allows the user to do the following,\n\n
     * /code
     * B = A;
     * /endcode
     *
     * where A and B are Distributed_Matrix.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESSES
     *
     * @param copy The reference to the Distributed_Matrix that needs to be copied.
     */
    Distributed_Matrix& operator= (const Distributed_Matrix &copy);
    /**
     * Sets the Distributed_Matrix equal to copy and returns the
     * the new Distributed_Matrix. This allows the user to do the following,\n\n
     * /code
     * B = A;
     * /endcode
     *
     * where A is a pointer to a distributed matrix.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESSES
     *
     * @param copy The pointer to the Distributed_Matrix that needs to be copied.
     */
    Distributed_Matrix& operator= (const Distributed_Matrix *copy);
    /**
     * Sets the Distributed_Matrix equal to copy and returns the
     * the new Distributed_Matrix. This allows the user to do the following,\n\n
     * /code
     * Distributed_Matrix A(B);
     * /endcode
     *
     * where A and B are Distributed_Matrix.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESSES
     *
     * @param copy The reference to the Distributed_Matrix that needs to be copied.
     */
    Distributed_Matrix(const Distributed_Matrix &copy);
    /**
     * Sets the Distributed_Matrix equal to copy and returns the
     * the new Distributed_Matrix. This allows the user to do the following,\n\n
     * /code
     * B = A;
     * /endcode
     *
     * where A is a pointer to a distributed matrix.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESSES
     *
     * @param copy The reference to the Distributed_Matrix that needs to be copied.
     */
    Distributed_Matrix(const Distributed_Matrix *copy);

    /**
     * Process 0 prints all of the tiles
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     */
    madness::Void print(const std::string str="") const __attribute__((weak_import)) ;

    /**
     * Process 0 gets and prints a specific tile.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param ti The tile row index.
     * @param tj The tile column index.
     */
    madness::Void print(const int &ti, const int &tj) const;
    /**
     * Process 0 gets and prints a specific tile.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param tij The tile index.
     */
    madness::Void print(const int &tij) const;

    /**
     * Clears the distributed matrix and deallocates the tiles.
     */
    void clear_matrix();

    /**
     * Returns the total number of tiles in the distributed matrix.
     */
    inline int ntiles() const {return ntiles_;}

    /**
     * Zero's the entire distributed matrix.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     */
    madness::Void zero();

    /**
     * Set the distributec matrix to the identity.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     */
    madness::Void identity();

    /**
     * Set the diagonal to zero.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     */
    madness::Void zero_diagonal();

    /**
     * This is the overloaded [] operator.  It returns a SubMatrix,
     * which the user can then use a second [] to access elements
     * and columns.
     *
     * @param row The row index.
     */
    SubMatrix<Distributed_Matrix> operator[] (const size_type &row)
    { return SubMatrix<Distributed_Matrix>(this, row); }

    /**
     * This is the overloaded [] operator.  It returns a SubMatrix,
     * which the user can then use a second [] to access an entire row.
     *
     * @param row A string that must be "*" to indicate the user is
     * accessing an entire row.
     */
    SubMatrix<Distributed_Matrix> operator[] (const std::string &row)
    {
        if (row == "*")
            return SubMatrix<Distributed_Matrix>(this, row);
        else
            throw PSIEXCEPTION("The only allowed strings are *./n");
    }

    /**
     * Returns a future for the requested value.  If the user needs
     * access to the underlying data add .get() to the end of this function.
     *
     * @param ti The tile row index.
     * @param tj The tile column index.
     */
    madness::Future<double> get_val(const int &ti, const int &tj);

    /**
     * Sets a specific value to val. This actually submits a task to the taskq to
     * set the value, so you may need to use a fence/sync afterward to make sure the
     * value is set before doing something else.
     *
     * @param ti The tile row index.
     * @param tj The tile column index.
     * @param val The value.  This can be used in a task with val being a future.  MADNESS
     * handles stripping the value out of the future.
     */

    madness::Void set_val(const int &ti, const int &tj, const double &val);

    /**
     * Set the name of the Distributed_Matrix
     *
     * @param name The name of the Distributed_Matrix
     */
    void set_name(const std::string &name) { name_ = name; }

    /**
     * Returns the tile size of the Distributed_Matrix
     */
    inline int tile_size() const {return tile_sz_;}
    /**
     * Returns the number of rows in the Distributed_Matrix
     */
    inline int nrows() const {return nrows_;}
    /**
     * Returns the number of columns in the Distributed_Matrix
     */
    inline int ncols() const {return ncols_;}
    /**
     * Returns the number of tile columns in the Distributed_Matrix
     */
    inline int tile_ncols() const {return tile_ncols_;}
    /**
     * Returns the number of tile rows in the Distributed_Matrix
     */
    inline int tile_rows() const {return tile_nrows_;}


    /**
     * Fill the Distributed_Matrix with a given value.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param val The value to fill the Distributed_Matrix.
     */
    madness::Void fill(const double &val);
    /**
     * The overloaded equal operator.  This fills the Distributed_Matrix
     * with a given value. Allows the user to do the folowing,\n\n
     * /code
     * A = 2.0;
     * /endcode
     *
     * where A is a Distributed_Matrix.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param val The value to fill the Distributed_Matrix.
     */
    madness::Void operator=(const double &val) { fill(val); }

    /**
     * The overloaded == operator.  This checks to see if the two Distributed_Matrix
     * are the same (i.e. same size and tile size).\n\n
     * /code
     * A == B;
     * /endcode
     *
     * where A and B are a Distributed_Matrix.
     *
     * Returns true if the two Distributed_Matrix are the same.
     *
     * @param rhs The rhs Distributed_Matrix of the == operator.
     */
    bool operator==(const Distributed_Matrix &rhs) const;
    /**
     * The overloaded == operator.  This checks to see if the two Distributed_Matrix
     * are the same (i.e. same size and tile size).\n\n
     * /code
     * A == B;
     * /endcode
     *
     * where A is a Distributed_Matrix and B is a pointer to a Distributed_Matrix.
     *
     * Returns true if the two Distributed_Matrix are the same.
     *
     * @param rhs A pointer to the rhs Distributed_Matrix of the == operator.
     */
    bool operator==(const Distributed_Matrix *rhs) const;
    /**
     * The overloaded != operator.  This checks to see if the two Distributed_Matrix
     * are NOT the same (i.e. different size and/or tile size).  This allows the user
     * to do the following,\n\n
     * /code
     * A != B;
     * /endcode
     *
     * where A is a Distributed_Matrix.
     *
     * Returns true if the two Distributed_Matrix are NOT the same.
     *
     * @param rhs The rhs Distributed_Matrix of the != operator.
     */
    bool operator!=(const Distributed_Matrix &rhs) const;
    /**
     * The overloaded != operator.  This checks to see if the two Distributed_Matrix
     * are NOT the same (i.e. different size and/or tile size). This allows the user
     * to do the following,\n\n
     * /code
     * A != B;
     * /endcode
     *
     * where A is a Distributed_Matrix and B is a pointer to a Distributed_Matrix.
     *
     * Returns true if the two Distributed_Matrix are NOT the same.
     *
     * @param rhs A pointer to the rhs Distributed_Matrix of the != operator.
     */
    bool operator!=(const Distributed_Matrix *rhs) const;

    /**
     * The overloaded += operator. Adds the rhs Distributed_Matrix to the lhs
     * Distributed_Matrix (i.e. *this). This allows the user to do the following,\n\n
     * /code
     * B += A;
     * /endcode
     *
     * where A and B are a Distributed_Matrix.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param rhs The rhs Distributed_Matrix of the += operator.
     */
    madness::Void operator +=(const Distributed_Matrix &rhs);
    /**
     * The overloaded += operator. Adds the rhs Distributed_Matrix to the lhs
     * Distributed_Matrix (i.e. *this). This allows the user to do the following,\n\n
     * /code
     * B += A;
     * /endcode
     *
     * where A is a Distributed_Matrix and B is a pointer to a Distributed_Matrix.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param rhs A pointer to the rhs Distributed_Matrix of the += operator.
     */
    madness::Void operator +=(const Distributed_Matrix *rhs);

    /**
     * The overloaded + operator. Adds the rhs Distributed_Matrix to the lhs
     * Distributed_Matrix (i.e. *this) and returns a Distributed_Matrix that
     * contains the result. This allows the user to do the following,\n\n
     * /code
     * C = B + A;
     * /endcode
     *
     * where A, B and C are a Distributed_Matrix.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param rhs A pointer to the rhs Distributed_Matrix of the += operator.
     */
    Distributed_Matrix operator +(const Distributed_Matrix &rhs) const;
    /**
     * The overloaded + operator. Adds the rhs Distributed_Matrix to the lhs
     * Distributed_Matrix (i.e. *this) and returns a Distributed_Matrix that
     * contains the result. This allows the user to do the following,\n\n
     * /code
     * C = B + A;
     * /endcode
     *
     * where A and B are a Distributed_Matrix and C is a pointer to a Distributed_Matrix.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param rhs A pointer to the rhs Distributed_Matrix of the += operator.
     */
    Distributed_Matrix operator +(const Distributed_Matrix *rhs) const;

    /**
     * Scale every element in the Distributed_Matrix by the given value.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param val The value to scale the Distributed_Matrix.
     */
    madness::Void scale(const double &val);
    /**
     * Scale every element in the Distributed_Matrix by the given value.
     * This allows the user to do the following,\n\n
     * /code
     * A *= 2.0;
     * /endcode
     *
     * where A is a Distributed_Matrix.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param val The value to scale the Distributed_Matrix.
     */

    madness::Void operator*=(const double &val);

    /**
     * Compute the trace of the distributed matrix.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param val The value to scale the Distributed_Matrix.
     */
    double trace() const;

    /**
     * Transpose the distributed matrix. Currently, the transpose of the Distributed_Matrix
     * replaces the object that calls transpose(). So that in the following,
     * /code
     * A.transpose();
     * /endcode
     *
     * A now contains the transpose of the original A.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     */
    Distributed_Matrix& transpose();

    /**
     * Return the dot product of the object that calls vector_dot and rhs.  Currently,
     * The calling object and rhs must be the same (i.e. same size and distribution).
     *
     * TODO: Need a general vector dot function.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     */
    double vector_dot(const Distributed_Matrix &rhs);

    /**
     * Set a row to a specific value.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param i The matrix row index.
     * @param val The value to set all elements in the row.
     */
    madness::Void set_row(const int &i, const double & val);

    /**
     * Set a column to a specific value.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param j The matrix column index.
     * @param val The value to set all elements in the column.
     */
    madness::Void set_col(const int &j, const double &val);

    /**
     * Set a specific element equal to a value.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param i The matrix row index.
     * @param j The matrix column index.
     * @param val The value to set the element to.
     */
    madness::Void set(const int &i, const int &j, const double &val);
    /**
     * Set a specific row equal to a value.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param i The matrix row index.
     * @param col A string that must be "*" to indicate that all elements in a
     *            are being set.
     * @param val The value to set the row to.
     */
    madness::Void set(const int &i, const std::string &col, const double &val);
    /**
     * Set a specific column equal to a value.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param i A string that must be "*" to indicate that all elements in a
     *            are being set.
     * @param j The matrix column index.
     * @param val The value to set the column to.
     */
    madness::Void set(const std::string &row, const int &j, const double &val);

    /**
     * Set a specific element equal to a value.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param i The matrix row index.
     * @param j The matrix column index.
     * @param val A future that will have a value to set the element to.
     */
    madness::Void set(const int &i, const int &j, madness::Future<double> val);
    /**
     * Set a specific row equal to a value.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param i The matrix row index.
     * @param col A string that must be "*" to indicate that all elements in a
     *            are being set.
     * @param val A future that will have a value to set the row to.
     */
    madness::Void set(const int &i, const std::string &col, madness::Future<double> val);
    /**
     * Set a specific column equal to a value.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param i A string that must be "*" to indicate that all elements in a
     *            are being set.
     * @param j The matrix column index.
     * @param val A future that will have a value to set the column to.
     */
    madness::Void set(const std::string &row, const int &j, madness::Future<double> val);


    /**
     * Copy a row (xi) from a Distributed_Matrix (X) into a
     * the row (yi) of the calling Distributed_Matrix (Y).
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param yi The matrix row index of Y (i.e. row to copy to).
     * @param X The Distributed_Matrix to copy from.
     * @param xi The matrix row index of X (i.e. row to copy from).
     */
    madness::Void copy_row(const int &yi, Distributed_Matrix &X, const int &xi);

    /**
     * Copy a column (xj) from a Distributed_Matrix (X) into a
     * the column (yj) of the calling Distributed_Matrix (Y).
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param yj The matrix column index of Y (i.e. column to copy to).
     * @param X The Distributed_Matrix to copy from.
     * @param xj The matrix column index of X (i.e. column to copy from).
     */
    madness::Void copy_col(const int &yj, Distributed_Matrix &X, const int &xj);

    /**
     * Copy elements from a Distributed_Matrix (X) into a the calling
     * Distributed_Matrix (Y). This is equivalent to a dcopy. It will try
     * to use the most efficient algorithm to copy (i.e. if the user is copying
     * a row, this will call copy_row()). However, if an unsusuall pattern is
     * being used to copy, it will default to a copy algorithm that should work
     * in all cases.
     *
     * TODO: This needs to be tested further and improved if possible.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param length The number of elements to copy.
     * @param X The Distributed_Matrix to copy from.
     * @param x0 The row matrix index for the starting point of X.
     * @param xj The column matrix index for the starting point of X.
     * @param x_inc The number of values to increment X.
     * @param y0 The row matrix index for the starting point of Y.
     * @param yj The column matrix index for the starting point of Y.
     * @param y_inc The number of values to increment Y.
     */
    madness::Void copy(const int &length, Distributed_Matrix &X,
                       const int xi, const int xj, const int &x_inc,
                       const int yi, const int yj, const int &y_inc);

//    madness::Void dswap(const int &length, Distributed_Matrix &X,
//                        const int &xrow, const int &xcol, const int &x_inc,
//                        const int &yrow, const int &ycol, const int &y_inc);

//    double ddot(const int &length, Distributed_Matrix &X,
//                const int x[2], const int &x_inc,
//                const int y[2], const int &y_inc);

    /**
     * Do a matrix-matrix multiplication. This allow the user to do the following,\n\n
     * /code
     * C = A * B;
     * /endcode
     *
     * where A is the lhs of the multiplication, B is the rhs of the multiplication,
     * and C is the resuls. A, B, and are Distributed_Matrix.
     *
     * NOTE: MUST BE CALLED BY ALL PROCESS
     *
     * @param b_mat The rhs (Distributed_Matrix) of the matrix-matrix multiplication.
     */
    Distributed_Matrix operator* (Distributed_Matrix &b_mat);


};

} // End of namespace psi

// serialization for double* (only use for local pointers)
#ifdef HAVE_MADNESS
namespace madness {  namespace archive {

    /// Serialize a psi Tile
    template <class Archive>
    struct ArchiveStoreImpl< Archive, double* > {
        static void store(const Archive &ar, const double* t) {
        }
    };

    /// Deserialize a psi Tile ... existing psi Tile is replaced
    template <class Archive>
    struct ArchiveLoadImpl< Archive, double* > {
        static void load(const Archive& ar, double* t) {
        }
    };


} }
#endif // End of HAVE_MADNESS

#endif // dist_mat_h
