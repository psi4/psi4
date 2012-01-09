#ifndef _psi_src_lib_libmints_process_grid_h_
#define _psi_src_lib_libmints_process_grid_h_

#if HAVE_MADNESS

#include <cstdio>
#include <string>
#include <cstring>
#include <vector>
#include <cstdarg>

#include <exception.h>
#include "tensor/tensor.h"

#include <libparallel/parallel.h>
//#include <boost/multi_array.hpp>

namespace boost {
template<class T> class shared_ptr;
// Forward declarations for boost.python used in the extract_subsets

namespace python {
       class tuple;
}}

namespace psi {

class PSIO;
//class process_grid;

template <size_t NumDims=2>
class process_grid
{
private:
    int ndims_;
    std::vector<int> dims_;
    madness::Tensor<int> pgrid_;

    /// Process ID
    int me_;
    /// Number of MPI Processes
    int nprocs_;
    /// Number of threads
    int nthreads_;
    /// Communicator Type
    std::string comm_;

    void parallel_init()
    {
        me_ = Communicator::world->me();
        nprocs_ = Communicator::world->nproc();
        nthreads_ = Communicator::world->nthread();
        comm_ = Communicator::world->communicator();
    }

    void fill()
    {
        int v=0;

        int *ptr = pgrid_.ptr();
        int *end = ptr + pgrid_.size();
        while (ptr < end) {
            *ptr++ = v%nprocs_;
            v++;
        }
    }

public:
    /// The default constructor sets only initializes the parallel stuff
    process_grid()
    {
        parallel_init();
    }


    /**
          * This constructor allows the user to set the dim sizes and initialize the process grid
          *
          * @param ndim A vector containing the dimension sizes for the process grid
          */
    process_grid(std::vector<int> &dims)
    {
        // Initialize the parallel stuff
        parallel_init();
        // initialize the process grid and multi_array_ref
        initialize(dims);
    }

    /**
      * This allows the user to initialize the dim sizes of the process grid
      *
      * @param ndim A vector containing the dimension sizes for the process grid
      */
    void initialize(std::vector<int> &dims) {

        pgrid_.clear();


        // This checks to make sure the input vector contains the
        // correct number of dimension sizes.  If not, it will add
        // dimension sizes of one for the unprovided dim sizes.
        if (dims.size() < NumDims) {
            unsigned int size = NumDims-dims.size();
            for (unsigned int i=0; i < size; i++)  {
                dims.push_back(1);
            }
        }
        dims_ = dims;

        // Compute the number of elements (size) in the process group
        unsigned int size = 1;
        for (unsigned int i = 0; i < NumDims; ++i)  {
            size *= dims[i];
        }

        // Check to make sure we have enough processors for the process group
        if (size > nprocs_) {
            throw PSIEXCEPTION("The number of processors must be equal to or "
                               "greater than the product of the process group "
                               "dimension sizes.\n");
        }
        else {
            pgrid_ = madness::Tensor<int>(dims[0], dims[1]);
            fill();
        }
    }

    /// Returns the number of elements in the process group
    unsigned int size() const { return pgrid_.size(); }
    /// Returns the size of a given dimension
    unsigned int dim_size(unsigned int i) const { return pgrid_.dim(i); }
    /// Returns the number of dimensions in the process group
    unsigned int ndims() const { return NumDims; }

    std::vector<int> dimension_grid() const { return dims_; }

    inline int operator ()(const int &i, const int &j) const
    {
        if (NumDims == 2) return pgrid_(i,j);
        else throw PSIEXCEPTION("Number of indices does not match NumDims.\n");
    }

    inline int operator ()(const int &i, const int &j,
                           const int &k) const
    {
        if (NumDims == 3) return pgrid_(i,j,k);
        else throw PSIEXCEPTION("Number of indices does not match NumDims.\n");
    }


    bool operator==(const process_grid &comp) const
    {
        if (this->ndims_ != comp.ndims_) return false;
        else if (this->dims_ != comp.dims_) return false;
        else return true;
    }
    bool operator==(const process_grid *comp) const
    {
        if (this->ndims_ != comp->ndims_) return false;
        else if (this->dims_ != comp->dims_) return false;
        else return true;
    }

    bool operator !=(const process_grid &comp) const
    {
        if (*this == comp) return false;
        else return true;
    }

    bool operator !=(const process_grid *comp) const
    {
        if (*this == comp) return false;
        else return true;
    }

    std::vector<int> process_row(const int &i)
    {
        std::vector<int> temp(pgrid_.dim(1), 0.0);

        for (int j=0; j < pgrid_.dim(1); j++) {
            temp[j] = pgrid_(i,j);
        }
        return temp;
    }

    std::vector<int> process_column(const int &j)
    {
        std::vector<int> temp(pgrid_.dim(0), 0.0);

        for (int i=0; i < pgrid_.dim(0); i++) {
            temp[i] = pgrid_(i,j);
        }
        return temp;
    }

};



//template <typename T, unsigned int NumDims>
//class process_grid_ref;

//extern FILE *outfile;

//template <typename T, unsigned int NumDims>
//class process_grid_ref : public boost::multi_array_ref<T,NumDims>
//{
//private:

//    typedef boost::detail::multi_array::extent_gen<NumDims> extents;
//    typedef boost::multi_array_ref<T,NumDims> array_interface;

//    /// The process grid
//    std::vector<T> data_;


//    /// Process ID
//    int me_;
//    /// Number of MPI Processes
//    int nprocs_;
//    /// Number of threads
//    int nthreads_;
//    /// Communicator Type
//    std::string comm_;

//    void parallel_init()
//    {
//        me_ = Communicator::world->me();
//        nprocs_ = Communicator::world->nproc();
//        nthreads_ = Communicator::world->nthread();
//        comm_ = Communicator::world->communicator();
//    }

//    void fill()
//    {
//        int v=0;
//        std::vector<int>::iterator it;
//        for (it = data_.begin(); it != data_.end(); it++, v++) {
//            *it = v%nprocs_;
//        }
//    }


//public:

//    /// The default constructor sets only initializes the parallel stuff
//    process_grid_ref() : array_interface(NULL, extents())
//    {
//        parallel_init();
//    }


//    /**
//      * This constructor allows the user to set the dim sizes and initialize the process grid
//      *
//      * @param ndim A vector containing the dimension sizes for the process grid
//      */
//    process_grid_ref(std::vector<T> &dims) : array_interface(NULL, extents())
//    {
//        // Initialize the parallel stuff
//        parallel_init();
//        // initialize the process grid and multi_array_ref
//        initialize(dims);
//    }

//    process_grid_ref(const process_grid_ref<T, NumDims> &pgrid) : array_interface(NULL, extents())
//    {
//        std::vector<int> dims;
//        for (unsigned int i=0; i < NumDims; i++) {
//            dims.push_back(pgrid.dim_size(i));
//        }
//        this->initialize(dims);
//    }


//    process_grid_ref<T, NumDims>& operator=(const process_grid_ref<T, NumDims> &pgrid)
//    {
//        std::vector<int> dims;
//        for (unsigned int i=0; i < NumDims; i++) {
//            dims.push_back(pgrid.dim_size(i));
//        }
//        this->initialize(dims);

//        return *this;
//    }

//    /**
//      * This allows the user to initialize the dim sizes of the process grid
//      *
//      * @param ndim A vector containing the dimension sizes for the process grid
//      */
//    void initialize(std::vector<int> &dims) {

//        data_.clear();

//        // This checks to make sure the input vector contains the
//        // correct number of dimension sizes.  If not, it will add
//        // dimension sizes of one for the unprovided dim sizes.
//        if (dims.size() < NumDims) {
//            unsigned int size = NumDims-dims.size();
//            for (unsigned int i=0; i < size; i++)  {
//                dims.push_back(1);
//            }
//        }

//        // Compute the number of elements (size) in the process group
//        unsigned int size = 1;
//        for (unsigned int i = 0; i < NumDims; ++i)  {
//            size *= dims[i];
//        }

//        // Check to make sure we have enough processors for the process group
//        if (size > nprocs_) {
//            throw PSIEXCEPTION("The number of processors must be equal to or "
//                               "greater than the product of the process group "
//                               "dimension sizes.\n");
//        }
//        else {
//            // clear and resize the array of process id's
//            data_.clear();
//            data_.resize(size, 0);

//            // This points the multi_array_ref to our array of process id's
//            set_base_ptr(&data_[0]);
//            // Set the number of elements in the multi_array_ref
//            this->num_elements_ = size;
//            // copy the vector of dim sizes to a boost::array
//            boost::array<unsigned int,NumDims> shape;
//            std::copy(dims.begin(), dims.end(), shape.begin());
//            // This sets the dimension sizes of the multi_array_ref
//            reshape(shape);
//            // Fill in the process id's. Currently round-robin scheme.
//            fill();
//        }
//    }

//    /// Returns the number of elements in the process group
//    unsigned int size() const { return data_.size(); }
//    /// Returns the size of a given dimension
//    unsigned int dim_size(unsigned int i) const { return this->shape()[i]; }
//    /// Returns the number of dimensions in the process group
//    unsigned int ndims() const { return NumDims; }

//};


//class process_grid
//{
//private:
//    int ndims_;
//    std::vector<int> dims_;

//    /// These are the defined process grids
//    process_grid_ref<int, 1> P1;
//    process_grid_ref<int, 2> P2;
//    process_grid_ref<int, 3> P3;
//    process_grid_ref<int, 4> P4;
//    process_grid_ref<int, 5> P5;
//    process_grid_ref<int, 6> P6;

//    void initialize(const int &ndims, std::vector<int> &grid)
//    {
//        ndims_ = ndims;
//        dims_ = grid;
//        if (Communicator::world->me() == 0)
//            std::cout << "ndims = " << ndims_ << std::endl;
//        switch (ndims_) {
//            case 0: throw PSIEXCEPTION("init(0): Process grid number of dimensions must be between 0-6");
//                break;
//            case 1: P1.initialize(grid);
//                break;
//            case 2: P2.initialize(grid);
//                break;
//            case 3: P3.initialize(grid);
//                break;
//            case 4: P4.initialize(grid);
//                break;
//            case 5: P5.initialize(grid);
//                break;
//            case 6: P6.initialize(grid);
//                break;
//            default: throw PSIEXCEPTION("init(default): Process grid number of dimensions must be between 0-6");
//                break;
//        }
//    }
//public:

//    process_grid() { }
//    process_grid(const process_grid &pgrid)
//    {
//        int ndims = pgrid.ndims();
//        std::vector<int> grid = pgrid.dimension_grid();
//        this->initialize(ndims, grid);
//    }

//    process_grid(const int &ndims, std::vector<int> &grid)
//    {
//        ndims_ = ndims;
//        dims_ = grid;
//        switch (ndims_) {
//            case 0: throw PSIEXCEPTION("process_grid(0): Process grid number of dimensions must be between 0-6");
//                break;
//            case 1: P1.initialize(grid);
//                break;
//            case 2: P2.initialize(grid);
//                break;
//            case 3: P3.initialize(grid);
//                break;
//            case 4: P4.initialize(grid);
//                break;
//            case 5: P5.initialize(grid);
//                break;
//            case 6: P6.initialize(grid);
//                break;
//            default: throw PSIEXCEPTION("process_grid(default): Process grid number of dimensions must be between 0-6");
//                break;
//        }
//    }

//    std::vector<int> dimension_grid() const { return dims_; }
//    int ndims() const { return ndims_; }

//    template <typename T>
//    T& grid()
//    {
//        switch (ndims_) {
//            case 0: throw PSIEXCEPTION("grid(0): Process grid number of dimensions must be between 0-6");
//                break;
//            case 1: return P1;
//                break;
//            case 2: return P2;
//                break;
//            case 3: return P3;
//                break;
//            case 4: return P4;
//                break;
//            case 5: return P5;
//                break;
//            case 6: return P6;
//                break;
//            default: throw PSIEXCEPTION("grid(default): Process grid number of dimensions must be between 0-6");
//                break;
//        }
//    }


//    inline int operator()(const int &i) const
//    {
//        return P1[i];
//    }

//    inline int operator()(const int &i, const int &j) const
//    {
//        return P2[i][j];
//    }

//    inline int operator()(const int &i, const int &j,
//                          const int &k) const
//    {
//        return P3[i][j][k];
//    }

//    inline int operator ()(const int &i, const int &j,
//                           const int &k, const int &l) const
//    {
//        return P4[i][j][k][l];
//    }

//    inline int operator ()(const int &i, const int &j,
//                           const int &k, const int &l,
//                           const int &m) const
//    {
//        return P5[i][j][k][l][m];
//    }

//    inline int operator ()(const int &i, const int &j,
//                           const int &k, const int &l,
//                           const int &m, const int &n) const
//    {
//        return P6[i][j][k][l][m][n];
//    }

//    int dim_size(const int &i) const
//    {
//        switch (ndims_) {
//            case 0: throw PSIEXCEPTION("dim_size(0): Process grid number of dimensions must be between 0-6");
//            case 1: return P1.dim_size(i);
//            case 2: return P2.dim_size(i);
//            case 3: return P3.dim_size(i);
//            case 4: return P4.dim_size(i);
//            case 5: return P5.dim_size(i);
//            case 6: return P6.dim_size(i);
//            default: throw PSIEXCEPTION("dim_size(default): Process grid number of dimensions must be between 0-6");
//        }
//    }

//    bool operator==(const process_grid &comp) const
//    {
//        if (this->ndims_ != comp.ndims_) return false;
//        else if (this->dims_ != comp.dims_) return false;
//        else return true;
//    }
//    bool operator==(const process_grid *comp) const
//    {
//        if (this->ndims_ != comp->ndims_) return false;
//        else if (this->dims_ != comp->dims_) return false;
//        else return true;
//    }

//    bool operator !=(const process_grid &comp) const
//    {
//        if (*this == comp) return false;
//        else return true;
//    }

//    bool operator !=(const process_grid *comp) const
//    {
//        if (*this == comp) return false;
//        else return true;
//    }

//    std::vector<int> process_row(const int &i)
//    {
//        std::vector<int> temp(P2.dim_size(1), 0.0);

//        for (int j=0; j < P2.dim_size(1); j++) {
//            temp[j] = P2[i][j];
//        }
//        return temp;
//    }

//    std::vector<int> process_column(const int &j)
//    {
//        std::vector<int> temp(P2.dim_size(0), 0.0);

//        for (int i=0; i < P2.dim_size(0); i++) {
//            temp[i] = P2[i][j];
//        }
//        return temp;
//    }

//};



} // End of psi4 namespace

#endif
#endif // PROCESS_GRID_H
