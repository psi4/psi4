#include "detail.h"
#include "libmints_wrapper.h"
#include "libdist_wrapper.h"
#include "elemental_wrapper.h"
    
namespace psi { namespace libmatrix {

// Define the default matrix type. In general, developers should just use this type.
#if !defined(HAVE_MPI) && !defined(HAVE_ELEMENTAL)
    typedef libmints_globals            matrix_globals;
    typedef libmints_matrix_wrapper     matrix;       // No mpi and no elemental
#elif defined(HAVE_MPI) && !defined(HAVE_ELEMENTAL)
    typedef libdist_globals             matrix_globals;
    typedef libdist_matrix_wrapper      matrix;       // mpi and no elemental
#elif defined(HAVE_MPI) && defined(HAVE_ELEMENTAL)
    typedef libelemental_globals        matrix_globals; 
    typedef libelemental_matrix_wrapper matrix;       // mpi and elemental
#endif

    // Provide typedefs of all wrappers to the developer. In some cases using a serial_matrix is preferred.
    typedef libmints_matrix_wrapper     serial_matrix;
    typedef libdist_matrix_wrapper      dist_matrix;
    typedef libelemental_matrix_wrapper elemental_matrix;

    // A templated version of create matrix. Developer must provide template type, e.g.,
    //    serial_matrix mat = create_specific_matrix<serial_matrix>(...);
    template <typename T>
    T create_specific_matrix(const std::string& name, const Dimension& m, const Dimension& n) {
        return detail::create_matrix<T>::create(name, m, n);
    }

    // Create matrix function using the default "matrix" typedef'ed above.
    matrix create_matrix(const std::string& name, const Dimension& m, const Dimension& n) {
        return create_specific_matrix<matrix>(name, m, n);
    }

}}

