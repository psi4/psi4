#include "libmints_wrapper.h"
#include "elemental_wrapper.h"

namespace psi { namespace libmatrix {
#if defined(HAVE_ELEMENTAL)
    std::string libelemental_globals::interface_name = "libelemental";
    elem::mpi::Comm libelemental_globals::mpi_comm;
    int libelemental_globals::rank = -1;
#endif
}} // end namespaces

