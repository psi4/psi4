#ifndef _psi_src_lib_libparallel_openmp_h_
#define _psi_src_lib_libparallel_openmp_h_

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {

void init_openmp();

}
#endif // _psi_src_lib_libparallel_openmp_h_
