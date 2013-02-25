#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <libparallel/openmp.h>

namespace psi {

boost::shared_ptr<worldcomm> WorldComm;

worldcomm* Init_Communicator(int &argc, char **argv) {
    return init_specific_communicator<worldcomm>(argc, argv);
}

void init_openmp() {
#ifdef _OPENMP
    // Disable nested threads in OpenMP
    omp_set_nested(0);
#endif

    // If no OMP thread variable is set, set nthreads to default to 1
    if (Process::environment("OMP_NUM_THREADS") == "")
        Process::environment.set_n_threads(1);
}

}
