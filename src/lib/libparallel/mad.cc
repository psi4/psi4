#include "parallel.h"
#include <cstring>

using namespace psi;
using namespace boost;


#if HAVE_MADNESS == 1

    MadCommunicator::MadCommunicator(boost::shared_ptr<madness::World> madness_world) :
    Communicator() {
        madworld_ = madness_world;
        me_ = madworld_->rank();
        nproc_ = madworld_->nproc();

        // set the number of threads in the thread pool
        int shift = 0;
        char *cnthread = getenv("MAD_NUM_THREADS");
        // MAD_NUM_THREADS is total no. of application threads whereas
        // POOL_NTHREAD is just the number in the pool (one less)
        if (cnthread) shift = 1;
        if (cnthread == 0) cnthread = getenv("POOL_NTHREAD");

        if (cnthread) {
            int result = sscanf(cnthread, "%d", &nthread_);
            //if (result != 1)
            //    MADNESS_EXCEPTION("POOL_NTHREAD is not an integer", result);
            nthread_ -= shift;
        }
        else {
#ifdef _SC_NPROCESSORS_CONF
            int ncpu = sysconf(_SC_NPROCESSORS_CONF);
      //if (ncpu <= 0) MADNESS_EXCEPTION("ThreadBase: set_affinity_pattern: sysconf(_SC_NPROCESSORS_CONF)", ncpu);
#elif defined(HC_NCPU)
            int mib[2]={CTL_HW,HW_NCPU}, ncpu;
            size_t len = sizeof(ncpu);
      //if (sysctl(mib, 2, &ncpu, &len, NULL, 0) != 0)
      //  MADNESS_EXCEPTION("ThreadBase: sysctl(CTL_HW,HW_NCPU) failed", 0);
      //std::cout << "NCPU " << ncpu << std::endl;
#else
            int ncpu=1;
#endif
            if (ncpu < 2) ncpu = 2;
            nthread_ = ncpu - 1; // One less than # physical processors
        }
    }

    MadCommunicator::~MadCommunicator() 
        {}

    void MadCommunicator::sync() {
        madworld_->gop.fence();
    }

    void MadCommunicator::barrier() {
        madworld_->mpi.Barrier();
    }

    void MadCommunicator::raw_send(int target, const void *data, int nbyte){
        madworld_->mpi.Isend(data, nbyte, MPI_BYTE, target, 0);
    }

    void MadCommunicator::raw_recv(int sender, void *data, int nbyte) {
        madworld_->mpi.Irecv(data, nbyte, MPI_BYTE, sender, 0);
    }

    void MadCommunicator::sum(double *data, int nelem, double *receive_buffer, int target) {
        madworld_->gop.sum(data, nelem);
    }

    void MadCommunicator::sum(unsigned int *data, int nelem, unsigned int *receive_buffer, int target) {
        madworld_->gop.sum(data, nelem);
    }

    void MadCommunicator::sum(int *data, int nelem, int *receive_buffer, int target) {
        madworld_->gop.sum(data, nelem);
    }

    void MadCommunicator::sum(char *data, int nelem, char *receive_buffer, int target) {
        madworld_->gop.sum(data, nelem);
    }
    
    void MadCommunicator::sum(long *data, int nelem, long *receive_buffer, int target) {
        madworld_->gop.sum(data, nelem);
    }

    void MadCommunicator::raw_bcast(void *data, int nbyte, int broadcaster) {
        madworld_->gop.broadcast(data, nbyte, broadcaster);
    }

    void MadCommunicator::print(FILE *out) const {
        if (me() == 0) {
            fprintf(out, "\n    Using MadCommunicator (Number of procs = %d)\n", nproc());
            fprintf(out, "                          (Number of threads in pool = %d)\n\n", nthread_);
        }
    }

#endif
