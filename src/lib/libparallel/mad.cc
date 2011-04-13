#include "parallel.h"
#include <cstring>
#include <pthread.h>

using namespace psi;
using namespace boost;


#if HAVE_MADNESS == 1

MadCommunicator::MadCommunicator(boost::shared_ptr<madness::World> madness_world, const std::string &communicator) :
Communicator(communicator), madworld_(madness_world) {
    me_ = madworld_->rank();
    nproc_ = madworld_->nproc();

    // set the number of threads in the thread pool
    int shift = 0;
    char *cnthread = getenv("MAD_NUM_THREADS");
    // MAD_NUM_THREADS is total no. of application threads whereas
    // POOL_NTHREAD is just the number in the pool (one less)
    if (cnthread) shift = -1;
    if (cnthread == 0) {
        cnthread = getenv("POOL_NTHREAD");
        shift = 1;
    }

    if (cnthread) {
        int result = sscanf(cnthread, "%d", &nthread_);
        nthread_ += shift;
    }
    else {
#ifdef _SC_NPROCESSORS_CONF
        int ncpu = sysconf(_SC_NPROCESSORS_CONF);
#elif defined(HC_NCPU)
        int mib[2]={CTL_HW,HW_NCPU}, ncpu;
        size_t len = sizeof(ncpu);
#else
        int ncpu=1;
#endif
        if (ncpu < 2) ncpu = 2;
        nthread_ = ncpu - 1; // One less than # physical processors
    }

    std::cout << "nthread = " << nthread_ << std::endl;
    thread_count=0;
    mutex_ = boost::shared_ptr<madness::Spinlock>(get_mutex());
    while (thread_count < nthread_) {
        madworld_->taskq.add(*this, &MadCommunicator::set_thread_id);
    }


    std::cout << "waiting at the sync" << std::endl;
    sync();
    std::cout << "Proc " << me_ << " finished getting thread ids.  They are:" << std::endl;

    std::map<pthread_t, int>::iterator it;
    for (it=thread_id_.begin(); it != thread_id_.end(); it++)
        std::cout << "  [" << (*it).first << ", " << (*it).second << "]" << std::endl;

    sync();
}

madness::Void MadCommunicator::set_thread_id() {

    mutex_->lock();

    pthread_t tid;
    std::map<pthread_t, int>::iterator it;

    tid = pthread_self();
    for (it=thread_id_.begin(); it != thread_id_.end(); it++) {
        if (tid == (*it).first) {
           mutex_->unlock();
           sleep(1);
           return madness::None;
        }
    }
    thread_id_.insert( std::pair<pthread_t, int>(tid,thread_count) );
    thread_count++;
    mutex_->unlock();

    return madness::None;
}

MadCommunicator::~MadCommunicator()
    {}

void MadCommunicator::sync() {
    madworld_->gop.fence();
}

void MadCommunicator::raw_send(const void* data, int nbyte, int target){
    madworld_->mpi.Isend(data, nbyte, MPI_BYTE, target, 0);
}

void MadCommunicator::raw_recv(void* data, int nbyte, int sender) {
    madworld_->mpi.Irecv(data, nbyte, MPI_BYTE, sender, 0);
}

void MadCommunicator::raw_bcast(void* data, int nbyte, int broadcaster) {
    madworld_->gop.broadcast(data, nbyte, broadcaster);
}


#define SUMMEMBER(T) \
void MadCommunicator::raw_sum(T *data, int nelem, T *receive_buffer, int target) { \
    madworld_->gop.sum(data, nelem); \
}

SUMMEMBER(double)
SUMMEMBER(unsigned int)
SUMMEMBER(int)
SUMMEMBER(char)
SUMMEMBER(long)

void MadCommunicator::print(FILE *out) const {
    if (me() == 0) {
        fprintf(out, "\n    Using MadCommunicator (Number of procs = %d)\n", nproc());
        fprintf(out, "                          (Number of threads in pool = %d)\n\n", nthread_);
    }
}

#endif
