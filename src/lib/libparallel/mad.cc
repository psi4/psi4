#include "parallel.h"
#include <cstring>
#include <pthread.h>

#include <psi4-dec.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace boost;


#ifdef HAVE_MADNESS

MadCommunicator::MadCommunicator(const int &argc, char **argv) :
    Communicator()
{

    // Madness is initialized in the communicator constructor
    // Initialize madness and create a madness world communicator
    madness::initialize(argc, argv);

    madworld_ = SharedMadWorld(new madness::World(MPI::COMM_WORLD));

    me_ = madworld_->rank();
    nproc_ = madworld_->nproc();
    nthread_ = madworld_->taskq.nthread();
    communicator_ = "MADNESS";

    thread_id_.clear();

    // If we have threads in the thread pool, this gets all of the pthread ids from madness, and
    // creates a map of pthread ids to a unique integer value (i.e. 0, 1, 2, ...)
    int i = 0;
    if (Process::environment("MAD_NUM_THREADS") != "1" && Process::environment("POOL_NTHREAD") != "0") {
        pthread_t *id = madworld_->taskq.thread_id();
        for (i=0; i < nthread_; i++) {
            thread_id_.insert( std::pair<pthread_t, int>(id[i], i) );
        }
        free(id);
    }
    // We need to add the master thread because it can perform tasks that are in the task queue
    thread_id_.insert( std::pair<pthread_t, int>(pthread_self(), i));
    // Add one to the total number of threads for the master thread
    nthread_++;

#ifdef _OPENMP
    // We need to disable OpenMP threads since MADNESS has its own.
    omp_set_num_threads(1);
#endif
}

MadCommunicator::MadCommunicator(const MadCommunicator &copy) :
    Communicator(copy),
    madworld_(copy.madworld_), thread_id_(copy.thread_id_)
{ }

MadCommunicator& MadCommunicator::operator =(const MadCommunicator& other)
{
    if (this != &other) {
        me_   = other.me_;
        nproc_ = other.nproc_;
        nthread_ = other.nthread_;
        madworld_ = other.madworld_;
        thread_id_ = other.thread_id_;
    }
    return *this;
}



MadCommunicator::~MadCommunicator()
{   }

void MadCommunicator::sync()
{
    madworld_->gop.fence();
}

void MadCommunicator::raw_send(const void* data, int nbyte, int target)
{
    madworld_->mpi.Isend(data, nbyte, MPI_BYTE, target, 0);
}

void MadCommunicator::raw_recv(void* data, int nbyte, int sender)
{
    madworld_->mpi.Irecv(data, nbyte, MPI_BYTE, sender, 0);
}

void MadCommunicator::raw_bcast(void* data, int nbyte, int broadcaster)
{
    madworld_->gop.broadcast(data, nbyte, broadcaster);
}


#define SUMMEMBER(T) \
void MadCommunicator::raw_sum(T *data, int nelem, T *receive_buffer, int target) \
{ \
    madworld_->gop.sum(data, nelem); \
}

SUMMEMBER(double)
SUMMEMBER(unsigned int)
SUMMEMBER(int)
SUMMEMBER(char)
SUMMEMBER(long)

void MadCommunicator::print(FILE *out) const
{
    if (me() == 0) {
        fprintf(out, "\n    Using MadCommunicator (Number of procs = %d)\n", nproc());
        fprintf(out, "                          (Number of threads in pool = %d)\n\n", nthread_);
    }
}

/*
void function_put() {   }
future function_get() {  }
arg = where and how_much


void MadCommunicator::raw_put(const void* data, int nbyte, int owner, int target)
{
    if (me_ == owner) {
        madworld_->am.send(target,function, args)
    }
}

void MadCommunicator::raw_get(const void* data, int nbyte, int owner, int target)
{
    if (me_ == target) {
        future = madworld_->am.send(owner,function, args)
    }
}
*/


#endif
