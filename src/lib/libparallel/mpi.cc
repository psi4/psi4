#include "parallel.h"
#include <cstring>

using namespace psi;
using namespace boost;

#if HAVE_MPI == 1

MPICommunicator::MPICommunicator(const int &argc, char **argv)
    : Communicator()
{
#if HAVE_MADNESS == 1
    putenv("MAD_NUM_THREADS=1");

    // Initialize madness and create a madness world communicator
    madness::initialize(argc, argv);

    madworld_ = SharedMadWorld(new madness::World(MPI::COMM_WORLD));
    comm_ = MPI::COMM_WORLD;
    me_ = madworld_->rank();
    nproc_ = madworld_->size();
#else
    MPI_Init(const_cast<int*>(&argc), &argv);

    comm_ = MPI_COMM_WORLD;
    MPI_Comm_rank(comm_, &me_);
    MPI_Comm_size(comm_, &nproc_);
#endif
    nthread_ = 1;
    communicator_ = "MPI";
}

MPICommunicator::MPICommunicator(const MPICommunicator &copy) :
    Communicator(), me_(copy.me_), nproc_(copy.nproc_),
    nthread_(copy.nthread_), communicator_(copy.communicator_),
    comm_(copy.comm_)
{ }

MPICommunicator::~MPICommunicator()
{ }

MPICommunicator& MPICommunicator::operator =(const MPICommunicator& other)
{
    if (this != &other) {
        me_   = other.me_;
        nproc_ = other.nproc_;
        nthread_ = other.nthread_;
        comm_ = other.comm_;
    }
    return *this;
}

void MPICommunicator::sync()
{
    MPI_Barrier(comm_);
}


void MPICommunicator::raw_send(const void* data, int nbyte, int target)
{
    MPI_Send(const_cast<void*>(data), nbyte, MPI_BYTE, target, 0, comm_);
}

void MPICommunicator::raw_recv(void* data, int nbyte, int sender)
{
    MPI_Status status;
    MPI_Recv(data, nbyte, MPI_BYTE, sender, 0, comm_, &status);
}

void MPICommunicator::raw_bcast(void* data, int nbyte, int broadcaster)
{
    MPI_Bcast(data, nbyte, MPI_BYTE, broadcaster, comm_);
}

#define SUMMEMBER(T, M) \
void MPICommunicator::raw_sum(T *data, int n, T *receive_buffer, int target) \
{ \
    bool alloc = false; \
    if (receive_buffer == NULL) { \
        alloc = true; \
        receive_buffer = new T[n]; \
    } \
 \
    if (target >= 0) \
        MPI_Reduce(static_cast<void*>(data), static_cast<void*>(receive_buffer), n, M, MPI_SUM, target, comm_); \
    else \
        MPI_Allreduce(static_cast<void*>(data), static_cast<void*>(receive_buffer), n, M, MPI_SUM, comm_); \
 \
    if (alloc) { \
        ::memcpy(static_cast<void*>(data), static_cast<void*>(receive_buffer), sizeof(T)*n); \
        delete[] receive_buffer; \
    } \
}

SUMMEMBER(double, MPI_DOUBLE)
SUMMEMBER(unsigned int, MPI_INT)
SUMMEMBER(int, MPI_INT)
SUMMEMBER(char, MPI_CHAR)
SUMMEMBER(long, MPI_LONG)

void MPICommunicator::print(FILE *out) const
{
    if (me() == 0) {
        fprintf(out, "\n    Using MPICommunicator (Number of processes = %d)\n\n", nproc());
    }
}

void MPICommunicator::finalize() {
#if HAVE_MADNESS == 1
    madness::finalize();
#else
    MPI_Finalize();
#endif
}


#endif
