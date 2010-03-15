#include "parallel.h"
#include <cstring>

using namespace psi;
using namespace boost;

#if HAVE_MPI == 1

MPICommunicator::MPICommunicator(MPI_Comm comm)
    : Communicator(), comm_(comm)
{
    MPI_Comm_rank(comm_, &me_);
    MPI_Comm_size(comm_, &nproc_);
}

MPICommunicator::MPICommunicator(const MPICommunicator &copy)
    : Communicator(), comm_(copy.comm_)
{
    me_ = copy.me_;
    nproc_ = copy.nproc_;
}

MPICommunicator::~MPICommunicator()
{
}

MPICommunicator& MPICommunicator::operator =(const MPICommunicator& other)
{
    if (this != &other) {
        comm_ = other.comm_;
        me_   = other.me_;
        nproc_ = other.nproc_;
    }
    return *this;
}

void MPICommunicator::sync()
{
    MPI_Barrier(comm_);
}

void MPICommunicator::raw_send(int target, const void* data, int nbyte)
{
    MPI_Send(const_cast<void*>(data), nbyte, MPI_BYTE, target, 0, comm_);
}

void MPICommunicator::raw_recv(int sender, void* data, int nbyte)
{
    MPI_Status status;
    MPI_Recv(data, nbyte, MPI_BYTE, sender, 0, comm_, &status);
}

void MPICommunicator::raw_bcast(void* data, int nbyte, int broadcaster)
{
    MPI_Bcast(data, nbyte, MPI_BYTE, broadcaster, comm_);
}

#define SUMMEMBER(T, M) \
void MPICommunicator::sum(T* data, int n, T* receive_buffer, int target) \
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
        fprintf(out, "\n  Using MPICommunicator (Number of processes = %d\n)\n\n", nproc());
    }
}

#endif
