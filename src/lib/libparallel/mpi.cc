#include "parallel.h"

using namespace psi;
using namespace boost;

MPICommunicator::MPICommunicator(MPI_Comm comm) : comm_(comm)
{
    MPI_Comm_rank(comm_, &me_);
    MPI_Comm_size(comm_, &nproc_);
}

MPICommunicator::~MPICommunicator()
{
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
    MPI_Bcast(data, nbyte, MPI_BYTE, me_, comm_);
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
    MPI_Reduce(static_cast<void*>(data), static_cast<void*>(receive_buffer), n, M, MPI_SUM, me_, comm_); \
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
