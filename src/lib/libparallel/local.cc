#include "parallel.h"
#include <psi4-dec.h>
#include <cstring>

using namespace psi;
using namespace boost;

LocalCommunicator::LocalCommunicator(const std::string &communicator)
    : Communicator(communicator)
{
    me_ = 0;
    nproc_ = 1;
}

LocalCommunicator::LocalCommunicator(const LocalCommunicator &copy)
    : Communicator(copy.communicator_), communicator_(copy.communicator_)
{
    me_ = 0;
    nproc_ = 1;
}

LocalCommunicator::~LocalCommunicator()
{

}

LocalCommunicator& LocalCommunicator::operator =(const LocalCommunicator& other)
{
    if (this != &other) {

    }
    return *this;
}

void LocalCommunicator::raw_send(const void* data, int nbyte, int target)
{

}

void LocalCommunicator::raw_recv(void* data, int nbyte, int sender)
{

}

void LocalCommunicator::raw_bcast(void* data, int nbyte, int broadcaster)
{

}

#define SUMMEMBER(T) \
void LocalCommunicator::raw_sum(T *data, int n, T *receive_buffer, int target) \
{ \
    if (receive_buffer != 0) { \
        ::memcpy(receive_buffer, data, sizeof(T) * n); \
    } \
}

SUMMEMBER(double)
SUMMEMBER(unsigned int)
SUMMEMBER(int)
SUMMEMBER(char)
SUMMEMBER(long)

void LocalCommunicator::sync()
{

}

void LocalCommunicator::print(FILE *out) const
{
    fprintf(out, "\n    Using LocalCommunicator (Number of processes = 1)\n\n");
}


