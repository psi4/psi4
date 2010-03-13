#include "parallel.h"
#include <cstring>

using namespace psi;
using namespace boost;

LocalCommunicator::LocalCommunicator()
    : Communicator()
{
    nproc_ = 1;
}

LocalCommunicator::LocalCommunicator(const LocalCommunicator &copy)
    : Communicator()
{
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

void LocalCommunicator::raw_send(int target, const void *data, int nbyte)
{

}

void LocalCommunicator::raw_recv(int sender, void *data, int nbyte)
{

}

void LocalCommunicator::sum(double *data, int n, double *receive_buffer, int target)
{
    if (target != -1) {
        ::memcpy(receive_buffer, data, sizeof(double) * n);
    }
}

void LocalCommunicator::sum(unsigned int *data, int n, unsigned int *receive_buffer, int target)
{
    if (target != -1) {
        ::memcpy(receive_buffer, data, sizeof(unsigned int) * n);
    }
}

void LocalCommunicator::sum(int *data, int n, int *receive_buffer, int target)
{
    if (target != -1) {
        ::memcpy(receive_buffer, data, sizeof(int) * n);
    }
}

void LocalCommunicator::sum(char *data, int n, char *receive_buffer, int target)
{
    if (target != -1) {
        ::memcpy(receive_buffer, data, sizeof(char) *n );
    }
}

void LocalCommunicator::sum(long *data, int n, long *receive_buffer, int target)
{
    if (target != -1) {
        ::memcpy(receive_buffer, data, sizeof(long) * n);
    }
}

void LocalCommunicator::sync()
{

}
