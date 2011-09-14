#include "parallel.h"
#include <psi4-dec.h>
#include <cstring>

using namespace psi;
using namespace boost;

LocalCommunicator::LocalCommunicator(const int &argc, char **argv)
    : Communicator()
{

#ifdef HAVE_MADNESS
    putenv("MAD_NUM_THREADS=1");

    // Initialize madness and create a madness world communicator
    madness::initialize(argc, argv);

    madworld_ = SharedMadWorld(new madness::World(MPI::COMM_WORLD));
#endif

    me_ = 0;
    nproc_ = 1;
    nthread_ = 1;
    communicator_ = "LOCAL";
}

LocalCommunicator::LocalCommunicator(const LocalCommunicator &copy) :
    Communicator(copy)
{ }

LocalCommunicator::~LocalCommunicator()
{

}

LocalCommunicator& LocalCommunicator::operator =(const LocalCommunicator& other)
{
    if (this != &other) {
        me_   = other.me_;
        nproc_ = other.nproc_;
        nthread_ = other.nthread_;
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

void LocalCommunicator::finalize() {
#ifdef HAVE_MADNESS
    madness::finalize();
#endif
}


