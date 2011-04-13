#include <cstdio>
#include <exception.h>

#include "parallel.h"
#include "serialize.h"

using namespace psi;
using namespace boost;

shared_ptr<Communicator> Communicator::world;

Communicator::Communicator(const std::string &communicator) : me_(0), nproc_(0), communicator_(communicator) {};

Communicator::~Communicator() {};

//void Communicator::sync() {};

//void Communicator::raw_send(const void *data, int nbyte, int target) {};

//void Communicator::raw_recv(void *data, int nbyte, int sender) {};

//void Communicator::raw_bcast(void* data, int nbyte, int broadcaster) {};

//#define SUMMEMBER(type) \
//void Communicator::raw_sum(type *data, int n, type *receive_buffer, int target) \
//{ \
//}

//SUMMEMBER(double)
//SUMMEMBER(unsigned int)
//SUMMEMBER(int)
//SUMMEMBER(char)
//SUMMEMBER(long)

