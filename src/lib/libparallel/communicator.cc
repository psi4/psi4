#include <cstdio>
#include <exception.h>

#include "parallel.h"
#include "serialize.h"

using namespace psi;
using namespace boost;

shared_ptr<Communicator> Communicator::world;

Communicator::Communicator() : me_(0), nproc_(0)
{
}

Communicator::~Communicator()
{
}

#define SENDMEMBER(type) \
void Communicator::send(int target, const type *data, int ndata) \
{ \
    raw_send(target, static_cast<const void*>(data), ndata * sizeof(type)); \
}

SENDMEMBER(double)
SENDMEMBER(unsigned int)
SENDMEMBER(int)
SENDMEMBER(char)
SENDMEMBER(long)

void Communicator::send(int target, double data)
{
    raw_send(target, const_cast<const double*>(&data), 1);
}

void Communicator::send(int target, int data)
{
    raw_send(target, const_cast<const int*>(&data), 1);
}

#define RECVMEMBER(type) \
void Communicator::recv(int sender, type *data, int ndata) \
{ \
    raw_recv(sender, static_cast<void*>(data), ndata * sizeof(type)); \
}

RECVMEMBER(double)
RECVMEMBER(unsigned int)
RECVMEMBER(int)
RECVMEMBER(char)
RECVMEMBER(long)

void Communicator::recv(int sender, double& data)
{
    raw_recv(sender, &data, 1);
}

void Communicator::recv(int sender, int& data)
{
    raw_recv(sender, &data, 1);
}

#define BCASTMEMBER(type) \
void Communicator::bcast(type *data, int ndata, int broadcaster) \
{ \
    raw_bcast(static_cast<void *>(data), ndata * sizeof(type), broadcaster); \
}

void Communicator::raw_bcast(void*, int, int)
{
    
}

BCASTMEMBER(double)
BCASTMEMBER(unsigned int)
BCASTMEMBER(int)
BCASTMEMBER(char)
BCASTMEMBER(long)

void Communicator::bcast(Serializable *data, int broadcaster)
{
    data->bcast(this, broadcaster);
}

void Communicator::bcast(double& data, int broadcaster)
{
    raw_bcast(&data, sizeof(double), broadcaster);
}

void Communicator::bcast(int& data, int broadcaster)
{
    raw_bcast(&data, sizeof(int), broadcaster);
}

#define SUMMEMBER(type) \
void Communicator::sum(type *data, int n, type *receive_buffer, int) \
{ \
    if (receive_buffer) \
        ::memcpy(static_cast<void *>(data), static_cast<void *>(receive_buffer), n * sizeof(type)); \
}

SUMMEMBER(double)
SUMMEMBER(unsigned int)
SUMMEMBER(int)
SUMMEMBER(char)
SUMMEMBER(long)

void Communicator::sum(double &data)
{
    sum(&data, 1);
}

void Communicator::sum(int &data)
{
    sum(&data, 1);
}
