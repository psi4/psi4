#include <cstdio>
#include <exception.h>

#include "parallel.h"
#include "serialize.h"

using namespace psi;
using namespace boost;

shared_ptr<Communicator> Communicator::world;

Communicator::Communicator() :
    me_(0), nproc_(0), nthread_(1), communicator_("Derived Communicator")
{ }

Communicator::~Communicator()
{ }

void Communicator::bcast(std::string& data, int broadcaster)
{
    size_t length = 0;

    if (broadcaster == me_)
        length = data.size();
    bcast(&length, 1, broadcaster);

    if (broadcaster != me_)
        data.resize(length, 0xff);

    char* tdata = const_cast<char*>(data.c_str());
    bcast(tdata, length, broadcaster);
    tdata[length] = '\0';
}
