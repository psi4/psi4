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

