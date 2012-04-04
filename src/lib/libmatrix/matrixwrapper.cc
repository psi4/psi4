#include <cstring>
#include <libmatrix/matrixwrapper.h>

using namespace psi;
using namespace boost;

namespace psi {

MatrixWrapper::MatrixWrapper()
{

    me_ = Communicator::world->me();
    nprocs_ = Communicator::world->nproc();
    nthreads_ = Communicator::world->nthread();
    comm_ = Communicator::world->communicator();

}

MatrixWrapper::~MatrixWrapper()
{ }

} // End namespace psi

