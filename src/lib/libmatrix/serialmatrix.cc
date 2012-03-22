#include <cstring>
#include <libmatrix/matrixwrapper.h>

using namespace psi;
using namespace boost;

namespace psi {

SerialMatrix::SerialMatrix(int rows, int cols, int tile_sz)
{
    serialmat_ = SharedMatrix(new Matrix(rows, cols));
}

SerialMatrix::~SerialMatrix()
{ }

void SerialMatrix::print()
{
    serialmat_->print();
}

void SerialMatrix::fill(double val)
{
    serialmat_->set(val);
}

void SerialMatrix::add(boost::shared_ptr<MatrixWrapper> rhs)
{
//    boost::shared_ptr<SerialMatrix> tmp = boost::shared_polymorphic_downcast<SerialMatrix>(rhs);
    serialmat_->add(rhs->serialmat_);
}






} // End of namespace PSI
