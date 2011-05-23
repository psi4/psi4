#include "integralparameters.h"
#include <libmints/vector.h>

using namespace psi;
using namespace boost;

CorrelationFactor::CorrelationFactor(unsigned int nparam)
    : IntegralParameters(nparam)
{
}

CorrelationFactor::CorrelationFactor(boost::shared_ptr<Vector> coeff, boost::shared_ptr<Vector> exponent)
{
    set_params(coeff, exponent);
}

CorrelationFactor::~CorrelationFactor()
{
    delete[] coeff_;
    delete[] exponent_;
}

void CorrelationFactor::set_params(boost::shared_ptr<Vector> coeff, boost::shared_ptr<Vector> exponent)
{
    int nparam = coeff->dim();
    if (nparam) {
        coeff_ = new double[nparam];
        exponent_ = new double[nparam];
        for (int i=0; i<nparam; ++i) {
            coeff_[i] = coeff->get(0, i);
            exponent_[i] = exponent->get(0, i);
        }
    }
}
