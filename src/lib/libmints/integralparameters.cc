#include <libmints/vector.h>
#include "integralparameters.h"

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

FittedSlaterCorrelationFactor::FittedSlaterCorrelationFactor(double exponent)
    : CorrelationFactor(6)
{
    // Perform the fit.
    SharedVector exps(new Vector(6));
    SharedVector coeffs(new Vector(6));

    slater_exponent_ = exponent;

    // The fitting coefficients
    coeffs->set(0, 0, -0.3144);
    coeffs->set(0, 1, -0.3037);
    coeffs->set(0, 2, -0.1681);
    coeffs->set(0, 3, -0.09811);
    coeffs->set(0, 4, -0.06024);
    coeffs->set(0, 5, -0.03726);

    // and the exponents
    exps->set(0, 0, 0.2209);
    exps->set(0, 1, 1.004);
    exps->set(0, 2, 3.622);
    exps->set(0, 3, 12.16);
    exps->set(0, 4, 45.87);
    exps->set(0, 5, 254.4);

    // They just need to be scaled
    double expsq = exponent * exponent;
    exps->scale(expsq);
    set_params(coeffs, exps);
}
