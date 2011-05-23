#ifndef LIBMINTS_INTEGRALPARAMETERS_H
#define LIBMINTS_INTEGRALPARAMETERS_H

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class Vector;

class IntegralParameters
{
private:
    unsigned int nparam_;
public:
    IntegralParameters(unsigned int nparam=0)
        : nparam_(nparam) { }
    virtual ~IntegralParameters() {}

    unsigned int nparam() const { return nparam_; }
};

class CorrelationFactor : public IntegralParameters
{
private:
    double *coeff_;
    double *exponent_;

public:
    CorrelationFactor(unsigned int nparam);
    CorrelationFactor(boost::shared_ptr<Vector> coeff,
                      boost::shared_ptr<Vector> exponent);
    virtual ~CorrelationFactor();

    void set_params(boost::shared_ptr<Vector> coeff,
                    boost::shared_ptr<Vector> exponent);
    double *exponent() const { return exponent_; }
    double *coeff() const { return coeff_; }
};

class FittedSlaterCorrelationFactor : public CorrelationFactor
{
private:
    double slater_exponent_;

public:
    FittedSlaterCorrelationFactor(double exponent);
    double exponent(){return slater_exponent_;}
};

}

#endif // INTEGRALPARAMETERS_H
