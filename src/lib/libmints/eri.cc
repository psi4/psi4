#include "mints.h"

using namespace boost;
using namespace psi;

/////////
// Normal two-electron repulsion integrals
/////////

ERI::ERI(const IntegralFactory *integral, int deriv, double schwarz)
    : TwoElectronInt(integral, deriv, schwarz)
{
    // The +1 is needed for derivatives to work.
    fjt_ = new Taylor_Fjt(basis1()->max_am() +
                          basis2()->max_am() +
                          basis3()->max_am() +
                          basis4()->max_am() +
                          deriv_+1, 1e-15);
}

ERI::~ERI()
{
    delete fjt_;
}

/////////
// F12
/////////

F12::F12(boost::shared_ptr<CorrelationFactor> cf, const IntegralFactory *integral, int deriv, double schwarz)
    : TwoElectronInt(integral, deriv, schwarz)
{
    fjt_ = new F12Fundamental(cf,
                              basis1()->max_am() +
                              basis2()->max_am() +
                              basis3()->max_am() +
                              basis4()->max_am() +
                              deriv_+1);
}

F12::~F12()
{
    delete fjt_;
}

/////////
// F12 squared
/////////

F12Squared::F12Squared(boost::shared_ptr<CorrelationFactor> cf, const IntegralFactory *integral, int deriv, double schwarz)
    : TwoElectronInt(integral, deriv, schwarz)
{
    fjt_ = new F12SquaredFundamental(cf,
                                     basis1()->max_am() +
                                     basis2()->max_am() +
                                     basis3()->max_am() +
                                     basis4()->max_am() +
                                     deriv_+1);
}

F12Squared::~F12Squared()
{
    delete fjt_;
}

/////////
// F12G12
/////////

F12G12::F12G12(boost::shared_ptr<CorrelationFactor> cf, const IntegralFactory *integral, int deriv, double schwarz)
    : TwoElectronInt(integral, deriv, schwarz)
{
    fjt_ = new F12G12Fundamental(cf,
                                 basis1()->max_am() +
                                 basis2()->max_am() +
                                 basis3()->max_am() +
                                 basis4()->max_am() +
                                 deriv_+1);
}

F12G12::~F12G12()
{
    delete fjt_;
}

/////////
// F12DoubleCommutator
/////////

F12DoubleCommutator::F12DoubleCommutator(boost::shared_ptr<CorrelationFactor> cf, const IntegralFactory *integral, int deriv, double schwarz)
    : TwoElectronInt(integral, deriv, schwarz)
{
    fjt_ = new F12DoubleCommutatorFundamental(cf,
                                              basis1()->max_am() +
                                              basis2()->max_am() +
                                              basis3()->max_am() +
                                              basis4()->max_am() +
                                              deriv_+1);
}

F12DoubleCommutator::~F12DoubleCommutator()
{
    delete fjt_;
}

/////////
// ErfERI 
/////////

ErfERI::ErfERI(double omega, const IntegralFactory *integral, int deriv, double schwarz)
    : TwoElectronInt(integral, deriv, schwarz)
{
    // The +1 is needed for derivatives to work.
    fjt_ = new ErfFundamental(omega, 
                          basis1()->max_am() +
                          basis2()->max_am() +
                          basis3()->max_am() +
                          basis4()->max_am() +
                          deriv_+1);
}

ErfERI::~ErfERI()
{
    delete fjt_;
}

void ErfERI::setOmega(double omega)
{
    (static_cast<ErfFundamental*>(fjt_))->setOmega(omega);
}

/////////
// ErfComplementERI 
/////////

ErfComplementERI::ErfComplementERI(double omega, const IntegralFactory *integral, int deriv, double schwarz)
    : TwoElectronInt(integral, deriv, schwarz)
{
    // The +1 is needed for derivatives to work.
    fjt_ = new ErfComplementFundamental(omega, 
                          basis1()->max_am() +
                          basis2()->max_am() +
                          basis3()->max_am() +
                          basis4()->max_am() +
                          deriv_+1);
}

ErfComplementERI::~ErfComplementERI()
{
    delete fjt_;
}

void ErfComplementERI::setOmega(double omega)
{
    (static_cast<ErfComplementFundamental*>(fjt_))->setOmega(omega);
}

