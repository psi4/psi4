#ifndef CCWAVE_H
#define CCWAVE_H

// Forward declarations
namespace boost {
template<class T> class shared_ptr;
}

namespace psi {
class Wavefunction;
class Options;
}

namespace psi { namespace ccenergy {

class CCEnergyWavefunction : public Wavefunction
{
public:
    CCEnergyWavefunction(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options);
    virtual ~CCEnergyWavefunction();
    virtual bool restricted() const { return reference_wavefunction_->restricted(); }

    double compute_energy();

private:
    void init();
};

}}

#endif // CCWAVE_H
