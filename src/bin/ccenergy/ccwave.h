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
    virtual bool same_a_b_orbs() const { return reference_wavefunction_->same_a_b_orbs(); }
    virtual bool same_a_b_dens() const { return reference_wavefunction_->same_a_b_dens(); }

    double compute_energy();

private:
    void init();
};

}}

#endif // CCWAVE_H
