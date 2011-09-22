#ifndef CIWAVE_H
#define CIWAVE_H

// Forward declarations
namespace boost {
template<class T> class shared_ptr;
}

namespace psi {
class Wavefunction;
class Options;
}

namespace psi { namespace detci {

class CIWavefunction : public Wavefunction
{
public:
    CIWavefunction(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options);
    virtual ~CIWavefunction();

    double compute_energy();

private:
    void init();
};

}}

#endif // CIWAVE_H

