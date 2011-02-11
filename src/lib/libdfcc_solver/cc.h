#ifndef CC_H
#define CC_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <psi4-dec.h>

using namespace psi;

namespace psi { namespace dfcc {

class CC : public Wavefunction {
private:

protected:

public:
    CC(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~CC();

    virtual double compute_energy()=0;
};

}}

#endif
