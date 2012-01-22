/*
 *  dfmp2_energy.h
 *
 *
 *
 */

#ifndef DFMP2ENERGY_H
#define DFMP2ENERGY_H

#include <libmints/wavefunction.h>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class Options;
class PSIO;
class Chkpt;

namespace dfmp2 {

class DFMP2Energy : public Wavefunction {
public:
    DFMP2Energy(Options & options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    virtual ~DFMP2Energy() {}

    double compute_E();
};

}}

#endif
