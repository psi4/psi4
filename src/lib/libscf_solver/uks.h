#ifndef __math_test_uks_h__
#define __math_test_uks_h__

#include <libpsio/psio.hpp>
#include "hf.h"
#include "uhf.h"

#include <psi4-dec.h>

using namespace psi;

namespace psi { namespace scf {

class UKS : public UHF {
protected:

public:
    UKS(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    UKS(Options& options, shared_ptr<PSIO> psio);
    virtual ~UKS();

    double compute_energy();
};

}}

#endif
