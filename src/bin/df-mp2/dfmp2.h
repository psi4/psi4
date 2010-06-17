/*
 *  dfmp2.h
 *  matrix
 *
 *
 */

#ifndef DFMP2_H
#define DFMP2_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <libdiis/diismanager.h>
#include <psi4-dec.h>

using namespace psi;

namespace psi { namespace  dfmp2 {

class DFMP2 : public Wavefunction {
protected:

public:
protected:


public:
    DFMP2(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    /// Compute energy for the iteration.
    double compute_E();

    virtual ~DFMP2();
};

}}

#endif
