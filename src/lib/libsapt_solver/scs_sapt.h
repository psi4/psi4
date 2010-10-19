/*
 *  Header file for SAPT0 objects 
 *  Created by Rob Parrish on 07/21/2010
 *
 */

#ifndef SCS_SAPT_H
#define SCS_SAPT_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <psi4-dec.h>

#include "sapt0.h"

using namespace psi;

namespace psi { namespace sapt {

class SCS_SAPT : public SAPT0 {
private:

protected:
    virtual void print_header();
    virtual double print_results();

    void scale_spin_components();

    virtual void exch_disp20();

public:
    SCS_SAPT(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~SCS_SAPT();

    virtual double compute_energy();

};

}}

#endif
