/*
 *  Header file for SAPT0 objects 
 *  Created by Rob Parrish on 07/21/2010
 *
 */

#ifndef SAPT2P3_H
#define SAPT2P3_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <psi4-dec.h>

#include "sapt2p.h"

using namespace psi;

namespace psi { namespace sapt {

class SAPT2p3 : public SAPT2p {
private:

protected:
    virtual void print_header();
    void disp30();
    void elst13();
    virtual double print_results();

    double disp30_1(int, const char *, int, const char *, int, const char *, 
      int, int, int, int);
    double disp30_2(int, const char *, int, const char *, const char *, int, 
      const char *, const char *, int, int, int, int);

public:
    SAPT2p3(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~SAPT2p3();

    virtual double compute_energy();

};

}}

#endif
