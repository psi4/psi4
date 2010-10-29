/*
 *  Header file for SAPT0 objects 
 *  Created by Rob Parrish on 07/21/2010
 *
 */

#ifndef SAPT3BN7_H
#define SAPT3BN7_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <psi4-dec.h>

#include "sapt3bn6.h"

using namespace psi;

namespace psi { namespace sapt {

class SAPT3BN7 : public SAPT3BN6 {
private:
    double disp211_T_0(int, char *, char *, char *, int, char *, char *, 
      char *, int, char *, double *, double *, double *, int, int, int, int, 
      int, int);

    double disp220_T_0(int, char *, char *, char *, int, char *, char *, 
      char *, int, char *, double *, double *, double *, int, int, int, int, 
      int, int);

protected:
    virtual void print_header();

    void disp211_T();
    void disp220_T();

    virtual double print_results();

public:
    SAPT3BN7(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~SAPT3BN7();

    virtual double compute_energy();

};

}}

#endif
