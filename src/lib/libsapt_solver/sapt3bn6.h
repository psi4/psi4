/*
 *  Header file for SAPT0 objects 
 *  Created by Rob Parrish on 07/21/2010
 *
 */

#ifndef SAPT3BN6_H
#define SAPT3BN6_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <psi4-dec.h>

#include "sapt3bn5.h"

using namespace psi;

namespace psi { namespace sapt {

class SAPT3BN6 : public SAPT3BN5 {
private:
    double disp3100_1(char *, char *, char *, char *, char, char *, int, 
      char *, int, char *, double *, double *, int, int, int, int);
    double disp3100_2(char *, char *, char *, int, int);
    double disp3100_3(char *, char, char *, char, char *, char, int, char *,
      int, int, int, int, int, int);
    double disp3100_4(char *, char *, char *, int, int, int, int);

    double disp211_D_1(char *, char *, double *, double *, int, int, int, int);

    double disp220_S_1(char *, char *, double *, int, int);

    double disp220_D_1(char *, char *, double *, int, int);

    double disp220_Q_1(char *, char *, char *, int, int, int, int);
    double disp220_Q_2(char *, char *, int, char *, int, int, int, int);
    double disp220_Q_3(char *, char *, int, char *, int, int, int, int);

protected:
    virtual void print_header();

    virtual double print_results();

public:
    SAPT3BN6(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~SAPT3BN6();

    void disp3100();
    void disp211_D();
    void disp220_S();
    void disp220_D();
    void disp220_Q();

    virtual double compute_energy();

};

}}

#endif
