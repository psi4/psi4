/*
 *  Header file for SAPT0 objects 
 *  Created by Rob Parrish on 07/21/2010
 *
 */

#ifndef SAPT2P_H
#define SAPT2P_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <psi4-dec.h>

#include "sapt2.h"

using namespace psi;

namespace psi { namespace sapt {

class SAPT2p : public SAPT2 {
private:
    virtual void theta_ar();
    virtual void theta_bs();

    double disp210();
    double disp201();
    double disp211();
    double disp22s(char *, char *, int, char *, char *, int, int);
    double disp220d();
    double disp202d();
    double disp22q_1(char *, char *, char *, int, int);
    double disp22q_2(char *, char *, char *, int, char *, int, int);
    double disp22q_3(char *, char, char, char *, char *, int, int, int, int);
    double disp22q_4(char *, char *, char *, int, int, int, int);
    double disp22q_5(char *, char *, char *, int, int, int, int);
    double disp220t(int, char *, char *, char *, int, char *, int, char *,
      char *, double *, double *, int, int, int, int, int, int);
    void fzn_triples(int, char *, char *, char *, int, char *, int, char *,
      char *, int, int, int, int, int, int); 
    void natural_orbitalify_triples(int, char *, char *, char *, int, char *,
      char *, char *, double *, double *, double **, double **, int, int, 
      int, int, int, int, int, int);
     double nat_orb_disp20();

protected:
    virtual void print_header();
    virtual void disp20();
    virtual void disp21();
    virtual void disp22sdq();
    virtual void disp22t();
    virtual double print_results();

public:
    SAPT2p(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~SAPT2p();

    virtual double compute_energy();

};

}}

#endif
