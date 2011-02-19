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
    double disp210();
    double disp201();
    double disp211();
    double disp22s(const char *, const char *, int, const char *, const char *,
      int, int);
    double disp220d();
    double disp202d();
    double disp22q_1(const char *, const char *, const char *, int, int);
    double disp22q_2(const char *, const char *, const char *, int, 
      const char *, int, int);
    double disp22q_3(const char *, const char, const char, const char *, 
      const char *, int, int, int, int);
    double disp22q_4(const char *, const char *, const char *, int, int, 
      int, int);
    double disp22q_5(const char *, const char *, const char *, int, int, 
      int, int);
    double disp220t(int, const char *, const char *, const char *, int, 
      const char *, int, const char *, const char *, double *, double *, 
      int, int, int, int, int, int);
    void fzn_triples(int, const char *, const char *, const char *, int, 
      const char *, int, const char *, const char *, int, int, int, int, 
      int, int); 
    void natural_orbitalify_triples(int, const char *, const char *, 
      const char *, int, const char *, const char *, const char *, double *, 
      double *, double **, double **, int, int, int, int, int, int, int, int);
    double nat_orb_disp20();

protected:
    virtual void print_header();
    virtual void disp20();
    virtual void disp21();
    virtual void disp22sdq();
    virtual void disp22t();
    virtual double print_results();

    virtual void theta_ar();
    virtual void theta_bs();

public:
    SAPT2p(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~SAPT2p();

    virtual double compute_energy();

};

}}

#endif
