/*
 *  Header file for SAPT0 objects 
 *  Created by Rob Parrish on 07/21/2010
 *
 */

#ifndef SAPT2_H
#define SAPT2_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <psi4-dec.h>

#include "sapt0.h"

using namespace psi;

namespace psi { namespace sapt {

class SAPT2 : public SAPT0 {
private:
    double exch110(const char *);
    double exch101(const char *);
    double exch111();
    double exch120_k2f();
    double exch102_k2f();
    double exch120_k11u_1();
    double exch120_k11u_2();
    double exch120_k11u_3();
    double exch120_k11u_4();
    double exch120_k11u_5();
    double exch120_k11u_6();
    double exch102_k11u_1();
    double exch102_k11u_2();
    double exch102_k11u_3();
    double exch102_k11u_4();
    double exch102_k11u_5();
    double exch102_k11u_6();
    double ind22_1(double **, double **, double **, const char *, int, 
      const char *, const char *, const char *, double *, int, int);
    double ind22_2(double **, double **, double **, const char *, int, int);
    double ind22_3(double **, double **, const char *, const char *, int, int);
    double ind22_4(double **, const char *, int, const char *, int, int);
    double ind22_5(double **, const char *, double *, int, int);
    double ind22_6(double **, const char *, int, const char *, const char *, 
      const char *, int, int);
    double ind22_7(double **, const char *, const char *, const char *, int, 
      const char *, const char *, const char *, int, const char *, int, int, 
      int, int);

    virtual double exch_ind20respA_B();
    virtual double exch_ind20respB_A();

protected:
    virtual void print_header();
    void elst12();
    void exch11();
    void exch12();
    void ind22();
    virtual double print_results();

public:
    SAPT2(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~SAPT2();

    virtual double compute_energy();

};

}}

#endif
