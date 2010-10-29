/*
 *  Header file for SAPT0 objects 
 *  Created by Rob Parrish on 07/21/2010
 *
 */

#ifndef SAPT0_H
#define SAPT0_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <psi4-dec.h>

#include "sapt2b.h"

using namespace psi;

namespace psi { namespace sapt {

class SAPT0 : public SAPT2B {
private:
    void H3(double **);
    void Q1(double **);
    void Q3(double **);
    void Q5(double **);
    void H1(double **);
    void Q7(double **);
    void Q10(double **);
    void Q11(double **);
    void H2(double **);
    void Q6(double **);
    void Q13(double **);
    void H4(double **);
    void Q2(double **);
    void Q14(double **);

    virtual double exch_ind20respA_B();
    virtual double exch_ind20respB_A();

protected:
    virtual void print_header();
    void elst10();
    void exch10();
    virtual void disp20();
    virtual void exch_disp20();
    void ind20();
    virtual void exch_ind20();
    virtual double print_results();

    virtual void theta_ar();
    virtual void theta_bs();

    double exch_disp_1();
    double exch_disp_2();
    double exch_disp_3();
    double exch_disp_4();
    double exch_disp_5();
    double exch_disp_6();
    double exch_disp_7();

public:
    SAPT0(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~SAPT0();

    virtual double compute_energy();

};

}}

#endif
