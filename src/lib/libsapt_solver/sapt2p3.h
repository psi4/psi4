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
    void ind30_amps(char *, char *, int, char *, int, char *, double **,
      double **, double **, double **, double **, double *, int, int, int,
      int);
    double exch_ind30_20();
    double exch_ind30_02();
    void ind_disp_ov(char *, char *, int, char *, char *, double **, double *,
      int, int);
    void ind_disp_ovov();
    double exch_ind_disp30_21();
    double exch_ind_disp30_12();
    double **disp30_amps(int, int, int, double *, double *, int, int, int, 
      int, int, int);
    void frzn_disp30_prep();
    void natural_orbitalify_disp30();
    double exch_disp30_20();
    double exch_disp30_02();
    double exch_disp30_22();

protected:
    virtual void print_header();
    virtual void exch_disp20();
    void ind30();
    void exch_ind30();
    void ind_disp30();
    void exch_ind_disp30();
    void disp30();
    void exch_disp30();
    void elst13();
    virtual double print_results();

public:
    SAPT2p3(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~SAPT2p3();

    virtual double compute_energy();

};

}}

#endif
