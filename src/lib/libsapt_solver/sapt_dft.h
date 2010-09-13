/*
 *  Header file for SAPT0 objects 
 *  Created by Rob Parrish on 07/21/2010
 *
 */

#ifndef SAPT_DFT_H
#define SAPT_DFT_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <psi4-dec.h>

#include "sapt0.h"

using namespace psi;

namespace psi { namespace sapt {

class SAPT_DFT : public SAPT0 {
private:
    void lr_ints();

    double **D_lambda_F(double, int, char *, int, char *, double *, int, int);
    double **DF_FDDS(double, int, char *, int, char *, double *, int, int);

protected:
    virtual void print_header();
    void df_disp20_chf();
    virtual double print_results();

public:
    SAPT_DFT(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~SAPT_DFT();

    virtual double compute_energy();

};

}}

#endif
