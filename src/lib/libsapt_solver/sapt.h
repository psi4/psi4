/*
 *  Header file for SAPT objects 
 *  Created by Rob Parrish on 07/21/2010
 *
 */

#ifndef SAPT_H
#define SAPT_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <psi4-dec.h>

#include "structs.h"
#define INDEX(i,j) ((i>=j) ? (calc_info_.ioff[i] + j) : (calc_info_.ioff[j] + i))


using namespace psi;

namespace psi { namespace sapt {


class SAPT : public Wavefunction {
protected:
    params params_;
    calcinfo calc_info_;
    results results_;
    shared_ptr<BasisSet> ribasis_;
    shared_ptr<BasisSet> zero_;

public:
    SAPT(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);

    void setup_sapt();
    void get_params();
    virtual void print_header();
    void get_ribasis();
    void get_calc_info();
    void cleanup_calc_info();
    void oetrans();
    void strans();
    void vtrans();
    void lr_ints();
    double** get_DF_ints(int filenum, char* label, int length);
    void zero_disk(int a, char* b, char* c, int, int);
    double CHF(int, char *, char *, char *, double **, double **, double *, int, int);
    void diis_update(double **, double **, double **, int, int);
    void A_mat(int, char *, char *, char *, double **, double **, int, int, int);
    double **W_ints(int, char *, double *, double **, int, int);


    virtual ~SAPT();
};

}}

#endif
