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

    shared_ptr<BasisSet> quadbasis_;
    double C_x_; // Slater's Constant

    double E_MP2C_int_;
    double E_MP2C_delta_;
    double E_MP2_int_;
    double E_UCHF_disp_;
    double E_TDDFT_disp_;

    double **X0_A_;
    double **X0_B_;
    double **XC_A_;
    double **XC_B_;
    double **WA_;
    double **WB_;
    double **J_;
    double **Jinv_;
    double **Q_;
    double **Qinv_;
    double* d_A_;
    double* d_B_;

    void free_arrays();
    void allocate_arrays();
    void get_quad_basis();

    void compute_J();
    void compute_W();
    void compute_Q();
    void compute_d();
    void compute_X_0(double omega);
    void compute_X_coup(double omega);
    double compute_UCHF_disp();
    double compute_TDDFT_disp();

protected:
    virtual void print_header();
    virtual void disp20();
    virtual void exch_disp20();
    virtual double print_results();

public:
    SAPT_DFT(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~SAPT_DFT();

    virtual double compute_energy();

};


}}

#endif
