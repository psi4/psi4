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

    double E_MP2C_int_;
    double E_MP2C_delta_;
    double E_MP2_int_;
    double E_UCHF_disp_;
    double E_TDDFT_disp_;

    double **X0_A_;
    double **X0_B_;
    double **XC_A_;
    double **XC_B_;
    double **D_A_;
    double **D_B_;
    double **W_A_;
    double **W_B_;
    double **S_;
    double **J_;
    double **Jinv_;

    void free_arrays();
    void allocate_arrays();

    void compute_S();
    void compute_D();
    void compute_J();
    void compute_W();
    void compute_X_0(double omega);
    void compute_X_coup(double omega);
    double compute_UCHF_disp();
    double compute_TDDFT_disp();

protected:
    virtual void print_header();
    virtual double print_results();

public:
    SAPT_DFT(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~SAPT_DFT();

    virtual double compute_energy();

};

class OmegaQuadrature {
    protected:
        int npoints_;
        int index_;
        double* omega_;
        double* w_;
    public:
        OmegaQuadrature(int npoints);
        virtual ~OmegaQuadrature();

        double getWeight() { return w_[index_]; }
        double getOmega() { return omega_[index_]; }
        void nextPoint() { index_++; }
        void reset() { index_ = 0; }
        bool isDone() { return index_ >= npoints_; }
};

}}

#endif
