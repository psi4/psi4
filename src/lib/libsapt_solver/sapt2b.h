/*
 *  Header file for SAPT objects 
 *  Created by Rob Parrish on 07/21/2010
 *
 */

#ifndef SAPT2B_H
#define SAPT2B_H

//#define _MKL

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <psi4-dec.h>

#include "structs.h"
#include "sapt.h"

#define INDEX(i,j) ((i>=j) ? (calc_info_.ioff[i] + j) : (calc_info_.ioff[j] + i))


using namespace psi;

namespace psi { namespace sapt {


class SAPT2B : public SAPT {
private:
    workflow workflow_;

    void get_workflow();
    void get_calc_info();
    void cleanup_calc_info();

    // Integral Functions
    virtual void df_ints();
    virtual void oetrans();
    virtual void w_ints();

    // Amplitude functions
    void t_arar(int, int, int);
    void t_bsbs(int, int, int);
    void t_arbs(int);
    void t2_arar(int);
    void t2_bsbs(int);
    void Y2(const char *, const char *, const char *, const char *, 
      const char *, int, const char *, const char *, const char *, double *, 
      int, int);
    double *t2_solver(int, const char *, const char *, int, const char *, 
      const char *, const char *, double *, int, int, int);
    double *t2_solver_natorbs(int, const char *, const char *, int, 
      const char *, const char *, const char *, double *, int, int, int, 
      const char *, double **, int);
    void g_arar();
    void g_bsbs();
    void Y3_ar();
    void Y3_bs();
    void Y3_1(double **, int, const char *, const char *, const char *, int, 
      const char *, int, int);
    void Y3_2(double **, int, const char *, const char *, const char *, int, 
      const char *, const char *, const char *, const char *, int, int);
    void Y3_3(double **, int, const char *, int, const char *, const char *, 
      int, int);
    void Y3_4(double **, int, const char *, const char *, int, const char *, 
      int, int);
    void Y3_5(double **, int, const char *, const char *, const char *, int, 
      const char *, const char *, int, int);
    void Y3_6(double **, int, const char *, const char *, const char *, int, 
      const char *, int, int);

    // Natural Orbital Functions
    void natural_orbitalify(const char *, const char *, double *, double **, 
      int, int, const char);

protected:
    calcinfo calc_info_;
    results results_;
    noinfo no_info_;

    virtual void print_header()=0;
    virtual double print_results()=0;

    void compute_amplitudes();
    void cphf_induction();

    double **get_diag_AA_ints(const int);
    double **get_diag_BB_ints(const int);
    double **get_AA_ints(const int);
    double **get_BB_ints(const int);
    double **get_AB_ints(const int);
    double **get_AS_ints(const int);
    double **get_RB_ints(const int);
    double **get_AR_ints(const int);
    double **get_BS_ints(const int);
    double **get_RR_ints(const int);
    double **get_SS_ints(const int);

public:
    SAPT2B(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~SAPT2B();

    virtual double compute_energy()=0;
};

}}

#endif
