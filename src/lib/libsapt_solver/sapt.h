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
private:
    workflow workflow_;

    void get_params();
    void get_ribasis();
    void get_calc_info();
    void cleanup_calc_info();

    // Integral Functions
    void df_ints();
    void oetrans();
    void w_ints();

    // CPHF Functions
    double **uchf_ind(double **, double *, int, int);
    double **cphf_ind(int, char *, char *, char *, double **, double **, 
      double *, int, int);
    void A_mat(int, char *, char *, char *, double **, double **, int, int,
      int);
    void diis_update(double **, double **, double **, int, int);

    // Amplitude functions
    void t_arar(int, int, int);
    void t_bsbs(int, int, int);
    void t_arbs(int);
    void t2_arar(int);
    void t2_bsbs(int);
    void Y2(char *, char *, char *, char *, char *, int, char *, char *, 
      char *, double *, int, int);
    double *t2_solver(char *, char *, int, char *, char *, char *, double *,
      int, int);
    void g_arar();
    void g_bsbs();

    // Natural Orbital Functions
    void natural_orbitalify(char *, char *, double *, double **, int, int,
      char);

protected:
    params params_;
    calcinfo calc_info_;
    results results_;
    noinfo no_info_;
    shared_ptr<BasisSet> ribasis_;
    shared_ptr<BasisSet> zero_;

    virtual void print_header()=0;
    virtual double print_results()=0;

    void compute_integrals();
    void compute_amplitudes();
    void cphf_induction();

    void zero_disk(int, char *, char *, int, int);
    double **read_IJKL(int, char *, int, int);
    void write_IJKL(double **, int, char *, int, int);
    double **get_DF_ints(int,char *,int);
    double **get_diag_AA_ints(int);
    double **get_diag_BB_ints(int);
    double **get_AA_ints(int);
    double **get_BB_ints(int);
    double **get_AB_ints(int);
    double **get_AS_ints(int);
    double **get_RB_ints(int);
    double **get_AR_ints(int);
    double **get_BS_ints(int);
    double **get_RR_ints(int);
    double **get_SS_ints(int);

public:
    SAPT(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~SAPT();

    virtual double compute_energy()=0;
};

}}

#endif
