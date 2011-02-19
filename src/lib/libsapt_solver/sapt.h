/*
 *  Header file for SAPT objects 
 *  Created by Rob Parrish on 07/21/2010
 *
 */

#ifndef SAPT_H
#define SAPT_H

//#define _MKL

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
    void get_params();

    // Integral Functions
    virtual void df_ints()=0;
    virtual void oetrans()=0;
    virtual void w_ints()=0;

protected:
    shared_ptr<BasisSet> ribasis_;
    shared_ptr<BasisSet> zero_;
    params params_;

    virtual void print_header()=0;
    virtual double print_results()=0;

    void get_ribasis();
    void compute_integrals();

    void zero_disk(int, char *, char *, int, int);
    double **read_IJKL(int, char *, int, int);
    void write_IJKL(double **, int, char *, int, int);
    double **get_DF_ints(int,char *,int);
    double **IJKL_ints(int, char *, int, int, char *, int);
    double **IJIJ_ints(int, char *, int);

    void MO_NO_ov_DF_trans(int, int, char *, char *, int, int, int, double **);
    void MO_NO_vv_DF_trans(int, int, char *, char *, int, int, int, double **);

    double **uchf_ind(double **, double *, int, int);
    double **cphf_ind(int, char *, char *, char *, double **, double **,
      double *, int, int);
    void A_mat(int, char *, char *, char *, double **, double **, int, int,
      int);
    void diis_update(double **, double **, double **, int, int);

public:
    SAPT(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~SAPT();

    virtual double compute_energy()=0;
};

}}

#endif
