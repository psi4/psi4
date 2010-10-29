/*
 *  Header file for SAPT objects 
 *  Created by Rob Parrish on 07/21/2010
 *
 */

#ifndef SAPT3B_H
#define SAPT3B_H

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


class SAPT3B : public SAPT {
private:
    void get_calc_info();
    void cleanup_calc_info();

    // Integral Functions
    virtual void df_ints();
    virtual void oetrans();
    virtual void w_ints();

    void build_w(double **, double **, double **, double **, int, char *, 
      char *, char *, int, char *, int, int, int);

    // Amplitude Functions
    void t2_amps(int, char *, int, char *, int, char *, double *, double *,
      int, int, int, int);
    void T_calc(char *, char, int,char *, char *, int, int, int, int);
    void theta_amps(int, char *, char *, int, char *, char *, char *,
      double *, int, int);
    void S_amps(char *, char *, int, char *, char *, int, int);
    void K1_amps(char *, int, char *, int, char *, char *, int, char *, char *,
      int, int, int, int, int, int);
    void K2_amps(char *, char *, char *, int, char *, int, char *, int, int,
      int, int);
    void K3_amps(char *, char *, int, char *, int, int);

protected:
    three_body_info calc_info_;
    three_body_results results_;

    virtual void print_header()=0;
    virtual double print_results()=0;

    void compute_amplitudes();
    void cphf_induction();

public:
    SAPT3B(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~SAPT3B();

    virtual double compute_energy()=0;
};

}}

#endif
