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

#include "sapt.h"

using namespace psi;

namespace psi { namespace sapt {

class SAPT0 : public SAPT {
protected:
    double exch_disp_1();
    double exch_disp_2();
    double exch_disp_3();
    double exch_disp_4();
    double exch_disp_5();
    double exch_disp_6();
    double exch_disp_7();

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
    
    double **get_AA_ints(int dress) ;
    double **get_diag_AA_ints(int dress) ;
    double **get_BB_ints(int dress) ;
    double **get_diag_BB_ints(int dress) ;
    double **get_AB_ints(int dress) ;
    double **get_AS_ints(int dress) ;
    double **get_RB_ints(int dress) ;
    double **get_AR_ints(int dress) ;
    double **get_BS_ints(int dress) ;

    void ao_df_ints();
    void ao_df_ints_restart();
public:
    SAPT0(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);

    double compute_energy();
    void elst10();
    void exch10();
    void disp20();
    void theta_ar();
    void theta_bs();
    void exch_disp20();
    void ind20resp();   
    void exch_ind20respA_B();   
    void exch_ind20respB_A();   
    void print_results(); 
    virtual ~SAPT0();
};

}}

#endif
