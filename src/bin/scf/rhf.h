/*
 *  rhf.h
 *  matrix
 *
 *  Created by Justin Turney on 4/10/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef RHF_H
#define RHF_H

#include <libpsio/psio.hpp>
#include "hf.h"

 using namespace psi;
 
class RHF : public HF {
protected:
    Matrix F_;
    Matrix C_;
    Matrix D_;
    Matrix Dold_;
    Matrix G_;
        
    Matrix** diis_F_;
    Matrix** diis_E_;
        
    int num_diis_vectors_;
    double **diis_B_;
    int current_diis_fock_;
    int diis_enabled_;
        
    int use_out_of_core_;
    double *pk_;
    
    void compute_multipole();
    
    void form_initialF();
    void form_C();
    void form_D();
    double compute_initial_E();
    double compute_E();
    
    void form_G(); // Out of core (i think there is a bug here)
    void form_G_from_PK(); // In core PK
    void form_G_from_direct_integrals(); // Computes all ERIs each iteration.
    
    void form_PK();
    void form_F();
    
    void find_occupation(Matrix&);
    void save_fock();
    void diis();
    void allocate_PK();
    
    bool test_convergency();
    void save_information();
    
    void common_init();
public:
    RHF(psi::PSIO *psio, psi::Chkpt *chkpt = 0);
    RHF(psi::PSIO &psio, psi::Chkpt &chkpt);
    virtual ~RHF();
    
    double compute_energy();
};

#endif
