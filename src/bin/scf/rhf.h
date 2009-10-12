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

namespace psi { namespace scf {
     
class RHF : public HF {
protected:
    SharedMatrix F_;
    SharedMatrix C_;
    SharedMatrix D_;
    SharedMatrix Dold_;
    SharedMatrix G_;
    SharedMatrix J_;
    SharedMatrix K_;
        
    std::vector<SharedMatrix> diis_F_;
    std::vector<SharedMatrix> diis_E_;
    double Drms_;
    
    int num_diis_vectors_;
    double **diis_B_;
    int current_diis_fock_;
    int diis_enabled_;
        
    int use_out_of_core_;
    double *pk_;
    
    int mind_; //minimum sorted compound index
    int* sieve_ind_; //permutation matrix for sort
    int* cut_ind_; //critical sorted index for schwarz sieve
    double* norm_;
    
    double **B_ia_P_; //Three Index tensor for DF-SCF
    
    void compute_multipole();
    
    void form_initialF();
    void form_C();
    void form_D();
    double compute_initial_E();
    double compute_E();
    
    void schwarz_sieve();
    void form_B(); 
    
    
    void form_G(); // Out of core (i think there is a bug here)
    void form_G_from_PK(); // In core PK
    void form_G_from_direct_integrals(); // Computes all ERIs each iteration.
    void form_G_from_direct_integrals_schwarz(); // Computes all ERIs  with schwarz seive each iteration.
    void form_G_from_RI(); //Uses two- and three- index integrals
    void form_G_from_J_and_K(double scale_K_by = 1.0); // Computes G from J and K
    void form_J_and_K();    // Computes J and K matrices from the ERIs
    
    void form_PK();
    void form_F();
    
    void find_occupation(SharedMatrix);
    void save_fock();
    void diis();
    void allocate_PK();
    
    bool test_convergency();
    void save_information();
    
    void common_init();
public:
    RHF(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~RHF();
    
    double compute_energy();
};

}}

#endif
