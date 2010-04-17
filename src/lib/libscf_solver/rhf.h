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

    boost::shared_ptr<TwoBodyInt> eri_;

    double Drms_;

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
    void form_G_from_RI(); //Uses two- and three- index integrals
    void form_G_from_J_and_K(double scale_K_by = 1.0); // Computes G from J and K
    void form_J_and_K();    // Computes J and K matrices from the ERIs

    void form_J_and_K_from_direct_integrals();
    void form_J_from_RI();
    void form_K_from_RI();

    void form_PK();
    void form_F();

    void save_fock();
    void diis();
    void allocate_PK();

    //DOWN FOR MAINTENANCE
    //void save_RHF_grid(Options &options, shared_ptr<BasisSet>  bas, SharedMatrix D, SharedMatrix C);    
    //double *getCartesianGridExtents(Options &options, shared_ptr<Molecule> m);
    //int* getCartesianGridResolution(Options &options);    
    
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
