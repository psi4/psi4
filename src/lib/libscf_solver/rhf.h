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
#include <libmints/mints.h>
#include "hf.h"

using namespace psi;

namespace psi { namespace scf {

class RHF : public HF {
protected:
    SharedMatrix F_;
//    SharedMatrix C_;
    SharedMatrix Lref_; //Formal local orbital coefficients (from last guess)
    SharedMatrix L_; //Propagated local orbital coefficients
    SharedMatrix D_;
    SharedMatrix Dold_;
    SharedMatrix G_;
    SharedMatrix J_;
    SharedMatrix K_;
    SharedVector orbital_e_;

    //Local K stuff (sorry)
    double** I_; //Lodwin Charges
    int* domain_atoms_;
    int** domain_shell_start_;
    int** domain_shell_length_;
    int** domain_fun_start_;
    int** domain_fun_length_;
    int** fit_shell_start_;
    int** fit_shell_length_;
    int** fit_fun_start_;
    int** fit_fun_length_;
    int* primary_shell_start_;
    int* primary_shell_length_;
    int* primary_fun_start_;
    int* primary_fun_length_;
    int* aux_shell_start_;
    int* aux_shell_length_;
    int* aux_fun_start_;
    int* aux_fun_length_;
    int* domain_pairs_;
    int* domain_size_;
    int* fit_size_;
    bool* domain_changed_;
    int** atom_domains_;
    int** old_atom_domains_;
    int max_fit_size_;
    int max_domain_size_;
    int max_domain_pairs_;

    boost::shared_ptr<TwoBodyInt> eri_;

    double Drms_;

    double *pk_;

    void compute_multipole();

    bool load_or_compute_initial_C();
    void form_C();
    void form_D();
    double compute_initial_E();
    double compute_E();

    void form_G(); // Out of core (i think there is a bug here)
    void form_G_from_PK(); // In core PK
    void form_G_from_direct_integrals(); // Computes all ERIs each iteration.
    void form_G_from_RI(); //Uses two- and three- index integrals
    void form_G_from_RI_local_K(); //Uses two- and three- index integrals
    void form_G_from_J_and_K(double scale_K_by = 1.0); // Computes G from J and K
    void form_J_and_K();    // Computes J and K matrices from the ERIs

    void form_J_and_K_from_direct_integrals();
    void form_J_from_RI();
    void form_K_from_RI();

    void fully_localize_mos();
    void propagate_local_mos();
    void localized_Lodwin_charges();

    //Zillions of baby index matrices    
    void form_domain_bookkeeping();
    void free_domain_bookkeeping();
    void form_domains();

    //Some stuff for Ed Hohenstein's SAPT code
    void save_sapt_info();    
    
    //SAD Guess and propagation
    void compute_SAD_guess();

    //Save Dual Basis
    void save_dual_basis_projection();

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
