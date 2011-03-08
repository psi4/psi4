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

#define CUSTOM_PK_CODE 0

using namespace psi;

namespace psi {

class TwoBodySOInt;
class PSIO;
class Chkpt;

namespace scf {

class RHF : public HF {
protected:
    SharedMatrix D_;
    SharedMatrix Dold_;
    SharedMatrix G_;
    SharedMatrix J_;
    SharedMatrix K_;

    double *pk_;

    void form_C();
    void form_D();
    double compute_initial_E();
    virtual double compute_E();

    // Form G routines
    double **G_vector_;                                // Used in form_G_from_PK to handle threading.
    void form_G_from_direct_integrals_parallel();      // Computes all ERIs in parallel each iteration

    //Some stuff for Ed Hohenstein's SAPT code
    // TODO: This must be removed for a conforming SCF module
    // The SAPT driver should save the three references and extract info from
    // That point
    void save_sapt_info();

    // PK specific stuff
    size_t pk_size_;
    size_t pk_pairs_;
    int *pk_symoffset_;

    //Save Dual Basis
    void save_dual_basis_projection();

#if CUSTOM_PK_CODE
    void allocate_PK();
    void form_PK();
    void form_G_from_PK();
#endif

    void form_F();
    virtual void form_G();

    void save_fock();
    bool diis();

    bool test_convergency();
    void save_information();

    void common_init();

    // Finalize memory/files
    virtual void finalize();

    void save_density_and_energy();

public:
    RHF(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    RHF(Options& options, shared_ptr<PSIO> psio);
    virtual ~RHF();

    double compute_energy_parallel();

    virtual SharedMatrix Da() const;

    virtual bool restricted() const { return true; }
};

}}

#endif
