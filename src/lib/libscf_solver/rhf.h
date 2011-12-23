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

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class TwoBodySOInt;
class PSIO;
class Chkpt;
class Matrix;
class Vector;

namespace scf {

class RHF : public HF {
protected:
    SharedMatrix D_;
    SharedMatrix Dold_;
    SharedMatrix G_;
    SharedMatrix J_;
    SharedMatrix K_;


    void form_C();
    void form_D();
    virtual void damp_update();
    double compute_initial_E();
    virtual double compute_E();

    //Some stuff for Ed Hohenstein's SAPT code
    // TODO: This must be removed for a conforming SCF module
    // The SAPT driver should save the three references and extract info from
    // That point
    void save_sapt_info();


// PK specific stuff
#if CUSTOM_PK_CODE
    // Form G routines
    double **G_vector_;                                // Used in form_G_from_PK to handle threading.
    double *pk_;
    size_t pk_size_;
    size_t pk_pairs_;
    int *pk_symoffset_;
    void allocate_PK();
    void form_PK();
    void form_G_from_PK();
#endif

    virtual void form_F();
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
    RHF(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    RHF(Options& options, boost::shared_ptr<PSIO> psio);
    virtual ~RHF();


    virtual SharedMatrix Da() const;

    virtual bool restricted() const { return true; }
};

}}

#endif
