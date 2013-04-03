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
    virtual void stability_analysis();

    //Some stuff for Ed Hohenstein's SAPT code
    // TODO: This must be removed for a conforming SCF module
    // The SAPT driver should save the three references and extract info from
    // That point
    void save_sapt_info();

    virtual void form_F();
    virtual void form_G();
    virtual void compute_orbital_gradient(bool save_fock);

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

    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }
};

}}

#endif
