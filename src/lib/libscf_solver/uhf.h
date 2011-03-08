#ifndef __math_test_uhf_h__
#define __math_test_uhf_h__

#include <libpsio/psio.hpp>
#include "hf.h"

#include <psi4-dec.h>

using namespace psi;

namespace psi { namespace scf {

class UHF : public HF {
protected:
    SharedMatrix Dt_, Dtold_;
    SharedMatrix Ga_, Gb_, J_, Ka_, Kb_;

//    void allocate_PK();
    void form_initialF();
    void form_C();
    void form_D();
    double compute_initial_E();
    virtual double compute_E();

    virtual void form_G();
//    void form_G_from_PK();
//    void form_PK();
    virtual void form_F();

    void save_fock();
    bool diis();

    bool test_convergency();
    void save_information();
    void compute_spin_contamination();

    void common_init();

    void save_density_and_energy();

    // Finalize memory/files
    virtual void finalize();

public:
    UHF(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    UHF(Options& options, shared_ptr<PSIO> psio);
    virtual ~UHF();

    virtual bool restricted() const { return false; }
};

}}

#endif
