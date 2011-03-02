#ifndef __math_test_uhf_h__
#define __math_test_uhf_h__

#include <libpsio/psio.hpp>
#include "hf.h"

#include <psi4-dec.h>

using namespace psi;

namespace psi { namespace scf {

class UHF : public HF {
protected:
    SharedMatrix Da_, Db_, Dt_, Dtold_;
    SharedMatrix Ga_, Gb_, Ka_, Kb_;

    double *p_jk_;
    double *p_k_;

    void allocate_PK();
    void form_initialF();
    void form_C();
    void form_D();
    double compute_initial_E();
    double compute_E();

    void form_G_from_PK();
    void form_PK();
    virtual void form_F();
    virtual bool load_or_compute_initial_C();

    void save_fock();
    void diis();

    bool test_convergency();
    void save_information();

    void common_init();
    // Finalize memory/files
    virtual void finalize();

public:
    UHF(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    UHF(Options& options, shared_ptr<PSIO> psio);
    virtual ~UHF();

    double compute_energy();

    virtual bool restricted() const { return false; }
};

}}

#endif
