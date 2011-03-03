#ifndef __rohf_psi_h__
#define __rohf_psi_h__

#include <vector>
#include <libpsio/psio.hpp>
#include "hf.h"

#include <psi4-dec.h>

using namespace psi;

namespace psi { namespace scf {

class ROHF : public HF {
protected:
    SharedMatrix S_;
    SharedMatrix Feff_;
    SharedMatrix Fc_;
    SharedMatrix Fo_;
    SharedMatrix Dc_;
    SharedMatrix Do_;
    SharedMatrix Dc_old_;
    SharedMatrix Do_old_;
    SharedMatrix Gc_;
    SharedMatrix Go_;

    int use_out_of_core_;
    double *pk_;
    double *k_;             // Used in formation of _Fo

    void form_initialF();
    void form_initial_C();
    void form_C();
    void form_D();
    double compute_initial_E();
    double compute_E();

    void form_G_from_PK();
    void form_PK();
    void form_F();

//    void find_occupation(SharedMatrix);
    void save_fock();
    bool diis();
    void allocate_PK();

    bool test_convergency();

    void save_information();
    // Finalize memory/files
    virtual void finalize();

    void common_init();
public:
    ROHF(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    ROHF(Options& options, shared_ptr<PSIO> psio);
    ~ROHF();

    double compute_energy();
};

}}

#endif
