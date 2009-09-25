#ifndef __rohf_psi_h__
#define __rohf_psi_h__

#include <libpsio/psio.hpp>
#include "hf.h"

using namespace psi;

class ROHF : public HF {
protected:
    RefMatrix S_;
    RefMatrix Fc_;
    RefMatrix Fo_;
    RefMatrix Feff_;
    RefMatrix C_;
    RefMatrix Dc_;
    RefMatrix Do_;
    RefMatrix Dc_old_;
    RefMatrix Do_old_;
    RefMatrix Gc_;
    RefMatrix Go_;
    RefVector epsilon_;
    
    Ref<RefMatrix, SimpleReferenceCount, StandardArrayPolicy> diis_F_;
    Ref<RefMatrix, SimpleReferenceCount, StandardArrayPolicy> diis_E_;

    int num_diis_vectors_;
    double **diis_B_;
    int current_diis_fock_;
    int diis_enabled_;

    int use_out_of_core_;
    double *pk_;
    double *k_;             // Used in formation of _Fo

    int charge_;
    int multiplicity_;
    
    void form_initialF();
    void initial_guess();
    void form_C();
    void form_D();
    double compute_initial_E();
    double compute_E();

    void form_G();
    void form_G_from_PK();
    void form_PK();
    void form_F();

    void find_occupation(RefMatrix &);
    void save_fock();
    void diis();
    void allocate_PK();
    
    bool test_convergency();

    void save_information();
    
    void common_init();
public:
    ROHF(psi::PSIO *psio, psi::Chkpt *chkpt = 0);
    ROHF(Ref<psi::PSIO> &psio, Ref<psi::Chkpt> &chkpt);
    ~ROHF();

    double compute_energy();    
};

#endif
