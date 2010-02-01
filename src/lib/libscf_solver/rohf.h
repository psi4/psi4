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
    SharedMatrix Fc_;
    SharedMatrix Fo_;
    SharedMatrix Feff_;
    SharedMatrix C_;
    SharedMatrix Dc_;
    SharedMatrix Do_;
    SharedMatrix Dc_old_;
    SharedMatrix Do_old_;
    SharedMatrix Gc_;
    SharedMatrix Go_;
    SharedVector epsilon_;
    
    std::vector<SharedMatrix> diis_F_;
    std::vector<SharedMatrix> diis_E_;

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
    void form_G_from_RI();
    void form_G_from_direct_integrals();
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
    ROHF(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    ~ROHF();

    double compute_energy();    
};

}}

#endif
