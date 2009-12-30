#ifndef __math_test_uhf_h__
#define __math_test_uhf_h__

#include <libpsio/psio.hpp>
#include "hf.h"

#include <psi4-dec.h>

using namespace psi;

namespace psi { namespace scf {
    
class UHF : public HF {
protected:
    SharedMatrix Fa_, Fb_;
    SharedMatrix Da_, Db_, Dt_, Dtold_;
    SharedMatrix Ca_, Cb_;
    SharedMatrix Ga_, Gb_;
    int use_out_of_core_;

    double *p_jk_;
    double *p_k_;
    double Drms_;
    
    // The alpha and beta external potentials
    double ***Va_;
    double ***Vb_;

    void allocate_PK();
    void form_initialF();
    void form_C();
    void form_D();
    double compute_initial_E();
    double compute_E();

    void form_G();
    void form_G_from_PK();
    void form_PK();
    virtual void form_F();

    bool test_convergency();
    void save_information();
	
    void common_init();
    
    void compute_multipole();
public:
    UHF(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~UHF();

    /// Add an external potential to the alpha and beta Fock matrices
    void set_external_potential(double ***Va, double ***Vb)
                               {Va_ = Va;  Vb_ = Vb;}

    double compute_energy();
};

}}

#endif
