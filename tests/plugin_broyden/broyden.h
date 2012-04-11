#ifndef SRC_LIB_BroydenRHF_H
#define SRC_LIB_BroydenRHF_H    

#include <libscf_solver/rhf.h>
#include <libmints/typedefs.h>

namespace psi{

class Options;

namespace scf{


class BroydenRHF : public RHF {

protected:

    // Current Broyden iteration (zeroed on restart) 
    int broyden_iteration_;
    // Broyden status message
    std::string broyden_status_;

    // Noccpi
    Dimension nocc_;
    // Nvirpi
    Dimension nvir_;

    // Initial eigenvalues, for precondition
    SharedVector eps_occ_;
    // Initial eigenvalues, for precondition
    SharedVector eps_vir_;

    // Initial C matrix
    SharedMatrix C0_;
    // Fia residual matrix
    SharedMatrix Fia_;
    // Xia parameter matrix
    SharedMatrix Xia_;
    // Fia residual matrix
    SharedMatrix Fia_old_;
    // Xia parameter matrix
    SharedMatrix Xia_old_;
    // V = exp(X) matrix
    SharedMatrix V_;

    // S deques
    std::deque<SharedMatrix> s_;
    // S 2-norms (squared)
    std::deque<double> s2_; 

    // Y deques
    std::deque<SharedMatrix> y_;
    // S'Y dot products
    std::deque<double> p_;

public:
    BroydenRHF(Options &options);
    virtual ~BroydenRHF();
    double compute_energy();
protected:
    void bfgs_seed();
    void bfgs_step();
    void dfp_seed();
    void dfp_step();
    void sr1_seed();
    void sr1_step();
    void broyden_seed();
    void broyden_step();
    void precondition_Fia();
    void build_Fia();
    void rotate_orbitals();
};


}} // Namespaces

#endif // Header guard
