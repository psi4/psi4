#ifndef __math_test_cuhf_h__
#define __math_test_cuhf_h__

#include <libpsio/psio.hpp>
#include "hf.h"

namespace boost {
template<class T> class shared_ptr;
}

namespace psi { namespace scf {

/*

    Constrained Unrestricted Hartree-Fock

    Reference: T. Tsuchimochi and G.E. Scuseria, J. Chem. Phys. 133, 
               141102 (2010)

    This is an alternative formulation of ROHF as a contrained UHF. A 
    Lagrangian constraint is placed on the usual UHF procedure to remove 
    the spin contamination. The result is an ROHF energy and semicanonical 
    ROHF orbitals. The need to pick coupling coefficients is removed. 
    Koopmans' theorem is valid for CUHF orbital energes.

    It is claimed that CUHF does not suffer from the convergence problems
    of certain ROHF implementations (not sure how PSI's ROHF code does).
    CUHF retains the UHF-like trait that Ca != Cb. Also, the converged CUHF 
    wavefunction yields the correct value for <S^2>, however, this is only 
    true at convergence. It is possible that this increased flexibility 
    improves convergence.

    -- EGH, August 15th, 2011

    TODO:

    Probably can't handle NSO != NMO right now, should either fix this code 
    or the transform functions from Matrix.

    Using the UHF form for the Lagrangian, this is probably correct, but 
    should be checked.

*/

class CUHF : public HF {
protected:
    SharedMatrix Dt_, Dtold_;
    SharedMatrix J_, Ka_, Kb_;
    // Contributions to the Fock matrix from charge and spin density
    SharedMatrix Fp_, Fm_;
    // Charge denisty and natural orbitals (eigenvectors of charge density)
    SharedMatrix Dp_, Cno_, Cno_temp_;
    // Natural orbital occupations
    SharedVector No_;

    void form_initialF();
    void form_C();
    void form_D();
    double compute_initial_E();
    virtual void stability_analysis();
    virtual double compute_E();

    virtual void form_G();
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
    CUHF(Options& options, boost::shared_ptr<PSIO> psio, 
        boost::shared_ptr<Chkpt> chkpt);
    CUHF(Options& options, boost::shared_ptr<PSIO> psio);
    virtual ~CUHF();

    virtual bool same_a_b_orbs() const { return false; }
    virtual bool same_a_b_dens() const { return false; }
};

}}

#endif
