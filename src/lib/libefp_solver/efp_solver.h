#ifndef EFP_SOLVER_H
#define EFP_SOLVER_H
/*
 * EFP header
 */

// TODO: add -I to compilation
//#include"/Users/loriab/linux/libefp/src/efp.h"
#include "../../../libefp/install/include/efp.h"

namespace psi{
  class Options;
}

namespace boost {
template<class T> class shared_ptr;
}

namespace psi{ namespace efp{


class EFP {
    // warning: options_ is pointer to current options object, and may not reflect
    // proper efp options outside of common_init()
    Options & options_;
    protected:
        char ** frag_name_;
        int nfrag_;
        struct efp * efp_;
        boost::shared_ptr<Molecule>molecule_;
        bool elst_enabled_, pol_enabled_, disp_enabled_, exch_enabled_, do_grad_;
        /// Initialize options
        void common_init();
    public:
        EFP(Options& options);
        ~EFP();
  
        /// Set geometry
        void SetGeometry();

        /// Compute energy and/or gradietn
        void Compute();

        /// Designate which atoms are qm atoms
        void SetQMAtoms();

        /// Returns EFP contribution to SCF energy
        double scf_energy_update();

        /// Returns EFP contribution to V
        boost::shared_ptr<Matrix> modify_Fock();
};

}}

#endif
