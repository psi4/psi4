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
  //class PSIO;
  class Molecule;
}

namespace boost {
template<class T> class shared_ptr;
}

namespace psi{ namespace efp{


class EFP {
    Options & options_;
    //boost::shared_ptr<PSIO> psio_;
    protected:
        char ** frag_name;
        int nfrag;
        struct efp * efp;
        boost::shared_ptr<Molecule>molecule;
        bool elst_enabled, pol_enabled, disp_enabled, exch_enabled, do_grad;
    public:
        EFP(Options& options);
        ~EFP();
  
        /// Initialize options
        void common_init();
        /// Set geometry
        void SetGeometry();
        /// Compute energy and/or gradietn
        void Compute();
};

}}

#endif
