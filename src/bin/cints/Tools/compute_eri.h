#ifndef _psi_src_bin_cints_Tools_compute_eri_h
#define _psi_src_bin_cints_Tools_compute_eri_h

/*! \file compute_eri.h
    \ingroup CINTS
*/
#include <libint/libint.h>

namespace psi { 
  namespace CINTS {
#ifndef __cplusplus
#error "Tools/compute_eri.h cannot be used in C programs"
#endif
    
    //! Returns the number of integrals computed
    int compute_eri(double* target, Libint_t* Libint, int& si, int& sj, int& sk, int& sl,
		    int& inc1, int& inc2, int& inc3, int& inc4, const bool do_not_permute);
  }
}
#endif
