/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_bin_detci_ci_tol_h
#define _psi_src_bin_detci_ci_tol_h

namespace psi { namespace detci {

#define S_MAX 0.999995 /* Max overlap for mitrush vects and collapsed vecs */
#define HD_MIN 1.0E-4  /* Minimum diagonal element of preconditioner       */
#define ZERO 1e-10
#define MPn_ZERO 1e-14
#define SA_NORM_TOL 1.0E-5 /* norm of schmidt orthogonalized vector */
#define MPn_NORM_TOL 1.0E-12

}} // namespace psi::detci

#endif // header guard
