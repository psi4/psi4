#ifndef _psi_include_chkpt_params_h_
#define _psi_include_chkpt_params_h_

// Up to k-functions (L=7) currently supported
const int MAXANGMOM = 11;
const double LINDEP_CUTOFF = 1.0e-6;
//#define MAXANGMOM 11          /* Up to k-functions (L=7) currently supported */
//#define LINDEP_CUTOFF 1E-6    /* Threshold below which basis functions are considered
//                                 linearly dependent (see canonical orthogonalization in Szabo in Ostlund) */

#endif /* header guard */
