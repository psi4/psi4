/*! \file
    \ingroup MP2R12
    \brief Enter brief description of file here 
*/
/* Struct for input parameters */

#ifndef _psi_src_bin_mp2r12_params_h
#define _psi_src_bin_mp2r12_params_h

namespace psi{
  namespace mp2r12{

  struct Params {
    int print_lvl;          /* Printing level */
    int c_limit;            /* Whether to use the limiting form for the B^{-1}*V */
    double tolerance;       /* Cutoff for reading in integrals */
    char *wfn;              /* The wavefunction */
    int keep_integrals;     /* Keep the integrals? */
  };

  } /* Namespace mp2r12 */
} /* Namespace psi */

#endif /* Header guard */
