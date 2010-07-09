/*! \defgroup CINTS cints: The Integral Computation Suite */

/*!
  \file
  \ingroup CINTS
  \brief Integral Computation Module
*/
#include<cstdio>
#include<cstdlib>
#include<cmath>

#include<libipv1/ip_lib.h>
#include<libciomr/libciomr.h>
#include<libchkpt/chkpt.h>
#include<libpsio/psio.h>
#include<psifiles.h>
#include<libint/libint.h>
#include<liboptions/liboptions.h>

#include"defines.h"
#include"global.h"
#include <stdexcept>
#include <iostream>
#include"small_fns.h"
#include"parsing.h"
#include"prints.h"
#include"molecule.h"
#include"basisset.h"
#include"symm.h"
#include"dcr.h"
#include"gto.h"
#ifdef INCLUDE_Default_Ints
#include"oe_ints.h"
#include"te_ints.h"
#endif
#ifdef INCLUDE_Fock
#include"fock.h"
#endif
#ifdef INCLUDE_Default_Deriv1
#include"deriv1.h"
#endif
#ifdef INCLUDE_Default_Deriv2
void deriv1_fock(void);
#include"deriv2.h"
#endif
#ifdef INCLUDE_OEProp_Ints
#include"oeprop_ints.h"
#endif
#ifdef INCLUDE_MP2
#include"mp2.h"
#endif
#ifdef INCLUDE_MkPT2
#include"mkpt2.h"
#endif
#ifdef INCLUDE_R12_Ints
#include"r12_te_ints.h"
#endif
#ifdef INCLUDE_MP2R12
#include"mp2r12.h"
#endif
#ifdef INCLUDE_CC
#include"cc.h"
#endif
#ifdef INCLUDE_GIAO_Deriv
#include"giao_deriv.h"
#endif
#include <Tools/cdsalc.h>

namespace psi {
  namespace cints {

  /*-------------------------------
    External functions declaration
    -------------------------------*/
    void init_globals();
    void check_max_am();

//! CINTS main procedure.
PsiReturnType cints(Options & options/*, int argc, char *argv[]*/)
{
  try {
    /*--- Local variables ---*/
    int i,j,k,l,m,count;

    init_globals();
    start_io();
    parsing(options);

    /*--- Parse the command line ---*/
    parsing_cmdline(options);
    print_intro();
    setup();

    /*--- Prepare data ---*/
    init_molecule();
    init_symmetry();
    init_basisset();
    check_max_am();
    init_dcr();
    init_gto();

    /* If need to compute derivatives over SO -- compute SALC data */
    if (UserOptions.make_deriv1 && UserOptions.symm_ints)
      init_cdsalc();

    /*--- Print out some stuff ---*/
    print_scalars();
    print_basisset();

    /*--- Compute the integrals ---*/
#ifdef INCLUDE_Default_Ints
    if (UserOptions.make_oei) {
      /* Molecule.Enuc = */ compute_enuc();
      chkpt_wt_enuc(Molecule.Enuc);
      oe_ints();
    }
    if (UserOptions.make_eri)
      te_ints();
#endif
#ifdef INCLUDE_Fock
    if (UserOptions.make_fock)
      fock();
#endif
    /*--- Compute the derivative integrals ---*/
#ifdef INCLUDE_Default_Deriv1
    if (UserOptions.make_deriv1)
      deriv1();
#endif

    /*--- Compute second derivative integrals ---*/
#ifdef INCLUDE_Default_Deriv2
    if (UserOptions.make_deriv2)
      deriv2();
#endif

#ifdef INCLUDE_OEProp_Ints
    if (UserOptions.make_oeprop)
      oeprop_ints();
#endif
#ifdef INCLUDE_MP2
    if (UserOptions.make_mp2)
      mp2();
#endif
#ifdef INCLUDE_MkPT2
    if (UserOptions.make_mkpt2_ints)
      run_mkpt2();
#endif
#ifdef INCLUDE_R12_Ints
    if (UserOptions.make_r12ints)
      r12_te_ints();
#endif
#ifdef INCLUDE_MP2R12
    if (UserOptions.make_mp2r12)
      mp2r12();
#endif
#ifdef INCLUDE_CC
    if (UserOptions.make_cc_bt2)
      direct_cc();
#endif
#ifdef INCLUDE_GIAO_Deriv
    if (UserOptions.make_giao_deriv)
      giao_deriv();
#endif

    /*--- Cleanup ---*/
    if (UserOptions.make_deriv1 && UserOptions.symm_ints)
      cleanup_cdsalc();
    cleanup_gto();
    cleanup_symmetry();
    cleanup_basisset();
    cleanup_molecule();

    stop_io();
  } catch (std::exception e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "cints failed due to errors\n";
    punt(const_cast<char*>(e.what()));
  }
  return(Success);
}

}}

