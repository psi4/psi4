#include <psi4-dec.h>
#include <libmints/mints.h>
#include <liboptions/liboptions.h>
#include <libplugin/plugin.h>
#include <libpsio/psio.hpp>
#include "dfadc.h"

namespace psi{ namespace plugin_dfadc {

DFADC::DFADC(): Wavefunction(Process::environment.options, _default_psio_lib_)
{   
    //
    // Original references for the ADC(2) theory:
    //
    // * J. Schirmer, PRA 26 (1982) 2395.
    // * A. B. Trofimov, G. Stelter and J. Schirmer, JCP 117 (2002) 6402.
    //
    // and so on.
    //
    // Remarks from the other groups:
    //
    // * C. Haettig, in NIC Series, edited by J. Grotendorst, S. Bleugel, and D. Marx 
    //   (John von Neumann Institute for Computing, Julich, Switzerland, 2006), Vol. 31, p. 245.
    // * C. Haettig and K. Hald, PCCP 4 (2002) 2111.
    // * C. Haettig, Adv. Quantum Chem. 50 (2005) 37.
    // * D. Kats, D. Usvyat and M. Schuetz, PRA 83 (2011) 062503.
    //
    init();
    
    conv_      = options_.get_int("NEWTON_CONV");
    norm_tol_  = options_.get_int("NORM_TOL");
    pole_max_  = options_.get_int("POLE_MAX");
    sem_max_   = options_.get_int("SEM_MAX");
    num_roots_ = options_.get_int("NUM_ROOTS");
    cutoff_    = options_.get_int("CUTOFF");
    
    if(options_["DIM_RITZ"].size() == 2) {
      int *dim_ritz = options_.get_int_array("DIM_RITZ");
      init_ritz_ = dim_ritz[0];
      max_ritz_  = dim_ritz[1];
    }
    else if(options_["DIM_RITZ"].size() == 0) {
      init_ritz_ = 1;
      max_ritz_  = 10;
    }
    else {
      fprintf(outfile, "Length of DIM_RITZ vector must be 2!\n");
      throw PsiException("dfadc comparison error DIM_RITZ vector", __FILE__, __LINE__);
    }
    
    fprintf(outfile, "\t==> Input Parameters <==\n\n");
    fprintf(outfile, "\tNEWTON_CONV = %3d, NORM_TOL  = %3d\n", conv_, norm_tol_);
    fprintf(outfile, "\tPOLE_MAX    = %3d, SEM_MAX   = %3d\n", pole_max_, sem_max_);
    fprintf(outfile, "\tNUM_ROOTS   = %3d, CUTOFF    = %3d\n\n", num_roots_, cutoff_);
    fprintf(outfile, "\t*NXS        = %3d\n", naocc_*navir_);
    fprintf(outfile, "\tINIT_RITZ   = %3d, MAX_RITZ  = %3d\n", init_ritz_*num_roots_, max_ritz_*num_roots_);
    
    bool do_residual = options_.get_bool("DO_RESIDUAL");
    if(!do_residual) {
      fprintf(outfile, "\n\t#WARNING: Convergence of residual vector is not sought in the ADC(2) level.\n\n");
    }
}

DFADC::~DFADC()
{
}
    
}} // End Namespaces
