#include <psi4-dec.h>
#include <libmints/mints.h>
#include <liboptions/liboptions.h>
#include <libplugin/plugin.h>
#include <libpsio/psio.hpp>
#include "dfadc.h"

namespace psi{ namespace plugin_dfadc {

DFADC::DFADC(): Wavefunction(Process::environment.options, _default_psio_lib_)
{   
    init();
    
    conv_      = options_.get_int("NEWTON_CONV");
    norm_tol_  = options_.get_int("NORM_TOL");
    pole_max_  = options_.get_int("POLE_MAX");
    sem_max_   = options_.get_int("SEM_MAX");
    num_roots_ = options_.get_int("NUM_ROOTS");
    cutoff_    = options_.get_int("CUTOFF");
    
    fprintf(outfile, "\t==> Input Parameters <==\n\n");
    fprintf(outfile, "\tNEWTON_CONV = %3d, NORM_TOL = %3d\n", conv_, norm_tol_);
    fprintf(outfile, "\tPOLE_MAX    = %3d, SEM_MAX  = %3d\n", pole_max_, sem_max_);
    fprintf(outfile, "\tNUM_ROOTS   = %3d, CUTOFF   = %3d\n\n", num_roots_, cutoff_);
    fprintf(outfile, "\t*NXS        = %3d\n", naocc_*navir_);

}

DFADC::~DFADC()
{
}
    
}} // End Namespaces
