//#include "psi4-dec.h"
#include <libmints/mints.h>
#include <libciomr/libciomr.h>
#include <liboptions/liboptions.h>
#include <libplugin/plugin.h>
#include <libchkpt/chkpt.h>
#include "dfadc.h"

INIT_PLUGIN

namespace psi{ namespace plugin_dfadc {

extern "C" int
read_options(std::string name, Options &options){
    if(name == "plugin_dfadc" || options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
        /*- The amount of memory available (in Mb) -*/
        options.add_int("MEMORY", 1000);
        /*- The Reference -*/
        options.add_str("REFERENCE", "");
        /*- The DF basis -*/
        options.add_str("RI_BASIS_ADC", "");
        /*- The convergence criterion for pole searching step -*/
        options.add_int("NEWTON_CONV", 7);
        /*- The maximum numbers of the pole searching iteration  -*/
        options.add_int("POLE_MAX", 20);
        /*- Maximu iteration number in simultaneous expansion method -*/
        options.add_int("SEM_MAX", 100);
        /*- The cutoff norm of residual vector in SEM step -*/
        options.add_int("NORM_TOL", 6);
        /*- Indicates whether three-virtual integral is processed by semi-direct algorithm -*/
        options.add_int("CUTOFF", 2);
        /*- Partial renormalization scheme for the ground state wavefunction -*/
        options.add_bool("DO_PR", false);
        /*- Whether convergence of residual norm is necessary or not -*/ 
        options.add_bool("DO_RESIDUAL", true);
        /*- Initial and the maximum dimension of the Ritz space -*/
        options.add("DIM_RITZ", new ArrayType());
    }
}

extern "C" PsiReturnType
plugin_dfadc(Options &options)
{
    tstart();

    fprintf(outfile, "\n");
    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\tDF-ADC(2): An Algebraic-Diagrammatic Construction Code\n ");
    fprintf(outfile, "\t                 with density-fitting algorithm\n");
    fprintf(outfile, "\t                      Masaaki Saitow\n");
    fprintf(outfile, "\t********************************************************\n");
    
    boost::shared_ptr<DFADC> dfadc (new DFADC());
    dfadc->compute_energy();
    
    tstop();
    
    fprintf(outfile, "\n");
    fprintf(outfile, "  ∩==\n"); 
    fprintf(outfile, "(: 3)))== kskkskkskksk         BOOOoooooOOOON!\n"); 
    fprintf(outfile, "  ∪==\n");
    // Gonzales Naitou
    fflush(outfile);

    return Success;
}
    
}} // End Namespaces
