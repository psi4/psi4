#include <libmints/mints.h>
#include <libciomr/libciomr.h>
#include <liboptions/liboptions.h>
#include <libchkpt/chkpt.h>
#include <libplugin/plugin.h>
#include "psi4-dec.h"
#include "adc.h"

INIT_PLUGIN

namespace psi{ namespace adc {

extern "C" int
read_options(std::string name, Options &options){
if(name == "ADC") {
    /*- The amount of information printed to the output file -*/
    options.add_int("PRINT", 1);
    /*- How to cache quantities within the DPD library -*/
    options.add_int("CACHELEV", 2);
    /*- The amount of memory available (in Mb) -*/
    options.add_int("MEMORY", 1000);
    /*- The Reference -*/
    options.add_str("REFERENCE", "");
    /*- The convergence criterion for pole searching step -*/
    options.add_int("NEWTON_CONV", 7);
    /*- The maximum numbers of the pole searching iteration  -*/
    options.add_int("POLE_MAX", 20);
    /*- Maximu iteration number in simultaneous expansion method -*/
    options.add_int("SEM_MAX", 30);
    /*- The cutoff norm of residual vector in SEM step -*/
    options.add_int("NORM_TOL", 6);
    /*- The poles per irrep vector -*/
    options.add("STATES_PER_IRREP", new ArrayType());
    /*- Partial renormalization scheme for the ground state wavefunction -*/
    options.add_bool("PR", false);
}
}
    
extern "C" PsiReturnType
adc(Options &options)
{
    tstart();

    fprintf(outfile, "\n");
    fprintf(outfile, "\t****************************************\n");
    fprintf(outfile, "\t                 A D C                  \n");
    fprintf(outfile, "\t An Algebraic-Diagrammatic Construction \n");
    fprintf(outfile, "\t based on direct-product decomposition  \n");
    fprintf(outfile, "\t             Masaaki Saitow             \n");
    fprintf(outfile, "\t****************************************\n");
    fflush(outfile);

    boost::shared_ptr<ADC> adc (new ADC());
    adc->compute_energy();
    
    tstop();
    
    fprintf(outfile, "\n");
    fprintf(outfile, "  ∩==\n"); 
    fprintf(outfile, "(: 3)))== kskkskkskksk         BOOOoooooOOOON!\n"); 
    fprintf(outfile, "  ∪==\n"); 
    fflush(outfile);

    return Success;
}
    
}} // End Namespaces
