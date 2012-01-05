#include <libmints/mints.h>
#include <libciomr/libciomr.h>
#include <liboptions/liboptions.h>
#include <libchkpt/chkpt.h>
#include "psi4-dec.h"
#include "adc.h"

namespace psi{ namespace adc {
    
PsiReturnType
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
