#include "psi4-dec.h"
#include <libmints/mints.h>
#include <liboptions/liboptions.h>
#include <libplugin/plugin.h>
#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include "apps.h"

INIT_PLUGIN

namespace psi{ namespace libfock_plugin {

extern "C" int
read_options(std::string name, Options &options){

    if(name == "PLUGIN_FOCK" || options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
        /*- The amount of debug information printed
            to the output file -*/
        options.add_int("DEBUG", 0);
        /*- Do singlet states? Default true
         -*/
        options.add_bool("DO_SINGLETS", true);
        /*- Do triplet states? Default true
         -*/
        options.add_bool("DO_TRIPLETS", true);
        /*- Minimum singles amplitude to print in 
            CIS analysis
         -*/
        options.add_double("CIS_AMPLITUDE_CUTOFF", 0.15);
        /*- Memory safety factor for allocating JK
        -*/
        options.add_double("CIS_MEM_SAFETY_FACTOR",0.75);
    } 
    if(name == "JK" || options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
        /*- The amount of debug information printed
            to the output file -*/
        options.add_int("DEBUG", 0);
        /*- The maximum number of integral threads (0 for omp_get_max_threads()) 
         -*/
        options.add_int("OMP_N_THREAD", 0);
        /*- The schwarz cutoff value 
         -*/
        options.add_double("SCHWARZ_CUTOFF", 1.0E-12);
        /*- The maximum reciprocal condition allowed in the fitting metric 
         -*/
        options.add_double("FITTING_CONDITION", 1.0E-12);
        /*- Fitting algorithm (0 for old, 1 for new)
         -*/
        options.add_int("FITTING_ALGORITHM", 0);
        /*- SCF Type 
         -*/
        options.add_str("SCF_TYPE", "DIRECT", "DIRECT DF GPUDF");
        /*- Auxiliary basis for SCF 
         -*/
        options.add_str("RI_BASIS_SCF", ""); 
    } 
    if(name == "SOLVER" || options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
        /*- The amount of debug information printed
            to the output file -*/
        options.add_int("DEBUG", 0);
        /*- Solver maximum iterations
         -*/
        options.add_int("SOLVER_MAXITER",100);
        /*- Solver convergence threshold (max 2-norm) 
         -*/
        options.add_double("SOLVER_CONVERGENCE",1.0E-6);
        /*- DL Solver number of roots 
         -*/
        options.add_int("SOLVER_N_ROOT",1);
        /*- DL Solver number of guesses
         -*/
        options.add_int("SOLVER_N_GUESS",1);
        /*- DL Solver number of subspace vectors to collapse to
         -*/
        options.add_int("SOLVER_MIN_SUBSPACE",2);
        /*- DL Solver maximum number of subspace vectors 
         -*/
        options.add_int("SOLVER_MAX_SUBSPACE",6);
        /*- DL Solver minimum corrector norm to add to subspace
         -*/
        options.add_double("SOLVER_NORM",1.0E-6);
    } 
}

extern "C" PsiReturnType
plugin_libfock(Options &options)
{
    tstart();

    boost::shared_ptr<RCIS> cis (new RCIS());
    cis->compute_energy();

    tstop();

    return Success;
}

}} // End Namespaces
