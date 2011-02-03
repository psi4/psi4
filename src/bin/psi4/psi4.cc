/*! \file psi4.cc
\defgroup PSI4 The new PSI4 driver.
*/

//  PSI4 driver
//  Justin Turney
//  Rollin King

#if HAVE_MPI == 1
#include <mpi.h>
#endif
#include <getopt.h>
#include <stdio.h>
#include <psiconfig.h>
#include <libciomr/libciomr.h>
#include <liboptions/liboptions.h>
#include <libparallel/parallel.h>
#include <libmints/wavefunction.h>
#include "script.h"
#include <physconst.h>
#include <psifiles.h>
#include <libyeti/yeti.h>
//#include <molecular_system.h>

#include <psi4-def.h>
#define MAIN
#include "psi4.h"
#include "script.h"

namespace psi {
    int psi_start(int argc, char *argv[]);
    int psi_stop(FILE* infile, FILE* outfile, char* psi_file_prefix);
    void print_version(FILE *);
    void set_memory(FILE *outfile);
    int psi4_driver();
    void psiclean(void);

    PSIO *psio = NULL;
}

// This is the ONLY main function in PSI
int main(int argc, char **argv, char **envp)
{
    using namespace psi;

    // Setup the environment
    Process::arguments.init(argc, argv);
    Process::environment.init(envp);
    Process::environment.set_memory(256000000L);

    std::string communicator = Process::environment("COMMUNICATOR");

#if HAVE_MPI == 1

    #if HAVE_MADNESS == 1
        if(communicator == "MADNESS") {
            madness::initialize(argc, argv);
            Communicator::world = shared_ptr<Communicator>(new MadCommunicator(
                    shared_ptr<madness::World>(new madness::World(MPI::COMM_WORLD)) ));
        }
        else {
    #endif
            // Initialize MPI
            // Since we allow multiple MPICommunicators to be created and MPI_Init needs to be called
            // as soon as possible, main takes care of call MPI_Init and not the communicator itself.
            MPI_Init(&argc, &argv);

            // Create the global world communicator
            Communicator::world = shared_ptr<Communicator>(new MPICommunicator(MPI_COMM_WORLD));

    #if HAVE_MADNESS == 1
        }
    #endif

#else
    // Initialize local communicator
    Communicator::world = shared_ptr<Communicator>(new LocalCommunicator);
#endif

    Wavefunction::initialize_singletons();

    // Create the scripting object
    Script::language = shared_ptr<Script>(new Python);
    // Create base objects in the scripting language and initialize the language
    Script::language->initialize();

    if(psi_start(argc, argv) == PSI_RETURN_FAILURE) return EXIT_FAILURE;

    if(!clean_only) print_version(outfile);

    set_memory(outfile);

    psio_init();

    if(clean_only) {
        psiclean();
        exit(EXIT_SUCCESS);
    }

    // Okay, we might only need to make this function call if we're using IPV1
    // TODO: Need to come up with a way to interface with optking
    // 
    Script::language->run(infile);

    // Shut things down:
    psi_stop(infile, outfile, psi_file_prefix);
    Script::language->finalize();

#if HAVE_MPI == 1
    #if HAVE_MADNESS == 1
        if(communicator == "MADNESS")
            madness::finalize();
    #endif

    if(communicator == "MPI")
        MPI_Finalize();
#endif

    // This needs to be changed to a return value from the processed script
    return EXIT_SUCCESS;
}
