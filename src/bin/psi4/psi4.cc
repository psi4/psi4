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
#include "script.h"
#include <physconst.h>

#include <molecular_system.h>

#include <psi4-def.h>
#define MAIN
#include "psi4.h"
#include "script.h"

namespace psi {
    int psi_start(int argc, char *argv[]);
    int psi_stop(FILE* infile, FILE* outfile, char* psi_file_prefix);
    void print_version(FILE *);
    void set_memory(FILE *outfile);
    int psi4_driver(Options & options, int argc, char *argv[]);
    void psiclean(void);

    int read_options(const std::string &name, Options & options);
    void read_atom_basis(char ** & atom_basis, int num_atoms);
    PSIO *psio = NULL;
    std::map<std::string, PsiReturnType(*)(Options &, int, char *[])> dispatch_table;
    Options options;
}

// This is the ONLY main function in PSI
int main(int argc, char **argv, char **envp)
{
    using namespace psi;

    // Setup the environment
    Process::arguments.init(argc, argv);
    Process::environment.init(envp);

#if HAVE_MPI == 1
    // Initialize MPI
    // Since we allow multiple MPICommunicators to be created and MPI_Init needs to be
    // as soon as possible, main takes care of call MPI_Init and not the communicator itself.
    MPI_Init(&argc, &argv);

    // Create the global world communicator
    Communicator::world = shared_ptr<Communicator>(new MPICommunicator(MPI_COMM_WORLD));
#else
    // Initialize local communicator
    Communicator::world = shared_ptr<Communicator>(new LocalCommunicator);
#endif

    // Create the scripting object
    Script::language = shared_ptr<Script>(new Python);
    // Create base objects in the scripting language and initialize the language
    Script::language->initialize();

    psi_start(argc, argv);

    if(Communicator::world->me()==0 && !clean_only) print_version(outfile);

    set_memory(outfile);

    psio_init();
    psio_ipv1_config();

    if(clean_only) {
        psiclean();
        exit(EXIT_SUCCESS);
    }

    // Okay, we might only need to make this function call if we're using IPV1
    if (!script) {
        psi4_driver(options, argc, argv);
    }
    else {
        Script::language->run(infile);
    }

    // Shut things down:
    psi_stop(infile, outfile, psi_file_prefix);
    Script::language->finalize();

#if HAVE_MPI == 1
    MPI_Finalize();
#endif

    // This needs to be changed to a return value from the processed script
    return EXIT_SUCCESS;
}
