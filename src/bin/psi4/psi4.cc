/*! \file psi4.cc
\defgroup PSI4 The new PSI4 driver.
*/

//  PSI4 driver
//  Justin Turney
//  Rollin King

#ifdef HAVE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#include <getopt.h>
#include <stdio.h>
#include <psiconfig.h>
#include <libciomr/libciomr.h>
#include <liboptions/liboptions.h>
#include <libparallel/parallel.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libqt/qt.h>
#include "script.h"
#include <physconst.h>
#include <psifiles.h>

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

#ifdef _OPENMP
    // Disable nested threads in OpenMP
    omp_set_nested(0);
#endif

    // If no OMP thread variable is set, set nthreads to default to 1
    if (Process::environment("OMP_NUM_THREADS") == "")
        Process::environment.set_n_threads(1);

    std::string communicator = Process::environment("COMMUNICATOR");

#ifdef HAVE_MPI && HAVE_MADNESS
    if (communicator == "MADNESS")
        Communicator::world = SharedComm(new MadCommunicator(argc, argv));
    else if (communicator == "LOCAL")
        Communicator::world = SharedComm(new LocalCommunicator(argc, argv));
    else
        throw FeatureNotImplemented(communicator, "Communicator", __FILE__, __LINE__);
#else
    if (communicator == "MADNESS")
        throw FeatureNotImplemented("PSI4 was not built with ", communicator, __FILE__, __LINE__);
    else if (communicator == "LOCAL")
        Communicator::world = SharedComm(new LocalCommunicator(argc, argv));
    else
        throw FeatureNotImplemented(communicator, "Communicator", __FILE__, __LINE__);
#endif

    // There is only one timer:
    timer_init();

    // There should only be one of these in Psi4
    Wavefunction::initialize_singletons();

    // Create the scripting object
    Script::language = boost::shared_ptr<Script>(new Python);
    // Create base objects in the scripting language and initialize the language
    Script::language->initialize();

    if(psi_start(argc, argv) == PSI_RETURN_FAILURE) return EXIT_FAILURE;

    if(!clean_only) print_version(outfile);

    // Set the default memory limit for Psi4
    set_memory(outfile);

    // Initialize the I/O library
    psio_init();

    // If the user ran 'psi4 -w' catch it and exit
    if(clean_only) {
        psiclean();
        exit(EXIT_SUCCESS);
    }

    // Finish linking Psi4 into Python, preprocess the input file,
    // and then run the input file.
    Script::language->run(infile);

    // Automatically clean scratch, unless the user asked for a messy run
    if(!messy) PSIOManager::shared_object()->psiclean();

    // Shut things down:
    Communicator::world->sync();
    // There is only one timer:
    timer_done();

    psi_stop(infile, outfile, psi_file_prefix);
    Script::language->finalize();

    Communicator::world->sync();
    Communicator::world->finalize();

    // This needs to be changed to a return value from the processed script
    return EXIT_SUCCESS;
}
