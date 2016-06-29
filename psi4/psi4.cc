/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file psi4.cc
\defgroup PSI4 The new PSI4 driver.
*/

//  PSI4 driver
//  Justin Turney
//  Rollin King

#include <getopt.h>
#include <stdio.h>

#include "psi4/src/lib/libciomr/libciomr.h"
#include "psi4/src/lib/liboptions/liboptions.h"
#include "psi4/src/lib/libparallel/parallel.h"
#include "psi4/src/lib/libpsio/psio.h"
#include "psi4/src/lib/libpsio/psio.hpp"
#include "psi4/src/lib/libmints/wavefunction.h"
#include "psi4/src/lib/libqt/qt.h"
#include "script.h"
#include "psi4/include/physconst.h"
#include "psi4/include/psifiles.h"

#include <psi4-def.h>

#ifdef HAVE_AMBIT
#include <ambit/tensor.h>
#endif

#define MAIN

#include "psi4.h"
#include "script.h"

#ifdef _OPENMP
#include <omp.h>
#endif
namespace psi {
int psi_start(int argc, char *argv[]);
int psi_stop(FILE *infile, std::string, char *psi_file_prefix);
void print_version(std::string OutFileRMR);
void set_memory(std::string OutFileRMR);
int psi4_driver();
void psiclean(void);

extern int read_options(const std::string& name, Options& options, bool suppress_printing = false);

PSIO *psio = NULL;

}

#if !defined(MAKE_PYTHON_MODULE)
// This is the ONLY main function in PSI
int main(int argc, char **argv)
{
    using namespace psi;

    // Initialize external Ambit library
#ifdef HAVE_AMBIT
    ambit::initialize(argc, argv);
#endif

    // Setup the environment
    Process::environment.initialize();   // grabs the environment from the global environ variable

    //The next five lines used to live in WorldComm, they are here now
#ifdef _OPENMP
    omp_set_nested(0);
#endif
    if (Process::environment("OMP_NUM_THREADS") == "")
        Process::environment.set_n_threads(1);

    // There is only one timer:
    timer_init();

    // There should only be one of these in Psi4
    Wavefunction::initialize_singletons();

    // Create the scripting object
    Script::language = boost::shared_ptr<Script>(new Python);
    // Create base objects in the scripting language and initialize the language
    Script::language->initialize();

    if (psi_start(argc, argv) == PSI_RETURN_FAILURE) return EXIT_FAILURE;

    if (!clean_only) print_version("outfile");

    // Set the default memory limit for Psi4
    set_memory("outfile");

    // Initialize the I/O library
    psio_init();

    // If the user ran 'psi4 -w' catch it and exit
    if (clean_only) {
        psiclean();
        exit(EXIT_SUCCESS);
    }

    // Setup global options
    Process::environment.options.set_read_globals(true);
    read_options("", Process::environment.options, true);
    Process::environment.options.set_read_globals(false);

    // Finish linking Psi4 into Python, preprocess the input file,
    // and then run the input file.
    Script::language->run(infile);

    // Automatically clean scratch, unless the user asked for a messy run
    if (!messy) PSIOManager::shared_object()->psiclean();

    // Shut things down:
    // There is only one timer:
    timer_done();

    psi_stop(infile, "outfile", psi_file_prefix);
    Script::language->finalize();


    Process::environment.legacy_wavefunction().reset();

#ifdef HAVE_AMBIT
    ambit::finalize();
#endif

    // This needs to be changed to a return value from the processed script
    return EXIT_SUCCESS;
}
#endif