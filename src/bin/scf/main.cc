#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>

#include <libmints/basisset.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/factory.h>
#include <libmints/symmetry.h>
#include "hfenergy.h"

using namespace psi;

extern "C" {
    char *gprgid();
}

FILE *infile = NULL,
     *outfile = NULL;
char *psi_file_prefix = NULL;

char *gprgid()
{
    const char *prgid = "SCF";
    return const_cast<char*>(prgid);
}

std::string to_string(const int val);   // In matrix.cpp
// {
//     std::stringstream strm;
//     strm <<  val;
//     return strm.str();
// }

int main (int argc, char * argv[]) 
{
    psi_start(&infile,&outfile,&psi_file_prefix,argc-1,argv+1,0);
    tstart(outfile);
    
    Wavefunction::initialize_singletons();
    
    // psio_init();
    // psio_ipv1_config();
    PSIO psio;
    psiopp_ipv1_config(&psio);
    
    // chkpt_init(PSIO_OPEN_OLD);
    Chkpt chkpt(psio, PSIO_OPEN_OLD);
    
    // Initialize the psi3 timer library.
    timer_init();

    // Compute the Hartree-Fock energy
    HFEnergy hf(psio, chkpt);
    double hf_energy = hf.compute_energy();
    
    // Shut down psi. 
    timer_done();
    tstop(outfile);
    psi_stop(infile, outfile, psi_file_prefix);
    
    return EXIT_SUCCESS;
}
