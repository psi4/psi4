#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <liboptions/liboptions.h>
#include <libparallel/parallel.h>
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>

#include "mints.h"

#include <psi4-dec.h>

using namespace psi;

namespace psi {
extern FILE *infile;
}

// Globals
size_t ioff[MAX_IOFF];
double df[MAX_DF];
double bc[MAX_BC][MAX_BC];
double fac[MAX_FAC];

Wavefunction::Wavefunction(Options & options, shared_ptr<PSIO> psio) :
    options_(options), psio_(psio)
{
    chkpt_ = shared_ptr<Chkpt>(new Chkpt(psio.get(), PSIO_OPEN_OLD));
    common_init();
}

Wavefunction::Wavefunction(Options & options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt) :
    options_(options), psio_(psio), chkpt_(chkpt)
{
    common_init();
}

Wavefunction::~Wavefunction()
{
}

void Wavefunction::common_init()
{

    Wavefunction::initialize_singletons();

    if (options_.get_bool("NO_INPUT") == false) {

        // Initialize the matrix factory
        factory_.init_with_chkpt(chkpt_);

        // Initialize the basis set object
        basisset_ = shared_ptr<BasisSet>(new BasisSet(chkpt_));
        
        // Basis set object has reference to initialized molecule, grab it
        molecule_ = basisset_->molecule();

        // Read in the memory requirements from input
        fndcor(&(memory_), infile, outfile);
    }
    else {
        // Take the molecule from the environment
        molecule_ = Process::environment.molecule();
        shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options_.get_str("BASIS_PATH")));
        basisset_ = BasisSet::construct(parser, molecule_, options_.get_str("BASIS"));

        int nbf[] = { basisset_->nbf() };
        factory_.init_with(1, nbf, nbf);

        memory_ = Process::environment.get_memory();    

        //fprintf(outfile,"  Using %ld bytes of core memory\n",memory_);
    }

    nso_ = basisset_->nbf();
    nmo_ = basisset_->nbf();
    for (int k = 0; k < 8; k++) {
        nsopi_[k] = 0; 
        nmopi_[k] = 0;
        doccpi_[k] = 0; 
        soccpi_[k] = 0; 
    }

    // Read in the debug flag
    debug_ = options_.get_int("DEBUG");

    // Read in energy convergence threshold
    int thresh = options_.get_int("E_CONVERGE");
    energy_threshold_ = pow(10.0, (double)-thresh);

    // Read in density convergence threshold
    thresh = options_.get_int("D_CONVERGE");
    density_threshold_ = pow(10.0, (double)-thresh);
}

double Wavefunction::compute_energy()
{
    return 0.0;
}

void Wavefunction::initialize_singletons()
{
    static bool done = false;

    if (done)
        return;

    ioff[0] = 0;
    for (size_t i=1; i<MAX_IOFF; ++i)
        ioff[i] = ioff[i-1] + i;

    df[0] = 1.0;
    df[1] = 1.0;
    df[2] = 1.0;
    for (int i=3; i<MAX_DF; ++i)
        df[i] = (i-1)*df[i-2];

    for (int i=0; i<MAX_BC; ++i)
        for (int j=0; j<=i; ++j)
            bc[i][j] = combinations(i, j);

    fac[0] = 1.0;
    for (int i=1; i<MAX_FAC; ++i)
        fac[i] = i*fac[i-1];

    done = true;
}

