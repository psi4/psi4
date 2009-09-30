#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>

#include "factory.h"
#include "wavefunction.h"

#include <psi4-dec.h>

using namespace psi;

namespace psi {
extern FILE *infile;
}

// Globals
int ioff[MAX_IOFF];
double df[MAX_DF];
double bc[MAX_BC][MAX_BC];
double fac[MAX_FAC];

Wavefunction::Wavefunction(shared_ptr<PSIO> psio) : psio_(psio)
{
    chkpt_ = shared_ptr<Chkpt>(new Chkpt(psio.get(), PSIO_OPEN_OLD));
    common_init();
}

Wavefunction::Wavefunction(shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt) : psio_(psio), chkpt_(chkpt)
{
    common_init();
}

Wavefunction::~Wavefunction()
{
}

void Wavefunction::common_init()
{
    Wavefunction::initialize_singletons();

    // Was a chkpt object sent to us?
    // if (chkpt_ == 0) {
    //     // No, create a checkpoint object
    //     chkpt_ = new Chkpt(psio_, PSIO_OPEN_OLD);
    // }
    
    // Initialize the matrix factory
    factory_.init_with_chkpt(chkpt_);

    // Initialize the basis set object
    basisset_ = shared_ptr<BasisSet>(new BasisSet(chkpt_));
    
    // Basis set object has reference to initialized molecule, grab it
    molecule_ = basisset_->molecule();
    
    // Read in the memory requirements from input
    fndcor(&(memory_), infile, outfile);
    
    // Read in the debug flag
    debug_ = 0;
    ip_data(const_cast<char*>("DEBUG"), const_cast<char*>("%d"), &(debug_), 0);
    
    // Read in energy convergence threshold
    int thresh = 8;
    ip_data(const_cast<char*>("E_CONVERGE"), const_cast<char*>("%d"), &(thresh), 0);
    energy_threshold_ = pow(10.0, (double)-thresh);
    
    // Read in density convergence threshold
    thresh = 8;
    ip_data(const_cast<char*>("D_CONVERGE"), const_cast<char*>("%d"), &(thresh), 0);
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
    for (int i=1; i<MAX_IOFF; ++i)
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

