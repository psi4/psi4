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
        throw PSIEXCEPTION("Wavefunction::common_init: You must set NO_INPUT = TRUE.");
    }
    else {
        // Take the molecule from the environment
        molecule_ = Process::environment.molecule();

        // Check the point group of the molecule. If it is not set, set it.
        if (!molecule_->point_group()) {
            shared_ptr<PointGroup> pg = molecule_->find_point_group();
            if (!molecule_->symmetry_from_input().empty()) {
                // If they're not the same
                if (molecule_->symmetry_from_input() != pg->symbol()) {
                    shared_ptr<PointGroup> user(new PointGroup(molecule_->symmetry_from_input().c_str()));

                    // Make sure user is subgroup of pg
                    CorrelationTable corrtable(pg, user);

                    // If we make it here user is good.
                    pg = user;
                }
            }
            molecule_->set_point_group(pg);
        }

        // Load in the basis set
        shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options_.get_str("BASIS_PATH")));
        basisset_ = BasisSet::construct(parser, molecule_, options_.get_str("BASIS"));

        // Create an SO basis...we need the point group for this part.
        shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset_, basisset_, basisset_, basisset_));
        sobasisset_ = shared_ptr<SOBasisSet>(new SOBasisSet(basisset_, integral));

        // Obtain the dimension object to initialize the factory.
        const Dimension dimension = sobasisset_->dimension();
        factory_.init_with(dimension, dimension);

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

shared_ptr<Molecule> Wavefunction::molecule() const
{
    return molecule_    ;
}

boost::shared_ptr<BasisSet> Wavefunction::basisset() const
{
    return basisset_;
}

boost::shared_ptr<SOBasisSet> Wavefunction::sobasisset() const
{
    return sobasisset_;
}
