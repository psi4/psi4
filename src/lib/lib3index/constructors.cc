#include "3index.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libmints/mints.h>

//MKL Header
#ifdef HAVE_MKL
#include <mkl.h>
#endif

//OpenMP Header
//_OPENMP is defined by the compiler if it exists
#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;
using namespace psi;

namespace psi {

ThreeIndexTensor::ThreeIndexTensor(shared_ptr<PSIO> psio, shared_ptr<BasisSet> primary) :
    psio_(psio), primary_basis_(primary)
{
    zero_ = BasisSet::zero_ao_basis_set();
    schwarz_cutoff_ = 1.0E-14;
    form_schwarz_sieve(schwarz_cutoff_);
    memory_ = 0L;
    nthread_ = 1;
    print_ = 0;
    nbf_ = primary_basis_->nbf();

    psio_->open(PSIF_3INDEX, PSIO_OPEN_OLD);

}
ThreeIndexTensor::~ThreeIndexTensor()
{
    free(schwarz_shells_);
    free(schwarz_funs_);
    free(schwarz_shells_reverse_);
    free(schwarz_funs_reverse_);
    free(schwarz_shell_vals_);
    free(schwarz_fun_vals_);
}
void ThreeIndexTensor::finalize()
{
    psio_->close(PSIF_3INDEX, 1);
}
DFTensor::DFTensor(shared_ptr<PSIO> psio, shared_ptr<BasisSet> primary, \
    shared_ptr<BasisSet> auxiliary, shared_ptr<BasisSet> poisson) :
    ThreeIndexTensor(psio, primary), auxiliary_basis_(auxiliary),
    poisson_basis_(poisson)
{
    poisson_ = true;
    common_init();
}
DFTensor::DFTensor(shared_ptr<PSIO> psio, shared_ptr<BasisSet> primary, \
    shared_ptr<BasisSet> auxiliary) :
    ThreeIndexTensor(psio, primary), auxiliary_basis_(auxiliary)
{
    poisson_ = false;
    common_init();
}
void DFTensor::common_init()
{
    ngaussian_ = auxiliary_basis_->nbf();
    if (!poisson_) {
        npoisson_ = 0;
    } else {
        npoisson_ = poisson_basis_->nbf();
    }
    naux_ = ngaussian_ + npoisson_;
    nfin_ = naux_;
}
DFTensor::~DFTensor()
{
}
int DFTensor::get_nfit()
{
    return nfin_;
}
int DFTensor::get_nraw()
{
    return naux_;
}
DFSCFTensor::DFSCFTensor(shared_ptr<PSIO> psio, shared_ptr<BasisSet> primary, \
    shared_ptr<BasisSet> auxiliary, shared_ptr<BasisSet> poisson) :
    DFTensor(psio, primary, auxiliary, poisson)
{
}
DFSCFTensor::DFSCFTensor(shared_ptr<PSIO> psio, shared_ptr<BasisSet> primary, \
    shared_ptr<BasisSet> auxiliary) :
    DFTensor(psio, primary, auxiliary)
{
}
DFSCFTensor::~DFSCFTensor()
{
}
CDTensor::CDTensor(shared_ptr<PSIO> psio, shared_ptr<BasisSet> primary, double delta) :
    ThreeIndexTensor(psio, primary), delta_(delta)
{
}
CDTensor::~CDTensor()
{
}
int CDTensor::get_nfit()
{
    return 0;
}
int CDTensor::get_nraw()
{
    return 0;
}
shared_ptr<DFTensor> DFTensor::bootstrap_DFTensor()
{
    shared_ptr<PSIO> psio(new PSIO());

    Options& options = Process::environment.options;
    shared_ptr<Molecule> molecule = Process::environment.molecule();

    if (molecule.get() == 0) {
        fprintf(outfile, "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    molecule->update_geometry();

    // Read in the basis set
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options.get_str("BASIS_PATH")));
    shared_ptr<BasisSet> basis = BasisSet::construct(parser, molecule, options.get_str("BASIS"));
    shared_ptr<BasisSet> auxiliary = BasisSet::construct(parser, molecule, options.get_str("RI_BASIS_SCF"));
    shared_ptr<BasisSet> poisson;
    bool use_poisson = false;
    if (options.get_str("POISSON_BASIS_SCF") != "") {
        use_poisson = true;
        poisson = BasisSet::construct(parser, molecule, options.get_str("POISSON_BASIS_SCF"));
    }

    shared_ptr<DFTensor> df;
    if (use_poisson)
        df = shared_ptr<DFTensor>(new DFTensor(psio, basis, auxiliary, poisson));
    else
        df = shared_ptr<DFTensor>(new DFTensor(psio, basis, auxiliary));

    return df;
}
void DFTensor::print(FILE* out)
{
    fprintf(out, "  ------------------------\n");
    fprintf(out, "  ==> DF Tensor Object <==\n");
    fprintf(out, "  ==>    Rob Parrish   <==\n");
    fprintf(out, "  ------------------------\n\n");

    fprintf(out, "  MOLECULE:\n\n");
    primary_basis_->molecule()->print();

    fprintf(out, "  PRIMARY BASIS:\n\n");
    primary_basis_->print(out);

    fprintf(out, "  AUXILIARY BASIS:\n\n");
    auxiliary_basis_->print(out);

    if (poisson_) {
        fprintf(out, "  POISSON BASIS:\n\n");
        poisson_basis_->print(out);
    }

    fflush(outfile);
}
void CDTensor::print(FILE* out)
{
    fprintf(out, "  ------------------------\n");
    fprintf(out, "  ==> CD Tensor Object <==\n");
    fprintf(out, "  ==>    Rob Parrish   <==\n");
    fprintf(out, "  ------------------------\n\n");

    fprintf(outfile, "  Delta = %8.5E\n\n", delta_);

    fprintf(out, "  MOLECULE:\n\n");
    primary_basis_->molecule()->print();

    fprintf(out, "  PRIMARY BASIS:\n\n");
    primary_basis_->print(out);

    fflush(outfile);
}

}
