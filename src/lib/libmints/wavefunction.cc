#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <liboptions/liboptions.h>
#include <libparallel/parallel.h>
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libparallel/parallel.h>

#include "mints.h"

#include <psi4-dec.h>

#include <boost/python.hpp>
#include <boost/python/call.hpp>

using namespace boost;
using namespace psi;

namespace psi {
extern FILE *infile;
}

// Globals
size_t ioff[MAX_IOFF];
double df[MAX_DF];
double bc[MAX_BC][MAX_BC];
double fac[MAX_FAC];

Wavefunction::Wavefunction(Options & options, boost::shared_ptr<PSIO> psio) :
    options_(options), psio_(psio)
{
    chkpt_ = boost::shared_ptr<Chkpt>(new Chkpt(psio.get(), PSIO_OPEN_OLD));
    common_init();
}

Wavefunction::Wavefunction(Options & options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt) :
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

    // Take the molecule from the environment
    molecule_ = Process::environment.molecule();

    // Load in the basis set
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    basisset_ = BasisSet::construct(parser, molecule_, "BASIS");

    // Check the point group of the molecule. If it is not set, set it.
    if (!molecule_->point_group()) {
        molecule_->set_point_group(molecule_->find_point_group());
    }

    // Create an SO basis...we need the point group for this part.
    integral_ = boost::shared_ptr<IntegralFactory>(new IntegralFactory(basisset_, basisset_, basisset_, basisset_));
    sobasisset_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(basisset_, integral_));

    // Obtain the dimension object to initialize the factory.
    const Dimension dimension = sobasisset_->dimension();
    factory_ = boost::shared_ptr<MatrixFactory>(new MatrixFactory);
    factory_->init_with(dimension, dimension);

    nirrep_ = dimension.n();

    // Initialize array that hold dimensionality information
    nsopi_    = Dimension(nirrep_, "SOs per irrep");
    nmopi_    = Dimension(nirrep_, "MOs per irrep");
    nalphapi_ = Dimension(nirrep_, "Alpha electrons per irrep");
    nbetapi_  = Dimension(nirrep_, "Beta electrons per irrep");
    doccpi_   = Dimension(nirrep_, "Doubly occupied orbitals per irrep");
    soccpi_   = Dimension(nirrep_, "Singly occupied orbitals per irrep");
    frzcpi_   = Dimension(nirrep_, "Frozen core orbitals per irrep");
    frzvpi_   = Dimension(nirrep_, "Frozen virtual orbitals per irrep");

    // Obtain memory amount from the environment
    memory_ = Process::environment.get_memory();

    nso_ = basisset_->nbf();
    nmo_ = basisset_->nbf();
    for (int k = 0; k < nirrep_; k++) {
        nsopi_[k] = dimension[k];
        nmopi_[k] = 0;
        doccpi_[k] = 0;
        soccpi_[k] = 0;
        nalphapi_[k] = 0;
        nbetapi_[k] = 0;
    }

    energy_ = 0.0;

    // Read in the debug flag
    debug_ = options_.get_int("DEBUG");

    // Read in energy convergence threshold
    int thresh = options_.get_int("E_CONVERGE");
    energy_threshold_ = pow(10.0, (double)-thresh);

    // Read in density convergence threshold
    thresh = options_.get_int("D_CONVERGE");
    density_threshold_ = pow(10.0, (double)-thresh);

    density_fitted_ = false;
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

boost::shared_ptr<Molecule> Wavefunction::molecule() const
{
    return molecule_;
}

boost::shared_ptr<PSIO> Wavefunction::psio() const
{
    return psio_;
}

Options& Wavefunction::options() const
{
    return options_;
}

boost::shared_ptr<IntegralFactory> Wavefunction::integral() const
{
    return integral_;
}

boost::shared_ptr<BasisSet> Wavefunction::basisset() const
{
    return basisset_;
}

boost::shared_ptr<SOBasisSet> Wavefunction::sobasisset() const
{
    return sobasisset_;
}

boost::shared_ptr<MatrixFactory> Wavefunction::matrix_factory() const
{
    return factory_;
}

boost::shared_ptr<Wavefunction> Wavefunction::reference_wavefunction() const
{
    return reference_wavefunction_;
}

void Wavefunction::set_reference_wavefunction(const boost::shared_ptr<Wavefunction> wfn)
{
    reference_wavefunction_ = wfn;
}

void Wavefunction::add_preiteration_callback(PyObject *pyobject)
{
    if (pyobject != Py_None)
        precallbacks_.push_back(pyobject);
}

void Wavefunction::add_postiteration_callback(PyObject *pyobject)
{
    if (pyobject != Py_None)
        postcallbacks_.push_back(pyobject);
}

void Wavefunction::call_preiteration_callbacks()
{
    std::vector<void*>::const_iterator iter;
    for (iter = precallbacks_.begin(); iter != precallbacks_.end(); ++iter) {
        if ((PyObject*)*iter != Py_None)
            boost::python::call<void>((PyObject*)*iter, Process::environment.reference_wavefunction());
    }
}
void Wavefunction::call_postiteration_callbacks()
{
    std::vector<void*>::const_iterator iter;
    for (iter = postcallbacks_.begin(); iter != postcallbacks_.end(); ++iter) {
        if ((PyObject*)*iter != Py_None)
            boost::python::call<void>((PyObject*)*iter, Process::environment.reference_wavefunction());
    }
}

SharedMatrix Wavefunction::Ca() const {
    if (!Ca_)
        if (!reference_wavefunction_)
            throw PSIEXCEPTION("Wavefunction::Ca: Unable to obtain MO coefficients.");
        else
            return reference_wavefunction_->Ca();

    return Ca_;
}

SharedMatrix Wavefunction::Cb() const {
    if (!Cb_)
        if (!reference_wavefunction_)
            throw PSIEXCEPTION("Wavefunction::Cb: Unable to obtain MO coefficients.");
        else
            return reference_wavefunction_->Cb();

    return Cb_;
}

SharedMatrix Wavefunction::Fa() const
{
    return Fa_;
}

SharedMatrix Wavefunction::Fb() const
{
    return Fb_;
}

boost::shared_ptr<Matrix> Wavefunction::Lagrangian() const
{
    return Lagrangian_;
}

boost::shared_ptr<Matrix> Wavefunction::tpdm_gradient_contribution() const
{
    throw PSIEXCEPTION("This type of wavefunction has not defined a TPDM gradient contribution!");
}

boost::shared_ptr<Vector> Wavefunction::epsilon_a() const
{
    return epsilon_a_;
}

boost::shared_ptr<Vector> Wavefunction::epsilon_b() const
{
    return epsilon_b_;
}

const SharedMatrix Wavefunction::Da() const
{
    return Da_;
}

SharedMatrix Wavefunction::Db() const
{
    return Db_;
}

SharedMatrix Wavefunction::X() const
{
    return Lagrangian_;
}

SharedMatrix Wavefunction::gradient() const
{
    return gradient_;
}

void Wavefunction::set_gradient(SharedMatrix& grad)
{
    gradient_ = grad;
}

void Wavefunction::save() const
{
}
