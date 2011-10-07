#include <libmints/mints.h>
#include <libqt/qt.h>
#include <libpsio/psio.hpp>
#include <psi4-dec.h>
#include "apps.h"
#include "jk.h"
#include "hamiltonian.h"
#include "solver.h"

#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;

namespace psi {

RBase::RBase() :
    Wavefunction(Process::environment.options,_default_psio_lib_) 
{
    common_init();
}
RBase::~RBase()
{
}
void RBase::common_init()
{
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");

    reference_wavefunction_ = Process::environment.reference_wavefunction();
    
    if (!reference_wavefunction_) {
        throw PSIEXCEPTION("RBase: Run SCF first");
    }

    if (!reference_wavefunction_->restricted()) {
        throw PSIEXCEPTION("RBase: Reference is not restricted");
    }
    
    Ca_ = reference_wavefunction_->Ca();
    epsilon_a_ = reference_wavefunction_->epsilon_a();

    Dimension naoccpi(Ca_->nirrep());
    Dimension navirpi(Ca_->nirrep());

    for (int h = 0; h < Ca_->nirrep(); ++h) {   
        naoccpi[h] = reference_wavefunction_->doccpi()[h] - reference_wavefunction_->frzcpi()[h];
        navirpi[h] = reference_wavefunction_->nmopi()[h] - reference_wavefunction_->doccpi()[h] - reference_wavefunction_->frzvpi()[h];
    }

    Cfocc_ = boost::shared_ptr<Matrix>(new Matrix("Cfocc", reference_wavefunction_->nsopi(), reference_wavefunction_->frzcpi()));
    Cfvir_ = boost::shared_ptr<Matrix>(new Matrix("Cfvir", reference_wavefunction_->nsopi(), reference_wavefunction_->frzvpi()));
    Caocc_ = boost::shared_ptr<Matrix>(new Matrix("Caocc", reference_wavefunction_->nsopi(), naoccpi));
    Cavir_ = boost::shared_ptr<Matrix>(new Matrix("Cavir", reference_wavefunction_->nsopi(), navirpi));

    eps_focc_ = boost::shared_ptr<Vector>(new Vector("eps_focc", reference_wavefunction_->frzcpi()));
    eps_fvir_ = boost::shared_ptr<Vector>(new Vector("eps_fvir", reference_wavefunction_->frzvpi()));
    eps_aocc_ = boost::shared_ptr<Vector>(new Vector("eps_aocc", naoccpi));
    eps_avir_ = boost::shared_ptr<Vector>(new Vector("eps_avir", navirpi));
    
    for (int h = 0; h < Ca_->nirrep(); ++h) {
        for (int i = 0; i < reference_wavefunction_->frzcpi()[h]; ++i) {
            eps_focc_->set(h,i,epsilon_a_->get(h,i));
            C_DCOPY(Ca_->rowspi()[h],&Ca_->pointer(h)[0][i], Ca_->colspi()[h], &Cfocc_->pointer(h)[0][i], Cfocc_->colspi()[h]);
        }
        for (int i = 0; i < naoccpi[h]; ++i) {
            eps_aocc_->set(h,i,epsilon_a_->get(h,i + reference_wavefunction_->frzcpi()[h]));
            C_DCOPY(Ca_->rowspi()[h],&Ca_->pointer(h)[0][i + reference_wavefunction_->frzcpi()[h]], Ca_->colspi()[h], &Caocc_->pointer(h)[0][i], Caocc_->colspi()[h]);
        }
        for (int i = 0; i < navirpi[h]; ++i) {
            eps_avir_->set(h,i,epsilon_a_->get(h,i + reference_wavefunction_->doccpi()[h]));
            C_DCOPY(Ca_->rowspi()[h],&Ca_->pointer(h)[0][i + reference_wavefunction_->doccpi()[h]], Ca_->colspi()[h], &Cavir_->pointer(h)[0][i], Cavir_->colspi()[h]);
        }
        for (int i = 0; i < reference_wavefunction_->frzvpi()[h]; ++i) {
            eps_fvir_->set(h,i,epsilon_a_->get(h,i + reference_wavefunction_->doccpi()[h] + navirpi[h]));
            C_DCOPY(Ca_->rowspi()[h],&Ca_->pointer(h)[0][i + reference_wavefunction_->doccpi()[h] + navirpi[h]], Ca_->colspi()[h], &Cfvir_->pointer(h)[0][i], Cfvir_->colspi()[h]);
        }
    } 
    
    if (debug_) {
        Cfocc_->print();
        Caocc_->print();
        Cavir_->print();
        Cfvir_->print();
        eps_focc_->print();
        eps_aocc_->print();
        eps_avir_->print();
        eps_fvir_->print();
    }
}

RCIS::RCIS()
{
}
RCIS::~RCIS()
{
}
void RCIS::print_header()
{
    fprintf(outfile, "\n");
    fprintf(outfile, "         ------------------------------------------------------------\n");
    fprintf(outfile, "                                      CIS                           \n"); 
    fprintf(outfile, "                                  Rob Parrish                       \n"); 
    fprintf(outfile, "         ------------------------------------------------------------\n\n");

    fprintf(outfile, " ==> Geometry <==\n\n");
    molecule_->print();
    fprintf(outfile, "  Nuclear repulsion = %20.15f\n", basisset_->molecule()->nuclear_repulsion_energy());
    fprintf(outfile, "  Reference energy  = %20.15f\n\n", reference_wavefunction_->reference_energy());

    fprintf(outfile, "  ==> Primary Basis <==\n\n");

    basisset_->print_by_level(outfile, print_);

    if (debug_ > 1) {
        fprintf(outfile, "  ==> Fock Matrix (MO Basis) <==\n\n");
        eps_aocc_->print();
        eps_avir_->print();
    }
}
double RCIS::compute_energy()
{
    print_header(); 
    
    boost::shared_ptr<JK> jk = JK::build_JK(options_, false);
    boost::shared_ptr<CISRHamiltonian> H(new CISRHamiltonian(jk, Caocc_,Cavir_,eps_aocc_,eps_avir_));
    boost::shared_ptr<DLRSolver> solver = DLRSolver::build_solver(options_,H);

    // Knobs
    H->set_print(print_);
    H->set_debug(debug_);

    jk->initialize();
    solver->initialize();

    solver->print_header();
    H->print_header();
    jk->print_header();

    // Singlets
    H->set_singlet(true);
    solver->solve();

    const std::vector<boost::shared_ptr<Vector> > singlets = solver->eigenvectors();    
    const std::vector<std::vector<double> > E_singlets = solver->eigenvalues();    
    singlets_.clear();
    E_singlets_.clear();
    for (int N = 0; N < singlets_.size(); ++N) {
        singlets_.push_back(singlets[N]);
        E_singlets_.push_back(E_singlets[N][0]);
    }
    
    // Triplets
    solver->initialize();
    H->set_singlet(false);
    solver->solve();
    
    const std::vector<boost::shared_ptr<Vector> > triplets = solver->eigenvectors();    
    const std::vector<std::vector<double> > E_triplets = solver->eigenvalues();    
    triplets_.clear();
    E_triplets_.clear();
    for (int N = 0; N < triplets_.size(); ++N) {
        triplets_.push_back(triplets[N]);
        E_triplets_.push_back(E_triplets[N][0]);
    }
    
    solver->finalize();
    jk->finalize();

    return 0.0; 
}



}
