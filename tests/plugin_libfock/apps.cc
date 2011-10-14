#include <libmints/mints.h>
#include <libqt/qt.h>
#include <libpsio/psio.hpp>
#include <psi4-dec.h>
#include "apps.h"
#include "jk.h"
#include "hamiltonian.h"
#include "solver.h"

#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;
using namespace boost;

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

    Eref_ = reference_wavefunction_->reference_energy();
    
    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));
    boost::shared_ptr<PetiteList> pet(new PetiteList(basisset_, integral));
    AO2USO_ = boost::shared_ptr<Matrix>(pet->aotoso());

    Ca_ = reference_wavefunction_->Ca();
    epsilon_a_ = reference_wavefunction_->epsilon_a();

    C_ = Ca_;

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
    fprintf(outfile, "  Reference energy  = %20.15f\n\n", Eref_);



    basisset_->print_by_level(outfile, print_);

    if (debug_ > 1) {
        fprintf(outfile, "  ==> Fock Matrix (MO Basis) <==\n\n");
        eps_aocc_->print();
        eps_avir_->print();
    }
}
void RCIS::print_wavefunctions()
{
    std::vector<boost::tuple<double, int, int, int> > states;
    for (int n = 0; n < E_singlets_.size(); ++n) {
        states.push_back(boost::tuple<double,int,int,int>(E_singlets_[n],n,1,singlets_[n]->symmetry()));
    }
    for (int n = 0; n < E_triplets_.size(); ++n) {
        states.push_back(boost::tuple<double,int,int,int>(E_triplets_[n],n,3,triplets_[n]->symmetry()));
    }

    std::sort(states.begin(), states.end());

    fprintf(outfile, "  ==> Excitation Energies <==\n\n");

    fprintf(outfile,"  -----------------------------------------------\n");
    fprintf(outfile,"  %5s %11s %14s %14s\n",
        "State", "Description", "dE (H)", "dE (eV)");
    fprintf(outfile,"  -----------------------------------------------\n");
    char** labels = basisset_->molecule()->irrep_labels();
    for (int i = 0; i < states.size(); i++) {
        double E = get<0>(states[i]);
        int    j = get<1>(states[i]);
        int    m = get<2>(states[i]);
        int    h = get<3>(states[i]);
        fprintf(outfile,"  %-5d %1s%-5d(%3s) %14.6E %14.6E\n",
            i + 1, (m == 1 ? "S" : "T"), j + 1, labels[h], E, 27.211396132 * E);
    }
    fprintf(outfile,"  -----------------------------------------------\n");
    fprintf(outfile, "\n");

    for(int h = 0; h < Caocc_->nirrep(); ++h) free(labels[h]); free(labels);


    print_transitions();

    if (debug_ > 1) {
        if (singlets_.size())
            fprintf(outfile, "  ==> Singlet States <==\n\n");
            for (int n = 0; n < singlets_.size(); n++) {
                singlets_[n]->print();
                Dmo(singlets_[n])->print();
                std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Vector> > N = Nmo(singlets_[n]);
                N.first->print();
                N.second->print();
            }

        if (triplets_.size())
            fprintf(outfile, "  ==> Triplet States <==\n\n");
            for (int n = 0; n < triplets_.size(); n++) {
                triplets_[n]->print();
                Dmo(triplets_[n])->print();
                std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Vector> > N = Nmo(triplets_[n]);
                N.first->print();
                N.second->print();
            }
    }
}
void RCIS::print_transitions()
{
    if (!print_) return;

    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));
    boost::shared_ptr<OneBodyAOInt> dipole(fact->ao_dipole());

    // Get dipole integrals
    std::vector<boost::shared_ptr<Matrix> > dipole_ints;
    int nso = basisset_->nbf();
    dipole_ints.push_back(boost::shared_ptr<Matrix>(new Matrix("Dipole X", nso, nso))); 
    dipole_ints.push_back(boost::shared_ptr<Matrix>(new Matrix("Dipole Y", nso, nso))); 
    dipole_ints.push_back(boost::shared_ptr<Matrix>(new Matrix("Dipole Z", nso, nso))); 
    dipole->compute(dipole_ints);

    if (debug_) {
        for (int i = 0; i < dipole_ints.size(); i++) {
            dipole_ints[i]->print();
        }
    }

    std::vector<boost::tuple<double, int, int, int> > states;
    for (int n = 0; n < E_singlets_.size(); ++n) {
        states.push_back(boost::tuple<double,int,int,int>(E_singlets_[n],n,1,singlets_[n]->symmetry()));
    }
    for (int n = 0; n < E_triplets_.size(); ++n) {
        states.push_back(boost::tuple<double,int,int,int>(E_triplets_[n],n,3,triplets_[n]->symmetry()));
    }

    std::sort(states.begin(), states.end());

    fprintf(outfile, "  ==> GS->XS Oscillator Strengths <==\n\n");

    fprintf(outfile,"  --------------------------------------------------------------------\n");
    fprintf(outfile,"  %5s %11s %11s %11s %11s %14s\n",
        "State", "Description", "mu_x", "mu_y", "mu_z", "f");
    fprintf(outfile,"  --------------------------------------------------------------------\n");
    char** labels = basisset_->molecule()->irrep_labels();
    for (int i = 0; i < states.size(); i++) {
        double E = get<0>(states[i]);
        int    j = get<1>(states[i]);
        int    m = get<2>(states[i]);
        int    h = get<3>(states[i]);

        boost::shared_ptr<Matrix> TD = TDao((m == 1 ? singlets_[j] : triplets_[j]));

        double mu[3];

        // Transition dipole elements
        mu[0] = TD->vector_dot(dipole_ints[0]);
        mu[1] = TD->vector_dot(dipole_ints[1]);
        mu[2] = TD->vector_dot(dipole_ints[2]);

        // Oscillator strength
        double f = 2.0 / 3.0 * E * (mu[0] * mu[0] + mu[1] * mu[1] + mu[2] * mu[2]);

        fprintf(outfile,"  %-5d %1s%-5d(%3s) %11.3E %11.3E %11.3E %14.6E\n",
            i + 1, (m == 1 ? "S" : "T"), j + 1, labels[h], mu[0],mu[1],mu[2],f);
    }
    fprintf(outfile,"  --------------------------------------------------------------------\n");
    fprintf(outfile, "\n");
    for(int h = 0; h < Caocc_->nirrep(); ++h) free(labels[h]); free(labels);
}
boost::shared_ptr<Matrix> RCIS::TDmo(boost::shared_ptr<Matrix> T1)
{
    return T1; // you idiot
}
boost::shared_ptr<Matrix> RCIS::TDso(boost::shared_ptr<Matrix> T1)
{
    boost::shared_ptr<Matrix> D(new Matrix("TDso", T1->nirrep(), C_->rowspi(), C_->rowspi(), T1->symmetry()));

    double* temp = new double[C_->max_nrow() * (ULI) T1->max_nrow()];
    
    int symm = T1->symmetry();
    for (int h = 0; h < T1->nirrep(); h++) {

        int nocc = T1->rowspi()[h];
        int nvir = T1->colspi()[h^symm];
        int nsoocc = Caocc_->rowspi()[h];
        int nsovir = Cavir_->rowspi()[h^symm];

        if (!nocc || !nvir || !nsoocc || !nsovir) continue;
        
        double** Dp = D->pointer(h);
        double** Tp = T1->pointer(h);
        double** Cop = Caocc_->pointer(h);
        double** Cvp = Cavir_->pointer(h^symm);

        C_DGEMM('N','T',nocc,nsovir,nvir,1.0,Tp[0],nvir,Cvp[0],nvir,0.0,temp,nsovir);
        C_DGEMM('N','N',nsoocc,nsovir,nocc,1.0,Cop[0],nocc,temp,nsovir,0.0,Dp[0],nsovir);       
 
    }

    delete[] temp;

    return D;
}
boost::shared_ptr<Matrix> RCIS::TDao(boost::shared_ptr<Matrix> T1)
{
    boost::shared_ptr<Matrix> D = TDso(T1);

    boost::shared_ptr<Matrix> D2(new Matrix("TDao", AO2USO_->nrow(), AO2USO_->nrow()));

    double* temp = new double[AO2USO_->max_nrow() * AO2USO_->max_ncol()];

    int symm = D->symmetry();
    for (int h = 0; h < D->nirrep(); h++) {

        int nsol = AO2USO_->colspi()[h];
        int nsor = AO2USO_->colspi()[h^symm];
        int nao = AO2USO_->rowspi()[h];

        if (!nao || !nsol || !nsor) continue;

        double** Dp = D->pointer(h);
        double** Ulp = AO2USO_->pointer(h);
        double** Urp = AO2USO_->pointer(h^symm);
        double** D2p = D2->pointer();

        C_DGEMM('N','N',nao,nsor,nsol,1.0,Ulp[0],nsol,Dp[0],nsor,0.0,temp,nsor); 
        C_DGEMM('N','T',nao,nao,nsor,1.0,temp,nsor,Urp[0],nsor,1.0,D2p[0],nao); 
    }

    delete[] temp;

    return D2;
    return D2;
}
boost::shared_ptr<Matrix> RCIS::Dmo(boost::shared_ptr<Matrix> T1, bool diff)
{
    boost::shared_ptr<Matrix> D(new Matrix("Dmo", T1->nirrep(), reference_wavefunction_->nmopi(), reference_wavefunction_->nmopi()));

    int symm = T1->symmetry();

    /// Reference occupation (if not difference density)
    if (!diff) {
        for (int h = 0; h < D->nirrep(); ++h) {
            int nmo = D->rowspi()[h];
            for (int i = 0; i < eps_focc_->dimpi()[h] + eps_aocc_->dimpi()[h]; i++) {
                D->set(h,i,i,1.0);
            }
        }
    }

    /// Depletion of occupied space 
    for (int h = 0; h < D->nirrep(); ++h) {
        int nmo = D->rowspi()[h];
        int naocc = T1->rowspi()[h];
        int navir = T1->colspi()[h^symm];
        int nfocc = eps_focc_->dimpi()[h];
        
        if (!nmo || !naocc || !navir) continue;

        double** Tp = T1->pointer(h);
        double** Dp = D->pointer(h);

        C_DGEMM('N','T',naocc,naocc,navir,-0.5,Tp[0],navir,Tp[0],navir,1.0,&Dp[nfocc][nfocc],nmo); 
    } 
    
    /// Accumulation of virtual space 
    for (int h = 0; h < D->nirrep(); ++h) {
        int nmo = D->rowspi()[h];
        int naocc = T1->rowspi()[h^symm];
        int navir = T1->colspi()[h];
        int nocc = eps_focc_->dimpi()[h] + eps_aocc_->dimpi()[h];
        
        if (!nmo || !naocc || !navir) continue;

        double** Tp = T1->pointer(h^symm);
        double** Dp = D->pointer(h);

        C_DGEMM('T','N',navir,navir,naocc,0.5,Tp[0],navir,Tp[0],navir,1.0,&Dp[nocc][nocc],nmo); 
    } 
    
    return D;
}
boost::shared_ptr<Matrix> RCIS::Dso(boost::shared_ptr<Matrix> T1, bool diff)
{
    boost::shared_ptr<Matrix> D = Dmo(T1,diff);
    boost::shared_ptr<Matrix> D2(new Matrix("Dso", C_->nirrep(), C_->rowspi(), C_->rowspi()));

    double* temp = new double[C_->max_nrow() * C_->max_ncol()];

    for (int h = 0; h < D->nirrep(); h++) {

        int nmo = C_->colspi()[h];
        int nso = C_->rowspi()[h];

        if (!nmo || !nso) continue;

        double** Dp = D->pointer(h);
        double** Cp = C_->pointer(h);
        double** D2p = D2->pointer(h);

        C_DGEMM('N','N',nso,nmo,nmo,1.0,Cp[0],nmo,Dp[0],nmo,0.0,temp,nmo); 
        C_DGEMM('N','T',nso,nso,nmo,1.0,temp,nmo,Cp[0],nmo,0.0,D2p[0],nso); 

    }

    delete[] temp;

    return D2;
}
boost::shared_ptr<Matrix> RCIS::Dao(boost::shared_ptr<Matrix> T1, bool diff)
{
    boost::shared_ptr<Matrix> D = Dso(T1,diff);
    boost::shared_ptr<Matrix> D2(new Matrix("Dao", AO2USO_->nrow(), AO2USO_->nrow()));

    double* temp = new double[AO2USO_->max_nrow() * AO2USO_->max_ncol()];

    for (int h = 0; h < D->nirrep(); h++) {

        int nso = AO2USO_->colspi()[h];
        int nao = AO2USO_->rowspi()[h];

        if (!nao || !nso) continue;

        double** Dp = D->pointer(h);
        double** Up = AO2USO_->pointer(h);
        double** D2p = D2->pointer();

        C_DGEMM('N','N',nao,nso,nso,1.0,Up[0],nso,Dp[0],nso,0.0,temp,nso); 
        C_DGEMM('N','T',nao,nao,nso,1.0,temp,nso,Up[0],nso,1.0,D2p[0],nao); 
    }

    delete[] temp;

    return D2;
}
std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Vector> > RCIS::Nmo(boost::shared_ptr<Matrix> T1, bool diff)
{
    boost::shared_ptr<Matrix> D = Dmo(T1, diff);
    boost::shared_ptr<Matrix> C(new Matrix("Nmo", D->nirrep(), D->rowspi(), D->rowspi()));
    boost::shared_ptr<Vector> O(new Vector("Occupation", D->nirrep(), D->rowspi()));
    
    D->diagonalize(C,O,Matrix::Descending);

    return make_pair(C,O);
}
std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Vector> > RCIS::Nso(boost::shared_ptr<Matrix> T1, bool diff)
{
    std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Vector> > pair = Nmo(T1,diff);
    boost::shared_ptr<Matrix> N = pair.first;
    boost::shared_ptr<Vector> O = pair.second;

    boost::shared_ptr<Matrix> N2(new Matrix("Nso", C_->nirrep(), C_->rowspi(), C_->colspi()));

    for (int h = 0; h < N->nirrep(); h++) {

        int nmo = C_->colspi()[h];
        int nso = C_->rowspi()[h];

        if (!nmo || !nso) continue;

        double** Np = N->pointer(h);
        double** Cp = C_->pointer(h);
        double** N2p = N2->pointer(h);

        C_DGEMM('N','N',nso,nmo,nmo,1.0,Cp[0],nmo,Np[0],nmo,0.0,N2p[0],nmo); 
    }
    return make_pair(N2,O);
}
std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Vector> > RCIS::Nao(boost::shared_ptr<Matrix> T1, bool diff)
{
    std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Vector> > pair = Nso(T1,diff);
    boost::shared_ptr<Matrix> N = pair.first;
    boost::shared_ptr<Vector> O = pair.second;

    boost::shared_ptr<Matrix> N2(new Matrix("Nso", C_->nrow(), C_->ncol()));
    boost::shared_ptr<Matrix> N3(new Matrix("Nso", C_->nrow(), C_->ncol()));
    boost::shared_ptr<Vector> O2(new Vector("Occupation", C_->ncol()));

    int offset = 0;
    std::vector<std::pair<double,int> > index;
    for (int h = 0; h < C_->nirrep(); h++) {

        int ncol = C_->ncol();
        int nmo = C_->colspi()[h];
        int nso = AO2USO_->colspi()[h];
        int nao = AO2USO_->rowspi()[h];

        if (!nmo || !nso || !nao) continue;

        for (int i = 0; i < nmo; i++) {
            index.push_back(make_pair(O->get(h,i),i+offset));
        }

        double** Np = N->pointer(h);
        double** Up = AO2USO_->pointer(h);
        double** N2p = N2->pointer(h);

        C_DGEMM('N','N',nao,nmo,nso,1.0,Up[0],nso,Np[0],nmo,0.0,&N2p[0][offset],ncol); 

        offset += nmo;
    }

    std::sort(index.begin(), index.end(), std::greater<std::pair<double,int> >());
    
    int nmo = N2->colspi()[0];
    int nao = N2->rowspi()[0]; 

    for (int i = 0; i < nmo; i++) {
        double occ = index[i].first;
        int ind    = index[i].second;
        O2->set(i,occ);

        C_DCOPY(nao, &(N2->pointer()[0][ind]), nmo, &(N3->pointer()[0][i]), nmo);
    }

    return make_pair(N3,O2);
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

    if (print_) {
        fprintf(outfile, "  ==> Singlets <==\n\n");
    }

    solver->solve();

    const std::vector<boost::shared_ptr<Vector> > singlets = solver->eigenvectors();    
    const std::vector<std::vector<double> > E_singlets = solver->eigenvalues();    
    singlets_.clear();
    E_singlets_.clear();
    for (int N = 0; N < singlets.size(); ++N) {
        std::vector<boost::shared_ptr<Matrix> > t = H->unpack(singlets[N]);
        for (int h = 0; h < Caocc_->nirrep(); h++) {
            singlets_.push_back(t[h]);
            E_singlets_.push_back(E_singlets[N][h]);
        }
    }
    
    // Triplets
    solver->initialize();
    H->set_singlet(false);

    if (print_) {
        fprintf(outfile, "  ==> Triplets <==\n\n");
    }

    solver->solve();
    
    const std::vector<boost::shared_ptr<Vector> > triplets = solver->eigenvectors();    
    const std::vector<std::vector<double> > E_triplets = solver->eigenvalues();    
    triplets_.clear();
    E_triplets_.clear();
    for (int N = 0; N < triplets.size(); ++N) {
        std::vector<boost::shared_ptr<Matrix> > t = H->unpack(triplets[N]);
        for (int h = 0; h < Caocc_->nirrep(); h++) {
            triplets_.push_back(t[h]);
            E_triplets_.push_back(E_triplets[N][h]);
        }
    }
    
    solver->finalize();
    jk->finalize();

    print_wavefunctions();

    return 0.0; 
}



}
