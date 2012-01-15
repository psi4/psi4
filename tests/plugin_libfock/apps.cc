#include <libmints/mints.h>
#include <libqt/qt.h>
#include <libpsio/psio.hpp>
#include <psi4-dec.h>
#include <physconst.h>
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
    reference_wavefunction_ = Process::environment.reference_wavefunction();
    
    if (!reference_wavefunction_) {
        throw PSIEXCEPTION("RBase: Run SCF first");
    }

    if (!reference_wavefunction_->same_a_b_dens()) {
        throw PSIEXCEPTION("RBase: Reference is not restricted");
    }
    
    copy(reference_wavefunction_);
    
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");
    Eref_ = reference_wavefunction_->reference_energy();

    Cfocc_ = Ca_subset("SO","FROZEN_OCC");
    Caocc_ = Ca_subset("SO","ACTIVE_OCC");
    Cavir_ = Ca_subset("SO","ACTIVE_VIR");
    Cfvir_ = Ca_subset("SO","FROZEN_VIR");
    eps_focc_ = epsilon_a_subset("SO","FROZEN_OCC");
    eps_aocc_ = epsilon_a_subset("SO","ACTIVE_OCC");
    eps_avir_ = epsilon_a_subset("SO","ACTIVE_VIR");
    eps_fvir_ = epsilon_a_subset("SO","FROZEN_VIR");

    std::vector<SharedMatrix> Cs;
    Cs.push_back(Cfocc_);
    Cs.push_back(Caocc_);
    Cs.push_back(Cavir_);
    Cs.push_back(Cfvir_);
    C_ = Matrix::horzcat(Cs);

    if (debug_) {
        AO2SO_->print();
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

RCPHF::RCPHF()
{
}
RCPHF::~RCPHF()
{
}
void RCPHF::print_header()
{
    fprintf(outfile, "\n");
    fprintf(outfile, "         ------------------------------------------------------------\n");
    fprintf(outfile, "                                     CPHF                           \n"); 
    fprintf(outfile, "                                  Rob Parrish                       \n"); 
    fprintf(outfile, "         ------------------------------------------------------------\n\n");

    fprintf(outfile, "  ==> Geometry <==\n\n");
    molecule_->print();
    fprintf(outfile, "  Nuclear repulsion = %20.15f\n", basisset_->molecule()->nuclear_repulsion_energy());
    fprintf(outfile, "  Reference energy  = %20.15f\n\n", Eref_);

    fprintf(outfile, "  ==> Basis Set <==\n\n");
    basisset_->print_by_level(outfile, print_);

    if (tasks_.size()) {
        fprintf(outfile, "  ==> Named Tasks <==\n\n");
        for (std::set<std::string>::const_iterator it = tasks_.begin(); it != tasks_.end(); ++it) {
            fprintf(outfile, "    %s\n", (*it).c_str());
        }
        fprintf(outfile, "\n");
    }

    if (debug_ > 1) {
        fprintf(outfile, "  ==> Fock Matrix (MO Basis) <==\n\n");
        eps_aocc_->print();
        eps_avir_->print();
    }
}
void RCPHF::add_task(const std::string& task)
{
    tasks_.insert(task);
}
void RCPHF::add_named_tasks()
{
    if (tasks_.count("POLARIZABILITY")) {
        add_polarizability();
    }
}
void RCPHF::analyze_named_tasks()
{
    if (tasks_.count("POLARIZABILITY")) {
        analyze_polarizability();
    }
}
void RCPHF::add_polarizability()
{
    OperatorSymmetry msymm(1, molecule_, integral_, factory_);
    std::vector<SharedMatrix> dipole = msymm.create_matrices("SO Dipole");
    boost::shared_ptr<OneBodySOInt> ints(integral_->so_dipole());
    ints->compute(dipole);
    
    for (int i = 0; i < dipole.size(); i++) {
        std::stringstream ss;
        ss << "Dipole Perturbation " << (i == 0 ? "X" : (i == 1 ? "Y" : "Z"));
        SharedMatrix B(new Matrix(ss.str(), Caocc_->colspi(), Cavir_->colspi(), dipole[i]->symmetry()));
        
        int symm = dipole[i]->symmetry();
        double* temp = new double[dipole[i]->max_nrow() * Cavir_->max_ncol()];

        for (int h = 0; h < B->nirrep(); h++) {
            int nsol = dipole[i]->rowspi()[h];
            int nsor = dipole[i]->colspi()[h^symm];
            int noccl = Caocc_->colspi()[h];
            int nvirr = Cavir_->colspi()[h^symm];

            if (!nsol || !nsor || !noccl || !nvirr) continue;

            double** dp = dipole[i]->pointer(h);
            double** bp = B->pointer(h);
            double** Clp = Caocc_->pointer(h);
            double** Crp = Cavir_->pointer(h^symm);
            
            C_DGEMM('N','N',nsol,nvirr,nsor,1.0,dp[0],nsor,Crp[0],nvirr,0.0,temp,nvirr);
            C_DGEMM('T','N',noccl,nvirr,nsol,1.0,Clp[0],noccl,temp,nvirr,0.0,bp[0],nvirr);
        }

        delete[] temp;
    
        b_.push_back(B); 
    }    
}
void RCPHF::analyze_polarizability()
{
    std::vector<SharedMatrix> u;
    std::vector<SharedMatrix> d;

    d.push_back(b_[b_.size() - 3]);
    d.push_back(b_[b_.size() - 2]);
    d.push_back(b_[b_.size() - 1]);

    u.push_back(x_[x_.size() - 3]);
    u.push_back(x_[x_.size() - 2]);
    u.push_back(x_[x_.size() - 1]);

    x_.pop_back();
    x_.pop_back();
    x_.pop_back();

    b_.pop_back();
    b_.pop_back();
    b_.pop_back();

    // Analysis
    SharedMatrix polarizability(new Matrix("CPHF Polarizability", 3, 3));
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            polarizability->set(0,i,j,-4.0 * (d[i]->symmetry() == u[j]->symmetry() ? d[i]->vector_dot(u[j]) : 0.0));
        }
    }

    polarizability->print();
}
double RCPHF::compute_energy()
{
    // Main CPHF Header
    print_header(); 

    // Add named tasks to the force vector list
    add_named_tasks();
    
    // Construct components
    boost::shared_ptr<JK> jk = JK::build_JK(options_, false);
    boost::shared_ptr<CPHFRHamiltonian> H(new CPHFRHamiltonian(jk, Caocc_,Cavir_,eps_aocc_,eps_avir_));
    boost::shared_ptr<CGRSolver> solver = CGRSolver::build_solver(options_,H);

    // Extra Knobs
    H->set_print(print_);
    H->set_debug(debug_);

    // Addition of force vectors
    std::vector<SharedVector>& bref = solver->b();
    std::vector<SharedVector> b = H->pack(b_);
    for (int i = 0; i < b.size(); i++) {
        bref.push_back(b[i]);
    }

    // Initialization/Memory
    solver->initialize();

    unsigned long int solver_memory = solver->memory_estimate();
    unsigned long int remaining_memory = memory_ / 8L - solver_memory;
    unsigned long int effective_memory = (unsigned long int)(options_.get_double("CPHF_MEM_SAFETY_FACTOR") * remaining_memory);
    jk->set_memory(effective_memory);
    jk->initialize();

    // Component Headers
    solver->print_header();
    H->print_header();
    jk->print_header();

    if (print_) {
        fprintf(outfile, "  ==> CPHF Iterations <==\n\n");
    }

    if (options_.get_bool("EXPLICIT_HAMILTONIAN")) {
        SharedMatrix A = H->explicit_hamiltonian();
        A->print();
    }

    if (debug_) {
        for (int Q = 0; Q < b_.size(); Q++) {
            b_[Q]->print();
        }
    }

    solver->solve();

    std::vector<SharedMatrix> x1 = H->unpack(solver->x());
    for (int Q = 0; Q < x1.size(); Q++) {
        x_.push_back(x1[Q]); 
    }

    if (debug_) {
        for (int Q = 0; Q < x_.size(); Q++) {
            x_[Q]->print();
        }
    }

    analyze_named_tasks();

    // Finalize solver/JK memory
    solver->finalize();
    jk->finalize();

    return 0.0; 
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

    fprintf(outfile, "  ==> Geometry <==\n\n");
    molecule_->print();
    fprintf(outfile, "  Nuclear repulsion = %20.15f\n", basisset_->molecule()->nuclear_repulsion_energy());
    fprintf(outfile, "  Reference energy  = %20.15f\n\n", Eref_);

    fprintf(outfile, "  ==> Basis Set <==\n\n");
    basisset_->print_by_level(outfile, print_);

    if (debug_ > 1) {
        fprintf(outfile, "  ==> Fock Matrix (MO Basis) <==\n\n");
        eps_aocc_->print();
        eps_avir_->print();
    }
}
void RCIS::sort_states()
{
    for (int n = 0; n < E_singlets_.size(); ++n) {
        states_.push_back(boost::tuple<double,int,int,int>(E_singlets_[n],n,1,singlets_[n]->symmetry()));
    }
    for (int n = 0; n < E_triplets_.size(); ++n) {
        states_.push_back(boost::tuple<double,int,int,int>(E_triplets_[n],n,3,triplets_[n]->symmetry()));
    }

    std::sort(states_.begin(), states_.end());
}
void RCIS::print_wavefunctions()
{
    fprintf(outfile, "  ==> Excitation Energies <==\n\n");

    fprintf(outfile,"  -----------------------------------------------\n");
    fprintf(outfile,"  %5s %11s %14s %14s\n",
        "State", "Description", "dE (H)", "dE (eV)");
    fprintf(outfile,"  -----------------------------------------------\n");
    char** labels = basisset_->molecule()->irrep_labels();
    for (int i = 0; i < states_.size(); i++) {
        double E = get<0>(states_[i]);
        int    j = get<1>(states_[i]);
        int    m = get<2>(states_[i]);
        int    h = get<3>(states_[i]);
        fprintf(outfile,"  %-5d %1s%-5d(%3s) %14.6E %14.6E\n",
            i + 1, (m == 1 ? "S" : "T"), j + 1, labels[h], E, _hartree2ev * E);
    }
    fprintf(outfile,"  -----------------------------------------------\n");
    fprintf(outfile, "\n");

    for(int h = 0; h < Caocc_->nirrep(); ++h) free(labels[h]); free(labels);


    if (debug_ > 1) {
        if (singlets_.size()) { 
            fprintf(outfile, "  ==> Singlet States <==\n\n");
            for (int n = 0; n < singlets_.size(); n++) {
                singlets_[n]->print();
                Dmo(singlets_[n])->print();
                Dao(singlets_[n])->print();
            }
        }

        if (triplets_.size()) {
            fprintf(outfile, "  ==> Triplet States <==\n\n");
            for (int n = 0; n < triplets_.size(); n++) {
                triplets_[n]->print();
                Dmo(triplets_[n])->print();
                Dao(triplets_[n])->print();
            }
        }
    }
}
void RCIS::print_amplitudes()
{
    if (!print_) return;

    double cutoff = options_.get_double("CIS_AMPLITUDE_CUTOFF");

    fprintf(outfile, "  ==> Significant Amplitudes <==\n\n");

    fprintf(outfile,"  --------------------------------------------------\n");
    fprintf(outfile,"  %5s %11s %20s %11s\n",
        "State", "Description", "Excitation", "Amplitude");
    fprintf(outfile,"  --------------------------------------------------\n");
    char** labels = basisset_->molecule()->irrep_labels();
    for (int i = 0; i < states_.size(); i++) {
        double E = get<0>(states_[i]);
        int    j = get<1>(states_[i]);
        int    m = get<2>(states_[i]);
        int    h = get<3>(states_[i]);

        SharedMatrix T = ((m == 1 ? singlets_[j] : triplets_[j]));
        int symm = T->symmetry();
    
        std::vector<boost::tuple<double,int,int,int,int> > amps;
        for (int h2 = 0; h2 < T->nirrep(); h2++) {
    
            int naocc = T->rowspi()[h2];
            int navir = T->colspi()[h2^symm];

            if (!naocc || !navir) continue;

            double** Tp = T->pointer(h2);

            for (int i2 = 0; i2 < naocc; i2++) {
                for (int a2 = 0; a2 < navir; a2++) {
                    if (fabs(Tp[i2][a2]) > cutoff) { 
                        int ival = i2 + Cfocc_->colspi()[h2];
                        int aval = a2 + Cfocc_->colspi()[h2^symm] + Caocc_->colspi()[h2^symm];
                        amps.push_back(boost::tuple<double,int,int,int,int>(Tp[i2][a2],ival,h2,aval,h2^symm)); 
                    } 
                }
            }
        }

        if (amps.size()) {
            std::sort(amps.begin(), amps.end());
            std::reverse(amps.begin(), amps.end());
            fprintf(outfile,"  %-5d %1s%-5d(%3s) %5d%-3s -> %5d%-3s %11.3E\n",
                i + 1, (m == 1 ? "S" : "T"), j + 1, labels[h], 
                get<1>(amps[0]) + 1, labels[get<2>(amps[0])],
                get<3>(amps[0]) + 1, labels[get<4>(amps[0])], 
                get<0>(amps[0]));
            for (int index = 1; index < amps.size(); index++) {
                fprintf(outfile,"                    %5d%-3s -> %5d%-3s %11.3E\n",
                    get<1>(amps[index]) + 1, labels[get<2>(amps[index])],
                    get<3>(amps[index]) + 1, labels[get<4>(amps[index])], 
                    get<0>(amps[index]));
            }
        } else {
            fprintf(outfile,"  %-5d %1s%-5d(%3s) %s\n",
                i + 1, (m == 1 ? "S" : "T"), j + 1, labels[h], "No Significant Amplitudes");
        }

        fprintf(outfile,"  --------------------------------------------------\n");
    }
    fprintf(outfile, "\n");
    for(int h = 0; h < Caocc_->nirrep(); ++h) free(labels[h]); free(labels);
}
void RCIS::print_transitions()
{
    if (!print_) return;

    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));
    boost::shared_ptr<OneBodyAOInt> dipole(fact->ao_dipole());

    // Get dipole integrals
    std::vector<SharedMatrix > dipole_ints;
    int nso = basisset_->nbf();
    dipole_ints.push_back(SharedMatrix(new Matrix("Dipole X", nso, nso))); 
    dipole_ints.push_back(SharedMatrix(new Matrix("Dipole Y", nso, nso))); 
    dipole_ints.push_back(SharedMatrix(new Matrix("Dipole Z", nso, nso))); 
    dipole->compute(dipole_ints);

    fprintf(outfile, "  ==> GS->XS Oscillator Strengths <==\n\n");

    fprintf(outfile,"  --------------------------------------------------------------------\n");
    fprintf(outfile,"  %5s %11s %11s %11s %11s %14s\n",
        "State", "Description", "mu_x", "mu_y", "mu_z", "f");
    fprintf(outfile,"  --------------------------------------------------------------------\n");
    char** labels = basisset_->molecule()->irrep_labels();
    for (int i = 0; i < states_.size(); i++) {

        double E = get<0>(states_[i]);
        int    j = get<1>(states_[i]);
        int    m = get<2>(states_[i]);
        int    h = get<3>(states_[i]);

        double mu[3];
        ::memset((void*) mu, '\0', 3*sizeof(double));

        // Only singlets have nonzero transition density, triplets are zero due to t_ia = - t_ia_bar
        if (m == 1) {

            SharedMatrix TD = TDao(singlets_[j]);

            // Transition dipole elements
            mu[0] = TD->vector_dot(dipole_ints[0]);
            mu[1] = TD->vector_dot(dipole_ints[1]);
            mu[2] = TD->vector_dot(dipole_ints[2]);

        }

        // Oscillator strength
        double f = 2.0 / 3.0 * E * (mu[0] * mu[0] + mu[1] * mu[1] + mu[2] * mu[2]);

        fprintf(outfile,"  %-5d %1s%-5d(%3s) %11.3E %11.3E %11.3E %14.6E\n",
            i + 1, (m == 1 ? "S" : "T"), j + 1, labels[h], mu[0],mu[1],mu[2],f);
    }
    fprintf(outfile,"  --------------------------------------------------------------------\n");
    fprintf(outfile, "\n");
    for(int h = 0; h < Caocc_->nirrep(); ++h) free(labels[h]); free(labels);
}
void RCIS::print_densities()
{
    for (int i = 0; i < options_["CIS_OPDM_STATES"].size(); i++) {
        int state = options_["CIS_OPDM_STATES"][i].to_integer();
        bool singlet = (state > 0);
        state = abs(state);

        boost::shared_ptr<Matrix> D = Dao((singlet ? singlets_[state] : triplets_[state]), false);
        std::stringstream s;
        s << (singlet ? "S" : "T") << state << "_D.dat";

        FILE* fh = fopen(s.str().c_str(),"w");
        fwrite((void*)D->pointer()[0],sizeof(double),nso_ * nso_,fh);
        fclose(fh);        
    }
    for (int i = 0; i < options_["CIS_DOPDM_STATES"].size(); i++) {
        int state = options_["CIS_DOPDM_STATES"][i].to_integer();
        bool singlet = (state > 0);
        state = abs(state);

        boost::shared_ptr<Matrix> D = Dao((singlet ? singlets_[state] : triplets_[state]), true);
        std::stringstream s;
        s << (singlet ? "S" : "T") << state << "_dD.dat";

        FILE* fh = fopen(s.str().c_str(),"w");
        fwrite((void*)D->pointer()[0],sizeof(double),nso_ * nso_,fh);
        fclose(fh);        
    }
    for (int i = 0; i < options_["CIS_TOPDM_STATES"].size(); i++) {
        int state = options_["CIS_TOPDM_STATES"][i].to_integer();
        bool singlet = (state > 0);
        state = abs(state);

        boost::shared_ptr<Matrix> D = TDao((singlet ? singlets_[state] : triplets_[state]), singlet);
        std::stringstream s;
        s << (singlet ? "S" : "T") << state << "_TD.dat";

        FILE* fh = fopen(s.str().c_str(),"w");
        fwrite((void*)D->pointer()[0],sizeof(double),nso_ * nso_,fh);
        fclose(fh);        
    }
    //for (int i = 0; i < options_["CIS_NO_STATES"].size(); i++) {
    //    int state = options_["CIS_NO_STATES"][i].to_integer();
    //    bool singlet = (state > 0);
    //    state = abs(state);

    //    std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Vector> > stuff = Nso(state,singlet);
    //    std::stringstream s;
    //    s << (singlet ? "S" : "T") << state << "_N.dat";

    //    FILE* fh = fopen(s.str().c_str(),"w");
    //    fwrite((void*)stuff.first->pointer()[0],sizeof(double),nso_ * nmo_,fh);
    //    fwrite((void*)stuff.second->pointer(),sizeof(double),nmo_,fh);
    //    fclose(fh);        
    //}
    //for (int i = 0; i < options_["CIS_AD_STATES"].size(); i++) {
    //    int state = options_["CIS_AD_STATES"][i].to_integer();
    //    bool singlet = (state > 0);
    //    state = abs(state);

    //    std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Matrix> > stuff = ADso(state,singlet);
    //    std::stringstream s;
    //    s << (singlet ? "S" : "T") << state << "_AD.dat";

    //    FILE* fh = fopen(s.str().c_str(),"w");
    //    fwrite((void*)stuff.first->pointer()[0],sizeof(double),nso_ * nso_,fh);
    //    fwrite((void*)stuff.second->pointer()[0],sizeof(double),nso_ * nso_,fh);
    //    fclose(fh);        
    //}
}
SharedMatrix RCIS::TDmo(SharedMatrix T1, bool singlet)
{
    SharedMatrix TD(T1->clone());
        
    TD->scale((singlet ? sqrt(2.0) : 0.0));
    TD->set_name("TDmo"); 

    return T1; 
}
SharedMatrix RCIS::TDso(SharedMatrix T1, bool singlet)
{
    SharedMatrix D(new Matrix("TDso", T1->nirrep(), C_->rowspi(), C_->rowspi(), T1->symmetry()));

    // Triplets are zero
    if (!singlet) return D;
    
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
        C_DGEMM('N','N',nsoocc,nsovir,nocc,sqrt(2.0),Cop[0],nocc,temp,nsovir,0.0,Dp[0],nsovir);       
 
    }

    delete[] temp;

    return D;
}
SharedMatrix RCIS::TDao(SharedMatrix T1, bool singlet)
{
    SharedMatrix D = TDso(T1, singlet);

    SharedMatrix D2(new Matrix("TDao", AO2SO_->rowspi()[0], AO2SO_->rowspi()[0]));

    double* temp = new double[AO2SO_->max_nrow() * AO2SO_->max_ncol()];

    int symm = D->symmetry();
    for (int h = 0; h < D->nirrep(); h++) {

        int nsol = AO2SO_->colspi()[h];
        int nsor = AO2SO_->colspi()[h^symm];
        int nao = AO2SO_->rowspi()[h];

        if (!nao || !nsol || !nsor) continue;

        double** Dp = D->pointer(h);
        double** Ulp = AO2SO_->pointer(h);
        double** Urp = AO2SO_->pointer(h^symm);
        double** D2p = D2->pointer();

        C_DGEMM('N','N',nao,nsor,nsol,1.0,Ulp[0],nsol,Dp[0],nsor,0.0,temp,nsor); 
        C_DGEMM('N','T',nao,nao,nsor,1.0,temp,nsor,Urp[0],nsor,1.0,D2p[0],nao); 
    }

    delete[] temp;

    return D2;
}
SharedMatrix RCIS::Dmo(SharedMatrix T1, bool diff)
{
    SharedMatrix D(new Matrix("Dmo", T1->nirrep(), reference_wavefunction_->nmopi(), reference_wavefunction_->nmopi()));

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
SharedMatrix RCIS::Dso(SharedMatrix T1, bool diff)
{
    SharedMatrix D = Dmo(T1,diff);
    SharedMatrix D2(new Matrix("Dso", C_->nirrep(), C_->rowspi(), C_->rowspi()));

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
SharedMatrix RCIS::Dao(SharedMatrix T1, bool diff)
{
    SharedMatrix D = Dso(T1,diff);
    SharedMatrix D2(new Matrix("Dao", AO2SO_->rowspi()[0], AO2SO_->rowspi()[0]));

    double* temp = new double[AO2SO_->max_nrow() * AO2SO_->max_ncol()];

    for (int h = 0; h < D->nirrep(); h++) {

        int nso = AO2SO_->colspi()[h];
        int nao = AO2SO_->rowspi()[h];

        if (!nao || !nso) continue;

        double** Dp = D->pointer(h);
        double** Up = AO2SO_->pointer(h);
        double** D2p = D2->pointer();

        C_DGEMM('N','N',nao,nso,nso,1.0,Up[0],nso,Dp[0],nso,0.0,temp,nso); 
        C_DGEMM('N','T',nao,nao,nso,1.0,temp,nso,Up[0],nso,1.0,D2p[0],nao); 
    }

    delete[] temp;

    return D2;
}
std::pair<SharedMatrix, boost::shared_ptr<Vector> > RCIS::Nmo(SharedMatrix T1, bool diff)
{
    SharedMatrix D = Dmo(T1, diff);
    SharedMatrix C(new Matrix("Nmo", D->nirrep(), D->rowspi(), D->rowspi()));
    boost::shared_ptr<Vector> O(new Vector("Occupation", D->nirrep(), D->rowspi()));
    
    D->diagonalize(C,O,Matrix::Descending);

    return make_pair(C,O);
}
std::pair<SharedMatrix, boost::shared_ptr<Vector> > RCIS::Nso(SharedMatrix T1, bool diff)
{
    std::pair<SharedMatrix, boost::shared_ptr<Vector> > pair = Nmo(T1,diff);
    SharedMatrix N = pair.first;
    boost::shared_ptr<Vector> O = pair.second;

    SharedMatrix N2(new Matrix("Nso", C_->nirrep(), C_->rowspi(), C_->colspi()));

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
std::pair<SharedMatrix, boost::shared_ptr<Vector> > RCIS::Nao(SharedMatrix T1, bool diff)
{
    std::pair<SharedMatrix, boost::shared_ptr<Vector> > pair = Nso(T1,diff);
    SharedMatrix N = pair.first;
    boost::shared_ptr<Vector> O = pair.second;

    SharedMatrix N2(new Matrix("Nso", C_->nrow(), C_->ncol()));
    SharedMatrix N3(new Matrix("Nso", C_->nrow(), C_->ncol()));
    boost::shared_ptr<Vector> O2(new Vector("Occupation", C_->ncol()));

    int offset = 0;
    std::vector<std::pair<double,int> > index;
    for (int h = 0; h < C_->nirrep(); h++) {

        int ncol = C_->ncol();
        int nmo = C_->colspi()[h];
        int nso = AO2SO_->colspi()[h];
        int nao = AO2SO_->rowspi()[h];

        if (!nmo || !nso || !nao) continue;

        for (int i = 0; i < nmo; i++) {
            index.push_back(make_pair(O->get(h,i),i+offset));
        }

        double** Np = N->pointer(h);
        double** Up = AO2SO_->pointer(h);
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
    // Main CIS Header
    print_header(); 
    
    // Construct components
    boost::shared_ptr<JK> jk = JK::build_JK(options_, false);
    boost::shared_ptr<CISRHamiltonian> H(new CISRHamiltonian(jk, Caocc_,Cavir_,eps_aocc_,eps_avir_));
    boost::shared_ptr<DLRSolver> solver = DLRSolver::build_solver(options_,H);

    // Extra Knobs
    H->set_print(print_);
    H->set_debug(debug_);

    // Initialization/Memory
    solver->initialize();
    unsigned long int solver_memory = solver->memory_estimate();
    unsigned long int remaining_memory = memory_ / 8L - solver_memory;
    unsigned long int effective_memory = (unsigned long int)(options_.get_double("CIS_MEM_SAFETY_FACTOR") * remaining_memory);
    jk->set_memory(effective_memory);
    jk->initialize();

    // Component Headers
    solver->print_header();
    H->print_header();
    jk->print_header();

    // Singlets
    if (options_.get_bool("DO_SINGLETS")) {

        H->set_singlet(true);

        if (print_) {
            fprintf(outfile, "  ==> Singlets <==\n\n");
        }

        if (options_.get_bool("EXPLICIT_HAMILTONIAN")) {
            SharedMatrix H1 = H->explicit_hamiltonian();
            H1->print();
            H->set_singlet(false);
            SharedMatrix H3 = H->explicit_hamiltonian();
            H3->print();
            return 0.0;
        }

        solver->solve();
        
        // Unpack 
        const std::vector<boost::shared_ptr<Vector> > singlets = solver->eigenvectors();    
        const std::vector<std::vector<double> > E_singlets = solver->eigenvalues();    

        std::vector<boost::shared_ptr<Matrix> > evec_temp;
        std::vector<std::pair<double, int> > eval_temp;

        for (int N = 0, index = 0; N < singlets.size(); ++N) {
            std::vector<SharedMatrix > t = H->unpack(singlets[N]);
            for (int h = 0; h < Caocc_->nirrep(); h++) {
                // Spurious zero eigenvalue due to not enough states
                if (N >= singlets[N]->dimpi()[h]) continue; 
                evec_temp.push_back(t[h]);
                eval_temp.push_back(make_pair(E_singlets[N][h], index));
                index++;
            }
        }

        std::sort(eval_temp.begin(), eval_temp.end());

        singlets_.clear();
        E_singlets_.clear();

        for (int i = 0; i < eval_temp.size(); i++) {
            E_singlets_.push_back(eval_temp[i].first);
            singlets_.push_back(evec_temp[eval_temp[i].second]);
        }
    }
    
    // Triplets
    if (options_.get_bool("DO_TRIPLETS")) {
        // Triplets
        solver->initialize();
        H->set_singlet(false);
        
        if (print_) {
            fprintf(outfile, "  ==> Triplets <==\n\n");
        }

        solver->solve();
        
        const std::vector<boost::shared_ptr<Vector> > triplets = solver->eigenvectors();    
        const std::vector<std::vector<double> > E_triplets = solver->eigenvalues();    

        std::vector<boost::shared_ptr<Matrix> > evec_temp;
        std::vector<std::pair<double, int> > eval_temp;

        for (int N = 0, index = 0; N < triplets.size(); ++N) {
            std::vector<SharedMatrix > t = H->unpack(triplets[N]);
            for (int h = 0; h < Caocc_->nirrep(); h++) {
                // Spurious zero eigenvalue due to not enough states
                if (N >= triplets[N]->dimpi()[h]) continue; 
                evec_temp.push_back(t[h]);
                eval_temp.push_back(make_pair(E_triplets[N][h], index));
                index++;
            }
        }

        std::sort(eval_temp.begin(), eval_temp.end());

        triplets_.clear();
        E_triplets_.clear();

        for (int i = 0; i < eval_temp.size(); i++) {
            E_triplets_.push_back(eval_temp[i].first);
            triplets_.push_back(evec_temp[eval_temp[i].second]);
        }

    }
    
    // Finalize solver/JK memory
    solver->finalize();
    jk->finalize();

    // Print wavefunctions and properties
    sort_states();
    print_wavefunctions();
    print_amplitudes();
    print_transitions();
    print_densities();

    return 0.0; 
}

RTDHF::RTDHF()
{
}
RTDHF::~RTDHF()
{
}
void RTDHF::print_header()
{
    fprintf(outfile, "\n");
    fprintf(outfile, "         ------------------------------------------------------------\n");
    fprintf(outfile, "                                      TDHF                           \n"); 
    fprintf(outfile, "                                  Rob Parrish                       \n"); 
    fprintf(outfile, "         ------------------------------------------------------------\n\n");

    fprintf(outfile, "  ==> Geometry <==\n\n");
    molecule_->print();
    fprintf(outfile, "  Nuclear repulsion = %20.15f\n", basisset_->molecule()->nuclear_repulsion_energy());
    fprintf(outfile, "  Reference energy  = %20.15f\n\n", Eref_);

    fprintf(outfile, "  ==> Basis Set <==\n\n");
    basisset_->print_by_level(outfile, print_);

    if (debug_ > 1) {
        fprintf(outfile, "  ==> Fock Matrix (MO Basis) <==\n\n");
        eps_aocc_->print();
        eps_avir_->print();
    }
}
double RTDHF::compute_energy()
{
    // Main TDHF Header
    print_header(); 
    
    // Construct components
    boost::shared_ptr<JK> jk = JK::build_JK(options_, false);
    boost::shared_ptr<TDHFRHamiltonian> H(new TDHFRHamiltonian(jk, Caocc_,Cavir_,eps_aocc_,eps_avir_));
    boost::shared_ptr<DLRXSolver> solver = DLRXSolver::build_solver(options_,H);

    // Extra Knobs
    H->set_print(print_);
    H->set_debug(debug_);

    // Initialization/Memory
    solver->initialize();
    unsigned long int solver_memory = solver->memory_estimate();
    unsigned long int remaining_memory = memory_ / 8L - solver_memory;
    unsigned long int effective_memory = (unsigned long int)(options_.get_double("TDHF_MEM_SAFETY_FACTOR") * remaining_memory);
    jk->set_memory(effective_memory);
    jk->initialize();

    // Component Headers
    solver->print_header();
    H->print_header();
    jk->print_header();

    // Singlets
    if (options_.get_bool("DO_SINGLETS")) {

        H->set_singlet(true);

        if (print_) {
            fprintf(outfile, "  ==> Singlets <==\n\n");
        }

        solver->solve();

        // TODO Unpack
        //const std::vector<boost::shared_ptr<Vector> > singlets = solver->eigenvectors();    
        //const std::vector<std::vector<double> > E_singlets = solver->eigenvalues();    
        //singlets_.clear();
        //E_singlets_.clear();
        //for (int N = 0; N < singlets.size(); ++N) {
        //    std::vector<SharedMatrix > t = H->unpack(singlets[N]);
        //    for (int h = 0; h < Caocc_->nirrep(); h++) {
        //        // Spurious zero eigenvalue due to not enough states
        //        if (N >= singlets[N]->dimpi()[h]) continue; 
        //        singlets_.push_back(t[h]);
        //        E_singlets_.push_back(E_singlets[N][h]);
        //    }
        //}
    }
    
    // Triplets
    if (options_.get_bool("DO_TRIPLETS")) {
        // Triplets
        solver->initialize();
        H->set_singlet(false);
        
        if (print_) {
            fprintf(outfile, "  ==> Triplets <==\n\n");
        }

        solver->solve();
        
        // TODO Unpack
        //const std::vector<boost::shared_ptr<Vector> > triplets = solver->eigenvectors();    
        //const std::vector<std::vector<double> > E_triplets = solver->eigenvalues();    
        //triplets_.clear();
        //E_triplets_.clear();
        //for (int N = 0; N < triplets.size(); ++N) {
        //    std::vector<SharedMatrix > t = H->unpack(triplets[N]);
        //    for (int h = 0; h < Caocc_->nirrep(); h++) {
        //        // Spurious zero eigenvalue due to not enough states
        //        if (N >= triplets[N]->dimpi()[h]) continue; 
        //        triplets_.push_back(t[h]);
        //        E_triplets_.push_back(E_triplets[N][h]);
        //    }
        //}
    }
    
    // Finalize solver/JK memory
    solver->finalize();
    jk->finalize();

    // TODO
    // Print wavefunctions and properties
    //print_wavefunctions();
    //print_amplitudes();
    //print_transitions();
    //print_densities();

    return 0.0; 
}

}
