/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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


#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"
#include "psi4/physconst.h"
#include "apps.h"
#include "jk.h"
#include "v.h"
#include "hamiltonian.h"
#include "solver.h"
#include "psi4/libscf_solver/hf.h"
#include "psi4/libscf_solver/rhf.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/multipolesymmetry.h"

#include <algorithm>
#include <tuple>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {

RBase::RBase(SharedWavefunction ref_wfn, Options& options, bool use_symmetry) :
    Wavefunction(options),
    use_symmetry_(use_symmetry)
{
    shallow_copy(ref_wfn);

    set_reference(ref_wfn);

    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");
    bench_ = options_.get_int("BENCH");
    convergence_ = options_.get_double("SOLVER_CONVERGENCE");
}
RBase::RBase(bool flag) :
    Wavefunction(Process::environment.options)
{
    psio_ = _default_psio_lib_;
    throw PSIEXCEPTION("DGAS: Lets not let RMP do dirty hacks!\n");
    outfile->Printf( "Dirty hack %s\n\n", (flag ? "true" : "false"));
}
RBase::~RBase()
{
    postiterations();
}
void RBase::set_reference(SharedWavefunction ref_wfn)
{
    reference_wavefunction_ = ref_wfn;

    if (!reference_wavefunction_->same_a_b_dens()) {
        throw PSIEXCEPTION("RBase: Reference is not restricted");
    }

    Eref_ = reference_wavefunction_->reference_energy();

    if(use_symmetry_) {
        Cocc_  = Ca_subset("SO","OCC");
        Cfocc_ = Ca_subset("SO","FROZEN_OCC");
        Caocc_ = Ca_subset("SO","ACTIVE_OCC");
        Cavir_ = Ca_subset("SO","ACTIVE_VIR");
        Cfvir_ = Ca_subset("SO","FROZEN_VIR");
        eps_focc_ = epsilon_a_subset("SO","FROZEN_OCC");
        eps_aocc_ = epsilon_a_subset("SO","ACTIVE_OCC");
        eps_avir_ = epsilon_a_subset("SO","ACTIVE_VIR");
        eps_fvir_ = epsilon_a_subset("SO","FROZEN_VIR");
    } else {
        Cocc_  = Ca_subset("AO","OCC");
        Cfocc_ = Ca_subset("AO","FROZEN_OCC");
        Caocc_ = Ca_subset("AO","ACTIVE_OCC");
        Cavir_ = Ca_subset("AO","ACTIVE_VIR");
        Cfvir_ = Ca_subset("AO","FROZEN_VIR");
        eps_focc_ = epsilon_a_subset("AO","FROZEN_OCC");
        eps_aocc_ = epsilon_a_subset("AO","ACTIVE_OCC");
        eps_avir_ = epsilon_a_subset("AO","ACTIVE_VIR");
        eps_fvir_ = epsilon_a_subset("AO","FROZEN_VIR");
    }

    std::vector<SharedMatrix> Cs;
    Cs.push_back(Cfocc_);
    Cs.push_back(Caocc_);
    Cs.push_back(Cavir_);
    Cs.push_back(Cfvir_);
    C_ = Matrix::horzcat(Cs);
}
void RBase::preiterations()
{
    if (!jk_) {
        if (options_.get_bool("SAVE_JK")) {
            jk_ = (static_cast<psi::scf::HF*>(reference_wavefunction_.get()))->jk();
            outfile->Printf("    Reusing JK object from SCF.\n\n");
        } else {
            if (options_.get_str("SCF_TYPE") == "DF"){
                jk_ = JK::build_JK(basisset_, get_basisset("DF_BASIS_SCF"), options_);
            } else {
                jk_ = JK::build_JK(basisset_, BasisSet::zero_ao_basis_set(), options_);
            }
            unsigned long int effective_memory = (unsigned long int)(0.125 * options_.get_double("CPHF_MEM_SAFETY_FACTOR") * memory_);
            jk_->set_memory(effective_memory);
            jk_->initialize();
        }
    }

    if (!v_) {
        if (options_.get_str("MODULE") == "RCPKS" || options_.get_str("MODULE") == "RTDA" || options_.get_str("MODULE") == "RTDDFT") {
            throw PSIEXCEPTION("V is not currently enabled in apps.cc");
            // v_ = VBase::build_V(basisset_, options_, "RK");
            // v_->initialize();
        }
    }
}
void RBase::postiterations()
{
    jk_.reset();
}

RCPHF::RCPHF(SharedWavefunction ref_wfn, Options& options, bool use_symmetry) :
    RBase(ref_wfn, options, use_symmetry)
{
}
RCPHF::~RCPHF()
{
}
void RCPHF::print_header()
{
    outfile->Printf( "\n");
    outfile->Printf( "         ------------------------------------------------------------\n");
    outfile->Printf( "                                     CPHF                           \n");
    outfile->Printf( "                                  Rob Parrish                       \n");
    outfile->Printf( "         ------------------------------------------------------------\n\n");

    outfile->Printf( "  ==> Geometry <==\n\n");
    molecule_->print();
    outfile->Printf( "  Nuclear repulsion = %20.15f\n", basisset_->molecule()->nuclear_repulsion_energy());
    outfile->Printf( "  Reference energy  = %20.15f\n\n", Eref_);

    outfile->Printf( "  ==> Basis Set <==\n\n");
    basisset_->print_by_level("outfile", print_);

    if (tasks_.size()) {
        outfile->Printf( "  ==> Named Tasks <==\n\n");
        for (std::set<std::string>::const_iterator it = tasks_.begin(); it != tasks_.end(); ++it) {
            outfile->Printf( "    %s\n", (*it).c_str());
        }
        outfile->Printf( "\n");
    }

    if (debug_ > 1) {
        outfile->Printf( "  ==> Fock Matrix (MO Basis) <==\n\n");
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
    std::shared_ptr<OneBodySOInt> ints(integral_->so_dipole());
    ints->compute(dipole);

    for (size_t i = 0; i < dipole.size(); i++) {
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

        std::stringstream ss2;
        ss2 << (i == 0 ? "MU_X" : (i == 1 ? "MU_Y" : "MU_Z"));
        b_[ss2.str()] = B;
    }
}
void RCPHF::analyze_polarizability()
{
    std::vector<SharedMatrix> u;
    std::vector<SharedMatrix> d;

    d.push_back(b_["MU_X"]);
    d.push_back(b_["MU_Y"]);
    d.push_back(b_["MU_Z"]);

    u.push_back(x_["MU_X"]);
    u.push_back(x_["MU_Y"]);
    u.push_back(x_["MU_Z"]);

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

    if (!jk_)
        preiterations();

    // Construct components
    std::shared_ptr<CPHFRHamiltonian> H(new CPHFRHamiltonian(jk_, Caocc_,Cavir_,eps_aocc_,eps_avir_));
    std::shared_ptr<CGRSolver> solver = CGRSolver::build_solver(options_,H);

    // Extra Knobs
    H->set_print(print_);
    H->set_debug(debug_);
    solver->set_convergence(convergence_);

    // Addition of force vectors
    std::vector<SharedVector>& bref = solver->b();
    std::map<std::string, SharedVector> b = H->pack(b_);
    for (std::map<std::string, SharedVector>::const_iterator it = b.begin();
        it != b.end(); ++it) {
        bref.push_back((*it).second);
    }

    // Initialization/Memory
    solver->initialize();

    // Component Headers
    solver->print_header();
    H->print_header();
    jk_->print_header();

    if (print_) {
        outfile->Printf( "  ==> CPHF Iterations <==\n\n");
    }

    if (options_.get_bool("EXPLICIT_HAMILTONIAN")) {
        SharedMatrix A = H->explicit_hamiltonian();
        A->print();
    }

    if (debug_) {
        for (std::map<std::string, SharedMatrix>::const_iterator it = b_.begin();
            it != b_.end(); ++it) {
            (*it).second->print();
        }
    }

    solver->solve();

    std::vector<SharedMatrix> x1 = H->unpack(solver->x());

    int index = 0;
    for (std::map<std::string, SharedMatrix>::const_iterator it = b_.begin();
        it != b_.end(); ++it) {
        x_[(*it).first] = x1[index++];
    }

    if (debug_) {
        for (std::map<std::string, SharedMatrix>::const_iterator it = x_.begin();
            it != x_.end(); ++it) {
            (*it).second->print();
        }
    }

    analyze_named_tasks();

    solver->finalize();

    return 0.0;
}

RCIS::RCIS(SharedWavefunction ref_wfn, Options& options) :
           RBase(ref_wfn, options)
{
}
RCIS::~RCIS()
{
}
void RCIS::print_header()
{
    outfile->Printf( "\n");
    outfile->Printf( "         ------------------------------------------------------------\n");
    outfile->Printf( "                                      CIS                           \n");
    outfile->Printf( "                                  Rob Parrish                       \n");
    outfile->Printf( "         ------------------------------------------------------------\n\n");

    outfile->Printf( "  ==> Geometry <==\n\n");
    molecule_->print();
    outfile->Printf( "  Nuclear repulsion = %20.15f\n", basisset_->molecule()->nuclear_repulsion_energy());
    outfile->Printf( "  Reference energy  = %20.15f\n\n", Eref_);

    outfile->Printf( "  ==> Basis Set <==\n\n");
    basisset_->print_by_level("outfile", print_);

    if (debug_ > 1) {
        outfile->Printf( "  ==> Fock Matrix (MO Basis) <==\n\n");
        eps_aocc_->print();
        eps_avir_->print();
    }
}
void RCIS::sort_states()
{
    for (size_t n = 0; n < E_singlets_.size(); ++n) {
        states_.push_back(std::tuple<double,int,int,int>(E_singlets_[n],n,1,singlets_[n]->symmetry()));
    }
    for (size_t n = 0; n < E_triplets_.size(); ++n) {
        states_.push_back(std::tuple<double,int,int,int>(E_triplets_[n],n,3,triplets_[n]->symmetry()));
    }

    std::sort(states_.begin(), states_.end());
}
void RCIS::print_wavefunctions()
{
    outfile->Printf( "  ==> Excitation Energies <==\n\n");

    outfile->Printf("  -----------------------------------------------\n");
    outfile->Printf("  %5s %11s %14s %14s\n",
        "State", "Description", "dE (H)", "dE (eV)");
    outfile->Printf("  -----------------------------------------------\n");
    char** labels = basisset_->molecule()->irrep_labels();
    for (size_t i = 0; i < states_.size(); i++) {
        double E = std::get<0>(states_[i]);
        int    j = std::get<1>(states_[i]);
        int    m = std::get<2>(states_[i]);
        int    h = std::get<3>(states_[i]);
        outfile->Printf("  %-5d %1s%-5d(%3s) %14.6E %14.6E\n",
            i + 1, (m == 1 ? "S" : "T"), j + 1, labels[h], E, pc_hartree2ev * E);
    }
    outfile->Printf("  -----------------------------------------------\n");
    outfile->Printf( "\n");

    for(int h = 0; h < Caocc_->nirrep(); ++h) free(labels[h]); free(labels);


    if (debug_ > 1) {
        if (singlets_.size()) {
            outfile->Printf( "  ==> Singlet States <==\n\n");
            for (size_t n = 0; n < singlets_.size(); n++) {
                singlets_[n]->print();
                Dmo(singlets_[n])->print();
                Dao(singlets_[n])->print();
            }
        }

        if (triplets_.size()) {
            outfile->Printf( "  ==> Triplet States <==\n\n");
            for (size_t n = 0; n < triplets_.size(); n++) {
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

    outfile->Printf( "  ==> Significant Amplitudes <==\n\n");

    outfile->Printf("  --------------------------------------------------\n");
    outfile->Printf("  %5s %11s %20s %11s\n",
        "State", "Description", "Excitation", "Amplitude");
    outfile->Printf("  --------------------------------------------------\n");
    char** labels = basisset_->molecule()->irrep_labels();
    for (size_t i = 0; i < states_.size(); i++) {
//        double E = std::get<0>(states_[i]);
        int    j = std::get<1>(states_[i]);
        int    m = std::get<2>(states_[i]);
        int    h = std::get<3>(states_[i]);

        SharedMatrix T = ((m == 1 ? singlets_[j] : triplets_[j]));
        int symm = T->symmetry();

        std::vector<std::tuple<double,int,int,int,int> > amps;
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
                        amps.push_back(std::tuple<double,int,int,int,int>(Tp[i2][a2],ival,h2,aval,h2^symm));
                    }
                }
            }
        }

        if (amps.size()) {
            std::sort(amps.begin(), amps.end());
            std::reverse(amps.begin(), amps.end());
            outfile->Printf("  %-5d %1s%-5d(%3s) %5d%-3s -> %5d%-3s %11.3E\n",
                i + 1, (m == 1 ? "S" : "T"), j + 1, labels[h],
                std::get<1>(amps[0]) + 1, labels[std::get<2>(amps[0])],
                std::get<3>(amps[0]) + 1, labels[std::get<4>(amps[0])],
                std::get<0>(amps[0]));
            for (size_t index = 1; index < amps.size(); index++) {
                outfile->Printf("                    %5d%-3s -> %5d%-3s %11.3E\n",
                    std::get<1>(amps[index]) + 1, labels[std::get<2>(amps[index])],
                    std::get<3>(amps[index]) + 1, labels[std::get<4>(amps[index])],
                    std::get<0>(amps[index]));
            }
        } else {
            outfile->Printf("  %-5d %1s%-5d(%3s) %s\n",
                i + 1, (m == 1 ? "S" : "T"), j + 1, labels[h], "No Significant Amplitudes");
        }

        outfile->Printf("  --------------------------------------------------\n");
    }
    outfile->Printf( "\n");
    for(int h = 0; h < Caocc_->nirrep(); ++h) free(labels[h]); free(labels);
}
void RCIS::print_transitions()
{
    if (!print_) return;

    std::shared_ptr<IntegralFactory> fact(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));
    std::shared_ptr<OneBodyAOInt> dipole(fact->ao_dipole());

    // Get dipole integrals
    std::vector<SharedMatrix > dipole_ints;
    int nso = basisset_->nbf();
    dipole_ints.push_back(SharedMatrix(new Matrix("Dipole X", nso, nso)));
    dipole_ints.push_back(SharedMatrix(new Matrix("Dipole Y", nso, nso)));
    dipole_ints.push_back(SharedMatrix(new Matrix("Dipole Z", nso, nso)));
    dipole->compute(dipole_ints);

    outfile->Printf( "  ==> GS->XS Oscillator Strengths <==\n\n");

    outfile->Printf("  --------------------------------------------------------------------\n");
    outfile->Printf("  %5s %11s %11s %11s %11s %14s\n",
        "State", "Description", "mu_x", "mu_y", "mu_z", "f");
    outfile->Printf("  --------------------------------------------------------------------\n");
    char** labels = basisset_->molecule()->irrep_labels();
    for (size_t i = 0; i < states_.size(); i++) {

        double E = std::get<0>(states_[i]);
        int    j = std::get<1>(states_[i]);
        int    m = std::get<2>(states_[i]);
        int    h = std::get<3>(states_[i]);

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

        outfile->Printf("  %-5d %1s%-5d(%3s) %11.3E %11.3E %11.3E %14.6E\n",
            i + 1, (m == 1 ? "S" : "T"), j + 1, labels[h], mu[0],mu[1],mu[2],f);
    }
    outfile->Printf("  --------------------------------------------------------------------\n");
    outfile->Printf( "\n");
    for(int h = 0; h < Caocc_->nirrep(); ++h) free(labels[h]); free(labels);
}
void RCIS::print_densities()
{
    for (unsigned int i = 0; i < options_["CIS_OPDM_STATES"].size(); i++) {
        int state = options_["CIS_OPDM_STATES"][i].to_integer();
        bool singlet = (state > 0);
        state = abs(state);

        std::shared_ptr<Matrix> D = Dao((singlet ? singlets_[state-1] : triplets_[state-1]), false);
        std::stringstream s;
        s << (singlet ? "S" : "T") << state << "_D.dat";

        FILE* fh = fopen(s.str().c_str(),"w");
        fwrite((void*)D->pointer()[0],sizeof(double),nso_ * nso_,fh);
        fclose(fh);
    }
    for (unsigned int i = 0; i < options_["CIS_DOPDM_STATES"].size(); i++) {
        int state = options_["CIS_DOPDM_STATES"][i].to_integer();
        bool singlet = (state > 0);
        state = abs(state);

        std::shared_ptr<Matrix> D = Dao((singlet ? singlets_[state-1] : triplets_[state-1]), true);
        std::stringstream s;
        s << (singlet ? "S" : "T") << state << "_dD.dat";

        FILE* fh = fopen(s.str().c_str(),"w");
        fwrite((void*)D->pointer()[0],sizeof(double),nso_ * nso_,fh);
        fclose(fh);
    }
    for (unsigned int i = 0; i < options_["CIS_TOPDM_STATES"].size(); i++) {
        int state = options_["CIS_TOPDM_STATES"][i].to_integer();
        bool singlet = (state > 0);
        state = abs(state);

        std::shared_ptr<Matrix> D = TDao((singlet ? singlets_[state-1] : triplets_[state-1]), singlet);
        std::stringstream s;
        s << (singlet ? "S" : "T") << state << "_TD.dat";

        FILE* fh = fopen(s.str().c_str(),"w");
        fwrite((void*)D->pointer()[0],sizeof(double),nso_ * nso_,fh);
        fclose(fh);
    }
    for (unsigned int i = 0; i < options_["CIS_NO_STATES"].size(); i++) {
        int state = options_["CIS_NO_STATES"][i].to_integer();
        bool singlet = (state > 0);
        state = abs(state);

        std::pair<std::shared_ptr<Matrix>, std::shared_ptr<Vector> > stuff = Nao((singlet ? singlets_[state-1] : triplets_[state-1]),false);
        std::stringstream s;
        s << (singlet ? "S" : "T") << state << "_N.dat";

        FILE* fh = fopen(s.str().c_str(),"w");
        fwrite((void*)stuff.first->pointer()[0],sizeof(double),nso_ * nmo_,fh);
        fwrite((void*)stuff.second->pointer(),sizeof(double),nmo_,fh);
        fclose(fh);
    }
    for (unsigned int i = 0; i < options_["CIS_AD_STATES"].size(); i++) {
        int state = options_["CIS_AD_STATES"][i].to_integer();
        bool singlet = (state > 0);
        state = abs(state);

        std::pair<std::shared_ptr<Matrix>, std::shared_ptr<Matrix> > stuff = ADao((singlet ? singlets_[state-1] : triplets_[state-1]));
        std::stringstream s;
        s << (singlet ? "S" : "T") << state << "_AD.dat";

        FILE* fh = fopen(s.str().c_str(),"w");
        fwrite((void*)stuff.first->pointer()[0],sizeof(double),nso_ * nso_,fh);
        fwrite((void*)stuff.second->pointer()[0],sizeof(double),nso_ * nso_,fh);
        fclose(fh);
    }
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
    SharedMatrix D(new Matrix("Dmo", reference_wavefunction_->nmopi(), reference_wavefunction_->nmopi()));

    int symm = T1->symmetry();

    /// Reference occupation (if not difference density)
    if (!diff) {
        for (int h = 0; h < D->nirrep(); ++h) {
//            int nmo = D->rowspi()[h];
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
std::pair<SharedMatrix, std::shared_ptr<Vector> > RCIS::Nmo(SharedMatrix T1, bool diff)
{
    SharedMatrix D = Dmo(T1, diff);
    SharedMatrix C(new Matrix("Nmo", D->nirrep(), D->rowspi(), D->rowspi()));
    std::shared_ptr<Vector> O(new Vector("Occupation", D->rowspi()));

    D->diagonalize(C,O,descending);

    return std::make_pair(C,O);
}
std::pair<SharedMatrix, std::shared_ptr<Vector> > RCIS::Nso(SharedMatrix T1, bool diff)
{
    std::pair<SharedMatrix, std::shared_ptr<Vector> > pair = Nmo(T1,diff);
    SharedMatrix N = pair.first;
    std::shared_ptr<Vector> O = pair.second;

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
    return std::make_pair(N2,O);
}
std::pair<SharedMatrix, std::shared_ptr<Vector> > RCIS::Nao(SharedMatrix T1, bool diff)
{
    std::pair<SharedMatrix, std::shared_ptr<Vector> > pair = Nso(T1,diff);
    SharedMatrix N = pair.first;
    std::shared_ptr<Vector> O = pair.second;

    SharedMatrix N2(new Matrix("Nso", C_->nrow(), C_->ncol()));
    SharedMatrix N3(new Matrix("Nso", C_->nrow(), C_->ncol()));
    std::shared_ptr<Vector> O2(new Vector("Occupation", C_->ncol()));

    int offset = 0;
    std::vector<std::pair<double,int> > index;
    for (int h = 0; h < C_->nirrep(); h++) {

        int ncol = C_->ncol();
        int nmo = C_->colspi()[h];
        int nso = AO2SO_->colspi()[h];
        int nao = AO2SO_->rowspi()[h];

        if (!nmo || !nso || !nao) continue;

        for (int i = 0; i < nmo; i++) {
            index.push_back(std::make_pair(O->get(h,i),i+offset));
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

    return std::make_pair(N3,O2);
}
std::pair<SharedMatrix, SharedMatrix> RCIS::ADmo(SharedMatrix T1)
{
    std::pair<SharedMatrix, SharedVector> nos = Nmo(T1, true);
    SharedMatrix N = nos.first;
    SharedVector f = nos.second;

    SharedMatrix A(new Matrix("A", N->rowspi(), N->rowspi()));
    SharedMatrix D(new Matrix("D", N->rowspi(), N->rowspi()));
    for (int h = 0; h < N->nirrep(); h++) {
        int nrow = N->rowspi()[h];
        int ncol = N->colspi()[h];
        if (!nrow || !ncol) continue;
        double** Np = N->pointer(h);
        double** Ap = A->pointer(h);
        double** Dp = D->pointer(h);
        double*  fp = f->pointer(h);

        int nA = 0;
        for (int i = 0; i < ncol; i++) {
            if (fp[i] < 0.0) {
                break;
            }
            nA++;
        }

        int nD = ncol - nA;

        // Attach
        for (int i = 0; i < nA; i++) {
            C_DSCAL(nrow,sqrt(fp[i]),&Np[0][i],ncol);
        }
        C_DGEMM('N','T',nrow,nrow,nA,1.0,&Np[0][0],ncol,&Np[0][0],ncol,0.0,Ap[0],nrow);

        // Detach
        for (int i = nA; i < ncol; i++) {
            C_DSCAL(nrow,sqrt(-fp[i]),&Np[0][i],ncol);
        }
        C_DGEMM('N','T',nrow,nrow,nD,1.0,&Np[0][nA],ncol,&Np[0][nA],ncol,0.0,Dp[0],nrow);

    }

    return std::make_pair(A,D);
}
std::pair<SharedMatrix, SharedMatrix> RCIS::ADso(SharedMatrix T1)
{
    std::pair<SharedMatrix, SharedVector> nos = Nso(T1, true);
    SharedMatrix N = nos.first;
    SharedVector f = nos.second;

    SharedMatrix A(new Matrix("A", N->rowspi(), N->rowspi()));
    SharedMatrix D(new Matrix("D", N->rowspi(), N->rowspi()));
    for (int h = 0; h < N->nirrep(); h++) {
        int nrow = N->rowspi()[h];
        int ncol = N->colspi()[h];
        if (!nrow || !ncol) continue;
        double** Np = N->pointer(h);
        double** Ap = A->pointer(h);
        double** Dp = D->pointer(h);
        double*  fp = f->pointer(h);

        int nA = 0;
        for (int i = 0; i < ncol; i++) {
            if (fp[i] < 0.0) {
                break;
            }
            nA++;
        }

        int nD = ncol - nA;

        // Attach
        for (int i = 0; i < nA; i++) {
            C_DSCAL(nrow,sqrt(fp[i]),&Np[0][i],ncol);
        }
        C_DGEMM('N','T',nrow,nrow,nA,1.0,&Np[0][0],ncol,&Np[0][0],ncol,0.0,Ap[0],nrow);

        // Detach
        for (int i = nA; i < ncol; i++) {
            C_DSCAL(nrow,sqrt(-fp[i]),&Np[0][i],ncol);
        }
        C_DGEMM('N','T',nrow,nrow,nD,1.0,&Np[0][nA],ncol,&Np[0][nA],ncol,0.0,Dp[0],nrow);

    }

    return std::make_pair(A,D);
}
std::pair<SharedMatrix, SharedMatrix> RCIS::ADao(SharedMatrix T1)
{
    std::pair<SharedMatrix, SharedVector> nos = Nao(T1, true);
    SharedMatrix N = nos.first;
    SharedVector f = nos.second;

    SharedMatrix A(new Matrix("A", N->rowspi(), N->rowspi()));
    SharedMatrix D(new Matrix("D", N->rowspi(), N->rowspi()));
    for (int h = 0; h < N->nirrep(); h++) {
        int nrow = N->rowspi()[h];
        int ncol = N->colspi()[h];
        if (!nrow || !ncol) continue;
        double** Np = N->pointer(h);
        double** Ap = A->pointer(h);
        double** Dp = D->pointer(h);
        double*  fp = f->pointer(h);

        int nA = 0;
        for (int i = 0; i < ncol; i++) {
            if (fp[i] < 0.0) {
                break;
            }
            nA++;
        }

        int nD = ncol - nA;

        // Attach
        for (int i = 0; i < nA; i++) {
            C_DSCAL(nrow,sqrt(fp[i]),&Np[0][i],ncol);
        }
        C_DGEMM('N','T',nrow,nrow,nA,1.0,&Np[0][0],ncol,&Np[0][0],ncol,0.0,Ap[0],nrow);

        // Detach
        for (int i = nA; i < ncol; i++) {
            C_DSCAL(nrow,sqrt(-fp[i]),&Np[0][i],ncol);
        }
        C_DGEMM('N','T',nrow,nrow,nD,1.0,&Np[0][nA],ncol,&Np[0][nA],ncol,0.0,Dp[0],nrow);

    }

    return std::make_pair(A,D);
}

double RCIS::compute_energy()
{
    // Main CIS Header
    print_header();

    if (!jk_)
        preiterations();

    // Construct components
    std::shared_ptr<CISRHamiltonian> H(new CISRHamiltonian(jk_, Caocc_,Cavir_,eps_aocc_,eps_avir_));
    std::shared_ptr<DLRSolver> solver;
    if (options_.get_str("SOLVER_TYPE") == "DL")
        solver = DLRSolver::build_solver(options_,H);
    else if (options_.get_str("SOLVER_TYPE") == "RAYLEIGH")
        solver = RayleighRSolver::build_solver(options_,H);

    // Extra Knobs
    H->set_print(print_);
    H->set_debug(debug_);
    H->set_bench(bench_);
    H->set_exact_diagonal(options_.get_bool("SOLVER_EXACT_DIAGONAL"));
    solver->set_convergence(convergence_);

    // Initialization/Memory
    solver->initialize();

    // Component Headers
    solver->print_header();
    H->print_header();
    jk_->print_header();

    // Singlets
    if (options_.get_bool("DO_SINGLETS")) {

        H->set_singlet(true);

        if (print_) {
            outfile->Printf( "  ==> Singlets <==\n\n");
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
        const std::vector<std::shared_ptr<Vector> > singlets = solver->eigenvectors();
        const std::vector<std::vector<double> > E_singlets = solver->eigenvalues();

        std::vector<std::shared_ptr<Matrix> > evec_temp;
        std::vector<std::pair<double, int> > eval_temp;

        for (size_t N = 0, index = 0; N < singlets.size(); ++N) {
            std::vector<SharedMatrix > t = H->unpack(singlets[N]);
            for (int h = 0; h < Caocc_->nirrep(); h++) {
                // Spurious zero eigenvalue due to not enough states
                if (N >= (size_t)singlets[N]->dimpi()[h]) continue;
                evec_temp.push_back(t[h]);
                eval_temp.push_back(std::make_pair(E_singlets[N][h], index));
                index++;
            }
        }

        std::sort(eval_temp.begin(), eval_temp.end());

        singlets_.clear();
        E_singlets_.clear();

        for (size_t i = 0; i < eval_temp.size(); i++) {
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
            outfile->Printf( "  ==> Triplets <==\n\n");
        }

        solver->solve();

        const std::vector<std::shared_ptr<Vector> > triplets = solver->eigenvectors();
        const std::vector<std::vector<double> > E_triplets = solver->eigenvalues();

        std::vector<std::shared_ptr<Matrix> > evec_temp;
        std::vector<std::pair<double, int> > eval_temp;

        for (size_t N = 0, index = 0; N < triplets.size(); ++N) {
            std::vector<SharedMatrix > t = H->unpack(triplets[N]);
            for (int h = 0; h < Caocc_->nirrep(); h++) {
                // Spurious zero eigenvalue due to not enough states
                if (N >= (size_t)triplets[N]->dimpi()[h]) continue;
                evec_temp.push_back(t[h]);
                eval_temp.push_back(std::make_pair(E_triplets[N][h], index));
                index++;
            }
        }

        std::sort(eval_temp.begin(), eval_temp.end());

        triplets_.clear();
        E_triplets_.clear();

        for (size_t i = 0; i < eval_temp.size(); i++) {
            E_triplets_.push_back(eval_temp[i].first);
            triplets_.push_back(evec_temp[eval_temp[i].second]);
        }

    }

    // Finalize solver
    solver->finalize();

    // Print wavefunctions and properties
    sort_states();
    print_wavefunctions();
    print_amplitudes();
    print_transitions();
    print_densities();

    return 0.0;
}

RTDHF::RTDHF(SharedWavefunction ref_wfn, Options& options) :
           RBase(ref_wfn, options)
{
}
RTDHF::~RTDHF()
{
}
void RTDHF::print_header()
{
    outfile->Printf( "\n");
    outfile->Printf( "         ------------------------------------------------------------\n");
    outfile->Printf( "                                      TDHF                           \n");
    outfile->Printf( "                                  Rob Parrish                       \n");
    outfile->Printf( "         ------------------------------------------------------------\n\n");

    outfile->Printf( "  ==> Geometry <==\n\n");
    molecule_->print();
    outfile->Printf( "  Nuclear repulsion = %20.15f\n", basisset_->molecule()->nuclear_repulsion_energy());
    outfile->Printf( "  Reference energy  = %20.15f\n\n", Eref_);

    outfile->Printf( "  ==> Basis Set <==\n\n");
    basisset_->print_by_level("outfile", print_);

    if (debug_ > 1) {
        outfile->Printf( "  ==> Fock Matrix (MO Basis) <==\n\n");
        eps_aocc_->print();
        eps_avir_->print();
    }
}
double RTDHF::compute_energy()
{
    // Main TDHF Header
    print_header();

    if (!jk_)
        preiterations();

    // Construct components
    std::shared_ptr<TDHFRHamiltonian> H(new TDHFRHamiltonian(jk_, Caocc_,Cavir_,eps_aocc_,eps_avir_));
    std::shared_ptr<DLRXSolver> solver = DLRXSolver::build_solver(options_,H);

    // Extra Knobs
    H->set_print(print_);
    H->set_debug(debug_);
    solver->set_convergence(convergence_);

    // Initialization
    solver->initialize();

    // Component Headers
    solver->print_header();
    H->print_header();
    jk_->print_header();

    // Singlets
    if (options_.get_bool("DO_SINGLETS")) {

        H->set_singlet(true);

        if (print_) {
            outfile->Printf( "  ==> Singlets <==\n\n");
        }

        solver->solve();

        // TODO Unpack
        //const std::vector<std::shared_ptr<Vector> > singlets = solver->eigenvectors();
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
            outfile->Printf( "  ==> Triplets <==\n\n");
        }

        solver->solve();

        // TODO Unpack
        //const std::vector<std::shared_ptr<Vector> > triplets = solver->eigenvectors();
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

    // Finalize solver
    solver->finalize();

    // TODO
    // Print wavefunctions and properties
    //print_wavefunctions();
    //print_amplitudes();
    //print_transitions();
    //print_densities();

    return 0.0;
}

RCPKS::RCPKS(SharedWavefunction ref_wfn, Options& options) :
             RCPHF(ref_wfn, options)
{
}
RCPKS::~RCPKS()
{
}
void RCPKS::print_header()
{
    outfile->Printf( "\n");
    outfile->Printf( "         ------------------------------------------------------------\n");
    outfile->Printf( "                                     CPKS                            \n");
    outfile->Printf( "                                  Rob Parrish                       \n");
    outfile->Printf( "         ------------------------------------------------------------\n\n");

    outfile->Printf( "  ==> Geometry <==\n\n");
    molecule_->print();
    outfile->Printf( "  Nuclear repulsion = %20.15f\n", basisset_->molecule()->nuclear_repulsion_energy());
    outfile->Printf( "  Reference energy  = %20.15f\n\n", Eref_);

    outfile->Printf( "  ==> Basis Set <==\n\n");
    basisset_->print_by_level("outfile", print_);
}
double RCPKS::compute_energy()
{
    // Main CPKS Header
    print_header();

    // Add named tasks to the force vector list
    add_named_tasks();

    if (!jk_ || !v_)
        preiterations();

    // Construct components
    std::shared_ptr<CPKSRHamiltonian> H(new CPKSRHamiltonian(jk_,v_,Cocc_,Caocc_,Cavir_,eps_aocc_,eps_avir_));
    std::shared_ptr<CGRSolver> solver = CGRSolver::build_solver(options_,H);

    // Extra Knobs
    H->set_print(print_);
    H->set_debug(debug_);
    solver->set_convergence(convergence_);

    // Addition of force vectors
    std::vector<SharedVector>& bref = solver->b();
    std::map<std::string, SharedVector> b = H->pack(b_);
    for (std::map<std::string, SharedVector>::const_iterator it = b.begin();
        it != b.end(); ++it) {
        bref.push_back((*it).second);
    }

    // Initialization/Memory
    solver->initialize();

    // Component Headers
    solver->print_header();
    H->print_header();
    jk_->print_header();

    if (print_) {
        outfile->Printf( "  ==> CPHF Iterations <==\n\n");
    }

    if (options_.get_bool("EXPLICIT_HAMILTONIAN")) {
        SharedMatrix A = H->explicit_hamiltonian();
        A->print();
    }

    if (debug_) {
        for (std::map<std::string, SharedMatrix>::const_iterator it = b_.begin();
            it != b_.end(); ++it) {
            (*it).second->print();
        }
    }

    solver->solve();

    std::vector<SharedMatrix> x1 = H->unpack(solver->x());

    int index = 0;
    for (std::map<std::string, SharedMatrix>::const_iterator it = b_.begin();
        it != b_.end(); ++it) {
        x_[(*it).first] = x1[index++];
    }

    if (debug_) {
        for (std::map<std::string, SharedMatrix>::const_iterator it = x_.begin();
            it != x_.end(); ++it) {
            (*it).second->print();
        }
    }

    analyze_named_tasks();

    solver->finalize();

    return 0.0;
}

RTDA::RTDA(SharedWavefunction ref_wfn, Options& options) : RCIS(ref_wfn, options)
{
}
RTDA::~RTDA()
{
}
void RTDA::print_header()
{
    outfile->Printf( "\n");
    outfile->Printf( "         ------------------------------------------------------------\n");
    outfile->Printf( "                                      TDA                            \n");
    outfile->Printf( "                                  Rob Parrish                       \n");
    outfile->Printf( "         ------------------------------------------------------------\n\n");

    outfile->Printf( "  ==> Geometry <==\n\n");
    molecule_->print();
    outfile->Printf( "  Nuclear repulsion = %20.15f\n", basisset_->molecule()->nuclear_repulsion_energy());
    outfile->Printf( "  Reference energy  = %20.15f\n\n", Eref_);

    outfile->Printf( "  ==> Basis Set <==\n\n");
    basisset_->print_by_level("outfile", print_);
}
double RTDA::compute_energy()
{
    // Main TDA Header
    print_header();

    if (!jk_ || !v_)
        preiterations();

    // Construct components
    std::shared_ptr<TDARHamiltonian> H(new TDARHamiltonian(jk_,v_,Cocc_,Caocc_,Cavir_,eps_aocc_,eps_avir_));
    std::shared_ptr<DLRSolver> solver = DLRSolver::build_solver(options_,H);

    // Extra Knobs
    H->set_print(print_);
    H->set_debug(debug_);

    // Initialization/Memory
    solver->initialize();
    solver->set_convergence(convergence_);

    // Component Headers
    solver->print_header();
    H->print_header();
    jk_->print_header();

    // Singlets
    if (options_.get_bool("DO_SINGLETS")) {

        H->set_singlet(true);

        if (print_) {
            outfile->Printf( "  ==> Singlets <==\n\n");
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
        const std::vector<std::shared_ptr<Vector> > singlets = solver->eigenvectors();
        const std::vector<std::vector<double> > E_singlets = solver->eigenvalues();

        std::vector<std::shared_ptr<Matrix> > evec_temp;
        std::vector<std::pair<double, int> > eval_temp;

        for (size_t N = 0, index = 0; N < singlets.size(); ++N) {
            std::vector<SharedMatrix > t = H->unpack(singlets[N]);
            for (int h = 0; h < Caocc_->nirrep(); h++) {
                // Spurious zero eigenvalue due to not enough states
                if (N >= (size_t)singlets[N]->dimpi()[h]) continue;
                evec_temp.push_back(t[h]);
                eval_temp.push_back(std::make_pair(E_singlets[N][h], index));
                index++;
            }
        }

        std::sort(eval_temp.begin(), eval_temp.end());

        singlets_.clear();
        E_singlets_.clear();

        for (size_t i = 0; i < eval_temp.size(); i++) {
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
            outfile->Printf( "  ==> Triplets <==\n\n");
        }

        solver->solve();

        const std::vector<std::shared_ptr<Vector> > triplets = solver->eigenvectors();
        const std::vector<std::vector<double> > E_triplets = solver->eigenvalues();

        std::vector<std::shared_ptr<Matrix> > evec_temp;
        std::vector<std::pair<double, int> > eval_temp;

        for (size_t N = 0, index = 0; N < triplets.size(); ++N) {
            std::vector<SharedMatrix > t = H->unpack(triplets[N]);
            for (int h = 0; h < Caocc_->nirrep(); h++) {
                // Spurious zero eigenvalue due to not enough states
                if (N >= (size_t)triplets[N]->dimpi()[h]) continue;
                evec_temp.push_back(t[h]);
                eval_temp.push_back(std::make_pair(E_triplets[N][h], index));
                index++;
            }
        }

        std::sort(eval_temp.begin(), eval_temp.end());

        triplets_.clear();
        E_triplets_.clear();

        for (size_t i = 0; i < eval_temp.size(); i++) {
            E_triplets_.push_back(eval_temp[i].first);
            triplets_.push_back(evec_temp[eval_temp[i].second]);
        }

    }

    // Finalize solver
    solver->finalize();

    // Print wavefunctions and properties
    sort_states();
    print_wavefunctions();
    print_amplitudes();
    print_transitions();
    print_densities();

    return 0.0;
}

RTDDFT::RTDDFT(SharedWavefunction ref_wfn, Options& options) :
              RTDHF(ref_wfn, options)
{
}
RTDDFT::~RTDDFT()
{
}
void RTDDFT::print_header()
{
    outfile->Printf( "\n");
    outfile->Printf( "         ------------------------------------------------------------\n");
    outfile->Printf( "                                     TDDFT                           \n");
    outfile->Printf( "                                  Rob Parrish                       \n");
    outfile->Printf( "         ------------------------------------------------------------\n\n");

    outfile->Printf( "  ==> Geometry <==\n\n");
    molecule_->print();
    outfile->Printf( "  Nuclear repulsion = %20.15f\n", basisset_->molecule()->nuclear_repulsion_energy());
    outfile->Printf( "  Reference energy  = %20.15f\n\n", Eref_);

    outfile->Printf( "  ==> Basis Set <==\n\n");
    basisset_->print_by_level("outfile", print_);
}
double RTDDFT::compute_energy()
{
    // Main TDDFT Header
    print_header();

    if (!jk_ || !v_)
        preiterations();

    // Construct components
    std::shared_ptr<TDDFTRHamiltonian> H(new TDDFTRHamiltonian(jk_,v_,Cocc_,Caocc_,Cavir_,eps_aocc_,eps_avir_));
    std::shared_ptr<DLRXSolver> solver = DLRXSolver::build_solver(options_,H);

    // Extra Knobs
    H->set_print(print_);
    H->set_debug(debug_);
    solver->set_convergence(convergence_);

    // Initialization
    solver->initialize();

    // Component Headers
    solver->print_header();
    H->print_header();
    jk_->print_header();

    // Singlets
    if (options_.get_bool("DO_SINGLETS")) {

        H->set_singlet(true);

        if (print_) {
            outfile->Printf( "  ==> Singlets <==\n\n");
        }

        solver->solve();

        // TODO Unpack
        //const std::vector<std::shared_ptr<Vector> > singlets = solver->eigenvectors();
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
            outfile->Printf( "  ==> Triplets <==\n\n");
        }

        solver->solve();

        // TODO Unpack
        //const std::vector<std::shared_ptr<Vector> > triplets = solver->eigenvectors();
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

    // Finalize solver
    solver->finalize();

    // TODO
    // Print wavefunctions and properties
    //print_wavefunctions();
    //print_amplitudes();
    //print_transitions();
    //print_densities();

    return 0.0;
}


}
