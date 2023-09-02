/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "apps.h"
#include "jk.h"
#include "v.h"
#include "hamiltonian.h"
#include "solver.h"

#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"
#include "psi4/physconst.h"
#include "psi4/libscf_solver/hf.h"
#include "psi4/libscf_solver/rhf.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/multipolesymmetry.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/process.h"

#include <algorithm>
#include <functional>
#include <tuple>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {

RBase::RBase(SharedWavefunction ref_wfn, Options& options, bool use_symmetry)
    : Wavefunction(options), use_symmetry_(use_symmetry) {
    shallow_copy(ref_wfn);

    set_reference(ref_wfn);

    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");
    bench_ = options_.get_int("BENCH");
    convergence_ = options_.get_double("SOLVER_CONVERGENCE");
}
RBase::RBase(bool flag) : Wavefunction(Process::environment.options) {
    psio_ = _default_psio_lib_;
    throw PSIEXCEPTION("DGAS: Lets not let RMP do dirty hacks!\n");
    outfile->Printf("Dirty hack %s\n\n", (flag ? "true" : "false"));
}
RBase::~RBase() { postiterations(); }
void RBase::set_reference(SharedWavefunction ref_wfn) {
    reference_wavefunction_ = ref_wfn;

    if (!reference_wavefunction_->same_a_b_dens()) {
        throw PSIEXCEPTION("RBase: Reference is not restricted");
    }

    Eref_ = reference_wavefunction_->energy();

    if (use_symmetry_) {
        Cocc_ = Ca_subset("SO", "OCC");
        Cfocc_ = Ca_subset("SO", "FROZEN_OCC");
        Caocc_ = Ca_subset("SO", "ACTIVE_OCC");
        Cavir_ = Ca_subset("SO", "ACTIVE_VIR");
        Cfvir_ = Ca_subset("SO", "FROZEN_VIR");
        eps_focc_ = epsilon_a_subset("SO", "FROZEN_OCC");
        eps_aocc_ = epsilon_a_subset("SO", "ACTIVE_OCC");
        eps_avir_ = epsilon_a_subset("SO", "ACTIVE_VIR");
        eps_fvir_ = epsilon_a_subset("SO", "FROZEN_VIR");
    } else {
        Cocc_ = Ca_subset("AO", "OCC");
        Cfocc_ = Ca_subset("AO", "FROZEN_OCC");
        Caocc_ = Ca_subset("AO", "ACTIVE_OCC");
        Cavir_ = Ca_subset("AO", "ACTIVE_VIR");
        Cfvir_ = Ca_subset("AO", "FROZEN_VIR");
        eps_focc_ = epsilon_a_subset("AO", "FROZEN_OCC");
        eps_aocc_ = epsilon_a_subset("AO", "ACTIVE_OCC");
        eps_avir_ = epsilon_a_subset("AO", "ACTIVE_VIR");
        eps_fvir_ = epsilon_a_subset("AO", "FROZEN_VIR");
    }

    std::vector<SharedMatrix> Cs;
    Cs.push_back(Cfocc_);
    Cs.push_back(Caocc_);
    Cs.push_back(Cavir_);
    Cs.push_back(Cfvir_);
    C_ = linalg::horzcat(Cs);
}
void RBase::preiterations() {
    if (!jk_) {
        if (options_.get_bool("SAVE_JK")) {
            jk_ = (static_cast<psi::scf::HF*>(reference_wavefunction_.get()))->jk();
            outfile->Printf("    Reusing JK object from SCF.\n\n");
        } else {
            size_t effective_memory = (size_t)(0.125 * options_.get_double("CPHF_MEM_SAFETY_FACTOR") * memory_);
            jk_ = JK::build_JK(basisset_, get_basisset("DF_BASIS_SCF"), options_, false, effective_memory);
            jk_->set_memory(effective_memory);
            jk_->initialize();
        }
    }

    if (!v_) {
        if (options_.get_str("MODULE") == "RCPKS" || options_.get_str("MODULE") == "RTDA" ||
            options_.get_str("MODULE") == "RTDDFT") {
            throw PSIEXCEPTION("V is not currently enabled in apps.cc");
            // v_ = VBase::build_V(basisset_, options_, "RK");
            // v_->initialize();
        }
    }
}
void RBase::postiterations() { jk_.reset(); }

RCPHF::RCPHF(SharedWavefunction ref_wfn, Options& options, bool use_symmetry) : RBase(ref_wfn, options, use_symmetry) {}
RCPHF::~RCPHF() {}
void RCPHF::print_header() {
    outfile->Printf("\n");
    outfile->Printf("         ------------------------------------------------------------\n");
    outfile->Printf("                                     CPHF                           \n");
    outfile->Printf("                                  Rob Parrish                       \n");
    outfile->Printf("         ------------------------------------------------------------\n\n");

    outfile->Printf("  ==> Geometry <==\n\n");
    molecule_->print();
    outfile->Printf("  Nuclear repulsion = %20.15f\n",
                    basisset_->molecule()->nuclear_repulsion_energy(dipole_field_strength_));
    outfile->Printf("  Reference energy  = %20.15f\n\n", Eref_);

    outfile->Printf("  ==> Basis Set <==\n\n");
    basisset_->print_by_level("outfile", print_);

    if (tasks_.size()) {
        outfile->Printf("  ==> Named Tasks <==\n\n");
        for (std::set<std::string>::const_iterator it = tasks_.begin(); it != tasks_.end(); ++it) {
            outfile->Printf("    %s\n", (*it).c_str());
        }
        outfile->Printf("\n");
    }

    if (debug_ > 1) {
        outfile->Printf("  ==> Fock Matrix (MO Basis) <==\n\n");
        eps_aocc_->print();
        eps_avir_->print();
    }
}
void RCPHF::add_task(const std::string& task) { tasks_.insert(task); }
void RCPHF::add_named_tasks() {
    if (tasks_.count("POLARIZABILITY")) {
        add_polarizability();
    }
}
void RCPHF::analyze_named_tasks() {
    if (tasks_.count("POLARIZABILITY")) {
        analyze_polarizability();
    }
}
void RCPHF::add_polarizability() {
    OperatorSymmetry msymm(1, molecule_, integral_, factory_);
    std::vector<SharedMatrix> dipole = msymm.create_matrices("SO Dipole");
    std::shared_ptr<OneBodySOInt> ints(integral_->so_dipole());
    ints->compute(dipole);

    for (size_t i = 0; i < dipole.size(); i++) {
        std::stringstream ss;
        ss << "Dipole Perturbation " << (i == 0 ? "X" : (i == 1 ? "Y" : "Z"));
        auto B = std::make_shared<Matrix>(ss.str(), Caocc_->colspi(), Cavir_->colspi(), dipole[i]->symmetry());

        int symm = dipole[i]->symmetry();
        double* temp = new double[dipole[i]->max_nrow() * Cavir_->max_ncol()];

        for (int h = 0; h < B->nirrep(); h++) {
            int nsol = dipole[i]->rowspi()[h];
            int nsor = dipole[i]->colspi()[h ^ symm];
            int noccl = Caocc_->colspi()[h];
            int nvirr = Cavir_->colspi()[h ^ symm];

            if (!nsol || !nsor || !noccl || !nvirr) continue;

            double** dp = dipole[i]->pointer(h);
            double** bp = B->pointer(h);
            double** Clp = Caocc_->pointer(h);
            double** Crp = Cavir_->pointer(h ^ symm);

            C_DGEMM('N', 'N', nsol, nvirr, nsor, 1.0, dp[0], nsor, Crp[0], nvirr, 0.0, temp, nvirr);
            C_DGEMM('T', 'N', noccl, nvirr, nsol, 1.0, Clp[0], noccl, temp, nvirr, 0.0, bp[0], nvirr);
        }

        delete[] temp;

        std::stringstream ss2;
        ss2 << (i == 0 ? "MU_X" : (i == 1 ? "MU_Y" : "MU_Z"));
        b_[ss2.str()] = B;
    }
}
void RCPHF::analyze_polarizability() {
    std::vector<SharedMatrix> u;
    std::vector<SharedMatrix> d;

    d.push_back(b_["MU_X"]);
    d.push_back(b_["MU_Y"]);
    d.push_back(b_["MU_Z"]);

    u.push_back(x_["MU_X"]);
    u.push_back(x_["MU_Y"]);
    u.push_back(x_["MU_Z"]);

    // Analysis
    auto polarizability = std::make_shared<Matrix>("CPHF Polarizability", 3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            polarizability->set(0, i, j, -4.0 * (d[i]->symmetry() == u[j]->symmetry() ? d[i]->vector_dot(u[j]) : 0.0));
        }
    }

    polarizability->print();
}
double RCPHF::compute_energy() {
    // Main CPHF Header
    print_header();

    // Add named tasks to the force vector list
    add_named_tasks();

    if (!jk_) preiterations();

    // Construct components
    auto H = std::make_shared<CPHFRHamiltonian>(jk_, Caocc_, Cavir_, eps_aocc_, eps_avir_);
    std::shared_ptr<CGRSolver> solver = CGRSolver::build_solver(options_, H);

    // Extra Knobs
    H->set_print(print_);
    H->set_debug(debug_);
    solver->set_convergence(convergence_);

    // Addition of force vectors
    std::vector<SharedVector>& bref = solver->b();
    std::map<std::string, SharedVector> b = H->pack(b_);
    for (std::map<std::string, SharedVector>::const_iterator it = b.begin(); it != b.end(); ++it) {
        bref.push_back((*it).second);
    }

    // Initialization/Memory
    solver->initialize();

    // Component Headers
    solver->print_header();
    H->print_header();
    jk_->print_header();

    if (print_) {
        outfile->Printf("  ==> CPHF Iterations <==\n\n");
    }

    if (options_.get_bool("EXPLICIT_HAMILTONIAN")) {
        SharedMatrix A = H->explicit_hamiltonian();
        A->print();
    }

    if (debug_) {
        for (std::map<std::string, SharedMatrix>::const_iterator it = b_.begin(); it != b_.end(); ++it) {
            (*it).second->print();
        }
    }

    solver->solve();

    std::vector<SharedMatrix> x1 = H->unpack(solver->x());

    int index = 0;
    for (std::map<std::string, SharedMatrix>::const_iterator it = b_.begin(); it != b_.end(); ++it) {
        x_[(*it).first] = x1[index++];
    }

    if (debug_) {
        for (std::map<std::string, SharedMatrix>::const_iterator it = x_.begin(); it != x_.end(); ++it) {
            (*it).second->print();
        }
    }

    analyze_named_tasks();

    solver->finalize();

    return 0.0;
}
}  // namespace psi
