/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include "solver.h"

#include <algorithm>
#include <cmath>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

#include "points.h"
#include "hamiltonian.h"
#include "jk.h"

namespace psi {

Solver::Solver() { common_init(); }
Solver::~Solver() {}
void Solver::common_init() {
    print_ = 2;
    debug_ = 0;
    bench_ = 0;
    // Unlimited default
    memory_ = 0L;
    converged_ = false;
    iteration_ = 0;
    convergence_ = 0.0;
    criteria_ = 1.0E-6;
    maxiter_ = 100;
    precondition_ = "JACOBI";
    name_ = "Solver";
}

RSolver::RSolver(std::shared_ptr<RHamiltonian> H) : Solver(), H_(H) { name_ = "RSolver"; }
RSolver::~RSolver() {}

USolver::USolver(std::shared_ptr<UHamiltonian> H) : Solver(), H_(H) { name_ = "USolver"; }
USolver::~USolver() {}

CGRSolver::CGRSolver(std::shared_ptr<RHamiltonian> H) : RSolver(H) {
    nguess_ = 0;
    name_ = "CGR";
}
CGRSolver::~CGRSolver() {}
std::shared_ptr<CGRSolver> CGRSolver::build_solver(Options& options, std::shared_ptr<RHamiltonian> H) {
    auto solver = std::make_shared<CGRSolver>(H);

    if (options["PRINT"].has_changed()) {
        solver->set_print(options.get_int("PRINT") + 1);
    }
    if (options["DEBUG"].has_changed()) {
        solver->set_debug(options.get_int("DEBUG"));
    }
    if (options["BENCH"].has_changed()) {
        solver->set_bench(options.get_int("BENCH"));
    }
    if (options["SOLVER_PRECONDITION"].has_changed()) {
        if (options.get_str("SOLVER_PRECONDITION") == "SUBSPACE") {
            outfile->Printf("  !!!Warning!!!\n");
            outfile->Printf("  The subspace preconditioner has been broken for some time and was removed in Psi4 1.4.\n");
            outfile->Printf("  Setting preconditioner to Jacobi instead.\n\n");
            solver->set_precondition("JACOBI");
        } else {
            solver->set_precondition(options.get_str("SOLVER_PRECONDITION"));
        }
    } else if (options["SOLVER_MAXITER"].has_changed()) {
        solver->set_maxiter(options.get_int("SOLVER_MAXITER"));
    }
    if (options["SOLVER_CONVERGENCE"].has_changed()) {
        solver->set_convergence(options.get_double("SOLVER_CONVERGENCE"));
    }
    if (options["SOLVER_N_GUESS"].has_changed()) {
        solver->set_nguess(options.get_int("SOLVER_N_GUESS"));
    }

    return solver;
}
void CGRSolver::print_header() const {
    if (print_) {
        outfile->Printf("  ==> CGRSolver (by Rob Parrish) <==\n\n");
        outfile->Printf("   Number of roots    = %9zu\n", b_.size());
        outfile->Printf("   Preconditioning    = %9s\n", precondition_.c_str());
        outfile->Printf("   Convergence cutoff = %9.0E\n", criteria_);
        outfile->Printf("   Maximum iterations = %9d\n\n", maxiter_);
    }
}
size_t CGRSolver::memory_estimate() {
    size_t dimension = 0L;
    if (!diag_) diag_ = H_->diagonal();
    for (int h = 0; h < diag_->nirrep(); h++) {
        dimension += diag_->dimpi()[h];
    }
    return (6L * b_.size()) * dimension;
}
void CGRSolver::initialize() {
    finalize();

    int nvec = b_.size();
    for (int N = 0; N < nvec; ++N) {
        std::stringstream xs;
        xs << "Solution Vector " << N + 1;
        x_.push_back(std::make_shared<Vector>(xs.str(), b_[0]->dimpi()));
        std::stringstream Aps;
        Aps << "Product Vector " << N + 1;
        Ap_.push_back(std::make_shared<Vector>(Aps.str(), b_[0]->dimpi()));
        std::stringstream zs;
        zs << "Z Vector " << N + 1;
        z_.push_back(std::make_shared<Vector>(zs.str(), b_[0]->dimpi()));
        std::stringstream rs;
        rs << "Residual Vector " << N + 1;
        r_.push_back(std::make_shared<Vector>(rs.str(), b_[0]->dimpi()));
        std::stringstream ps;
        ps << "Conjugate Vector " << N + 1;
        p_.push_back(std::make_shared<Vector>(ps.str(), b_[0]->dimpi()));
        alpha_.push_back(0.0);
        beta_.push_back(0.0);
        r_nrm2_.push_back(0.0);
        z_r_.push_back(0.0);
        r_converged_.push_back(false);
    }

    diag_ = H_->diagonal();
}
void CGRSolver::solve() {
    iteration_ = 0;
    converged_ = false;
    nconverged_ = 0;
    convergence_ = 0.0;

    if (print_ > 1) {
        outfile->Printf("  => Iterations <=\n\n");
        outfile->Printf("  %10s %4s %10s %10s %11s\n", "", "Iter", "Converged", "Remaining", "Residual");
    }

    setup();
    guess();
    products_x();
    residual();
    update_z();
    update_p();

    do {
        iteration_++;

        products_p();
        alpha();
        update_x();
        update_r();
        check_convergence();
        if (print_) {
            outfile->Printf("  %-10s %4d %10d %10zu %11.3E\n", name_.c_str(), iteration_, nconverged_,
                            b_.size() - nconverged_, convergence_);
        }
        update_z();
        beta();
        update_p();

    } while (iteration_ < maxiter_ && !converged_);

    if (print_ > 1) {
        outfile->Printf("\n");
        if (!converged_) {
            outfile->Printf("    %sSolver did not converge.\n\n", name_.c_str());
        } else {
            outfile->Printf("    %sSolver converged.\n\n", name_.c_str());
        }
    }
}
void CGRSolver::finalize() {
    Ap_.clear();
    z_.clear();
    r_.clear();
    p_.clear();
    alpha_.clear();
    beta_.clear();
    r_nrm2_.clear();
    z_r_.clear();
    r_converged_.clear();
    diag_.reset();
}
void CGRSolver::setup() {
    if (shifts_.size() == 0) {
        shifts_ = std::vector<std::vector<double>>(diag_->nirrep());
        for (int h = 0; h < diag_->nirrep(); h++) {
            shifts_[h] = std::vector<double>(b_.size(), 0.0);
        }
    }

}
void CGRSolver::guess() {
    for (size_t N = 0; N < b_.size(); ++N) {
        for (int h = 0; h < b_[N]->nirrep(); ++h) {
            int n = b_[N]->dimpi()[h];
            if (!n) continue;
            auto bp = b_[N]->pointer();
            auto xp = x_[N]->pointer();
            auto dp = diag_->pointer();
            if (precondition_ == "JACOBI") {
                double lambda = shifts_[h][N];
                for (int i = 0; i < n; ++i) {
                    xp[i] = bp[i] / (dp[i] - lambda);
                }
            } else {
                for (int i = 0; i < n; ++i) {
                    xp[i] = bp[i];
                }
            }
        }
    }

    if (debug_) {
        diag_->print();
        for (size_t N = 0; N < b_.size(); N++) {
            x_[N]->print();
            b_[N]->print();
        }
    }
}
void CGRSolver::residual() {
    for (size_t N = 0; N < b_.size(); ++N) {
        r_[N]->copy(*Ap_[N]);
        r_[N]->scale(-1.0);
        r_[N]->add(*b_[N]);
    }

    if (debug_) {
        outfile->Printf("  > Residuals x <\n\n");
        for (size_t N = 0; N < r_.size(); N++) {
            r_[N]->print();
        }
    }
}
void CGRSolver::products_x() {
    H_->product(x_, Ap_);

    for (int h = 0; h < diag_->nirrep(); h++) {
        for (size_t i = 0; i < x_.size(); i++) {
            if (shifts_[h][i] != 0.0) {
                double lambda = shifts_[h][i];
                C_DAXPY(diag_->dimpi()[h], -lambda, x_[i]->pointer(h), 1, Ap_[i]->pointer(h), 1);
            }
        }
    }

    if (debug_) {
        outfile->Printf("  > Products x <\n\n");
        for (size_t N = 0; N < Ap_.size(); N++) {
            Ap_[N]->print();
        }
    }
}
void CGRSolver::products_p() {
    std::vector<std::shared_ptr<Vector>> p;
    std::vector<std::shared_ptr<Vector>> Ap;

    for (size_t N = 0; N < b_.size(); ++N) {
        if (r_converged_[N]) continue;
        p.push_back(p_[N]);
        Ap.push_back(Ap_[N]);
    }

    H_->product(p, Ap);

    for (size_t N = 0; N < b_.size(); ++N) {
        if (r_converged_[N]) continue;
        for (int h = 0; h < diag_->nirrep(); h++) {
            if (shifts_[h][N] != 0.0) {
                double lambda = shifts_[h][N];
                C_DAXPY(diag_->dimpi()[h], -lambda, p_[N]->pointer(h), 1, Ap_[N]->pointer(h), 1);
            }
        }
    }

    if (debug_) {
        outfile->Printf("  > Products p <\n\n");
        for (size_t N = 0; N < Ap_.size(); N++) {
            Ap_[N]->print();
        }
    }
}
void CGRSolver::alpha() {
    for (size_t N = 0; N < b_.size(); ++N) {
        if (r_converged_[N]) continue;
        z_r_[N] = 0.0;
        double p_Ap = 0.0;
        for (int h = 0; h < b_[N]->nirrep(); ++h) {
            int n = b_[N]->dimpi()[h];
            if (!n) continue;
            double* rp = r_[N]->pointer(h);
            double* zp = z_[N]->pointer(h);
            double* pp = p_[N]->pointer(h);
            double* App = Ap_[N]->pointer(h);
            z_r_[N] += C_DDOT(n, rp, 1, zp, 1);
            p_Ap += C_DDOT(n, pp, 1, App, 1);
        }
        alpha_[N] = z_r_[N] / p_Ap;
    }

    if (debug_) {
        outfile->Printf("  > Alpha <\n\n");
        for (size_t N = 0; N < alpha_.size(); N++) {
            outfile->Printf("Alpha %zu = %24.16E\n", N + 1, alpha_[N]);
        }
    }
}
void CGRSolver::update_x() {
    for (size_t N = 0; N < b_.size(); ++N) {
        if (r_converged_[N]) continue;
        for (int h = 0; h < b_[N]->nirrep(); ++h) {
            int n = b_[N]->dimpi()[h];
            if (!n) continue;
            double* xp = x_[N]->pointer(h);
            double* pp = p_[N]->pointer(h);
            C_DAXPY(n, alpha_[N], pp, 1, xp, 1);
        }
    }

    if (debug_) {
        outfile->Printf("  > Update x <\n\n");
        for (size_t N = 0; N < x_.size(); N++) {
            x_[N]->print();
        }
    }
}
void CGRSolver::update_r() {
    for (size_t N = 0; N < b_.size(); ++N) {
        if (r_converged_[N]) continue;
        for (int h = 0; h < b_[N]->nirrep(); ++h) {
            int n = b_[N]->dimpi()[h];
            if (!n) continue;
            double* rp = r_[N]->pointer(h);
            double* App = Ap_[N]->pointer(h);
            C_DAXPY(n, -alpha_[N], App, 1, rp, 1);
        }
    }

    if (debug_) {
        outfile->Printf("  > Update r <\n\n");
        for (size_t N = 0; N < r_.size(); N++) {
            r_[N]->print();
        }
    }
}
void CGRSolver::check_convergence() {
    convergence_ = 0.0;
    for (size_t N = 0; N < b_.size(); ++N) {
        if (r_converged_[N]) continue;
        double R2 = 0.0;
        double B2 = 0.0;
        for (int h = 0; h < b_[N]->nirrep(); ++h) {
            int n = b_[N]->dimpi()[h];
            if (!n) continue;
            double* rp = r_[N]->pointer(h);
            double* bp = b_[N]->pointer(h);
            B2 += C_DDOT(n, bp, 1, bp, 1);
            R2 += C_DDOT(n, rp, 1, rp, 1);
        }
        r_nrm2_[N] = sqrt(R2 / B2);
        if (convergence_ < r_nrm2_[N]) {
            convergence_ = r_nrm2_[N];
        }
        if (r_nrm2_[N] < criteria_) {
            r_converged_[N] = true;
            nconverged_++;
        }
    }
    if ((size_t)nconverged_ == b_.size()) {
        converged_ = true;
    }
}
void CGRSolver::update_z() {
    for (size_t N = 0; N < b_.size(); ++N) {
        if (r_converged_[N]) continue;
        for (int h = 0; h < b_[N]->nirrep(); ++h) {
            int n = b_[N]->dimpi()[h];
            if (!n) continue;
            double* zp = z_[N]->pointer();
            double* rp = r_[N]->pointer();
            double* dp = diag_->pointer();
            if (precondition_ == "JACOBI") {
                double lambda = shifts_[h][N];
                for (int i = 0; i < n; ++i) {
                    zp[i] = rp[i] / (dp[i] - lambda);
                }
            } else {
                for (int i = 0; i < n; ++i) {
                    zp[i] = rp[i];
                }
            }
        }
    }

    if (debug_) {
        outfile->Printf("  > Update z <\n\n");
        for (size_t N = 0; N < z_.size(); N++) {
            z_[N]->print();
        }
    }
}
void CGRSolver::beta() {
    for (size_t N = 0; N < b_.size(); ++N) {
        if (r_converged_[N]) continue;
        double zr = 0.0;
        for (int h = 0; h < b_[N]->nirrep(); ++h) {
            int n = b_[N]->dimpi()[h];
            if (!n) continue;
            double* rp = r_[N]->pointer();
            double* zp = z_[N]->pointer();
            zr += C_DDOT(n, rp, 1, zp, 1);
        }
        beta_[N] = zr / z_r_[N];
    }

    if (debug_) {
        outfile->Printf("  > Beta <\n\n");
        for (size_t N = 0; N < beta_.size(); N++) {
            outfile->Printf("Beta %zu = %24.16E\n", N + 1, beta_[N]);
        }
    }
}
void CGRSolver::update_p() {
    for (size_t N = 0; N < b_.size(); ++N) {
        if (r_converged_[N]) continue;
        p_[N]->scale(beta_[N]);
        p_[N]->add(*z_[N]);
    }

    if (debug_) {
        outfile->Printf("  > Update p <\n\n");
        for (size_t N = 0; N < p_.size(); N++) {
            p_[N]->print();
        }
    }
}

DLUSolver::DLUSolver(std::shared_ptr<UHamiltonian> H)
    : USolver(H),
      nroot_(1),
      norm_(1.E-6),
      max_subspace_(6),
      min_subspace_(2),
      nguess_(1),
      nsubspace_(0),
      nconverged_(0) {
    name_ = "DLU";
}

DLUSolver::~DLUSolver(){};

std::shared_ptr<DLUSolver> DLUSolver::build_solver(Options& options, std::shared_ptr<UHamiltonian> H) {
    auto solver = std::make_shared<DLUSolver>(H);

    if (options["PRINT"].has_changed()) {
        solver->set_print(options.get_int("PRINT") + 1);
    }
    if (options["DEBUG"].has_changed()) {
        solver->set_debug(options.get_int("DEBUG"));
    }
    if (options["BENCH"].has_changed()) {
        solver->set_bench(options.get_int("BENCH"));
    }
    if (options["SOLVER_MAXITER"].has_changed()) {
        solver->set_maxiter(options.get_int("SOLVER_MAXITER"));
    }
    if (options["SOLVER_CONVERGENCE"].has_changed()) {
        solver->set_convergence(options.get_double("SOLVER_CONVERGENCE"));
    }
    if (options["SOLVER_N_ROOT"].has_changed()) {
        solver->set_nroot(options.get_int("SOLVER_N_ROOT"));
    }
    if (options["SOLVER_N_GUESS"].has_changed()) {
        solver->set_nguess(options.get_int("SOLVER_N_GUESS"));
    } else {
        solver->set_nguess(3);
    }
    if (options["SOLVER_MIN_SUBSPACE"].has_changed()) {
        solver->set_min_subspace(options.get_int("SOLVER_MIN_SUBSPACE"));
    }
    if (options["SOLVER_MAX_SUBSPACE"].has_changed()) {
        solver->set_max_subspace(options.get_int("SOLVER_MAX_SUBSPACE"));
    } else {
        solver->set_max_subspace(12);
    }
    if (options["SOLVER_NORM"].has_changed()) {
        solver->set_norm(options.get_double("SOLVER_NORM"));
    }
    if (options["SOLVER_PRECONDITION"].has_changed()) {
        solver->set_precondition(options.get_str("SOLVER_PRECONDITION"));
    }

    return solver;
}

void DLUSolver::print_header() const {
    if (print_) {
        outfile->Printf("  ==> DLUSolver (by Jerome Gonthier) <== \n");
        outfile->Printf(" ==> (Basically a copy-paste from R.Parrish DLR Solver) <== \n\n");
        outfile->Printf("   Number of roots         = %11d\n", nroot_);
        outfile->Printf("   Number of guess vectors = %11d\n", nguess_);
        outfile->Printf("   Maximum subspace size   = %11d\n", max_subspace_);
        outfile->Printf("   Minimum subspace size   = %11d\n", min_subspace_);
        outfile->Printf("   Subspace expansion norm = %11.0E\n", norm_);
        outfile->Printf("   Convergence cutoff      = %11.0E\n", criteria_);
        outfile->Printf("   Maximum iterations      = %11d\n", maxiter_);
        outfile->Printf("   Preconditioning         = %11s\n\n", precondition_.c_str());
    }
}

/*// Commented implementation below is from DLR solver but this function
// is never called in the DLU solver.
size_t DLUSolver::memory_estimate()
{
    size_t dimension = 0L;
    if (!diag_) diag_ = H_->diagonal();
    for (int h = 0; h < diag_->nirrep(); h++) {
        dimension += diag_->dimpi()[h];
    }
    return (2L * max_subspace_ + 3L * nroot_ + 1L) * dimension;
}
*/

void DLUSolver::initialize() {
    finalize();

    c_.clear();
    E_.clear();

    // Diagonal is given to us as an alpha,beta pair
    // that we have to concatenate.

    diag_components = H_->diagonal();

    diag_ = contract_pair(diag_components);

    // We get the dimension of the smallest irrep

    int nirrep = diag_->nirrep();
    const Dimension& dim = diag_->dimpi();
    int mindim = dim[0];

    for (int symm = 1; symm < nirrep; ++symm) {
        if (dim[symm] < mindim) mindim = dim[symm];
    }

    // The maximum subspace dimension should never be larger than
    // the dimension of the smallest irrep. Otherwise, cols. and lines
    // of zeroes are added to the matrix to diagonalize.

    int new_sub = mindim - nroot_;

    if (max_subspace_ > new_sub) {
        outfile->Printf("  SOLVER_MAX_SUBSPACE should not be larger than the dimension \n");
        outfile->Printf("  of the smallest irrep - SOLVER_N_ROOT.\n");
        outfile->Printf("  Setting SOLVER_MAX_SUBSPACE to %4i.\n\n", new_sub);
        max_subspace_ = new_sub;
    }
}

// Contract an alpha/beta pair into a new vector. Each irrep is separately contracted.
std::shared_ptr<Vector> DLUSolver::contract_pair(
    std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector>> components) {
    int nirrepa = components.first->nirrep();
    int nirrepb = components.second->nirrep();
    if (nirrepa != nirrepb) {
        throw PSIEXCEPTION("Alpha and Beta should have same number of irreps.\n");
    }

    const Dimension& dima = components.first->dimpi();
    const Dimension& dimb = components.second->dimpi();
    Dimension dims(nirrepa);
    for (int symm = 0; symm < nirrepa; ++symm) {
        dims[symm] = dima[symm] + dimb[symm];
    }

    auto vec = std::make_shared<Vector>("UStab Alpha + Beta", dims);

    double val = 0;
    for (int symm = 0; symm < nirrepa; ++symm) {
        for (int i = 0; i < dima[symm]; ++i) {
            val = components.first->get(symm, i);
            vec->set(symm, i, val);
        }
        for (int i = 0; i < dimb[symm]; ++i) {
            val = components.second->get(symm, i);
            vec->set(symm, dima[symm] + i, val);
        }
    }

    return vec;
}

// Contract an alpha/beta pair into a result vector. It must have the dimension of
// alpah + beta in each irrep.
void DLUSolver::contract_pair(std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector>> components,
                              std::shared_ptr<Vector> result) {
    int nirrepa = components.first->nirrep();
    int nirrepb = components.second->nirrep();

    if (nirrepa != nirrepb) {
        throw PSIEXCEPTION("Alpha and Beta should have same number of irreps.\n");
    }

    const Dimension& dima = components.first->dimpi();
    const Dimension& dimb = components.second->dimpi();
    const Dimension& dims = result->dimpi();
    for (int symm = 0; symm < nirrepa; ++symm) {
        if (dims[symm] != dima[symm] + dimb[symm]) {
            throw PSIEXCEPTION("Result vector dimpi should be the sum of alpha and beta.\n");
        }
    }

    double val = 0;
    for (int symm = 0; symm < nirrepa; ++symm) {
        for (int i = 0; i < dima[symm]; ++i) {
            val = components.first->get(symm, i);
            result->set(symm, i, val);
        }
        for (int i = 0; i < dimb[symm]; ++i) {
            val = components.second->get(symm, i);
            result->set(symm, dima[symm] + i, val);
        }
    }
}

// Expand a vector into a pair of alpha/beta, created in the routine.
std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector>> DLUSolver::expand_pair(std::shared_ptr<Vector> vec) {
    int nirrepa = diag_components.first->nirrep();
    int nirrepb = diag_components.second->nirrep();
    int nirrep = vec->nirrep();

    if ((nirrep != nirrepa) || (nirrep != nirrepb)) {
        throw PSIEXCEPTION("Full vector irrep does not correspond to alpha or beta.\n");
    }

    const Dimension& dima = diag_components.first->dimpi();
    const Dimension& dimb = diag_components.second->dimpi();
    const Dimension& dims = vec->dimpi();

    for (int symm = 0; symm < nirrep; ++symm) {
        if (dims[symm] != dima[symm] + dimb[symm]) {
            throw PSIEXCEPTION("Wrong irrep dimension of input vector.\n");
        }
    }

    auto pairalpha = std::make_shared<Vector>("UStab Alpha", dima);
    auto pairbeta = std::make_shared<Vector>("UStab Beta", dimb);

    double val = 0;
    for (int symm = 0; symm < nirrep; ++symm) {
        for (int i = 0; i < dima[symm]; ++i) {
            val = vec->get(symm, i);
            pairalpha->set(symm, i, val);
        }
        for (int i = 0; i < dimb[symm]; ++i) {
            val = vec->get(symm, dima[symm] + i);
            pairbeta->set(symm, i, val);
        }
    }

    return std::make_pair(pairalpha, pairbeta);
}

// Expand a vector into an alpha/beta pair result, whose dimensions must match the sum
// of alpha and beta in the input vector for each irrep.
void DLUSolver::expand_pair(std::shared_ptr<Vector> vec,
                            std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector>> result) {
    int nirrepa = result.first->nirrep();
    int nirrepb = result.second->nirrep();
    int nirrep = vec->nirrep();

    if ((nirrep != nirrepa) || (nirrep != nirrepb)) {
        throw PSIEXCEPTION("Full vector irrep does not correspond to alpha or beta.\n");
    }

    const Dimension& dima = result.first->dimpi();
    const Dimension& dimb = result.second->dimpi();
    const Dimension& dims = vec->dimpi();

    for (int symm = 0; symm < nirrep; ++symm) {
        if (dims[symm] != dima[symm] + dimb[symm]) {
            throw PSIEXCEPTION("Wrong irrep dimension of input vector.\n");
        }
    }

    double val = 0;
    for (int symm = 0; symm < nirrep; ++symm) {
        for (int i = 0; i < dima[symm]; ++i) {
            val = vec->get(symm, i);
            result.first->set(symm, i, val);
        }
        for (int i = 0; i < dimb[symm]; ++i) {
            val = vec->get(symm, dima[symm] + i);
            result.second->set(symm, i, val);
        }
    }
}

void DLUSolver::solve() {
    iteration_ = 0;
    converged_ = false;
    nconverged_ = 0;
    convergence_ = 0.0;

    if (print_ > 1) {
        outfile->Printf("  => Iterations <=\n\n");
        outfile->Printf("  %10s %4s %10s %10s %11s\n", "", "Iter", "Converged", "Subspace", "Residual");
    }

    // Compute the first set of sigma vectors
    guess();
    sigma();

    do {
        iteration_++;

        // Compute subspace Hamiltonian
        subspaceHamiltonian();
        // Diagonalize subspace Hamiltonian
        subspaceDiagonalization();
        // Find eigenvectors
        eigenvecs();
        // Find eigenvalues
        eigenvals();
        // Find residuals, update convergence
        residuals();

        if (print_) {
            outfile->Printf("  %-10s %4d %10d %10d %11.3E\n", name_.c_str(), iteration_, nconverged_, nsubspace_,
                            convergence_);
        }

        // Check for convergence
        if (converged_ || iteration_ >= maxiter_) break;

        // Find delta correctors
        correctors();
        // Collapse subspace if needed
        subspaceCollapse();
        // Orthogonalize/add significant correctors
        subspaceExpansion();
        // Compute new sigma vectors
        sigma();

    } while (true);

    if (print_ > 1) {
        outfile->Printf("\n");
        if (!converged_) {
            outfile->Printf("    %sSolver did not converge.\n\n", name_.c_str());
        } else {
            outfile->Printf("    %sSolver converged.\n\n", name_.c_str());
        }
    }
}

void DLUSolver::finalize() {
    b_.clear();
    s_.clear();
    G_.reset();
    a_.reset();
    l_.reset();
    r_.clear();
    n_.clear();
    d_.clear();
    diag_.reset();
}

void DLUSolver::guess() {
    // Find the nguess strongest diagonals in the A matrix
    Dimension rank(diag_->nirrep());
    A_inds_.clear();
    A_inds_.resize(diag_->nirrep());
    for (int h = 0; h < diag_->nirrep(); ++h) {
        int n = diag_->dimpi()[h];
        if (!n) continue;

        std::vector<std::pair<double, int>> d;
        for (int i = 0; i < n; ++i) {
            d.push_back(std::make_pair(diag_->get(h, i), i));
        }
        std::sort(d.begin(), d.end());

        int r = 0;
        for (int i = 0; (i < nguess_) && (i < n); ++i) {
            A_inds_[h].push_back(d[i].second);
            r++;
        }
        rank[h] = r;
    }

    // Preconditioner submatrix and Guess Hamiltonian
    A_ = std::make_shared<Matrix>("A_IJ (Preconditioner)", rank, rank);
    for (int i = 0; i < nguess_; i += max_subspace_) {
        b_.clear();
        s_.clear();
        int n = (max_subspace_ > (nguess_ - i) ? (nguess_ - i) : max_subspace_);
        for (int j = 0; j < n; j++) {
            size_t k = i + j;
            b_.push_back(std::make_shared<Vector>("Delta Guess", diag_->dimpi()));
            for (int h = 0; h < diag_->nirrep(); h++) {
                if (k >= A_inds_[h].size()) continue;
                b_[j]->set(h, A_inds_[h][k], 1.0);
            }
        }

        // Low rank!
        sigma();

        for (int j = 0; j < n; j++) {
            size_t k = i + j;
            for (int h = 0; h < diag_->nirrep(); h++) {
                if (k >= A_inds_[h].size()) continue;
                double** Ap = A_->pointer(h);
                double* sp = s_[j]->pointer(h);
                for (int l = 0; l < rank[h]; l++) {
                    Ap[k][l] = sp[A_inds_[h][l]];
                }
            }
        }
    }

    s_.clear();
    b_.clear();

    SharedMatrix A2(A_->clone());
    auto U = std::make_shared<Matrix>("U", rank, rank);
    auto L = std::make_shared<Vector>("L", rank);

    A2->diagonalize(U, L);

    for (int i = 0; i < nroot_; i++) {
        std::stringstream ss;
        ss << "Guess " << i;
        b_.push_back(std::make_shared<Vector>(ss.str(), diag_->dimpi()));
        for (int h = 0; h < diag_->nirrep(); h++) {
            double** Up = U->pointer(h);
            double* bp = b_[i]->pointer(h);
            for (int j = 0; j < rank[h]; j++) {
                bp[A_inds_[h][j]] = Up[j][i];
            }
        }
    }

    nsubspace_ = nroot_;

    if (debug_) {
        outfile->Printf("   > Guess <\n\n");
        diag_->print();
        A_->print();
        for (size_t i = 0; i < b_.size(); ++i) {
            b_[i]->print();
        }
    }
}
void DLUSolver::sigma() {
    int n = b_.size() - s_.size();
    int offset = s_.size();
    for (int i = 0; i < n; i++) {
        std::stringstream s;
        s << "Sigma Vector " << (i + offset);
        s_.push_back(std::make_shared<Vector>(s.str(), diag_->dimpi()));
    }

    std::vector<std::shared_ptr<Vector>> x;
    std::vector<std::shared_ptr<Vector>> b;

    for (int i = offset; i < offset + n; i++) {
        x.push_back(b_[i]);
        b.push_back(s_[i]);
    }

    std::vector<std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector>>> xpair;
    std::vector<std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector>>> bpair;

    // Get the big concatenated alpha/beta vectors into individual pair components.

    for (int i = 0; i < n; i++) {
        xpair.push_back(expand_pair(x[i]));
        bpair.push_back(expand_pair(b[i]));
    }

    H_->product(xpair, bpair);

    // Concatenate the pair back into single vectors

    for (int i = 0; i < n; i++) {
        contract_pair(xpair[i], x[i]);
        contract_pair(bpair[i], b[i]);
    }

    if (debug_) {
        outfile->Printf("   > Sigma <\n\n");
        for (size_t i = 0; i < s_.size(); i++) {
            s_[i]->print();
        }
    }
}
void DLUSolver::subspaceHamiltonian() {
    int n = s_.size();
    int nirrep = diag_->nirrep();
    int* npi = new int[nirrep];
    for (int h = 0; h < nirrep; ++h) {
        npi[h] = n;
    }

    G_ = std::make_shared<Matrix>("Subspace Hamiltonian", nirrep, npi, npi);
    delete[] npi;

    for (int h = 0; h < nirrep; ++h) {
        int dimension = diag_->dimpi()[h];

        if (!dimension) continue;

        double** Gp = G_->pointer(h);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                Gp[i][j] = Gp[j][i] = C_DDOT(dimension, b_[i]->pointer(h), 1, s_[j]->pointer(h), 1);
            }
        }
    }

    if (debug_) {
        outfile->Printf("   > SubspaceHamiltonian <\n\n");
        G_->print();
    }
}
void DLUSolver::subspaceDiagonalization() {
    int n = s_.size();
    int nirrep = diag_->nirrep();
    Dimension npi(nirrep);
    for (int h = 0; h < nirrep; ++h) {
        npi[h] = n;
    }

    SharedMatrix G2(G_->clone());
    a_ = std::make_shared<Matrix>("Subspace Eigenvectors", npi, npi);
    l_ = std::make_shared<Vector>("Subspace Eigenvalues", npi);

    G2->diagonalize(a_, l_);

    // Resort to remove false zeros for cases with too small of irreps
    for (int h = 0; h < nirrep; h++) {
        int dim = diag_->dimpi()[h];

        int nfalse = n - dim;

        if (nfalse <= 0) continue;

        double** ap = a_->pointer(h);
        double* lp = l_->pointer(h);

        for (int m = 0; m < n - nfalse; m++) {
            lp[m] = lp[m + nfalse];
            C_DCOPY(n, &ap[0][m + nfalse], n, &ap[0][m], n);
        }

        for (int m = 0; m < nfalse; m++) {
            lp[n - m - 1] = 0;
            C_DSCAL(n, 0.0, &ap[0][n - m - 1], n);
        }
    }

    if (debug_) {
        outfile->Printf("   > SubspaceDiagonalize <\n\n");
        a_->print();
        l_->print();
    }
}

void DLUSolver::eigenvecs() {
    if (c_.size() != (size_t)nroot_) {
        c_.clear();
        for (int m = 0; m < nroot_; ++m) {
            std::stringstream s;
            s << "Eigenvector " << m;
            auto c = std::make_shared<Vector>(s.str().c_str(), diag_->dimpi());
            c_.push_back(c);
        }
    }

    for (int h = 0; h < diag_->nirrep(); ++h) {
        int dimension = diag_->dimpi()[h];

        if (!dimension) continue;

        double** ap = a_->pointer(h);
        for (int m = 0; m < nroot_; m++) {
            double* cp = c_[m]->pointer(h);
            ::memset((void*)cp, '\0', dimension * sizeof(double));
            for (size_t i = 0; i < b_.size(); i++) {
                double* bp = b_[i]->pointer(h);
                C_DAXPY(dimension, ap[i][m], bp, 1, cp, 1);
            }
        }
    }

    if (debug_) {
        outfile->Printf("   > Eigenvectors <\n\n");
        for (size_t m = 0; m < c_.size(); m++) {
            c_[m]->print();
        }
    }
}

void DLUSolver::eigenvals() {
    E_.clear();
    E_.resize(nroot_);

    for (int h = 0; h < diag_->nirrep(); ++h) {
        for (int k = 0; k < nroot_; k++) {
            E_[k].push_back(l_->get(h, k));
        }
    }

    if (debug_) {
        outfile->Printf("   > Eigenvalues <\n\n");
        for (size_t m = 0; m < E_.size(); m++) {
            for (size_t h = 0; h < E_[0].size(); ++h) {
                outfile->Printf("    Eigenvalue %zu, Irrep %d = %24.16E\n", m, h, E_[m][h]);
            }
        }
        outfile->Printf("\n");
    }
}

void DLUSolver::residuals() {
    n_.resize(nroot_);
    nconverged_ = 0;

    if (r_.size() != (size_t)nroot_) {
        r_.clear();
        for (int k = 0; k < nroot_; ++k) {
            // Residual k
            std::stringstream s;
            s << "Residual Vector " << k;
            r_.push_back(std::make_shared<Vector>(s.str().c_str(), diag_->dimpi()));
        }
    }

    for (int k = 0; k < nroot_; k++) {
        double R2 = 0.0;
        double S2 = 0.0;

        for (int h = 0; h < diag_->nirrep(); ++h) {
            int dimension = diag_->dimpi()[h];
            if (!dimension) continue;

            double** ap = a_->pointer(h);
            double* lp = l_->pointer(h);
            double* rp = r_[k]->pointer(h);
            double* cp = c_[k]->pointer(h);

            ::memset((void*)rp, '\0', dimension * sizeof(double));

            for (size_t i = 0; i < b_.size(); i++) {
                double* sp = s_[i]->pointer(h);
                C_DAXPY(dimension, ap[i][k], sp, 1, rp, 1);
            }
            S2 += C_DDOT(dimension, rp, 1, rp, 1);

            C_DAXPY(dimension, -lp[k], cp, 1, rp, 1);

            R2 += C_DDOT(dimension, rp, 1, rp, 1);
        }

        // Residual norm k
        double rnorm = sqrt(R2 / S2);
        n_[k] = rnorm;
        if (rnorm < criteria_) {
            nconverged_++;
        }
    }

    // Global convergence check
    convergence_ = 0.0;
    for (int k = 0; k < nroot_; k++) {
        if (convergence_ < n_[k]) convergence_ = n_[k];
    }

    if (nconverged_ == nroot_) converged_ = true;
    if (debug_) {
        outfile->Printf("   > Residuals <\n\n");
        for (size_t i = 0; i < r_.size(); i++) {
            r_[i]->print();
        }
        for (size_t i = 0; i < n_.size(); i++) {
            outfile->Printf("    Residual %zu = %24.16E\n", i, n_[i]);
        }
        outfile->Printf("\n");
        outfile->Printf("    %d of %d roots converged, we are %s\n\n", nconverged_, nroot_,
                        (converged_ ? "converged" : "not converged"));
    }
}

void DLUSolver::correctors() {
    // Only add correctors for roots that are not converged
    d_.clear();

    for (int k = 0; k < nroot_; k++) {
        // Do not attempt to add a corrector if root is already converged
        if (n_[k] < criteria_) continue;

        std::stringstream s;
        s << "Corrector Vector " << k;
        auto d = std::make_shared<Vector>(s.str().c_str(), diag_->dimpi());

        for (int h = 0; h < diag_->nirrep(); ++h) {
            int dimension = diag_->dimpi()[h];
            if (!dimension) continue;

            double* hp = diag_->pointer(h);
            double lambda = E_[k][h];
            double* dp = d->pointer(h);
            double* rp = r_[k]->pointer(h);

            if (precondition_ == "SUBSPACE") {
                for (int m = 0; m < dimension; m++) {
                    dp[m] = rp[m] / (hp[m] - lambda);
                }

                // Cannot subspace precondition on the first iteration
                if (iteration_ > 1) {
                    int rank = A_inds_[h].size();
                    auto A2 = std::make_shared<Matrix>("A2", rank, rank);
                    double** A2p = A2->pointer();
                    double** Ap = A_->pointer(h);
                    ::memcpy((void*)A2p[0], (void*)Ap[0], sizeof(double) * rank * rank);
                    for (int i = 0; i < rank; i++) {
                        A2p[i][i] -= lambda;
                    }

                    int* ipiv = new int[rank];
                    int info = C_DGETRF(rank, rank, A2p[0], rank, ipiv);
                    // Only apply the improved preconditioner if nonsingular
                    if (!info) {
                        double* v = new double[rank];
                        for (int i = 0; i < rank; i++) {
                            v[i] = rp[A_inds_[h][i]];
                        }
                        C_DGETRS('N', rank, 1, A2p[0], rank, ipiv, v, rank);
                        for (int i = 0; i < rank; i++) {
                            dp[A_inds_[h][i]] = v[i];
                        }
                        delete[] v;
                    }
                    delete[] ipiv;
                }

            } else if (precondition_ == "JACOBI") {
                for (int m = 0; m < dimension; m++) {
                    dp[m] = rp[m] / (lambda - hp[m]);
                }
            } else {
                C_DCOPY(dimension, rp, 1, dp, 1);
            }

            // Substitute r for this vector, if norm is bad
            double norm = sqrt(C_DDOT(dimension, dp, 1, dp, 1));
            if (norm != norm || std::isinf(norm)) {
                C_DCOPY(dimension, rp, 1, dp, 1);
                norm = sqrt(C_DDOT(dimension, dp, 1, dp, 1));
            }

            double scale = 1.0 / norm;
            if (scale != scale || std::isinf(scale)) {
                scale = 0.0;
            }

            // Normalize the correctors
            C_DSCAL(dimension, scale, dp, 1);
        }

        d_.push_back(d);
    }
    if (debug_) {
        outfile->Printf("   > Correctors <\n\n");
        for (size_t i = 0; i < d_.size(); i++) {
            d_[i]->print();
        }
    }
}

void DLUSolver::subspaceExpansion() {
    if (debug_) {
        outfile->Printf("   > SubspaceExpansion <\n\n");
    }

    // Which vectors are significant?
    std::vector<bool> sig(d_.size());
    for (size_t i = 0; i < d_.size(); ++i) {
        sig[i] = false;
    }

    // Orthonormalize d_ via Modified Gram-Schmidt
    for (int h = 0; h < diag_->nirrep(); ++h) {
        int dimension = diag_->dimpi()[h];
        if (!dimension) continue;

        // Remove the projection of d on b from b
        for (size_t i = 0; i < d_.size(); ++i) {
            for (size_t j = 0; j < b_.size(); ++j) {
                double* dp = d_[i]->pointer(h);
                double* bp = b_[j]->pointer(h);

                double r_ji = C_DDOT(dimension, dp, 1, bp, 1);
                C_DAXPY(dimension, -r_ji, bp, 1, dp, 1);
            }
        }

        // Remove the self-projection of d on d from d
        for (size_t i = 0; i < d_.size(); ++i) {
            double* dip = d_[i]->pointer(h);
            double r_ii = sqrt(C_DDOT(dimension, dip, 1, dip, 1));
            C_DSCAL(dimension, (r_ii > norm_ ? 1.0 / r_ii : 0.0), dip, 1);
            for (size_t j = i + 1; j < d_.size(); ++j) {
                double* djp = d_[j]->pointer(h);
                double r_ij = C_DDOT(dimension, djp, 1, dip, 1);
                C_DAXPY(dimension, -r_ij, dip, 1, djp, 1);
            }
            if (r_ii > norm_) {
                sig[i] = sig[i] | true;
            }
        }
    }

    // Add significant vectors
    for (size_t i = 0; i < d_.size(); ++i) {
        if (sig[i]) {
            b_.push_back(d_[i]);
        }
    }

    nsubspace_ = b_.size();

    if (debug_) {
        outfile->Printf("Final subspace after addition\n\n");
        for (size_t i = 0; i < b_.size(); i++) {
            b_[i]->print();
        }
    }
}

void DLUSolver::subspaceCollapse() {
    if (nsubspace_ <= max_subspace_) return;

    std::vector<std::shared_ptr<Vector>> s2;
    std::vector<std::shared_ptr<Vector>> b2;

    for (int k = 0; k < min_subspace_; ++k) {
        std::stringstream bs;
        bs << "Subspace Vector " << k;
        b2.push_back(std::make_shared<Vector>(bs.str(), diag_->dimpi()));
        std::stringstream ss;
        ss << "Sigma Vector " << k;
        s2.push_back(std::make_shared<Vector>(ss.str(), diag_->dimpi()));
    }

    int n = a_->rowspi()[0];
    for (int k = 0; k < min_subspace_; ++k) {
        for (int h = 0; h < diag_->nirrep(); ++h) {
            int dimension = diag_->dimpi()[h];
            if (!dimension) continue;

            double** ap = a_->pointer(h);
            double* b2p = b2[k]->pointer(h);
            double* s2p = s2[k]->pointer(h);

            for (int i = 0; i < n; ++i) {
                double* bp = b_[i]->pointer(h);
                double* sp = s_[i]->pointer(h);
                C_DAXPY(dimension, ap[i][k], sp, 1, s2p, 1);
                C_DAXPY(dimension, ap[i][k], bp, 1, b2p, 1);
            }
        }
    }

    s_ = s2;
    b_ = b2;
    nsubspace_ = b_.size();

    if (debug_) {
        outfile->Printf("   > SubspaceCollapse <\n\n");
        for (size_t i = 0; i < b_.size(); i++) {
            b_[i]->print();
        }
        for (size_t i = 0; i < s_.size(); i++) {
            s_[i]->print();
        }
    }
}
}  // namespace psi
