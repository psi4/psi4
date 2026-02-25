/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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
            outfile->Printf(
                "  The subspace preconditioner has been broken for some time and was removed in Psi4 1.4.\n");
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
}  // namespace psi
