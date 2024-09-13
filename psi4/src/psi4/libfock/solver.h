/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#ifndef SOLVER_H
#define SOLVER_H

#include <psi4/libmints/typedefs.h>
#include <psi4/liboptions/liboptions.h>

#include <vector>
#include <string>

namespace psi {

class Vector;
class RHamiltonian;

class Solver {
    // => BASE CLASSES <= //

   protected:
    /// Print flag, defaults to 1
    int print_;
    /// Debug flag, defaults to 0
    int debug_;
    /// Bench flag, defaults to 0
    int bench_;
    /// Name of solver (set by subclasses)
    std::string name_;
    /// Memory available, in doubles, defaults to 0 => Unlimited storage
    size_t memory_;

    /// Convergence criteria, defaults to 1.0E-6
    double criteria_;
    /// Maximum number of iterations, defaults to 100
    int maxiter_;
    /// Converged or not?
    bool converged_;
    /// Convergence measure at this iteration
    double convergence_;
    /// Current iteration count
    int iteration_;
    /// Preconditioner type
    std::string precondition_;

    /// Common initialization
    void common_init();

   public:
    // => Constructors < = //

    /// Default Constructor
    Solver();
    /// Destructor
    virtual ~Solver();

    // => Knobs <= //

    /// Set precondition type (specific to solver type)
    void set_precondition(const std::string& precondition) { precondition_ = precondition; }
    /// Set maximum vector storage space (defaults to 0 MB => Unlimited storage)
    void set_memory(size_t memory) { memory_ = memory; }
    /// Set maximum number of iterations (defaults to 100)
    void set_maxiter(int maxiter) { maxiter_ = maxiter; }
    /// Set convergence criteria (defaults to 1.0E-6)
    void set_convergence(double criteria) { criteria_ = criteria; }
    /// Print flag (defaults to 1)
    void set_print(int print) { print_ = print; }
    /// Debug flag (defaults to 0)
    void set_debug(int debug) { debug_ = debug; }
    /// Bench flag (defaults to 0)
    void set_bench(int bench) { bench_ = bench; }

    // => Accessors <= //

    /// What iteration is it?
    int iteration() const { return iteration_; }
    /// Did the solver converge?
    bool converged() const { return converged_; }
    /**
    * Print header information regarding Solver
    * type on output file
    */
    virtual void print_header() const = 0;
    /// Estimate of maximum memory usage (in doubles)
    virtual size_t memory_estimate() = 0;

    // => Computers <= //

    /**
     * Warm it all up, everything you got.
     * Cmon you apes, you wanna live forever?!
     */
    virtual void initialize() = 0;
    /**
     * Perform iterative solution,
     * based on current state
     */
    virtual void solve() = 0;
    /**
     * Method to clear off overhead memory
     * without destroying the object
     */
    virtual void finalize() = 0;
};

class RSolver : public Solver {
   protected:
    /// Reference to underlying RHamiltonian
    std::shared_ptr<RHamiltonian> H_;

   public:
    // => Constructors < = //

    RSolver(std::shared_ptr<RHamiltonian> H);
    /// Destructor
    ~RSolver() override;

    // => Accessors <= //

    /**
    * Knob to swap out a Hamiltonian object
    * @param H new RHamiltonian object
    */
    void set_Hamiltonian(std::shared_ptr<RHamiltonian> H) { H_ = H; }

    /**
    * Pointer to the Hamiltonian object
    * @return current RHamiltonian object
    */
    std::shared_ptr<RHamiltonian> H() const { return H_; }
};

// => APPLIED CLASSES <= //

// Conjugate gradients with a spin-restricted Hamiltonian.
// Solves the equations Hx=b.
class CGRSolver : public RSolver {
   protected:
    /// Force vectors
    std::vector<std::shared_ptr<Vector> > b_;
    /// Solution vectors
    std::vector<std::shared_ptr<Vector> > x_;
    /// Product vectors
    std::vector<std::shared_ptr<Vector> > Ap_;
    /// M^-1 x
    std::vector<std::shared_ptr<Vector> > z_;
    /// Residual vectors
    std::vector<std::shared_ptr<Vector> > r_;
    /// Conjugate directions
    std::vector<std::shared_ptr<Vector> > p_;
    /// Alpha values
    std::vector<double> alpha_;
    /// Beta values
    std::vector<double> beta_;
    /// Residual norm
    std::vector<double> r_nrm2_;
    /// z'r
    std::vector<double> z_r_;
    /// Which vectors are converged?
    std::vector<bool> r_converged_;
    /// Number of converged vectors
    int nconverged_;

    /// Diagonal M, for guess and Jacobi preconditioning
    std::shared_ptr<Vector> diag_;
    /// A subspace matrix, for preconditioning
    SharedMatrix A_;
    /// A subspace indices
    std::vector<std::vector<int> > A_inds_;
    /// Shifts [to solve (A-mI)]; Outer vector indexes irreps, inner indexes basis vectors of that irrep
    std::vector<std::vector<double> > shifts_;
    /// Number of guess vectors to use for subspace preconditioner
    int nguess_;

    /// Initializes shifts_ to 0
    void setup();
    void guess();
    void residual();
    /// Write (H_-shifts_)*x_ to Ap_
    void products_x();
    void products_p();
    void alpha();
    void update_x();
    void update_r();
    void check_convergence();
    void update_z();
    void beta();
    void update_p();

   public:
    CGRSolver(std::shared_ptr<RHamiltonian> H);
    ~CGRSolver() override;

    /// Static constructor, uses Options object
    static std::shared_ptr<CGRSolver> build_solver(Options& options, std::shared_ptr<RHamiltonian> H);

    std::vector<std::shared_ptr<Vector> >& x() { return x_; }
    std::vector<std::shared_ptr<Vector> >& b() { return b_; }

    void print_header() const override;
    size_t memory_estimate() override;
    void initialize() override;
    void solve() override;
    void finalize() override;

    void set_shifts(const std::vector<std::vector<double> >& shifts) { shifts_ = shifts; }
    void set_A(SharedMatrix A, const std::vector<std::vector<int> > inds) {
        A_ = A;
        A_inds_ = inds;
    }
    void set_nguess(int nguess) { nguess_ = nguess; }
};
}
#endif
