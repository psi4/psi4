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

#ifndef SOLVER_H
#define SOLVER_H

#include <psi4/libmints/typedefs.h>
#include <psi4/liboptions/liboptions.h>

#include <vector>
#include <string>

namespace psi {

class Vector;
class RHamiltonian;
class UHamiltonian;

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

class USolver : public Solver {
   protected:
    /// Reference to underlying UHamiltonian
    std::shared_ptr<UHamiltonian> H_;

   public:
    // => Constructors < = //

    /// Reference to underlying RHamiltonian
    USolver(std::shared_ptr<UHamiltonian> H);
    /// Destructor
    ~USolver() override;

    // => Accessors <= //

    /**
    * Knob to swap out a Hamiltonian object
    * @param H new UHamiltonian object
    */
    void set_Hamiltonian(std::shared_ptr<UHamiltonian> H) { H_ = H; }

    /**
    * Pointer to the Hamiltonian object
    * @return current UHamiltonian object
    */
    std::shared_ptr<UHamiltonian> H() const { return H_; }
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

// Class for solving UHF stability analysis.
class DLUSolver : public USolver {
   protected:
    // => Control parameters <= //

    /// Number of desired roots
    int nroot_;
    /// Required norm for subspace expansion
    double norm_;
    /// Maximum allowed subspace size
    int max_subspace_;
    /// Minimum allowed subspace size (after collapse)
    int min_subspace_;
    /// Number of guess vectors to build
    int nguess_;

    // => Iteration values <= //

    /// The current subspace size
    int nsubspace_;
    /// The number of converged roots
    int nconverged_;

    // => State values <= //

    /// Current eigenvectors (nroots)
    std::vector<std::shared_ptr<Vector> > c_;
    /// Current eigenvalues (nroots)
    std::vector<std::vector<double> > E_;
    /// B vectors (nsubspace)
    std::vector<std::shared_ptr<Vector> > b_;
    /// Sigma vectors (nsubspace)
    std::vector<std::shared_ptr<Vector> > s_;
    /// Delta Subspace Hamiltonian (preconditioner)
    SharedMatrix A_;
    /// Delta Subspace indices
    std::vector<std::vector<int> > A_inds_;
    /// G_ij Subspace Hamiltonian (nsubspace x nsubspace)
    SharedMatrix G_;
    /// Subspace eigenvectors (nsubspace x nsubspace)
    SharedMatrix a_;
    /// Subspace eigenvalues (nsubspace)
    std::shared_ptr<Vector> l_;
    /// Residual vectors (nroots)
    std::vector<std::shared_ptr<Vector> > r_;
    /// Residual vector 2-norms (nroots)
    std::vector<double> n_;
    /// Correction vectors (nroots)
    std::vector<std::shared_ptr<Vector> > d_;
    /// Diagonal of Hamiltonian
    std::shared_ptr<Vector> diag_;
    /// Diagonal components of UHamiltonian
    std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> > diag_components;

    // => Run routines <= //

    // Guess, based on diagonal
    void guess();
    // Compute sigma vectors for the given set of b
    void sigma();
    // Compute subspace Hamiltonian
    void subspaceHamiltonian();
    // Diagonalize subspace Hamiltonian
    void subspaceDiagonalization();
    // Find eigenvectors
    void eigenvecs();
    // Find eigenvalues
    void eigenvals();
    // Find residuals, update convergence
    void residuals();
    // Find correctors
    virtual void correctors();
    // Orthogonalize/add significant correctors
    void subspaceExpansion();
    // Collapse subspace if needed
    void subspaceCollapse();

   public:
    // => Constructors <= //

    /// Constructor
    DLUSolver(std::shared_ptr<UHamiltonian> H);
    /// Destructor
    ~DLUSolver() override;

    /// Static constructor, uses Options object
    static std::shared_ptr<DLUSolver> build_solver(Options& options, std::shared_ptr<UHamiltonian> H);

    // => Required Methods <= //

    void print_header() const override;
    size_t memory_estimate() override { return 0; };
    void initialize() override;
    void solve() override;
    void finalize() override;

    // => Accessors <= //

    /// Eigenvectors, by state
    const std::vector<std::shared_ptr<Vector> >& eigenvectors() const { return c_; }
    /// Eigenvalues, by state/irrep
    const std::vector<std::vector<double> >& eigenvalues() const { return E_; }
    /// Convert an alpha/beta pair into a single vector
    std::shared_ptr<Vector> contract_pair(std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> > components);
    void contract_pair(std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> > components,
                       std::shared_ptr<Vector> result);
    /// Convert a single vector into an alpha/beta pair
    std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> > expand_pair(std::shared_ptr<Vector> vec);
    void expand_pair(std::shared_ptr<Vector> vec, std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> > result);

    // => Knobs <= //

    /// Set number of roots (defaults to 1)
    void set_nroot(int nroot) { nroot_ = nroot; }
    /// Set maximum subspace size (defaults to 6)
    void set_max_subspace(double max_subspace) { max_subspace_ = max_subspace; }
    /// Set minimum subspace size, for collapse (defaults to 2)
    void set_min_subspace(double min_subspace) { min_subspace_ = min_subspace; }
    /// Set number of guesses (defaults to 1)
    void set_nguess(int nguess) { nguess_ = nguess; }
    /// Set norm critera for adding vectors to subspace (defaults to 1.0E-6)
    void set_norm(double norm) { norm_ = norm; }
};
}
#endif
