/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#ifndef _psi_src_lib_libmints_local_h_
#define _psi_src_lib_libmints_local_h_

#include <string>
#include <vector>
#include <memory>

#include "psi4/pragma.h"

namespace psi {
class Options;

class Matrix;

class BasisSet;

class PSI_API Localizer {
   protected:
    // => Parameters <= //

    /// Print flag
    int print_;
    /// Debug flug
    int debug_;
    /// Bench flag
    int bench_;

    /// Relative convergence criteria
    double convergence_;
    /// Maximum orbital-gradient convergence criterion
    double gradient_convergence_;
    /// Maximum number of iterations
    int maxiter_;
    /// Use augmented-Hessian steps after the Jacobi startup sweeps
    bool use_augmented_hessian_;
    /// First iteration on which an augmented-Hessian step may be attempted
    int augmented_hessian_start_;
    /// Largest number of independent rotations for which the dense Hessian is formed
    int augmented_hessian_max_rotations_;
    /// Maximum Davidson subspace dimension for matrix-free augmented-Hessian solves
    int augmented_hessian_max_subspace_;
    /// Euclidean trust radius for an augmented-Hessian rotation step
    double augmented_hessian_trust_radius_;
    /// Positive-curvature threshold used to identify localization saddles
    double saddle_tolerance_;

    /// Primary orbital basis set
    std::shared_ptr<BasisSet> primary_;
    /// Delocalized Orbitals
    std::shared_ptr<Matrix> C_;

    // => Targets <= //

    /// Localized Orbitals
    std::shared_ptr<Matrix> L_;
    /// MO -> LO transformation
    std::shared_ptr<Matrix> U_;
    /// Did the algorithm converge?
    bool converged_;

    /// Set defaults
    void common_init();
    /// Apply common localization controls from an Options object
    void configure(Options &options);

    /// Maximize sum_K sum_p [(U^T A_K U)_pp]^power for symmetric orbital-space operators A_K.
    /// This is the common form of the Boys, generalized Pipek--Mezey, and IBO objectives.
    /// Optional range boundaries restrict rotations to independent orbital subspaces.
    void localize_matrix_objective(std::vector<std::shared_ptr<Matrix>> operators, const std::string &label,
                                   int power = 2, const std::vector<int> &ranges = {});

   public:
    // => Constructors <= //

    Localizer(std::shared_ptr<BasisSet> primary_, std::shared_ptr<Matrix> C);

    virtual ~Localizer();

    static std::shared_ptr<Localizer> build(const std::string &type, std::shared_ptr<BasisSet> primary,
                                            std::shared_ptr<Matrix> C);

    static std::shared_ptr<Localizer> build(const std::string &type, std::shared_ptr<BasisSet> primary,
                                            std::shared_ptr<Matrix> C, Options &options);

    static std::shared_ptr<Localizer> build(std::shared_ptr<BasisSet> primary, std::shared_ptr<Matrix> C,
                                            Options &options);

    // => Computers <= //

    /// Print out the localization algorithm and parameters
    virtual void print_header() const = 0;

    /// Perform the localization algorithm
    virtual void localize() = 0;

    /// Given a Fock matrix in the original basis (usually diagonal), produce an ordered copy in the local basis, and
    /// reorder L and U
    std::shared_ptr<Matrix> fock_update(std::shared_ptr<Matrix> F_orig,
                                        const std::vector<int> &ranges = {});

    // => Accessors <= //

    std::shared_ptr<Matrix> L() const { return L_; }

    std::shared_ptr<Matrix> U() const { return U_; }

    bool converged() const { return converged_; }

    // => Knobs <= //

    void set_print(int print) { print_ = print; }

    void set_debug(int debug) { debug_ = debug; }

    void set_bench(int bench) { bench_ = bench; }

    void set_convergence(double convergence) { convergence_ = convergence; }

    void set_gradient_convergence(double convergence) { gradient_convergence_ = convergence; }

    void set_maxiter(int maxiter) { maxiter_ = maxiter; }

    void set_use_augmented_hessian(bool enabled) { use_augmented_hessian_ = enabled; }

    void set_augmented_hessian_start(int iteration) { augmented_hessian_start_ = iteration; }

    void set_augmented_hessian_max_rotations(int rotations) { augmented_hessian_max_rotations_ = rotations; }

    void set_augmented_hessian_max_subspace(int vectors) { augmented_hessian_max_subspace_ = vectors; }

    void set_augmented_hessian_trust_radius(double radius) { augmented_hessian_trust_radius_ = radius; }

    void set_saddle_tolerance(double tolerance) { saddle_tolerance_ = tolerance; }
};

class PSI_API BoysLocalizer : public Localizer {
   protected:
    /// Set defaults
    void common_init();

   public:
    BoysLocalizer(std::shared_ptr<BasisSet> primary, std::shared_ptr<Matrix> C);

    ~BoysLocalizer() override;

    void print_header() const override;

    void localize() override;
};

class PSI_API PMLocalizer : public Localizer {
   protected:
    /// Set defaults
    void common_init();
    /// Symmetric atomic population operators in the input-orbital basis
    std::vector<std::shared_ptr<Matrix>> population_matrices_;
    /// Label for the population partition used by the generalized PM functional
    std::string population_method_;
    /// Even power used in the generalized PM objective (normally 2; IBO conventionally uses 4)
    int power_;

   public:
    PMLocalizer(std::shared_ptr<BasisSet> primary, std::shared_ptr<Matrix> C);

    PMLocalizer(std::shared_ptr<BasisSet> primary, std::shared_ptr<Matrix> C,
                const std::vector<std::shared_ptr<Matrix>> &population_matrices,
                const std::string &population_method = "USER");

    ~PMLocalizer() override;

    void print_header() const override;

    void localize() override;

    void set_power(int power) { power_ = power; }
};

/**
 * Intrinsic-bond-orbital localizer.
 *
 * This libmints implementation deliberately coexists with the historical
 * FISAPT IBOLocalizer2 while that client is migrated and regression-tested.
 * It constructs Knizia intrinsic atomic orbitals (IAOs), converts their
 * atom-resolved projectors into generalized Pipek--Mezey population
 * operators, and delegates the orbital optimization (including augmented-
 * Hessian saddle control) to Localizer::localize_matrix_objective.
 */
class PSI_API IBOLocalizer : public PMLocalizer {
   protected:
    /// Minimal basis used to construct IAOs
    std::shared_ptr<BasisSet> minao_;
    /// Complete occupied space defining the determinant projector. This may
    /// include frozen occupied orbitals not present in C_.
    std::shared_ptr<Matrix> C_reference_;
    /// Include basis functions on ghost centers in the IAO partition
    bool use_ghosts_;
    /// Eigenvalue cutoff used in IAO metric inverse square roots
    double condition_;
    /// Boundaries of independently localized orbital blocks
    std::vector<int> ranges_;

    /// Maps compact, non-ghosted atoms/IAOs to the complete molecule/minimal basis
    std::vector<int> true_atoms_;
    std::vector<int> true_iaos_;
    std::vector<int> iaos_to_atoms_;

    /// AO overlap, orthonormal IAOs in the primary basis, and final orbital populations
    std::shared_ptr<Matrix> S_;
    std::shared_ptr<Matrix> A_;
    std::shared_ptr<Matrix> Q_;

    void common_init();
    void build_iaos();
    void build_population_matrices();
    void update_orbital_charges();
    std::shared_ptr<Matrix> orbital_charges(const std::shared_ptr<Matrix> &orbitals) const;

   public:
    IBOLocalizer(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> minao,
                 std::shared_ptr<Matrix> C);
    IBOLocalizer(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> minao,
                 std::shared_ptr<Matrix> C, std::shared_ptr<Matrix> C_reference);
    ~IBOLocalizer() override;

    static std::shared_ptr<IBOLocalizer> build(std::shared_ptr<BasisSet> primary,
                                               std::shared_ptr<BasisSet> minao,
                                               std::shared_ptr<Matrix> C);
    static std::shared_ptr<IBOLocalizer> build(std::shared_ptr<BasisSet> primary,
                                               std::shared_ptr<BasisSet> minao,
                                               std::shared_ptr<Matrix> C,
                                               std::shared_ptr<Matrix> C_reference);
    static std::shared_ptr<IBOLocalizer> build(std::shared_ptr<BasisSet> primary,
                                               std::shared_ptr<BasisSet> minao,
                                               std::shared_ptr<Matrix> C, Options &options);
    static std::shared_ptr<IBOLocalizer> build(std::shared_ptr<BasisSet> primary,
                                               std::shared_ptr<BasisSet> minao,
                                               std::shared_ptr<Matrix> C,
                                               std::shared_ptr<Matrix> C_reference, Options &options);

    void print_header() const override;
    void localize() override;
    void print_charges(double scale = 2.0);

    std::shared_ptr<Matrix> A() const { return A_; }
    std::shared_ptr<Matrix> Q() const;

    void set_use_ghosts(bool use_ghosts) {
        if (use_ghosts_ != use_ghosts) {
            use_ghosts_ = use_ghosts;
            S_.reset();
            A_.reset();
            Q_.reset();
            population_matrices_.clear();
            L_.reset();
            U_.reset();
            converged_ = false;
        }
    }
    void set_condition(double condition) {
        if (condition_ != condition) {
            condition_ = condition;
            S_.reset();
            A_.reset();
            Q_.reset();
            population_matrices_.clear();
            L_.reset();
            U_.reset();
            converged_ = false;
        }
    }
    void set_ranges(const std::vector<int> &ranges) {
        ranges_ = ranges;
        L_.reset();
        U_.reset();
        Q_.reset();
        converged_ = false;
    }
};

class PSI_API ERLocalizer : public Localizer {
   protected:
    /// Set defaults
    void common_init();
    /// Auxiliary basis set
    std::shared_ptr<BasisSet> auxiliary_;
    /// Reusable AO collocation factor, x^I_mu
    std::shared_ptr<Matrix> x_ao_;
    /// Reusable THC coupling factor, Z^IJ
    std::shared_ptr<Matrix> Z_;

   public:
    ERLocalizer(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, std::shared_ptr<Matrix> C);

    /// Construct from a previously computed AO THC factorization.  The factors
    /// are orbital independent and may be shared by successive localizations.
    ERLocalizer(std::shared_ptr<BasisSet> primary, std::shared_ptr<Matrix> C, std::shared_ptr<Matrix> x_ao,
                std::shared_ptr<Matrix> Z);

    ~ERLocalizer() override;

    void print_header() const override;

    void localize() override;
};

}  // Namespace psi

#endif
