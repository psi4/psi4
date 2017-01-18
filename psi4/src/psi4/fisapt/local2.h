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

#ifndef FISAPT_LOCAL2_H
#define FISAPT_LOCAL2_H

#include <vector>
#include <memory>

namespace psi {

class Matrix;
class BasisSet;

namespace fisapt {

class IBOLocalizer2 {

protected:

    // => Overall Parameters <= //

    /// Print flag
    int print_;
    /// Debug flug
    int debug_;
    /// Bench flag
    int bench_;

    /// Relative convergence criteria
    double convergence_;
    /// Maximum number of iterations
    int maxiter_;

    // => IAO Parameters <= //

    /// Use ghost IAOs?
    bool use_ghosts_;
    /// IAO localization power (4 or 2)
    int power_;
    /// Metric condition for IAO
    double condition_;

    /// Occupied orbitals, in primary basis
    std::shared_ptr<Matrix> C_;
    /// Primary orbital basis set
    std::shared_ptr<BasisSet> primary_;
    /// MinAO orbital baiss set
    std::shared_ptr<BasisSet> minao_;

    // => Stars Parameters <= //

    /// Do stars treatment?
    bool use_stars_;
    /// Charge completeness for two-center orbitals
    double stars_completeness_;
    /// List of centers for stars
    std::vector<int> stars_;

    // => IAO Data <= //

    /// Map from non-ghosted to full atoms: true_atoms[ind_true] = ind_full
    std::vector<int> true_atoms_;
    /// Map from non-ghosted IAOs to full IAOs: true_iaos[ind_true] = ind_full
    std::vector<int> true_iaos_;
    /// Map from non-ghosted IAOs to non-ghosted atoms
    std::vector<int> iaos_to_atoms_;

    /// Overlap matrix in full basis
    std::shared_ptr<Matrix> S_;
    /// Non-ghosted IAOs in full basis
    std::shared_ptr<Matrix> A_;


    /// Set defaults
    void common_init();

    /// Build the IAOs
    void build_iaos();
    /// Localization task (returns U and L)
    static std::map<std::string, std::shared_ptr<Matrix> > localize_task(
        std::shared_ptr<Matrix> L,                        // Matrix of <i|m> [nocc x nmin]
        const std::vector<std::vector<int> >& minao_inds,   // List of minao indices per active center
        const std::vector<std::pair<int, int> >& rot_inds,  // List of allowed rotations (unique)
        double convergence,                                 // Convergence criterion
        int maxiter,                                        // Maximum number of iterations
        int power                                           // Localization metric power
        );
    /// Energy-ordered local orbital permutation [nmo(local) x nmo(ordered)]
    static std::shared_ptr<Matrix> reorder_orbitals(
        std::shared_ptr<Matrix> F,
        const std::vector<int>& ranges);
    /// Orbital atomic charges (natom x nmo)
    std::shared_ptr<Matrix> orbital_charges(
        std::shared_ptr<Matrix> L
        );


public:

    // => Constructors <= //

    IBOLocalizer2(
        std::shared_ptr<BasisSet> primary,
        std::shared_ptr<BasisSet> minao,
        std::shared_ptr<Matrix> C);

    virtual ~IBOLocalizer2();

    /// Build IBO with defaults from Options object (including MINAO_BASIS)
    static std::shared_ptr<IBOLocalizer2> build(
        std::shared_ptr<BasisSet> primary,
        std::shared_ptr<BasisSet> minao,
        std::shared_ptr<Matrix> C,
        Options& options);

    // => Computers <= //

    /// Print out the localization algorithm and parameters
    virtual void print_header() const;

    /// Localize the orbitals, returns the matrices L [nbf x nmo], U [nmo(dlocal) x nmo(local)], and F [nmo x nmo]
    std::map<std::string, std::shared_ptr<Matrix> > localize(
        std::shared_ptr<Matrix> Cocc,              // Orbitals to localize [nbf x nmo], must live in C above
        std::shared_ptr<Matrix> Focc,              // Fock matrix of orbitals to localize [nmo x nmo]
        const std::vector<int>& ranges = std::vector<int>() // [0, nfocc, nocc] will separately localize core and valence
        );
    /// Print the charges
    void print_charges(double scale = 2.0);

    // => Knobs <= //

    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }
    void set_bench(int bench) { bench_ = bench; }
    void set_convergence(double convergence) { convergence_ = convergence; }
    void set_maxiter(int maxiter) { maxiter_ = maxiter; }
    void set_use_ghosts(bool use_ghosts) { use_ghosts_ = use_ghosts; }
    void set_condition(double condition) { condition_ = condition; }
    void set_power(double power) { power_ = power; }
    void set_use_stars(bool use_stars) { use_stars_ = use_stars; }
    void set_stars_completeness(double stars_completeness) { stars_completeness_ = stars_completeness; }
    void set_stars(const std::vector<int>& stars) { stars_ = stars; }

};

} // Namespace fisapt

} // Namespace psi

#endif
