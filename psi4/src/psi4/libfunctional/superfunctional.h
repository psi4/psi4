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

#ifndef SUPERFUNCTIONAL_H
#define SUPERFUNCTIONAL_H

#include "psi4/libmints/typedefs.h"
#include "psi4/pragma.h"
#include <map>
#include <vector>
#include <cstdlib>
#include <string>
#include <optional>
namespace psi {

class Functional;
class BlockOPoints;

/**
 * SuperFunctional: High-level semilocal DFA object
 *
 * This class directly holds semilocal DFA objects
 * and sketches the definition of the generalized Kohn-Sham functional
 * via alpha and omega parameters for long-range and global hybrid
 * exact exchange mixing, and alpha and omega parameters for long-range
 * and global MP2-like correlation.
 *
 * A DFT functional is defined as:
 *
 * E_XC = E_X + E_C + E(-D)
 * E_X  = (1-\alpha_x) E_X^DFA [\omega_x] + \alpha_x E_X^HF + (1-\alpha_x) E_X^HF,LR
 * E_C  = E_C^DFA [\omega_c] + \alpha_c E_C^MP2
 *
 *
 * E(-D) is an empirical dispersion correction of some form
 *
 **/
class SuperFunctional {
   protected:
    // => Meta-Data <= //
    std::string name_;
    std::string description_;
    std::string citation_;
    std::string xclib_description_;
    bool locked_;

    // => Exchange-side DFA functionals <= //
    std::vector<std::shared_ptr<Functional>> x_functionals_;
    double x_alpha_;
    double x_beta_;
    double x_omega_;

    // => Correlation-side DFA functionals <= //
    std::vector<std::shared_ptr<Functional>> c_functionals_;
    double c_alpha_;
    double c_ss_alpha_;
    double c_os_alpha_;
    double c_omega_;

    // => Asymptotic corrections <= //
    bool needs_grac_;
    std::shared_ptr<Functional> grac_x_functional_;
    std::shared_ptr<Functional> grac_c_functional_;
    double grac_shift_;
    double grac_alpha_;
    double grac_beta_;

    // => VV10 parameters corrections <= //
    bool needs_vv10_;
    double vv10_b_;
    double vv10_c_;
    double vv10_beta_;

    // => Functional values and partials <= //
    bool libxc_xc_func_;
    int max_points_;
    int deriv_;
    // For the functional itself. Each vector elt. is the contribution of a
    // particular point to an integral. Populated by compute_functional. Only fields
    // needed for the particular functional, given LSDA/GGA/Meta status and spin-polarization
    // are populated. Some fields change definition depending on spin polarization.
    // LSDA:
    //   always:
    //      V_RHO_A : ∂/∂ρ_α
    //   if is_unpolarized():
    //      V_RHO_B : ∂/∂ρ_β
    //      V_RHO_A_RHO_A : ∂^2/∂ρ_α^2
    //      V_RHO_A_RHO_B : ∂^2/∂ρ_α∂ρ_β
    //      V_RHO_B_RHO_B : ∂^2/∂ρ_β^2
    //   else:
    //      V_RHO_A_RHO_A : (∂^2/∂ρ_α^2 + ∂^2/∂ρ_α∂ρ_β) / 2
    //
    // GGA:
    //   if is_unpolarized():
    //      V_GAMMA_AA : ∂/∂γ_αα
    //      V_GAMMA_AB : ∂/∂γ_αβ
    //      V_GAMMA_BB : ∂/∂γ_ββ
    //      V_RHO_A_GAMMA_AA : ∂^2/∂ρ_α∂γ_αα
    //      V_RHO_A_GAMMA_AB : ∂^2/∂ρ_α∂γ_αβ
    //      V_RHO_A_GAMMA_BB : ∂^2/∂ρ_α∂γ_ββ
    //      V_RHO_B_GAMMA_AA : ∂^2/∂ρ_β∂γ_αα
    //      V_RHO_B_GAMMA_AB : ∂^2/∂ρ_β∂γ_αβ
    //      V_RHO_B_GAMMA_BB : ∂^2/∂ρ_β∂γ_ββ
    //      V_GAMMA_AA_GAMMA_AA : ∂^2/∂γ_αα^2
    //      V_GAMMA_AA_GAMMA_AB : ∂^2/∂γ_αα∂γ_αβ
    //      V_GAMMA_AA_GAMMA_BB : ∂^2/∂γ_αα∂γ_ββ
    //      V_GAMMA_AB_GAMMA_AB : ∂^2/∂γ_αβ^2
    //      V_GAMMA_AB_GAMMA_BB : ∂^2/∂γ_αβ∂γ_ββ
    //      V_GAMMA_BB_GAMMA_BB : ∂^2/∂γ_ββ^2
    //   else:
    //      V_GAMMA_AA : (∂/∂γ_αα + ∂/∂γ_αβ + ∂/∂γ_ββ) / 2
    //      V_RHO_A_GAMMA_AA : (∂/∂ρ_α (∂/∂γ_αα + ∂/∂γ_αβ + ∂/∂γ_ββ)) / 4
    //      V_GAMMA_AA_GAMMA_AA : (∂^2/∂γ_αα^2 + ∂^2/∂γ_αα∂γ_αβ + ∂^2/∂γ_αα∂γ_ββ + ∂^2/∂γ_αβ^2) / 8
    //        ...note that for closed shells, ∂^2/∂γ_αβ^2 = 2 ∂^2/∂γ_αα∂γ_αβ
    std::map<std::string, SharedVector> values_;
    // For GRAC
    std::map<std::string, SharedVector> ac_values_;
    // For VV10
    std::map<std::string, SharedVector> vv_values_;

    // => Other LibXC settings
    double density_tolerance_;

    // Set up a null Superfunctional
    void common_init();

    // Check if we can edit this Superfunctional
    void can_edit();

   public:
    // => Constructors (Use the factory constructor, or really know what's up) <= //

    SuperFunctional();
    virtual ~SuperFunctional();

    // Build a blank superfunctional
    static std::shared_ptr<SuperFunctional> blank();
    static std::shared_ptr<SuperFunctional> XC_build(std::string name, bool unpolarized, const std::optional<std::map<std::string, double>>& );

    // Clones a *polarized*, complete superfunctional. Used, e.g., in spin-symmetry-
    // breaking eigenvectors of the MO hessian or linear response eigenproblem.
    std::shared_ptr<SuperFunctional> build_polarized();
    // Builds a worker version of the superfunctional
    std::shared_ptr<SuperFunctional> build_worker();

    // Allocate values (MUST be called after adding new functionals to the superfunctional)
    void allocate();

    // => Computers <= //

    // Populates values_
    // If not spin polarized, singlet controls whether singlet or triplet is asked for.
    std::map<std::string, SharedVector>& compute_functional(const std::map<std::string, SharedVector>& vals,
                                                            int npoints = -1, bool singlet = true);
    void test_functional(SharedVector rho_a, SharedVector rho_b, SharedVector gamma_aa, SharedVector gamma_ab,
                         SharedVector gamma_bb, SharedVector tau_a, SharedVector tau_b);

    // Compute the cache data for VV10 dispersion
    std::map<std::string, SharedVector> compute_vv10_cache(const std::map<std::string, SharedVector>& vals,
                                                           std::shared_ptr<BlockOPoints> block, double rho_thresh,
                                                           int npoints = -1, bool internal = false);

    // Copmutes the Cache data for VV10 dispersion
    double compute_vv10_kernel(const std::map<std::string, SharedVector>& vals,
                               const std::vector<std::map<std::string, SharedVector>>& vv10_cache,
                               std::shared_ptr<BlockOPoints> block, int npoints = -1, bool do_grad = false);

    // => Input/Output <= //

    std::map<std::string, SharedVector>& values() { return values_; }
    SharedVector value(const std::string& key) { return values_[key]; }
    SharedVector vv_value(const std::string& key) { return vv_values_[key]; }

    std::vector<std::shared_ptr<Functional>>& x_functionals() { return x_functionals_; }
    std::vector<std::shared_ptr<Functional>>& c_functionals() { return c_functionals_; }
    std::shared_ptr<Functional> grac_x_functional() { return grac_x_functional_; }
    std::shared_ptr<Functional> grac_c_functional() { return grac_c_functional_; }

    std::shared_ptr<Functional> x_functional(const std::string& name);
    std::shared_ptr<Functional> c_functional(const std::string& name);
    void add_x_functional(std::shared_ptr<Functional> fun);
    void add_c_functional(std::shared_ptr<Functional> fun);
    void set_grac_x_functional(std::shared_ptr<Functional> fun) {
        needs_grac_ = true;
        grac_x_functional_ = fun;
    }
    void set_grac_c_functional(std::shared_ptr<Functional> fun) {
        needs_grac_ = true;
        grac_c_functional_ = fun;
    }

    // => Setters <= //

    void set_lock(bool locked) { locked_ = locked; }
    void set_do_vv10(bool do_vv10) { needs_vv10_ = do_vv10; }
    void set_name(const std::string& name) { name_ = name; }
    void set_description(const std::string& description) { description_ = description; }
    void set_citation(const std::string& citation) { citation_ = citation; }
    void set_xclib_description(const std::string& description) { xclib_description_ = description; }

    void set_max_points(int max_points) { max_points_ = max_points; }
    void set_deriv(int deriv) { deriv_ = deriv; }

    void set_x_omega(double omega);
    void set_c_omega(double omega);
    void set_x_alpha(double alpha);
    void set_x_beta(double beta);
    void set_c_alpha(double alpha);
    void set_c_ss_alpha(double alpha);
    void set_c_os_alpha(double alpha);
    void set_vv10_b(double b);
    void set_vv10_c(double c);
    void set_grac_shift(double grac_shift);
    void set_grac_alpha(double grac_alpha);
    void set_grac_beta(double grac_beta);
    void set_density_tolerance(double cut);
    void print_density_threshold(std::string out_fname = "outfile", int print = 1) const;
    void py_print_density_threshold() const { print_density_threshold("outfile", 1); }
    // => Accessors <= //

    std::string name() const { return name_; }
    std::string description() const { return description_; }
    std::string citation() const { return citation_; }
    std::string xclib_description() const { return xclib_description_; }

    int ansatz() const;
    int max_points() const { return max_points_; }
    int deriv() const { return deriv_; }

    double x_omega() const { return x_omega_; }
    double c_omega() const { return c_omega_; }
    double x_alpha() const { return x_alpha_; }
    double x_beta() const { return x_beta_; }
    double c_alpha() const { return c_alpha_; }
    double c_ss_alpha() const { return c_ss_alpha_; }
    double c_os_alpha() const { return c_os_alpha_; }
    double vv10_b() const { return vv10_b_; }
    double vv10_c() const { return vv10_c_; }
    double grac_shift() const { return grac_shift_; }
    double grac_alpha() const { return grac_alpha_; }
    double grac_beta() const { return grac_beta_; }
    double density_tolerance() const { return density_tolerance_; }

    bool needs_xc() const { return ((c_functionals_.size() + x_functionals_.size()) > 0); }
    bool needs_vv10() const { return needs_vv10_; };
    bool needs_grac() const { return needs_grac_; };
    bool PSI_API is_unpolarized() const;
    bool PSI_API is_meta() const;
    bool PSI_API is_gga() const;
    bool is_x_lrc() const { return x_omega_ != 0.0; }
    bool is_c_lrc() const { return c_omega_ != 0.0; }
    bool is_x_hybrid() const { return x_alpha_ != 0.0; }
    bool is_c_hybrid() const { return c_alpha_ != 0.0; }
    bool is_c_scs_hybrid() const { return c_os_alpha_ != 0.0 || c_ss_alpha_ != 0.0; }
    bool is_libxc_func() const { return libxc_xc_func_; }

    // => Utility <= //

    void print(std::string out_fname = "outfile", int print = 1) const;
    void py_print() const { print("outfile", 1); }
    void py_print_detail(int level) const { print("outfile", level); }
};
}  // namespace psi

#endif
