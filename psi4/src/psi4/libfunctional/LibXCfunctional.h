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

#ifndef LibXC_FUNCTIONAL_H
#define LibXC_FUNCTIONAL_H

#include "psi4/libfunctional/functional.h"
#include "psi4/libmints/typedefs.h"

#include <xc.h>

#include <map>

struct xc_func_type;

/**
 * LibXC functional wrapper
 **/
namespace psi {

class LibXCFunctional : public Functional {
    // Wrapper to the LibXC library

   private:
    std::string xc_func_name_;
    std::unique_ptr<xc_func_type> xc_functional_;
    int func_id_;
    bool user_omega_;
    bool exc_; // Can we compute functional at a point?
    bool vxc_; // Can we compute first derivative at a point?
    bool fxc_; // Can we compute second derivative at a point?

    // **ONLY** Used to pass information up the chain.
    // Is GGA?
    bool gga_;
    // Is Meta?
    bool meta_;
    // Is LRC?
    bool lrc_;
    // Unpolarized? (restricted)
    bool unpolarized_;

    // Exchange
    double global_exch_;
    double lr_exch_;

    // Needs vv10
    bool needs_vv10_;
    double vv10_b_;
    double vv10_c_;

    // User defined tweakers
    // * Libxc needs all set at once as list c. v5.1.0, but store as richer map anyways
    std::map<std::string, double> user_tweakers_;

   public:
    LibXCFunctional(std::string xc_name, bool unpolarized);
    ~LibXCFunctional() override;

    void compute_functional(const std::map<std::string, SharedVector>& in,
                            const std::map<std::string, SharedVector>& out, int npoints, int deriv) override;

    // Clones a *polarized*, complete functional. Used, e.g., in spin-symmetry-
    // breaking eigenvectors of the MO hessian or linear response eigenproblem.
    std::shared_ptr<Functional> build_polarized() override;
    // Clones a *worker* for the functional. This is not a complete functional
    std::shared_ptr<Functional> build_worker() override;

    // Setters and getters
    void set_gga(bool gga) { gga_ = gga; }
    void set_meta(bool meta) { meta_ = meta; }
    void set_density_cutoff(double cut) override;
    void set_omega(double omega);
    void set_tweak(std::map<std::string, double>, bool);
    std::vector<std::tuple<std::string, int, double>> get_mix_data();

    bool is_meta() const override { return meta_; }
    bool is_gga() const override { return gga_; }
    bool is_lrc() const override { return lrc_; }
    bool is_unpolarized() const override { return unpolarized_; }

    // Make queries to libxc
    std::map<std::string, double> query_libxc(const std::string& functional);

    // Only used to pass information up the chain
    double global_exchange() { return global_exch_; }
    double lr_exchange() { return lr_exch_; }
    double needs_vv10() { return needs_vv10_; }
    double vv10_b() { return vv10_b_; }
    double vv10_c() { return vv10_c_; }
    double density_cutoff() const override { return xc_functional_->dens_threshold; }
    
    // Information about what type of LibXC functional this is
    int kind() override { return xc_functional_->info->kind; }

    // Get libxc provenance stamp
    static std::string xclib_description();
};
}  // namespace psi

#endif
