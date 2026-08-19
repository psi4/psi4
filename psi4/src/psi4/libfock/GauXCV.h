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

#ifndef GAUXCV_H
#define GAUXCV_H

#include "v.h"

#ifdef USING_gauxc
#include <memory>
#include <optional>

#include <Eigen/Core>
#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>

#include "gauxc_interface.h"
#endif

namespace psi {

#ifdef USING_gauxc
/// Holds the entire GauXC state and all conversions between the Psi4 and
/// GauXC conventions. Shared by GauXCRV and GauXCUV.
class GauXCEngine {
   public:
    using matrix_type = Eigen::MatrixXd;

    /// Builds grid, load balancer and integrator. Throws on unsupported cases.
    GauXCEngine(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options,
                bool polarized);

    /// Psi4 density -> GauXC convention (permutation, sph->cart). No spin scaling.
    matrix_type to_gauxc_density(SharedMatrix D) const;
    /// GauXC potential -> Psi4 convention, written into V_out.
    void from_gauxc_potential(const matrix_type& V_gauxc, SharedMatrix V_out) const;

    GauXC::XCIntegrator<matrix_type>& integrator() { return *integrator_; }
    bool use_gpu() const { return use_gpu_; }
    bool force_cartesian() const { return force_cartesian_; }

   private:
    std::shared_ptr<BasisSet> primary_;
    Options& options_;

    bool use_gpu_ = false;
    bool force_cartesian_ = false;
    SharedMatrix cartao_to_ao_matrix_ = nullptr;
    std::optional<Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> > permutation_matrix_;
    std::shared_ptr<GauXC::XCIntegrator<matrix_type> > integrator_;
    std::unique_ptr<GauXC::XCIntegratorFactory<matrix_type> > integrator_factory_;
    std::unique_ptr<GauXC::RuntimeEnvironment> runtime_;
};
#endif

/// RKS XC quadrature via the GauXC library. Inherits all derivative paths from RV.
class GauXCRV : public RV {
   protected:
#ifdef USING_gauxc
    /// `mutable` because the const-qualified print_header also needs the engine.
    mutable std::unique_ptr<GauXCEngine> engine_;
    /// Builds the engine on first use, then caches it.
    GauXCEngine& engine() const;
#endif

   public:
    GauXCRV(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options);
    ~GauXCRV() override;

    void compute_V(std::vector<SharedMatrix> ret) override;
    void print_header() const override;
};

/// UKS XC quadrature via the GauXC library. Inherits all derivative paths from UV.
class GauXCUV : public UV {
   protected:
#ifdef USING_gauxc
    mutable std::unique_ptr<GauXCEngine> engine_;
    GauXCEngine& engine() const;
#endif

   public:
    GauXCUV(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options);
    ~GauXCUV() override;

    void compute_V(std::vector<SharedMatrix> ret) override;
    void print_header() const override;
};

}  // namespace psi
#endif
