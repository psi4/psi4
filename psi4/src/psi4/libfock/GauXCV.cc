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

#include "GauXCV.h"

#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/psi4-dec.h"

#ifdef USING_gauxc
#include <cmath>

#include "gauxc_interface.h"

#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libqt/qt.h"

#include "psi4/libfunctional/functional.h"

#include <exchcxx/xc_functional.hpp>

#include <gauxc/load_balancer.hpp>
#include <gauxc/molecular_weights.hpp>
#include <gauxc/molgrid/defaults.hpp>
#include <gauxc/util/environment.hpp>
#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/impl.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>
#endif

namespace psi {

#ifdef USING_gauxc
GauXCEngine::GauXCEngine(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary,
                         Options& options, bool polarized)
    : primary_(primary), options_(options) {
    use_gpu_ = options_.get_bool("DFT_V_USE_GPU");
    const auto ex = use_gpu_ ? GauXC::ExecutionSpace::Device : GauXC::ExecutionSpace::Host;
    runtime_ = gauxc_interface::make_runtime_environment(use_gpu_, 0.01 * options_.get_int("DFT_V_GPU_MEM"));

    const int gauxc_max_am = GauXC::gauxc_max_am(ex, GauXC::SupportedAlg::XC);
    if (primary_->max_am() > gauxc_max_am) {
        throw PSIEXCEPTION("Basis set has higher-AM shells (Max AM = " + std::to_string(primary_->max_am()) +
                           ") than this GauXC instance supports (Max AM = " + std::to_string(gauxc_max_am) +
                           "). Use a smaller basis or a GauXC install with higher AM support.");
    }

    force_cartesian_ = use_gpu_ && primary_->has_puream();
    if (force_cartesian_) {
        MintsHelper helper(primary_, options_, 0);
        cartao_to_ao_matrix_ = helper.cartao_to_ao_transform();
        if (cartao_to_ao_matrix_->nirrep() != 1) {
            throw PSIEXCEPTION(
                "GPU execution of the GauXC XC quadrature with a spherical basis currently requires C1 symmetry.");
        }
    }

#if psi4_SHGSHELL_ORDERING == LIBINT_SHGSHELL_ORDERING_GAUSSIAN
    if (primary_->has_puream()) {
        permutation_matrix_ = gauxc_interface::generate_permutation_matrix(primary_, gauxc_max_am);
    }
#endif

    auto gauxc_mol = gauxc_interface::psi4_to_gauxc_molecule(primary_->molecule());
    auto gauxc_basis = gauxc_interface::psi4_to_gauxc_basisset<double>(primary_, 1.0e-10, force_cartesian_);

    auto gauxc_grid = GauXC::MolGridFactory::create_default_molgrid(
        gauxc_mol, gauxc_interface::to_gauxc_pruning_scheme(options_.get_str("DFT_PRUNING_SCHEME")),
        GauXC::BatchSize(512), gauxc_interface::to_gauxc_radial_scheme(options_.get_str("DFT_RADIAL_SCHEME")),
        GauXC::RadialSize(options_.get_int("DFT_RADIAL_POINTS")),
        GauXC::AngularSize(options_.get_int("DFT_SPHERICAL_POINTS")));

    GauXC::LoadBalancerFactory lb_factory(ex, "Default");
    auto load_balancer = lb_factory.get_instance(*runtime_, gauxc_mol, gauxc_grid, gauxc_basis);

    GauXC::MolecularWeightsSettings mw_settings;
    mw_settings.weight_alg = gauxc_interface::to_gauxc_weight_scheme(options_.get_str("DFT_NUCLEAR_SCHEME"));
    GauXC::MolecularWeightsFactory mw_factory(ex, "Default", mw_settings);
    mw_factory.get_instance().modify_weights(load_balancer);

    // Mirror Psi4's own SuperFunctional decomposition instead of matching a
    // composite name: every ingredient Psi4 uses carries its libxc kernel name
    // ("XC_LDA_X", ...), and ExchCXX can construct kernels from exactly those.
    const auto spin = polarized ? ExchCXX::Spin::Polarized : ExchCXX::Spin::Unpolarized;
    std::vector<std::pair<double, ExchCXX::XCKernel> > kernels;
    auto add_kernels = [&](const std::vector<std::shared_ptr<Functional> >& fs) {
        for (const auto& f : fs) {
            kernels.emplace_back(f->alpha(), ExchCXX::XCKernel(ExchCXX::libxc_name_string(f->name()), spin));
        }
    };
    add_kernels(functional->x_functionals());
    add_kernels(functional->c_functionals());
    if (kernels.empty()) {
        throw PSIEXCEPTION("GauXC XC quadrature: SuperFunctional carries no exchange-correlation ingredients.");
    }
    if (functional->needs_vv10()) {
        throw PSIEXCEPTION(
            "DFT_V_ALGORITHM=GAUXC does not support the VV10 non-local correlation kernel. "
            "Use DFT_V_ALGORITHM=INTERNAL for this functional.");
    }
    if (functional->needs_grac()) {
        throw PSIEXCEPTION(
            "DFT_V_ALGORITHM=GAUXC does not support GRAC asymptotic corrections. "
            "Use DFT_V_ALGORITHM=INTERNAL for this functional.");
    }

    ExchCXX::HybCoeffs hyb;
    hyb.alpha = functional->x_alpha();
    if (functional->is_x_lrc()) {
        hyb.beta = functional->x_beta();
        hyb.omega = functional->x_omega();
    }
    ExchCXX::XCFunctional func(kernels, hyb);

    integrator_factory_ =
        std::make_unique<GauXC::XCIntegratorFactory<matrix_type> >(ex, "Replicated", "Default", "Default", "Default");
    integrator_ = integrator_factory_->get_shared_instance(func, load_balancer);
}

GauXCEngine::matrix_type GauXCEngine::to_gauxc_density(SharedMatrix D) const {
    SharedMatrix work = D;
    if (force_cartesian_) {
        work = std::make_shared<Matrix>();
        work->transform(D, cartao_to_ao_matrix_);
    }
    matrix_type P = work->eigen_map();
    if (permutation_matrix_.has_value()) {
        const auto& perm = permutation_matrix_.value();
        P = perm * P * perm.transpose();
    }
    return P;
}

void GauXCEngine::from_gauxc_potential(const matrix_type& V_gauxc, SharedMatrix V_out) const {
    matrix_type V = V_gauxc;
    if (permutation_matrix_.has_value()) {
        const auto& perm = permutation_matrix_.value();
        V = perm.transpose() * V * perm;
    }
    if (force_cartesian_) {
        auto buffer = std::make_shared<Matrix>(V.rows(), V.cols());
        buffer->eigen_map() = V;
        V_out->back_transform(buffer, cartao_to_ao_matrix_);
    } else {
        V_out->eigen_map() = V;
    }
}

GauXCEngine& GauXCRV::engine() const {
    if (!engine_) {
        engine_ = std::make_unique<GauXCEngine>(functional_, primary_, options_, /*polarized=*/false);
    }
    return *engine_;
}

GauXCEngine& GauXCUV::engine() const {
    if (!engine_) {
        engine_ = std::make_unique<GauXCEngine>(functional_, primary_, options_, /*polarized=*/true);
    }
    return *engine_;
}
#endif

GauXCRV::GauXCRV(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options)
    : RV(functional, primary, options) {}

GauXCRV::~GauXCRV() {}

void GauXCRV::compute_V(std::vector<SharedMatrix> ret) {
#ifdef USING_gauxc
    timer_on("GauXCRV: Form V");
    if ((D_AO_.size() != 1) || (ret.size() != 1)) {
        throw PSIEXCEPTION("GauXCRV: RKS should have only one D/V Matrix");
    }

    // Psi4 hands over the alpha density (see v.h) and GauXC's RKS entry point
    // expects exactly that -- it forms the total density internally. Verified
    // empirically: passing 2*D scales E_xc by 2^(4/3) for LDA exchange.
    GauXCEngine::matrix_type P = engine().to_gauxc_density(D_AO_[0]);
    auto [exc, vxc] = engine().integrator().eval_exc_vxc(P);

    auto V_AO = std::make_shared<Matrix>("V AO Temp", nbf_, nbf_);
    engine().from_gauxc_potential(vxc, V_AO);

    if (AO2USO_) {
        ret[0]->apply_symmetry(V_AO, AO2USO_);
    } else {
        ret[0]->copy(V_AO);
    }

    quad_values_["FUNCTIONAL"] = exc;
    quad_values_["VV10"] = 0.0;
    quad_values_["RHO_A"] = 0.0;
    quad_values_["RHO_AX"] = 0.0;
    quad_values_["RHO_AY"] = 0.0;
    quad_values_["RHO_AZ"] = 0.0;
    quad_values_["RHO_B"] = 0.0;
    quad_values_["RHO_BX"] = 0.0;
    quad_values_["RHO_BY"] = 0.0;
    quad_values_["RHO_BZ"] = 0.0;

    if (std::isnan(quad_values_["FUNCTIONAL"])) {
        throw PSIEXCEPTION("GauXCRV: Integrated DFT functional to get NaN.");
    }
    timer_off("GauXCRV: Form V");
#else
    throw PSIEXCEPTION("DFT_V_ALGORITHM=GAUXC requires Psi4 built with ENABLE_gauxc=ON.");
#endif
}

void GauXCRV::print_header() const { RV::print_header(); }

GauXCUV::GauXCUV(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options)
    : UV(functional, primary, options) {}

GauXCUV::~GauXCUV() {}

void GauXCUV::compute_V(std::vector<SharedMatrix> ret) {
#ifdef USING_gauxc
    timer_on("GauXCUV: Form V");
    if ((D_AO_.size() != 2) || (ret.size() != 2)) {
        throw PSIEXCEPTION("GauXCUV: UKS should have two D/V Matrices");
    }
    if (functional_->needs_grac()) {
        throw PSIEXCEPTION("GauXCUV: UKS cannot compute GRAC corrections.");
    }

    // GauXC's UKS entry point takes the scalar and the spin-magnetization
    // density, Ps = Da + Db and Pz = Da - Db. See GauXC's reference host work
    // driver, which forms rho_+/- as 0.5*(rho_s +/- rho_z).
    GauXCEngine::matrix_type Da = engine().to_gauxc_density(D_AO_[0]);
    GauXCEngine::matrix_type Db = engine().to_gauxc_density(D_AO_[1]);
    GauXCEngine::matrix_type Ps = Da + Db;
    GauXCEngine::matrix_type Pz = Da - Db;

    auto [exc, vxc_s, vxc_z] = engine().integrator().eval_exc_vxc(Ps, Pz);

    // Chain rule inverse of the density combination above.
    GauXCEngine::matrix_type Va = vxc_s + vxc_z;
    GauXCEngine::matrix_type Vb = vxc_s - vxc_z;

    auto Va_AO = std::make_shared<Matrix>("Va AO Temp", nbf_, nbf_);
    auto Vb_AO = std::make_shared<Matrix>("Vb AO Temp", nbf_, nbf_);
    engine().from_gauxc_potential(Va, Va_AO);
    engine().from_gauxc_potential(Vb, Vb_AO);

    if (AO2USO_) {
        ret[0]->apply_symmetry(Va_AO, AO2USO_);
        ret[1]->apply_symmetry(Vb_AO, AO2USO_);
    } else {
        ret[0]->copy(Va_AO);
        ret[1]->copy(Vb_AO);
    }

    quad_values_["FUNCTIONAL"] = exc;
    quad_values_["VV10"] = 0.0;
    quad_values_["RHO_A"] = 0.0;
    quad_values_["RHO_AX"] = 0.0;
    quad_values_["RHO_AY"] = 0.0;
    quad_values_["RHO_AZ"] = 0.0;
    quad_values_["RHO_B"] = 0.0;
    quad_values_["RHO_BX"] = 0.0;
    quad_values_["RHO_BY"] = 0.0;
    quad_values_["RHO_BZ"] = 0.0;

    if (std::isnan(quad_values_["FUNCTIONAL"])) {
        throw PSIEXCEPTION("GauXCUV: Integrated DFT functional to get NaN.");
    }
    timer_off("GauXCUV: Form V");
#else
    throw PSIEXCEPTION("DFT_V_ALGORITHM=GAUXC requires Psi4 built with ENABLE_gauxc=ON.");
#endif
}

void GauXCUV::print_header() const { UV::print_header(); }

}  // namespace psi
