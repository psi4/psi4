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

#ifdef USING_gauxc

#include "gauxc_int.h"

#include "psi4/libfunctional/LibXCfunctional.h"
#include "psi4/libfunctional/superfunctional.h"

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"

#include "psi4/liboptions/liboptions.h"

#include <gauxc/molecular_weights.hpp>
#include <gauxc/molgrid/defaults.hpp>

namespace psi {

void GauRV::initialize() {
    // TODO: Allow for Device execspace, depending on flags. This will add GPU support.
    const auto gauxc_execspace = GauXC::ExecutionSpace::Host;
    GauXC::LoadBalancerFactory lb_factory(gauxc_execspace, options_.get_str("GAUXC_LOAD_BALANCER_KERNEL"));
    auto rt = std::make_unique<GauXC::RuntimeEnvironment>( GAUXC_MPI_CODE(MPI_COMM_WORLD) );
    auto gauxc_mol = primary_->molecule()->to_gauxc_molecule();

    std::unordered_map<std::string, GauXC::PruningScheme> pruning_scheme_map = {
      {"ROBUST", GauXC::PruningScheme::Robust},
      {"TREUTLER", GauXC::PruningScheme::Treutler},
      {"NONE", GauXC::PruningScheme::Unpruned}
    };
    auto pruning_scheme = options_.get_str("GAUXC_PRUNING_SCHEME");

    // generate map for radial quadrature schemes 
    std::unordered_map<std::string, GauXC::RadialQuad> radial_scheme_map = { 
      {"TREUTLER", GauXC::RadialQuad::TreutlerAhlrichs},
      {"MURA", GauXC::RadialQuad::MuraKnowles},
      {"EM", GauXC::RadialQuad::MurrayHandyLaming}
    };
    auto radial_scheme = options_.get_str("GAUXC_RADIAL_SCHEME");
    auto grid_batch_size = options_.get_int("GAUXC_GRID_BATCH_SIZE");
    auto radial_points = options_.get_int("GAUXC_RADIAL_POINTS");
    auto spherical_points = options_.get_int("GAUXC_SPHERICAL_POINTS");

    auto gauxc_grid = GauXC::MolGridFactory::create_default_molgrid(
        gauxc_mol, 
        pruning_scheme_map[pruning_scheme],
        GauXC::BatchSize(grid_batch_size), // TODO: Value taken from MPQC. Experimental.
        radial_scheme_map[radial_scheme], 
        GauXC::RadialSize(radial_points),
        GauXC::AngularSize(spherical_points)
    );
    auto gauxc_primary = primary_->to_gauxc_basisset<double>(1.0e-10, false); // TODO: Allow customization
    auto load_balancer = lb_factory.get_shared_instance(*rt, gauxc_mol, gauxc_grid, gauxc_primary);
    GauXC::MolecularWeightsFactory mw_factory(gauxc_execspace, "Default",
                                          GauXC::MolecularWeightsSettings{});
    auto mw = mw_factory.get_instance();
    mw.modify_weights(*load_balancer);

    // TODO: Allow for more options here. This is quick-and-dirty.
    GauXC::XCIntegratorFactory<Eigen::MatrixXd> integrator_factory(
          gauxc_execspace, "Replicated", "Default", "Default", "Default");
    std::vector<std::pair<double, ExchCXX::XCKernel>> init_vec;

    for (const auto functionalComponent: functional_->x_functionals()) {
        auto xcfunc = dynamic_pointer_cast<LibXCFunctional>(functionalComponent);
        if (xcfunc == nullptr) {
            throw PSIEXCEPTION("GauXC integration requires LibXC functionals.");
        }
        auto name = xcfunc->name();
        auto alpha = xcfunc->alpha();
        // TODO: Update Spin.
        init_vec.emplace_back(alpha, ExchCXX::XCKernel(ExchCXX::libxc_name_string(name), ExchCXX::Spin::Unpolarized));
    }
    for (const auto functionalComponent: functional_->c_functionals()) {
        auto xcfunc = dynamic_pointer_cast<LibXCFunctional>(functionalComponent);
        if (xcfunc == nullptr) {
            throw PSIEXCEPTION("GauXC integration requires LibXC functionals.");
        }
        auto name = xcfunc->name();
        auto alpha = xcfunc->alpha();
        init_vec.emplace_back(alpha, ExchCXX::XCKernel(ExchCXX::libxc_name_string(name), ExchCXX::Spin::Unpolarized));
    }
    auto gxc_func = std::make_shared<GauXC::functional_type>(init_vec);
    integrator_ =
          integrator_factory.get_shared_instance(gxc_func, load_balancer);
}

std::map<std::string, double> GauRV::compute_V(std::vector<SharedMatrix> ret) {
    Eigen::MatrixXd eigen_d = D_AO_[0]->eigen_map();
#if psi4_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD
    if (primary_->has_puream()) {
        auto permuter = primary_->generate_permutation_to_cca();
        eigen_d = permuter * eigen_d * permuter.transpose();
    }
#endif
    auto [e_xc, v_xc] = integrator_->eval_exc_vxc(eigen_d);
#if psi4_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD
    if (primary_->has_puream()) {
        auto permuter = primary_->generate_permutation_to_cca();
        v_xc = permuter.transpose() * v_xc * permuter;
    }
#endif

    // Set the result
    auto ao_result = std::make_shared<Matrix>(v_xc);
    if (AO2USO_) {
        ret[0]->apply_symmetry(ao_result, AO2USO_);
    } else {
        ret[0]->copy(ao_result);
    }

    std::map<std::string, double> quad_values;
    quad_values["VV10"] = 0.0;
    quad_values["FUNCTIONAL"] = e_xc;
    quad_values["RHO_A"] = 0.0;
    quad_values["RHO_AX"] = 0.0;
    quad_values["RHO_AY"] = 0.0;
    quad_values["RHO_AZ"] = 0.0;
    quad_values["RHO_B"] = 0.0;
    quad_values["RHO_BX"] = 0.0;
    quad_values["RHO_BY"] = 0.0;
    quad_values["RHO_BZ"] = 0.0;
    return quad_values;
}


} // namespace psi

#endif
