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
#include "psi4/libmints/eigen_interface.h"
#include "psi4/libmints/molecule.h"

#include "psi4/liboptions/liboptions.h"

#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include "psi4/libqt/qt.h"

#include <gauxc/molecular_weights.hpp>
#include <gauxc/molgrid/defaults.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>

namespace psi {

void GauXCBase::initialize() {
    use_gpu_ = options_.get_bool("GAUXC_USE_GPU");
#ifndef GAUXC_HAS_DEVICE
    if (use_gpu_) {
        throw PSIEXCEPTION("GAUXC_USE_GPU was requested, but GauXC was not built with device support.");
    }
#endif
    const auto gauxc_execspace = use_gpu_ ? GauXC::ExecutionSpace::Device : GauXC::ExecutionSpace::Host;
    GauXC::LoadBalancerFactory lb_factory(gauxc_execspace, options_.get_str("GAUXC_LOAD_BALANCER_KERNEL"));
    std::unique_ptr<GauXC::RuntimeEnvironment> rt = nullptr;
 #ifdef GAUXC_HAS_DEVICE 
    if (use_gpu_) {
        rt = std::make_unique<GauXC::DeviceRuntimeEnvironment>( GAUXC_MPI_CODE(MPI_COMM_WORLD,) 0.01*options_.get_int("GAUXC_GPU_MEM"));
    } else {
        rt = std::make_unique<GauXC::RuntimeEnvironment>( GAUXC_MPI_CODE(MPI_COMM_WORLD) );
    }
#else
        rt = std::make_unique<GauXC::RuntimeEnvironment>( GAUXC_MPI_CODE(MPI_COMM_WORLD) );
#endif
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
        init_vec.emplace_back(alpha, ExchCXX::XCKernel(ExchCXX::libxc_name_string(name), this->spin()));
    }
    for (const auto functionalComponent: functional_->c_functionals()) {
        auto xcfunc = dynamic_pointer_cast<LibXCFunctional>(functionalComponent);
        if (xcfunc == nullptr) {
            throw PSIEXCEPTION("GauXC integration requires LibXC functionals.");
        }
        auto name = xcfunc->name();
        auto alpha = xcfunc->alpha();
        init_vec.emplace_back(alpha, ExchCXX::XCKernel(ExchCXX::libxc_name_string(name), this->spin()));
    }
    auto gxc_func = std::make_shared<GauXC::functional_type>(init_vec);
    integrator_ =
          integrator_factory.get_shared_instance(gxc_func, load_balancer);
}

void GauXCBase::print_header() const {
    outfile->Printf("  ===> DFT Potential <===\n\n");
    functional_->print("outfile", options_.get_int("PRINT"));
    outfile->Printf("   => GauXC Quadrature <=\n\n");
    outfile->Printf("    Execution Space        = %14s\n", use_gpu_ ? "Device (GPU)" : "Host (CPU)");
    outfile->Printf("    Radial Scheme          = %14s\n", options_.get_str("GAUXC_RADIAL_SCHEME").c_str());
    outfile->Printf("    Pruning Scheme         = %14s\n", options_.get_str("GAUXC_PRUNING_SCHEME").c_str());
    outfile->Printf("\n");
    outfile->Printf("    Radial Points          = %14d\n", options_.get_int("GAUXC_RADIAL_POINTS"));
    outfile->Printf("    Spherical Points       = %14d\n", options_.get_int("GAUXC_SPHERICAL_POINTS"));
    outfile->Printf("    Grid Batch Size        = %14d\n", options_.get_int("GAUXC_GRID_BATCH_SIZE"));
    outfile->Printf("\n");

}

std::map<std::string, double> GauRV::compute_V(std::vector<SharedMatrix> ret) {
    timer_on("GauRV: Form V");
    Eigen::MatrixXd eigen_d = psi::linalg::eigen_map(*D_AO_[0]);
#if psi4_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD
    if (primary_->has_puream()) {
        auto permuter = linalg::generate_permutation_to_cca(*primary_);
        eigen_d = permuter * eigen_d * permuter.transpose();
    }
#endif
    auto [e_xc, v_xc] = integrator_->eval_exc_vxc(eigen_d);
#if psi4_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD
    if (primary_->has_puream()) {
        auto permuter = linalg::generate_permutation_to_cca(*primary_);
        v_xc = permuter.transpose() * v_xc * permuter;
    }
#endif

    // Set the result
    auto ao_result = psi::linalg::matrix_from_eigen(v_xc);
    if (AO2USO_) {
        (*ret[0]).apply_symmetry(ao_result, *AO2USO_);
    } else {
        (*ret[0]).copy(ao_result);
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
    timer_off("GauRV: Form V");
    return quad_values;
}

std::map<std::string, double> GauUV::compute_V(std::vector<SharedMatrix> ret) {
    timer_on("GauUV: Form V");
    auto Ds = D_AO_[0]->clone();
    Ds->add(D_AO_[1]);
    auto Dz = D_AO_[0]->clone();
    Dz->subtract(D_AO_[1]);
    Eigen::MatrixXd eigen_ds = psi::linalg::eigen_map(*Ds);
    Eigen::MatrixXd eigen_dz = psi::linalg::eigen_map(*Dz);
#if psi4_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD
    if (primary_->has_puream()) {
        auto permuter = linalg::generate_permutation_to_cca(*primary_);
        eigen_ds = permuter * eigen_ds * permuter.transpose();
        eigen_dz = permuter * eigen_dz * permuter.transpose();
    }
#endif
    auto [e_xc, v_s, v_z] = integrator_->eval_exc_vxc(eigen_ds, eigen_dz);
    auto v_a = static_cast<Eigen::MatrixXd>(v_s + v_z);
    auto v_b = static_cast<Eigen::MatrixXd>(v_s - v_z);
#if psi4_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD
    if (primary_->has_puream()) {
        auto permuter = linalg::generate_permutation_to_cca(*primary_);
        v_a = permuter.transpose() * v_a * permuter;
        v_b = permuter.transpose() * v_b * permuter;
    }
#endif

    // Set the result
    auto ao_a = psi::linalg::matrix_from_eigen(v_a);
    auto ao_b = psi::linalg::matrix_from_eigen(v_b);
    if (AO2USO_) {
        ret[0]->apply_symmetry(ao_a, *AO2USO_);
        ret[1]->apply_symmetry(ao_b, *AO2USO_);
    } else {
        ret[0]->copy(ao_a);
        ret[1]->copy(ao_b);
    }

    quad_values_["VV10"] = 0.0;
    quad_values_["FUNCTIONAL"] = e_xc;
    quad_values_["RHO_A"] = 0.0;
    quad_values_["RHO_AX"] = 0.0;
    quad_values_["RHO_AY"] = 0.0;
    quad_values_["RHO_AZ"] = 0.0;
    quad_values_["RHO_B"] = 0.0;
    quad_values_["RHO_BX"] = 0.0;
    quad_values_["RHO_BY"] = 0.0;
    quad_values_["RHO_BZ"] = 0.0;

    timer_off("GauUV: Form V");
    return quad_values_;
}

SharedMatrix GauRV::compute_gradient() {
    Eigen::MatrixXd eigen_d = psi::linalg::eigen_map(*D_AO_[0]);
#if psi4_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD
    if (primary_->has_puream()) {
        auto permuter = linalg::generate_permutation_to_cca(*primary_);
        eigen_d = permuter * eigen_d * permuter.transpose();
    }
#endif
    auto grad_xc = integrator_->eval_exc_grad(eigen_d);

    int natom = primary_->molecule()->natom();
    auto G = std::make_shared<Matrix>("XC Gradient", natom, 3);

    G->copy_from(grad_xc.data());

    quad_values_["FUNCTIONAL"] = 0.0;
    quad_values_["RHO_A"] = 0.0;
    quad_values_["RHO_AX"] = 0.0;
    quad_values_["RHO_AY"] = 0.0;
    quad_values_["RHO_AZ"] = 0.0;
    quad_values_["RHO_B"] = 0.0;
    quad_values_["RHO_BX"] = 0.0;
    quad_values_["RHO_BY"] = 0.0;
    quad_values_["RHO_BZ"] = 0.0;

    return G;
}

SharedMatrix GauUV::compute_gradient() {
    auto Ds = D_AO_[0]->clone();
    Ds->add(D_AO_[1]);
    auto Dz = D_AO_[0]->clone();
    Dz->subtract(D_AO_[1]);
    Eigen::MatrixXd eigen_ds = psi::linalg::eigen_map(*Ds);
    Eigen::MatrixXd eigen_dz = psi::linalg::eigen_map(*Dz);
#if psi4_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD
    if (primary_->has_puream()) {
        auto permuter = linalg::generate_permutation_to_cca(*primary_);
        eigen_ds = permuter * eigen_ds * permuter.transpose();
        eigen_dz = permuter * eigen_dz * permuter.transpose();
    }
#endif
    auto grad_xc = integrator_->eval_exc_grad(eigen_ds, eigen_dz);

    int natom = primary_->molecule()->natom();
    auto G = std::make_shared<Matrix>("XC Gradient", natom, 3);

    G->copy_from(grad_xc.data());

    quad_values_["FUNCTIONAL"] = 0.0;
    quad_values_["RHO_A"] = 0.0;
    quad_values_["RHO_AX"] = 0.0;
    quad_values_["RHO_AY"] = 0.0;
    quad_values_["RHO_AZ"] = 0.0;
    quad_values_["RHO_B"] = 0.0;
    quad_values_["RHO_BX"] = 0.0;
    quad_values_["RHO_BY"] = 0.0;
    quad_values_["RHO_BZ"] = 0.0;

    return G;
}

void GauRV::compute_Vx(const std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret) {
    timer_on("GauRV: Form Vx");

    // Validation
    if (D_AO_.size() != 1) {
        throw PSIEXCEPTION("Vx: RKS should have only one D Matrix");
    }
    if ((Dx.size() != ret.size())) {
        throw PSIEXCEPTION("Vx: RKS input and output sizes should be the same.");
    }
    if (Dx.size() == 0) {
        throw PSIEXCEPTION("Vx: Can't compute with matrix-vector prodcuts with no vectors.");
    }

    Eigen::MatrixXd eigen_d = psi::linalg::eigen_map(*D_AO_[0]);
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> permuter;
#if psi4_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD
    if (primary_->has_puream()) {
        permuter = linalg::generate_permutation_to_cca(*primary_);
        eigen_d = permuter * eigen_d * permuter.transpose();
    }
#endif

    for (size_t i = 0; i < Dx.size(); i++) {
        Eigen::MatrixXd eigen_dx = psi::linalg::eigen_map(*Dx[i]);
#if psi4_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD
        if (primary_->has_puream()) {
            eigen_dx = permuter * eigen_dx * permuter.transpose();
        }
#endif
	GauXC::IntegratorSettingsEXC_GRAD set;
	set.include_weight_derivatives = false;
        auto vx = integrator_->eval_fxc_contraction(eigen_d, eigen_dx, set);
#if psi4_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD
        if (primary_->has_puream()) {
            vx = permuter.transpose() * vx * permuter;
        }
#endif
        ret[i]->copy_from(vx.data());
        ret[i]->scale(0.5);
    }
    timer_off("GauRV: Form Vx");
}

void GauUV::compute_Vx(const std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret) {
    timer_on("GauUV: Form Vx");

    // Validation
    if (D_AO_.size() != 2) {
        throw PSIEXCEPTION("Vx: UKS should have only one D Matrix");
    }
    if ((Dx.size() != ret.size())) {
        throw PSIEXCEPTION("Vx: UKS input and output sizes should be the same.");
    }
    if (Dx.size() == 0) {
        throw PSIEXCEPTION("Vx: Can't compute with matrix-vector prodcuts with no vectors.");
    }
    if (Dx.size() % 2) {
        throw PSIEXCEPTION("Vx: Need an even number of pseudo-densities (alpha, beta pairs).");
    }
    auto num_pairs = Dx.size() / 2;

    auto Ds = D_AO_[0]->clone();
    Ds->add(D_AO_[1]);
    auto Dz = D_AO_[0]->clone();
    Dz->subtract(D_AO_[1]);
    Eigen::MatrixXd eigen_ds = psi::linalg::eigen_map(*Ds);
    Eigen::MatrixXd eigen_dz = psi::linalg::eigen_map(*Dz);
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> permuter;
#if psi4_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD
    if (primary_->has_puream()) {
        permuter = linalg::generate_permutation_to_cca(*primary_);
        eigen_ds = permuter * eigen_ds * permuter.transpose();
        eigen_dz = permuter * eigen_dz * permuter.transpose();
    }
#endif

    for (size_t i = 0; i < num_pairs; i++) {
        auto Ds = Dx[2*i]->clone();
        Ds->add(Dx[2*i+1]);
        auto Dz = Dx[2*i]->clone();
        Dz->subtract(Dx[2*i+1]);
        Eigen::MatrixXd eigen_dxs = psi::linalg::eigen_map(*Ds);
        Eigen::MatrixXd eigen_dxz = psi::linalg::eigen_map(*Dz);
#if psi4_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD
        if (primary_->has_puream()) {
            eigen_dxs = permuter * eigen_dxs * permuter.transpose();
            eigen_dxz = permuter * eigen_dxz * permuter.transpose();
        }
#endif
        auto [vxs, vxz] = integrator_->eval_fxc_contraction(eigen_ds, eigen_dz, eigen_dxs, eigen_dxz);
        auto vx_a = static_cast<Eigen::MatrixXd>(vxs + vxz);
        auto vx_b = static_cast<Eigen::MatrixXd>(vxs - vxz);
#if psi4_SHGSHELL_ORDERING != LIBINT_SHGSHELL_ORDERING_STANDARD
        if (primary_->has_puream()) {
            vx_a = permuter.transpose() * vx_a * permuter;
            vx_b = permuter.transpose() * vx_b * permuter;
        }
#endif

        // Set the result
        auto ao_a = psi::linalg::matrix_from_eigen(vx_a);
        auto ao_b = psi::linalg::matrix_from_eigen(vx_b);
        if (AO2USO_) {
            (*ret[2*i]).apply_symmetry(ao_a, *AO2USO_);
            (*ret[2*i+1]).apply_symmetry(ao_b, *AO2USO_);
        } else {
            (*ret[2*i]).copy(ao_a);
            (*ret[2*i+1]).copy(ao_b);
        }
    }
    timer_off("GauUV: Form Vx");
}

} // namespace psi

#endif
