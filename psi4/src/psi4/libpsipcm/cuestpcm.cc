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

 #ifdef USING_cuEST

 #include <utility> 
 #include <map>
 #include <string>
 #include "psi4/libmints/basisset.h"
 #include "psi4/libmints/matrix.h"
 #include "psi4/libmints/molecule.h"
 #include "psi4/libmints/mintshelper.h"
 #include "psi4/liboptions/liboptions.h"
 #include "psi4/libpsipcm/cuestpcm.h"
 #include "psi4/physconst.h"
 #include "cuest.h"
 #include "psi4/libfock/cuESTCommon.h"
 extern cuestHandle_t cuest_handle;
 
namespace psi {
namespace detail {
    /// Optimized zeta values for given Lebedev grid size
    std::map<int, double> pcm_optimized_zetas_ = {
        {14, 4.865},
        {26, 4.855},
        {50, 4.893},
        {110, 4.901},
        {194, 4.903},
        {302, 4.905},
        {434, 4.906},
        {590, 4.905},
        {770, 4.899},
        {974, 4.907},
        {1202, 4.907}      
    };
    /// Bondi radii for atoms, in Angstrom
    std::map<std::string, double> pcm_bondi_radii_ = {
        {"H", 1.10},
        {"HE", 1.40},
        {"LI", 1.81},
        {"BE", 1.53},
        {"B", 1.92},
        {"C", 1.70},
        {"N", 1.55},
        {"O", 1.52},
        {"F", 1.47},
        {"NE", 1.54},
        {"NA", 2.27},
        {"MG", 1.73},
        {"AL", 1.84},
        {"SI", 2.10},
        {"P", 1.80},
        {"S", 1.80},
        {"CL", 1.75},
        {"AR", 1.88},
        {"K", 2.75},
        {"CA", 2.31},
        {"SC", 2.15},
        {"TI", 2.11},
        {"V", 2.07},
        {"CR", 2.06},
        {"MN", 2.05},
        {"FE", 2.04},
        {"CO", 2.00},
        {"NI", 1.97},
        {"CU", 1.96},
        {"ZN", 2.01},
        {"GA", 1.87},
        {"GE", 2.11},
        {"AS", 1.85},
        {"SE", 1.90},
        {"BR", 1.83},
        {"KR", 2.02},
        {"RB", 3.03},
        {"SR", 2.49},
        {"Y", 2.32},
        {"ZR", 2.23},
        {"NB", 2.18},
        {"MO", 2.17},
        {"TC", 2.16},
        {"RU", 2.13},
        {"RH", 2.10},
        {"PD", 2.10},
        {"AG", 2.11},
        {"CD", 2.18},
        {"IN", 1.93},
        {"SN", 2.17},
        {"SB", 2.06},
        {"TE", 2.06},
        {"I", 1.98},
        {"XE", 2.16},
        {"CS", 3.43},
        {"BA", 2.68},
        {"LA", 2.43},
        {"CE", 2.42},
        {"PR", 2.40},
        {"ND", 2.39},
        {"TB", 2.33},
        {"DY", 2.31},
        {"HO", 2.30},
        {"ER", 2.29},
        {"TM", 2.27},
        {"YB", 2.26},
        {"LU", 2.24},
        {"TL", 1.96},
        {"PB", 2.02},
        {"BI", 2.07},
        {"PO", 1.97},
        {"AT", 2.02},
        {"RN", 2.20},
        {"FR", 3.48},
        {"RA", 2.83},
        {"AC", 2.47},
        {"TH", 2.45},
        {"PA", 2.43},
        {"U", 2.41},
        {"NP", 2.39},
        {"PU", 2.43},
        {"AM", 2.44},
        {"CM", 2.45},
        {"BK", 2.44},
        {"CF", 2.45},
        {"ES", 2.45},
        {"FM", 2.45},
        {"MD", 2.46},
        {"NO", 2.46},
        {"LR", 2.46},    
    };
} // namespace detail

cuestPCM::cuestPCM(const Options &options, const std::shared_ptr<MintsHelper> &mintshelper) :
    mintshelper_(mintshelper) {
    // Create the integral plan
    const std::shared_ptr<Molecule> molecule = mintshelper->basisset()->molecule();
    size_t natom = molecule->natom();
    size_t nbf = mintshelper->basisset()->nbf();

    cuestPCMIntPlanParameters_t pcm_int_params;
    CHECK_CUEST(cuestParametersCreate(CUEST_PCMINTPLAN_PARAMETERS, reinterpret_cast<void**>(&pcm_int_params)));
    int heavy_atom_spherical_points = options.get_int("CUEST_PCM_HEAVY_ATOM_SPHERICAL_POINTS");
    double heavy_atom_zeta = 0.0;
    try {
        heavy_atom_zeta = detail::pcm_optimized_zetas_[heavy_atom_spherical_points];
    } catch (...) {
        throw PSIEXCEPTION("Invalid PCM_HEAVY_ATOM_SPHERICAL_POINTS; must be a valid Lebedev grid size.");
    }

    int hydrogen_atom_spherical_points = options.get_int("CUEST_PCM_HYDROGEN_ATOM_SPHERICAL_POINTS");
    double hydrogen_atom_zeta = 0.0;
    try {
        hydrogen_atom_zeta = detail::pcm_optimized_zetas_[hydrogen_atom_spherical_points];
    } catch (...) {
        throw PSIEXCEPTION("Invalid PCM_HYDROGEN_ATOM_SPHERICAL_POINTS; must be a valid Lebedev grid size.");
    }

    std::vector<double> radii(natom);
    std::vector<double> zetas(natom);
    std::vector<double> charges(natom);
    std::vector<uint64_t> grid_sizes(natom);
    for (size_t i = 0; i < natom; i++) {
        std::string symbol = molecule->symbol(i);
        std::transform(symbol.begin(), symbol.end(), symbol.begin(),
            [](unsigned char c){ return std::toupper(c); });
        if (detail::pcm_bondi_radii_.find(symbol) != detail::pcm_bondi_radii_.end()) {
            double bondi_radius = options.get_double("CUEST_PCM_BONDI_RADII_SCALE") * detail::pcm_bondi_radii_[symbol];
            radii[i] = bondi_radius / pc_bohr2angstroms;
        } else {
            throw PSIEXCEPTION("Invalid atom symbol " + symbol + "; must be a valid atom symbol.");
        }
        zetas[i] = symbol == "H" ? hydrogen_atom_zeta : heavy_atom_zeta;
        charges[i] = molecule->Z(i);
        grid_sizes[i] = symbol == "H" ? hydrogen_atom_spherical_points : heavy_atom_spherical_points;
    }
    double pcm_epsilon = options.get_double("CUEST_PCM_DIELECTRIC");

    cuestPCMIntPlanParameters_t pcm_integral_plan_parameters;
    CHECK_CUEST(cuestParametersCreate(CUEST_PCMINTPLAN_PARAMETERS, reinterpret_cast<void**>(&pcm_integral_plan_parameters)));
    double pcm_x = options.get_double("CUEST_PCM_X_PREFACTOR");
    CHECK_CUEST(cuestParametersConfigure(CUEST_PCMINTPLAN_PARAMETERS, pcm_integral_plan_parameters, CUEST_PCMINTPLAN_PARAMETERS_X_PREFACTOR, &pcm_x, sizeof(pcm_x)));
    double pcm_cutoff = options.get_double("CUEST_PCM_CUTOFF");
    CHECK_CUEST(cuestParametersConfigure(CUEST_PCMINTPLAN_PARAMETERS, pcm_integral_plan_parameters, CUEST_PCMINTPLAN_PARAMETERS_CUTOFF, &pcm_cutoff, sizeof(pcm_cutoff)));
    pcm_convergence_ = options.get_double("CUEST_PCM_CONVERGENCE");
    cuestWorkspaceDescriptor_t temporary_workspace_descriptor{};
    cuestWorkspaceDescriptor_t persistent_workspace_descriptor{};
    CHECK_CUEST(cuestPCMIntPlanCreateWorkspaceQuery(
        cuest_handle,
        mintshelper_->cuest_oeint_plan(),
        pcm_integral_plan_parameters,
        &persistent_workspace_descriptor,
        &temporary_workspace_descriptor,
        grid_sizes.data(),
        pcm_epsilon,
        zetas.data(),
        radii.data(),
        charges.data(),
        &pcm_integral_plan_));
    pcm_integral_ws_ptr_ = cuest_common::allocateWorkspace(&persistent_workspace_descriptor);
    cuestWorkspace_t* temporary_workspace = cuest_common::allocateWorkspace(&temporary_workspace_descriptor);
    CHECK_CUEST(cuestPCMIntPlanCreate(
        cuest_handle,
        mintshelper_->cuest_oeint_plan(),
        pcm_integral_plan_parameters,
        pcm_integral_ws_ptr_,
        temporary_workspace,
        grid_sizes.data(),
        pcm_epsilon,
        zetas.data(),
        radii.data(),
        charges.data(),
        &pcm_integral_plan_));
    cuest_common::freeWorkspace(temporary_workspace);
    CHECK_CUEST(cuestParametersDestroy(CUEST_PCMINTPLAN_PARAMETERS, pcm_integral_plan_parameters));

    uint64_t npoints = 0;
    CHECK_CUEST(cuestQuery(
        cuest_handle,
        CUEST_PCMINTPLAN,
        pcm_integral_plan_,
        CUEST_PCMINTPLAN_NUM_POINT,
        &npoints,
        sizeof(npoints)
        ));

    cudaError_t err = cudaMalloc((void**)&d_q1_matrix_, npoints * sizeof(double));
    if (err != cudaSuccess) {
        throw PSIEXCEPTION("cudaMalloc failed in cuestPCM constructor");
    }
    cudaMemset(d_q1_matrix_, 0, npoints * sizeof(double));
    err = cudaMalloc((void**)&d_q2_matrix_, npoints * sizeof(double));
    if (err != cudaSuccess) {
        throw PSIEXCEPTION("cudaMalloc failed in cuestPCM constructor");
    }
    cudaMemset(d_q2_matrix_, 0, npoints * sizeof(double));
    err = cudaMalloc((void**)&d_V_matrix_, nbf * nbf * sizeof(double));
    if (err != cudaSuccess) {
        throw PSIEXCEPTION("cudaMalloc failed in cuestPCM constructor");
    }
    err = cudaMalloc((void**)&d_D_matrix_, nbf * nbf * sizeof(double));
    if (err != cudaSuccess) {
        throw PSIEXCEPTION("cudaMalloc failed in cuestPCM constructor");
    }
}

SharedMatrix cuestPCM::compute_PCM_gradient(const SharedMatrix &D) {
    int natom = mintshelper_->basisset()->molecule()->natom();
    SharedMatrix grad = std::make_shared<Matrix>("PCM Gradient", natom, 3);
    double *d_V_gradient_matrix = nullptr;
    cudaError_t err = cudaMalloc((void**)&d_V_gradient_matrix, natom * 3 * sizeof(double));
    if (err != cudaSuccess) {
        throw PSIEXCEPTION("cudaMalloc failed in cuestPCM::compute_PCM_gradient");
    }

    cuestPCMResults_t pcm_results;
    CHECK_CUEST(cuestResultsCreate(CUEST_PCM_RESULTS, &pcm_results));

    // Copy the density matrix to the device
    int nbf = mintshelper_->basisset()->nbf();
    err = cudaMemcpy(d_D_matrix_, D->pointer()[0], nbf * nbf * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        throw PSIEXCEPTION("cudaMemcpy failed in cuestPCM::compute_PCM_terms");
    }
    cuestPCMDerivativeComputeParameters_t pcm_derivative_compute_parameters;
    CHECK_CUEST(cuestParametersCreate(CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS, &pcm_derivative_compute_parameters));
    CHECK_CUEST(cuestParametersConfigure(CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS, pcm_derivative_compute_parameters, CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS_CONVERGENCE_THRESHOLD, &pcm_convergence_, sizeof(pcm_convergence_)));
    cuestWorkspaceDescriptor_t temporary_workspace_descriptor{};

    CHECK_CUEST(cuestPCMDerivativeComputeWorkspaceQuery(
        cuest_handle,
        pcm_integral_plan_,
        pcm_derivative_compute_parameters,
        &temporary_workspace_descriptor,
        d_D_matrix_,
        d_q1_matrix_,
        d_q2_matrix_,
        pcm_results,
        d_V_gradient_matrix
    ));
    cuestWorkspace_t* temporary_workspace = cuest_common::allocateWorkspace(&temporary_workspace_descriptor);
    CHECK_CUEST(cuestPCMDerivativeCompute(
        cuest_handle,
        pcm_integral_plan_,
        pcm_derivative_compute_parameters,
        temporary_workspace,
        d_D_matrix_,
        d_q1_matrix_,
        d_q2_matrix_,
        pcm_results,
        d_V_gradient_matrix
    ));
    cuest_common::freeWorkspace(temporary_workspace);
    CHECK_CUEST(cuestParametersDestroy(CUEST_PCMDERIVATIVECOMPUTE_PARAMETERS, pcm_derivative_compute_parameters));
    CHECK_CUEST(cuestResultsDestroy(CUEST_PCM_RESULTS, pcm_results));

    // Copy the gradient matrix to the host
    err = cudaMemcpy(grad->get_pointer(0), d_V_gradient_matrix, natom * 3 * sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        throw PSIEXCEPTION("cudaMemcpy failed in cuestPCM::compute_PCM_gradient");
    }
    err = cudaFree(d_V_gradient_matrix);
    if (err != cudaSuccess) {
        throw PSIEXCEPTION("cudaFree failed in cuestPCM::compute_PCM_gradient");
    }

    return grad;
}

std::pair<double, SharedMatrix> cuestPCM::compute_PCM_terms(const SharedMatrix &D) {
    cuestPCMResults_t pcm_results;
    CHECK_CUEST(cuestResultsCreate(CUEST_PCM_RESULTS, &pcm_results));

    // Copy the density matrix to the device
    int nbf = mintshelper_->basisset()->nbf();
    cudaError_t err = cudaMemcpy(d_D_matrix_, D->pointer()[0], nbf * nbf * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        throw PSIEXCEPTION("cudaMemcpy failed in cuestPCM::compute_PCM_terms");
    }
    cuestPCMPotentialComputeParameters_t pcm_potential_compute_parameters;
    CHECK_CUEST(cuestParametersCreate(CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS, &pcm_potential_compute_parameters));
    CHECK_CUEST(cuestParametersConfigure(CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS, pcm_potential_compute_parameters, CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS_CONVERGENCE_THRESHOLD, &pcm_convergence_, sizeof(pcm_convergence_)));
    cuestWorkspaceDescriptor_t temporary_workspace_descriptor{};

    CHECK_CUEST(cuestPCMPotentialComputeWorkspaceQuery(
        cuest_handle,
        pcm_integral_plan_,
        pcm_potential_compute_parameters,
        &temporary_workspace_descriptor,
        d_D_matrix_,
        d_q1_matrix_,
        d_q2_matrix_,
        pcm_results,
        d_V_matrix_
    ));
    cuestWorkspace_t* temporary_workspace = cuest_common::allocateWorkspace(&temporary_workspace_descriptor);
    CHECK_CUEST(cuestPCMPotentialCompute(
        cuest_handle,
        pcm_integral_plan_,
        pcm_potential_compute_parameters,
        temporary_workspace,
        d_D_matrix_,
        d_q1_matrix_,
        d_q2_matrix_,
        pcm_results,
        d_V_matrix_
    ));

    // The output q values should be used as the guess on the next iteration
    std::swap(d_q1_matrix_, d_q2_matrix_);
    cuest_common::freeWorkspace(temporary_workspace);
    CHECK_CUEST(cuestParametersDestroy(CUEST_PCMPOTENTIALCOMPUTE_PARAMETERS, pcm_potential_compute_parameters));
    int32_t converged = 0;
    CHECK_CUEST(cuestResultsQuery(CUEST_PCM_RESULTS, pcm_results, CUEST_PCMRESULT_CONVERGED, &converged, sizeof(converged)));
    if (!converged) {
        throw PSIEXCEPTION("PCM potential computation did not converge.");
    }
    double Epcm = 0.0;
    CHECK_CUEST(cuestResultsQuery(CUEST_PCM_RESULTS, pcm_results, CUEST_PCMRESULT_PCM_DIELECTRIC_ENERGY, &Epcm, sizeof(Epcm)));
    if (!converged) {
        throw PSIEXCEPTION("PCM potential computation did not converge.");
    }
    CHECK_CUEST(cuestResultsDestroy(CUEST_PCM_RESULTS, pcm_results));

    // Copy the potential matrix to the host
    SharedMatrix V = D->clone();
    err = cudaMemcpy(V->pointer()[0], d_V_matrix_, nbf * nbf * sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        throw PSIEXCEPTION("cudaMemcpy failed in cuestPCM::compute_PCM_terms");
    }
    return {Epcm, V};
}


cuestPCM::~cuestPCM() {
    if (pcm_integral_plan_ != nullptr) {
        CHECK_CUEST(cuestPCMIntPlanDestroy(pcm_integral_plan_));
        pcm_integral_plan_ = nullptr;
    }
    if (pcm_integral_ws_ptr_ != nullptr) {
        cuest_common::freeWorkspace(pcm_integral_ws_ptr_);
        pcm_integral_ws_ptr_ = nullptr;
    }
    if (d_q1_matrix_ != nullptr) {
        cudaFree(d_q1_matrix_);
        d_q1_matrix_ = nullptr;
    }
    if (d_q2_matrix_ != nullptr) {
        cudaFree(d_q2_matrix_);
        d_q2_matrix_ = nullptr;
    }
    if (d_D_matrix_ != nullptr) {
        cudaFree(d_D_matrix_);
        d_D_matrix_ = nullptr;
    }
    if (d_V_matrix_ != nullptr) {
        cudaFree(d_V_matrix_);
        d_V_matrix_ = nullptr;
    }
}

SharedMatrix cuestPCM::compute_V(const SharedMatrix &D) {
    return nullptr;
}

} // namespace psi

#endif // USING_cuEST