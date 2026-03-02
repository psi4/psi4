#include "cuESTJK.h"

#ifdef USING_cuEST

#include "cuESTCommon.h"

#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/integral.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/psi4-dec.h"

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace psi {

cuESTJK::cuESTJK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options)
    : JK(primary),
      options_(options),
      auxiliary_(auxiliary),
      condition_(1.0E-12),
      cuest_primary_basis_(nullptr),
      cuest_auxiliary_basis_(nullptr),
      cuest_pair_list_(nullptr),
      cuest_df_plan_(nullptr),
      plan_built_(false) {
    primary_persistent_ws_ = {};
    auxiliary_persistent_ws_ = {};
    pair_persistent_ws_ = {};
    df_persistent_ws_ = {};
    compute_temp_ws_ = {};
    exchange_max_ws_desc_ = {};
}

cuESTJK::~cuESTJK() {
    destroy_cuest_objects();
}

void cuESTJK::allocate_workspace(cuestWorkspaceDescriptor_t& desc, cuestWorkspace_t& ws) {
    cuest_common::alloc_workspace(desc, ws);
}

void cuESTJK::free_workspace(cuestWorkspace_t& ws) {
    cuest_common::free_workspace(ws);
}

cuestAOBasis_t cuESTJK::build_cuest_basis(std::shared_ptr<BasisSet> basis,
                                           std::vector<cuestAOShell_t>& shells_out,
                                           cuestWorkspace_t& persistent_ws) {
    return cuest_common::build_cuest_basis(basis, shells_out, persistent_ws);
}

void cuESTJK::destroy_cuest_objects() {
    if (cuest_df_plan_) {
        cuestDFIntPlanDestroy(cuest_df_plan_);
        cuest_df_plan_ = nullptr;
    }
    if (cuest_pair_list_) {
        cuestAOPairListDestroy(cuest_pair_list_);
        cuest_pair_list_ = nullptr;
    }
    if (cuest_auxiliary_basis_) {
        cuestAOBasisDestroy(cuest_auxiliary_basis_);
        cuest_auxiliary_basis_ = nullptr;
    }
    if (cuest_primary_basis_) {
        cuestAOBasisDestroy(cuest_primary_basis_);
        cuest_primary_basis_ = nullptr;
    }

    for (auto& s : cuest_primary_shells_) cuestAOShellDestroy(s);
    cuest_primary_shells_.clear();
    for (auto& s : cuest_auxiliary_shells_) cuestAOShellDestroy(s);
    cuest_auxiliary_shells_.clear();

    free_workspace(df_persistent_ws_);
    free_workspace(pair_persistent_ws_);
    free_workspace(auxiliary_persistent_ws_);
    free_workspace(primary_persistent_ws_);
    free_workspace(compute_temp_ws_);

    cuestParametersDestroy(CUEST_DFCOULOMBCOMPUTE_PARAMETERS, coulomb_compute_params_);
    cuestParametersDestroy(CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS, exchange_compute_params_);

    plan_built_ = false;
}

size_t cuESTJK::memory_estimate() {
    return 0;
}

void cuESTJK::print_header() const {
    if (print_) {
        outfile->Printf("  ==> cuESTJK: GPU-Accelerated Density-Fitted J/K Matrices <==\n\n");
        outfile->Printf("    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) outfile->Printf("    Omega:             %11.3E\n", omega_);
        outfile->Printf("    Fitting Coverage:  %11.0E\n", condition_);
        outfile->Printf("\n");
    }
}

void cuESTJK::preiterations() {
    if (plan_built_) return;

    cuest_primary_basis_ = build_cuest_basis(primary_, cuest_primary_shells_, primary_persistent_ws_);
    cuest_auxiliary_basis_ = build_cuest_basis(auxiliary_, cuest_auxiliary_shells_, auxiliary_persistent_ws_);

    auto mol = primary_->molecule();
    int natom = mol->natom();

    std::vector<double> xyz(natom * 3);
    for (int A = 0; A < natom; A++) {
        xyz[3 * A + 0] = mol->x(A);
        xyz[3 * A + 1] = mol->y(A);
        xyz[3 * A + 2] = mol->z(A);
    }

    cuestAOPairListParameters_t pair_params;
    CHECK_CUEST(cuestParametersCreate(CUEST_AOPAIRLIST_PARAMETERS, reinterpret_cast<void**>(&pair_params)));

    cuestWorkspaceDescriptor_t pair_persistent_desc = {}, pair_temp_desc = {};
    CHECK_CUEST(cuestAOPairListCreateWorkspaceQuery(
        cuest_handle, cuest_primary_basis_, static_cast<uint64_t>(natom), xyz.data(),
        cutoff_, pair_params, &pair_persistent_desc, &pair_temp_desc, nullptr));

    allocate_workspace(pair_persistent_desc, pair_persistent_ws_);

    cuestWorkspace_t pair_temp_ws = {};
    allocate_workspace(pair_temp_desc, pair_temp_ws);

    CHECK_CUEST(cuestAOPairListCreate(
        cuest_handle, cuest_primary_basis_, static_cast<uint64_t>(natom), xyz.data(),
        cutoff_, pair_params,
        &pair_persistent_ws_, &pair_temp_ws, &cuest_pair_list_));

    free_workspace(pair_temp_ws);
    cuestParametersDestroy(CUEST_AOPAIRLIST_PARAMETERS, pair_params);

    cuestDFIntPlanParameters_t df_params;
    CHECK_CUEST(cuestParametersCreate(CUEST_DFINTPLAN_PARAMETERS, reinterpret_cast<void**>(&df_params)));

    CHECK_CUEST(cuestParametersConfigure(CUEST_DFINTPLAN_PARAMETERS, df_params,
        CUEST_DFINTPLAN_PARAMETERS_FITTING_CUTOFF, &condition_, sizeof(double)));

    cuestWorkspaceDescriptor_t df_persistent_desc = {}, df_temp_desc = {};
    CHECK_CUEST(cuestDFIntPlanCreateWorkspaceQuery(
        cuest_handle, cuest_primary_basis_, cuest_auxiliary_basis_,
        cuest_pair_list_, df_params,
        &df_persistent_desc, &df_temp_desc, nullptr));

    allocate_workspace(df_persistent_desc, df_persistent_ws_);

    cuestWorkspace_t df_temp_ws = {};
    allocate_workspace(df_temp_desc, df_temp_ws);

    CHECK_CUEST(cuestDFIntPlanCreate(
        cuest_handle, cuest_primary_basis_, cuest_auxiliary_basis_,
        cuest_pair_list_, df_params,
        &df_persistent_ws_, &df_temp_ws, &cuest_df_plan_));

    free_workspace(df_temp_ws);
    cuestParametersDestroy(CUEST_DFINTPLAN_PARAMETERS, df_params);

    exchange_max_ws_desc_.hostBufferSizeInBytes = 0;
    exchange_max_ws_desc_.deviceBufferSizeInBytes = static_cast<size_t>(2) * 1024 * 1024 * 1024;

    CHECK_CUEST(cuestParametersCreate(
        CUEST_DFCOULOMBCOMPUTE_PARAMETERS,
        &coulomb_compute_params_));
    CHECK_CUEST(cuestParametersCreate(
        CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS,
        &exchange_compute_params_));
    plan_built_ = true;
}

void cuESTJK::compute_JK() {
    using clock = std::chrono::high_resolution_clock;
    auto t_total_start = clock::now();

    int nbf = primary_->nbf();
    size_t nbf2_bytes = static_cast<size_t>(nbf) * nbf * sizeof(double);

    auto t0 = clock::now();
    double* d_D = nullptr;
    double* d_J = nullptr;
    double* d_K = nullptr;
    double* d_C = nullptr;

    cudaMalloc(reinterpret_cast<void**>(&d_D), nbf2_bytes);
    cudaMalloc(reinterpret_cast<void**>(&d_J), nbf2_bytes);
    cudaMalloc(reinterpret_cast<void**>(&d_K), nbf2_bytes);
    auto t_alloc = clock::now();

    cuestWorkspaceDescriptor_t j_temp_desc = {};
    cuestWorkspaceDescriptor_t k_temp_desc = {};
    size_t max_host = 0, max_device = 0;

    if (do_J_) {
        CHECK_CUEST(cuestDFCoulombComputeWorkspaceQuery(
            cuest_handle,
            cuest_df_plan_,
            coulomb_compute_params_,
            &j_temp_desc,
            d_D,
            d_J));
        max_host = std::max(max_host, j_temp_desc.hostBufferSizeInBytes);
        max_device = std::max(max_device, j_temp_desc.deviceBufferSizeInBytes);
    }

    if (do_K_) {
        for (size_t N = 0; N < D_ao_.size(); N++) {
            int nocc = C_left_ao_[N]->ncol();
            if (nocc == 0) continue;

            CHECK_CUEST(cuestDFSymmetricExchangeComputeWorkspaceQuery(
                cuest_handle,
                cuest_df_plan_,
                exchange_compute_params_,
                &exchange_max_ws_desc_,
                &k_temp_desc,
                static_cast<uint64_t>(nocc),
                nullptr,
                d_K));
            max_host = std::max(max_host, k_temp_desc.hostBufferSizeInBytes);
            max_device = std::max(max_device, k_temp_desc.deviceBufferSizeInBytes);
        }
    }

    cuestWorkspaceDescriptor_t total_desc;
    total_desc.hostBufferSizeInBytes = max_host;
    total_desc.deviceBufferSizeInBytes = max_device;

    free_workspace(compute_temp_ws_);
    allocate_workspace(total_desc, compute_temp_ws_);
    auto t_ws = clock::now();

    double ms_J_compute = 0.0, ms_K_compute = 0.0;
    double ms_memcpy_h2d = 0.0, ms_memcpy_d2h = 0.0;
    double ms_transpose = 0.0;

    for (size_t N = 0; N < D_ao_.size(); N++) {
        if (do_J_) {
            auto tc0 = clock::now();
            cudaMemcpy(d_D, D_ao_[N]->get_pointer(), nbf2_bytes, cudaMemcpyHostToDevice);
            auto tc1 = clock::now();
            CHECK_CUEST(cuestDFCoulombCompute(
                cuest_handle,
                cuest_df_plan_,
                coulomb_compute_params_,
                &compute_temp_ws_,
                d_D,
                d_J));
            cudaDeviceSynchronize();
            auto tc2 = clock::now();
            cudaMemcpy(J_ao_[N]->get_pointer(), d_J, nbf2_bytes, cudaMemcpyDeviceToHost);
            auto tc3 = clock::now();
            ms_memcpy_h2d += std::chrono::duration<double, std::milli>(tc1 - tc0).count();
            ms_J_compute += std::chrono::duration<double, std::milli>(tc2 - tc1).count();
            ms_memcpy_d2h += std::chrono::duration<double, std::milli>(tc3 - tc2).count();
        }

        if (do_K_) {
            int nocc = C_left_ao_[N]->ncol();
            if (nocc > 0) {
                size_t c_bytes = static_cast<size_t>(nocc) * nbf * sizeof(double);

                auto tt0 = clock::now();
                std::vector<double> C_row_major(nocc * nbf);
                double** Cp = C_left_ao_[N]->pointer();
                for (int i = 0; i < nocc; i++) {
                    for (int mu = 0; mu < nbf; mu++) {
                        C_row_major[i * nbf + mu] = Cp[mu][i];
                    }
                }
                auto tt1 = clock::now();
                ms_transpose += std::chrono::duration<double, std::milli>(tt1 - tt0).count();

                cudaFree(d_C);
                cudaMalloc(reinterpret_cast<void**>(&d_C), c_bytes);

                auto tk0 = clock::now();
                cudaMemcpy(d_C, C_row_major.data(), c_bytes, cudaMemcpyHostToDevice);
                auto tk1 = clock::now();

                CHECK_CUEST(cuestDFSymmetricExchangeCompute(
                    cuest_handle,
                    cuest_df_plan_,
                    exchange_compute_params_,
                    &exchange_max_ws_desc_,
                    &compute_temp_ws_,
                    static_cast<uint64_t>(nocc),
                    d_C,
                    d_K));
                cudaDeviceSynchronize();
                auto tk2 = clock::now();

                cudaMemcpy(K_ao_[N]->get_pointer(), d_K, nbf2_bytes, cudaMemcpyDeviceToHost);
                auto tk3 = clock::now();

                ms_memcpy_h2d += std::chrono::duration<double, std::milli>(tk1 - tk0).count();
                ms_K_compute += std::chrono::duration<double, std::milli>(tk2 - tk1).count();
                ms_memcpy_d2h += std::chrono::duration<double, std::milli>(tk3 - tk2).count();
            }
        }
    }

    auto t_free_start = clock::now();
    cudaFree(d_D);
    cudaFree(d_J);
    cudaFree(d_K);
    if (d_C) cudaFree(d_C);
    auto t_free_end = clock::now();

    double ms_alloc = std::chrono::duration<double, std::milli>(t_alloc - t0).count();
    double ms_wsquery = std::chrono::duration<double, std::milli>(t_ws - t_alloc).count();
    double ms_free = std::chrono::duration<double, std::milli>(t_free_end - t_free_start).count();
    double ms_total = std::chrono::duration<double, std::milli>(t_free_end - t_total_start).count();

    outfile->Printf("    cuESTJK compute_JK: total=%7.2fms | alloc=%5.2fms ws=%5.2fms "
                     "J=%6.2fms K=%6.2fms memcpy(H2D)=%5.2fms memcpy(D2H)=%5.2fms "
                     "transpose=%5.2fms free=%5.2fms\n",
                     ms_total, ms_alloc, ms_wsquery, ms_J_compute, ms_K_compute,
                     ms_memcpy_h2d, ms_memcpy_d2h, ms_transpose, ms_free);
}

void cuESTJK::postiterations() {

    destroy_cuest_objects();
}

} // namespace psi

#endif // USING_cuEST
