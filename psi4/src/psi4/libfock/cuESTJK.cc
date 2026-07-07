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
// The interface to cuEST was contributed by NVIDIA under the following terms:
// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: LGPL-3.0-only

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

#include <cublas_v2.h>
#include <cusolverDn.h>

extern cublasHandle_t cublas_handle;

namespace psi {

cuESTJK::cuESTJK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options)
    : JK(primary),
      options_(options),
      auxiliary_(auxiliary),
      condition_(1.0E-12),
      pq_threshold_(1.0E-14),
      cuest_primary_basis_(nullptr),
      cuest_auxiliary_basis_(nullptr),
      cuest_pair_list_(nullptr),
      cuest_df_plan_(nullptr),
      cuest_coulomb_compute_params_(nullptr),
      cuest_exchange_compute_params_(nullptr),
      initialized_(false) 
{
    cuest_common::ensure_cuest_initialized();

    cuest_primary_basis_ = primary->cuest_basis();
    cuest_auxiliary_basis_ = auxiliary->cuest_basis();

    condition_ = options.get_double("DF_FITTING_CONDITION");
    pq_threshold_ = options.get_double("INTS_TOLERANCE");
}

cuESTJK::~cuESTJK() {
    destroy_cuest_objects();
}

void cuESTJK::preiterations()
{
    auto mol = primary_->molecule();
    size_t natom = mol->natom();

    cuestWorkspaceDescriptor_t* persistentWorkspaceDescriptor = (cuestWorkspaceDescriptor_t*) malloc(sizeof(cuestWorkspaceDescriptor_t));
    cuestWorkspaceDescriptor_t* temporaryWorkspaceDescriptor = (cuestWorkspaceDescriptor_t*) malloc(sizeof(cuestWorkspaceDescriptor_t));

    // AO Pair List
    cuestAOPairListParameters_t pair_params;
    CHECK_CUEST(cuestParametersCreate(CUEST_AOPAIRLIST_PARAMETERS, reinterpret_cast<void**>(&pair_params)));

    CHECK_CUEST(cuestAOPairListCreateWorkspaceQuery(cuest_handle, cuest_primary_basis_,
        static_cast<uint64_t>(natom), mol->geometry().pointer()[0], pq_threshold_, pair_params, persistentWorkspaceDescriptor, temporaryWorkspaceDescriptor, nullptr));

    cuest_pair_list_ws_ptr_ = cuest_common::allocateWorkspace(persistentWorkspaceDescriptor);
    cuestWorkspace_t* temporaryPairListWorkspace = cuest_common::allocateWorkspace(temporaryWorkspaceDescriptor);

    CHECK_CUEST(cuestAOPairListCreate(cuest_handle, cuest_primary_basis_,
        static_cast<uint64_t>(natom), mol->geometry().pointer()[0], pq_threshold_, pair_params, cuest_pair_list_ws_ptr_, temporaryPairListWorkspace, &cuest_pair_list_));

    cuest_common::freeWorkspace(temporaryPairListWorkspace);
    CHECK_CUEST(cuestParametersDestroy(CUEST_AOPAIRLIST_PARAMETERS, pair_params));

    // DF Int Plan List
    cuestDFIntPlanParameters_t dfint_params;
    CHECK_CUEST(cuestParametersCreate(CUEST_DFINTPLAN_PARAMETERS, reinterpret_cast<void**>(&dfint_params)));

    CHECK_CUEST(cuestDFIntPlanCreateWorkspaceQuery(cuest_handle,
        cuest_primary_basis_, cuest_auxiliary_basis_, cuest_pair_list_,
        dfint_params, persistentWorkspaceDescriptor, temporaryWorkspaceDescriptor, nullptr));

    cuest_dfint_plan_ws_ptr_ = cuest_common::allocateWorkspace(persistentWorkspaceDescriptor);;
    cuestWorkspace_t* temporaryDFIntPlanWorkspace = cuest_common::allocateWorkspace(temporaryWorkspaceDescriptor);;

    CHECK_CUEST(cuestDFIntPlanCreate(cuest_handle,
        cuest_primary_basis_, cuest_auxiliary_basis_, cuest_pair_list_,
        dfint_params, cuest_dfint_plan_ws_ptr_, temporaryDFIntPlanWorkspace, &cuest_df_plan_));

    cuest_common::freeWorkspace(temporaryDFIntPlanWorkspace);
    CHECK_CUEST(cuestParametersDestroy(CUEST_DFINTPLAN_PARAMETERS, dfint_params));

    free(persistentWorkspaceDescriptor);
    free(temporaryWorkspaceDescriptor);

    // Create J/K parameters
    CHECK_CUEST(cuestParametersCreate(CUEST_DFCOULOMBCOMPUTE_PARAMETERS, reinterpret_cast<void**>(&cuest_coulomb_compute_params_)));
    CHECK_CUEST(cuestParametersCreate(CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS, reinterpret_cast<void**>(&cuest_exchange_compute_params_)));

    initialized_ = true;
}

void cuESTJK::destroy_cuest_objects() {
    if (cuest_df_plan_) {
        CHECK_CUEST(cuestDFIntPlanDestroy(cuest_df_plan_));
        cuest_df_plan_ = nullptr;
    }
    if (cuest_pair_list_) {
        CHECK_CUEST(cuestAOPairListDestroy(cuest_pair_list_));
        cuest_pair_list_ = nullptr;
    }
    if (cuest_pair_list_ws_ptr_) {
        cuest_common::freeWorkspace(cuest_pair_list_ws_ptr_);
        cuest_pair_list_ws_ptr_ = nullptr;
    }
    if (cuest_dfint_plan_ws_ptr_) {
        cuest_common::freeWorkspace(cuest_dfint_plan_ws_ptr_);
        cuest_dfint_plan_ws_ptr_ = nullptr;
    }
    if (cuest_coulomb_compute_params_) {
        CHECK_CUEST(cuestParametersDestroy(CUEST_DFCOULOMBCOMPUTE_PARAMETERS, cuest_coulomb_compute_params_));
        cuest_coulomb_compute_params_ = nullptr;
    }
    if (cuest_exchange_compute_params_) {
        CHECK_CUEST(cuestParametersDestroy(CUEST_DFSYMMETRICEXCHANGECOMPUTE_PARAMETERS, cuest_exchange_compute_params_));
        cuest_exchange_compute_params_ = nullptr;
    }
}

size_t cuESTJK::memory_estimate() {
    return 0;
}

void cuESTJK::print_header() const {
    if (print_) {
        outfile->Printf("  ==> cuESTJK: GPU-Accelerated Density-Fitted J/K Matrices <==\n\n");
        outfile->Printf("    J tasked:              %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:              %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    wK tasked:             %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) 
        outfile->Printf("    Omega:                 %11.3E\n", omega_);
        outfile->Printf("    Pseudoinverse cutoff:  %11.1E\n", condition_);
        outfile->Printf("    Threshold PQ:          %11.1E\n", pq_threshold_);
        outfile->Printf("\n");
    }
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
 
    cuestWorkspaceDescriptor_t* j_temp_desc = (cuestWorkspaceDescriptor_t*) malloc(sizeof(cuestWorkspaceDescriptor_t));
    cuestWorkspaceDescriptor_t* k_temp_desc = (cuestWorkspaceDescriptor_t*) malloc(sizeof(cuestWorkspaceDescriptor_t));

    cuestWorkspaceDescriptor_t* variableBufferSize = (cuestWorkspaceDescriptor_t*) malloc(sizeof(cuestWorkspaceDescriptor_t));
    variableBufferSize->hostBufferSizeInBytes = 0;
    variableBufferSize->deviceBufferSizeInBytes = 8000000000;
    
    size_t max_host = 0, max_device = 0;
 
    if (do_J_) {
        CHECK_CUEST(cuestDFCoulombComputeWorkspaceQuery(
            cuest_handle,
            cuest_df_plan_,
            cuest_coulomb_compute_params_,
            j_temp_desc,
            d_D,
            d_J));
        max_host = std::max(max_host, j_temp_desc->hostBufferSizeInBytes);
        max_device = std::max(max_device, j_temp_desc->deviceBufferSizeInBytes);
    }
 
    if (do_K_) {
        for (size_t N = 0; N < D_ao_.size(); N++) {
            int nocc = C_left_ao_[N]->ncol();
            if (nocc == 0) continue;
 
            CHECK_CUEST(cuestDFSymmetricExchangeComputeWorkspaceQuery(
                cuest_handle,
                cuest_df_plan_,
                cuest_exchange_compute_params_,
                variableBufferSize,
                k_temp_desc,
                static_cast<uint64_t>(nocc),
                nullptr,
                d_K));
            max_host = std::max(max_host, k_temp_desc->hostBufferSizeInBytes);
            max_device = std::max(max_device, k_temp_desc->deviceBufferSizeInBytes);
        }
    }
 
    cuestWorkspaceDescriptor_t* total_desc = (cuestWorkspaceDescriptor_t*) malloc(sizeof(cuestWorkspaceDescriptor_t));
    total_desc->hostBufferSizeInBytes = max_host;
    total_desc->deviceBufferSizeInBytes = max_device;
 
    cuestWorkspace_t* temporaryJKWorkspace = cuest_common::allocateWorkspace(total_desc);

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
                cuest_coulomb_compute_params_,
                temporaryJKWorkspace,
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
                    cuest_exchange_compute_params_,
                    variableBufferSize,
                    temporaryJKWorkspace,
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
 
    free(j_temp_desc);
    free(k_temp_desc);
    free(total_desc);
    free(variableBufferSize);

    cuest_common::freeWorkspace(temporaryJKWorkspace);

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

void cuESTJK::postiterations() 
{
    destroy_cuest_objects();
}

} // namespace psi

#endif // USING_cuEST
