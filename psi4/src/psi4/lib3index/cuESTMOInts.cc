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

#include "cuESTMOInts.h"

#ifdef USING_cuEST

#include <cuda_runtime.h>

#include "psi4/libfock/cuESTCommon.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/psi4-dec.h"

namespace psi {

cuESTMOInts::cuESTMOInts(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options, bool integral_direct) 
    : primary_(primary), 
      auxiliary_(auxiliary), 
      pq_threshold_(1.0E-14), 
      integral_direct_(integral_direct ? 1 : 0) 
{
    if (!primary_ || !auxiliary_) {
        throw PSIEXCEPTION("cuESTMOInts: primary and auxiliary basis sets must be non-null.");
    }
    if (!auxiliary_->has_puream()) {
        throw PSIEXCEPTION("cuESTMOInts requires a spherical-harmonic auxiliary basis set; '" + auxiliary_->name() + "' is Cartesian.");
    }
    
    pq_threshold_ = options.get_double("INTS_TOLERANCE");
    common_init();
}

cuESTMOInts::~cuESTMOInts() {
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
}

void cuESTMOInts::common_init() {
    // basis handles
    cuest_primary_basis_ = primary_->cuest_basis();
    cuest_auxiliary_basis_ = auxiliary_->cuest_basis();

    cuestWorkspaceDescriptor_t persistentWorkspaceDescriptor;
    cuestWorkspaceDescriptor_t temporaryWorkspaceDescriptor;

    // AO Pair List
    auto mol = primary_->molecule();
    size_t natom = mol->natom();

    cuestAOPairListParameters_t pair_params;
    CHECK_CUEST(cuestParametersCreate(
        CUEST_AOPAIRLIST_PARAMETERS,
        reinterpret_cast<void**>(&pair_params)));

    CHECK_CUEST(cuestAOPairListCreateWorkspaceQuery(
        cuest_handle,
        cuest_primary_basis_,
        static_cast<uint64_t>(natom),
        mol->geometry().pointer()[0],
        pq_threshold_,
        pair_params,
        &persistentWorkspaceDescriptor,
        &temporaryWorkspaceDescriptor,
        nullptr));

    outfile->Printf("    cuEST AO pair list workspace: %8.3f GiB persistent, %8.3f GiB temporary\n",
                    persistentWorkspaceDescriptor.deviceBufferSizeInBytes / (1024.0 * 1024.0 * 1024.0),
                    temporaryWorkspaceDescriptor.deviceBufferSizeInBytes / (1024.0 * 1024.0 * 1024.0));

    cuest_pair_list_ws_ptr_ = cuest_common::allocateWorkspace(&persistentWorkspaceDescriptor);
    cuestWorkspace_t* temporaryPairListWorkspace = cuest_common::allocateWorkspace(&temporaryWorkspaceDescriptor);

    CHECK_CUEST(cuestAOPairListCreate(
        cuest_handle,
        cuest_primary_basis_,
        static_cast<uint64_t>(natom),
        mol->geometry().pointer()[0],
        pq_threshold_,
        pair_params,
        cuest_pair_list_ws_ptr_,
        temporaryPairListWorkspace,
        &cuest_pair_list_));

    cuest_common::freeWorkspace(temporaryPairListWorkspace);
    CHECK_CUEST(cuestParametersDestroy(
        CUEST_AOPAIRLIST_PARAMETERS,
        pair_params));

    cuestDFIntPlanParameters_t dfint_params;
    CHECK_CUEST(cuestParametersCreate(
        CUEST_DFINTPLAN_PARAMETERS,
        reinterpret_cast<void**>(&dfint_params)));

    int32_t integral_direct = integral_direct_;
    CHECK_CUEST(cuestParametersConfigure(
        CUEST_DFINTPLAN_PARAMETERS,
        dfint_params,
        CUEST_DFINTPLAN_PARAMETERS_THREE_INDEX_INTEGRAL_DIRECT,
        &integral_direct,
        sizeof(int32_t)));

    CHECK_CUEST(cuestDFIntPlanCreateWorkspaceQuery(
        cuest_handle,
        cuest_primary_basis_,
        cuest_auxiliary_basis_,
        cuest_pair_list_,
        dfint_params,
        &persistentWorkspaceDescriptor,
        &temporaryWorkspaceDescriptor,
        nullptr));

    outfile->Printf("    cuEST DF plan workspace:      %8.3f GiB persistent, %8.3f GiB temporary\n",
                    persistentWorkspaceDescriptor.deviceBufferSizeInBytes / (1024.0 * 1024.0 * 1024.0),
                    temporaryWorkspaceDescriptor.deviceBufferSizeInBytes / (1024.0 * 1024.0 * 1024.0));

    cuest_dfint_plan_ws_ptr_ = cuest_common::allocateWorkspace(&persistentWorkspaceDescriptor);
    cuestWorkspace_t* temporaryDFIntPlanWorkspace = cuest_common::allocateWorkspace(&temporaryWorkspaceDescriptor);

    CHECK_CUEST(cuestDFIntPlanCreate(
        cuest_handle,
        cuest_primary_basis_,
        cuest_auxiliary_basis_,
        cuest_pair_list_,
        dfint_params,
        cuest_dfint_plan_ws_ptr_,
        temporaryDFIntPlanWorkspace,
        &cuest_df_plan_));

    cuest_common::freeWorkspace(temporaryDFIntPlanWorkspace);
    CHECK_CUEST(cuestParametersDestroy(
        CUEST_DFINTPLAN_PARAMETERS,
        dfint_params));
}

void cuESTMOInts::compute_3c_ints(SharedMatrix Cleft, SharedMatrix Cright, double* d_out) {
    if (Cleft->nirrep() != 1 || Cright->nirrep() != 1) {
        throw PSIEXCEPTION("cuESTMOInts::compute_3c_ints: coefficient matrices must be C1 (AO basis).");
    }

    const size_t nao = primary_->nbf();
    const uint64_t nleft = static_cast<uint64_t>(Cleft->colspi()[0]);
    const uint64_t nright = static_cast<uint64_t>(Cright->colspi()[0]);

    if (static_cast<size_t>(Cleft->rowspi()[0]) != nao || static_cast<size_t>(Cright->rowspi()[0]) != nao) {
        throw PSIEXCEPTION("cuESTMOInts::compute_3c_ints: coefficient matrices must have nbf rows.");
    }
    if (nleft == 0 || nright == 0) {
        throw PSIEXCEPTION("cuESTMOInts::compute_3c_ints: cuEST requires at least one left and one right orbital.");
    }

    // cuEST needs Cleft & Cright as [nmo, nao]
    SharedMatrix CleftT = Cleft->transpose();
    SharedMatrix CrightT = Cright->transpose();

    double* d_Cleft = nullptr;
    double* d_Cright = nullptr;

    // TODO: add memory estimate for these calls
    CHECK_CUDA(cudaMalloc((void**)&d_Cleft, sizeof(double) * nleft * nao));
    CHECK_CUDA(cudaMalloc((void**)&d_Cright, sizeof(double) * nright * nao));
    CHECK_CUDA(cudaMemcpy(d_Cleft, CleftT->pointer()[0], sizeof(double) * nleft * nao, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_Cright, CrightT->pointer()[0], sizeof(double) * nright * nao, cudaMemcpyHostToDevice));

    // presently there are no configurable parameters for MO ints
    cuestDFMOIntegralsComputeParameters_t params;
    CHECK_CUEST(cuestParametersCreate(CUEST_DFMOINTEGRALSCOMPUTE_PARAMETERS, reinterpret_cast<void**>(&params)));

    // internal scratch for the transform, 2 GB is a reasonable default.
    // this does not limit total (on device) memory usage
    cuestWorkspaceDescriptor_t variableBufferSize;
    variableBufferSize.hostBufferSizeInBytes = 0;
    variableBufferSize.deviceBufferSizeInBytes = 2000000000;

    cuestWorkspaceDescriptor_t temporaryWorkspaceDescriptor;

    // a single (left, right) pair -> a single [naux, nleft, nright] tensor
    const uint64_t num_matrices = 1;

    CHECK_CUEST(cuestDFMOIntegralsComputeWorkspaceQuery(
        cuest_handle, 
        cuest_df_plan_, 
        params, 
        &variableBufferSize,
        &temporaryWorkspaceDescriptor, 
        num_matrices, 
        &nleft, 
        &nright,
        d_Cleft, 
        d_Cright, 
        d_out));

    outfile->Printf("    cuEST AO->MO workspace:       %8.3f GiB temporary\n",
                    temporaryWorkspaceDescriptor.deviceBufferSizeInBytes / (1024.0 * 1024.0 * 1024.0));

    cuestWorkspace_t* temporaryDFMOWorkspace = cuest_common::allocateWorkspace(&temporaryWorkspaceDescriptor);

    CHECK_CUEST(cuestDFMOIntegralsCompute(
        cuest_handle, 
        cuest_df_plan_, 
        params, 
        &variableBufferSize,
        temporaryDFMOWorkspace, 
        num_matrices, 
        &nleft, 
        &nright, 
        d_Cleft, 
        d_Cright,
        d_out));

    cuest_common::freeWorkspace(temporaryDFMOWorkspace);
    CHECK_CUEST(cuestParametersDestroy(CUEST_DFMOINTEGRALSCOMPUTE_PARAMETERS, params));

    CHECK_CUDA(cudaFree(d_Cleft));
    CHECK_CUDA(cudaFree(d_Cright));
}

}  // namespace psi

#endif  // USING_cuEST
