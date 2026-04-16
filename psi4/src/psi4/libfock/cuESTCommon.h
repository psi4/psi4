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

#ifndef CUEST_COMMON_H
#define CUEST_COMMON_H

#ifdef USING_cuEST

#include <cuest.h>
#include <cuda_runtime.h>

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libpsi4util/exception.h"

#include <cstdlib>
#include <memory>
#include <sstream>
#include <vector>

extern cuestHandle_t cuest_handle;

namespace psi {
namespace cuest_common {

inline void check_cuest(cuestStatus_t status, const char* func) {
    if (status != CUEST_STATUS_SUCCESS) {
        std::ostringstream msg;
        msg << "cuEST error in " << func << " (status code " << static_cast<int>(status) << ")";
        throw PsiException(msg.str(), __FILE__, __LINE__);
    }
}

#define CHECK_CUEST(call) ::psi::cuest_common::check_cuest((call), #call)

inline cuestWorkspace_t* allocateWorkspace(const cuestWorkspaceDescriptor_t* workspaceDescriptor)
{
    /* Check that a valid workspace descriptor has been provided. */
    if (workspaceDescriptor == NULL) {
        fprintf(stderr, "Invalid argument: workspaceDescriptor must not be NULL\n");
        exit(EXIT_FAILURE);
    }

    /* Allocate the workspace structure. */
    cuestWorkspace_t* workspace = (cuestWorkspace_t*) malloc(sizeof(cuestWorkspace_t));
    if (!workspace) {
        fprintf(stderr, "Failed to allocate cuestWorkspace_t struct\n");
        exit(EXIT_FAILURE);
    }

    /* Set the length of the host and device buffers. */
    workspace->hostBufferSizeInBytes   = workspaceDescriptor->hostBufferSizeInBytes;
    workspace->deviceBufferSizeInBytes = workspaceDescriptor->deviceBufferSizeInBytes;
    workspace->hostBuffer              = (uintptr_t) NULL;
    workspace->deviceBuffer            = (uintptr_t) NULL;

    /* Allocate the host buffer if the size is non-zero. */
    if (workspace->hostBufferSizeInBytes) {
        void* hostPtr = (void*) malloc(workspace->hostBufferSizeInBytes);
        if (!hostPtr) {
            fprintf(stderr, "Failed to allocate host buffer\n");
            free(workspace);
            exit(EXIT_FAILURE);
        }
        workspace->hostBuffer = (uintptr_t) hostPtr;
    }

    /* Allocate the device buffer if the size is non-zero. */
    if (workspace->deviceBufferSizeInBytes) {
        void* devicePtr = NULL;
        cudaError_t err = cudaMalloc(&devicePtr, workspace->deviceBufferSizeInBytes);
        if (err != cudaSuccess) {
            fprintf(stderr, "Failed to allocate device buffer: %s\n", cudaGetErrorString(err));
            if (workspace->hostBuffer) free((void*) workspace->hostBuffer);
            free(workspace);
            exit(EXIT_FAILURE);
        }
        workspace->deviceBuffer = (uintptr_t) devicePtr;
    }

    /* Return the fully allocated workspace. */
    return workspace;
}

inline void freeWorkspace(cuestWorkspace_t* workspace)
{
    if (!workspace) return;

    /* Frss the host buffer if it is not NULL. */
    if (workspace->hostBuffer) {
        free((void*) workspace->hostBuffer);
    }

    /* Frss the device buffer if it is not NULL. */
    if (workspace->deviceBuffer) {
        cudaError_t err = cudaFree((void*) workspace->deviceBuffer);
        if (err != cudaSuccess) {
            fprintf(stderr, "cudaFree failed: %s\n", cudaGetErrorString(err));
            free(workspace);
            exit(EXIT_FAILURE);
        }
    }

    /* Free the workspace structure. */
    free(workspace);
}

}  // namespace cuest_common
}  // namespace psi

#endif  // USING_cuEST
#endif  // CUEST_COMMON_H
