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
