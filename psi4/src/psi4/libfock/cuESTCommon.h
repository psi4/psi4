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

inline void alloc_workspace(cuestWorkspaceDescriptor_t& desc, cuestWorkspace_t& ws) {
    ws = {};
    if (desc.hostBufferSizeInBytes > 0) {
        ws.hostBuffer = reinterpret_cast<uintptr_t>(malloc(desc.hostBufferSizeInBytes));
        ws.hostBufferSizeInBytes = desc.hostBufferSizeInBytes;
    }
    if (desc.deviceBufferSizeInBytes > 0) {
        void* dev_ptr = nullptr;
        cudaMalloc(&dev_ptr, desc.deviceBufferSizeInBytes);
        ws.deviceBuffer = reinterpret_cast<uintptr_t>(dev_ptr);
        ws.deviceBufferSizeInBytes = desc.deviceBufferSizeInBytes;
    }
}

inline void free_workspace(cuestWorkspace_t& ws) {
    if (ws.hostBuffer) {
        free(reinterpret_cast<void*>(ws.hostBuffer));
        ws.hostBuffer = 0;
        ws.hostBufferSizeInBytes = 0;
    }
    if (ws.deviceBuffer) {
        cudaFree(reinterpret_cast<void*>(ws.deviceBuffer));
        ws.deviceBuffer = 0;
        ws.deviceBufferSizeInBytes = 0;
    }
}

inline cuestAOBasis_t build_cuest_basis(std::shared_ptr<BasisSet> basis,
                                        cuestWorkspace_t& persistent_ws) {
    auto mol = basis->molecule();
    int natom = mol->natom();

    cuestAOShellParameters_t shell_params;
    CHECK_CUEST(cuestParametersCreate(CUEST_AOSHELL_PARAMETERS, reinterpret_cast<void**>(&shell_params)));

    std::vector<cuestAOShell_t> shells_out;
    std::vector<uint64_t> shells_per_atom(natom);

    for (int A = 0; A < natom; A++) {
        int nshell_on_atom = basis->nshell_on_center(A);
        shells_per_atom[A] = static_cast<uint64_t>(nshell_on_atom);
        for (int Q = 0; Q < nshell_on_atom; Q++) {
            int shell_idx = basis->shell_on_center(A, Q);
            const auto& shell = basis->shell(shell_idx);
            cuestAOShell_t cuest_shell;
            CHECK_CUEST(cuestAOShellCreate(cuest_handle, shell.is_pure() ? 1 : 0,
                                           static_cast<uint64_t>(shell.am()),
                                           static_cast<uint64_t>(shell.nprimitive()),
                                           shell.exps(), shell.coefs(), shell_params, &cuest_shell));
            shells_out.push_back(cuest_shell);
        }
    }
    cuestParametersDestroy(CUEST_AOSHELL_PARAMETERS, shell_params);

    cuestAOBasisParameters_t basis_params;
    CHECK_CUEST(cuestParametersCreate(CUEST_AOBASIS_PARAMETERS, reinterpret_cast<void**>(&basis_params)));

    cuestWorkspaceDescriptor_t persistent_desc = {}, temp_desc = {};
    CHECK_CUEST(cuestAOBasisCreateWorkspaceQuery(cuest_handle, static_cast<uint64_t>(natom),
        shells_per_atom.data(), shells_out.data(), basis_params, &persistent_desc, &temp_desc, nullptr));

    alloc_workspace(persistent_desc, persistent_ws);
    cuestWorkspace_t temp_ws = {};
    alloc_workspace(temp_desc, temp_ws);

    cuestAOBasis_t cuest_basis;
    CHECK_CUEST(cuestAOBasisCreate(cuest_handle, static_cast<uint64_t>(natom),
        shells_per_atom.data(), shells_out.data(), basis_params, &persistent_ws, &temp_ws, &cuest_basis));

    free_workspace(temp_ws);
    cuestParametersDestroy(CUEST_AOBASIS_PARAMETERS, basis_params);

    for (auto& s : shells_out) cuestAOShellDestroy(s);
    shells_out.clear();

    return cuest_basis;
}

inline cuestAOPairList_t build_cuest_pairlist(cuestAOBasis_t basis, int natom,
                                              const double* xyz, double cutoff,
                                              cuestWorkspace_t& persistent_ws) {
    cuestAOPairListParameters_t pair_params;
    CHECK_CUEST(cuestParametersCreate(CUEST_AOPAIRLIST_PARAMETERS, reinterpret_cast<void**>(&pair_params)));

    cuestWorkspaceDescriptor_t p_desc = {}, t_desc = {};
    CHECK_CUEST(cuestAOPairListCreateWorkspaceQuery(cuest_handle, basis,
        static_cast<uint64_t>(natom), xyz, cutoff, pair_params, &p_desc, &t_desc, nullptr));

    alloc_workspace(p_desc, persistent_ws);
    cuestWorkspace_t temp_ws = {};
    alloc_workspace(t_desc, temp_ws);

    cuestAOPairList_t pair_list;
    CHECK_CUEST(cuestAOPairListCreate(cuest_handle, basis,
        static_cast<uint64_t>(natom), xyz, cutoff, pair_params, &persistent_ws, &temp_ws, &pair_list));

    free_workspace(temp_ws);
    cuestParametersDestroy(CUEST_AOPAIRLIST_PARAMETERS, pair_params);
    return pair_list;
}

}  // namespace cuest_common
}  // namespace psi

#endif  // USING_cuEST
#endif  // CUEST_COMMON_H
