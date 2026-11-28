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

#ifdef USING_gauxc

#include "gauxc_device.h"

#include <sstream>

#include <gauxc/gauxc_config.hpp>

#if defined(GAUXC_HAS_CUDA) && defined(PSI4_HAS_CUDA_RUNTIME)
#include <cuda_runtime.h>
#include "psi4/libpsi4util/cuda_diagnostics.h"
#endif

namespace psi {
namespace gauxc_device {

std::string check(const std::string& option) {
#ifndef GAUXC_HAS_DEVICE
    throw PSIEXCEPTION(option + " was requested, but GauXC was not built with device support.");
#elif defined(GAUXC_HAS_CUDA) && defined(PSI4_HAS_CUDA_RUNTIME)
    // Before handing the GPU to GauXC, confirm there is one and identify it
    using namespace cuda_diagnostics;
    const std::string requested = option + " requested, but ";
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess) {
        throw PSIEXCEPTION(cuda_failure_message(requested + "CUDA device discovery failed", err, "", "GauXC"));
    }
    if (device_count == 0) {
        throw PSIEXCEPTION(requested + "no CUDA-capable GPU was found. Set " + option +
                           " false to run GauXC on the CPU.");
    }
    int device_id = 0;
    err = cudaGetDevice(&device_id);
    if (err != cudaSuccess) {
        throw PSIEXCEPTION(cuda_failure_message(requested + "cudaGetDevice failed", err, "", "GauXC"));
    }
    cudaDeviceProp props;
    err = cudaGetDeviceProperties(&props, device_id);
    if (err != cudaSuccess) {
        throw PSIEXCEPTION(cuda_failure_message(requested + "cudaGetDeviceProperties failed", err, "", "GauXC"));
    }
    std::ostringstream desc;
    desc << describe_cuda_device(device_id, props) << ", compute capability " << props.major << "." << props.minor;
    return desc.str();
#else
    return std::string();
#endif
}

// GauXC's own device errors (e.g., CUBLAS_STATUS_ARCH_MISMATCH, no kernel image for device) read like
// internal faults. Which architectures a GauXC binary carries isn't exported, so explain the likely causes.
std::string failure_message(const std::string& what, const std::string& option, const std::string& device_desc) {
    std::ostringstream msg;
    msg << "GauXC failed to run on the GPU";
    if (!device_desc.empty()) msg << ": " << device_desc;
    msg << ".\n\n" << what << "\n"
        << "  This is almost always the GPU or the GauXC build, not the input. GauXC runs only on GPU\n"
        << "  architectures it was compiled for (list them with `cuobjdump --list-elf libgauxc.so`;\n"
        << "  CUDA 13 cannot target compute capability below 7.5). Run on a supported GPU, rebuild GauXC\n"
        << "  with this compute capability in CMAKE_CUDA_ARCHITECTURES, or set " << option << " false to\n"
        << "  run GauXC on the CPU.\n";
    return msg.str();
}

}  // namespace gauxc_device
}  // namespace psi

#endif  // USING_gauxc
