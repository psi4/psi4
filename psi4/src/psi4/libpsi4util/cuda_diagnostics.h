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

#ifndef _psi_src_lib_libpsi4util_cuda_diagnostics_h_
#define _psi_src_lib_libpsi4util_cuda_diagnostics_h_

/*! \file cuda_diagnostics.h
 *  \brief Readable explanations of CUDA failures, shared by every GPU-using component.
 *
 *  Only available when psi4 links the CUDA runtime (PSI4_HAS_CUDA_RUNTIME), e.g., through
 *  cuEST or a device-enabled GauXC.
 */

#ifdef PSI4_HAS_CUDA_RUNTIME

#include <string>

#include <cuda_runtime.h>

namespace psi {
namespace cuda_diagnostics {

/*! For CUDA status codes that mean "this GPU or this node is broken/misprovisioned" rather
 *  than "the program did something wrong", a diagnosis plus the command that confirms it.
 *  Empty string for codes that are not known environmental faults.
 */
std::string cuda_unhealthy_environment_hint(cudaError_t err);

/*! Text for a failed CUDA call. `what` is the ordinary human-readable description; when the
 *  status code is a known environmental fault, the health diagnosis is appended. `component`
 *  (e.g., "cuEST", "GauXC") is named alongside Psi4 as not being at fault. Returns rather than
 *  throws so that PSIEXCEPTION still records the caller's __FILE__/__LINE__.
 */
std::string cuda_failure_message(const std::string& what, cudaError_t err, const std::string& device_desc = "",
                                 const std::string& component = "");

/*! Identify the device well enough to quote in a ticket: logical index, model name, and the
 *  PCI bus ID that nvidia-smi and the kernel log also report.
 */
std::string describe_cuda_device(int device_id, const cudaDeviceProp& props);

}  // namespace cuda_diagnostics
}  // namespace psi

#endif  // PSI4_HAS_CUDA_RUNTIME

#endif  // _psi_src_lib_libpsi4util_cuda_diagnostics_h_
