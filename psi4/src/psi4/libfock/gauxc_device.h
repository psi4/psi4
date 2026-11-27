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

#ifndef LIBFOCK_GAUXC_DEVICE_H
#define LIBFOCK_GAUXC_DEVICE_H

/*! \file gauxc_device.h
 *  \brief Readable failures for GauXC GPU requests, shared by the DFT (GAUXC_USE_GPU) and
 *  sn-LinK (SNLINK_USE_GPU) interfaces.
 */

#ifdef USING_gauxc

#include <exception>
#include <string>
#include <utility>

#include "psi4/libpsi4util/exception.h"

namespace psi {
namespace gauxc_device {

/*! For a GPU request made through `option` (e.g., "GAUXC_USE_GPU"), throws a readable PsiException
 *  if GauXC was built without device support or no usable CUDA device is present. Otherwise returns
 *  a description of the device for output (empty when Psi4 can't query CUDA itself).
 */
std::string check(const std::string& option);

/*! Explains a GauXC failure on the GPU: names the device and the likely causes (a GauXC build
 *  lacking this GPU's architecture), and how to fall back to the CPU through `option`.
 */
std::string failure_message(const std::string& what, const std::string& option, const std::string& device_desc);

/*! Runs `step` (GauXC device setup or first kernel launches) and returns its result. When `use_gpu`,
 *  a GauXC failure is rethrown as a PsiException explained by failure_message().
 */
template <typename Step>
decltype(auto) run(bool use_gpu, const std::string& option, const std::string& device_desc, Step&& step) {
    if (!use_gpu) return std::forward<Step>(step)();
    try {
        return std::forward<Step>(step)();
    } catch (const PsiException&) {
        throw;
    } catch (const std::exception& e) {
        throw PSIEXCEPTION(failure_message(e.what(), option, device_desc));
    }
}

}  // namespace gauxc_device
}  // namespace psi

#endif  // USING_gauxc

#endif  // LIBFOCK_GAUXC_DEVICE_H
