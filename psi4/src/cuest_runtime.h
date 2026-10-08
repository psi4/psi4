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
// The interface to cuEST was contributed by NVIDIA under the following terms:
// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef _psi_src_cuest_runtime_h_
#define _psi_src_cuest_runtime_h_

/*! \file cuest_runtime.h
 *  \brief Psi4-module-lifetime interface to the cuEST/CUDA runtime.
 *
 *  This header is safe to include unconditionally: in a build without cuEST it
 *  still declares the whole interface, with the entry points inlined away to
 *  nothing.  Callers therefore never need an `#ifdef USING_cuEST` of their own.
 *
 *  The implementation lives in cuest_runtime.cc, which is only added to the
 *  `core` target when cuEST is enabled.  Everything to do with GPU discovery,
 *  compute-capability and device-health checking, and CUDA/cuBLAS/cuSOLVER/
 *  cuEST handle lifetime is there rather than in core.cc.
 */

namespace psi {
namespace cuest_runtime {

/*! Whether Psi4 was *built* with cuEST support.
 *
 *  Compile-time only.  A build with cuEST says nothing about whether the
 *  machine actually has a usable GPU -- that is not known until the first
 *  cuEST-requiring computation triggers initialization, which is why this is
 *  the only query offered here.
 */
constexpr bool built_with_cuest() {
#ifdef USING_cuEST
    return true;
#else
    return false;
#endif
}

#ifdef USING_cuEST

/*! Release the cuEST/CUDA handles held for the lifetime of the psi4 module.
 *
 *  Idempotent and safe on a machine that never brought cuEST up (no GPU, no
 *  cuEST-requiring computation, or a failed initialization): initialization is
 *  lazy, so in all those cases there is nothing to release and this returns
 *  immediately.  Called from psi4_python_module_finalize().
 */
void shutdown();

#else

inline void shutdown() {}

#endif  // USING_cuEST

}  // namespace cuest_runtime
}  // namespace psi

#endif  // _psi_src_cuest_runtime_h_
