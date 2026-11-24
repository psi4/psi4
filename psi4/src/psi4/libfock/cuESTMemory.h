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

#ifndef CUEST_MEMORY_H
#define CUEST_MEMORY_H

#ifdef USING_cuEST

#include <cstddef>

#include <cuda_runtime.h>

#include "psi4/libfock/cuESTCommon.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/psi4-dec.h"

namespace psi {
namespace cuest_common {

/**
 * Device-memory report table for cuEST workflows (all sizes in GiB).
 *
 * print_memory_header(): prints the device total and the column labels.
 * print_memory_row():    prints one step. Call it BEFORE allocating that step's buffers.
 *                        "In use" is the projected device usage while the step runs:
 *                        memory currently in use (cudaMemGetInfo) + persistent + temporary.
 */

constexpr double bytes_per_gib = 1024.0 * 1024.0 * 1024.0;

inline void print_memory_header() {
    size_t free_bytes = 0, total_bytes = 0;
    CHECK_CUDA(cudaMemGetInfo(&free_bytes, &total_bytes));

    outfile->Printf("  ==> cuEST Device Memory [GiB] <==\n\n");
    outfile->Printf("    Device total: %.3f\n\n", total_bytes / bytes_per_gib);
    outfile->Printf("    %-20s %12s %12s %12s\n", "Step", "Persistent", "Temporary", "In use");
    outfile->Printf("    %-20s %12s %12s %12s\n", "--------------------", "----------", "----------", "----------");
}

inline void print_memory_row(const char* label, size_t persistent_bytes, size_t temporary_bytes) {
    size_t free_bytes = 0, total_bytes = 0;
    CHECK_CUDA(cudaMemGetInfo(&free_bytes, &total_bytes));
    const size_t in_use_bytes = (total_bytes - free_bytes) + persistent_bytes + temporary_bytes;

    outfile->Printf("    %-20s %12.3f %12.3f %12.3f\n", label, persistent_bytes / bytes_per_gib,
                    temporary_bytes / bytes_per_gib, in_use_bytes / bytes_per_gib);
}

}  // namespace cuest_common
}  // namespace psi

#endif  // USING_cuEST
#endif  // CUEST_MEMORY_H
