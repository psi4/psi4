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

/*!
 ** \file
 ** \ingroup PSIO
 */
#include <cstring>
#include "psi4/libpsio/psio.h"

namespace psi {

/// @brief Decode the current value of the errno variable to obtain a human-readable error message from the operating
/// system. This can be used to figure out why a preceeding system call (eg. lseek) has failed. Callers should save the
/// value of errno into a local variable **immediately** after the call that should be checked returns, because there
/// are a lot of things in C++ that can potentially fail and overwrite the global errno with a new error code.
/// @param errno_in : the error code from the OS
/// @return Human-readable error message from the OS. May or may not be in English.
std::string decode_errno(const int errno_in) { return std::string(std::strerror(errno_in)); }
}  // namespace psi
