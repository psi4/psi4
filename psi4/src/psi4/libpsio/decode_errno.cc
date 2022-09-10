/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include <cerrno>
#include <cstring>
#include "psi4/libpsio/psio.h"

namespace psi {

/*!
 ** \ingroup PSIO
 **
 ** decode_errno(): Decode the current value of the errno variable to obtain a human-readable error message from the
 ** operating system. This can be used to figure out why a preceeding system call (eg. lseek) has failed. The error
 ** message is returned as a convenient std::string.
 */
std::string decode_errno() { return std::string(std::strerror(errno)); }
}  // namespace psi
