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
#include <cstring>
#include "psi4/libpsio/psio.h"

namespace psi {
/// @brief
/// @param context
/// @param unit
/// @param errno_in
/// @return
std::string psio_read_err_msg(const std::string& context, const size_t unit, const int errno_in) {
    std::string errmsg = "READ failed. Error description from the OS: " + decode_errno(errno_in);
    errmsg += '\n' + context + ", unit ";
    errmsg += std::to_string(unit) + ".\n";
}

/// @brief
/// @param context
/// @param unit
/// @return
std::string psio_read_err_msg_some(const std::string& context, const size_t unit) {
    std::string errmsg = "READ failed. Only some of the bytes were read!";
    errmsg += '\n' + context + ", unit ";
    errmsg += std::to_string(unit) + ".\n";
}
}  // namespace psi
