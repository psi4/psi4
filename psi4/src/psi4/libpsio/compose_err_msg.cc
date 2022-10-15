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
#include <optional>
#include "psi4/libpsio/psio.h"

namespace psi {
/// @brief Composes an error message explaining that an IO system call has failed, with the error message from the OS
/// (based on errno), some context from the caller of this function and the unit number. Callers should save the
/// value of errno into a local variable **immediately** after the IO call returns, as there are a lot of things in C++
/// that can potentially fail and overwrite the global errno with a new error code.
/// @param beginning : Beginning of the error message. Should not end with a newline.
/// @param context : Some information about the context of the IO call that has failed. Should not end with a newline.
/// @param unit : Unit number that the IO failed on
/// @param errno_in : The value of errno after the failed IO call
/// @return String explaining the error. Ends with a newline.
std::string psio_compose_err_msg(const std::string& beginning, const std::string& context, const size_t unit,
                                 const std::optional<int> errno_in /* = std::nullopt*/) {
    std::string errmsg = beginning;
    if (errno_in.has_value()) errmsg += " Error description from the OS: " + decode_errno(errno_in.value());
    errmsg += '\n' + context + ", unit " + std::to_string(unit) + ".\n";
    return errmsg;
}
}  // namespace psi
