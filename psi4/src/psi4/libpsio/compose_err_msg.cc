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
/// @return 
std::string psio_compose_err_msg_core(const std::string& context, const size_t unit) {
    std::string errmsg =  context + ", unit ";
    errmsg += std::to_string(unit) + ".\n";
    return errmsg;
}

/// @brief Composes an error message explaining that an lseek system call has failed, with the error message from the OS
/// (based on errno), some context from the caller of this function and the unit number. Callers should save the
/// value of errno into a local variable **immediately** after the lseek call that should be checked returns, because
/// there are a lot of things in C++ that can potentially fail and overwrite the global errno with a new error code.
/// @param context : Some information about the context of the lseek call that has failed
/// @param unit : Unit number that that was attempted to be lseek'd
/// @param errno_in : The value of errno after the failed lseek call
/// @return String explaining the error
std::string psio_lseek_err_msg(const std::string& context, const size_t unit, const int errno_in) {
    std::string errmsg = "LSEEK failed. Error description from the OS: " + decode_errno(errno_in) + '\n';
    errmsg += psio_compose_err_msg_core(context, unit);
    return errmsg;
}

/// @brief
/// @param context
/// @param unit
/// @param errno_in
/// @return
std::string psio_read_err_msg(const std::string& context, const size_t unit, const int errno_in) {
    std::string errmsg = "READ failed. Error description from the OS: " + decode_errno(errno_in) + '\n';
    errmsg += psio_compose_err_msg_core(context, unit);
    return errmsg;
}

/// @brief
/// @param context
/// @param unit
/// @return
std::string psio_read_err_msg_some(const std::string& context, const size_t unit) {
    std::string errmsg = "READ failed. Only some of the bytes were read!" + '\n';
    errmsg += psio_compose_err_msg_core(context, unit);
    return errmsg;
}

/// @brief
/// @param context
/// @param unit
/// @param errno_in
/// @return
std::string psio_write_err_msg(const std::string& context, const size_t unit, const int errno_in) {
    std::string errmsg = "WRITE failed. Error description from the OS: " + decode_errno(errno_in) + '\n';
    errmsg += psio_compose_err_msg_core(context, unit);
    return errmsg;
}

/// @brief
/// @param context
/// @param unit
/// @return
std::string psio_write_err_msg_some(const std::string& context, const size_t unit) {
    std::string errmsg = "WRITE failed. Only some of the bytes were written! Maybe the disk is full?" + '\n';
    errmsg += psio_compose_err_msg_core(context, unit);
    return errmsg;
}
}  // namespace psi
