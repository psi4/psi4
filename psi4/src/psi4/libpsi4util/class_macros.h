/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#ifndef _psi_src_lib_libpsi4util_class_macros_h_
#define _psi_src_lib_libpsi4util_class_macros_h_

// This macro may be used in a class to declare
// a private variable with a set and get function
#define READ_WRITE_VAR(type, name)                  \
   public:                                          \
    type get_name() const { return name; }          \
    void set_name(type _value_) { name = _value_; } \
                                                    \
   private:                                         \
    type name;

// This macro may be used in a class to declare
// a private variable with a get function
#define READ_VAR(type, name)               \
   public:                                 \
    type get_name() const { return name; } \
                                           \
   private:                                \
    type name;

#endif  // _psi_src_lib_libpsi4util_class_macros_h_
