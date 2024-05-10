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

#ifndef psi_include_psi4_dec_h
#define psi_include_psi4_dec_h

#include "psi4/pragma.h"
#include <string>
#include <memory>

/// This is all defined in python.cc initialize
namespace psi {

class PsiOutStream;
extern PSI_API std::shared_ptr<PsiOutStream> outfile;
extern std::string outfile_name;

extern char *psi_file_prefix;
extern std::string restart_id;  // Does not have a default

enum PsiReturnType { Success, Failure, Balk, EndLoop };

// Very useful regex for matching floating point numbers
#define NUMBER "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

// Hard zero
#define PSI_ZERO 1.0E-14

void die_if_not_converged();

}  // End namespace psi
#endif
