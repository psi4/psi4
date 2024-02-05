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

#pragma once

#include <functional>
#include <map>
#include <memory>
#include <string>

#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
class Wavefunction;
class Matrix;
struct dpdbuf4;
namespace fcidump {
/*!  \fn void fcidump_tei_helper(int nirrep, bool restricted, std::map<std::string, int> DPD_info, double
 * ints_tolerance,
 *                  std::string fname = "INTDUMP")
 *  \brief Write integrals to file in FCIDUMP format
 *  \param[in] nirrep number of irreps
 *  \param[in] bool whether RHF or UHF
 *  \param[in] DPD_info DPD instance and MO spaces IDs
 *  \param[in] ints_tolerance tolerance for integrals to be written to file
 *  \param[in] fname name of the FCIDUMP file
 */
void fcidump_tei_helper(int nirrep, bool restricted, std::map<std::string, int> DPD_info, double ints_tolerance,
                        std::string fname = "INTDUMP");

namespace detail {
using OrbitalIndexing = std::function<int(const int)>;

void write_tei_to_disk(std::shared_ptr<PsiOutStream> intdump, int nirrep, dpdbuf4& K, double ints_tolerance,
                       OrbitalIndexing indx1, OrbitalIndexing indx2);
}  // End namespace detail
}  // End namespace fcidump
}  // End namespace psi
