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

#ifndef DPD_MOSPACE_H
#define DPD_MOSPACE_H

#include <vector>
#include <string>
#include "psi4/libmints/dimension.h"

namespace psi {

std::vector<std::string> dpd_split(const std::string &indices);

class DPDMOSpace {
   protected:
    // name of the space
    char label_;
    // list of allowed orbital-index labels
    std::vector<std::string> indices_;
    // number of irreps
    int nIrrep_;
    // number of orbitals
    int nOrb_;
    // number of orbitals per irrep
    std::vector<int> orbPI_;
    // irrep of each orbital
    std::vector<int> orbSym_;

   public:
    DPDMOSpace(const char label, const std::string &indices, std::vector<int> orbspi);
    DPDMOSpace(const char label, const std::string &indices, Dimension orbspi);
    DPDMOSpace();
    ~DPDMOSpace();

    char label() { return label_; }
    std::vector<std::string> indices() { return indices_; }
    int nIrrep() { return nIrrep_; }
    int nOrb() { return nOrb_; }
    std::vector<int> orbPI() { return orbPI_; }
    std::vector<int> orbSym() { return orbSym_; }

    void print();
    bool operator==(const char *c);
    bool operator==(const std::string &c);
    bool operator==(DPDMOSpace &lhs);
    friend bool operator==(const char *c, const DPDMOSpace &rhs);
    friend bool operator==(const std::string &c, const DPDMOSpace &rhs);
    std::vector<std::string> overlap(DPDMOSpace &rhs);
};

}  // namespace psi

#endif  // DPD_MOSPACE_H
