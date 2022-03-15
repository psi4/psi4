/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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

//
// corrtab.h
//
// Additional modifications made by Justin Turney <jturney@ccqc.uga.edu>
// for use in PSI4.
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef _math_symmetry_corrtab_h
#define _math_symmetry_corrtab_h

#include <iostream>

#include "psi4/libmints/pointgrp.h"

namespace psi {

// //////////////////////////////////////////////////////////////////

/** The CorrelationTable class provides a correlation
    table between two point groups.
*/
class PSI_API CorrelationTable {
   private:
    std::shared_ptr<PointGroup> group_;
    std::shared_ptr<PointGroup> subgroup_;

    int n_;
    int subn_;
    int* ngamma_;
    int** gamma_;

    void clear();

   public:
    CorrelationTable();

    /// Create a correlation table for the two groups.
    CorrelationTable(const std::shared_ptr<PointGroup>& group, const std::shared_ptr<PointGroup>& subgroup);

    ~CorrelationTable();

    /// Returns the higher order point group.
    std::shared_ptr<PointGroup> group() const { return group_; }
    /// Returns the lower order point group.
    std::shared_ptr<PointGroup> subgroup() const { return subgroup_; }

    /** Initalize the correlation table.  Returns 0 for success and nonzero
        for failure.  This will fail if the subgroup is not really a subgroup
        of group. */
    int initialize_table(const std::shared_ptr<PointGroup>& group, const std::shared_ptr<PointGroup>& subgroup);

    /// Converts error codes from initialize_table into a text string.
    const char* error(int errcod);

    /// Returns the number of irreps in the high order group.
    int n() const { return n_; }
    /// Returns the number of irreps in the subgroup.
    int subn() const { return subn_; }
    /// Returns the degeneracy of the irrep.
    int degen(int igamma) const;
    /// Returns the degeneracy of the subgroup irrep.
    int subdegen(int igamma) const;
    /// Returns the number of irreps in the low order group that an irrep
    // from the high order group can be reduced to.
    int ngamma(int igamma) const { return ngamma_[igamma]; }
    /** Returns the irreps in the low order group that an irrep from the
        high order group can be reduced to. */
    int gamma(int igamma, int i) const { return gamma_[igamma][i]; }

    // void print(std::ostream &o=ExEnv::out0()) const;
};

}  // namespace psi

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
