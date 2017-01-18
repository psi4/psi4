/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef THREE_INDEX_FITTER
#define THREE_INDEX_FITTER

#include "psi4/psi4-dec.h"


namespace psi {

class BasisSet;
class Matrix;
class Vector;

/*!
 * Small utility class to compute d_A = J_AB^{-1} (Q|mn) D_mn coefficients for QM/MM
 */
class DFChargeFitter {

protected:
    /// Print flag (defaults to 1)
    int print_;
    /// Debug flag (defaults to 0)
    int debug_;

    /// Target coefficients
    SharedVector d_;
    /// Driving density
    SharedMatrix D_;

    /// Primary Basis Set
    std::shared_ptr<BasisSet> primary_;
    /// Auxiliary Basis Set
    std::shared_ptr<BasisSet> auxiliary_;

public:

    DFChargeFitter();
    ~DFChargeFitter();

    SharedVector fit();

    void setD(SharedMatrix D) { D_ = D; }
    void setPrimary(std::shared_ptr<BasisSet> primary) { primary_ = primary; }
    void setAuxiliary(std::shared_ptr<BasisSet> auxiliary) { auxiliary_ = auxiliary; }

    SharedVector d() const { return d_; }

    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }
};

} // Namespace psi
#endif
