/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

#include <stdexcept>
#include "psi4/libciomr/libciomr.h"
#include "psi4/physconst.h"
#include "psi4/libmints/giao_potential.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))

;
using namespace psi;

GiaoPotentialInt::GiaoPotentialInt(std::vector<SphericalTransform>& st, std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2, int deriv) :
    OneBodyAOInt(st, bs1, bs2, deriv), overlap_recur_(bs1->max_am()+2, bs2->max_am()+2)
{
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    int maxnao1 = (maxam1+1)*(maxam1+2)/2;
    int maxnao2 = (maxam2+1)*(maxam2+2)/2;

    // deriv_ == 0 -> first order GIAO derivatives 
    // Increase buffer size to handle x, y, and z components
    if (deriv_ == 0) {
        buffer_ = new double[3*maxnao1*maxnao2];
        set_chunks(3);
    }

    if (deriv_ > 0) {
        throw std::runtime_error("GiaoPotentialInt: does not support any derivatives");
    }
}

GiaoPotentialInt::~GiaoPotentialInt()
{
    delete[] buffer_;
}

 /*  The engine only supports segmented basis sets */

void GiaoPotentialInt::compute_pair(const GaussianShell& s1, const GaussianShell& s2)
{
}
