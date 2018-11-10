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

#include "integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include <cmath>
#include <cctype>
#include <cstdio>
#include "psi4/psifiles.h"
#include "mospace.h"

using namespace psi;

/**
 * Transform the two-electron integrals from the SO to the MO basis in the spaces specified
 *
 * @param s1 - the MO space for the first index
 * @param s2 - the MO space for the second index
 * @param s3 - the MO space for the third index
 * @param s4 - the MO space for the fourth index
 */
void IntegralTransform::transform_tei(const std::shared_ptr<MOSpace> s1, const std::shared_ptr<MOSpace> s2,
                                      const std::shared_ptr<MOSpace> s3, const std::shared_ptr<MOSpace> s4,
                                      HalfTrans ht) {
    check_initialized();
    // Only do the first half if the "make" flag is set
    if (ht == HalfTrans::MakeAndKeep || ht == HalfTrans::MakeAndNuke) transform_tei_first_half(s1, s2);

    if (ht == HalfTrans::ReadAndNuke || ht == HalfTrans::MakeAndNuke) {
        keepHtInts_ = false;
    } else {
        keepHtInts_ = true;
    }
    transform_tei_second_half(s1, s2, s3, s4);
}
