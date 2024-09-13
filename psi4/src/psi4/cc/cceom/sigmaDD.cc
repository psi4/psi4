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

/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cceom {

void FDD(int i, int C_irr);
void WabefDD(int i, int C_irr);
void WmnijDD(int i, int C_irr);
void WmbejDD(int i, int C_irr);
void WmnefDD(int i, int C_irr);

/* This function computes the H-bar doubles-doubles block contribution
to a Sigma vector stored at Sigma plus 'i' */

void sigmaDD(int i, int C_irr) {
    timer_on("FDD");
    FDD(i, C_irr);
    timer_off("FDD");
    timer_on("WmnijDD");
    WmnijDD(i, C_irr);
    timer_off("WmnijDD");
    timer_on("WabefDD");
    WabefDD(i, C_irr);
    timer_off("WabefDD");
    timer_on("WmbejDD");
    WmbejDD(i, C_irr);
    timer_off("WmbejDD");
    timer_on("WmnefDD");
    WmnefDD(i, C_irr);
    timer_off("WmnefDD");

    return;
}

}  // namespace cceom
}  // namespace psi
