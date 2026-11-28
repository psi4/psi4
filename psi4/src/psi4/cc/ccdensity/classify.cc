/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libiwl/iwl_writer.h"
#include "psi4/libciomr/libciomr.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#include "globals.h"

namespace psi {
namespace ccdensity {

void classify(int p, int q, int r, int s, double value, IWLWriter &ABuf, IWLWriter &BBuf, IWLWriter &CBuf,
              IWLWriter &DBuf, IWLWriter &EBuf, IWLWriter &FBuf) {
    int *occ, *vir, *socc;
    int *cc_occ, *cc_vir;
    int dirac = 1;
    int soccs;

    /* NB that integrals involving frozen orbitals are KEPT in this code */

    occ = frozen.occ;
    vir = frozen.vir;
    socc = frozen.socc;
    cc_occ = frozen.allcc_occ;
    cc_vir = frozen.allcc_vir;

    soccs = socc[p] + socc[q] + socc[r] + socc[s];

    /* A (oo|oo) integrals */
    if ((occ[p] && occ[q] && occ[r] && occ[s]))
        ABuf.write(cc_occ[p], cc_occ[q], cc_occ[r], cc_occ[s], value, dirac);

    /* B (vv|vv) integrals */
    if ((vir[p] && vir[q] && vir[r] && vir[s]))
        BBuf.write(cc_vir[p], cc_vir[q], cc_vir[r], cc_vir[s], value, dirac);

    /* C (oo|vv) integrals */
    if (soccs > 1) {
        if ((occ[p] && occ[q] && vir[r] && vir[s]))
            CBuf.write(cc_occ[p], cc_occ[q], cc_vir[r], cc_vir[s], value, dirac);
        if ((occ[r] && occ[s] && vir[p] && vir[q]))
            CBuf.write(cc_occ[r], cc_occ[s], cc_vir[p], cc_vir[q], value, dirac);
    } else if ((occ[p] && occ[q] && vir[r] && vir[s]))
        CBuf.write(cc_occ[p], cc_occ[q], cc_vir[r], cc_vir[s], value, dirac);
    else if ((occ[r] && occ[s] && vir[p] && vir[q]))
        CBuf.write(cc_occ[r], cc_occ[s], cc_vir[p], cc_vir[q], value, dirac);

    /* D (ov|ov) integrals */
    if (soccs > 1) {
        if ((occ[p] && vir[q] && occ[r] && vir[s]))
            DBuf.write(cc_occ[p], cc_vir[q], cc_occ[r], cc_vir[s], value, dirac);
        if ((occ[q] && vir[p] && occ[r] && vir[s]))
            DBuf.write(cc_occ[q], cc_vir[p], cc_occ[r], cc_vir[s], value, dirac);
        if ((occ[p] && vir[q] && occ[s] && vir[r]))
            DBuf.write(cc_occ[p], cc_vir[q], cc_occ[s], cc_vir[r], value, dirac);
        if ((occ[q] && vir[p] && occ[s] && vir[r]))
            DBuf.write(cc_occ[q], cc_vir[p], cc_occ[s], cc_vir[r], value, dirac);
    } else if ((occ[p] && vir[q] && occ[r] && vir[s]))
        DBuf.write(cc_occ[p], cc_vir[q], cc_occ[r], cc_vir[s], value, dirac);
    else if ((occ[q] && vir[p] && occ[r] && vir[s]))
        DBuf.write(cc_occ[q], cc_vir[p], cc_occ[r], cc_vir[s], value, dirac);
    else if ((occ[p] && vir[q] && occ[s] && vir[r]))
        DBuf.write(cc_occ[p], cc_vir[q], cc_occ[s], cc_vir[r], value, dirac);
    else if ((occ[q] && vir[p] && occ[s] && vir[r]))
        DBuf.write(cc_occ[q], cc_vir[p], cc_occ[s], cc_vir[r], value, dirac);

    /* E (vo|oo) integrals */
    if (soccs > 1) {
        if ((vir[p] && occ[q] && occ[r] && occ[s]))
            EBuf.write(cc_vir[p], cc_occ[q], cc_occ[r], cc_occ[s], value, dirac);
        if ((vir[q] && occ[p] && occ[r] && occ[s]))
            EBuf.write(cc_vir[q], cc_occ[p], cc_occ[r], cc_occ[s], value, dirac);
        if ((vir[r] && occ[s] && occ[p] && occ[q]))
            EBuf.write(cc_vir[r], cc_occ[s], cc_occ[p], cc_occ[q], value, dirac);
        if ((vir[s] && occ[r] && occ[p] && occ[q]))
            EBuf.write(cc_vir[s], cc_occ[r], cc_occ[p], cc_occ[q], value, dirac);
    } else if ((vir[p] && occ[q] && occ[r] && occ[s]))
        EBuf.write(cc_vir[p], cc_occ[q], cc_occ[r], cc_occ[s], value, dirac);
    else if ((vir[p] && occ[q] && occ[s] && occ[r]))
        EBuf.write(cc_vir[q], cc_occ[p], cc_occ[r], cc_occ[s], value, dirac);
    else if ((vir[r] && occ[s] && occ[p] && occ[q]))
        EBuf.write(cc_vir[r], cc_occ[s], cc_occ[p], cc_occ[q], value, dirac);
    else if ((vir[r] && occ[s] && occ[q] && occ[p]))
        EBuf.write(cc_vir[s], cc_occ[r], cc_occ[p], cc_occ[q], value, dirac);

    /* F (ov|vv) integrals */
    if (soccs > 1) {
        if ((occ[p] && vir[q] && vir[r] && vir[s]))
            FBuf.write(cc_occ[p], cc_vir[q], cc_vir[r], cc_vir[s], value, dirac);
        if ((occ[q] && vir[p] && vir[r] && vir[s]))
            FBuf.write(cc_occ[q], cc_vir[p], cc_vir[r], cc_vir[s], value, dirac);
        if ((occ[r] && vir[s] && vir[p] && vir[q]))
            FBuf.write(cc_occ[r], cc_vir[s], cc_vir[p], cc_vir[q], value, dirac);
        if ((occ[s] && vir[r] && vir[p] && vir[q]))
            FBuf.write(cc_occ[s], cc_vir[r], cc_vir[p], cc_vir[q], value, dirac);
    } else if ((occ[p] && vir[q] && vir[r] && vir[s]))
        FBuf.write(cc_occ[p], cc_vir[q], cc_vir[r], cc_vir[s], value, dirac);
    else if ((occ[q] && vir[p] && vir[r] && vir[s]))
        FBuf.write(cc_occ[q], cc_vir[p], cc_vir[r], cc_vir[s], value, dirac);
    else if ((occ[r] && vir[s] && vir[p] && vir[q]))
        FBuf.write(cc_occ[r], cc_vir[s], cc_vir[p], cc_vir[q], value, dirac);
    else if ((occ[s] && vir[r] && vir[p] && vir[q]))
        FBuf.write(cc_occ[s], cc_vir[r], cc_vir[p], cc_vir[q], value, dirac);
}

}  // namespace ccdensity
}  // namespace psi
