/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "cclambda.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cclambda {

void DL2(const struct L_Params& L_params);
void cc2_faeL2(int L_irr);
void cc2_fmiL2(int L_irr);
void WijmnL2(int L_irr);
void WefabL2(int L_irr);
void WejabL2(int L_irr);
void WijmbL2(int L_irr);
void L1FL2(int L_irr);
void dijabL2(int L_irr);

void BL2_AO(int L_irr);

void CCLambdaWavefunction::cc2_L2_build(const struct L_Params& L_params) {
    int L_irr;
    L_irr = L_params.irrep;

    DL2(L_params);
    if (params.print & 2) status("<ij||ab> -> L2", "outfile");

#ifdef EOM_DEBUG
    check_sum("DL2", L_irr);
#endif

    cc2_faeL2(L_irr);

#ifdef EOM_DEBUG
    check_sum("FaeL2", L_irr);
#endif

    cc2_fmiL2(L_irr);

#ifdef EOM_DEBUG
    check_sum("FmiL2", L_irr);
#endif
    if (params.print & 2) status("F -> L2", "outfile");

    WijmbL2(L_irr);

#ifdef EOM_DEBUG
    check_sum("WmnieL2", L_irr);
#endif
    if (params.print & 2) status("Wmnie -> L2", "outfile");

    WejabL2(L_irr);

#ifdef EOM_DEBUG
    check_sum("WejabL2", L_irr);
#endif
    if (params.print & 2) status("Wamef -> L2", "outfile");

    L1FL2(L_irr);

#ifdef EOM_DEBUG
    check_sum("L1FL2", L_irr);
#endif
    if (params.print & 2) status("L1*F -> L2", "outfile");

    dijabL2(L_irr);

#ifdef EOM_DEBUG
    check_sum("after D2s", L_irr);
#endif
    if (params.print & 2) status("L2 amplitudes", "outfile");
}

}  // namespace cclambda
}  // namespace psi
