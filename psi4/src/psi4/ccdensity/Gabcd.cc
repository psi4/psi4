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

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void Gabcd(void)
{
  dpdbuf4 G, L, T;
  int G_irr;
  G_irr = params.G_irr;
  bool T2_L2_V = true;

  /*  T2 * L2 * V is absent in CC2 Lagrangian */
  if (params.wfn == "CC2" && params.dertype ==1) T2_L2_V = false;

  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 5, 5, 5, 5, 0, "GAbCd");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 0, 5, 0, 5, 0, "LIjAb");
    if (T2_L2_V) 
       global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    else 
       global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "t1_IjAb");
    global_dpd_->contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    if (params.ground)
      global_dpd_->buf4_symm(&G);
    global_dpd_->buf4_close(&G);
  }
  else if(params.ref == 1) { /** RHF/ROHF **/

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 2, 7, 2, 7, 0, "LIJAB");
    if (T2_L2_V) 
       global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    else 
       global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "t1_IJAB");
    global_dpd_->contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    if (params.ground)
      global_dpd_->buf4_symm(&G);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "Gabcd");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 2, 7, 2, 7, 0, "Lijab");
    if (T2_L2_V) 
       global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    else 
       global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "t1_ijab");
    global_dpd_->contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    if (params.ground)
      global_dpd_->buf4_symm(&G);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 5, 5, 5, 5, 0, "GAbCd");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 0, 5, 0, 5, 0, "LIjAb");
    if (T2_L2_V) 
       global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    else 
       global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "t1_IjAb");
    global_dpd_->contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    if (params.ground)
      global_dpd_->buf4_symm(&G);
    global_dpd_->buf4_close(&G);
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 2, 7, 2, 7, 0, "LIJAB");
    if (T2_L2_V)
       global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    else 
       global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "t1_IJAB");
    global_dpd_->contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_symm(&G);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 17, 17, 17, 17, 0, "Gabcd");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 12, 17, 12, 17, 0, "Lijab");
    if (T2_L2_V)
       global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    else 
       global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "t1_ijab");
    global_dpd_->contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_symm(&G);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 28, 28, 28, 28, 0, "GAbCd");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 22, 28, 22, 28, 0, "LIjAb");
    if (T2_L2_V) 
       global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    else 
       global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "t1_IjAb");
    global_dpd_->contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_symm(&G);
    global_dpd_->buf4_close(&G);
  }
}

}} // namespace psi::ccdensity
