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

void zero_onepdm(struct RHO_Params rho_params)
{
  dpdfile2 D;
  int G_irr;
  G_irr = params.G_irr;

  if ( (params.ref == 0) || (params.ref == 1) ) {
    global_dpd_->file2_init(&D, PSIF_CC_OEI, G_irr, 0, 0, rho_params.DIJ_lbl);
    global_dpd_->file2_scm(&D, 0.0);
    global_dpd_->file2_close(&D);
    global_dpd_->file2_init(&D, PSIF_CC_OEI, G_irr, 0, 0, rho_params.Dij_lbl);
    global_dpd_->file2_scm(&D, 0.0);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, G_irr, 1, 1, rho_params.DAB_lbl);
    global_dpd_->file2_scm(&D, 0.0);
    global_dpd_->file2_close(&D);
    global_dpd_->file2_init(&D, PSIF_CC_OEI, G_irr, 1, 1, rho_params.Dab_lbl);
    global_dpd_->file2_scm(&D, 0.0);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, G_irr, 0, 1, rho_params.DIA_lbl);
    global_dpd_->file2_scm(&D, 0.0);
    global_dpd_->file2_close(&D);
    global_dpd_->file2_init(&D, PSIF_CC_OEI, G_irr, 0, 1, rho_params.Dia_lbl);
    global_dpd_->file2_scm(&D, 0.0);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, G_irr, 0, 1, rho_params.DAI_lbl);
    global_dpd_->file2_scm(&D, 0.0);
    global_dpd_->file2_close(&D);
    global_dpd_->file2_init(&D, PSIF_CC_OEI, G_irr, 0, 1, rho_params.Dai_lbl);
    global_dpd_->file2_scm(&D, 0.0);
    global_dpd_->file2_close(&D);
  }
  else if (params.ref == 2) {
    global_dpd_->file2_init(&D, PSIF_CC_OEI, G_irr, 0, 0, rho_params.DIJ_lbl);
    global_dpd_->file2_scm(&D, 0.0);
    global_dpd_->file2_close(&D);
    global_dpd_->file2_init(&D, PSIF_CC_OEI, G_irr, 2, 2, rho_params.Dij_lbl);
    global_dpd_->file2_scm(&D, 0.0);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, G_irr, 1, 1, rho_params.DAB_lbl);
    global_dpd_->file2_scm(&D, 0.0);
    global_dpd_->file2_close(&D);
    global_dpd_->file2_init(&D, PSIF_CC_OEI, G_irr, 3, 3, rho_params.Dab_lbl);
    global_dpd_->file2_scm(&D, 0.0);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, G_irr, 0, 1, rho_params.DIA_lbl);
    global_dpd_->file2_scm(&D, 0.0);
    global_dpd_->file2_close(&D);
    global_dpd_->file2_init(&D, PSIF_CC_OEI, G_irr, 2, 3, rho_params.Dia_lbl);
    global_dpd_->file2_scm(&D, 0.0);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, G_irr, 0, 1, rho_params.DAI_lbl);
    global_dpd_->file2_scm(&D, 0.0);
    global_dpd_->file2_close(&D);
    global_dpd_->file2_init(&D, PSIF_CC_OEI, G_irr, 2, 3, rho_params.Dai_lbl);
    global_dpd_->file2_scm(&D, 0.0);
    global_dpd_->file2_close(&D);
  }
}

void zero_twopdm(void)
{
  dpdbuf4 G;
  int G_irr;
  G_irr = params.G_irr;

  if ( (params.ref == 0) || (params.ref == 1) ) {
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "GIJKL");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "Gijkl");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 0, 0, 0, 0, 0, "GIjKl");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "Gabcd");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 5, 5, 5, 5, 0, "GAbCd");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 10, 2, 10, 0, "GIJKA");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 10, 2, 10, 0, "Gijka");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 0, 10, 0, 10, 0, "GIjKa");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 0, 10, 0, 10, 0, "GiJkA");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIBJA");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "Gibja");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbJa");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBjA");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbjA");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBJa");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 11, 7, 11, 7, 0, "GCIAB");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 11, 7, 11, 7, 0, "Gciab");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 11, 5, 11, 5, 0, "GCiAb");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 11, 5, 11, 5, 0, "GcIaB");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "GIJAB");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "Gijab");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 0, 5, 0, 5, 0, "GIjAb");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
  }
  else if (params.ref == 2) {
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "GIJKL");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 12, 12, 12, 12, 0, "Gijkl");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 22, 22, 22, 22, 0, "GIjKl");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 17, 17, 17, 17, 0, "Gabcd");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 28, 28, 28, 28, 0, "GAbCd");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 20, 2, 20, 0, "GIJKA");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 12, 30, 12, 30, 0, "Gijka");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 22, 24, 22, 24, 0, "GIjKa");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 23, 27, 23, 27, 0, "GiJkA");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 20, 20, 20, 20, 0, "GIBJA");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 30, 30, 30, 30, 0, "Gibja");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 24, 24, 24, 24, 0, "GIbJa");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 27, 27, 27, 27, 0, "GiBjA");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 24, 27, 24, 27, 0, "GIbjA");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 27, 24, 27, 24, 0, "GiBJa");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 21, 7, 21, 7, 0, "GCIAB");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 31, 17, 31, 17, 0, "Gciab");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 26, 28, 26, 28, 0, "GCiAb");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 25, 29, 25, 29, 0, "GcIaB");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "GIJAB");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 12, 17, 12, 17, 0, "Gijab");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 24, 28, 24, 28, 0, "GIjAb");
    global_dpd_->buf4_scm(&G, 0.0);
    global_dpd_->buf4_close(&G);
  }
}

}} // namespace psi::ccdensity
