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

void G_norm(void) {
  dpdfile2 G1;
  dpdbuf4 G;
  double value, value1, dot_IA, dot_ia, dot_AI, dot_ai;
  int G_irr = 0;

  outfile->Printf("Calculating overlaps of CC_OEI\n");
  global_dpd_->file2_init(&G1, PSIF_CC_OEI, G_irr, 0, 1, "DIA");
  dot_IA = global_dpd_->file2_dot_self(&G1);
  global_dpd_->file2_close(&G1);
  global_dpd_->file2_init(&G1, PSIF_CC_OEI, G_irr, 0, 1, "Dia");
  dot_ia = global_dpd_->file2_dot_self(&G1);
  global_dpd_->file2_close(&G1);
  global_dpd_->file2_init(&G1, PSIF_CC_OEI, G_irr, 0, 1, "DAI");
  dot_AI = global_dpd_->file2_dot_self(&G1);
  global_dpd_->file2_close(&G1);
  global_dpd_->file2_init(&G1, PSIF_CC_OEI, G_irr, 0, 1, "Dai");
  dot_ai = global_dpd_->file2_dot_self(&G1);
  global_dpd_->file2_close(&G1);
  /*
  outfile->Printf("<DIA|DIA> = %15.10lf\n", dot_IA);
  outfile->Printf("<Dia|Dia> = %15.10lf\n", dot_ia);
  outfile->Printf("<DAI|DAI> = %15.10lf\n", dot_AI);
  outfile->Printf("<Dai|Dai> = %15.10lf\n", dot_ai);
  */
  outfile->Printf("\t<Dpq|Dqp>     = %15.10lf\n", dot_IA+dot_ia+dot_AI+dot_ai);

  outfile->Printf("Calculating overlaps of CC_GAMMA\n");

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "GIJKL");
  value = global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "Gijkl");
  value += global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 0, 0, 0, 0, 0, "GIjKl");
  value += global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  outfile->Printf("\t<Gijkl|Gijkl> = %15.10lf\n", value);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 10, 2, 10, 0, "GIJKA");
  value = global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 10, 2, 10, 0, "Gijka");
  value += global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  value += global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  value += global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  outfile->Printf("\t<Gijka|Gijka> = %15.10lf\n",value);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "GIJAB");
  value = global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "Gijab");
  value += global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  value += global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  outfile->Printf("\t<Gijab|Gijab> = %15.10lf\n", value);


  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIBJA");
  value = global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "Gibja");
  value += global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbJa");
  value += global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBjA");
  value += global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbjA");
  value += global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBJa");
  value += global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  outfile->Printf("\t<Gibja|Gibja> = %15.10lf\n",value);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  value = global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 11, 7, 11, 7, 0, "Gciab");
  value += global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  value += global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  value += global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  outfile->Printf("\t<Gciab|Gciab> = %15.10lf\n",value);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
  value = global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "Gabcd");
  value += global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 5, 5, 5, 5, 0, "GAbCd");
  value += global_dpd_->buf4_dot_self(&G);
  global_dpd_->buf4_close(&G);
  outfile->Printf("\t<Gabcd|Gabcd> = %15.10lf\n", value);

  return;
}

}} // namespace psi::ccdensity
