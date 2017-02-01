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

    void Gijkl(void)
    {
      dpdbuf4 V, G;
      int G_irr;
      G_irr = params.G_irr;

      if(params.ref == 0) { /** RHF **/
	global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 0, 0, 0, 0, 0, "VMnIj");
	global_dpd_->buf4_copy(&V, PSIF_CC_GAMMA, "GIjKl");
	global_dpd_->buf4_close(&V);
	if (params.ground) {
	  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 0, 0, 0, 0, 0, "GIjKl");
	  global_dpd_->buf4_symm(&G);
	  global_dpd_->buf4_close(&G);
	}
      }
      else if(params.ref == 1) { /** ROHF **/

	global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 2, 2, 2, 2, 0, "VMNIJ");
	global_dpd_->buf4_copy(&V, PSIF_CC_GAMMA, "GIJKL");
	global_dpd_->buf4_close(&V);
	if (params.ground) {
	  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "GIJKL");
	  global_dpd_->buf4_symm(&G);
	  global_dpd_->buf4_close(&G);
	}

	global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 2, 2, 2, 2, 0, "Vmnij");
	global_dpd_->buf4_copy(&V, PSIF_CC_GAMMA, "Gijkl");
	global_dpd_->buf4_close(&V);
	if (params.ground) {
	  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "Gijkl");
	  global_dpd_->buf4_symm(&G);
	  global_dpd_->buf4_close(&G);
	}

	global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 0, 0, 0, 0, 0, "VMnIj");
	global_dpd_->buf4_copy(&V, PSIF_CC_GAMMA, "GIjKl");
	global_dpd_->buf4_close(&V);
	if (params.ground) {
	  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 0, 0, 0, 0, 0, "GIjKl");
	  global_dpd_->buf4_symm(&G);
	  global_dpd_->buf4_close(&G);
	}
      }
      else if(params.ref == 2) { /** UHF **/

	global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 2, 2, 2, 2, 0, "VMNIJ");
	global_dpd_->buf4_copy(&V, PSIF_CC_GAMMA, "GIJKL");
	global_dpd_->buf4_close(&V);
	if (params.ground) {
	  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "GIJKL");
	  global_dpd_->buf4_symm(&G);
	  global_dpd_->buf4_close(&G);
	}

	global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 12, 12, 12, 12, 0, "Vmnij");
	global_dpd_->buf4_copy(&V, PSIF_CC_GAMMA, "Gijkl");
	global_dpd_->buf4_close(&V);
	if (params.ground) {
	  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 12, 12, 12, 12, 0, "Gijkl");
	  global_dpd_->buf4_symm(&G);
	  global_dpd_->buf4_close(&G);
	}

	global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 22, 22, 22, 22, 0, "VMnIj");
	global_dpd_->buf4_copy(&V, PSIF_CC_GAMMA, "GIjKl");
	global_dpd_->buf4_close(&V);
	if (params.ground) {
	  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 22, 22, 22, 22, 0, "GIjKl");
	  global_dpd_->buf4_symm(&G);
	  global_dpd_->buf4_close(&G);
	}
      }
    }

  }} // namespace psi::ccdensity
