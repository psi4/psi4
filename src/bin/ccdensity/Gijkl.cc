/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
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
	dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 0, 0, 0, 0, 0, "VMnIj");
	dpd_->buf4_copy(&V, PSIF_CC_GAMMA, "GIjKl");
	dpd_->buf4_close(&V);
	if (params.ground) {
	  dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 0, 0, 0, 0, 0, "GIjKl");
	  dpd_->buf4_symm(&G);
	  dpd_->buf4_close(&G);
	}
      }
      else if(params.ref == 1) { /** ROHF **/

	dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 2, 2, 2, 2, 0, "VMNIJ");
	dpd_->buf4_copy(&V, PSIF_CC_GAMMA, "GIJKL");
	dpd_->buf4_close(&V);
	if (params.ground) {
	  dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "GIJKL");
	  dpd_->buf4_symm(&G);
	  dpd_->buf4_close(&G);
	}

	dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 2, 2, 2, 2, 0, "Vmnij");
	dpd_->buf4_copy(&V, PSIF_CC_GAMMA, "Gijkl");
	dpd_->buf4_close(&V);
	if (params.ground) {
	  dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "Gijkl");
	  dpd_->buf4_symm(&G);
	  dpd_->buf4_close(&G);
	}

	dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 0, 0, 0, 0, 0, "VMnIj");
	dpd_->buf4_copy(&V, PSIF_CC_GAMMA, "GIjKl");
	dpd_->buf4_close(&V);
	if (params.ground) {
	  dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 0, 0, 0, 0, 0, "GIjKl");
	  dpd_->buf4_symm(&G);
	  dpd_->buf4_close(&G);
	}
      }
      else if(params.ref == 2) { /** UHF **/

	dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 2, 2, 2, 2, 0, "VMNIJ");
	dpd_->buf4_copy(&V, PSIF_CC_GAMMA, "GIJKL");
	dpd_->buf4_close(&V);
	if (params.ground) {
	  dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 2, 2, 2, 2, 0, "GIJKL");
	  dpd_->buf4_symm(&G);
	  dpd_->buf4_close(&G);
	}

	dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 12, 12, 12, 12, 0, "Vmnij");
	dpd_->buf4_copy(&V, PSIF_CC_GAMMA, "Gijkl");
	dpd_->buf4_close(&V);
	if (params.ground) {
	  dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 12, 12, 12, 12, 0, "Gijkl");
	  dpd_->buf4_symm(&G);
	  dpd_->buf4_close(&G);
	}

	dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 22, 22, 22, 22, 0, "VMnIj");
	dpd_->buf4_copy(&V, PSIF_CC_GAMMA, "GIjKl");
	dpd_->buf4_close(&V);
	if (params.ground) {
	  dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 22, 22, 22, 22, 0, "GIjKl");
	  dpd_->buf4_symm(&G);
	  dpd_->buf4_close(&G);
	}
      }
    }

  }} // namespace psi::ccdensity
