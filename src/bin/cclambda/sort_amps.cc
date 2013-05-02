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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void sort_amps(int L_irr)
{
  dpdbuf4 L2;

  if(params.ref == 0) {
    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_scmcopy(&L2, PSIF_CC_LAMBDA, "2 LIjAb - LIjBa", 2);
    dpd_buf4_sort_axpy(&L2, PSIF_CC_LAMBDA, pqsr, 0, 5, "2 LIjAb - LIjBa", -1);
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 10, 10, "LIAjb");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, psqr, 10, 10, "LIbjA");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_buf4_scmcopy(&L2, PSIF_CC_LAMBDA, "2 LIAjb - LIbjA", 2);
    dpd_buf4_sort_axpy(&L2, PSIF_CC_LAMBDA, psrq, 10, 10, "2 LIAjb - LIbjA", -1);
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, qpsr, 0, 5, "LiJaB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 10, 10, "LiaJB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, rqps, 10, 10, "LjAIb");
    dpd_buf4_close(&L2);
  }
  
  if(params.ref == 1) { /** ROHF **/
    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 10, 10, "LIAjb");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, psqr, 10, 10, "LIbjA");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, qpsr, 0, 5, "LiJaB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 10, 10, "LiaJB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, rqps, 10, 10, "LjAIb");
    dpd_buf4_close(&L2);
    
    /* Build L2IAJB List */
    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 10, 10, "LIAJB");
    dpd_buf4_close(&L2);
    /* Build L2iajb List */
    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "Lijab");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 10, 10, "Liajb");
    dpd_buf4_close(&L2);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, qpsr, 23, 29, "LiJaB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 20, 20, "LIAJB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 10, 15, 12, 17, 0, "Lijab");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 30, 30, "Liajb");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 20, 30, "LIAjb");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, prqs, 30, 20, "LiaJB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, psrq, 24, 27, "LIbjA");
    dpd_buf4_sort(&L2, PSIF_CC_LAMBDA, rqps, 27, 24, "LjAIb");
    dpd_buf4_close(&L2);
  }

}


}} // namespace psi::cclambda
