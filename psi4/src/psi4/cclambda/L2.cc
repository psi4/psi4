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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

/* L2_build():

*/

void DL2(struct L_Params L_params);
void FaeL2(int L_irr);
void FmiL2(int L_irr);
void WijmnL2(int L_irr);
void WefabL2(int L_irr);
void WejabL2(int L_irr);
void WijmbL2(int L_irr);
void L1FL2(int L_irr);
void WmbejL2(int L_irr);
void GaeL2(int L_irr);
void GmiL2(int L_irr);
void dijabL2(int L_irr);

void BL2_AO(int L_irr);
void status(const char *, std::string);

void L2_build(struct L_Params L_params) {
  dpdbuf4 L2;
  int L_irr;
  L_irr = L_params.irrep;

  DL2(L_params);
  if(params.print & 2) status("<ij||ab> -> L2", "outfile");

#ifdef EOM_DEBUG
  check_sum("DL2", L_irr);
#endif

  WijmnL2(L_irr);
#ifdef EOM_DEBUG
  check_sum("WijmnL2", L_irr);
#endif
  if(params.print & 2) status("Wmnij -> L2", "outfile");

  WefabL2(L_irr);
#ifdef EOM_DEBUG
  check_sum("WefabL2", L_irr);
#endif
  if(params.print & 2) status("Wabef -> L2", "outfile");

  WejabL2(L_irr);
#ifdef EOM_DEBUG
  check_sum("WejabL2", L_irr);
#endif
  if(params.print & 2) status("Wamef -> L2", "outfile");

  WijmbL2(L_irr);
#ifdef EOM_DEBUG
  check_sum("WijmbL2", L_irr);
#endif
  if(params.print & 2) status("Wmnie -> L2", "outfile");

  GaeL2(L_irr);
#ifdef EOM_DEBUG
  check_sum("GaeL2", L_irr);
#endif

  GmiL2(L_irr);
#ifdef EOM_DEBUG
  check_sum("GmiL2", L_irr);
#endif
  if(params.print & 2) status("G -> L2", "outfile");

  /* For RHF-CCSD response calculations, save all the above
     contributions to the L2 residual for use in the ccresponse code
     (specifically, HX1Y1 and LHX1Y1). */
  if(params.ref == 0 && params.dertype == 3) {
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_copy(&L2, PSIF_CC_LAMPS, "LHX1Y1 Residual I");
    global_dpd_->buf4_close(&L2);
  }

  FaeL2(L_irr);
#ifdef EOM_DEBUG
  check_sum("FaeL2", L_irr);
#endif

  FmiL2(L_irr);
#ifdef EOM_DEBUG
  check_sum("FmiL2", L_irr);
#endif
  if(params.print & 2) status("F -> L2", "outfile");

  WmbejL2(L_irr);
#ifdef EOM_DEBUG
  check_sum("WmbejL2", L_irr);
#endif
  if(params.print & 2) status("Wmbej -> L2", "outfile");

  if(!params.sekino) L1FL2(L_irr); /* should be dropped for Sekino-Bartlett modelIII approach */
#ifdef EOM_DEBUG
  check_sum("L1FL2", L_irr);
#endif
  if(params.print & 2) status("L1*F -> L2", "outfile");

  dijabL2(L_irr);
#ifdef EOM_DEBUG
  check_sum("after D2s", L_irr);
#endif
  if(params.print & 2) status("L2 amplitudes", "outfile");
}


}} // namespace psi::cclambda
