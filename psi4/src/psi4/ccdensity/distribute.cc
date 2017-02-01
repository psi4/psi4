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
#include "psi4/libiwl/iwl.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void classify(int p, int q, int r, int s, double value,
	      struct iwlbuf *ABuf, struct iwlbuf *BBuf,
	      struct iwlbuf *CBuf, struct iwlbuf *DBuf,
	      struct iwlbuf *EBuf, struct iwlbuf *FBuf);

void distribute(void)
{
  double tolerance;
  struct iwlbuf InBuf;
  struct iwlbuf ABuf, BBuf, CBuf, DBuf, EBuf, FBuf;
  int lastbuf;
  Value *valptr;
  Label *lblptr;
  int idx, p, q, r, s;
  double value;

  tolerance = params.tolerance;

  iwl_buf_init(&InBuf, PSIF_MO_TEI, tolerance, 1, 1);
  iwl_buf_init(&ABuf, 90, tolerance, 0, 0);
  iwl_buf_init(&BBuf, 91, tolerance, 0, 0);
  iwl_buf_init(&CBuf, 92, tolerance, 0, 0);
  iwl_buf_init(&DBuf, 93, tolerance, 0, 0);
  iwl_buf_init(&EBuf, 94, tolerance, 0, 0);
  iwl_buf_init(&FBuf, 95, tolerance, 0, 0);

  /* Run through the buffer that's already available */
  lblptr = InBuf.labels;
  valptr = InBuf.values;
  lastbuf = InBuf.lastbuf;

  for (idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
    p = (int) lblptr[idx++];
    q = (int) lblptr[idx++];
    r = (int) lblptr[idx++];
    s = (int) lblptr[idx++];

    value = (double) valptr[InBuf.idx];

    /* Check integral into each class */
    classify(p,q,r,s,value,&ABuf,&BBuf,&CBuf,&DBuf,&EBuf,&FBuf);

/*    outfile->Printf( "(%d %d|%d %d) = %20.10lf\n", p, q, r, s, value);  */

    } /* end loop through current buffer */

  /* Now run through the rest of the buffers in the file */
  while (!lastbuf) {
    iwl_buf_fetch(&InBuf);
    lastbuf = InBuf.lastbuf;

    for (idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
      p = (int) lblptr[idx++];
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];

      value = (double) valptr[InBuf.idx];

      /* Check integral into each class */
      classify(p,q,r,s,value,&ABuf,&BBuf,&CBuf,&DBuf,&EBuf,&FBuf);

/*      outfile->Printf( "(%d %d|%d %d) = %20.10lf\n", p, q, r, s, value); */

      } /* end loop through current buffer */
    } /* end loop over reading buffers */


  iwl_buf_close(&InBuf, 1);
  iwl_buf_flush(&ABuf, 1);
  iwl_buf_flush(&BBuf, 1);
  iwl_buf_flush(&CBuf, 1);
  iwl_buf_flush(&DBuf, 1);
  iwl_buf_flush(&EBuf, 1);
  iwl_buf_flush(&FBuf, 1);
  iwl_buf_close(&ABuf, 1);
  iwl_buf_close(&BBuf, 1);
  iwl_buf_close(&CBuf, 1);
  iwl_buf_close(&DBuf, 1);
  iwl_buf_close(&EBuf, 1);
  iwl_buf_close(&FBuf, 1);


}


}} // namespace psi::ccdensity
