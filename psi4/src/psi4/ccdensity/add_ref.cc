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
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void add_ref(struct iwlbuf *OutBuf)
{
  int i,j;
  int nfzc, nclsd, nopen;

  nfzc = moinfo.nfzc;
  nclsd = moinfo.nclsd;
  nopen = moinfo.nopen;

  /*** One-electron component ***/

  for(i=0; i < (nfzc + nclsd); i++)
      moinfo.opdm[i][i] += 2.0;

  for(i=nfzc + nclsd; i < (nfzc + nclsd + nopen); i++)
      moinfo.opdm[i][i] += 1.0;


  /*** Two-electron component ***/

  /* docc-docc */
  for(i=0; i < (nfzc + nclsd); i++) {
      iwl_buf_wrt_val(OutBuf, i, i, i, i, 1.0, 0, "outfile", 0);
      for(j=0; j < i; j++) {
	  iwl_buf_wrt_val(OutBuf, i, i, j, j, 2.0, 0, "outfile", 0);
	  iwl_buf_wrt_val(OutBuf, i, j, j, i,-1.0, 0, "outfile", 0);
	}
    }

  /* socc-docc && socc-socc*/
  for(i=(nfzc + nclsd); i < (nfzc + nclsd + nopen); i++) {
      for(j=0; j < (nfzc + nclsd); j++) {
	  iwl_buf_wrt_val(OutBuf, i, i, j, j, 1.0, 0, "outfile", 0);
	  iwl_buf_wrt_val(OutBuf, i, j, j, i,-0.5, 0, "outfile", 0);
	}
      for(j=(nfzc + nclsd); j < i; j++) {
	  iwl_buf_wrt_val(OutBuf, i, i, j, j, 0.5, 0, "outfile", 0);
	  iwl_buf_wrt_val(OutBuf, i, j, j, i,-0.5, 0, "outfile", 0);
	}
    }
}

}} // namespace psi::ccdensity
