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

/*!
  \file
  \ingroup IWL
*/
#include <cstdio>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "iwl.h"
#include "iwl.hpp"

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

namespace psi {

void IWL::read_two(PSIO *psio, int itap, double *ints, int *ioff, int norbs,
    int nfzc, int nfzv, int printflg, std::string out)
{
    IWL Buf(psio, itap, 0.0, 1, 1);
    if ((nfzc == 0) && (nfzv == 0))
        Buf.read_all(ints, ioff, ioff, 0, ioff, printflg, out);
    else
        Buf.read_all_active(ints, ioff, ioff, 0, ioff, nfzc, norbs-nfzv-1, printflg, out);
}

/*!
** iwl_rdtwo(): read two electron ints from the given file.
** The "iwl" stands for "integrals with labels," and this is the proposed
** new standard for storing two-electron integrals and their (absolute)
** orbital labels.
**
**    \param itap     = unit to read from
**    \param ints     = two electron integrals (already allocated)
**    \param ioff     = the old ioff array for lexical ordering
**    \param norbs    = number of orbitals
**    \param nfzc     = number of frozen core orbitals
**    \param nfzv     = number of frozen virtual orbitals
**    \param printflg = print integrals as they're read
**    \param out  = output file pointer
**
** David Sherrill, 1995
** \ingroup IWL
*/
void iwl_rdtwo(int itap, double *ints, int *ioff, int norbs,
      int nfzc, int nfzv, int printflg, std::string out)
{
  struct iwlbuf Buf;

  iwl_buf_init(&Buf, itap, 0.0, 1, 1);
  if ((nfzc == 0) && (nfzv == 0))
    iwl_buf_rd_all(&Buf, ints, ioff, ioff, 0, ioff, printflg, out);
  else
    iwl_buf_rd_all_act(&Buf, ints, ioff, ioff, 0, ioff, nfzc, norbs-nfzv-1,
                       printflg, out);
  iwl_buf_close(&Buf, 1);
}

}
