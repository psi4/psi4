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

namespace psi {

void IWL::write_two(PSIO *psio, int itap, int nbfso, double *ints, int *ioff,
    double toler, int printflg, std::string out)
{
    IWL Buf(psio, itap, toler, 0, 0);
    Buf.write_all(nbfso, ints, ioff, printflg, out);
    Buf.flush(1);
}

/*!
** iwl_wrttwo()
**
** Write two electron ints to output in lexical order
** The "iwl" stands for "integrals with labels," and this is the proposed
** new standard for storing two-electron integrals and their (absolute)
** orbital labels.  This function closes the output file when finished.
**
**    \param itap     = unit to write to
**    \param nbfso    = number of basis functions in symmetry orbitals
**    \param ints     = two electron integrals
**    \param ioff     = the old ioff array for lexical ordering
**    \param printflg = print flag (1 or 0)
**    \param out  =  output file
**
** Revised 6/27/96 by CDS
** \ingroup IWL
*/
void iwl_wrttwo(int itap, int nbfso, double *ints, int *ioff, double toler,
                int printflg, std::string out)
{
  struct iwlbuf Buf;

  iwl_buf_init(&Buf, itap, toler, 0, 0);
  iwl_buf_wrt_all(&Buf, nbfso, ints, ioff, printflg, out);
  iwl_buf_flush(&Buf, 1);
  iwl_buf_close(&Buf, 1);

}

}
