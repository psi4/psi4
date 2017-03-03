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
  \brief Filter out unneeded frozen core/virt integrals
  \ingroup QT
*/

namespace psi {
	
/*!
** filter(): Filter out undesired (frozen core/virt) integrals
**
** Given a lower-triangle array of integrals in the full
** space of orbitals as well as numbers of frozen core and virtual
** orbitals, this function returns a list of integrals involving only
** active orbitals.
**
** TDC, June 2001
**
** Note: Based on the code written by CDS in the original
** iwl_rd_one_all_act() function in LIBIWL.
**
** \ingroup QT
*/

void filter(double *input, double *output, int *ioff, int norbs, 
            int nfzc, int nfzv)
{
  int i, j, ij, IJ;
  int nact;

  nact = norbs - nfzc - nfzv;

  for(i=0,ij=0; i < nact; i++) {
    for(j=0; j <= i; j++,ij++) {
      IJ = ioff[i+nfzc] + (j + nfzc);
      output[ij] = input[IJ];
    }
  }
}

}
