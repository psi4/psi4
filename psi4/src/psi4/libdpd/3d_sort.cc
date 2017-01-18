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
    \ingroup DPD
    \brief Enter brief description of file here
*/

/*! \defgroup DPD libdpd: The Direct-Product Decomposition Library */

/* sort_3d(): Sorts a 3-index array stored as a two-by-one index
** array, (ab,c), into any desired ordering.
**
** Arguments:
**   double ***Win, ***Wout:  The input and output arrays,
**   respectively.  The leading dimension is the number of irreps in
**   the point group.  The memory for these obviously must be
**   provided by the caller.
**
**   int nirreps: The number of irreps in the point group.
**
**   int *rowtot, **rowidx, ***roworb: The integer lookup arrays for,
**   respectively, the total number of rows per irrep, the row index for
**   a given pair of orbitals, and the orbitals associated with a
**   given row index.  These generally come from the original DPD
**   buf4 whose elements are being sorted here.
**
**   int *asym, *bsym: Integer lookup arrays for the orbital
**   symmetries of the a and b compound indices on the row.
**
**   int *aoff, *boff: Integer lookup arrays for the orbital offsets
**   by irrep of the a and b compound indices on the row.
**
**   int *cpi, *coff: Integer lookup arrays for the number of
**   orbitals per irrep and the orbital offset by irrep for the c
**   column index.
**
**   int **rowidx_out: Integer lookup array for the row index for a
**   given pair of orbitals for the target ordering.
**
**   enum pattern index: One of acb, cab, cba, bac, bca.
**
**   int sum: Boolean to indicate if the contents of Wout should be
**   overwritten or accumulated.
**
**  Originally written by TDC for (T) correction code, June 2001.
**  Moved to libdpd for use by (T) and CC3 codes by TDC, July 2004.
*/

#include <cstdio>
#include "dpd.h"
#include "psi4/psi4-dec.h"
namespace psi {

void DPD::sort_3d(double ***Win, double ***Wout, int nirreps, int h, int *rowtot, int **rowidx,
                  int ***roworb, int *asym, int *bsym, int *aoff, int *boff,
                  int *cpi, int *coff, int **rowidx_out, enum pattern index, int sum)
{
    int Ga, Gb, Gc;
    int Gab, Gac, Gca, Gcb, Gbc, Gba;
    int A, B, C, a, b, c;
    int ab, ac, ca, cb, bc, ba;

    switch(index) {

    case abc:
        outfile->Printf( "\ndpd_3d_sort: abc pattern is invalid.\n");
        dpd_error("3d_sort", "outfile");
        break;

    case acb:
        for(Gab=0; Gab < nirreps; Gab++) {
            Gc = h ^ Gab;

            for(ab=0; ab < rowtot[Gab]; ab++) {

                A = roworb[Gab][ab][0];
                B = roworb[Gab][ab][1];

                Ga = asym[A]; Gb = bsym[B];
                Gac = Ga ^ Gc;

                b = B - boff[Gb];

                for(c=0; c < cpi[Gc]; c++) {
                    C = coff[Gc] + c;

                    ac = rowidx_out[A][C];

                    if(sum) Wout[Gac][ac][b] += Win[Gab][ab][c];
                    else Wout[Gac][ac][b] = Win[Gab][ab][c];
                }
            }
        }

        break;

    case cab:
        for(Gab=0; Gab < nirreps; Gab++) {
            Gc = h ^ Gab;

            for(ab=0; ab < rowtot[Gab]; ab++) {

                A = roworb[Gab][ab][0];
                B = roworb[Gab][ab][1];

                Ga = asym[A]; Gb = bsym[B];
                Gca = Ga ^ Gc;

                b = B - boff[Gb];

                for(c=0; c < cpi[Gc]; c++) {
                    C = coff[Gc] + c;

                    ca = rowidx_out[C][A];

                    if(sum) Wout[Gca][ca][b] += Win[Gab][ab][c];
                    else Wout[Gca][ca][b] = Win[Gab][ab][c];
                }
            }
        }

        break;

    case cba:
        for(Gab=0; Gab < nirreps; Gab++) {
            Gc = h ^ Gab;

            for(ab=0; ab < rowtot[Gab]; ab++) {

                A = roworb[Gab][ab][0];
                B = roworb[Gab][ab][1];

                Ga = asym[A]; Gb = bsym[B];
                a = A - aoff[Ga];

                Gcb = Gc ^ Gb;

                for(c=0; c < cpi[Gc]; c++) {
                    C = coff[Gc] + c;

                    cb = rowidx_out[C][B];

                    if(sum) Wout[Gcb][cb][a] += Win[Gab][ab][c];
                    else Wout[Gcb][cb][a] = Win[Gab][ab][c];
                }
            }
        }

        break;

    case bca:
        for(Gab=0; Gab < nirreps; Gab++) {
            Gc = h ^ Gab;

            for(ab=0; ab < rowtot[Gab]; ab++) {

                A = roworb[Gab][ab][0];
                B = roworb[Gab][ab][1];

                Ga = asym[A]; Gb = bsym[B];
                a = A - aoff[Ga];

                Gbc = Gb ^ Gc;

                for(c=0; c < cpi[Gc]; c++) {
                    C = coff[Gc] + c;

                    bc = rowidx_out[B][C];

                    if(sum) Wout[Gbc][bc][a] += Win[Gab][ab][c];
                    else Wout[Gbc][bc][a] = Win[Gab][ab][c];
                }
            }
        }

        break;

    case bac:
        for(Gab=0; Gab < nirreps; Gab++) {
            Gc = h ^ Gab;
            Gba = Gab;

            for(ab=0; ab < rowtot[Gab]; ab++) {

                A = roworb[Gab][ab][0];
                B = roworb[Gab][ab][1];

                ba = rowidx_out[B][A];

                for(c=0; c < cpi[Gc]; c++) {
                    C = coff[Gc] + c;

                    if(sum) Wout[Gba][ba][c] += Win[Gab][ab][c];
                    else Wout[Gba][ba][c] = Win[Gab][ab][c];
                }
            }
        }

        break;

    }
}

}
