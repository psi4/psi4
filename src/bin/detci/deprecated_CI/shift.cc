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
**  \ingroup DETCI
**  \brief Shifts SCF eigenvalues corresponding to SOCC orbitals for ZAPTn 
**
** Shifts SCF eigenvalues corresponding to SOCC orbitals for ZAPTn
** Alpha eigenvalues used for occupied SOCC (sigma+), 
** beta for unoccupied (sigma-)
**
** Done so that existing HD_AVG = ORB_ENER MPn machinery can be used unaltered.
**
** Steven E. Wheeler
** CCQC/UCLA
** 12/23/2007
*/

#include <cstdio>
#include <cstdlib>
#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "structs.h"
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define EXTERN
#include "globals.h"

namespace psi { namespace detci {

void zapt_shift(double *TEI, int nirreps, int nmo, int *doccpi, int *soccpi, 
   int *orbspi, int *frzdoccpi, int *reorder)
{
    int h1, h2;
    int x, y, i, j;
    int offset, offset2;
    int ij, ijij;
    int docc, socc;
    int docc2, socc2;
    int totfzc;

    for(h1=0,totfzc=0; h1<nirreps; h1++)
        totfzc += frzdoccpi[h1];

    for(h1 = 0,offset = 0; h1 < nirreps; h1++) {
        if(h1>0) offset += orbspi[h1-1];
        docc = doccpi[h1];
        socc = soccpi[h1];
        for(x = offset+docc; x<offset+docc+socc; x++)
            for(h2 = 0,offset2 = 0; h2 < nirreps; h2++) {
                if(h2>0) offset2 += orbspi[h2-1];
                docc2 = doccpi[h2];
                socc2 = soccpi[h2];
                for(y=offset2+docc2;y<offset2+docc2+socc2;y++) {
                    i = reorder[x] - totfzc;
                    j = reorder[y] - totfzc;
                    ij = INDEX(i,j);
                    ijij = ioff[ij] + ij;
                    CalcInfo.scfeigvala[i+totfzc] -= 0.5*TEI[ijij];
                    CalcInfo.scfeigvalb[i+totfzc] += 0.5*TEI[ijij];
                }
        }
    }
}

}} // namespace psi::detci

