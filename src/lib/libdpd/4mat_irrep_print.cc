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
    \ingroup DPD
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "dpd.h"
#include "psi4-dec.h"
namespace psi {

int DPD::mat4_irrep_print(double **matrix, dpdparams4 *Params,
                     int block, int my_irrep, FILE *outfile)
{
    div_t fraction;
    int i,j,r_irrep;
    int rows, cols, cols_per_page, num_pages, last_page, page, first_col;

    cols_per_page = 5;

    r_irrep = block^my_irrep;

    rows = Params->rowtot[block];
    cols = Params->coltot[r_irrep];

    /* Determine the number of cols_per_page groups */
    fraction = div(cols,cols_per_page);
    num_pages = fraction.quot;  /* Number of complete column groups */
    last_page = fraction.rem;  /* Number of columns in last group */

    /* Loop over the complete column groups */
    for(page=0; page < num_pages; page++) {
        first_col = page*cols_per_page;

        psi::fprintf(outfile,"\n           ");
        for(i=first_col; i < first_col+cols_per_page; i++)
            psi::fprintf(outfile,"              %5d",i);

        psi::fprintf(outfile,"\n               ");
        for(i=first_col; i < first_col+cols_per_page; i++)
            psi::fprintf(outfile,"          (%3d,%3d)",
                    Params->colorb[r_irrep][i][0], Params->colorb[r_irrep][i][1]);

        psi::fprintf (outfile,"\n");
        for(i=0; i < rows; i++) {
            psi::fprintf(outfile,"\n%5d  (%3d,%3d)",i,
                    Params->roworb[block][i][0], Params->roworb[block][i][1]);

            for(j=first_col; j < first_col+cols_per_page; j++)
                psi::fprintf (outfile,"%19.15f",matrix[i][j]);
        }

        psi::fprintf (outfile,"\n");
    }

    /* Now print the remaining columns */
    if(last_page) {
        first_col = page*cols_per_page;

        psi::fprintf(outfile,"\n           ");
        for(i=first_col; i < first_col+last_page; i++)
            psi::fprintf(outfile,"              %5d",i);

        psi::fprintf(outfile,"\n               ");
        for(i=first_col; i < first_col+last_page; i++)
            psi::fprintf(outfile,"          (%3d,%3d)",
                    Params->colorb[r_irrep][i][0], Params->colorb[r_irrep][i][1]);

        psi::fprintf (outfile,"\n");
        for(i=0; i < rows; i++) {
            psi::fprintf(outfile,"\n%5d  (%3d,%3d)",i,
                    Params->roworb[block][i][0], Params->roworb[block][i][1]);

            for(j=first_col; j < first_col+last_page; j++)
                psi::fprintf (outfile,"%19.15f",matrix[i][j]);
        }

        psi::fprintf (outfile,"\n");
    }

    return 0;

}

}

