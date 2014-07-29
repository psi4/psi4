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

namespace psi {

int DPD::file2_mat_print(dpdfile2 *File, FILE *outfile)
{
    div_t fraction;
    int i,j;
    int rows, cols, cols_per_page, num_pages, last_page, page, first_col;
    int h, my_irrep;
    dpdparams2 *Params;

    Params = File->params;
    my_irrep = File->my_irrep;

    cols_per_page = 9;

    for(h=0; h < File->params->nirreps; h++) {

        psi::fprintf(outfile, "\n\tFile %3d DPD File2: %s\n", File->filenum,
                File->label);
        psi::fprintf(outfile,   "\tMatrix for Irrep %1d\n", h);
        psi::fprintf(outfile,   "\t----------------------------------------\n");

        rows = Params->rowtot[h];
        cols = Params->coltot[h^my_irrep];

        /* Determine the number of cols_per_page groups */
        fraction = div(cols,cols_per_page);
        num_pages = fraction.quot;  /* Number of complete column groups */
        last_page = fraction.rem;  /* Number of columns in last group */

        /* Loop over the complete column groups */
        for(page=0; page < num_pages; page++) {
            first_col = page*cols_per_page;

            psi::fprintf(outfile,"\n            ");
            for(i=first_col; i < first_col+cols_per_page; i++)
                psi::fprintf(outfile,"         %5d     ",i);

            psi::fprintf(outfile,"\n            ");
            for(i=first_col; i < first_col+cols_per_page; i++)
                psi::fprintf(outfile,"          (%3d)    ",
                        Params->colorb[h^my_irrep][i]);

            psi::fprintf (outfile,"\n");
            for(i=0; i < rows; i++) {
                psi::fprintf(outfile,"\n%5d  (%3d)",i, Params->roworb[h][i]);

                for(j=first_col; j < first_col+cols_per_page; j++)
                    psi::fprintf (outfile,"%19.15f",File->matrix[h][i][j]);
            }

            psi::fprintf (outfile,"\n");
        }

        /* Now print the remaining columns */
        if(last_page) {
            first_col = page*cols_per_page;

            psi::fprintf(outfile,"\n            ");
            for(i=first_col; i < first_col+last_page; i++)
                psi::fprintf(outfile,"         %5d     ",i);

            psi::fprintf(outfile,"\n            ");
            for(i=first_col; i < first_col+last_page; i++)
                psi::fprintf(outfile,"          (%3d)    ",
                        Params->colorb[h^my_irrep][i]);

            psi::fprintf (outfile,"\n");
            for(i=0; i < rows; i++) {
                psi::fprintf(outfile,"\n%5d  (%3d)",i, Params->roworb[h][i]);

                for(j=first_col; j < first_col+last_page; j++)
                    psi::fprintf (outfile,"%19.15f", File->matrix[h][i][j]);
            }
            psi::fprintf (outfile,"\n");
        }
    }

    return 0;
}

}
