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
#include <cstdio>
#include <cstdlib>
#include "dpd.h"
#include "psi4/psi4-dec.h"
#include "psi4/libparallel/ParallelPrinter.h"
namespace psi {

int DPD::mat4_irrep_print(double **matrix, dpdparams4 *Params,
                     int block, int my_irrep, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
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

        outfile->Printf("\n           ");
        for(i=first_col; i < first_col+cols_per_page; i++)
            outfile->Printf("              %5d",i);

        outfile->Printf("\n               ");
        for(i=first_col; i < first_col+cols_per_page; i++)
            outfile->Printf("          (%3d,%3d)",
                    Params->colorb[r_irrep][i][0], Params->colorb[r_irrep][i][1]);

        outfile->Printf("\n");
        for(i=0; i < rows; i++) {
            outfile->Printf("\n%5d  (%3d,%3d)",i,
                    Params->roworb[block][i][0], Params->roworb[block][i][1]);

            for(j=first_col; j < first_col+cols_per_page; j++)
                outfile->Printf("%19.15f",matrix[i][j]);
        }

        outfile->Printf("\n");
    }

    /* Now print the remaining columns */
    if(last_page) {
        first_col = page*cols_per_page;

        outfile->Printf("\n           ");
        for(i=first_col; i < first_col+last_page; i++)
            outfile->Printf("              %5d",i);

        outfile->Printf("\n               ");
        for(i=first_col; i < first_col+last_page; i++)
            outfile->Printf("          (%3d,%3d)",
                    Params->colorb[r_irrep][i][0], Params->colorb[r_irrep][i][1]);

        outfile->Printf("\n");
        for(i=0; i < rows; i++) {
            outfile->Printf("\n%5d  (%3d,%3d)",i,
                    Params->roworb[block][i][0], Params->roworb[block][i][1]);

            for(j=first_col; j < first_col+last_page; j++)
                outfile->Printf("%19.15f",matrix[i][j]);
        }

        outfile->Printf("\n");
    }

    return 0;

}

}
