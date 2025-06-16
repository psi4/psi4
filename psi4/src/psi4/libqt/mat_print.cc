/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include <cstdio>
#include <cstdlib>
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
/*!
** \file
** \brief Print a matrix to a file in a formatted style
** \ingroup QT
*/

namespace psi {

/*!
** mat_print(): Prints a matrix to a file in a formatted style
**
** \param matrix  = matrix to print
** \param rows    = number of rows
** \param cols    = number of columns
** \param outfile = output file pointer for printing
**
** \ingroup QT
*/
void mat_print(double **matrix, int rows, int cols, std::string out) {
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    div_t fraction;
    int i, j;
    int cols_per_page, num_pages, last_page, page, first_col;

    cols_per_page = 5;

    /* Determine the number of cols_per_page groups */
    fraction = div(cols, cols_per_page);
    num_pages = fraction.quot; /* Number of complete column groups */
    last_page = fraction.rem;  /* Number of columns in last group */

    /* Loop over the complete column groups */
    for (page = 0; page < num_pages; page++) {
        first_col = page * cols_per_page;

        outfile->Printf("\n      ");
        for (i = first_col; i < first_col + cols_per_page; i++) printer->Printf("         %5d        ", i);

        printer->Printf("\n");
        for (i = 0; i < rows; i++) {
            printer->Printf("\n%5d ", i);

            for (j = first_col; j < first_col + cols_per_page; j++) printer->Printf("%22.15f", matrix[i][j]);
        }

        printer->Printf("\n");
    }

    /* Now print the remaining columns */
    if (last_page) {
        first_col = page * cols_per_page;

        printer->Printf("\n      ");
        for (i = first_col; i < first_col + last_page; i++) printer->Printf("         %5d        ", i);

        printer->Printf("\n");
        for (i = 0; i < rows; i++) {
            printer->Printf("\n%5d ", i);

            for (j = first_col; j < first_col + last_page; j++) printer->Printf("%22.15f", matrix[i][j]);
        }

        printer->Printf("\n");
    }
}

}  // namespace psi
