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
#include "dpd.h"

namespace psi {

/* dpd_file2_print(): Prints out data for all irreps of a two-index dpdfile.
**
** Arguments:
**   struct dpdfile2 *File: A pointer to the dpdfile to be printed.
**   FILE *outfile: The formatted output file stream.
*/

int DPD::file2_print(dpdfile2 *File, FILE *outfile)
{
    int i, my_irrep;
    dpdparams2 *Params;

    my_irrep = File->my_irrep;
    Params = File->params;

    psi::fprintf(outfile, "\n\tDPD File2: %s\n", File->label);
    psi::fprintf(outfile,   "\tDPD Parameters:\n");
    psi::fprintf(outfile,   "\t------------------\n");
    psi::fprintf(outfile,   "\tpnum = %d   qnum = %d   irrep = %d \n",
            Params->pnum, Params->qnum, File->my_irrep);
    psi::fprintf(outfile,   "\tIrreps = %1d\n\n", Params->nirreps);
    psi::fprintf(outfile, "\t   Row and column dimensions for DPD Block:\n");
    psi::fprintf(outfile, "\t   ----------------------------------------\n");
    for(i=0; i < Params->nirreps; i++)
        psi::fprintf(outfile,   "\t   Irrep: %1d row = %5d\tcol = %5d\n", i,
                Params->rowtot[i], Params->coltot[i^my_irrep]);
    fflush(outfile);

    file2_mat_init(File);
    file2_mat_rd(File);
    file2_mat_print(File, outfile);
    file2_mat_close(File);

    return 0;

}

}
