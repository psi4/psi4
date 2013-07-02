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

int DPD::file4_print(dpdfile4 *File, FILE *outfile)
{
    int i, h, my_irrep;
    dpdparams4 *Params;

    my_irrep = File->my_irrep;
    Params = File->params;

    fprintf(outfile, "\n\tDPD File4: %s\n", File->label);
    fprintf(outfile, "\n\tDPD Parameters:\n");
    fprintf(outfile,   "\t---------------\n");
    fprintf(outfile,   "\tpqnum = %d   rsnum = %d\n",
            Params->pqnum, Params->rsnum);
    fprintf(outfile, "\t   Row and column dimensions for DPD Block:\n");
    fprintf(outfile, "\t   ----------------------------------------\n");
    for(i=0; i < Params->nirreps; i++)
        fprintf(outfile,   "\t   Irrep: %1d row = %5d\tcol = %5d\n", i,
                Params->rowtot[i], Params->coltot[i^my_irrep]);
    fflush(outfile);

    for(h=0; h < File->params->nirreps; h++) {
        fprintf(outfile, "\n\tFile %3d DPD File4: %s\n", File->filenum,
                File->label);
        fprintf(outfile,   "\tMatrix for Irrep %1d\n", h);
        fprintf(outfile,   "\t----------------------------------------\n");
        file4_mat_irrep_init(File, h);
        file4_mat_irrep_rd(File, h);
        mat4_irrep_print(File->matrix[h], File->params, h, my_irrep, outfile);
        file4_mat_irrep_close(File, h);
    }

    return 0;

}

}
