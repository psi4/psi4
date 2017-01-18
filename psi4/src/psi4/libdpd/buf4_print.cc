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
#include "dpd.h"
#include "psi4/psi4-dec.h"
#include "psi4/libparallel/ParallelPrinter.h"
namespace psi {

/* dpd_buf4_print(): Prints out data for all irreps of a dpd
** four-index buffer.
**
** Arguments:
**   dpdbuf4 *Buf: A pointer to the dpdbuf to be printed.
**   std::string OutFileRMR: The formatted output file stream.
**   int print_data: 0 = print buf4 parameters only; 1 = print matrices
*/

int DPD::buf4_print(dpdbuf4 *Buf, std::string out, int print_data)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
   int h, i, all_buf_irrep;
    dpdparams4 *Params;

    all_buf_irrep = Buf->file.my_irrep;
    Params = Buf->params;

    outfile->Printf( "\n\tDPD Buf4 for file4: %s\n", Buf->file.label);
    outfile->Printf( "\n\tDPD Parameters:\n");
    outfile->Printf(   "\t---------------\n");
    outfile->Printf(   "\tpqnum = %d   rsnum = %d\n",
            Params->pqnum, Params->rsnum);
    outfile->Printf( "\t   Row and column dimensions for DPD Block:\n");
    outfile->Printf( "\t   ----------------------------------------\n");
    for(i=0; i < Params->nirreps; i++)
        outfile->Printf(   "\t   Irrep: %1d row = %5d\tcol = %5d\n", i,
                Params->rowtot[i], Params->coltot[i^all_buf_irrep]);


    if(print_data) {
        for(h=0; h < Buf->params->nirreps; h++) {
            outfile->Printf( "\n\tFile %3d DPD Buf4: %s\n", Buf->file.filenum, Buf->file.label);
            outfile->Printf(   "\tMatrix for Irrep %1d\n", h);
            outfile->Printf(   "\t----------------------------------------\n");
            buf4_mat_irrep_init(Buf, h);
            buf4_mat_irrep_rd(Buf, h);
            mat4_irrep_print(Buf->matrix[h], Buf->params, h, all_buf_irrep, "outfile");
            buf4_mat_irrep_close(Buf, h);
        }
        outfile->Printf("\tTrace = %20.14f\n", buf4_trace(Buf));
    }

    return 0;

}

}
