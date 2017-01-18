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
#include "psi4/libparallel/ParallelPrinter.h"
namespace psi {

/* dpd_file2_print(): Prints out data for all irreps of a two-index dpdfile.
**
** Arguments:
**   struct dpdfile2 *File: A pointer to the dpdfile to be printed.
**   std::string OutFileRMR: The formatted output file stream.
*/

int DPD::file2_print(dpdfile2 *File, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
   int i, my_irrep;
    dpdparams2 *Params;

    my_irrep = File->my_irrep;
    Params = File->params;

    printer->Printf( "\n\tDPD File2: %s\n", File->label);
    printer->Printf(   "\tDPD Parameters:\n");
    printer->Printf(   "\t------------------\n");
    printer->Printf(   "\tpnum = %d   qnum = %d   irrep = %d \n",
            Params->pnum, Params->qnum, File->my_irrep);
    printer->Printf(   "\tIrreps = %1d\n\n", Params->nirreps);
    printer->Printf( "\t   Row and column dimensions for DPD Block:\n");
    printer->Printf( "\t   ----------------------------------------\n");
    for(i=0; i < Params->nirreps; i++)
        printer->Printf(   "\t   Irrep: %1d row = %5d\tcol = %5d\n", i,
                Params->rowtot[i], Params->coltot[i^my_irrep]);


    file2_mat_init(File);
    file2_mat_rd(File);
    file2_mat_print(File, "outfile");
    file2_mat_close(File);

    return 0;

}

}
