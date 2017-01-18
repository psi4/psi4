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

#include "psi4/psi4-dec.h"
#include "psi4/libparallel/ParallelPrinter.h"
#include "psi4/libqt/qt.h"
#include <cmath>
#include "adc.h"

namespace psi{ namespace adc{

// This source file is based on the sorting function in ccenergy, cceom and cis.

struct onestack{
    double value;
    int i;
    int a;
};

//void onestack_insert(struct onestack *stack, double value, int i, int a, int level, int stacklen);

void
ADCWfn::amps_write(dpdfile2 *B, int length, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
   struct onestack *t1stack;
    int Gia = B->my_irrep;

    t1stack = (struct onestack*)malloc(length*sizeof(struct onestack));
    for(int m = 0; m < length; m++) { t1stack[m].value = 0; t1stack[m].i = 0; t1stack[m].a = 0; }

    global_dpd_->file2_mat_init(B);
    global_dpd_->file2_mat_rd(B);

    int numt1 = 0;
    for(int h = 0;h < nirrep_;h++){
        numt1 += B->params->rowtot[h] * B->params->coltot[h^Gia];
        for(int i = 0;i < B->params->rowtot[h];i++){
            int I = B->params->roworb[h][i];
            for(int a = 0;a < B->params->coltot[h^Gia];a++){
                int A = B->params->colorb[h^Gia][a];
                double value = B->matrix[h][i][a];
                for(int m = 0;m < length;m++){
                    if((fabs(value)-fabs(t1stack[m].value)) > 1e-12){
                        onestack_insert(t1stack, value, I, A, m, length);
                        break;
                    }
                }
            }
        }
    }
    global_dpd_->file2_mat_close(B);

    for(int m = 0;m < ((numt1 < length) ? numt1 : length);m++){
        if(fabs(t1stack[m].value) > 1e-6){
            printer->Printf( "\t        %3d %3d %20.10f\n", t1stack[m].i, t1stack[m].a, t1stack[m].value);
        }
    }
    free(t1stack);
}

void
ADCWfn::onestack_insert(struct onestack *stack, double value, int i, int a, int level, int stacklen)
{
    struct onestack temp;

    temp = stack[level];

    stack[level].value = value;
    stack[level].i = i;
    stack[level].a = a;

    value = temp.value;
    i = temp.i;
    a = temp.a;

    for(int l = level;l < stacklen-1;l++){
        temp = stack[l+1];

        stack[l+1].value = value;
        stack[l+1].i = i;
        stack[l+1].a = a;

        value = temp.value;
        i = temp.i;
        a = temp.a;
    }
}

}} // End Namespaces
