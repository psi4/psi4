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
#include "psi4/libciomr/libciomr.h"
#include "dpd.h"

namespace psi {

DPD::~DPD()
{
    int h,i,j,k,cnt;
    /*  dpd_file2_cache_print(stdout); */
    file2_cache_close();
    /*  dpd_file4_cache_print(stdout);*/
    file4_cache_close();

    if(params4)
        for(i=0; i < num_pairs; i++)
            for(j=0; j < num_pairs; j++)
                free_int_matrix(params4[i][j].start13);

    if(orboff){
        for(i=0; i < num_subspaces; i++)
            free(orboff[i]);
        free(orboff);
    }

    if(pairidx && pairorb){
        for(i=0; i < num_subspaces; i++) {
            for(j=0; j < 5; j++) {
                free_int_matrix(pairidx[5*i+j]);
                for(k=0; k < nirreps; k++)
                    if(pairtot[5*i+j][k])
                        free_int_matrix(pairorb[5*i+j][k]);
                free(pairorb[5*i+j]);
            }
        }
        for(i=0,cnt=5*num_subspaces; i < num_subspaces; i++) {
            for(j=i+1; j < num_subspaces; j++,cnt+=2) {
                free_int_matrix(pairidx[cnt]);
                free_int_matrix(pairidx[cnt+1]);
                for(k=0; k < nirreps; k++) {
                    if(pairtot[cnt][k])
                        free_int_matrix(pairorb[cnt][k]);
                    if(pairtot[cnt+1][k])
                        free_int_matrix(pairorb[cnt+1][k]);
                }
                free(pairorb[cnt]);
                free(pairorb[cnt+1]);
            }
        }
        free(pairidx);
        free(pairorb);
    }

    if(orbs2 && orbidx2){
        for(i=0; i < num_subspaces; i++) {
            free(orbidx2[i]);
            for(j=0; j < nirreps; j++) {
                if(orbspi[i][j])
                    free(orbs2[i][j]);
            }
            free(orbs2[i]);
        }
        free(orbidx2);
        free(orbs2);
    }

    if(orbspi and orbsym){
        for(i=0; i < num_subspaces; i++) {
            free(orbspi[i]);
            free(orbsym[i]);
        }
        free(orbspi);
        free(orbsym);
    }

    if(pairtot)
        free_int_matrix(pairtot);

    if(numorbs)
        free(numorbs);

    if(params4){
        for(i=0; i < num_pairs; i++)
            free(params4[i]);
        free(params4);
    }
    if(params2){
        for(i=0; i < num_subspaces; i++)
            free(params2[i]);
        free(params2);
    }

    /*
    printf("memory = %d; memfree = %d\n",
    dpd_main.memory, dpd_main.memfree);
  */
}

}
