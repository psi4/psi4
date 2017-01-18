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
#include <cstdarg>
#include <vector>
#include "psi4/libciomr/libciomr.h"
#include "dpd.h"
#include "psi4/libpsi4util/exception.h"

namespace psi {

DPD *global_dpd_ = NULL;
int dpd_default = 0;
DPD* dpd_list[2] = {NULL, NULL};
dpd_gbl dpd_main;

struct dpdpair{
    int *left_orbspi;
    int *right_orbspi;
    int *left_orbsym;
    int *right_orbsym;
    int *left_orboff;
    int *right_orboff;
    int permlr;
    int ler;
};

int dpd_set_default(int dpd_num)
{
  dpd_default = dpd_num;
  global_dpd_ = dpd_list[dpd_num];
  return 0;
}

extern int dpd_init(int dpd_num, int nirreps, long int memory, int cachetype,
            int *cachefiles, int **cachelist, dpd_file4_cache_entry *priority,
            int num_subspaces, std::vector<int*> &spaceArrays)
{
    if(dpd_list[dpd_num])
        throw PSIEXCEPTION("Attempting to initilize new DPD instance before the old one was freed.");
    dpd_list[dpd_num] = new DPD(dpd_num, nirreps, memory, cachetype, cachefiles, cachelist,
                                priority, num_subspaces, spaceArrays);
    dpd_default = dpd_num;
    global_dpd_ = dpd_list[dpd_num];
    return 0;
}

extern int dpd_close(int dpd_num)
{
    if(dpd_list[dpd_num] == 0)
        throw PSIEXCEPTION("Attempting to close a non-existent DPD instance.");
    delete dpd_list[dpd_num];
    dpd_list[dpd_num] = 0;

    return 0;
}

extern long int dpd_memfree(void)
{
  return dpd_main.memory - (dpd_main.memused -
                dpd_main.memcache +
                dpd_main.memlocked);
}

extern void dpd_memset(long int memory)
{
  dpd_main.memory = memory;
}

DPD::DPD():
    nirreps(0),
    num_subspaces(0),
    num_pairs(0),
    numorbs(0),
    orboff(0),
    pairtot(0),
    orbspi(0),
    orbsym(0),
    orbidx2(0),
    pairidx(0),
    orbs2(0),
    pairorb(0),
    params2(0),
    params4(0)
{}

DPD::DPD(int dpd_num, int nirreps, long int memory, int cachetype,
         int *cachefiles, int **cachelist, dpd_file4_cache_entry *priority,
         int num_subspaces, std::vector<int*> &spaceArrays)
{
    init(dpd_num, nirreps, memory, cachetype, cachefiles, cachelist,
             priority, num_subspaces, spaceArrays);
}

/* Another constructor to use new DPDMOSpaces for determining pair indices directly from strings */
DPD::DPD(int dpd_num, int nirreps, long int memory, int cachetype,
             int *cachefiles, int **cachelist,
             dpd_file4_cache_entry *priority, int num_subspaces, std::vector<DPDMOSpace> &spaces)
{
  std::vector<int*> spaceArrays;
  int *tmparray;

  for(int i=0; i < num_subspaces; i++) {
    tmparray = init_int_array(nirreps);
    for(int j=0; j < spaces[i].nIrrep(); j++) tmparray[j] = spaces[i].orbPI()[j];
    spaceArrays.push_back(tmparray);

    tmparray = init_int_array(spaces[i].nOrb());
    for(int j=0; j < spaces[i].nOrb(); j++) tmparray[j] = spaces[i].orbSym()[j];
    spaceArrays.push_back(tmparray);

    moSpaces.push_back(spaces[i]);
  }

  init(dpd_num, nirreps, memory, cachetype, cachefiles, cachelist, priority, num_subspaces, spaceArrays);
}

/* This is the original function call, but is now just a wrapper to the same function
 * that takes the spaces in a vector instead of using variable argument lists */
int DPD::init(int dpd_num, int nirreps, long int memory, int cachetype,
             int *cachefiles, int **cachelist,
             dpd_file4_cache_entry *priority, int num_subspaces, ...)
{
    std::vector<int*> spaceArrays;
    va_list ap;
    int *tmparray;

    va_start(ap, num_subspaces);
    for(int i=0; i < num_subspaces; i++) {
        tmparray = va_arg(ap, int *);
        spaceArrays.push_back(tmparray);
        tmparray = va_arg(ap, int *);
        spaceArrays.push_back(tmparray);
    }
    va_end(ap);
    return init(dpd_num, nirreps, memory, cachetype, cachefiles,
                    cachelist, priority, num_subspaces, spaceArrays);
}

/* This is the original function code, but modified to take a vector of the orbital
 * space information arrays, rather than a variable argument list; the former is
 * easier to construct for some code that creates an arbitrary number of spaces */
int DPD::init(int dpd_num_in, int nirreps_in, long int memory_in, int cachetype_in,
                  int *cachefiles_in, int **cachelist_in, dpd_file4_cache_entry *priority_in,
                  int num_subspaces_in, std::vector<int*> &spaceArrays_in)
{
    int h,h0,h1,cnt,***dp,l_irrep,r_irrep,p,q;
    int i,j,k,l,*count,offset1,offset2;
    int *tmparray;
    dpdpair *pairs;
    int nump, nrows, Gp, offset;

    nirreps = nirreps_in;
    num_subspaces = num_subspaces_in;

    dpd_main.memory = memory_in/sizeof(double);  /* Available memory in doubles */
    dpd_main.memused = 0; /* At first... */
    dpd_main.memcache = 0; /* At first... */
    dpd_main.memlocked = 0; /* At first... */

    dpd_main.cachetype = cachetype_in;
    dpd_main.cachelist = cachelist_in;
    dpd_main.cachefiles = cachefiles_in;
    dpd_main.file4_cache_priority = priority_in;

    /* Construct binary direct product array */
    dp = (int ***) malloc(nirreps * sizeof(int **));
    for(h=0; h < nirreps; h++) {
        dp[h] = init_int_matrix(nirreps,2);
        cnt=0;
        for(h0=0; h0 < nirreps; h0++) {
            for(h1=0; h1 < nirreps; h1++) {
                if((h0^h1)==h) {
                    dp[h][cnt][0] = h0;
                    dp[h][cnt++][1] = h1;
                }
            }
        }
    }

    /* Grab the irrep population and orbital symmetry arrays from the arg list */
    orbspi = (int **) malloc(sizeof(int *) * num_subspaces);
    orbsym = (int **) malloc(sizeof(int *) * num_subspaces);
    numorbs = (int *) malloc(num_subspaces * sizeof(int));
    for(i=0; i < num_subspaces; i++) {
        orbspi[i] = (int *) malloc(sizeof(int) * nirreps);
        tmparray = spaceArrays_in[2*i];
        for(j=0; j < nirreps; j++) orbspi[i][j] = tmparray[j];

        /* Compute the number of orbitals in this subspace */
        numorbs[i] = 0;
        for(h=0; h < nirreps; h++)
            numorbs[i] += orbspi[i][h];

        orbsym[i] = (int *) malloc(sizeof(int) * numorbs[i]);
        tmparray = spaceArrays_in[2*i+1];
        for(j=0; j < numorbs[i]; j++) orbsym[i][j] = tmparray[j];
    }

    /* Compute the orbital offset arrays */
    orboff = (int **) malloc(num_subspaces * sizeof(int *));
    for(i=0; i < num_subspaces; i++) {
        orboff[i] = init_int_array(nirreps);
        for(j=1; j < nirreps; j++)
            orboff[i][j] = orboff[i][j-1] + orbspi[i][j-1];
    }

    /* Compute the number of bra or ket index combinations */
    num_pairs = (num_subspaces * (num_subspaces - 1)) + (5 * num_subspaces);

    /* Set up the pair structs for later use */
    pairs = (dpdpair *) malloc(num_pairs * sizeof(dpdpair));

    /* Build the row/column dimension arrays */
    pairtot = init_int_matrix(num_pairs, nirreps);
    /* Loop over the groups of five "diagonal" pairs */
    for(i=0; i < num_subspaces; i++) {

        pairs[5*i].left_orbspi = orbspi[i];
        pairs[5*i].left_orbsym = orbsym[i];
        pairs[5*i].left_orboff = orboff[i];
        pairs[5*i].right_orbspi = orbspi[i];
        pairs[5*i].right_orbsym = orbsym[i];
        pairs[5*i].right_orboff = orboff[i];
        pairs[5*i].permlr = 0;
        pairs[5*i].ler = 0;

        pairs[5*i+1].left_orbspi = orbspi[i];
        pairs[5*i+1].left_orbsym = orbsym[i];
        pairs[5*i+1].left_orboff = orboff[i];
        pairs[5*i+1].right_orbspi = orbspi[i];
        pairs[5*i+1].right_orbsym = orbsym[i];
        pairs[5*i+1].right_orboff = orboff[i];
        pairs[5*i+1].permlr = 1;
        pairs[5*i+1].ler = 0;

        pairs[5*i+2].left_orbspi = orbspi[i];
        pairs[5*i+2].left_orbsym = orbsym[i];
        pairs[5*i+2].left_orboff = orboff[i];
        pairs[5*i+2].right_orbspi = orbspi[i];
        pairs[5*i+2].right_orbsym = orbsym[i];
        pairs[5*i+2].right_orboff = orboff[i];
        pairs[5*i+2].permlr = -1;
        pairs[5*i+2].ler = 0;

        pairs[5*i+3].left_orbspi = orbspi[i];
        pairs[5*i+3].left_orbsym = orbsym[i];
        pairs[5*i+3].left_orboff = orboff[i];
        pairs[5*i+3].right_orbspi = orbspi[i];
        pairs[5*i+3].right_orbsym = orbsym[i];
        pairs[5*i+3].right_orboff = orboff[i];
        pairs[5*i+3].permlr = 1;
        pairs[5*i+3].ler = 1;

        pairs[5*i+4].left_orbspi = orbspi[i];
        pairs[5*i+4].left_orbsym = orbsym[i];
        pairs[5*i+4].left_orboff = orboff[i];
        pairs[5*i+4].right_orbspi = orbspi[i];
        pairs[5*i+4].right_orbsym = orbsym[i];
        pairs[5*i+4].right_orboff = orboff[i];
        pairs[5*i+4].permlr = -1;
        pairs[5*i+4].ler = 1;

        for(j=0; j < nirreps; j++)
            for(k=0; k < nirreps; k++) {
                l_irrep = dp[j][k][0]; r_irrep = dp[j][k][1];

                /* orbspi,orbspi */
                pairtot[5*i][j] += orbspi[i][l_irrep] * orbspi[i][r_irrep];

                if(l_irrep > r_irrep) {
                    /* orbspi < orbspi, +1 */
                    pairtot[5*i+1][j] += orbspi[i][l_irrep] * orbspi[i][r_irrep];
                    /* orbspi < orbspi, -1 */
                    pairtot[5*i+2][j] += orbspi[i][l_irrep] * orbspi[i][r_irrep];
                    /* orbspi <= orbspi, +1 */
                    pairtot[5*i+3][j] += orbspi[i][l_irrep] * orbspi[i][r_irrep];
                    /* orbspi <= orbspi, -1 */
                    pairtot[5*i+4][j] += orbspi[i][l_irrep] * orbspi[i][r_irrep];
                }
                else if(l_irrep == r_irrep) {
                    /* orbspi < orbspi, +1 */
                    pairtot[5*i+1][j] +=
                            (orbspi[i][l_irrep] * (orbspi[i][l_irrep]-1))/2;
                    /* orbspi < orbspi, -1 */
                    pairtot[5*i+2][j] +=
                            (orbspi[i][l_irrep] * (orbspi[i][l_irrep]-1))/2;
                    /* orbspi <= orbspi, +1 */
                    pairtot[5*i+3][j] +=
                            (orbspi[i][l_irrep] * (orbspi[i][l_irrep]+1))/2;
                    /* orbspi <= orbspi, -1 */
                    pairtot[5*i+4][j] +=
                            (orbspi[i][l_irrep] * (orbspi[i][l_irrep]+1))/2;
                }
            }
    }

    /* Loop over the remaining "off diagonal" pairs */
    for(i=0,cnt=5*num_subspaces; i < num_subspaces; i++)
        for(j=i+1; j < num_subspaces; j++,cnt+=2) {

            pairs[cnt].left_orbspi = orbspi[i];
            pairs[cnt].left_orbsym = orbsym[i];
            pairs[cnt].left_orboff = orboff[i];
            pairs[cnt].right_orbspi = orbspi[j];
            pairs[cnt].right_orbsym = orbsym[j];
            pairs[cnt].right_orboff = orboff[j];
            pairs[cnt].permlr = 0;
            pairs[cnt].ler = 0;

            pairs[cnt+1].left_orbspi = orbspi[j];
            pairs[cnt+1].left_orbsym = orbsym[j];
            pairs[cnt+1].left_orboff = orboff[j];
            pairs[cnt+1].right_orbspi = orbspi[i];
            pairs[cnt+1].right_orbsym = orbsym[i];
            pairs[cnt+1].right_orboff = orboff[i];
            pairs[cnt+1].permlr = 0;
            pairs[cnt+1].ler = 0;

            for(k=0; k < nirreps; k++)
                for(l=0; l < nirreps; l++) {

                    l_irrep = dp[k][l][0]; r_irrep = dp[k][l][1];

                    /* orbspi[i],orbspi[j] */
                    pairtot[cnt][k] += orbspi[i][l_irrep] * orbspi[j][r_irrep];
                    /* orbspi[j],orbspi[i] */
                    pairtot[cnt+1][k] += orbspi[j][l_irrep] * orbspi[i][r_irrep];

                }
        }

    /* Temporary check until I'm sure I'm doing this right */
    if(num_pairs != cnt) { printf("Error in dpd_init()!\n"); exit(PSI_RETURN_FAILURE); }

    /* Build the row/column index lookup arrays */
    pairidx = (int ***) malloc(num_pairs * sizeof(int **));
    pairorb = (int ****) malloc(num_pairs * sizeof(int ***));
    count = init_int_array(nirreps);
    /* Loop over the groups of five "diagonal" pairs */
    for(i=0; i < num_subspaces; i++) {

        for(l=0; l < 5; l++) {
            pairidx[5*i+l] = init_int_matrix(numorbs[i],numorbs[i]);
            for(j=0; j < numorbs[i]; j++)
                for(k=0; k < numorbs[i]; k++)
                    pairidx[5*i+l][j][k] = -1;

            pairorb[5*i+l] = (int ***) malloc(nirreps * sizeof(int **));
            for(j=0; j < nirreps; j++) {
                pairorb[5*i+l][j] =
                        pairtot[5*i+l][j] ? init_int_matrix(pairtot[5*i+l][j],2) : NULL;
                for(k=0; k < pairtot[5*i+l][j]; k++) {
                    pairorb[5*i+l][j][k][0] = -1;
                    pairorb[5*i+l][j][k][1] = -1;
                }
            }
        }

        zero_int_array(count,nirreps);

        /* orbspi[i],orbspi[i] */
        for(j=0; j < nirreps; j++)
            for(k=0; k < nirreps; k++) {
                h0 = dp[j][k][0]; h1 = dp[j][k][1];
                offset1 = orboff[i][h0];  offset2 = orboff[i][h1];
                for(p=0; p < orbspi[i][h0]; p++)
                    for(q=0; q < orbspi[i][h1]; q++) {
                        pairorb[5*i][j][count[j]][0] = p+offset1;
                        pairorb[5*i][j][count[j]][1] = q+offset2;
                        pairidx[5*i][p+offset1][q+offset2] = count[j]++;
                    }
            }

        zero_int_array(count, nirreps);

        /* orbspi[i] < orbspi[i], +1 */
        for(j=0; j < nirreps; j++)
            for(k=0; k < nirreps; k++) {
                h0 = dp[j][k][0]; h1 = dp[j][k][1];
                offset1 = orboff[i][h0]; offset2 = orboff[i][h1];
                if(h0 == h1) {
                    for(p=0; p < orbspi[i][h0]; p++)
                        for(q=0; q < p; q++) {
                            pairorb[5*i+1][j][count[j]][0] = p+offset1;
                            pairorb[5*i+1][j][count[j]][1] = q+offset2;
                            pairidx[5*i+1][p+offset1][q+offset2] = count[j];
                            pairidx[5*i+1][q+offset2][p+offset1] = count[j]++;
                        }
                }
                else if(h0 > h1) {
                    for(p=0; p < orbspi[i][h0]; p++)
                        for(q=0; q < orbspi[i][h1]; q++) {
                            pairorb[5*i+1][j][count[j]][0] = p+offset1;
                            pairorb[5*i+1][j][count[j]][1] = q+offset2;
                            pairidx[5*i+1][p+offset1][q+offset2] = count[j];
                            pairidx[5*i+1][q+offset2][p+offset1] = count[j]++;
                        }
                }
            }

        zero_int_array(count, nirreps);

        /* orbspi[i] < orbspi[i], -1 */
        for(j=0; j < nirreps; j++)
            for(k=0; k < nirreps; k++) {
                h0 = dp[j][k][0]; h1 = dp[j][k][1];
                offset1 = orboff[i][h0]; offset2 = orboff[i][h1];
                if(h0 == h1) {
                    for(p=0; p < orbspi[i][h0]; p++)
                        for(q=0; q < p; q++) {
                            pairorb[5*i+2][j][count[j]][0] = p+offset1;
                            pairorb[5*i+2][j][count[j]][1] = q+offset2;
                            pairidx[5*i+2][p+offset1][q+offset2] = count[j];
                            pairidx[5*i+2][q+offset2][p+offset1] = count[j]++;
                        }
                }
                else if(h0 > h1) {
                    for(p=0; p < orbspi[i][h0]; p++)
                        for(q=0; q < orbspi[i][h1]; q++) {
                            pairorb[5*i+2][j][count[j]][0] = p+offset1;
                            pairorb[5*i+2][j][count[j]][1] = q+offset2;
                            pairidx[5*i+2][p+offset1][q+offset2] = count[j];
                            pairidx[5*i+2][q+offset2][p+offset1] = count[j]++;
                        }
                }
            }

        zero_int_array(count, nirreps);

        /* orbspi[i] <= orbspi[i], +1 */
        for(j=0; j < nirreps; j++)
            for(k=0; k < nirreps; k++) {
                h0 = dp[j][k][0]; h1 = dp[j][k][1];
                offset1 = orboff[i][h0]; offset2 = orboff[i][h1];
                if(h0 == h1) {
                    for(p=0; p < orbspi[i][h0]; p++)
                        for(q=0; q <= p; q++) {
                            pairorb[5*i+3][j][count[j]][0] = p+offset1;
                            pairorb[5*i+3][j][count[j]][1] = q+offset2;
                            pairidx[5*i+3][p+offset1][q+offset2] = count[j];
                            pairidx[5*i+3][q+offset2][p+offset1] = count[j]++;
                        }
                }
                else if(h0 > h1) {
                    for(p=0; p < orbspi[i][h0]; p++)
                        for(q=0; q < orbspi[i][h1]; q++) {
                            pairorb[5*i+3][j][count[j]][0] = p+offset1;
                            pairorb[5*i+3][j][count[j]][1] = q+offset2;
                            pairidx[5*i+3][p+offset1][q+offset2] = count[j];
                            pairidx[5*i+3][q+offset2][p+offset1] = count[j]++;
                        }
                }
            }

        zero_int_array(count, nirreps);

        /* orbspi[i] <= orbspi[i], -1 */
        for(j=0; j < nirreps; j++)
            for(k=0; k < nirreps; k++) {
                h0 = dp[j][k][0]; h1 = dp[j][k][1];
                offset1 = orboff[i][h0]; offset2 = orboff[i][h1];
                if(h0 == h1) {
                    for(p=0; p < orbspi[i][h0]; p++)
                        for(q=0; q <= p; q++) {
                            pairorb[5*i+4][j][count[j]][0] = p+offset1;
                            pairorb[5*i+4][j][count[j]][1] = q+offset2;
                            pairidx[5*i+4][p+offset1][q+offset2] = count[j];
                            pairidx[5*i+4][q+offset2][p+offset1] = count[j]++;
                        }
                }
                else if(h0 > h1) {
                    for(p=0; p < orbspi[i][h0]; p++)
                        for(q=0; q < orbspi[i][h1]; q++) {
                            pairorb[5*i+4][j][count[j]][0] = p+offset1;
                            pairorb[5*i+4][j][count[j]][1] = q+offset2;
                            pairidx[5*i+4][p+offset1][q+offset2] = count[j];
                            pairidx[5*i+4][q+offset2][p+offset1] = count[j]++;
                        }
                }
            }

    }

    /* Loop over the remaining "off diagonal" pairs */
    for(i=0,cnt=5*num_subspaces; i < num_subspaces; i++) {
        for(j=i+1; j < num_subspaces; j++,cnt+=2) {

            pairidx[cnt] = init_int_matrix(numorbs[i],numorbs[j]);
            for(k=0; k < numorbs[i]; k++)
                for(l=0; l < numorbs[j]; l++)
                    pairidx[cnt][k][l] = -1;
            pairidx[cnt+1] = init_int_matrix(numorbs[j],numorbs[i]);
            for(k=0; k < numorbs[j]; k++)
                for(l=0; l < numorbs[i]; l++)
                    pairidx[cnt+1][k][l] = -1;

            pairorb[cnt] = (int ***) malloc(nirreps * sizeof(int **));
            pairorb[cnt+1] = (int ***) malloc(nirreps * sizeof(int **));
            for(k=0; k < nirreps; k++) {
                pairorb[cnt][k] =
                        pairtot[cnt][k] ? init_int_matrix(pairtot[cnt][k],2) : NULL;
                pairorb[cnt+1][k] =
                        pairtot[cnt+1][k] ? init_int_matrix(pairtot[cnt+1][k],2) : NULL;
                for(l=0; l < pairtot[cnt][k]; l++) {
                    pairorb[cnt][k][l][0] = -1;
                    pairorb[cnt][k][l][1] = -1;
                }
                for(l=0; l < pairtot[cnt+1][k]; l++) {
                    pairorb[cnt+1][k][l][0] = -1;
                    pairorb[cnt+1][k][l][1] = -1;
                }

            }

            zero_int_array(count, nirreps);

            for(k=0; k < nirreps; k++)
                for(l=0; l < nirreps; l++) {
                    h0 = dp[k][l][0]; h1 = dp[k][l][1];
                    offset1 = orboff[i][h0];  offset2 = orboff[j][h1];
                    for(p=0; p < orbspi[i][h0]; p++)
                        for(q=0; q < orbspi[j][h1]; q++) {
                            pairorb[cnt][k][count[k]][0] = p+offset1;
                            pairorb[cnt][k][count[k]][1] = q+offset2;
                            pairidx[cnt][p+offset1][q+offset2] = count[k]++;
                        }
                }

            zero_int_array(count, nirreps);

            for(k=0; k < nirreps; k++)
                for(l=0; l < nirreps; l++) {
                    h0 = dp[k][l][0]; h1 = dp[k][l][1];
                    offset1 = orboff[j][h0];  offset2 = orboff[i][h1];
                    for(p=0; p < orbspi[j][h0]; p++)
                        for(q=0; q < orbspi[i][h1]; q++) {
                            pairorb[cnt+1][k][count[k]][0] = p+offset1;
                            pairorb[cnt+1][k][count[k]][1] = q+offset2;
                            pairidx[cnt+1][p+offset1][q+offset2] = count[k]++;
                        }
                }
        }
    }

    /* Temporary check until I'm sure I'm doing this right */
    if(num_pairs != cnt) { printf("Error in dpd_init()!\n"); exit(PSI_RETURN_FAILURE); }

    /* Now generate the global list of DPD parameters */
    params4 = (dpdparams4 **) malloc(num_pairs*sizeof(dpdparams4 *));
    for(i=0; i < num_pairs; i++)
        params4[i] = (dpdparams4 *) malloc(num_pairs*sizeof(dpdparams4));

    for(i=0; i < num_pairs; i++) {
        for(j=0; j < num_pairs; j++) {
            params4[i][j].nirreps = nirreps_in;

            params4[i][j].pqnum = i;
            params4[i][j].rsnum = j;

            params4[i][j].rowtot = pairtot[i];
            params4[i][j].coltot = pairtot[j];

            params4[i][j].rowidx = pairidx[i];
            params4[i][j].colidx = pairidx[j];

            params4[i][j].roworb = pairorb[i];
            params4[i][j].colorb = pairorb[j];

            params4[i][j].ppi = pairs[i].left_orbspi;
            params4[i][j].qpi = pairs[i].right_orbspi;

            params4[i][j].rpi = pairs[j].left_orbspi;
            params4[i][j].spi = pairs[j].right_orbspi;

            params4[i][j].psym = pairs[i].left_orbsym;
            params4[i][j].qsym = pairs[i].right_orbsym;

            params4[i][j].rsym = pairs[j].left_orbsym;
            params4[i][j].ssym = pairs[j].right_orbsym;

            params4[i][j].poff = pairs[i].left_orboff;
            params4[i][j].qoff = pairs[i].right_orboff;

            params4[i][j].roff = pairs[j].left_orboff;
            params4[i][j].soff = pairs[j].right_orboff;

            params4[i][j].perm_pq = pairs[i].permlr;
            params4[i][j].perm_rs = pairs[j].permlr;
            params4[i][j].peq = pairs[i].ler;
            params4[i][j].res = pairs[j].ler;

        }
    }

    /* generate the start13 lookup array */
    for(i=0; i < num_pairs; i++) {
        for(j=0; j < num_pairs; j++) {

            for(h=0,nump=0; h < nirreps; h++)
                nump += params4[i][j].ppi[h];

            params4[i][j].start13 = init_int_matrix(nirreps, nump);

            for(h=0; h < nirreps; h++) { /* h = Gamma_pq */
                for(p=0; p < nump; p++)
                    params4[i][j].start13[h][p] = -1; /* error checking */
                nrows = 0;
                for(Gp=0; Gp < nirreps; Gp++) { /* Gamma_p */
                    for(p=0; p < params4[i][j].ppi[Gp]; p++) {
                        offset = params4[i][j].poff[Gp];
                        if(params4[i][j].qpi[Gp^h])
                            params4[i][j].start13[h][offset + p] = nrows;
                        nrows += params4[i][j].qpi[Gp^h];
                    } /* p */
                } /* Gp */
            } /* h */

        } /* j */
    } /* i */

    /* Now generate the global list of one-electron DPD parameters */
    orbidx2 = (int **) malloc(num_subspaces*sizeof(int *));
    for(i=0; i < num_subspaces; i++) {
        orbidx2[i] = init_int_array(numorbs[i]);
        for(j=0; j < numorbs[i]; j++)
            orbidx2[i][j] = -1;
    }

    orbs2 = (int ***) malloc(num_subspaces*sizeof(int **));
    for(i=0;i < num_subspaces; i++) {
        orbs2[i] = (int **) malloc(nirreps*sizeof(int *));
        for(j=0; j < nirreps; j++) {
            orbs2[i][j] = orbspi[i][j] ? init_int_array(orbspi[i][j]) : NULL;
            for(k=0; k < orbspi[i][j]; k++)
                orbs2[i][j][k] = -1;
        }
    }

    for(i=0; i < num_subspaces; i++) {
        zero_int_array(count, nirreps);
        for(j=0; j < nirreps; j++) {
            offset1 = orboff[i][j];
            for(p=0; p < orbspi[i][j]; p++) {
                orbs2[i][j][count[j]] = p+offset1;
                orbidx2[i][p+offset1] = count[j]++;
            }
        }
    }

    params2 = (dpdparams2 **) malloc(num_subspaces*sizeof(dpdparams2 *));
    for(i=0; i < num_subspaces; i++)
        params2[i] = (dpdparams2 *) malloc(num_subspaces*sizeof(dpdparams2));

    for(i=0,cnt=0; i < num_subspaces; i++) {
        for(j=0; j < num_subspaces; j++,cnt++) {
            params2[i][j].nirreps = nirreps;

            params2[i][j].pnum = i;
            params2[i][j].qnum = j;

            params2[i][j].rowtot = orbspi[i];
            params2[i][j].coltot = orbspi[j];

            params2[i][j].rowidx = orbidx2[i];
            params2[i][j].colidx = orbidx2[j];

            params2[i][j].roworb = orbs2[i];
            params2[i][j].colorb = orbs2[j];

            params2[i][j].ppi = orbspi[i];
            params2[i][j].qpi = orbspi[j];

            params2[i][j].poff = orboff[i];
            params2[i][j].qoff = orboff[j];

            params2[i][j].psym = orbsym[i];
            params2[i][j].qsym = orbsym[j];
        }
    }

    free(count);

    free(pairs);

    for(h=0; h < nirreps; h++)
        free_int_matrix(dp[h]);
    free(dp);

    /* Set the default DPD set to the current one */
//          dpd_set_default(dpd_num_in);

    /* Init the Cache Linked Lists */
    file2_cache_init();
    file4_cache_init();

    return 0;
}


}
