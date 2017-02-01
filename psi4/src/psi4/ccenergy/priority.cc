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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "psi4/libciomr/libciomr.h"
#include "psi4/psifiles.h"
#include "psi4/libdpd/dpd.h"
#include "ccwave.h"
namespace psi { namespace ccenergy {



void CCEnergyWavefunction::init_priority_list(void)
{

    if(!cache_priority_list_)
        throw PSIEXCEPTION("Cache priority list must be allocated before calling *CCEnergyWavefunction::init_priority_list().");
    strcpy(cache_priority_list_[  0].label, "D <ij||ab> (i>j,a>b)");
    cache_priority_list_[  0].filenum = 105;
    cache_priority_list_[  0].irrep = 0;
    cache_priority_list_[  0].pqnum =  2;
    cache_priority_list_[  0].rsnum =  7;
    cache_priority_list_[  0].priority = 231;
    cache_priority_list_[  0].next = &(cache_priority_list_[  0+1]);
    cache_priority_list_[  0].last = NULL;

    strcpy(cache_priority_list_[  1].label, "tIJAB");
    cache_priority_list_[  1].filenum = 109;
    cache_priority_list_[  1].irrep = 0;
    cache_priority_list_[  1].pqnum =  2;
    cache_priority_list_[  1].rsnum =  7;
    cache_priority_list_[  1].priority = 1163;
    cache_priority_list_[  1].next = &(cache_priority_list_[  1+1]);
    cache_priority_list_[  1].last = &(cache_priority_list_[  1-1]);

    strcpy(cache_priority_list_[  2].label, "tijab");
    cache_priority_list_[  2].filenum = 109;
    cache_priority_list_[  2].irrep = 0;
    cache_priority_list_[  2].pqnum =  2;
    cache_priority_list_[  2].rsnum =  7;
    cache_priority_list_[  2].priority = 1163;
    cache_priority_list_[  2].next = &(cache_priority_list_[  2+1]);
    cache_priority_list_[  2].last = &(cache_priority_list_[  2-1]);

    strcpy(cache_priority_list_[  3].label, "dIJAB");
    cache_priority_list_[  3].filenum = 108;
    cache_priority_list_[  3].irrep = 0;
    cache_priority_list_[  3].pqnum =  1;
    cache_priority_list_[  3].rsnum =  6;
    cache_priority_list_[  3].priority = 77;
    cache_priority_list_[  3].next = &(cache_priority_list_[  3+1]);
    cache_priority_list_[  3].last = &(cache_priority_list_[  3-1]);

    strcpy(cache_priority_list_[  4].label, "dijab");
    cache_priority_list_[  4].filenum = 108;
    cache_priority_list_[  4].irrep = 0;
    cache_priority_list_[  4].pqnum =  1;
    cache_priority_list_[  4].rsnum =  6;
    cache_priority_list_[  4].priority = 77;
    cache_priority_list_[  4].next = &(cache_priority_list_[  4+1]);
    cache_priority_list_[  4].last = &(cache_priority_list_[  4-1]);

    strcpy(cache_priority_list_[  5].label, "D <ij|ab>");
    cache_priority_list_[  5].filenum = 105;
    cache_priority_list_[  5].irrep = 0;
    cache_priority_list_[  5].pqnum =  0;
    cache_priority_list_[  5].rsnum =  5;
    cache_priority_list_[  5].priority = 535;
    cache_priority_list_[  5].next = &(cache_priority_list_[  5+1]);
    cache_priority_list_[  5].last = &(cache_priority_list_[  5-1]);

    strcpy(cache_priority_list_[  6].label, "tIjAb");
    cache_priority_list_[  6].filenum = 109;
    cache_priority_list_[  6].irrep = 0;
    cache_priority_list_[  6].pqnum =  0;
    cache_priority_list_[  6].rsnum =  5;
    cache_priority_list_[  6].priority = 1241;
    cache_priority_list_[  6].next = &(cache_priority_list_[  6+1]);
    cache_priority_list_[  6].last = &(cache_priority_list_[  6-1]);

    strcpy(cache_priority_list_[  7].label, "dIjAb");
    cache_priority_list_[  7].filenum = 108;
    cache_priority_list_[  7].irrep = 0;
    cache_priority_list_[  7].pqnum =  0;
    cache_priority_list_[  7].rsnum =  5;
    cache_priority_list_[  7].priority = 77;
    cache_priority_list_[  7].next = &(cache_priority_list_[  7+1]);
    cache_priority_list_[  7].last = &(cache_priority_list_[  7-1]);

    strcpy(cache_priority_list_[  8].label, "tauIJAB");
    cache_priority_list_[  8].filenum = 109;
    cache_priority_list_[  8].irrep = 0;
    cache_priority_list_[  8].pqnum =  2;
    cache_priority_list_[  8].rsnum =  7;
    cache_priority_list_[  8].priority = 1161;
    cache_priority_list_[  8].next = &(cache_priority_list_[  8+1]);
    cache_priority_list_[  8].last = &(cache_priority_list_[  8-1]);

    strcpy(cache_priority_list_[  9].label, "tauijab");
    cache_priority_list_[  9].filenum = 109;
    cache_priority_list_[  9].irrep = 0;
    cache_priority_list_[  9].pqnum =  2;
    cache_priority_list_[  9].rsnum =  7;
    cache_priority_list_[  9].priority = 1161;
    cache_priority_list_[  9].next = &(cache_priority_list_[  9+1]);
    cache_priority_list_[  9].last = &(cache_priority_list_[  9-1]);

    strcpy(cache_priority_list_[ 10].label, "tauIjAb");
    cache_priority_list_[ 10].filenum = 109;
    cache_priority_list_[ 10].irrep = 0;
    cache_priority_list_[ 10].pqnum =  0;
    cache_priority_list_[ 10].rsnum =  5;
    cache_priority_list_[ 10].priority = 1161;
    cache_priority_list_[ 10].next = &(cache_priority_list_[ 10+1]);
    cache_priority_list_[ 10].last = &(cache_priority_list_[ 10-1]);

    strcpy(cache_priority_list_[ 11].label, "tauIjbA");
    cache_priority_list_[ 11].filenum = 109;
    cache_priority_list_[ 11].irrep = 0;
    cache_priority_list_[ 11].pqnum =  0;
    cache_priority_list_[ 11].rsnum =  5;
    cache_priority_list_[ 11].priority = 543;
    cache_priority_list_[ 11].next = &(cache_priority_list_[ 11+1]);
    cache_priority_list_[ 11].last = &(cache_priority_list_[ 11-1]);

    strcpy(cache_priority_list_[ 12].label, "tauiJaB");
    cache_priority_list_[ 12].filenum = 109;
    cache_priority_list_[ 12].irrep = 0;
    cache_priority_list_[ 12].pqnum =  0;
    cache_priority_list_[ 12].rsnum =  5;
    cache_priority_list_[ 12].priority = 389;
    cache_priority_list_[ 12].next = &(cache_priority_list_[ 12+1]);
    cache_priority_list_[ 12].last = &(cache_priority_list_[ 12-1]);

    strcpy(cache_priority_list_[ 13].label, "tautIJAB");
    cache_priority_list_[ 13].filenum = 109;
    cache_priority_list_[ 13].irrep = 0;
    cache_priority_list_[ 13].pqnum =  2;
    cache_priority_list_[ 13].rsnum =  7;
    cache_priority_list_[ 13].priority = 931;
    cache_priority_list_[ 13].next = &(cache_priority_list_[ 13+1]);
    cache_priority_list_[ 13].last = &(cache_priority_list_[ 13-1]);

    strcpy(cache_priority_list_[ 14].label, "tautijab");
    cache_priority_list_[ 14].filenum = 109;
    cache_priority_list_[ 14].irrep = 0;
    cache_priority_list_[ 14].pqnum =  2;
    cache_priority_list_[ 14].rsnum =  7;
    cache_priority_list_[ 14].priority = 931;
    cache_priority_list_[ 14].next = &(cache_priority_list_[ 14+1]);
    cache_priority_list_[ 14].last = &(cache_priority_list_[ 14-1]);

    strcpy(cache_priority_list_[ 15].label, "tautIjAb");
    cache_priority_list_[ 15].filenum = 109;
    cache_priority_list_[ 15].irrep = 0;
    cache_priority_list_[ 15].pqnum =  0;
    cache_priority_list_[ 15].rsnum =  5;
    cache_priority_list_[ 15].priority = 931;
    cache_priority_list_[ 15].next = &(cache_priority_list_[ 15+1]);
    cache_priority_list_[ 15].last = &(cache_priority_list_[ 15-1]);

    strcpy(cache_priority_list_[ 16].label, "tiJaB");
    cache_priority_list_[ 16].filenum = 109;
    cache_priority_list_[ 16].irrep = 0;
    cache_priority_list_[ 16].pqnum =  0;
    cache_priority_list_[ 16].rsnum =  5;
    cache_priority_list_[ 16].priority = 541;
    cache_priority_list_[ 16].next = &(cache_priority_list_[ 16+1]);
    cache_priority_list_[ 16].last = &(cache_priority_list_[ 16-1]);

    strcpy(cache_priority_list_[ 17].label, "tIAJB");
    cache_priority_list_[ 17].filenum = 109;
    cache_priority_list_[ 17].irrep = 0;
    cache_priority_list_[ 17].pqnum = 10;
    cache_priority_list_[ 17].rsnum = 10;
    cache_priority_list_[ 17].priority = 693;
    cache_priority_list_[ 17].next = &(cache_priority_list_[ 17+1]);
    cache_priority_list_[ 17].last = &(cache_priority_list_[ 17-1]);

    strcpy(cache_priority_list_[ 18].label, "tiajb");
    cache_priority_list_[ 18].filenum = 109;
    cache_priority_list_[ 18].irrep = 0;
    cache_priority_list_[ 18].pqnum = 10;
    cache_priority_list_[ 18].rsnum = 10;
    cache_priority_list_[ 18].priority = 693;
    cache_priority_list_[ 18].next = &(cache_priority_list_[ 18+1]);
    cache_priority_list_[ 18].last = &(cache_priority_list_[ 18-1]);

    strcpy(cache_priority_list_[ 19].label, "tIAjb");
    cache_priority_list_[ 19].filenum = 109;
    cache_priority_list_[ 19].irrep = 0;
    cache_priority_list_[ 19].pqnum = 10;
    cache_priority_list_[ 19].rsnum = 10;
    cache_priority_list_[ 19].priority = 849;
    cache_priority_list_[ 19].next = &(cache_priority_list_[ 19+1]);
    cache_priority_list_[ 19].last = &(cache_priority_list_[ 19-1]);

    strcpy(cache_priority_list_[ 20].label, "tiaJB");
    cache_priority_list_[ 20].filenum = 109;
    cache_priority_list_[ 20].irrep = 0;
    cache_priority_list_[ 20].pqnum = 10;
    cache_priority_list_[ 20].rsnum = 10;
    cache_priority_list_[ 20].priority = 693;
    cache_priority_list_[ 20].next = &(cache_priority_list_[ 20+1]);
    cache_priority_list_[ 20].last = &(cache_priority_list_[ 20-1]);

    strcpy(cache_priority_list_[ 21].label, "tIbjA");
    cache_priority_list_[ 21].filenum = 109;
    cache_priority_list_[ 21].irrep = 0;
    cache_priority_list_[ 21].pqnum = 10;
    cache_priority_list_[ 21].rsnum = 10;
    cache_priority_list_[ 21].priority = 619;
    cache_priority_list_[ 21].next = &(cache_priority_list_[ 21+1]);
    cache_priority_list_[ 21].last = &(cache_priority_list_[ 21-1]);

    strcpy(cache_priority_list_[ 22].label, "tjAIb");
    cache_priority_list_[ 22].filenum = 109;
    cache_priority_list_[ 22].irrep = 0;
    cache_priority_list_[ 22].pqnum = 10;
    cache_priority_list_[ 22].rsnum = 10;
    cache_priority_list_[ 22].priority = 541;
    cache_priority_list_[ 22].next = &(cache_priority_list_[ 22+1]);
    cache_priority_list_[ 22].last = &(cache_priority_list_[ 22-1]);

    strcpy(cache_priority_list_[ 23].label, "C <ia||jb> (ia,bj)");
    cache_priority_list_[ 23].filenum = 104;
    cache_priority_list_[ 23].irrep = 0;
    cache_priority_list_[ 23].pqnum = 10;
    cache_priority_list_[ 23].rsnum = 11;
    cache_priority_list_[ 23].priority = 75;
    cache_priority_list_[ 23].next = &(cache_priority_list_[ 23+1]);
    cache_priority_list_[ 23].last = &(cache_priority_list_[ 23-1]);

    strcpy(cache_priority_list_[ 24].label, "WMBEJ");
    cache_priority_list_[ 24].filenum = 115;
    cache_priority_list_[ 24].irrep = 0;
    cache_priority_list_[ 24].pqnum = 10;
    cache_priority_list_[ 24].rsnum = 11;
    cache_priority_list_[ 24].priority = 1823;
    cache_priority_list_[ 24].next = &(cache_priority_list_[ 24+1]);
    cache_priority_list_[ 24].last = &(cache_priority_list_[ 24-1]);

    strcpy(cache_priority_list_[ 25].label, "Wmbej");
    cache_priority_list_[ 25].filenum = 115;
    cache_priority_list_[ 25].irrep = 0;
    cache_priority_list_[ 25].pqnum = 10;
    cache_priority_list_[ 25].rsnum = 11;
    cache_priority_list_[ 25].priority = 1823;
    cache_priority_list_[ 25].next = &(cache_priority_list_[ 25+1]);
    cache_priority_list_[ 25].last = &(cache_priority_list_[ 25-1]);

    strcpy(cache_priority_list_[ 26].label, "C <ia|jb>");
    cache_priority_list_[ 26].filenum = 104;
    cache_priority_list_[ 26].irrep = 0;
    cache_priority_list_[ 26].pqnum = 10;
    cache_priority_list_[ 26].rsnum = 10;
    cache_priority_list_[ 26].priority = 227;
    cache_priority_list_[ 26].next = &(cache_priority_list_[ 26+1]);
    cache_priority_list_[ 26].last = &(cache_priority_list_[ 26-1]);

    strcpy(cache_priority_list_[ 27].label, "WmBEj");
    cache_priority_list_[ 27].filenum = 115;
    cache_priority_list_[ 27].irrep = 0;
    cache_priority_list_[ 27].pqnum = 10;
    cache_priority_list_[ 27].rsnum = 10;
    cache_priority_list_[ 27].priority = 1823;
    cache_priority_list_[ 27].next = &(cache_priority_list_[ 27+1]);
    cache_priority_list_[ 27].last = &(cache_priority_list_[ 27-1]);

    strcpy(cache_priority_list_[ 28].label, "WMbeJ");
    cache_priority_list_[ 28].filenum = 115;
    cache_priority_list_[ 28].irrep = 0;
    cache_priority_list_[ 28].pqnum = 10;
    cache_priority_list_[ 28].rsnum = 10;
    cache_priority_list_[ 28].priority = 1823;
    cache_priority_list_[ 28].next = &(cache_priority_list_[ 28+1]);
    cache_priority_list_[ 28].last = &(cache_priority_list_[ 28-1]);

    strcpy(cache_priority_list_[ 29].label, "D <ij|ab> (ib,aj)");
    cache_priority_list_[ 29].filenum = 105;
    cache_priority_list_[ 29].irrep = 0;
    cache_priority_list_[ 29].pqnum = 10;
    cache_priority_list_[ 29].rsnum = 11;
    cache_priority_list_[ 29].priority = 227;
    cache_priority_list_[ 29].next = &(cache_priority_list_[ 29+1]);
    cache_priority_list_[ 29].last = &(cache_priority_list_[ 29-1]);

    strcpy(cache_priority_list_[ 30].label, "WMbEj");
    cache_priority_list_[ 30].filenum = 115;
    cache_priority_list_[ 30].irrep = 0;
    cache_priority_list_[ 30].pqnum = 10;
    cache_priority_list_[ 30].rsnum = 11;
    cache_priority_list_[ 30].priority = 1823;
    cache_priority_list_[ 30].next = &(cache_priority_list_[ 30+1]);
    cache_priority_list_[ 30].last = &(cache_priority_list_[ 30-1]);

    strcpy(cache_priority_list_[ 31].label, "WmBeJ");
    cache_priority_list_[ 31].filenum = 115;
    cache_priority_list_[ 31].irrep = 0;
    cache_priority_list_[ 31].pqnum = 10;
    cache_priority_list_[ 31].rsnum = 11;
    cache_priority_list_[ 31].priority = 1823;
    cache_priority_list_[ 31].next = &(cache_priority_list_[ 31+1]);
    cache_priority_list_[ 31].last = &(cache_priority_list_[ 31-1]);

    strcpy(cache_priority_list_[ 32].label, "F <ia||bc> (ia,b>c)");
    cache_priority_list_[ 32].filenum = 107;
    cache_priority_list_[ 32].irrep = 0;
    cache_priority_list_[ 32].pqnum = 10;
    cache_priority_list_[ 32].rsnum =  7;
    cache_priority_list_[ 32].priority = 455;
    cache_priority_list_[ 32].next = &(cache_priority_list_[ 32+1]);
    cache_priority_list_[ 32].last = &(cache_priority_list_[ 32-1]);

    strcpy(cache_priority_list_[ 33].label, "F <ia|bc>");
    cache_priority_list_[ 33].filenum = 107;
    cache_priority_list_[ 33].irrep = 0;
    cache_priority_list_[ 33].pqnum = 10;
    cache_priority_list_[ 33].rsnum =  5;
    cache_priority_list_[ 33].priority = 379;
    cache_priority_list_[ 33].next = &(cache_priority_list_[ 33+1]);
    cache_priority_list_[ 33].last = &(cache_priority_list_[ 33-1]);

    strcpy(cache_priority_list_[ 34].label, "E <ij||ka> (i>j,ak)");
    cache_priority_list_[ 34].filenum = 106;
    cache_priority_list_[ 34].irrep = 0;
    cache_priority_list_[ 34].pqnum =  2;
    cache_priority_list_[ 34].rsnum = 11;
    cache_priority_list_[ 34].priority = 75;
    cache_priority_list_[ 34].next = &(cache_priority_list_[ 34+1]);
    cache_priority_list_[ 34].last = &(cache_priority_list_[ 34-1]);

    strcpy(cache_priority_list_[ 35].label, "E <ai|jk>");
    cache_priority_list_[ 35].filenum = 106;
    cache_priority_list_[ 35].irrep = 0;
    cache_priority_list_[ 35].pqnum = 11;
    cache_priority_list_[ 35].rsnum =  0;
    cache_priority_list_[ 35].priority = 759;
    cache_priority_list_[ 35].next = &(cache_priority_list_[ 35+1]);
    cache_priority_list_[ 35].last = &(cache_priority_list_[ 35-1]);

    strcpy(cache_priority_list_[ 36].label, "E <ij|ka>");
    cache_priority_list_[ 36].filenum = 106;
    cache_priority_list_[ 36].irrep = 0;
    cache_priority_list_[ 36].pqnum =  0;
    cache_priority_list_[ 36].rsnum = 10;
    cache_priority_list_[ 36].priority = 151;
    cache_priority_list_[ 36].next = &(cache_priority_list_[ 36+1]);
    cache_priority_list_[ 36].last = &(cache_priority_list_[ 36-1]);

    strcpy(cache_priority_list_[ 37].label, "WMBEJ");
    cache_priority_list_[ 37].filenum = 111;
    cache_priority_list_[ 37].irrep = 0;
    cache_priority_list_[ 37].pqnum = 10;
    cache_priority_list_[ 37].rsnum = 10;
    cache_priority_list_[ 37].priority = 1899;
    cache_priority_list_[ 37].next = &(cache_priority_list_[ 37+1]);
    cache_priority_list_[ 37].last = &(cache_priority_list_[ 37-1]);

    strcpy(cache_priority_list_[ 38].label, "Wmbej");
    cache_priority_list_[ 38].filenum = 111;
    cache_priority_list_[ 38].irrep = 0;
    cache_priority_list_[ 38].pqnum = 10;
    cache_priority_list_[ 38].rsnum = 10;
    cache_priority_list_[ 38].priority = 1899;
    cache_priority_list_[ 38].next = &(cache_priority_list_[ 38+1]);
    cache_priority_list_[ 38].last = &(cache_priority_list_[ 38-1]);

    strcpy(cache_priority_list_[ 39].label, "WMbEj");
    cache_priority_list_[ 39].filenum = 111;
    cache_priority_list_[ 39].irrep = 0;
    cache_priority_list_[ 39].pqnum = 10;
    cache_priority_list_[ 39].rsnum = 10;
    cache_priority_list_[ 39].priority = 1899;
    cache_priority_list_[ 39].next = &(cache_priority_list_[ 39+1]);
    cache_priority_list_[ 39].last = &(cache_priority_list_[ 39-1]);

    strcpy(cache_priority_list_[ 40].label, "WmBeJ");
    cache_priority_list_[ 40].filenum = 111;
    cache_priority_list_[ 40].irrep = 0;
    cache_priority_list_[ 40].pqnum = 10;
    cache_priority_list_[ 40].rsnum = 10;
    cache_priority_list_[ 40].priority = 1899;
    cache_priority_list_[ 40].next = &(cache_priority_list_[ 40+1]);
    cache_priority_list_[ 40].last = &(cache_priority_list_[ 40-1]);

    strcpy(cache_priority_list_[ 41].label, "WMbeJ");
    cache_priority_list_[ 41].filenum = 111;
    cache_priority_list_[ 41].irrep = 0;
    cache_priority_list_[ 41].pqnum = 10;
    cache_priority_list_[ 41].rsnum = 10;
    cache_priority_list_[ 41].priority = 1519;
    cache_priority_list_[ 41].next = &(cache_priority_list_[ 41+1]);
    cache_priority_list_[ 41].last = &(cache_priority_list_[ 41-1]);

    strcpy(cache_priority_list_[ 42].label, "WmBEj");
    cache_priority_list_[ 42].filenum = 111;
    cache_priority_list_[ 42].irrep = 0;
    cache_priority_list_[ 42].pqnum = 10;
    cache_priority_list_[ 42].rsnum = 10;
    cache_priority_list_[ 42].priority = 1519;
    cache_priority_list_[ 42].next = &(cache_priority_list_[ 42+1]);
    cache_priority_list_[ 42].last = &(cache_priority_list_[ 42-1]);

    strcpy(cache_priority_list_[ 43].label, "D <ij||ab> (ia,jb)");
    cache_priority_list_[ 43].filenum = 105;
    cache_priority_list_[ 43].irrep = 0;
    cache_priority_list_[ 43].pqnum = 10;
    cache_priority_list_[ 43].rsnum = 10;
    cache_priority_list_[ 43].priority = 303;
    cache_priority_list_[ 43].next = &(cache_priority_list_[ 43+1]);
    cache_priority_list_[ 43].last = &(cache_priority_list_[ 43-1]);

    strcpy(cache_priority_list_[ 44].label, "D <ij|ab> (ia,jb)");
    cache_priority_list_[ 44].filenum = 105;
    cache_priority_list_[ 44].irrep = 0;
    cache_priority_list_[ 44].pqnum = 10;
    cache_priority_list_[ 44].rsnum = 10;
    cache_priority_list_[ 44].priority = 303;
    cache_priority_list_[ 44].next = &(cache_priority_list_[ 44+1]);
    cache_priority_list_[ 44].last = &(cache_priority_list_[ 44-1]);

    strcpy(cache_priority_list_[ 45].label, "Y (ME,JN)");
    cache_priority_list_[ 45].filenum = 115;
    cache_priority_list_[ 45].irrep = 0;
    cache_priority_list_[ 45].pqnum = 10;
    cache_priority_list_[ 45].rsnum =  0;
    cache_priority_list_[ 45].priority = 4103;
    cache_priority_list_[ 45].next = &(cache_priority_list_[ 45+1]);
    cache_priority_list_[ 45].last = &(cache_priority_list_[ 45-1]);

    strcpy(cache_priority_list_[ 46].label, "D <ij||ab> (ia,bj)");
    cache_priority_list_[ 46].filenum = 105;
    cache_priority_list_[ 46].irrep = 0;
    cache_priority_list_[ 46].pqnum = 10;
    cache_priority_list_[ 46].rsnum = 11;
    cache_priority_list_[ 46].priority = 151;
    cache_priority_list_[ 46].next = &(cache_priority_list_[ 46+1]);
    cache_priority_list_[ 46].last = &(cache_priority_list_[ 46-1]);

    strcpy(cache_priority_list_[ 47].label, "D <ij|ab> (ia,bj)");
    cache_priority_list_[ 47].filenum = 105;
    cache_priority_list_[ 47].irrep = 0;
    cache_priority_list_[ 47].pqnum = 10;
    cache_priority_list_[ 47].rsnum = 11;
    cache_priority_list_[ 47].priority = 151;
    cache_priority_list_[ 47].next = &(cache_priority_list_[ 47+1]);
    cache_priority_list_[ 47].last = &(cache_priority_list_[ 47-1]);

    strcpy(cache_priority_list_[ 48].label, "D <ij|ab> (ib,ja)");
    cache_priority_list_[ 48].filenum = 105;
    cache_priority_list_[ 48].irrep = 0;
    cache_priority_list_[ 48].pqnum = 10;
    cache_priority_list_[ 48].rsnum = 10;
    cache_priority_list_[ 48].priority = 303;
    cache_priority_list_[ 48].next = &(cache_priority_list_[ 48+1]);
    cache_priority_list_[ 48].last = &(cache_priority_list_[ 48-1]);

    strcpy(cache_priority_list_[ 49].label, "D <ij||ab>");
    cache_priority_list_[ 49].filenum = 105;
    cache_priority_list_[ 49].irrep = 0;
    cache_priority_list_[ 49].pqnum =  0;
    cache_priority_list_[ 49].rsnum =  5;
    cache_priority_list_[ 49].priority = 75;
    cache_priority_list_[ 49].next = &(cache_priority_list_[ 49+1]);
    cache_priority_list_[ 49].last = &(cache_priority_list_[ 49-1]);

    strcpy(cache_priority_list_[ 50].label, "D <ij||ab> (i>j,ab)");
    cache_priority_list_[ 50].filenum = 105;
    cache_priority_list_[ 50].irrep = 0;
    cache_priority_list_[ 50].pqnum =  2;
    cache_priority_list_[ 50].rsnum =  5;
    cache_priority_list_[ 50].priority = 75;
    cache_priority_list_[ 50].next = &(cache_priority_list_[ 50+1]);
    cache_priority_list_[ 50].last = &(cache_priority_list_[ 50-1]);

    strcpy(cache_priority_list_[ 51].label, "D <ij||ab> (ij,a>b)");
    cache_priority_list_[ 51].filenum = 105;
    cache_priority_list_[ 51].irrep = 0;
    cache_priority_list_[ 51].pqnum =  0;
    cache_priority_list_[ 51].rsnum =  7;
    cache_priority_list_[ 51].priority = 75;
    cache_priority_list_[ 51].next = &(cache_priority_list_[ 51+1]);
    cache_priority_list_[ 51].last = &(cache_priority_list_[ 51-1]);

    strcpy(cache_priority_list_[ 52].label, "C <ia||jb>");
    cache_priority_list_[ 52].filenum = 104;
    cache_priority_list_[ 52].irrep = 0;
    cache_priority_list_[ 52].pqnum = 10;
    cache_priority_list_[ 52].rsnum = 10;
    cache_priority_list_[ 52].priority = 227;
    cache_priority_list_[ 52].next = &(cache_priority_list_[ 52+1]);
    cache_priority_list_[ 52].last = &(cache_priority_list_[ 52-1]);

    strcpy(cache_priority_list_[ 53].label, "A <ij|kl>");
    cache_priority_list_[ 53].filenum = 102;
    cache_priority_list_[ 53].irrep = 0;
    cache_priority_list_[ 53].pqnum =  0;
    cache_priority_list_[ 53].rsnum =  0;
    cache_priority_list_[ 53].priority = 151;
    cache_priority_list_[ 53].next = &(cache_priority_list_[ 53+1]);
    cache_priority_list_[ 53].last = &(cache_priority_list_[ 53-1]);

    strcpy(cache_priority_list_[ 54].label, "WMNIJ");
    cache_priority_list_[ 54].filenum = 111;
    cache_priority_list_[ 54].irrep = 0;
    cache_priority_list_[ 54].pqnum =  2;
    cache_priority_list_[ 54].rsnum =  2;
    cache_priority_list_[ 54].priority = 8891;
    cache_priority_list_[ 54].next = &(cache_priority_list_[ 54+1]);
    cache_priority_list_[ 54].last = &(cache_priority_list_[ 54-1]);

    strcpy(cache_priority_list_[ 55].label, "Wmnij");
    cache_priority_list_[ 55].filenum = 111;
    cache_priority_list_[ 55].irrep = 0;
    cache_priority_list_[ 55].pqnum =  2;
    cache_priority_list_[ 55].rsnum =  2;
    cache_priority_list_[ 55].priority = 8891;
    cache_priority_list_[ 55].next = &(cache_priority_list_[ 55+1]);
    cache_priority_list_[ 55].last = &(cache_priority_list_[ 55-1]);

    strcpy(cache_priority_list_[ 56].label, "WMnIj");
    cache_priority_list_[ 56].filenum = 111;
    cache_priority_list_[ 56].irrep = 0;
    cache_priority_list_[ 56].pqnum =  0;
    cache_priority_list_[ 56].rsnum =  0;
    cache_priority_list_[ 56].priority = 2127;
    cache_priority_list_[ 56].next = &(cache_priority_list_[ 56+1]);
    cache_priority_list_[ 56].last = &(cache_priority_list_[ 56-1]);

    strcpy(cache_priority_list_[ 57].label, "E <ij||ka> (i>j,ka)");
    cache_priority_list_[ 57].filenum = 106;
    cache_priority_list_[ 57].irrep = 0;
    cache_priority_list_[ 57].pqnum =  2;
    cache_priority_list_[ 57].rsnum = 10;
    cache_priority_list_[ 57].priority = 75;
    cache_priority_list_[ 57].next = &(cache_priority_list_[ 57+1]);
    cache_priority_list_[ 57].last = &(cache_priority_list_[ 57-1]);

    strcpy(cache_priority_list_[ 58].label, "W (MN,IJ)");
    cache_priority_list_[ 58].filenum = 115;
    cache_priority_list_[ 58].irrep = 0;
    cache_priority_list_[ 58].pqnum =  2;
    cache_priority_list_[ 58].rsnum =  0;
    cache_priority_list_[ 58].priority = 2583;
    cache_priority_list_[ 58].next = &(cache_priority_list_[ 58+1]);
    cache_priority_list_[ 58].last = &(cache_priority_list_[ 58-1]);

    strcpy(cache_priority_list_[ 59].label, "ZIJMA");
    cache_priority_list_[ 59].filenum = 113;
    cache_priority_list_[ 59].irrep = 0;
    cache_priority_list_[ 59].pqnum =  2;
    cache_priority_list_[ 59].rsnum = 10;
    cache_priority_list_[ 59].priority = 455;
    cache_priority_list_[ 59].next = &(cache_priority_list_[ 59+1]);
    cache_priority_list_[ 59].last = &(cache_priority_list_[ 59-1]);

    strcpy(cache_priority_list_[ 60].label, "Zijma");
    cache_priority_list_[ 60].filenum = 113;
    cache_priority_list_[ 60].irrep = 0;
    cache_priority_list_[ 60].pqnum =  2;
    cache_priority_list_[ 60].rsnum = 10;
    cache_priority_list_[ 60].priority = 455;
    cache_priority_list_[ 60].next = &(cache_priority_list_[ 60+1]);
    cache_priority_list_[ 60].last = &(cache_priority_list_[ 60-1]);

    strcpy(cache_priority_list_[ 61].label, "ZIjMa");
    cache_priority_list_[ 61].filenum = 113;
    cache_priority_list_[ 61].irrep = 0;
    cache_priority_list_[ 61].pqnum =  0;
    cache_priority_list_[ 61].rsnum = 10;
    cache_priority_list_[ 61].priority = 455;
    cache_priority_list_[ 61].next = &(cache_priority_list_[ 61+1]);
    cache_priority_list_[ 61].last = &(cache_priority_list_[ 61-1]);

    strcpy(cache_priority_list_[ 62].label, "ZIjmA");
    cache_priority_list_[ 62].filenum = 113;
    cache_priority_list_[ 62].irrep = 0;
    cache_priority_list_[ 62].pqnum =  0;
    cache_priority_list_[ 62].rsnum = 10;
    cache_priority_list_[ 62].priority = 379;
    cache_priority_list_[ 62].next = &(cache_priority_list_[ 62+1]);
    cache_priority_list_[ 62].last = &(cache_priority_list_[ 62-1]);

    strcpy(cache_priority_list_[ 63].label, "ZIJAM");
    cache_priority_list_[ 63].filenum = 113;
    cache_priority_list_[ 63].irrep = 0;
    cache_priority_list_[ 63].pqnum =  2;
    cache_priority_list_[ 63].rsnum = 11;
    cache_priority_list_[ 63].priority = 455;
    cache_priority_list_[ 63].next = &(cache_priority_list_[ 63+1]);
    cache_priority_list_[ 63].last = &(cache_priority_list_[ 63-1]);

    strcpy(cache_priority_list_[ 64].label, "Zijam");
    cache_priority_list_[ 64].filenum = 113;
    cache_priority_list_[ 64].irrep = 0;
    cache_priority_list_[ 64].pqnum =  2;
    cache_priority_list_[ 64].rsnum = 11;
    cache_priority_list_[ 64].priority = 455;
    cache_priority_list_[ 64].next = &(cache_priority_list_[ 64+1]);
    cache_priority_list_[ 64].last = &(cache_priority_list_[ 64-1]);

    strcpy(cache_priority_list_[ 65].label, "ZIjAm");
    cache_priority_list_[ 65].filenum = 113;
    cache_priority_list_[ 65].irrep = 0;
    cache_priority_list_[ 65].pqnum =  0;
    cache_priority_list_[ 65].rsnum = 11;
    cache_priority_list_[ 65].priority = 455;
    cache_priority_list_[ 65].next = &(cache_priority_list_[ 65+1]);
    cache_priority_list_[ 65].last = &(cache_priority_list_[ 65-1]);

    strcpy(cache_priority_list_[ 66].label, "New tIJAB");
    cache_priority_list_[ 66].filenum = 109;
    cache_priority_list_[ 66].irrep = 0;
    cache_priority_list_[ 66].pqnum =  2;
    cache_priority_list_[ 66].rsnum =  7;
    cache_priority_list_[ 66].priority = 58571;
    cache_priority_list_[ 66].next = &(cache_priority_list_[ 66+1]);
    cache_priority_list_[ 66].last = &(cache_priority_list_[ 66-1]);

    strcpy(cache_priority_list_[ 67].label, "New tijab");
    cache_priority_list_[ 67].filenum = 109;
    cache_priority_list_[ 67].irrep = 0;
    cache_priority_list_[ 67].pqnum =  2;
    cache_priority_list_[ 67].rsnum =  7;
    cache_priority_list_[ 67].priority = 58571;
    cache_priority_list_[ 67].next = &(cache_priority_list_[ 67+1]);
    cache_priority_list_[ 67].last = &(cache_priority_list_[ 67-1]);

    strcpy(cache_priority_list_[ 68].label, "New tIjAb");
    cache_priority_list_[ 68].filenum = 109;
    cache_priority_list_[ 68].irrep = 0;
    cache_priority_list_[ 68].pqnum =  0;
    cache_priority_list_[ 68].rsnum =  5;
    cache_priority_list_[ 68].priority = 10843;
    cache_priority_list_[ 68].next = &(cache_priority_list_[ 68+1]);
    cache_priority_list_[ 68].last = &(cache_priority_list_[ 68-1]);

    strcpy(cache_priority_list_[ 69].label, "T (I>J,AB)");
    cache_priority_list_[ 69].filenum = 115;
    cache_priority_list_[ 69].irrep = 0;
    cache_priority_list_[ 69].pqnum =  2;
    cache_priority_list_[ 69].rsnum =  5;
    cache_priority_list_[ 69].priority = 7295;
    cache_priority_list_[ 69].next = &(cache_priority_list_[ 69+1]);
    cache_priority_list_[ 69].last = &(cache_priority_list_[ 69-1]);

    strcpy(cache_priority_list_[ 70].label, "T (IJ,A>B)");
    cache_priority_list_[ 70].filenum = 115;
    cache_priority_list_[ 70].irrep = 0;
    cache_priority_list_[ 70].pqnum =  0;
    cache_priority_list_[ 70].rsnum =  7;
    cache_priority_list_[ 70].priority = 4711;
    cache_priority_list_[ 70].next = &(cache_priority_list_[ 70+1]);
    cache_priority_list_[ 70].last = &(cache_priority_list_[ 70-1]);

    strcpy(cache_priority_list_[ 71].label, "B <ab||cd> (a>b,c>d)");
    cache_priority_list_[ 71].filenum = 103;
    cache_priority_list_[ 71].irrep = 0;
    cache_priority_list_[ 71].pqnum =  7;
    cache_priority_list_[ 71].rsnum =  7;
    cache_priority_list_[ 71].priority = 75;
    cache_priority_list_[ 71].next = &(cache_priority_list_[ 71+1]);
    cache_priority_list_[ 71].last = &(cache_priority_list_[ 71-1]);

    strcpy(cache_priority_list_[ 72].label, "B <ab|cd>");
    cache_priority_list_[ 72].filenum = 103;
    cache_priority_list_[ 72].irrep = 0;
    cache_priority_list_[ 72].pqnum =  5;
    cache_priority_list_[ 72].rsnum =  5;
    cache_priority_list_[ 72].priority = 75;
    cache_priority_list_[ 72].next = &(cache_priority_list_[ 72+1]);
    cache_priority_list_[ 72].last = &(cache_priority_list_[ 72-1]);

    strcpy(cache_priority_list_[ 73].label, "Z(ab,ij)");
    cache_priority_list_[ 73].filenum = 115;
    cache_priority_list_[ 73].irrep = 0;
    cache_priority_list_[ 73].pqnum =  7;
    cache_priority_list_[ 73].rsnum =  2;
    cache_priority_list_[ 73].priority = 683;
    cache_priority_list_[ 73].next = &(cache_priority_list_[ 73+1]);
    cache_priority_list_[ 73].last = &(cache_priority_list_[ 73-1]);

    strcpy(cache_priority_list_[ 74].label, "Z(ij,ab)");
    cache_priority_list_[ 74].filenum = 115;
    cache_priority_list_[ 74].irrep = 0;
    cache_priority_list_[ 74].pqnum =  2;
    cache_priority_list_[ 74].rsnum =  7;
    cache_priority_list_[ 74].priority = 911;
    cache_priority_list_[ 74].next = &(cache_priority_list_[ 74+1]);
    cache_priority_list_[ 74].last = &(cache_priority_list_[ 74-1]);

    strcpy(cache_priority_list_[ 75].label, "Z(Ab,Ij)");
    cache_priority_list_[ 75].filenum = 115;
    cache_priority_list_[ 75].irrep = 0;
    cache_priority_list_[ 75].pqnum =  5;
    cache_priority_list_[ 75].rsnum =  0;
    cache_priority_list_[ 75].priority = 379;
    cache_priority_list_[ 75].next = &(cache_priority_list_[ 75+1]);
    cache_priority_list_[ 75].last = &(cache_priority_list_[ 75-1]);

    strcpy(cache_priority_list_[ 76].label, "Z(Ij,Ab)");
    cache_priority_list_[ 76].filenum = 115;
    cache_priority_list_[ 76].irrep = 0;
    cache_priority_list_[ 76].pqnum =  0;
    cache_priority_list_[ 76].rsnum =  5;
    cache_priority_list_[ 76].priority = 455;
    cache_priority_list_[ 76].next = &(cache_priority_list_[ 76+1]);
    cache_priority_list_[ 76].last = &(cache_priority_list_[ 76-1]);

    strcpy(cache_priority_list_[ 77].label, "T (JI,A>B)");
    cache_priority_list_[ 77].filenum = 115;
    cache_priority_list_[ 77].irrep = 0;
    cache_priority_list_[ 77].pqnum =  0;
    cache_priority_list_[ 77].rsnum =  7;
    cache_priority_list_[ 77].priority = 911;
    cache_priority_list_[ 77].next = &(cache_priority_list_[ 77+1]);
    cache_priority_list_[ 77].last = &(cache_priority_list_[ 77-1]);

    strcpy(cache_priority_list_[ 78].label, "F <ai|bc>");
    cache_priority_list_[ 78].filenum = 107;
    cache_priority_list_[ 78].irrep = 0;
    cache_priority_list_[ 78].pqnum = 11;
    cache_priority_list_[ 78].rsnum =  5;
    cache_priority_list_[ 78].priority = 75;
    cache_priority_list_[ 78].next = &(cache_priority_list_[ 78+1]);
    cache_priority_list_[ 78].last = &(cache_priority_list_[ 78-1]);

    strcpy(cache_priority_list_[ 79].label, "T (I>J,BA)");
    cache_priority_list_[ 79].filenum = 115;
    cache_priority_list_[ 79].irrep = 0;
    cache_priority_list_[ 79].pqnum =  2;
    cache_priority_list_[ 79].rsnum =  5;
    cache_priority_list_[ 79].priority = 911;
    cache_priority_list_[ 79].next = &(cache_priority_list_[ 79+1]);
    cache_priority_list_[ 79].last = &(cache_priority_list_[ 79-1]);

    strcpy(cache_priority_list_[ 80].label, "E <ia|jk>");
    cache_priority_list_[ 80].filenum = 106;
    cache_priority_list_[ 80].irrep = 0;
    cache_priority_list_[ 80].pqnum = 10;
    cache_priority_list_[ 80].rsnum =  0;
    cache_priority_list_[ 80].priority = 75;
    cache_priority_list_[ 80].next = &(cache_priority_list_[ 80+1]);
    cache_priority_list_[ 80].last = &(cache_priority_list_[ 80-1]);

    strcpy(cache_priority_list_[ 81].label, "T2 (IA,JB)");
    cache_priority_list_[ 81].filenum = 115;
    cache_priority_list_[ 81].irrep = 0;
    cache_priority_list_[ 81].pqnum = 10;
    cache_priority_list_[ 81].rsnum = 10;
    cache_priority_list_[ 81].priority = 683;
    cache_priority_list_[ 81].next = &(cache_priority_list_[ 81+1]);
    cache_priority_list_[ 81].last = &(cache_priority_list_[ 81-1]);

    strcpy(cache_priority_list_[ 82].label, "T2 (IJ,AB)");
    cache_priority_list_[ 82].filenum = 115;
    cache_priority_list_[ 82].irrep = 0;
    cache_priority_list_[ 82].pqnum =  0;
    cache_priority_list_[ 82].rsnum =  5;
    cache_priority_list_[ 82].priority = 4103;
    cache_priority_list_[ 82].next = &(cache_priority_list_[ 82+1]);
    cache_priority_list_[ 82].last = &(cache_priority_list_[ 82-1]);

    strcpy(cache_priority_list_[ 83].label, "T2 (JI,AB)");
    cache_priority_list_[ 83].filenum = 115;
    cache_priority_list_[ 83].irrep = 0;
    cache_priority_list_[ 83].pqnum =  0;
    cache_priority_list_[ 83].rsnum =  5;
    cache_priority_list_[ 83].priority = 1367;
    cache_priority_list_[ 83].next = &(cache_priority_list_[ 83+1]);
    cache_priority_list_[ 83].last = &(cache_priority_list_[ 83-1]);

    strcpy(cache_priority_list_[ 84].label, "T2 (IJ,BA)");
    cache_priority_list_[ 84].filenum = 115;
    cache_priority_list_[ 84].irrep = 0;
    cache_priority_list_[ 84].pqnum =  0;
    cache_priority_list_[ 84].rsnum =  5;
    cache_priority_list_[ 84].priority = 1367;
    cache_priority_list_[ 84].next = &(cache_priority_list_[ 84+1]);
    cache_priority_list_[ 84].last = &(cache_priority_list_[ 84-1]);

    strcpy(cache_priority_list_[ 85].label, "T2 (JI,BA)");
    cache_priority_list_[ 85].filenum = 115;
    cache_priority_list_[ 85].irrep = 0;
    cache_priority_list_[ 85].pqnum =  0;
    cache_priority_list_[ 85].rsnum =  5;
    cache_priority_list_[ 85].priority = 1367;
    cache_priority_list_[ 85].next = &(cache_priority_list_[ 85+1]);
    cache_priority_list_[ 85].last = &(cache_priority_list_[ 85-1]);

    strcpy(cache_priority_list_[ 86].label, "T2 (ia,jb)");
    cache_priority_list_[ 86].filenum = 115;
    cache_priority_list_[ 86].irrep = 0;
    cache_priority_list_[ 86].pqnum = 10;
    cache_priority_list_[ 86].rsnum = 10;
    cache_priority_list_[ 86].priority = 683;
    cache_priority_list_[ 86].next = &(cache_priority_list_[ 86+1]);
    cache_priority_list_[ 86].last = &(cache_priority_list_[ 86-1]);

    strcpy(cache_priority_list_[ 87].label, "T2 (ij,ab)");
    cache_priority_list_[ 87].filenum = 115;
    cache_priority_list_[ 87].irrep = 0;
    cache_priority_list_[ 87].pqnum =  0;
    cache_priority_list_[ 87].rsnum =  5;
    cache_priority_list_[ 87].priority = 1367;
    cache_priority_list_[ 87].next = &(cache_priority_list_[ 87+1]);
    cache_priority_list_[ 87].last = &(cache_priority_list_[ 87-1]);

    strcpy(cache_priority_list_[ 88].label, "T2 (ji,ab)");
    cache_priority_list_[ 88].filenum = 115;
    cache_priority_list_[ 88].irrep = 0;
    cache_priority_list_[ 88].pqnum =  0;
    cache_priority_list_[ 88].rsnum =  5;
    cache_priority_list_[ 88].priority = 455;
    cache_priority_list_[ 88].next = &(cache_priority_list_[ 88+1]);
    cache_priority_list_[ 88].last = &(cache_priority_list_[ 88-1]);

    strcpy(cache_priority_list_[ 89].label, "T2 (ij,ba)");
    cache_priority_list_[ 89].filenum = 115;
    cache_priority_list_[ 89].irrep = 0;
    cache_priority_list_[ 89].pqnum =  0;
    cache_priority_list_[ 89].rsnum =  5;
    cache_priority_list_[ 89].priority = 455;
    cache_priority_list_[ 89].next = &(cache_priority_list_[ 89+1]);
    cache_priority_list_[ 89].last = &(cache_priority_list_[ 89-1]);

    strcpy(cache_priority_list_[ 90].label, "T2 (ji,ba)");
    cache_priority_list_[ 90].filenum = 115;
    cache_priority_list_[ 90].irrep = 0;
    cache_priority_list_[ 90].pqnum =  0;
    cache_priority_list_[ 90].rsnum =  5;
    cache_priority_list_[ 90].priority = 455;
    cache_priority_list_[ 90].next = &(cache_priority_list_[ 90+1]);
    cache_priority_list_[ 90].last = &(cache_priority_list_[ 90-1]);

    strcpy(cache_priority_list_[ 91].label, "T2 (IA,jb)");
    cache_priority_list_[ 91].filenum = 115;
    cache_priority_list_[ 91].irrep = 0;
    cache_priority_list_[ 91].pqnum = 10;
    cache_priority_list_[ 91].rsnum = 10;
    cache_priority_list_[ 91].priority = 1291;
    cache_priority_list_[ 91].next = &(cache_priority_list_[ 91+1]);
    cache_priority_list_[ 91].last = &(cache_priority_list_[ 91-1]);

    strcpy(cache_priority_list_[ 92].label, "T2 (Ij,Ab) 1");
    cache_priority_list_[ 92].filenum = 115;
    cache_priority_list_[ 92].irrep = 0;
    cache_priority_list_[ 92].pqnum =  0;
    cache_priority_list_[ 92].rsnum =  5;
    cache_priority_list_[ 92].priority = 455;
    cache_priority_list_[ 92].next = &(cache_priority_list_[ 92+1]);
    cache_priority_list_[ 92].last = &(cache_priority_list_[ 92-1]);

    strcpy(cache_priority_list_[ 93].label, "T2 (Ib,jA)");
    cache_priority_list_[ 93].filenum = 115;
    cache_priority_list_[ 93].irrep = 0;
    cache_priority_list_[ 93].pqnum = 10;
    cache_priority_list_[ 93].rsnum = 10;
    cache_priority_list_[ 93].priority = 683;
    cache_priority_list_[ 93].next = &(cache_priority_list_[ 93+1]);
    cache_priority_list_[ 93].last = &(cache_priority_list_[ 93-1]);

    strcpy(cache_priority_list_[ 94].label, "T2 (Ij,Ab) 2");
    cache_priority_list_[ 94].filenum = 115;
    cache_priority_list_[ 94].irrep = 0;
    cache_priority_list_[ 94].pqnum =  0;
    cache_priority_list_[ 94].rsnum =  5;
    cache_priority_list_[ 94].priority = 455;
    cache_priority_list_[ 94].next = &(cache_priority_list_[ 94+1]);
    cache_priority_list_[ 94].last = &(cache_priority_list_[ 94-1]);

    strcpy(cache_priority_list_[ 95].label, "Y (MB,JI)");
    cache_priority_list_[ 95].filenum = 115;
    cache_priority_list_[ 95].irrep = 0;
    cache_priority_list_[ 95].pqnum = 10;
    cache_priority_list_[ 95].rsnum =  0;
    cache_priority_list_[ 95].priority = 1367;
    cache_priority_list_[ 95].next = &(cache_priority_list_[ 95+1]);
    cache_priority_list_[ 95].last = &(cache_priority_list_[ 95-1]);

    strcpy(cache_priority_list_[ 96].label, "T2 (AB,JI)");
    cache_priority_list_[ 96].filenum = 115;
    cache_priority_list_[ 96].irrep = 0;
    cache_priority_list_[ 96].pqnum =  5;
    cache_priority_list_[ 96].rsnum =  0;
    cache_priority_list_[ 96].priority = 1367;
    cache_priority_list_[ 96].next = &(cache_priority_list_[ 96+1]);
    cache_priority_list_[ 96].last = &(cache_priority_list_[ 96-1]);

    strcpy(cache_priority_list_[ 97].label, "Y (mA,jI)");
    cache_priority_list_[ 97].filenum = 115;
    cache_priority_list_[ 97].irrep = 0;
    cache_priority_list_[ 97].pqnum = 10;
    cache_priority_list_[ 97].rsnum =  0;
    cache_priority_list_[ 97].priority = 683;
    cache_priority_list_[ 97].next = &(cache_priority_list_[ 97+1]);
    cache_priority_list_[ 97].last = &(cache_priority_list_[ 97-1]);

    strcpy(cache_priority_list_[ 98].label, "T2 (bA,jI)");
    cache_priority_list_[ 98].filenum = 115;
    cache_priority_list_[ 98].irrep = 0;
    cache_priority_list_[ 98].pqnum =  5;
    cache_priority_list_[ 98].rsnum =  0;
    cache_priority_list_[ 98].priority = 683;
    cache_priority_list_[ 98].next = &(cache_priority_list_[ 98+1]);
    cache_priority_list_[ 98].last = &(cache_priority_list_[ 98-1]);

    strcpy(cache_priority_list_[ 99].label, "T2 (Ij,Ab)");
    cache_priority_list_[ 99].filenum = 115;
    cache_priority_list_[ 99].irrep = 0;
    cache_priority_list_[ 99].pqnum =  0;
    cache_priority_list_[ 99].rsnum =  5;
    cache_priority_list_[ 99].priority = 1823;
    cache_priority_list_[ 99].next = &(cache_priority_list_[ 99+1]);
    cache_priority_list_[ 99].last = &(cache_priority_list_[ 99-1]);

    strcpy(cache_priority_list_[100].label, "Y (Mb,Ij)");
    cache_priority_list_[100].filenum = 115;
    cache_priority_list_[100].irrep = 0;
    cache_priority_list_[100].pqnum = 10;
    cache_priority_list_[100].rsnum =  0;
    cache_priority_list_[100].priority = 683;
    cache_priority_list_[100].next = &(cache_priority_list_[100+1]);
    cache_priority_list_[100].last = &(cache_priority_list_[100-1]);

    strcpy(cache_priority_list_[101].label, "T2 (Ab,Ij)");
    cache_priority_list_[101].filenum = 115;
    cache_priority_list_[101].irrep = 0;
    cache_priority_list_[101].pqnum =  5;
    cache_priority_list_[101].rsnum =  0;
    cache_priority_list_[101].priority = 683;
    cache_priority_list_[101].next = &(cache_priority_list_[101+1]);
    cache_priority_list_[101].last = &(cache_priority_list_[101-1]);

    strcpy(cache_priority_list_[102].label, "Y(Mb,jI)");
    cache_priority_list_[102].filenum = 115;
    cache_priority_list_[102].irrep = 0;
    cache_priority_list_[102].pqnum = 10;
    cache_priority_list_[102].rsnum =  0;
    cache_priority_list_[102].priority = 683;
    cache_priority_list_[102].next = &(cache_priority_list_[102+1]);
    cache_priority_list_[102].last = &(cache_priority_list_[102-1]);

    strcpy(cache_priority_list_[103].label, "T2 (Ab,jI)");
    cache_priority_list_[103].filenum = 115;
    cache_priority_list_[103].irrep = 0;
    cache_priority_list_[103].pqnum =  5;
    cache_priority_list_[103].rsnum =  0;
    cache_priority_list_[103].priority = 683;
    cache_priority_list_[103].next = &(cache_priority_list_[103+1]);
    cache_priority_list_[103].last = &(cache_priority_list_[103-1]);

    strcpy(cache_priority_list_[104].label, "Y(mA,Ij)");
    cache_priority_list_[104].filenum = 115;
    cache_priority_list_[104].irrep = 0;
    cache_priority_list_[104].pqnum = 10;
    cache_priority_list_[104].rsnum =  0;
    cache_priority_list_[104].priority = 683;
    cache_priority_list_[104].next = &(cache_priority_list_[104+1]);
    cache_priority_list_[104].last = &(cache_priority_list_[104-1]);

    strcpy(cache_priority_list_[105].label, "T2 (bA,Ij)");
    cache_priority_list_[105].filenum = 115;
    cache_priority_list_[105].irrep = 0;
    cache_priority_list_[105].pqnum =  5;
    cache_priority_list_[105].rsnum =  0;
    cache_priority_list_[105].priority = 683;
    cache_priority_list_[105].next = &(cache_priority_list_[105+1]);
    cache_priority_list_[105].last = &(cache_priority_list_[105-1]);

    strcpy(cache_priority_list_[106].label, "T2(IJ,AB) DIIS");
    cache_priority_list_[106].filenum = 115;
    cache_priority_list_[106].irrep = 0;
    cache_priority_list_[106].pqnum =  2;
    cache_priority_list_[106].rsnum =  7;
    cache_priority_list_[106].priority = 73;
    cache_priority_list_[106].next = &(cache_priority_list_[106+1]);
    cache_priority_list_[106].last = &(cache_priority_list_[106-1]);

    strcpy(cache_priority_list_[107].label, "T2(Ij,Ab) DIIS");
    cache_priority_list_[107].filenum = 115;
    cache_priority_list_[107].irrep = 0;
    cache_priority_list_[107].pqnum =  0;
    cache_priority_list_[107].rsnum =  5;
    cache_priority_list_[107].priority = 73;
    cache_priority_list_[107].next = &(cache_priority_list_[107+1]);
    cache_priority_list_[107].last = &(cache_priority_list_[107-1]);

    /* Hand-coded entries follow */

    strcpy(cache_priority_list_[108].label, "X(5,0)");
    cache_priority_list_[108].filenum = 115;
    cache_priority_list_[108].irrep = 0;
    cache_priority_list_[108].pqnum =  5;
    cache_priority_list_[108].rsnum =  0;
    cache_priority_list_[108].priority = 99999;
    cache_priority_list_[108].next = &(cache_priority_list_[109]);
    cache_priority_list_[108].last = &(cache_priority_list_[107]);

    strcpy(cache_priority_list_[109].label, "X(0,5) 1");
    cache_priority_list_[109].filenum = 115;
    cache_priority_list_[109].irrep = 0;
    cache_priority_list_[109].pqnum =  0;
    cache_priority_list_[109].rsnum =  5;
    cache_priority_list_[109].priority = 99999;
    cache_priority_list_[109].next = &(cache_priority_list_[110]);
    cache_priority_list_[109].last = &(cache_priority_list_[108]);

    strcpy(cache_priority_list_[110].label, "X(0,5) 2");
    cache_priority_list_[110].filenum = 115;
    cache_priority_list_[110].irrep = 0;
    cache_priority_list_[110].pqnum =  0;
    cache_priority_list_[110].rsnum =  5;
    cache_priority_list_[110].priority = 99999;
    cache_priority_list_[110].next = &(cache_priority_list_[111]);
    cache_priority_list_[110].last = &(cache_priority_list_[109]);

    strcpy(cache_priority_list_[111].label, "X(0,5) 3");
    cache_priority_list_[111].filenum = 115;
    cache_priority_list_[111].irrep = 0;
    cache_priority_list_[111].pqnum =  0;
    cache_priority_list_[111].rsnum =  5;
    cache_priority_list_[111].priority = 99999;
    cache_priority_list_[111].next = &(cache_priority_list_[112]);
    cache_priority_list_[111].last = &(cache_priority_list_[110]);

    strcpy(cache_priority_list_[112].label, "X(0,5) 4");
    cache_priority_list_[112].filenum = 115;
    cache_priority_list_[112].irrep = 0;
    cache_priority_list_[112].pqnum =  0;
    cache_priority_list_[112].rsnum =  5;
    cache_priority_list_[112].priority = 99999;
    cache_priority_list_[112].next = NULL;
    cache_priority_list_[112].last = &(cache_priority_list_[111]);
}
}} // namespace psi::ccenergy
