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
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include <libdpd/dpd.h>
#include "ccwave.h"
namespace psi { namespace ccenergy {



dpd_file4_cache_entry *CCEnergyWavefunction::priority_list(void)
{

  strcpy(list_[  0].label, "D <ij||ab> (i>j,a>b)");
  list_[  0].filenum = 105;
  list_[  0].irrep = 0;
  list_[  0].pqnum =  2;
  list_[  0].rsnum =  7;
  list_[  0].priority = 231;
  list_[  0].next = &(list_[  0+1]);
  list_[  0].last = NULL;

  strcpy(list_[  1].label, "tIJAB");
  list_[  1].filenum = 109;
  list_[  1].irrep = 0;
  list_[  1].pqnum =  2;
  list_[  1].rsnum =  7;
  list_[  1].priority = 1163;
  list_[  1].next = &(list_[  1+1]);
  list_[  1].last = &(list_[  1-1]);

  strcpy(list_[  2].label, "tijab");
  list_[  2].filenum = 109;
  list_[  2].irrep = 0;
  list_[  2].pqnum =  2;
  list_[  2].rsnum =  7;
  list_[  2].priority = 1163;
  list_[  2].next = &(list_[  2+1]);
  list_[  2].last = &(list_[  2-1]);

  strcpy(list_[  3].label, "dIJAB");
  list_[  3].filenum = 108;
  list_[  3].irrep = 0;
  list_[  3].pqnum =  1;
  list_[  3].rsnum =  6;
  list_[  3].priority = 77;
  list_[  3].next = &(list_[  3+1]);
  list_[  3].last = &(list_[  3-1]);

  strcpy(list_[  4].label, "dijab");
  list_[  4].filenum = 108;
  list_[  4].irrep = 0;
  list_[  4].pqnum =  1;
  list_[  4].rsnum =  6;
  list_[  4].priority = 77;
  list_[  4].next = &(list_[  4+1]);
  list_[  4].last = &(list_[  4-1]);

  strcpy(list_[  5].label, "D <ij|ab>");
  list_[  5].filenum = 105;
  list_[  5].irrep = 0;
  list_[  5].pqnum =  0;
  list_[  5].rsnum =  5;
  list_[  5].priority = 535;
  list_[  5].next = &(list_[  5+1]);
  list_[  5].last = &(list_[  5-1]);

  strcpy(list_[  6].label, "tIjAb");
  list_[  6].filenum = 109;
  list_[  6].irrep = 0;
  list_[  6].pqnum =  0;
  list_[  6].rsnum =  5;
  list_[  6].priority = 1241;
  list_[  6].next = &(list_[  6+1]);
  list_[  6].last = &(list_[  6-1]);

  strcpy(list_[  7].label, "dIjAb");
  list_[  7].filenum = 108;
  list_[  7].irrep = 0;
  list_[  7].pqnum =  0;
  list_[  7].rsnum =  5;
  list_[  7].priority = 77;
  list_[  7].next = &(list_[  7+1]);
  list_[  7].last = &(list_[  7-1]);

  strcpy(list_[  8].label, "tauIJAB");
  list_[  8].filenum = 109;
  list_[  8].irrep = 0;
  list_[  8].pqnum =  2;
  list_[  8].rsnum =  7;
  list_[  8].priority = 1161;
  list_[  8].next = &(list_[  8+1]);
  list_[  8].last = &(list_[  8-1]);

  strcpy(list_[  9].label, "tauijab");
  list_[  9].filenum = 109;
  list_[  9].irrep = 0;
  list_[  9].pqnum =  2;
  list_[  9].rsnum =  7;
  list_[  9].priority = 1161;
  list_[  9].next = &(list_[  9+1]);
  list_[  9].last = &(list_[  9-1]);

  strcpy(list_[ 10].label, "tauIjAb");
  list_[ 10].filenum = 109;
  list_[ 10].irrep = 0;
  list_[ 10].pqnum =  0;
  list_[ 10].rsnum =  5;
  list_[ 10].priority = 1161;
  list_[ 10].next = &(list_[ 10+1]);
  list_[ 10].last = &(list_[ 10-1]);

  strcpy(list_[ 11].label, "tauIjbA");
  list_[ 11].filenum = 109;
  list_[ 11].irrep = 0;
  list_[ 11].pqnum =  0;
  list_[ 11].rsnum =  5;
  list_[ 11].priority = 543;
  list_[ 11].next = &(list_[ 11+1]);
  list_[ 11].last = &(list_[ 11-1]);

  strcpy(list_[ 12].label, "tauiJaB");
  list_[ 12].filenum = 109;
  list_[ 12].irrep = 0;
  list_[ 12].pqnum =  0;
  list_[ 12].rsnum =  5;
  list_[ 12].priority = 389;
  list_[ 12].next = &(list_[ 12+1]);
  list_[ 12].last = &(list_[ 12-1]);

  strcpy(list_[ 13].label, "tautIJAB");
  list_[ 13].filenum = 109;
  list_[ 13].irrep = 0;
  list_[ 13].pqnum =  2;
  list_[ 13].rsnum =  7;
  list_[ 13].priority = 931;
  list_[ 13].next = &(list_[ 13+1]);
  list_[ 13].last = &(list_[ 13-1]);

  strcpy(list_[ 14].label, "tautijab");
  list_[ 14].filenum = 109;
  list_[ 14].irrep = 0;
  list_[ 14].pqnum =  2;
  list_[ 14].rsnum =  7;
  list_[ 14].priority = 931;
  list_[ 14].next = &(list_[ 14+1]);
  list_[ 14].last = &(list_[ 14-1]);

  strcpy(list_[ 15].label, "tautIjAb");
  list_[ 15].filenum = 109;
  list_[ 15].irrep = 0;
  list_[ 15].pqnum =  0;
  list_[ 15].rsnum =  5;
  list_[ 15].priority = 931;
  list_[ 15].next = &(list_[ 15+1]);
  list_[ 15].last = &(list_[ 15-1]);

  strcpy(list_[ 16].label, "tiJaB");
  list_[ 16].filenum = 109;
  list_[ 16].irrep = 0;
  list_[ 16].pqnum =  0;
  list_[ 16].rsnum =  5;
  list_[ 16].priority = 541;
  list_[ 16].next = &(list_[ 16+1]);
  list_[ 16].last = &(list_[ 16-1]);

  strcpy(list_[ 17].label, "tIAJB");
  list_[ 17].filenum = 109;
  list_[ 17].irrep = 0;
  list_[ 17].pqnum = 10;
  list_[ 17].rsnum = 10;
  list_[ 17].priority = 693;
  list_[ 17].next = &(list_[ 17+1]);
  list_[ 17].last = &(list_[ 17-1]);

  strcpy(list_[ 18].label, "tiajb");
  list_[ 18].filenum = 109;
  list_[ 18].irrep = 0;
  list_[ 18].pqnum = 10;
  list_[ 18].rsnum = 10;
  list_[ 18].priority = 693;
  list_[ 18].next = &(list_[ 18+1]);
  list_[ 18].last = &(list_[ 18-1]);

  strcpy(list_[ 19].label, "tIAjb");
  list_[ 19].filenum = 109;
  list_[ 19].irrep = 0;
  list_[ 19].pqnum = 10;
  list_[ 19].rsnum = 10;
  list_[ 19].priority = 849;
  list_[ 19].next = &(list_[ 19+1]);
  list_[ 19].last = &(list_[ 19-1]);

  strcpy(list_[ 20].label, "tiaJB");
  list_[ 20].filenum = 109;
  list_[ 20].irrep = 0;
  list_[ 20].pqnum = 10;
  list_[ 20].rsnum = 10;
  list_[ 20].priority = 693;
  list_[ 20].next = &(list_[ 20+1]);
  list_[ 20].last = &(list_[ 20-1]);

  strcpy(list_[ 21].label, "tIbjA");
  list_[ 21].filenum = 109;
  list_[ 21].irrep = 0;
  list_[ 21].pqnum = 10;
  list_[ 21].rsnum = 10;
  list_[ 21].priority = 619;
  list_[ 21].next = &(list_[ 21+1]);
  list_[ 21].last = &(list_[ 21-1]);

  strcpy(list_[ 22].label, "tjAIb");
  list_[ 22].filenum = 109;
  list_[ 22].irrep = 0;
  list_[ 22].pqnum = 10;
  list_[ 22].rsnum = 10;
  list_[ 22].priority = 541;
  list_[ 22].next = &(list_[ 22+1]);
  list_[ 22].last = &(list_[ 22-1]);

  strcpy(list_[ 23].label, "C <ia||jb> (ia,bj)");
  list_[ 23].filenum = 104;
  list_[ 23].irrep = 0;
  list_[ 23].pqnum = 10;
  list_[ 23].rsnum = 11;
  list_[ 23].priority = 75;
  list_[ 23].next = &(list_[ 23+1]);
  list_[ 23].last = &(list_[ 23-1]);

  strcpy(list_[ 24].label, "WMBEJ");
  list_[ 24].filenum = 115;
  list_[ 24].irrep = 0;
  list_[ 24].pqnum = 10;
  list_[ 24].rsnum = 11;
  list_[ 24].priority = 1823;
  list_[ 24].next = &(list_[ 24+1]);
  list_[ 24].last = &(list_[ 24-1]);

  strcpy(list_[ 25].label, "Wmbej");
  list_[ 25].filenum = 115;
  list_[ 25].irrep = 0;
  list_[ 25].pqnum = 10;
  list_[ 25].rsnum = 11;
  list_[ 25].priority = 1823;
  list_[ 25].next = &(list_[ 25+1]);
  list_[ 25].last = &(list_[ 25-1]);

  strcpy(list_[ 26].label, "C <ia|jb>");
  list_[ 26].filenum = 104;
  list_[ 26].irrep = 0;
  list_[ 26].pqnum = 10;
  list_[ 26].rsnum = 10;
  list_[ 26].priority = 227;
  list_[ 26].next = &(list_[ 26+1]);
  list_[ 26].last = &(list_[ 26-1]);

  strcpy(list_[ 27].label, "WmBEj");
  list_[ 27].filenum = 115;
  list_[ 27].irrep = 0;
  list_[ 27].pqnum = 10;
  list_[ 27].rsnum = 10;
  list_[ 27].priority = 1823;
  list_[ 27].next = &(list_[ 27+1]);
  list_[ 27].last = &(list_[ 27-1]);

  strcpy(list_[ 28].label, "WMbeJ");
  list_[ 28].filenum = 115;
  list_[ 28].irrep = 0;
  list_[ 28].pqnum = 10;
  list_[ 28].rsnum = 10;
  list_[ 28].priority = 1823;
  list_[ 28].next = &(list_[ 28+1]);
  list_[ 28].last = &(list_[ 28-1]);

  strcpy(list_[ 29].label, "D <ij|ab> (ib,aj)");
  list_[ 29].filenum = 105;
  list_[ 29].irrep = 0;
  list_[ 29].pqnum = 10;
  list_[ 29].rsnum = 11;
  list_[ 29].priority = 227;
  list_[ 29].next = &(list_[ 29+1]);
  list_[ 29].last = &(list_[ 29-1]);

  strcpy(list_[ 30].label, "WMbEj");
  list_[ 30].filenum = 115;
  list_[ 30].irrep = 0;
  list_[ 30].pqnum = 10;
  list_[ 30].rsnum = 11;
  list_[ 30].priority = 1823;
  list_[ 30].next = &(list_[ 30+1]);
  list_[ 30].last = &(list_[ 30-1]);

  strcpy(list_[ 31].label, "WmBeJ");
  list_[ 31].filenum = 115;
  list_[ 31].irrep = 0;
  list_[ 31].pqnum = 10;
  list_[ 31].rsnum = 11;
  list_[ 31].priority = 1823;
  list_[ 31].next = &(list_[ 31+1]);
  list_[ 31].last = &(list_[ 31-1]);

  strcpy(list_[ 32].label, "F <ia||bc> (ia,b>c)");
  list_[ 32].filenum = 107;
  list_[ 32].irrep = 0;
  list_[ 32].pqnum = 10;
  list_[ 32].rsnum =  7;
  list_[ 32].priority = 455;
  list_[ 32].next = &(list_[ 32+1]);
  list_[ 32].last = &(list_[ 32-1]);

  strcpy(list_[ 33].label, "F <ia|bc>");
  list_[ 33].filenum = 107;
  list_[ 33].irrep = 0;
  list_[ 33].pqnum = 10;
  list_[ 33].rsnum =  5;
  list_[ 33].priority = 379;
  list_[ 33].next = &(list_[ 33+1]);
  list_[ 33].last = &(list_[ 33-1]);

  strcpy(list_[ 34].label, "E <ij||ka> (i>j,ak)");
  list_[ 34].filenum = 106;
  list_[ 34].irrep = 0;
  list_[ 34].pqnum =  2;
  list_[ 34].rsnum = 11;
  list_[ 34].priority = 75;
  list_[ 34].next = &(list_[ 34+1]);
  list_[ 34].last = &(list_[ 34-1]);

  strcpy(list_[ 35].label, "E <ai|jk>");
  list_[ 35].filenum = 106;
  list_[ 35].irrep = 0;
  list_[ 35].pqnum = 11;
  list_[ 35].rsnum =  0;
  list_[ 35].priority = 759;
  list_[ 35].next = &(list_[ 35+1]);
  list_[ 35].last = &(list_[ 35-1]);

  strcpy(list_[ 36].label, "E <ij|ka>");
  list_[ 36].filenum = 106;
  list_[ 36].irrep = 0;
  list_[ 36].pqnum =  0;
  list_[ 36].rsnum = 10;
  list_[ 36].priority = 151;
  list_[ 36].next = &(list_[ 36+1]);
  list_[ 36].last = &(list_[ 36-1]);

  strcpy(list_[ 37].label, "WMBEJ");
  list_[ 37].filenum = 111;
  list_[ 37].irrep = 0;
  list_[ 37].pqnum = 10;
  list_[ 37].rsnum = 10;
  list_[ 37].priority = 1899;
  list_[ 37].next = &(list_[ 37+1]);
  list_[ 37].last = &(list_[ 37-1]);

  strcpy(list_[ 38].label, "Wmbej");
  list_[ 38].filenum = 111;
  list_[ 38].irrep = 0;
  list_[ 38].pqnum = 10;
  list_[ 38].rsnum = 10;
  list_[ 38].priority = 1899;
  list_[ 38].next = &(list_[ 38+1]);
  list_[ 38].last = &(list_[ 38-1]);

  strcpy(list_[ 39].label, "WMbEj");
  list_[ 39].filenum = 111;
  list_[ 39].irrep = 0;
  list_[ 39].pqnum = 10;
  list_[ 39].rsnum = 10;
  list_[ 39].priority = 1899;
  list_[ 39].next = &(list_[ 39+1]);
  list_[ 39].last = &(list_[ 39-1]);

  strcpy(list_[ 40].label, "WmBeJ");
  list_[ 40].filenum = 111;
  list_[ 40].irrep = 0;
  list_[ 40].pqnum = 10;
  list_[ 40].rsnum = 10;
  list_[ 40].priority = 1899;
  list_[ 40].next = &(list_[ 40+1]);
  list_[ 40].last = &(list_[ 40-1]);

  strcpy(list_[ 41].label, "WMbeJ");
  list_[ 41].filenum = 111;
  list_[ 41].irrep = 0;
  list_[ 41].pqnum = 10;
  list_[ 41].rsnum = 10;
  list_[ 41].priority = 1519;
  list_[ 41].next = &(list_[ 41+1]);
  list_[ 41].last = &(list_[ 41-1]);

  strcpy(list_[ 42].label, "WmBEj");
  list_[ 42].filenum = 111;
  list_[ 42].irrep = 0;
  list_[ 42].pqnum = 10;
  list_[ 42].rsnum = 10;
  list_[ 42].priority = 1519;
  list_[ 42].next = &(list_[ 42+1]);
  list_[ 42].last = &(list_[ 42-1]);

  strcpy(list_[ 43].label, "D <ij||ab> (ia,jb)");
  list_[ 43].filenum = 105;
  list_[ 43].irrep = 0;
  list_[ 43].pqnum = 10;
  list_[ 43].rsnum = 10;
  list_[ 43].priority = 303;
  list_[ 43].next = &(list_[ 43+1]);
  list_[ 43].last = &(list_[ 43-1]);

  strcpy(list_[ 44].label, "D <ij|ab> (ia,jb)");
  list_[ 44].filenum = 105;
  list_[ 44].irrep = 0;
  list_[ 44].pqnum = 10;
  list_[ 44].rsnum = 10;
  list_[ 44].priority = 303;
  list_[ 44].next = &(list_[ 44+1]);
  list_[ 44].last = &(list_[ 44-1]);

  strcpy(list_[ 45].label, "Y (ME,JN)");
  list_[ 45].filenum = 115;
  list_[ 45].irrep = 0;
  list_[ 45].pqnum = 10;
  list_[ 45].rsnum =  0;
  list_[ 45].priority = 4103;
  list_[ 45].next = &(list_[ 45+1]);
  list_[ 45].last = &(list_[ 45-1]);

  strcpy(list_[ 46].label, "D <ij||ab> (ia,bj)");
  list_[ 46].filenum = 105;
  list_[ 46].irrep = 0;
  list_[ 46].pqnum = 10;
  list_[ 46].rsnum = 11;
  list_[ 46].priority = 151;
  list_[ 46].next = &(list_[ 46+1]);
  list_[ 46].last = &(list_[ 46-1]);

  strcpy(list_[ 47].label, "D <ij|ab> (ia,bj)");
  list_[ 47].filenum = 105;
  list_[ 47].irrep = 0;
  list_[ 47].pqnum = 10;
  list_[ 47].rsnum = 11;
  list_[ 47].priority = 151;
  list_[ 47].next = &(list_[ 47+1]);
  list_[ 47].last = &(list_[ 47-1]);

  strcpy(list_[ 48].label, "D <ij|ab> (ib,ja)");
  list_[ 48].filenum = 105;
  list_[ 48].irrep = 0;
  list_[ 48].pqnum = 10;
  list_[ 48].rsnum = 10;
  list_[ 48].priority = 303;
  list_[ 48].next = &(list_[ 48+1]);
  list_[ 48].last = &(list_[ 48-1]);

  strcpy(list_[ 49].label, "D <ij||ab>");
  list_[ 49].filenum = 105;
  list_[ 49].irrep = 0;
  list_[ 49].pqnum =  0;
  list_[ 49].rsnum =  5;
  list_[ 49].priority = 75;
  list_[ 49].next = &(list_[ 49+1]);
  list_[ 49].last = &(list_[ 49-1]);

  strcpy(list_[ 50].label, "D <ij||ab> (i>j,ab)");
  list_[ 50].filenum = 105;
  list_[ 50].irrep = 0;
  list_[ 50].pqnum =  2;
  list_[ 50].rsnum =  5;
  list_[ 50].priority = 75;
  list_[ 50].next = &(list_[ 50+1]);
  list_[ 50].last = &(list_[ 50-1]);

  strcpy(list_[ 51].label, "D <ij||ab> (ij,a>b)");
  list_[ 51].filenum = 105;
  list_[ 51].irrep = 0;
  list_[ 51].pqnum =  0;
  list_[ 51].rsnum =  7;
  list_[ 51].priority = 75;
  list_[ 51].next = &(list_[ 51+1]);
  list_[ 51].last = &(list_[ 51-1]);

  strcpy(list_[ 52].label, "C <ia||jb>");
  list_[ 52].filenum = 104;
  list_[ 52].irrep = 0;
  list_[ 52].pqnum = 10;
  list_[ 52].rsnum = 10;
  list_[ 52].priority = 227;
  list_[ 52].next = &(list_[ 52+1]);
  list_[ 52].last = &(list_[ 52-1]);

  strcpy(list_[ 53].label, "A <ij|kl>");
  list_[ 53].filenum = 102;
  list_[ 53].irrep = 0;
  list_[ 53].pqnum =  0;
  list_[ 53].rsnum =  0;
  list_[ 53].priority = 151;
  list_[ 53].next = &(list_[ 53+1]);
  list_[ 53].last = &(list_[ 53-1]);

  strcpy(list_[ 54].label, "WMNIJ");
  list_[ 54].filenum = 111;
  list_[ 54].irrep = 0;
  list_[ 54].pqnum =  2;
  list_[ 54].rsnum =  2;
  list_[ 54].priority = 8891;
  list_[ 54].next = &(list_[ 54+1]);
  list_[ 54].last = &(list_[ 54-1]);

  strcpy(list_[ 55].label, "Wmnij");
  list_[ 55].filenum = 111;
  list_[ 55].irrep = 0;
  list_[ 55].pqnum =  2;
  list_[ 55].rsnum =  2;
  list_[ 55].priority = 8891;
  list_[ 55].next = &(list_[ 55+1]);
  list_[ 55].last = &(list_[ 55-1]);

  strcpy(list_[ 56].label, "WMnIj");
  list_[ 56].filenum = 111;
  list_[ 56].irrep = 0;
  list_[ 56].pqnum =  0;
  list_[ 56].rsnum =  0;
  list_[ 56].priority = 2127;
  list_[ 56].next = &(list_[ 56+1]);
  list_[ 56].last = &(list_[ 56-1]);

  strcpy(list_[ 57].label, "E <ij||ka> (i>j,ka)");
  list_[ 57].filenum = 106;
  list_[ 57].irrep = 0;
  list_[ 57].pqnum =  2;
  list_[ 57].rsnum = 10;
  list_[ 57].priority = 75;
  list_[ 57].next = &(list_[ 57+1]);
  list_[ 57].last = &(list_[ 57-1]);

  strcpy(list_[ 58].label, "W (MN,IJ)");
  list_[ 58].filenum = 115;
  list_[ 58].irrep = 0;
  list_[ 58].pqnum =  2;
  list_[ 58].rsnum =  0;
  list_[ 58].priority = 2583;
  list_[ 58].next = &(list_[ 58+1]);
  list_[ 58].last = &(list_[ 58-1]);

  strcpy(list_[ 59].label, "ZIJMA");
  list_[ 59].filenum = 113;
  list_[ 59].irrep = 0;
  list_[ 59].pqnum =  2;
  list_[ 59].rsnum = 10;
  list_[ 59].priority = 455;
  list_[ 59].next = &(list_[ 59+1]);
  list_[ 59].last = &(list_[ 59-1]);

  strcpy(list_[ 60].label, "Zijma");
  list_[ 60].filenum = 113;
  list_[ 60].irrep = 0;
  list_[ 60].pqnum =  2;
  list_[ 60].rsnum = 10;
  list_[ 60].priority = 455;
  list_[ 60].next = &(list_[ 60+1]);
  list_[ 60].last = &(list_[ 60-1]);

  strcpy(list_[ 61].label, "ZIjMa");
  list_[ 61].filenum = 113;
  list_[ 61].irrep = 0;
  list_[ 61].pqnum =  0;
  list_[ 61].rsnum = 10;
  list_[ 61].priority = 455;
  list_[ 61].next = &(list_[ 61+1]);
  list_[ 61].last = &(list_[ 61-1]);

  strcpy(list_[ 62].label, "ZIjmA");
  list_[ 62].filenum = 113;
  list_[ 62].irrep = 0;
  list_[ 62].pqnum =  0;
  list_[ 62].rsnum = 10;
  list_[ 62].priority = 379;
  list_[ 62].next = &(list_[ 62+1]);
  list_[ 62].last = &(list_[ 62-1]);

  strcpy(list_[ 63].label, "ZIJAM");
  list_[ 63].filenum = 113;
  list_[ 63].irrep = 0;
  list_[ 63].pqnum =  2;
  list_[ 63].rsnum = 11;
  list_[ 63].priority = 455;
  list_[ 63].next = &(list_[ 63+1]);
  list_[ 63].last = &(list_[ 63-1]);

  strcpy(list_[ 64].label, "Zijam");
  list_[ 64].filenum = 113;
  list_[ 64].irrep = 0;
  list_[ 64].pqnum =  2;
  list_[ 64].rsnum = 11;
  list_[ 64].priority = 455;
  list_[ 64].next = &(list_[ 64+1]);
  list_[ 64].last = &(list_[ 64-1]);

  strcpy(list_[ 65].label, "ZIjAm");
  list_[ 65].filenum = 113;
  list_[ 65].irrep = 0;
  list_[ 65].pqnum =  0;
  list_[ 65].rsnum = 11;
  list_[ 65].priority = 455;
  list_[ 65].next = &(list_[ 65+1]);
  list_[ 65].last = &(list_[ 65-1]);

  strcpy(list_[ 66].label, "New tIJAB");
  list_[ 66].filenum = 109;
  list_[ 66].irrep = 0;
  list_[ 66].pqnum =  2;
  list_[ 66].rsnum =  7;
  list_[ 66].priority = 58571;
  list_[ 66].next = &(list_[ 66+1]);
  list_[ 66].last = &(list_[ 66-1]);

  strcpy(list_[ 67].label, "New tijab");
  list_[ 67].filenum = 109;
  list_[ 67].irrep = 0;
  list_[ 67].pqnum =  2;
  list_[ 67].rsnum =  7;
  list_[ 67].priority = 58571;
  list_[ 67].next = &(list_[ 67+1]);
  list_[ 67].last = &(list_[ 67-1]);

  strcpy(list_[ 68].label, "New tIjAb");
  list_[ 68].filenum = 109;
  list_[ 68].irrep = 0;
  list_[ 68].pqnum =  0;
  list_[ 68].rsnum =  5;
  list_[ 68].priority = 10843;
  list_[ 68].next = &(list_[ 68+1]);
  list_[ 68].last = &(list_[ 68-1]);

  strcpy(list_[ 69].label, "T (I>J,AB)");
  list_[ 69].filenum = 115;
  list_[ 69].irrep = 0;
  list_[ 69].pqnum =  2;
  list_[ 69].rsnum =  5;
  list_[ 69].priority = 7295;
  list_[ 69].next = &(list_[ 69+1]);
  list_[ 69].last = &(list_[ 69-1]);

  strcpy(list_[ 70].label, "T (IJ,A>B)");
  list_[ 70].filenum = 115;
  list_[ 70].irrep = 0;
  list_[ 70].pqnum =  0;
  list_[ 70].rsnum =  7;
  list_[ 70].priority = 4711;
  list_[ 70].next = &(list_[ 70+1]);
  list_[ 70].last = &(list_[ 70-1]);

  strcpy(list_[ 71].label, "B <ab||cd> (a>b,c>d)");
  list_[ 71].filenum = 103;
  list_[ 71].irrep = 0;
  list_[ 71].pqnum =  7;
  list_[ 71].rsnum =  7;
  list_[ 71].priority = 75;
  list_[ 71].next = &(list_[ 71+1]);
  list_[ 71].last = &(list_[ 71-1]);

  strcpy(list_[ 72].label, "B <ab|cd>");
  list_[ 72].filenum = 103;
  list_[ 72].irrep = 0;
  list_[ 72].pqnum =  5;
  list_[ 72].rsnum =  5;
  list_[ 72].priority = 75;
  list_[ 72].next = &(list_[ 72+1]);
  list_[ 72].last = &(list_[ 72-1]);

  strcpy(list_[ 73].label, "Z(ab,ij)");
  list_[ 73].filenum = 115;
  list_[ 73].irrep = 0;
  list_[ 73].pqnum =  7;
  list_[ 73].rsnum =  2;
  list_[ 73].priority = 683;
  list_[ 73].next = &(list_[ 73+1]);
  list_[ 73].last = &(list_[ 73-1]);

  strcpy(list_[ 74].label, "Z(ij,ab)");
  list_[ 74].filenum = 115;
  list_[ 74].irrep = 0;
  list_[ 74].pqnum =  2;
  list_[ 74].rsnum =  7;
  list_[ 74].priority = 911;
  list_[ 74].next = &(list_[ 74+1]);
  list_[ 74].last = &(list_[ 74-1]);

  strcpy(list_[ 75].label, "Z(Ab,Ij)");
  list_[ 75].filenum = 115;
  list_[ 75].irrep = 0;
  list_[ 75].pqnum =  5;
  list_[ 75].rsnum =  0;
  list_[ 75].priority = 379;
  list_[ 75].next = &(list_[ 75+1]);
  list_[ 75].last = &(list_[ 75-1]);

  strcpy(list_[ 76].label, "Z(Ij,Ab)");
  list_[ 76].filenum = 115;
  list_[ 76].irrep = 0;
  list_[ 76].pqnum =  0;
  list_[ 76].rsnum =  5;
  list_[ 76].priority = 455;
  list_[ 76].next = &(list_[ 76+1]);
  list_[ 76].last = &(list_[ 76-1]);

  strcpy(list_[ 77].label, "T (JI,A>B)");
  list_[ 77].filenum = 115;
  list_[ 77].irrep = 0;
  list_[ 77].pqnum =  0;
  list_[ 77].rsnum =  7;
  list_[ 77].priority = 911;
  list_[ 77].next = &(list_[ 77+1]);
  list_[ 77].last = &(list_[ 77-1]);

  strcpy(list_[ 78].label, "F <ai|bc>");
  list_[ 78].filenum = 107;
  list_[ 78].irrep = 0;
  list_[ 78].pqnum = 11;
  list_[ 78].rsnum =  5;
  list_[ 78].priority = 75;
  list_[ 78].next = &(list_[ 78+1]);
  list_[ 78].last = &(list_[ 78-1]);

  strcpy(list_[ 79].label, "T (I>J,BA)");
  list_[ 79].filenum = 115;
  list_[ 79].irrep = 0;
  list_[ 79].pqnum =  2;
  list_[ 79].rsnum =  5;
  list_[ 79].priority = 911;
  list_[ 79].next = &(list_[ 79+1]);
  list_[ 79].last = &(list_[ 79-1]);

  strcpy(list_[ 80].label, "E <ia|jk>");
  list_[ 80].filenum = 106;
  list_[ 80].irrep = 0;
  list_[ 80].pqnum = 10;
  list_[ 80].rsnum =  0;
  list_[ 80].priority = 75;
  list_[ 80].next = &(list_[ 80+1]);
  list_[ 80].last = &(list_[ 80-1]);

  strcpy(list_[ 81].label, "T2 (IA,JB)");
  list_[ 81].filenum = 115;
  list_[ 81].irrep = 0;
  list_[ 81].pqnum = 10;
  list_[ 81].rsnum = 10;
  list_[ 81].priority = 683;
  list_[ 81].next = &(list_[ 81+1]);
  list_[ 81].last = &(list_[ 81-1]);

  strcpy(list_[ 82].label, "T2 (IJ,AB)");
  list_[ 82].filenum = 115;
  list_[ 82].irrep = 0;
  list_[ 82].pqnum =  0;
  list_[ 82].rsnum =  5;
  list_[ 82].priority = 4103;
  list_[ 82].next = &(list_[ 82+1]);
  list_[ 82].last = &(list_[ 82-1]);

  strcpy(list_[ 83].label, "T2 (JI,AB)");
  list_[ 83].filenum = 115;
  list_[ 83].irrep = 0;
  list_[ 83].pqnum =  0;
  list_[ 83].rsnum =  5;
  list_[ 83].priority = 1367;
  list_[ 83].next = &(list_[ 83+1]);
  list_[ 83].last = &(list_[ 83-1]);

  strcpy(list_[ 84].label, "T2 (IJ,BA)");
  list_[ 84].filenum = 115;
  list_[ 84].irrep = 0;
  list_[ 84].pqnum =  0;
  list_[ 84].rsnum =  5;
  list_[ 84].priority = 1367;
  list_[ 84].next = &(list_[ 84+1]);
  list_[ 84].last = &(list_[ 84-1]);

  strcpy(list_[ 85].label, "T2 (JI,BA)");
  list_[ 85].filenum = 115;
  list_[ 85].irrep = 0;
  list_[ 85].pqnum =  0;
  list_[ 85].rsnum =  5;
  list_[ 85].priority = 1367;
  list_[ 85].next = &(list_[ 85+1]);
  list_[ 85].last = &(list_[ 85-1]);

  strcpy(list_[ 86].label, "T2 (ia,jb)");
  list_[ 86].filenum = 115;
  list_[ 86].irrep = 0;
  list_[ 86].pqnum = 10;
  list_[ 86].rsnum = 10;
  list_[ 86].priority = 683;
  list_[ 86].next = &(list_[ 86+1]);
  list_[ 86].last = &(list_[ 86-1]);

  strcpy(list_[ 87].label, "T2 (ij,ab)");
  list_[ 87].filenum = 115;
  list_[ 87].irrep = 0;
  list_[ 87].pqnum =  0;
  list_[ 87].rsnum =  5;
  list_[ 87].priority = 1367;
  list_[ 87].next = &(list_[ 87+1]);
  list_[ 87].last = &(list_[ 87-1]);

  strcpy(list_[ 88].label, "T2 (ji,ab)");
  list_[ 88].filenum = 115;
  list_[ 88].irrep = 0;
  list_[ 88].pqnum =  0;
  list_[ 88].rsnum =  5;
  list_[ 88].priority = 455;
  list_[ 88].next = &(list_[ 88+1]);
  list_[ 88].last = &(list_[ 88-1]);

  strcpy(list_[ 89].label, "T2 (ij,ba)");
  list_[ 89].filenum = 115;
  list_[ 89].irrep = 0;
  list_[ 89].pqnum =  0;
  list_[ 89].rsnum =  5;
  list_[ 89].priority = 455;
  list_[ 89].next = &(list_[ 89+1]);
  list_[ 89].last = &(list_[ 89-1]);

  strcpy(list_[ 90].label, "T2 (ji,ba)");
  list_[ 90].filenum = 115;
  list_[ 90].irrep = 0;
  list_[ 90].pqnum =  0;
  list_[ 90].rsnum =  5;
  list_[ 90].priority = 455;
  list_[ 90].next = &(list_[ 90+1]);
  list_[ 90].last = &(list_[ 90-1]);

  strcpy(list_[ 91].label, "T2 (IA,jb)");
  list_[ 91].filenum = 115;
  list_[ 91].irrep = 0;
  list_[ 91].pqnum = 10;
  list_[ 91].rsnum = 10;
  list_[ 91].priority = 1291;
  list_[ 91].next = &(list_[ 91+1]);
  list_[ 91].last = &(list_[ 91-1]);

  strcpy(list_[ 92].label, "T2 (Ij,Ab) 1");
  list_[ 92].filenum = 115;
  list_[ 92].irrep = 0;
  list_[ 92].pqnum =  0;
  list_[ 92].rsnum =  5;
  list_[ 92].priority = 455;
  list_[ 92].next = &(list_[ 92+1]);
  list_[ 92].last = &(list_[ 92-1]);

  strcpy(list_[ 93].label, "T2 (Ib,jA)");
  list_[ 93].filenum = 115;
  list_[ 93].irrep = 0;
  list_[ 93].pqnum = 10;
  list_[ 93].rsnum = 10;
  list_[ 93].priority = 683;
  list_[ 93].next = &(list_[ 93+1]);
  list_[ 93].last = &(list_[ 93-1]);

  strcpy(list_[ 94].label, "T2 (Ij,Ab) 2");
  list_[ 94].filenum = 115;
  list_[ 94].irrep = 0;
  list_[ 94].pqnum =  0;
  list_[ 94].rsnum =  5;
  list_[ 94].priority = 455;
  list_[ 94].next = &(list_[ 94+1]);
  list_[ 94].last = &(list_[ 94-1]);

  strcpy(list_[ 95].label, "Y (MB,JI)");
  list_[ 95].filenum = 115;
  list_[ 95].irrep = 0;
  list_[ 95].pqnum = 10;
  list_[ 95].rsnum =  0;
  list_[ 95].priority = 1367;
  list_[ 95].next = &(list_[ 95+1]);
  list_[ 95].last = &(list_[ 95-1]);

  strcpy(list_[ 96].label, "T2 (AB,JI)");
  list_[ 96].filenum = 115;
  list_[ 96].irrep = 0;
  list_[ 96].pqnum =  5;
  list_[ 96].rsnum =  0;
  list_[ 96].priority = 1367;
  list_[ 96].next = &(list_[ 96+1]);
  list_[ 96].last = &(list_[ 96-1]);

  strcpy(list_[ 97].label, "Y (mA,jI)");
  list_[ 97].filenum = 115;
  list_[ 97].irrep = 0;
  list_[ 97].pqnum = 10;
  list_[ 97].rsnum =  0;
  list_[ 97].priority = 683;
  list_[ 97].next = &(list_[ 97+1]);
  list_[ 97].last = &(list_[ 97-1]);

  strcpy(list_[ 98].label, "T2 (bA,jI)");
  list_[ 98].filenum = 115;
  list_[ 98].irrep = 0;
  list_[ 98].pqnum =  5;
  list_[ 98].rsnum =  0;
  list_[ 98].priority = 683;
  list_[ 98].next = &(list_[ 98+1]);
  list_[ 98].last = &(list_[ 98-1]);

  strcpy(list_[ 99].label, "T2 (Ij,Ab)");
  list_[ 99].filenum = 115;
  list_[ 99].irrep = 0;
  list_[ 99].pqnum =  0;
  list_[ 99].rsnum =  5;
  list_[ 99].priority = 1823;
  list_[ 99].next = &(list_[ 99+1]);
  list_[ 99].last = &(list_[ 99-1]);

  strcpy(list_[100].label, "Y (Mb,Ij)");
  list_[100].filenum = 115;
  list_[100].irrep = 0;
  list_[100].pqnum = 10;
  list_[100].rsnum =  0;
  list_[100].priority = 683;
  list_[100].next = &(list_[100+1]);
  list_[100].last = &(list_[100-1]);

  strcpy(list_[101].label, "T2 (Ab,Ij)");
  list_[101].filenum = 115;
  list_[101].irrep = 0;
  list_[101].pqnum =  5;
  list_[101].rsnum =  0;
  list_[101].priority = 683;
  list_[101].next = &(list_[101+1]);
  list_[101].last = &(list_[101-1]);

  strcpy(list_[102].label, "Y(Mb,jI)");
  list_[102].filenum = 115;
  list_[102].irrep = 0;
  list_[102].pqnum = 10;
  list_[102].rsnum =  0;
  list_[102].priority = 683;
  list_[102].next = &(list_[102+1]);
  list_[102].last = &(list_[102-1]);

  strcpy(list_[103].label, "T2 (Ab,jI)");
  list_[103].filenum = 115;
  list_[103].irrep = 0;
  list_[103].pqnum =  5;
  list_[103].rsnum =  0;
  list_[103].priority = 683;
  list_[103].next = &(list_[103+1]);
  list_[103].last = &(list_[103-1]);

  strcpy(list_[104].label, "Y(mA,Ij)");
  list_[104].filenum = 115;
  list_[104].irrep = 0;
  list_[104].pqnum = 10;
  list_[104].rsnum =  0;
  list_[104].priority = 683;
  list_[104].next = &(list_[104+1]);
  list_[104].last = &(list_[104-1]);

  strcpy(list_[105].label, "T2 (bA,Ij)");
  list_[105].filenum = 115;
  list_[105].irrep = 0;
  list_[105].pqnum =  5;
  list_[105].rsnum =  0;
  list_[105].priority = 683;
  list_[105].next = &(list_[105+1]);
  list_[105].last = &(list_[105-1]);

  strcpy(list_[106].label, "T2(IJ,AB) DIIS");
  list_[106].filenum = 115;
  list_[106].irrep = 0;
  list_[106].pqnum =  2;
  list_[106].rsnum =  7;
  list_[106].priority = 73;
  list_[106].next = &(list_[106+1]);
  list_[106].last = &(list_[106-1]);

  strcpy(list_[107].label, "T2(Ij,Ab) DIIS");
  list_[107].filenum = 115;
  list_[107].irrep = 0;
  list_[107].pqnum =  0;
  list_[107].rsnum =  5;
  list_[107].priority = 73;
  list_[107].next = &(list_[107+1]);
  list_[107].last = &(list_[107-1]);

  /* Hand-coded entries follow */

  strcpy(list_[108].label, "X(5,0)");
  list_[108].filenum = 115;
  list_[108].irrep = 0;
  list_[108].pqnum =  5;
  list_[108].rsnum =  0;
  list_[108].priority = 99999;
  list_[108].next = &(list_[109]);
  list_[108].last = &(list_[107]);

  strcpy(list_[109].label, "X(0,5) 1");
  list_[109].filenum = 115;
  list_[109].irrep = 0;
  list_[109].pqnum =  0;
  list_[109].rsnum =  5;
  list_[109].priority = 99999;
  list_[109].next = &(list_[110]);
  list_[109].last = &(list_[108]);

  strcpy(list_[110].label, "X(0,5) 2");
  list_[110].filenum = 115;
  list_[110].irrep = 0;
  list_[110].pqnum =  0;
  list_[110].rsnum =  5;
  list_[110].priority = 99999;
  list_[110].next = &(list_[111]);
  list_[110].last = &(list_[109]);

  strcpy(list_[111].label, "X(0,5) 3");
  list_[111].filenum = 115;
  list_[111].irrep = 0;
  list_[111].pqnum =  0;
  list_[111].rsnum =  5;
  list_[111].priority = 99999;
  list_[111].next = &(list_[112]);
  list_[111].last = &(list_[110]);

  strcpy(list_[112].label, "X(0,5) 4");
  list_[112].filenum = 115;
  list_[112].irrep = 0;
  list_[112].pqnum =  0;
  list_[112].rsnum =  5;
  list_[112].priority = 99999;
  list_[112].next = NULL;
  list_[112].last = &(list_[111]);

  return(&(list_[0]));
}
}} // namespace psi::ccenergy
