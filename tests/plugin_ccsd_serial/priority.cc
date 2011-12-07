/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <ccfiles.h>
#include <libdpd/dpd.h>

namespace psi { namespace ccenergy {
void hmm();

#define NUM_ENTRIES 113
struct dpd_file4_cache_entry list[NUM_ENTRIES];

struct dpd_file4_cache_entry *priority_list(void)
{
  extern struct dpd_file4_cache_entry list[NUM_ENTRIES];

  strcpy(list[  0].label, "D <ij||ab> (i>j,a>b)");
  list[  0].filenum = 105;
  list[  0].irrep = 0;
  list[  0].pqnum =  2;
  list[  0].rsnum =  7;
  list[  0].priority = 231;
  list[  0].next = &(list[  0+1]);
  list[  0].last = NULL;

  strcpy(list[  1].label, "tIJAB");
  list[  1].filenum = 109;
  list[  1].irrep = 0;
  list[  1].pqnum =  2;
  list[  1].rsnum =  7;
  list[  1].priority = 1163;
  list[  1].next = &(list[  1+1]);
  list[  1].last = &(list[  1-1]);

  strcpy(list[  2].label, "tijab");
  list[  2].filenum = 109;
  list[  2].irrep = 0;
  list[  2].pqnum =  2;
  list[  2].rsnum =  7;
  list[  2].priority = 1163;
  list[  2].next = &(list[  2+1]);
  list[  2].last = &(list[  2-1]);

  strcpy(list[  3].label, "dIJAB");
  list[  3].filenum = 108;
  list[  3].irrep = 0;
  list[  3].pqnum =  1;
  list[  3].rsnum =  6;
  list[  3].priority = 77;
  list[  3].next = &(list[  3+1]);
  list[  3].last = &(list[  3-1]);

  strcpy(list[  4].label, "dijab");
  list[  4].filenum = 108;
  list[  4].irrep = 0;
  list[  4].pqnum =  1;
  list[  4].rsnum =  6;
  list[  4].priority = 77;
  list[  4].next = &(list[  4+1]);
  list[  4].last = &(list[  4-1]);

  strcpy(list[  5].label, "D <ij|ab>");
  list[  5].filenum = 105;
  list[  5].irrep = 0;
  list[  5].pqnum =  0;
  list[  5].rsnum =  5;
  list[  5].priority = 535;
  list[  5].next = &(list[  5+1]);
  list[  5].last = &(list[  5-1]);

  strcpy(list[  6].label, "tIjAb");
  list[  6].filenum = 109;
  list[  6].irrep = 0;
  list[  6].pqnum =  0;
  list[  6].rsnum =  5;
  list[  6].priority = 1241;
  list[  6].next = &(list[  6+1]);
  list[  6].last = &(list[  6-1]);

  strcpy(list[  7].label, "dIjAb");
  list[  7].filenum = 108;
  list[  7].irrep = 0;
  list[  7].pqnum =  0;
  list[  7].rsnum =  5;
  list[  7].priority = 77;
  list[  7].next = &(list[  7+1]);
  list[  7].last = &(list[  7-1]);

  strcpy(list[  8].label, "tauIJAB");
  list[  8].filenum = 109;
  list[  8].irrep = 0;
  list[  8].pqnum =  2;
  list[  8].rsnum =  7;
  list[  8].priority = 1161;
  list[  8].next = &(list[  8+1]);
  list[  8].last = &(list[  8-1]);

  strcpy(list[  9].label, "tauijab");
  list[  9].filenum = 109;
  list[  9].irrep = 0;
  list[  9].pqnum =  2;
  list[  9].rsnum =  7;
  list[  9].priority = 1161;
  list[  9].next = &(list[  9+1]);
  list[  9].last = &(list[  9-1]);

  strcpy(list[ 10].label, "tauIjAb");
  list[ 10].filenum = 109;
  list[ 10].irrep = 0;
  list[ 10].pqnum =  0;
  list[ 10].rsnum =  5;
  list[ 10].priority = 1161;
  list[ 10].next = &(list[ 10+1]);
  list[ 10].last = &(list[ 10-1]);

  strcpy(list[ 11].label, "tauIjbA");
  list[ 11].filenum = 109;
  list[ 11].irrep = 0;
  list[ 11].pqnum =  0;
  list[ 11].rsnum =  5;
  list[ 11].priority = 543;
  list[ 11].next = &(list[ 11+1]);
  list[ 11].last = &(list[ 11-1]);

  strcpy(list[ 12].label, "tauiJaB");
  list[ 12].filenum = 109;
  list[ 12].irrep = 0;
  list[ 12].pqnum =  0;
  list[ 12].rsnum =  5;
  list[ 12].priority = 389;
  list[ 12].next = &(list[ 12+1]);
  list[ 12].last = &(list[ 12-1]);

  strcpy(list[ 13].label, "tautIJAB");
  list[ 13].filenum = 109;
  list[ 13].irrep = 0;
  list[ 13].pqnum =  2;
  list[ 13].rsnum =  7;
  list[ 13].priority = 931;
  list[ 13].next = &(list[ 13+1]);
  list[ 13].last = &(list[ 13-1]);

  strcpy(list[ 14].label, "tautijab");
  list[ 14].filenum = 109;
  list[ 14].irrep = 0;
  list[ 14].pqnum =  2;
  list[ 14].rsnum =  7;
  list[ 14].priority = 931;
  list[ 14].next = &(list[ 14+1]);
  list[ 14].last = &(list[ 14-1]);

  strcpy(list[ 15].label, "tautIjAb");
  list[ 15].filenum = 109;
  list[ 15].irrep = 0;
  list[ 15].pqnum =  0;
  list[ 15].rsnum =  5;
  list[ 15].priority = 931;
  list[ 15].next = &(list[ 15+1]);
  list[ 15].last = &(list[ 15-1]);

  strcpy(list[ 16].label, "tiJaB");
  list[ 16].filenum = 109;
  list[ 16].irrep = 0;
  list[ 16].pqnum =  0;
  list[ 16].rsnum =  5;
  list[ 16].priority = 541;
  list[ 16].next = &(list[ 16+1]);
  list[ 16].last = &(list[ 16-1]);

  strcpy(list[ 17].label, "tIAJB");
  list[ 17].filenum = 109;
  list[ 17].irrep = 0;
  list[ 17].pqnum = 10;
  list[ 17].rsnum = 10;
  list[ 17].priority = 693;
  list[ 17].next = &(list[ 17+1]);
  list[ 17].last = &(list[ 17-1]);

  strcpy(list[ 18].label, "tiajb");
  list[ 18].filenum = 109;
  list[ 18].irrep = 0;
  list[ 18].pqnum = 10;
  list[ 18].rsnum = 10;
  list[ 18].priority = 693;
  list[ 18].next = &(list[ 18+1]);
  list[ 18].last = &(list[ 18-1]);

  strcpy(list[ 19].label, "tIAjb");
  list[ 19].filenum = 109;
  list[ 19].irrep = 0;
  list[ 19].pqnum = 10;
  list[ 19].rsnum = 10;
  list[ 19].priority = 849;
  list[ 19].next = &(list[ 19+1]);
  list[ 19].last = &(list[ 19-1]);

  strcpy(list[ 20].label, "tiaJB");
  list[ 20].filenum = 109;
  list[ 20].irrep = 0;
  list[ 20].pqnum = 10;
  list[ 20].rsnum = 10;
  list[ 20].priority = 693;
  list[ 20].next = &(list[ 20+1]);
  list[ 20].last = &(list[ 20-1]);

  strcpy(list[ 21].label, "tIbjA");
  list[ 21].filenum = 109;
  list[ 21].irrep = 0;
  list[ 21].pqnum = 10;
  list[ 21].rsnum = 10;
  list[ 21].priority = 619;
  list[ 21].next = &(list[ 21+1]);
  list[ 21].last = &(list[ 21-1]);

  strcpy(list[ 22].label, "tjAIb");
  list[ 22].filenum = 109;
  list[ 22].irrep = 0;
  list[ 22].pqnum = 10;
  list[ 22].rsnum = 10;
  list[ 22].priority = 541;
  list[ 22].next = &(list[ 22+1]);
  list[ 22].last = &(list[ 22-1]);

  strcpy(list[ 23].label, "C <ia||jb> (ia,bj)");
  list[ 23].filenum = 104;
  list[ 23].irrep = 0;
  list[ 23].pqnum = 10;
  list[ 23].rsnum = 11;
  list[ 23].priority = 75;
  list[ 23].next = &(list[ 23+1]);
  list[ 23].last = &(list[ 23-1]);

  strcpy(list[ 24].label, "WMBEJ");
  list[ 24].filenum = 115;
  list[ 24].irrep = 0;
  list[ 24].pqnum = 10;
  list[ 24].rsnum = 11;
  list[ 24].priority = 1823;
  list[ 24].next = &(list[ 24+1]);
  list[ 24].last = &(list[ 24-1]);

  strcpy(list[ 25].label, "Wmbej");
  list[ 25].filenum = 115;
  list[ 25].irrep = 0;
  list[ 25].pqnum = 10;
  list[ 25].rsnum = 11;
  list[ 25].priority = 1823;
  list[ 25].next = &(list[ 25+1]);
  list[ 25].last = &(list[ 25-1]);

  strcpy(list[ 26].label, "C <ia|jb>");
  list[ 26].filenum = 104;
  list[ 26].irrep = 0;
  list[ 26].pqnum = 10;
  list[ 26].rsnum = 10;
  list[ 26].priority = 227;
  list[ 26].next = &(list[ 26+1]);
  list[ 26].last = &(list[ 26-1]);

  strcpy(list[ 27].label, "WmBEj");
  list[ 27].filenum = 115;
  list[ 27].irrep = 0;
  list[ 27].pqnum = 10;
  list[ 27].rsnum = 10;
  list[ 27].priority = 1823;
  list[ 27].next = &(list[ 27+1]);
  list[ 27].last = &(list[ 27-1]);

  strcpy(list[ 28].label, "WMbeJ");
  list[ 28].filenum = 115;
  list[ 28].irrep = 0;
  list[ 28].pqnum = 10;
  list[ 28].rsnum = 10;
  list[ 28].priority = 1823;
  list[ 28].next = &(list[ 28+1]);
  list[ 28].last = &(list[ 28-1]);

  strcpy(list[ 29].label, "D <ij|ab> (ib,aj)");
  list[ 29].filenum = 105;
  list[ 29].irrep = 0;
  list[ 29].pqnum = 10;
  list[ 29].rsnum = 11;
  list[ 29].priority = 227;
  list[ 29].next = &(list[ 29+1]);
  list[ 29].last = &(list[ 29-1]);

  strcpy(list[ 30].label, "WMbEj");
  list[ 30].filenum = 115;
  list[ 30].irrep = 0;
  list[ 30].pqnum = 10;
  list[ 30].rsnum = 11;
  list[ 30].priority = 1823;
  list[ 30].next = &(list[ 30+1]);
  list[ 30].last = &(list[ 30-1]);

  strcpy(list[ 31].label, "WmBeJ");
  list[ 31].filenum = 115;
  list[ 31].irrep = 0;
  list[ 31].pqnum = 10;
  list[ 31].rsnum = 11;
  list[ 31].priority = 1823;
  list[ 31].next = &(list[ 31+1]);
  list[ 31].last = &(list[ 31-1]);

  strcpy(list[ 32].label, "F <ia||bc> (ia,b>c)");
  list[ 32].filenum = 107;
  list[ 32].irrep = 0;
  list[ 32].pqnum = 10;
  list[ 32].rsnum =  7;
  list[ 32].priority = 455;
  list[ 32].next = &(list[ 32+1]);
  list[ 32].last = &(list[ 32-1]);

  strcpy(list[ 33].label, "F <ia|bc>");
  list[ 33].filenum = 107;
  list[ 33].irrep = 0;
  list[ 33].pqnum = 10;
  list[ 33].rsnum =  5;
  list[ 33].priority = 379;
  list[ 33].next = &(list[ 33+1]);
  list[ 33].last = &(list[ 33-1]);

  strcpy(list[ 34].label, "E <ij||ka> (i>j,ak)");
  list[ 34].filenum = 106;
  list[ 34].irrep = 0;
  list[ 34].pqnum =  2;
  list[ 34].rsnum = 11;
  list[ 34].priority = 75;
  list[ 34].next = &(list[ 34+1]);
  list[ 34].last = &(list[ 34-1]);

  strcpy(list[ 35].label, "E <ai|jk>");
  list[ 35].filenum = 106;
  list[ 35].irrep = 0;
  list[ 35].pqnum = 11;
  list[ 35].rsnum =  0;
  list[ 35].priority = 759;
  list[ 35].next = &(list[ 35+1]);
  list[ 35].last = &(list[ 35-1]);

  strcpy(list[ 36].label, "E <ij|ka>");
  list[ 36].filenum = 106;
  list[ 36].irrep = 0;
  list[ 36].pqnum =  0;
  list[ 36].rsnum = 10;
  list[ 36].priority = 151;
  list[ 36].next = &(list[ 36+1]);
  list[ 36].last = &(list[ 36-1]);

  strcpy(list[ 37].label, "WMBEJ");
  list[ 37].filenum = 111;
  list[ 37].irrep = 0;
  list[ 37].pqnum = 10;
  list[ 37].rsnum = 10;
  list[ 37].priority = 1899;
  list[ 37].next = &(list[ 37+1]);
  list[ 37].last = &(list[ 37-1]);

  strcpy(list[ 38].label, "Wmbej");
  list[ 38].filenum = 111;
  list[ 38].irrep = 0;
  list[ 38].pqnum = 10;
  list[ 38].rsnum = 10;
  list[ 38].priority = 1899;
  list[ 38].next = &(list[ 38+1]);
  list[ 38].last = &(list[ 38-1]);

  strcpy(list[ 39].label, "WMbEj");
  list[ 39].filenum = 111;
  list[ 39].irrep = 0;
  list[ 39].pqnum = 10;
  list[ 39].rsnum = 10;
  list[ 39].priority = 1899;
  list[ 39].next = &(list[ 39+1]);
  list[ 39].last = &(list[ 39-1]);

  strcpy(list[ 40].label, "WmBeJ");
  list[ 40].filenum = 111;
  list[ 40].irrep = 0;
  list[ 40].pqnum = 10;
  list[ 40].rsnum = 10;
  list[ 40].priority = 1899;
  list[ 40].next = &(list[ 40+1]);
  list[ 40].last = &(list[ 40-1]);

  strcpy(list[ 41].label, "WMbeJ");
  list[ 41].filenum = 111;
  list[ 41].irrep = 0;
  list[ 41].pqnum = 10;
  list[ 41].rsnum = 10;
  list[ 41].priority = 1519;
  list[ 41].next = &(list[ 41+1]);
  list[ 41].last = &(list[ 41-1]);

  strcpy(list[ 42].label, "WmBEj");
  list[ 42].filenum = 111;
  list[ 42].irrep = 0;
  list[ 42].pqnum = 10;
  list[ 42].rsnum = 10;
  list[ 42].priority = 1519;
  list[ 42].next = &(list[ 42+1]);
  list[ 42].last = &(list[ 42-1]);

  strcpy(list[ 43].label, "D <ij||ab> (ia,jb)");
  list[ 43].filenum = 105;
  list[ 43].irrep = 0;
  list[ 43].pqnum = 10;
  list[ 43].rsnum = 10;
  list[ 43].priority = 303;
  list[ 43].next = &(list[ 43+1]);
  list[ 43].last = &(list[ 43-1]);

  strcpy(list[ 44].label, "D <ij|ab> (ia,jb)");
  list[ 44].filenum = 105;
  list[ 44].irrep = 0;
  list[ 44].pqnum = 10;
  list[ 44].rsnum = 10;
  list[ 44].priority = 303;
  list[ 44].next = &(list[ 44+1]);
  list[ 44].last = &(list[ 44-1]);

  strcpy(list[ 45].label, "Y (ME,JN)");
  list[ 45].filenum = 115;
  list[ 45].irrep = 0;
  list[ 45].pqnum = 10;
  list[ 45].rsnum =  0;
  list[ 45].priority = 4103;
  list[ 45].next = &(list[ 45+1]);
  list[ 45].last = &(list[ 45-1]);

  strcpy(list[ 46].label, "D <ij||ab> (ia,bj)");
  list[ 46].filenum = 105;
  list[ 46].irrep = 0;
  list[ 46].pqnum = 10;
  list[ 46].rsnum = 11;
  list[ 46].priority = 151;
  list[ 46].next = &(list[ 46+1]);
  list[ 46].last = &(list[ 46-1]);

  strcpy(list[ 47].label, "D <ij|ab> (ia,bj)");
  list[ 47].filenum = 105;
  list[ 47].irrep = 0;
  list[ 47].pqnum = 10;
  list[ 47].rsnum = 11;
  list[ 47].priority = 151;
  list[ 47].next = &(list[ 47+1]);
  list[ 47].last = &(list[ 47-1]);

  strcpy(list[ 48].label, "D <ij|ab> (ib,ja)");
  list[ 48].filenum = 105;
  list[ 48].irrep = 0;
  list[ 48].pqnum = 10;
  list[ 48].rsnum = 10;
  list[ 48].priority = 303;
  list[ 48].next = &(list[ 48+1]);
  list[ 48].last = &(list[ 48-1]);

  strcpy(list[ 49].label, "D <ij||ab>");
  list[ 49].filenum = 105;
  list[ 49].irrep = 0;
  list[ 49].pqnum =  0;
  list[ 49].rsnum =  5;
  list[ 49].priority = 75;
  list[ 49].next = &(list[ 49+1]);
  list[ 49].last = &(list[ 49-1]);

  strcpy(list[ 50].label, "D <ij||ab> (i>j,ab)");
  list[ 50].filenum = 105;
  list[ 50].irrep = 0;
  list[ 50].pqnum =  2;
  list[ 50].rsnum =  5;
  list[ 50].priority = 75;
  list[ 50].next = &(list[ 50+1]);
  list[ 50].last = &(list[ 50-1]);

  strcpy(list[ 51].label, "D <ij||ab> (ij,a>b)");
  list[ 51].filenum = 105;
  list[ 51].irrep = 0;
  list[ 51].pqnum =  0;
  list[ 51].rsnum =  7;
  list[ 51].priority = 75;
  list[ 51].next = &(list[ 51+1]);
  list[ 51].last = &(list[ 51-1]);

  strcpy(list[ 52].label, "C <ia||jb>");
  list[ 52].filenum = 104;
  list[ 52].irrep = 0;
  list[ 52].pqnum = 10;
  list[ 52].rsnum = 10;
  list[ 52].priority = 227;
  list[ 52].next = &(list[ 52+1]);
  list[ 52].last = &(list[ 52-1]);

  strcpy(list[ 53].label, "A <ij|kl>");
  list[ 53].filenum = 102;
  list[ 53].irrep = 0;
  list[ 53].pqnum =  0;
  list[ 53].rsnum =  0;
  list[ 53].priority = 151;
  list[ 53].next = &(list[ 53+1]);
  list[ 53].last = &(list[ 53-1]);

  strcpy(list[ 54].label, "WMNIJ");
  list[ 54].filenum = 111;
  list[ 54].irrep = 0;
  list[ 54].pqnum =  2;
  list[ 54].rsnum =  2;
  list[ 54].priority = 8891;
  list[ 54].next = &(list[ 54+1]);
  list[ 54].last = &(list[ 54-1]);

  strcpy(list[ 55].label, "Wmnij");
  list[ 55].filenum = 111;
  list[ 55].irrep = 0;
  list[ 55].pqnum =  2;
  list[ 55].rsnum =  2;
  list[ 55].priority = 8891;
  list[ 55].next = &(list[ 55+1]);
  list[ 55].last = &(list[ 55-1]);

  strcpy(list[ 56].label, "WMnIj");
  list[ 56].filenum = 111;
  list[ 56].irrep = 0;
  list[ 56].pqnum =  0;
  list[ 56].rsnum =  0;
  list[ 56].priority = 2127;
  list[ 56].next = &(list[ 56+1]);
  list[ 56].last = &(list[ 56-1]);

  strcpy(list[ 57].label, "E <ij||ka> (i>j,ka)");
  list[ 57].filenum = 106;
  list[ 57].irrep = 0;
  list[ 57].pqnum =  2;
  list[ 57].rsnum = 10;
  list[ 57].priority = 75;
  list[ 57].next = &(list[ 57+1]);
  list[ 57].last = &(list[ 57-1]);

  strcpy(list[ 58].label, "W (MN,IJ)");
  list[ 58].filenum = 115;
  list[ 58].irrep = 0;
  list[ 58].pqnum =  2;
  list[ 58].rsnum =  0;
  list[ 58].priority = 2583;
  list[ 58].next = &(list[ 58+1]);
  list[ 58].last = &(list[ 58-1]);

  strcpy(list[ 59].label, "ZIJMA");
  list[ 59].filenum = 113;
  list[ 59].irrep = 0;
  list[ 59].pqnum =  2;
  list[ 59].rsnum = 10;
  list[ 59].priority = 455;
  list[ 59].next = &(list[ 59+1]);
  list[ 59].last = &(list[ 59-1]);

  strcpy(list[ 60].label, "Zijma");
  list[ 60].filenum = 113;
  list[ 60].irrep = 0;
  list[ 60].pqnum =  2;
  list[ 60].rsnum = 10;
  list[ 60].priority = 455;
  list[ 60].next = &(list[ 60+1]);
  list[ 60].last = &(list[ 60-1]);

  strcpy(list[ 61].label, "ZIjMa");
  list[ 61].filenum = 113;
  list[ 61].irrep = 0;
  list[ 61].pqnum =  0;
  list[ 61].rsnum = 10;
  list[ 61].priority = 455;
  list[ 61].next = &(list[ 61+1]);
  list[ 61].last = &(list[ 61-1]);

  strcpy(list[ 62].label, "ZIjmA");
  list[ 62].filenum = 113;
  list[ 62].irrep = 0;
  list[ 62].pqnum =  0;
  list[ 62].rsnum = 10;
  list[ 62].priority = 379;
  list[ 62].next = &(list[ 62+1]);
  list[ 62].last = &(list[ 62-1]);

  strcpy(list[ 63].label, "ZIJAM");
  list[ 63].filenum = 113;
  list[ 63].irrep = 0;
  list[ 63].pqnum =  2;
  list[ 63].rsnum = 11;
  list[ 63].priority = 455;
  list[ 63].next = &(list[ 63+1]);
  list[ 63].last = &(list[ 63-1]);

  strcpy(list[ 64].label, "Zijam");
  list[ 64].filenum = 113;
  list[ 64].irrep = 0;
  list[ 64].pqnum =  2;
  list[ 64].rsnum = 11;
  list[ 64].priority = 455;
  list[ 64].next = &(list[ 64+1]);
  list[ 64].last = &(list[ 64-1]);

  strcpy(list[ 65].label, "ZIjAm");
  list[ 65].filenum = 113;
  list[ 65].irrep = 0;
  list[ 65].pqnum =  0;
  list[ 65].rsnum = 11;
  list[ 65].priority = 455;
  list[ 65].next = &(list[ 65+1]);
  list[ 65].last = &(list[ 65-1]);

  strcpy(list[ 66].label, "New tIJAB");
  list[ 66].filenum = 109;
  list[ 66].irrep = 0;
  list[ 66].pqnum =  2;
  list[ 66].rsnum =  7;
  list[ 66].priority = 58571;
  list[ 66].next = &(list[ 66+1]);
  list[ 66].last = &(list[ 66-1]);

  strcpy(list[ 67].label, "New tijab");
  list[ 67].filenum = 109;
  list[ 67].irrep = 0;
  list[ 67].pqnum =  2;
  list[ 67].rsnum =  7;
  list[ 67].priority = 58571;
  list[ 67].next = &(list[ 67+1]);
  list[ 67].last = &(list[ 67-1]);

  strcpy(list[ 68].label, "New tIjAb");
  list[ 68].filenum = 109;
  list[ 68].irrep = 0;
  list[ 68].pqnum =  0;
  list[ 68].rsnum =  5;
  list[ 68].priority = 10843;
  list[ 68].next = &(list[ 68+1]);
  list[ 68].last = &(list[ 68-1]);

  strcpy(list[ 69].label, "T (I>J,AB)");
  list[ 69].filenum = 115;
  list[ 69].irrep = 0;
  list[ 69].pqnum =  2;
  list[ 69].rsnum =  5;
  list[ 69].priority = 7295;
  list[ 69].next = &(list[ 69+1]);
  list[ 69].last = &(list[ 69-1]);

  strcpy(list[ 70].label, "T (IJ,A>B)");
  list[ 70].filenum = 115;
  list[ 70].irrep = 0;
  list[ 70].pqnum =  0;
  list[ 70].rsnum =  7;
  list[ 70].priority = 4711;
  list[ 70].next = &(list[ 70+1]);
  list[ 70].last = &(list[ 70-1]);

  strcpy(list[ 71].label, "B <ab||cd> (a>b,c>d)");
  list[ 71].filenum = 103;
  list[ 71].irrep = 0;
  list[ 71].pqnum =  7;
  list[ 71].rsnum =  7;
  list[ 71].priority = 75;
  list[ 71].next = &(list[ 71+1]);
  list[ 71].last = &(list[ 71-1]);

  strcpy(list[ 72].label, "B <ab|cd>");
  list[ 72].filenum = 103;
  list[ 72].irrep = 0;
  list[ 72].pqnum =  5;
  list[ 72].rsnum =  5;
  list[ 72].priority = 75;
  list[ 72].next = &(list[ 72+1]);
  list[ 72].last = &(list[ 72-1]);

  strcpy(list[ 73].label, "Z(ab,ij)");
  list[ 73].filenum = 115;
  list[ 73].irrep = 0;
  list[ 73].pqnum =  7;
  list[ 73].rsnum =  2;
  list[ 73].priority = 683;
  list[ 73].next = &(list[ 73+1]);
  list[ 73].last = &(list[ 73-1]);

  strcpy(list[ 74].label, "Z(ij,ab)");
  list[ 74].filenum = 115;
  list[ 74].irrep = 0;
  list[ 74].pqnum =  2;
  list[ 74].rsnum =  7;
  list[ 74].priority = 911;
  list[ 74].next = &(list[ 74+1]);
  list[ 74].last = &(list[ 74-1]);

  strcpy(list[ 75].label, "Z(Ab,Ij)");
  list[ 75].filenum = 115;
  list[ 75].irrep = 0;
  list[ 75].pqnum =  5;
  list[ 75].rsnum =  0;
  list[ 75].priority = 379;
  list[ 75].next = &(list[ 75+1]);
  list[ 75].last = &(list[ 75-1]);

  strcpy(list[ 76].label, "Z(Ij,Ab)");
  list[ 76].filenum = 115;
  list[ 76].irrep = 0;
  list[ 76].pqnum =  0;
  list[ 76].rsnum =  5;
  list[ 76].priority = 455;
  list[ 76].next = &(list[ 76+1]);
  list[ 76].last = &(list[ 76-1]);

  strcpy(list[ 77].label, "T (JI,A>B)");
  list[ 77].filenum = 115;
  list[ 77].irrep = 0;
  list[ 77].pqnum =  0;
  list[ 77].rsnum =  7;
  list[ 77].priority = 911;
  list[ 77].next = &(list[ 77+1]);
  list[ 77].last = &(list[ 77-1]);

  strcpy(list[ 78].label, "F <ai|bc>");
  list[ 78].filenum = 107;
  list[ 78].irrep = 0;
  list[ 78].pqnum = 11;
  list[ 78].rsnum =  5;
  list[ 78].priority = 75;
  list[ 78].next = &(list[ 78+1]);
  list[ 78].last = &(list[ 78-1]);

  strcpy(list[ 79].label, "T (I>J,BA)");
  list[ 79].filenum = 115;
  list[ 79].irrep = 0;
  list[ 79].pqnum =  2;
  list[ 79].rsnum =  5;
  list[ 79].priority = 911;
  list[ 79].next = &(list[ 79+1]);
  list[ 79].last = &(list[ 79-1]);

  strcpy(list[ 80].label, "E <ia|jk>");
  list[ 80].filenum = 106;
  list[ 80].irrep = 0;
  list[ 80].pqnum = 10;
  list[ 80].rsnum =  0;
  list[ 80].priority = 75;
  list[ 80].next = &(list[ 80+1]);
  list[ 80].last = &(list[ 80-1]);

  strcpy(list[ 81].label, "T2 (IA,JB)");
  list[ 81].filenum = 115;
  list[ 81].irrep = 0;
  list[ 81].pqnum = 10;
  list[ 81].rsnum = 10;
  list[ 81].priority = 683;
  list[ 81].next = &(list[ 81+1]);
  list[ 81].last = &(list[ 81-1]);

  strcpy(list[ 82].label, "T2 (IJ,AB)");
  list[ 82].filenum = 115;
  list[ 82].irrep = 0;
  list[ 82].pqnum =  0;
  list[ 82].rsnum =  5;
  list[ 82].priority = 4103;
  list[ 82].next = &(list[ 82+1]);
  list[ 82].last = &(list[ 82-1]);

  strcpy(list[ 83].label, "T2 (JI,AB)");
  list[ 83].filenum = 115;
  list[ 83].irrep = 0;
  list[ 83].pqnum =  0;
  list[ 83].rsnum =  5;
  list[ 83].priority = 1367;
  list[ 83].next = &(list[ 83+1]);
  list[ 83].last = &(list[ 83-1]);

  strcpy(list[ 84].label, "T2 (IJ,BA)");
  list[ 84].filenum = 115;
  list[ 84].irrep = 0;
  list[ 84].pqnum =  0;
  list[ 84].rsnum =  5;
  list[ 84].priority = 1367;
  list[ 84].next = &(list[ 84+1]);
  list[ 84].last = &(list[ 84-1]);

  strcpy(list[ 85].label, "T2 (JI,BA)");
  list[ 85].filenum = 115;
  list[ 85].irrep = 0;
  list[ 85].pqnum =  0;
  list[ 85].rsnum =  5;
  list[ 85].priority = 1367;
  list[ 85].next = &(list[ 85+1]);
  list[ 85].last = &(list[ 85-1]);

  strcpy(list[ 86].label, "T2 (ia,jb)");
  list[ 86].filenum = 115;
  list[ 86].irrep = 0;
  list[ 86].pqnum = 10;
  list[ 86].rsnum = 10;
  list[ 86].priority = 683;
  list[ 86].next = &(list[ 86+1]);
  list[ 86].last = &(list[ 86-1]);

  strcpy(list[ 87].label, "T2 (ij,ab)");
  list[ 87].filenum = 115;
  list[ 87].irrep = 0;
  list[ 87].pqnum =  0;
  list[ 87].rsnum =  5;
  list[ 87].priority = 1367;
  list[ 87].next = &(list[ 87+1]);
  list[ 87].last = &(list[ 87-1]);

  strcpy(list[ 88].label, "T2 (ji,ab)");
  list[ 88].filenum = 115;
  list[ 88].irrep = 0;
  list[ 88].pqnum =  0;
  list[ 88].rsnum =  5;
  list[ 88].priority = 455;
  list[ 88].next = &(list[ 88+1]);
  list[ 88].last = &(list[ 88-1]);

  strcpy(list[ 89].label, "T2 (ij,ba)");
  list[ 89].filenum = 115;
  list[ 89].irrep = 0;
  list[ 89].pqnum =  0;
  list[ 89].rsnum =  5;
  list[ 89].priority = 455;
  list[ 89].next = &(list[ 89+1]);
  list[ 89].last = &(list[ 89-1]);

  strcpy(list[ 90].label, "T2 (ji,ba)");
  list[ 90].filenum = 115;
  list[ 90].irrep = 0;
  list[ 90].pqnum =  0;
  list[ 90].rsnum =  5;
  list[ 90].priority = 455;
  list[ 90].next = &(list[ 90+1]);
  list[ 90].last = &(list[ 90-1]);

  strcpy(list[ 91].label, "T2 (IA,jb)");
  list[ 91].filenum = 115;
  list[ 91].irrep = 0;
  list[ 91].pqnum = 10;
  list[ 91].rsnum = 10;
  list[ 91].priority = 1291;
  list[ 91].next = &(list[ 91+1]);
  list[ 91].last = &(list[ 91-1]);

  strcpy(list[ 92].label, "T2 (Ij,Ab) 1");
  list[ 92].filenum = 115;
  list[ 92].irrep = 0;
  list[ 92].pqnum =  0;
  list[ 92].rsnum =  5;
  list[ 92].priority = 455;
  list[ 92].next = &(list[ 92+1]);
  list[ 92].last = &(list[ 92-1]);

  strcpy(list[ 93].label, "T2 (Ib,jA)");
  list[ 93].filenum = 115;
  list[ 93].irrep = 0;
  list[ 93].pqnum = 10;
  list[ 93].rsnum = 10;
  list[ 93].priority = 683;
  list[ 93].next = &(list[ 93+1]);
  list[ 93].last = &(list[ 93-1]);

  strcpy(list[ 94].label, "T2 (Ij,Ab) 2");
  list[ 94].filenum = 115;
  list[ 94].irrep = 0;
  list[ 94].pqnum =  0;
  list[ 94].rsnum =  5;
  list[ 94].priority = 455;
  list[ 94].next = &(list[ 94+1]);
  list[ 94].last = &(list[ 94-1]);

  strcpy(list[ 95].label, "Y (MB,JI)");
  list[ 95].filenum = 115;
  list[ 95].irrep = 0;
  list[ 95].pqnum = 10;
  list[ 95].rsnum =  0;
  list[ 95].priority = 1367;
  list[ 95].next = &(list[ 95+1]);
  list[ 95].last = &(list[ 95-1]);

  strcpy(list[ 96].label, "T2 (AB,JI)");
  list[ 96].filenum = 115;
  list[ 96].irrep = 0;
  list[ 96].pqnum =  5;
  list[ 96].rsnum =  0;
  list[ 96].priority = 1367;
  list[ 96].next = &(list[ 96+1]);
  list[ 96].last = &(list[ 96-1]);

  strcpy(list[ 97].label, "Y (mA,jI)");
  list[ 97].filenum = 115;
  list[ 97].irrep = 0;
  list[ 97].pqnum = 10;
  list[ 97].rsnum =  0;
  list[ 97].priority = 683;
  list[ 97].next = &(list[ 97+1]);
  list[ 97].last = &(list[ 97-1]);

  strcpy(list[ 98].label, "T2 (bA,jI)");
  list[ 98].filenum = 115;
  list[ 98].irrep = 0;
  list[ 98].pqnum =  5;
  list[ 98].rsnum =  0;
  list[ 98].priority = 683;
  list[ 98].next = &(list[ 98+1]);
  list[ 98].last = &(list[ 98-1]);

  strcpy(list[ 99].label, "T2 (Ij,Ab)");
  list[ 99].filenum = 115;
  list[ 99].irrep = 0;
  list[ 99].pqnum =  0;
  list[ 99].rsnum =  5;
  list[ 99].priority = 1823;
  list[ 99].next = &(list[ 99+1]);
  list[ 99].last = &(list[ 99-1]);

  strcpy(list[100].label, "Y (Mb,Ij)");
  list[100].filenum = 115;
  list[100].irrep = 0;
  list[100].pqnum = 10;
  list[100].rsnum =  0;
  list[100].priority = 683;
  list[100].next = &(list[100+1]);
  list[100].last = &(list[100-1]);

  strcpy(list[101].label, "T2 (Ab,Ij)");
  list[101].filenum = 115;
  list[101].irrep = 0;
  list[101].pqnum =  5;
  list[101].rsnum =  0;
  list[101].priority = 683;
  list[101].next = &(list[101+1]);
  list[101].last = &(list[101-1]);

  strcpy(list[102].label, "Y(Mb,jI)");
  list[102].filenum = 115;
  list[102].irrep = 0;
  list[102].pqnum = 10;
  list[102].rsnum =  0;
  list[102].priority = 683;
  list[102].next = &(list[102+1]);
  list[102].last = &(list[102-1]);

  strcpy(list[103].label, "T2 (Ab,jI)");
  list[103].filenum = 115;
  list[103].irrep = 0;
  list[103].pqnum =  5;
  list[103].rsnum =  0;
  list[103].priority = 683;
  list[103].next = &(list[103+1]);
  list[103].last = &(list[103-1]);

  strcpy(list[104].label, "Y(mA,Ij)");
  list[104].filenum = 115;
  list[104].irrep = 0;
  list[104].pqnum = 10;
  list[104].rsnum =  0;
  list[104].priority = 683;
  list[104].next = &(list[104+1]);
  list[104].last = &(list[104-1]);

  strcpy(list[105].label, "T2 (bA,Ij)");
  list[105].filenum = 115;
  list[105].irrep = 0;
  list[105].pqnum =  5;
  list[105].rsnum =  0;
  list[105].priority = 683;
  list[105].next = &(list[105+1]);
  list[105].last = &(list[105-1]);

  strcpy(list[106].label, "T2(IJ,AB) DIIS");
  list[106].filenum = 115;
  list[106].irrep = 0;
  list[106].pqnum =  2;
  list[106].rsnum =  7;
  list[106].priority = 73;
  list[106].next = &(list[106+1]);
  list[106].last = &(list[106-1]);

  strcpy(list[107].label, "T2(Ij,Ab) DIIS");
  list[107].filenum = 115;
  list[107].irrep = 0;
  list[107].pqnum =  0;
  list[107].rsnum =  5;
  list[107].priority = 73;
  list[107].next = &(list[107+1]);
  list[107].last = &(list[107-1]);

  /* Hand-coded entries follow */

  strcpy(list[108].label, "X(5,0)");
  list[108].filenum = 115;
  list[108].irrep = 0;
  list[108].pqnum =  5;
  list[108].rsnum =  0;
  list[108].priority = 99999;
  list[108].next = &(list[109]);
  list[108].last = &(list[107]);

  strcpy(list[109].label, "X(0,5) 1");
  list[109].filenum = 115;
  list[109].irrep = 0;
  list[109].pqnum =  0;
  list[109].rsnum =  5;
  list[109].priority = 99999;
  list[109].next = &(list[110]);
  list[109].last = &(list[108]);

  strcpy(list[110].label, "X(0,5) 2");
  list[110].filenum = 115;
  list[110].irrep = 0;
  list[110].pqnum =  0;
  list[110].rsnum =  5;
  list[110].priority = 99999;
  list[110].next = &(list[111]);
  list[110].last = &(list[109]);

  strcpy(list[111].label, "X(0,5) 3");
  list[111].filenum = 115;
  list[111].irrep = 0;
  list[111].pqnum =  0;
  list[111].rsnum =  5;
  list[111].priority = 99999;
  list[111].next = &(list[112]);
  list[111].last = &(list[110]);

  strcpy(list[112].label, "X(0,5) 4");
  list[112].filenum = 115;
  list[112].irrep = 0;
  list[112].pqnum =  0;
  list[112].rsnum =  5;
  list[112].priority = 99999;
  list[112].next = NULL;
  list[112].last = &(list[111]);

  return(&(list[0]));
}
void hmm(){
  dpdbuf4 newt2;
  dpd_buf4_init(&newt2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_mat_irrep_init(&newt2, 0);
  dpd_buf4_mat_irrep_rd(&newt2, 0);
  for (long int a=0; a<19; a++){
      for (long int b=0; b<19; b++){
          for (long int i=0; i<5; i++){
              for (long int j=0; j<5; j++){
                  printf("elsewhere %5li %5li %5li %5li %20.12lf\n",a,b,i,j,newt2.matrix[0][i*5+j][a*19+b]);fflush(stdout);
              }
          }
      }
  }
  dpd_buf4_mat_irrep_close(&newt2, 0);
  dpd_buf4_close(&newt2);

}
}} // namespace psi::ccenergy
