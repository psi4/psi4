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

/*!
  \file
  \ingroup DETCI
  \brief Contains code to do block-to-block single replacement lists
*/

#include <cstdio>
#include <cstdlib>
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/detci/structs.h"

namespace psi { namespace detci {

/* DEFINES */
#define MAX_EL 30
#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

/* PROTOS */
extern int subgr_lex_addr(struct level *head, int *occs, int nel, int norb);

void b2brepl(unsigned char **occs, int *Jcnt, int **Jij, int **Joij,
   int **Jridx, signed char **Jsgn, struct olsen_graph *Graph,
   int Ilist, int Jlist, int len, struct calcinfo *Cinfo);
void b2bgen1(unsigned char **occs, int *Jcnt, int **Jij, int **Joij,
   int **Jridx, signed char **Jsgn, struct level *subgr_head,
   int len, int ijsym, int nel, int ras1_lvl, int ras3_lvl, int ras4_lvl,
   struct calcinfo *Cinfo);
void b2bgen2(unsigned char **occs, int *Jcnt, int **Jij, int **Joij,
   int **Jridx, signed char **Jsgn, struct level *subgr_head,
   int up, int down, int len, int ijsym, int nel, int ras1_lvl, int ras3_lvl,
   int ras4_lvl, struct calcinfo *Cinfo);


/*
** b2brepl
**
** Generate block to block single replacements on the fly.
** Generates replacements from strings I in list Ilist to strings J in
**    list Jlist.
**
** Parameters:
**    occs     = array of occupied orbitals for each string I
**    Jcnt     = array to hold number of replacements for each string I
**    Jij      = matrix of ij's to all J's for each I (Jij[I][J])
**    Joij     = matrix of olsen ij's to all J's for each I (Joij[I][J])
**    Jridx    = matrix of relative indices of resultant strings J,
**               i.e. J'[J] = Jridx[I][J] where J' is the proper address of J
**    Jsgn     = matrix of signs as above
**    Graph    = Olsen Graph for the relevant strings
**    Ilist    = list number for I's
**    Jlist    = list number for J's
**    len      = length of occs array (how many I strings)
**
** David Sherrill
** August 1995
*/
void b2brepl(unsigned char **occs, int *Jcnt, int **Jij, int **Joij,
      int **Jridx, signed char **Jsgn, struct olsen_graph *Graph,
      int Ilist, int Jlist, int len, struct calcinfo *Cinfo)
{
   int I_n1, I_n2, I_n3, I_n4;
   int J_n1, J_n2, J_n3, J_n4;
   int D_n1, D_n2, D_n3, D_n4;
   int Ilist_ir, Jlist_ir;
   int nel,Icode,Jcode,ijsym,up,down;
   struct level *subgr_head;

   /* zero out Jcnt so there's no mistake */
   zero_int_array(Jcnt, len);

   /* get pointer to subgraph */
   Icode = Ilist % Graph->subgr_per_irrep;
   Jcode = Jlist % Graph->subgr_per_irrep;
   Ilist_ir = Ilist / Graph->subgr_per_irrep;
   Jlist_ir = Jlist / Graph->subgr_per_irrep;
   subgr_head = Graph->sg[Jlist_ir][Jcode].lvl;
   /* first figure out how many electrons in RAS I, II, III, IV for
    * each of the blocks
    */
   nel = Graph->num_el_expl;
   I_n1 = Graph->encode[0][Icode];
   I_n3 = Graph->encode[1][Icode];
   I_n4 = Graph->encode[2][Icode];
   I_n2 = nel - I_n1 - I_n3 - I_n4;
   J_n1 = Graph->encode[0][Jcode];
   J_n3 = Graph->encode[1][Jcode];
   J_n4 = Graph->encode[2][Jcode];
   J_n2 = nel - J_n1 - J_n3 - J_n4;
   if (I_n1 < 0 || I_n2 < 0 || I_n3 < 0 || J_n1 < 0 || J_n2 < 0 || J_n3 < 0
          || I_n4 < 0 || J_n4 < 0) {
     outfile->Printf("b2brepl: got less than 1 electron in a partition\n");
      return;
      }

   /* now figure out the differences */
   D_n1 = J_n1 - I_n1;
   D_n2 = J_n2 - I_n2;
   D_n3 = J_n3 - I_n3;
   D_n4 = J_n4 - I_n4;

   /* are these ok? */
   if (abs(D_n1) + abs(D_n2) + abs(D_n3) + abs(D_n4) > 2)
      return;

   /* get ijsym */
   ijsym = Ilist_ir ^ Jlist_ir;

   /* figure out the case */
   if (D_n1 == 0 && D_n2 == 0 && D_n3 == 0 && D_n4 == 0) {
      b2bgen1(occs,Jcnt,Jij,Joij,Jridx,Jsgn,subgr_head,len,ijsym,nel,
         Graph->ras1_lvl, Graph->ras3_lvl, Graph->ras4_lvl, Cinfo);
      }
   else { /* figure out which is 1 and which is -1 */
      if (D_n1 == 1) up = 0;
      else if (D_n2 == 1) up = 1;
      else if (D_n3 == 1) up = 2;
      else if (D_n4 == 1) up = 3;
      if (D_n1 == -1) down = 0;
      else if (D_n2 == -1) down = 1;
      else if (D_n3 == -1) down = 2;
      else if (D_n4 == -1) down = 3;
      b2bgen2(occs,Jcnt,Jij,Joij,Jridx,Jsgn,subgr_head,up,down,len,ijsym,nel,
         Graph->ras1_lvl, Graph->ras3_lvl, Graph->ras4_lvl, Cinfo);
      }
}



/*
** b2bgen1: Generate all single replacements going to another block
**    in which the number of electrons in each RAS partition must
**    remain constant (i.e. staying in same code, maybe irrep changes)
**
*/
void b2bgen1(unsigned char **occs, int *Jcnt, int **Jij, int **Joij,
      int **Jridx, signed char **Jsgn, struct level *subgr_head,
      int len, int ijsym, int nel, int ras1_lvl, int ras3_lvl, int ras4_lvl,
      struct calcinfo *Cinfo)
{
   int I;
   int O[MAX_EL], T[MAX_EL], ras_occs[4][MAX_EL], ecnt[4];
   int i,j,k,l,m,ij,oij,orb,r1cnt,r2cnt,r3cnt,r4cnt;
   int isym,jsym;
   int cnt,ridx,norb;
   signed char sgn;
   int ras,hole,part,abshole,hops,iused;
   int **ras_orbs[4], **ras_opi;

   for (i=0; i<4; i++) ras_orbs[i] = Cinfo->ras_orbs[i];
   ras_opi = Cinfo->ras_opi;
   norb = Cinfo->num_ci_orbs;

   /* loop over strings */
   for (I=0; I<len; I++) {

      cnt = 0;  /* how many singlerepls found for this string */

      /* divide up electrons into their subspaces */
      r1cnt = r2cnt = r3cnt = r4cnt = 0;
      for (i=0; i<nel; i++) {
         orb = (int) occs[I][i];
         O[i] = orb;
         if (orb <= ras1_lvl) ras_occs[0][r1cnt++] = orb;
         else if (orb >= ras3_lvl && orb < ras4_lvl) ras_occs[2][r3cnt++] = orb;
         else if (orb >= ras4_lvl) ras_occs[3][r4cnt++] = orb;
         else ras_occs[1][r2cnt++] = orb;
         }

      ecnt[0] = r1cnt;
      ecnt[1] = r2cnt;
      ecnt[2] = r3cnt;
      ecnt[3] = r4cnt;

      /* do diagonals first */
      if (ijsym == 0) {
         ridx = subgr_lex_addr(subgr_head, O, nel, norb);
         if (ridx < 0) {
           outfile->Printf("b2bgen1: invalid string index = %d\n", ridx);
            continue;
            }
         for (k=0; k<nel; k++) {
            i = O[k];
            ij = ioff[i] + i;
            oij = i * norb + i;
            Jij[I][cnt] = ij;
            Joij[I][cnt] = oij;
            Jsgn[I][cnt] = 1;
            Jridx[I][cnt] = ridx;
            cnt++;
            }
         }

      /* Done with diagonals.  Now loop over RAS subspaces */
      for (ras=0; ras<4; ras++) {
         /* loop over excited electrons/orbitals */
         for (hole=0; hole<ecnt[ras]; hole++) {

            abshole = hole;
            for (k=0; k<ras; k++) abshole += ecnt[k];

            j = ras_occs[ras][hole];
            if (j < Cinfo->num_expl_cor_orbs) continue;
            jsym = Cinfo->orbsym[j + Cinfo->num_drc_orbs];
            isym = ijsym ^ jsym;
            for (part=0; part<Cinfo->ras_opi[ras][isym]; part++) {
               i = ras_orbs[ras][isym][part];

               /* make sure i is not occupied already */
               for (k=0,l=0; k<ecnt[ras]; k++) {
                  if (i == ras_occs[ras][k]) {l=1; break;}
                  }
               if (l==1) continue;

               /* arrange the occupied list in order */
               for (k=0,l=0; k<ras; k++) {
                  for (m=0; m<ecnt[k]; m++) {
                     T[l++] = ras_occs[k][m];
                     }
                  }
               for (k=0,iused=0; k<ecnt[ras]; k++) {
                  if (!iused && i < ras_occs[ras][k]) {
                     iused=1;
                     hops = l;
                     T[l++] = i;
                     }
                  if (k != hole) {
                     T[l++] = ras_occs[ras][k];
                     }
                  }
               if (!iused) {
                  hops = l;
                  T[l++] = i;
                  }
               for (k=ras+1; k<4; k++) {
                  for (m=0; m<ecnt[k]; m++) {
                     T[l++] = ras_occs[k][m];
                     }
                  }
               ridx = subgr_lex_addr(subgr_head, T, nel, norb);
               /* wouldn't be found if i is actually occupied */
               if (ridx < 0) {
                 outfile->Printf("b2bgen1: invalid string index = %d\n", ridx);
                  continue;
                  }

               /* get the sign */
               sgn = ((abshole + hops) % 2) ? -1 : 1;
               ij = ioff[MAX0(i,j)] + MIN0(i,j);
               oij = i * norb + j;

               /* store these results and increment the counter */
               Jij[I][cnt] = ij;
               Joij[I][cnt] = oij;
               Jsgn[I][cnt] = sgn;
               Jridx[I][cnt] = ridx;
               cnt++;

               } /* end loop over particles */
            } /* end loop over holes */
         } /* end loop over RAS subspaces */
      Jcnt[I] = cnt;
      } /* end loop over strings I */

}


void b2bgen2(unsigned char **occs, int *Jcnt, int **Jij, int **Joij,
      int **Jridx, signed char **Jsgn, struct level *subgr_head,
      int up, int down, int len, int ijsym, int nel, int ras1_lvl,
      int ras3_lvl, int ras4_lvl, struct calcinfo *Cinfo)
{
   int I;
   int O[MAX_EL], T[MAX_EL], ras_occs[4][MAX_EL], ecnt[4];
   int i,j,k,l,ij,oij,orb,r1cnt,r2cnt,r3cnt,r4cnt;
   int isym,jsym;
   int cnt,ridx,norb;
   signed char sgn;
   int hole,part,abshole,hops,iused;
   int **ras_orbs, *ras_opi, *ras_occs_excite, *ras_occs_virt;

   norb = Cinfo->num_ci_orbs;

   /* loop over strings */
   for (I=0; I<len; I++) {

      cnt = 0;  /* how many singlerepls found for this string */

      /* divide up electrons into their subspaces */
      r1cnt = r2cnt = r3cnt = r4cnt = 0;
      for (i=0; i<nel; i++) {
         orb = (int) occs[I][i];
         O[i] = orb;
         if (orb <= ras1_lvl) ras_occs[0][r1cnt++] = orb;
         else if (orb>=ras3_lvl && orb<ras4_lvl) ras_occs[2][r3cnt++] = orb;
         else if (orb >= ras4_lvl) ras_occs[3][r4cnt++] = orb;
         else ras_occs[1][r2cnt++] = orb;
         }

      ecnt[0] = r1cnt;
      ecnt[1] = r2cnt;
      ecnt[2] = r3cnt;
      ecnt[3] = r4cnt;
      ras_occs_excite = ras_occs[down];
      ras_occs_virt = ras_occs[up];
      ras_opi = Cinfo->ras_opi[up];
      ras_orbs = Cinfo->ras_orbs[up];

      for (hole=0; hole<ecnt[down]; hole++) {

         abshole = hole;
         for (k=0; k<down; k++) abshole += ecnt[k];

         j = ras_occs_excite[hole];
         if (j < Cinfo->num_expl_cor_orbs) continue;
         jsym = Cinfo->orbsym[j + Cinfo->num_drc_orbs];
         isym = ijsym ^ jsym;
         for (part=0; part<Cinfo->ras_opi[up][isym]; part++) {
            i = ras_orbs[isym][part];

            /* make sure i is not occupied already */
            for (k=0,l=0; k<ecnt[up]; k++) {
               if (i == ras_occs_virt[k]) {l=1; break;}
               }
            if (l==1) continue;

            /* arrange the occupied list in order */
            for (k=0,iused=0,l=0; k<nel; k++) {
               if (!iused && i < O[k]) { iused=1; hops = l; T[l++] = i; }
               if (k != abshole) T[l++] = O[k];
               }
            if (!iused) { hops = l; T[l++] = i; }

            ridx = subgr_lex_addr(subgr_head, T, nel, norb);
            /* wouldn't be found if i is actually occupied */
            if (ridx < 0) {
              outfile->Printf("b2bgen2: invalid string index = %d\n", ridx);
               continue;
               }

            /* get the sign */
            sgn = ((abshole + hops) % 2) ? -1 : 1;

            ij = ioff[MAX0(i,j)] + MIN0(i,j);
            oij = i * norb + j;

            /* store these results and increment the counter */
            Jij[I][cnt] = ij;
            Joij[I][cnt] = oij;
            Jsgn[I][cnt] = sgn;
            Jridx[I][cnt] = ridx;
            cnt++;

            } /* end loop over particles */
         } /* end loop over holes */

      Jcnt[I] = cnt;
      } /* end loop over strings I */

}


void b2brepl_test(unsigned char ***occs, int *Jcnt, int **Jij, int **Joij,
      int **Jridx, signed char **Jsgn, struct olsen_graph *Graph, struct calcinfo *Cinfo)
{
   int i, j, nirreps, ncodes;
   int Iirrep, Icode, Ilistnum;
   int Jirrep, Jcode, Jlistnum;
   struct stringgraph *Isubgraph, *Jsubgraph;

   nirreps = Graph->nirreps;
   ncodes = Graph->subgr_per_irrep;

   outfile->Printf("\nTesting block to block single-replacements b2brepl()\n");
   for (Iirrep=0,Ilistnum=0; Iirrep<nirreps; Iirrep++) {
      for (Icode=0; Icode<ncodes; Icode++,Ilistnum++) {
         Isubgraph = Graph->sg[Iirrep] + Icode;
         if (!Isubgraph->num_strings) continue;
         for (Jirrep=0,Jlistnum=0; Jirrep<nirreps; Jirrep++) {
            for (Jcode=0; Jcode<ncodes; Jcode++,Jlistnum++) {
               Jsubgraph = Graph->sg[Jirrep] + Jcode;
               if (!Jsubgraph->num_strings) continue;

               b2brepl(occs[Ilistnum], Jcnt, Jij, Joij, Jridx, Jsgn,
                  Graph, Ilistnum, Jlistnum, Isubgraph->num_strings, Cinfo);

               for (i=0; i<Isubgraph->num_strings; i++) {
                  outfile->Printf( "\nString %4d (",i);
                     for (j=0; j<Graph->num_el_expl; j++) {
                        outfile->Printf( "%2d ", (int) occs[Ilistnum][i][j]);
                        }
                  outfile->Printf( ")\n   Links:\n") ;
                  for (j=0; j<Jcnt[i]; j++) {
                     outfile->Printf( "   %3d [%3d] %c (%2d %3d)\n",
                        Jij[i][j], Joij[i][j], (Jsgn[i][j]==1) ? '+' : '-',
                        Jlistnum, Jridx[i][j]);
                     }
                  }

               } /* end loop over Jcodes */
            } /* end loop over Jirrep */
         } /* end loop over Icodes */
      } /* end loop over Iirrep */
}

}} // namespace psi::detci
