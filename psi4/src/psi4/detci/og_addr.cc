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
    \ingroup DETCI
    \brief Enter brief description of file here
*/

/*
** OG_ADDR.C: Code to calculate lexical addresses of strings for
**    olsen graph structures.
**
** C. David Sherrill, June 1995
**
*/

#include <cstdio>
#include "psi4/psi4-dec.h"
#include "psi4/detci/structs.h"

namespace psi { namespace detci {

/* FUNCTION PROTOTYPES */
int subgr_lex_addr(struct level *head, int *occs, int nel, int norb);

/*
** subgr_lex_addr():  Function takes a pointer to the head of a subgraph
**    and a list of occupied orbitals, and returns the lexical
**    address (within the given subgraph) of the walk.
**
** Parameters:
**    head  =  pointer to first level of subgraph (struct level *)
**    occs  =  integer array containing a list of occupied orbitals
**    nel   =  number of electrons in the walk (or length of occs[])
**    norb  =  number of orbitals in subgraph
**
** Returns: the lexical address of the walk, OR -1 if not found
*/
int subgr_lex_addr(struct level *head, int *occs, int nel, int norb)
{
   int i=0, j=1, c=0;
   int addr=0;
   struct level *curr;

   curr = head;

   while (i < norb) {
      if (c < nel && occs[c] == i) {
         addr += curr->y[j-1];
         j = curr->k[1][j-1];
         c++;
         }
      else {
         j = curr->k[0][j-1];
         }
      if (j == 0) {
         outfile->Printf( "(subgr_lex_addr): Impossible walk!\n");
         return(-1);
         }
      i++;
      curr++;
      }

   return(addr);
}


/*
** og_lex_addr():  Function determines the lexical address for a given
**    string.  Also returns the id for the string list containing the
**    given string.
**
** Parameters:
**    Graph    =  pointer to olsen graph
**    occs     =  array holding orbital numbers for occupied orbitals
**    nel      =  number of explicit electrons (i.e. length of occs array)
**    listnum  =  ptr to hold id of the list containing the given string
**
** Returns: the relative index within the given list, OR -1 if not found
**
*/
int og_lex_addr(struct olsen_graph *Graph, int *occs, int nel,
      int *listnum)
{
   int i,j,irrep,code;
   int inras1 = 0, inras2 = 0, inras3 = 0, inras4 = 0;
   int addr;
   int *orbsym;
   struct stringgraph *subgraph;

   irrep = Graph->drc_sym;
   orbsym = Graph->orbsym + Graph->num_drc_orbs;

   for (i=0; i<nel; i++) {
      j = occs[i];
      irrep ^= orbsym[j];
      if (j <= Graph->ras1_lvl) inras1++;
      else if (j >= Graph->ras3_lvl && j < Graph->ras4_lvl) inras3++;
      else if (j >= Graph->ras4_lvl) inras4++;
      else inras2++;
      }
   // inras1 += Graph->num_drc_orbs; // CDS 4/15
   inras1 -= Graph->ras1_min;
   if (inras1 < 0) return(-1);
   if (inras3 > Graph->ras3_max) return(-1);
   if (inras4 > Graph->ras4_max) return(-1);
   code = Graph->decode[inras1][inras3][inras4];
   if (code < 0) return(-1);

   subgraph = Graph->sg[irrep] + code;
   if (subgraph->num_strings < 1) return(-1);

   *listnum = irrep * Graph->subgr_per_irrep + code;

   addr = subgr_lex_addr(subgraph->lvl, occs, nel, Graph->num_orb);
   return(addr);

}



/*
** str_abs2rel(): Function returns the relative index and string list
**    number corresponding to the given absolute string number.
**
** Parameters:
**    absidx  = the absolute index
**    relidx  = ptr to hold the relative index
**    listnum = ptr to hold the code for the list holding this string
**    Graph   = olsen graph
*/
void str_abs2rel(int absidx, int *relidx, int *listnum,
      struct olsen_graph *Graph)
{
   int tot=0;
   int irrep, code;

   for (irrep=0; irrep<Graph->nirreps; irrep++) {
      if (tot + Graph->str_per_irrep[irrep] > absidx) break;
      else tot += Graph->str_per_irrep[irrep];
      }

   for (code=0; code<Graph->subgr_per_irrep; code++) {
      if (tot + Graph->sg[irrep][code].num_strings > absidx) break;
      else tot += Graph->sg[irrep][code].num_strings;
      }
   *relidx = absidx - tot;
   *listnum = irrep * Graph->subgr_per_irrep + code;
}


/*
** str_rel2abs(): Function returns the absolute index for a string
**    when given a list number and a relative index within that list.
**
** Parameters:
**    relidx  = the relative index
**    listnum = the code for the list holding this string
**    Graph   = olsen graph
**
** Returns: the absolute index
*/
int str_rel2abs(int relidx, int listnum, struct olsen_graph *Graph)
{
   int tot=0;

   tot = Graph->list_offset[listnum] + relidx;

   return(tot);
}

}} // namespace psi::detci
