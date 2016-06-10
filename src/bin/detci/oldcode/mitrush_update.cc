/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

/* need to fix to add phase.  should probably move most of code
** to a mitrush_blk_upd() function in C
*/
void mitrush_update(CIvect &Cc, CIvect &Cl, CIvect &Sc, CIvect &Sl,
      double norm, double acur, double alast)
{
   double *ccptr, *clptr, *scptr, *slptr;
   int i,al,bl,ai,bi,blk;
   int cac,cbc,len;

   if (Cc.icore == 1 && Cl.icore == 1 && Sc.icore == 1 && Sl.icore == 1) {
      if (Cc.Ms0) {
         for (blk=0; blk < Cc.num_blocks; blk++) {
            cac = Cc.Ia_code[blk];
            cbc = Cc.Ib_code[blk];
            if (cac < cbc) continue;
            if (cac != cbc) {
               len = Cc.Ia_size[blk] * Cc.Ib_size[blk];
               ccptr = Cc.blocks[blk][0];
               clptr = Cl.blocks[blk][0];
               scptr = Sc.blocks[blk][0];
               slptr = Sl.blocks[blk][0];
               for (i=0; i<len; i++) {
                  ccptr[i] = norm * (acur * ccptr[i] + alast * clptr[i]);
                  scptr[i] = norm * (acur * scptr[i] + alast * slptr[i]);
                  }
               for (i=0; i<H0block.size; i++) {
                  al = H0block.alplist[i];
                  bl = H0block.betlist[i];
                  ai = H0block.alpidx[i];
                  bi = H0block.betidx[i];
                  if (al == cac && bl == cbc) {
                     H0block.c0b[i] = Cc.blocks[blk][ai][bi];
                     H0block.s0b[i] = Sc.blocks[blk][ai][bi];
                     } 
                  else if (al == cbc && bl == cac) {
                     H0block.c0b[i] = Cc.blocks[blk][bi][ai];
                     H0block.s0b[i] = Sc.blocks[blk][bi][ai];
                     }
                  }
               }
            if (cac == cbc) {
               len = Cc.Ia_size[blk] * Cc.Ib_size[blk];
               ccptr = Cc.blocks[blk][0];
               clptr = Cl.blocks[blk][0];
               scptr = Sc.blocks[blk][0];
               slptr = Sl.blocks[blk][0];
               for (i=0; i<len; i++) {
                  ccptr[i] = norm * (acur * ccptr[i] + alast * clptr[i]);
                  scptr[i] = norm * (acur * scptr[i] + alast * slptr[i]);
                  }
               for (i=0; i<H0block.size; i++) {
                  al = H0block.alplist[i];
                  bl = H0block.betlist[i];
                  ai = H0block.alpidx[i];
                  bi = H0block.betidx[i];
                  if (al == cac && bl == cbc) {
                     H0block.c0b[i] = Cc.blocks[blk][ai][bi];
                     H0block.s0b[i] = Sc.blocks[blk][ai][bi];
                     }
                  }
               }
            } /* end loop over blocks */
         } /* end Ms0 case */

      else { /* not Ms = 0 */
         for (blk=0; blk < Cc.num_blocks; blk++) {
            cac = Cc.Ia_code[blk];
            cbc = Cc.Ib_code[blk];
            len = Cc.Ia_size[blk] * Cc.Ib_size[blk];
            ccptr = Cc.blocks[blk][0];
            clptr = Cl.blocks[blk][0];
            scptr = Sc.blocks[blk][0];
            slptr = Sl.blocks[blk][0];
            for (i=0; i<len; i++) {
               ccptr[i] = norm * (acur * ccptr[i] + alast * clptr[i]);
               scptr[i] = norm * (acur * scptr[i] + alast * slptr[i]);
               }
            for (i=0; i<H0block.size; i++) {
               al = H0block.alplist[i];
               bl = H0block.betlist[i];
               ai = H0block.alpidx[i];
               bi = H0block.betidx[i];
               if (al == cac && bl == cbc) {
                  H0block.c0b[i] = Cc.blocks[blk][ai][bi];
                  H0block.s0b[i] = Sc.blocks[blk][ai][bi];
                  }
               }
            } /* end loop over blocks */
         } /* end Ms != 0 case */

      } /* end all in-core option */
   else {
      printf("(mitrush_update): unrecognized icore option\n");
      return;
      }
}