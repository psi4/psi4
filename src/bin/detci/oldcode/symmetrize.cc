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

/*
** CIvect::symmetrize(): This function symmetrizes the CI vector
**    to maintain the appropriate spin symmetry.  It is the same
**    as the symnorm function except that it does no normalization.
**    Assume that this function is called only if Ms=0
**
** Parameters:
**    phase = the exponent S in C(Ia,Ib) = (-1)^S * C(Ib,Ia)
**
*/
void CIvect::symmetrize(double phase)
{
   int i,j;
   int blk,buf,irrep,ac,bc,len,upper;
   double **mat,*arr;
   double phase;

   if (icore == 1) {

      read(cur_vect, 0);
      for (blk=0; blk<num_blocks; blk++) {
         ac = Ia_code[blk];
         bc = Ib_code[blk];
         mat = blocks[blk];
         if (ac == bc) { /* diagonal block */
            for (i=0; i<Ia_size[blk]; i++) {
               for (j=0; j<i; j++) {
                  mat[j][i] = mat[i][j] * phase;
                  }
               }
            }
         if (ac > bc) { /* off-diagonal block */
            upper = decode[bc][ac];
            if (upper >= 0) {
               for (i=0; i<Ia_size[blk]; i++) {
                  for (j=0; j<Ib_size[blk]; j++) {
                     blocks[upper][j][i] = mat[i][j] * phase;
                     }
                  }
               }
            }
         } /* end loop over blocks */

      write(cur_vect, 0); 

      } /* end icore == 1 */


   else if (icore == 2) { /* irrep at a time */

      for (buf=0; buf<buf_per_vect; buf++) {
         irrep = buf2blk[buf];
         if (buf_offdiag[buf]) continue;

         /* diagonal irrep, symmetrize and normalize */
         read(cur_vect, buf);
         for (blk=first_ablk[irrep]; blk<=last_ablk[irrep]; blk++) {
            ac = Ia_code[blk];
            bc = Ib_code[blk];
            mat = blocks[blk];
            if (ac == bc) { /* diagonal block */
               for (i=0; i<Ia_size[blk]; i++) {
                  for (j=0; j<i; j++) {
                     mat[j][i] = mat[i][j] * phase;
                     }
                  }
               }
            if (ac > bc) { /* off-diagonal block in lower triangle */
               upper = decode[bc][ac];
               if (upper >= 0) {
                  for (i=0; i<Ia_size[blk]; i++) {
                     for (j=0; j<Ib_size[blk]; j++) {
                        blocks[upper][j][i] = mat[i][j] * phase;
                        }
                     }
                  }
               }
            } /* end loop over blocks */
         write(cur_vect, buf);
         } /* end loop over buffers */
      } /* end icore==2 */

   else if (icore==0) { /* one RAS block at a time */
      for (buf=0; buf<buf_per_vect; buf++) {
         blk = buf2blk[buf];
         ac = Ia_code[blk];
         bc = Ib_code[blk];
         mat = blocks[blk];
         if (ac == bc) { /* diagonal block */
            read(cur_vect, buf); /* just do the ones we need to */
            for (i=0; i<Ia_size[blk]; i++) {
               for (j=0; j<i; j++) {
                  mat[j][i] = mat[i][j] * phase;
                  }
               }
            write(cur_vect, buf);
            }
         } /* end loop over buffers */

      } /* end case icore==0 */


   else {
      printf("(CIvect::symmetrize): Unrecognized icore option\n");
      return;
      }

}

