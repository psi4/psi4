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
** no M_s = 0 simplifications here yet: assumes everything available
*/
void olsen_update(CIvect &Cn, CIvect &Cc, CIvect &Sc, CIvect &Hd, 
      double E, double E_est, double *norm, double *ovrlap)
{

   int i,blk,cac,cbc,len;
   double tn,to;

   /* if everything is in core */
   if (Cn.icore == 1 && Cc.icore == 1 && Sc.icore == 1 && Hd.icore == 1) {
      *norm = 0.0;
      *ovrlap = 0.0;
      if (Cn.Ms0) {
         for (blk = 0; blk < Cn.num_blocks; blk++) {
            cac = Cn.Ia_code[blk];
            cbc = Cn.Ib_code[blk];
            if (cac < cbc) continue;
            if (cac != cbc) {
               len = Cn.Ia_size[blk] * Cn.Ib_size[blk];
               H0block_ols_upd(E, E_est, 1.00,
                  Cn.blocks[blk], Cc.blocks[blk], Sc.blocks[blk],
                  Hd.blocks[blk], &tn, &to,   
                  Cn.Ia_code[blk], Cn.Ib_code[blk], len, 1);
               *norm += 2.0 * tn; 
               *ovrlap += 2.0 * to;
               } 
            else {
               len = Cn.Ia_size[blk];
               H0block_diag_ols_upd(E, E_est, 1.000,
                  Cn.blocks[blk], Cc.blocks[blk],
                  Sc.blocks[blk], Hd.blocks[blk], &tn, &to, 
                  Cn.Ia_code[blk], Cn.Ib_code[blk], len);
               *norm += tn;
               *ovrlap += to;
               }
            }
         *norm = sqrt(1.0 / *norm);      
         }
      else {
         for (blk = 0; blk < Cn.num_blocks; blk++) {
            cac = Cn.Ia_code[blk];
            cbc = Cn.Ib_code[blk];
            len = Cn.Ia_size[blk] * Cn.Ib_size[blk];
            H0block_ols_upd(E, E_est, 1.00,
               Cn.blocks[blk], Cc.blocks[blk], Sc.blocks[blk],
               Hd.blocks[blk], &tn, &to,
               Cn.Ia_code[blk], Cn.Ib_code[blk], len, 0);
            *norm += tn;
            *ovrlap += to;
            }
         *norm = sqrt(1.0 / *norm);
         }

      } /* end all in-core option */
   else {
      printf("(olsen_update): Unavailable icore option\n");
      return;
      }

}         

