/*
** H0block_ols_upd()
**
** Parameters:
**    Ms0  =  1 if lower triangle blocks of C are used to update upper 
**            triangle parts of H0block.
*/
void H0block_ols_upd(double E, double E_est, double phase,  
      double **Cn, double **Cc, double **S, double **Hd,
      double *norm, double *ovrlap, int ac, int bc,
      int len, int Ms0)
{
   int i;
   int al,bl;
   double tval,tval2,c;
   double *ccptr, *sptr, *cnptr, *hdptr;
   double nx = 0.0, ox = 0.0;

   ccptr = Cc[0];
   sptr = S[0];
   cnptr = Cn[0];
   hdptr = Hd[0];

   for (i=0; i<len; i++) {
      tval = hdptr[i] - E;
      if (fabs(tval) < HD_MIN) tval = HD_MIN ; /* prevent /0 */
      tval = 1.0 / tval;
      tval2 = ccptr[i] + tval * (E_est * ccptr[i] - sptr[i]);
      cnptr[i] = tval2;
      nx += tval2 * tval2;
      ox += tval2 * ccptr[i];
      }         


   for (i=0; i<H0block.size; i++) {
      if (ac != H0block.alplist[i] || bc != H0block.betlist[i]) continue;
      al = H0block.alpidx[i];
      bl = H0block.betidx[i];
      tval = Hd[al][bl] - E;
      if (fabs(tval) < HD_MIN) tval = HD_MIN ; /* prevent /0 */
      tval = 1.0 / tval;
      c = Cc[al][bl];
      /* tval2 = c + tval * (E_est * c - S[al][bl]); */
      tval2 = Cn[al][bl];
      nx -= tval2 * tval2;
      ox -= tval2 * c;
      tval = c + E_est * H0block.c0bp[i];
      tval -= H0block.s0bp[i];
      Cn[al][bl] = tval;
      H0block.c0b[i] = tval; /* this should gather all of c0b. Norm later */
      nx += tval * tval;
      ox += tval * c;
      }

   /* get upper triangle of H0block */
   for (i=0; i<H0block.size; i++) {
      if (ac != H0block.betlist[i] || bc != H0block.alplist[i]) continue;
      al = H0block.alpidx[i];
      bl = H0block.betidx[i];
      H0block.c0b[i] = phase * Cn[bl][al]; 
      }

   *norm = nx;
   *ovrlap = ox;
}


/* call for Ms=0 cases.  If not Ms=0, then call H0block_ols_upd */
void H0block_diag_ols_upd(double E, double E_est, double phase, 
      double **Cn, double **Cc, double **S, double **Hd,
      double *norm, double *ovrlap, int ac, int bc, int len)
{
   int i,j;
   int al,bl;
   double tval,tval2,c;
   double nx = 0.0, ox = 0.0;

   for (i=0; i<len; i++) {
      tval = Hd[i][i] - E;
      if (fabs(tval) < HD_MIN) tval = HD_MIN ; /* prevent /0 */
      tval = 1.0 / tval;
      tval2 = Cc[i][i] + tval * (E_est * Cc[i][i] - S[i][i]);
      Cn[i][i] = tval2;
      nx += tval2 * tval2;
      ox += tval2 * Cc[i][i];
      
      for (j=0; j<i; j++) {
         tval = Hd[i][j] - E;
         if (fabs(tval) < HD_MIN) tval = HD_MIN ; /* prevent /0 */
         tval = 1.0 / tval;
         tval2 = Cc[i][j] + tval * (E_est * Cc[i][j] - S[i][j]);
         Cn[i][j] = tval2;
         nx += 2.0 * tval2 * tval2;
         ox += 2.0 * tval2 * Cc[i][j];
         }         
       }


   /* the code below assumes Ms=0 and that h0block members come in pairs
    * if off-diagonal
    */
   for (i=0; i<H0block.size; i++) {
      if (ac != H0block.alplist[i] || bc != H0block.betlist[i]) continue;
      al = H0block.alpidx[i];
      bl = H0block.betidx[i];
      c = Cc[al][bl];
      if (al > bl) {
         tval2 = Cn[al][bl];
         nx -= 2.0 * tval2 * tval2;
         ox -= 2.0 * tval2 * c;         
         tval = c + E_est * H0block.c0bp[i];
         tval -= H0block.s0bp[i];
         Cn[al][bl] = tval;
         nx += 2.0 * tval * tval;
         ox += 2.0 * tval * c;
         H0block.c0b[i] = tval; /* this should gather all of c0b. Norm later */
         }
      else if (al == bl) {
         tval2 = Cn[al][bl];
         nx -= tval2 * tval2;
         ox -= tval2 * c;         
         tval = c + E_est * H0block.c0bp[i];
         tval -= H0block.s0bp[i];
         Cn[al][bl] = tval;
         nx += tval * tval;
         ox += tval * c;
         H0block.c0b[i] = tval; /* this should gather all of c0b. Norm later */
         }
      }

   for (i=0; i<H0block.size; i++) {
      if (ac != H0block.alplist[i] || bc != H0block.betlist[i]) continue;
      al = H0block.alpidx[i];
      bl = H0block.betidx[i];
      if (al < bl) {
         H0block.c0b[i] = phase * Cn[bl][al];
         }
      }

   *norm = nx;
   *ovrlap = ox;
}

