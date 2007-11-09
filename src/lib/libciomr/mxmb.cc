/*!
** \file mxmb.cc
** \ingroup (CIOMR)
*/
 
#include "includes.h"

extern "C" {

  static void mxmbol(double **a, int ia, int ja, double **b, int ib, int jb, 
                     double **c, int ic, int jc, int nrow, int nlnk, int ncol) {
    abort();
  }
  
/*!
** mxmb: multiplies two rectangular matrices together.  If in=1 (n=a,b,c),
** then normal multiply.  If jn=1 (n=a,b,c) then multiply the transpose
** of matrix n.
** \ingroup (CIOMR)
*/
void mxmb(double **a, int ia, int ja, double **b, int ib, int jb, 
          double **c, int ic, int jc, int nrow, int nlnk, int ncol)
   {
      if (ic == 1) {
         if (ia == 1) {
            if (ib == 1) {
               mmult(a,0,b,0,c,0,nrow,nlnk,ncol,0);
               }
            else {
               if (jb == 1) {
                  mmult(a,0,b,1,c,0,nrow,nlnk,ncol,0);
                  }
               else {
                  mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
                  }
               }
            }
         else {
            if (ja == 1) {
               if (ib == 1) {
                  mmult(a,1,b,0,c,0,nrow,nlnk,ncol,0);
                  }
               else {
                  if (jb == 1) {
                     mmult(a,1,b,1,c,0,nrow,nlnk,ncol,0);
                     }
                  else {
                     mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
                     }
                  }
               }
            else {
               mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
               }
            }
         }
      else {
         if (jc == 1) {
            if (ia == 1) {
               if (ib == 1) {
                  mmult(a,0,b,0,c,1,nrow,nlnk,ncol,0);
                  }
               else {
                  if (jb == 1) {
                     mmult(a,0,b,1,c,1,nrow,nlnk,ncol,0);
                     }
                  else {
                     mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
                     }
                  }
               }
            else {
               if (ja == 1) {
                  if (ib == 1) {
                     mmult(a,1,b,0,c,1,nrow,nlnk,ncol,0);
                     }
                  else {
                     if (jb == 1) {
                        mmult(a,1,b,1,c,1,nrow,nlnk,ncol,0);
                        }
                     else {
                        mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
                        }
                     }
                  }
               else {
                  mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
                  }
               }
            }
         else {
            mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
            }
         }
      }

} /* extern "C" */
