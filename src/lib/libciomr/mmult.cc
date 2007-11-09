/*!
   \file mmult.cc
   \ingroup (CIOMR)
*/

#include <math.h>
#include <libciomr/libciomr.h>

extern "C" {

static int keep_nr=0;
static int keep_nl=0;
static int keep_nc=0;
static double **aa,**bb;

/*!
**                                                             
** mmult():
** a reasonably fast matrix multiply (at least on the DEC3100) 
** written by ETS                                              
**                                                             
** AF,BF,and CF are fortran arrays                             
**                                                             
** ta,tb and tc indicate whether the corresponding arrays are  
**              to be converted to their transpose             
**                                                             
** nr,nl,nc are the number of rows,links,and columns in the    
**          final matrices to be multiplied together           
**          if ta=0 AF should have the dimensions nr x nl      
**          if ta=1 AF should have the dimensions nl x nr      
**          if tb=0 BF should have the dimensions nl x nc      
**          if tb=1 BF should have the dimensions nc x nl      
**          if tc=0 CF should have the dimensions nr x nc      
**          if tc=1 CF should have the dimensions nc x nr      
**                                                             
** add is 1 if this matrix is to be added to the one passed    
**        in as CF, 0 otherwise                                
**
** \ingroup (CIOMR)
*/
void mmult(double **AF, int ta, double **BF, int tb, double **CF, int tc,
	   int nr, int nl, int nc, int add)
{
   int odd_nr,odd_nc,odd_nl;
   int i,j,k,ij;
   double t00,t01,t10,t11;
   double **a,**b;
   double *att,*bt;
   double *at1,*bt1;

   if(!aa) {
      aa = (double **) init_matrix(nr,nl);
      bb = (double **) init_matrix(nc,nl);
      keep_nr = nr;
      keep_nl = nl;
      keep_nc = nc;
      }

   if(nl > keep_nl) {
      free_matrix(aa,keep_nr);
      free_matrix(bb,keep_nc);
      keep_nl = nl;
      keep_nr = (nr > keep_nr) ? nr : keep_nr;
      keep_nc = (nc > keep_nc) ? nc : keep_nc;
      aa = (double **) init_matrix(keep_nr,keep_nl);
      bb = (double **) init_matrix(keep_nc,keep_nl);
      }
   if(nr > keep_nr) {
      free_matrix(aa,keep_nr);
      keep_nr = nr;
      aa = (double **) init_matrix(keep_nr,keep_nl);
      }
   if(nc > keep_nc) {
      free_matrix(bb,keep_nc);
      keep_nc = nc;
      bb = (double **) init_matrix(keep_nc,keep_nl);
      }

   odd_nr = (nr)%2;
   odd_nc = (nc)%2;
   odd_nl = (nl)%2;

   a=aa;
   if(ta)
      for(i=0; i < nr ; i++)
         for(j=0; j < nl ; j++)
            a[i][j] = AF[j][i];
   else
      a=AF;

   b=bb;
   if(tb)
      b=BF;
   else
      for(i=0; i < nc ; i++)
         for(j=0; j < nl ; j++)
            b[i][j] = BF[j][i];
      
   for(j=0; j < nc-1 ; j+=2) {
      for(i=0; i < nr-1 ; i+=2) {
         att=a[i]; bt=b[j];
         at1=a[i+1]; bt1=b[j+1];
         if(add) {
            if(tc) {
               t00 = CF[j][i];
               t01 = CF[j+1][i];
               t10 = CF[j][i+1];
               t11 = CF[j+1][i+1];
               }
            else {
               t00 = CF[i][j];
               t01 = CF[i][j+1];
               t10 = CF[i+1][j];
               t11 = CF[i+1][j+1];
               }
            }
         else t00=t01=t10=t11=0.0;
         for(k=nl; k ; k--,att++,bt++,at1++,bt1++) {
            t00 += *att * *bt;
            t01 += *att * *bt1;
            t10 += *at1 * *bt;
            t11 += *at1 * *bt1;
            }
         if(tc) {
            CF[j][i]=t00;
            CF[j+1][i]=t01;
            CF[j][i+1]=t10;
            CF[j+1][i+1]=t11;
            }
         else {
            CF[i][j]=t00;
            CF[i][j+1]=t01;
            CF[i+1][j]=t10;
            CF[i+1][j+1]=t11;
            }
         }
      if(odd_nr) {
         att=a[i]; bt=b[j];
         bt1=b[j+1];
         if(add) {
            if(tc) {
               t00 = CF[j][i];
               t01 = CF[j+1][i];
               }
            else {
               t00 = CF[i][j];
               t01 = CF[i][j+1];
               }
            }
         else t00=t01=0.0;
         for(k= nl; k ; k--,att++,bt++,bt1++) {
            t00 += *att * *bt;
            t01 += *att * *bt1;
            }
         if(tc) {
            CF[j][i]=t00;
            CF[j+1][i]=t01;
            }
         else {
            CF[i][j]=t00;
            CF[i][j+1]=t01;
            }
         }
      }
   if(odd_nc) {
      for(i=0; i < nr-1 ; i+=2) {
         att=a[i]; bt=b[j];
         at1=a[i+1];
         if(add) {
            if(tc) {
               t00 = CF[j][i];
               t10 = CF[j][i+1];
               }
            else {
               t00 = CF[i][j];
               t10 = CF[i+1][j];
               }
            }
         else t00=t10=0.0;
         for(k= nl; k ; k--,att++,bt++,at1++) {
            t00 += *att * *bt;
            t10 += *at1 * *bt;
            }
         if(tc) {
            CF[j][i]=t00;
            CF[j][i+1]=t10;
            }
         else {
            CF[i][j]=t00;
            CF[i+1][j]=t10;
            }
         }
      if(odd_nr) {
         att=a[i]; bt=b[j];
         if(add)
            t00 = (tc) ? CF[j][i] : CF[i][j];
         else t00=0.0;
         for(k=nl; k ; k--,att++,bt++)
            t00 += *att * *bt;
         if(tc) CF[j][i]=t00;
         else CF[i][j]=t00;
         }
      }
   }

} /* extern "C" */
