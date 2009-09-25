/*
** GENCIC.CC
**
** This file contains code for the General CI (GENCI) module
**
** Matt Leininger
** Center for Computational Quantum Chemistry, UGA
** 1996
*/


#define EXTERN

/*** INCLUDES ***/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include "structs.h"
#include "globals.h"
#include "genci.h"

namespace psi { namespace detci {

/*** DEFINES ***/
#define NORM_TOL 1.0E-5
#define INDEX(i,j) ( (i>j) ? (ioff[(i)] + (j)): (ioff[(j)] + (i)) )
#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

/*
** SCHMIDT_ADDOC(): Assume A is a orthogonal matrix.  
** This function Gram-Schmidt
** orthogonalizes a new vector v and adds it to matrix A.  A must contain
** a free row pointer for a new row.  Don't add orthogonalized v' if 
** norm(v') < NORM_TOL.
**
** David Sherrill, Feb 1994
** Matt Leininger, June 1996 Rewrote for out-of-core 
**
** Arguments:
**    buffer4 = buffer to store d or f vectors 
**    buffer5 = buffer to store b vectors 
**    buf_size = size of buffers
**    extra_buf = size of the last extra buffer of a vector
**    num_buf = number of buffers of size buf_size in a vector length N
**    f_index = 
**    N = dimension (length) of b, d, or f vectors 
**
** Returns: 1 if a vector is added to A, 0 otherwise
*/

int schmidt_addoc(double *buffer4, double *buffer5, int buf_size, int extra_buf,
                int num_buf, PSI_FPTR d_index, int N, int L,
                int b_file, int d_file)
{
   double *dotval, normval = 0.0, tval;
   int i, j, I;
   PSI_FPTR b_index = 0, jnk = 0, f_index;

   /* determine array of dot products b.f */
   f_index = d_index; 
   dotval = init_array(L); 
   for (i=0; i<num_buf; i++) { /* loop over num_buf's in b and f */
      wreadw(d_file, (char *) buffer4, sizeof(double)*buf_size, 
             d_index, &d_index); 
      b_index = i*buf_size*sizeof(double); 
      for (j=0; j<L; j++) { /* loop over b vectors */
         wreadw(b_file, (char *) buffer5, sizeof(double)*buf_size, 
                b_index, &jnk); 
         dot_arr(buffer4, buffer5, buf_size, &tval);
         dotval[j] += tval; 
         b_index += N*sizeof(double); 
         } /* end j loop */
      } /* end i loop */
   if (extra_buf != 0) {
   b_index = num_buf*buf_size*sizeof(double); 
   wreadw(d_file, (char *) buffer4, sizeof(double)*extra_buf, 
          d_index, &d_index); 
   for (i=0; i<L; i++) {
      wreadw(b_file, (char *) buffer5, sizeof(double)*extra_buf, 
             b_index, &jnk);
      dot_arr(buffer4, buffer5, extra_buf, &tval);
      dotval[i] += tval; 
      b_index += N*sizeof(double);
      } /* end i loop */
     } /* end extra_buf */

   /* schmidt orthronormalize f to the set of b vectors */
   d_index -= N*sizeof(double); 
   for (i=0; i<num_buf; i++) {
      wreadw(d_file, (char *) buffer4, sizeof(double)*buf_size, 
             d_index, &d_index);
      b_index = i*buf_size*sizeof(double); 
      for (j=0; j<L; j++) { /* loop over b vectors */
         wreadw(b_file, (char *) buffer5, sizeof(double)*buf_size, 
                b_index, &jnk);
         b_index += N*sizeof(double);
         for (I=0; I<buf_size; I++) 
            buffer4[I] -= dotval[j]*buffer5[I]; 
         } /* end j loop */
      dot_arr(buffer4, buffer4, buf_size, &tval); 
      normval += tval; 
      wwritw(d_file, (char *) buffer4, sizeof(double)*buf_size, 
             f_index, &f_index); 
      } /* end i loop */
   if (extra_buf != 0) {
     b_index = num_buf*buf_size*sizeof(double);
     wreadw(d_file, (char *) buffer4, sizeof(double)*extra_buf, 
            d_index, &d_index);
     for (i=0; i<L; i++) {
        wreadw(b_file, (char *) buffer5, sizeof(double)*extra_buf,
               b_index, &jnk);
        b_index += N*sizeof(double);
        for (I=0; I<extra_buf; I++)
           buffer4[I] -= dotval[i]*buffer5[I];
        } /* end i loop */
     dot_arr(buffer4, buffer4, extra_buf, &tval); 
     normval += tval; 
     wwritw(d_file, (char *) buffer4, sizeof(double)*extra_buf, 
            f_index, &f_index); 
     } /* end extra_buf */
   normval = sqrt(normval); 
   /* fprintf(outfile,"\n normval = %15.10f\n", normval); */ 
 
  /* test norm of new b vector */
  free(dotval); 
  if (normval < NORM_TOL)
    return(0); 
  else {
      f_index -= N*sizeof(double); 
      b_index = N*L*sizeof(double); 
      for (i=0; i<num_buf; i++) {
         wreadw(d_file, (char *) buffer4, sizeof(double)*buf_size, 
                f_index, &f_index);
         for (I=0; I<buf_size; I++)  
            buffer4[I] /= normval; 
         wwritw(b_file, (char *) buffer4, sizeof(double)*buf_size, 
                b_index, &b_index); 
         } /* end i loop */
      if (extra_buf != 0) {
        wreadw(d_file, (char *) buffer4, sizeof(double)*extra_buf, 
               f_index, &f_index);
        for (I=0; I<extra_buf; I++)  
           buffer4[I] /= normval; 
        wwritw(b_file, (char *) buffer4, sizeof(double)*extra_buf, 
               b_index, &b_index); 
        } /* end extra_buf */

      return(1); 
      } /* end else */

}


/*
** NORMALIZE.C : Normalize a set of vectors
**
** Assume we're normalizing the ROWS
**
** David Sherrill, Feb 1994
** Matt Leininger, June 1996 modified for out-of-core vectors
*/
void v_normalize(double *A, PSI_FPTR index, int buf_size, 
               int extra_buf, int num_buf, int d_file)
{
   double normval = 0.0, tval;
   register int i, I;
   PSI_FPTR jnk = 0, index2; 

   index2 = index; 
  
   /* divide each row by the square root of its norm */
   for (i=0; i<num_buf; i++) {
      wreadw(d_file, (char *) A, sizeof(double)*buf_size, index, &index); 
      dot_arr(A, A, buf_size, &tval);
      normval += tval;
      } /* end i loop */
   if (extra_buf != 0) {
     wreadw(d_file, (char *) A, sizeof(double)*extra_buf, index, &index); 
     dot_arr(A, A, extra_buf, &tval);
     normval += tval;
     } /* end extra_buf */
   normval = sqrt(normval); 

   index = index2; 
   for (i=0; i<num_buf; i++) {
      wreadw(d_file, (char *) A, sizeof(double)*buf_size, index, &index);
      for (I=0; I<buf_size; I++) 
         A[I] /= normval;  
      wwritw(d_file, (char *) A, sizeof(double)*buf_size, index2, &index2); 
      } /* end i loop */
   if (extra_buf != 0) {
     wreadw(d_file, (char *) A, sizeof(double)*extra_buf, index, &index);
     for (I=0; I<extra_buf; I++)
        A[I] /= normval;
     wwritw(d_file, (char *) A, sizeof(double)*extra_buf, index2, &index2);  
     } /* end extra_buf */

}

/*
** V_SCHMIDT(): Assume A is a orthogonal matrix.  This function Gram-Schmidt
** orthogonalizes a set of given vectors.
** David Sherrill, Feb 1994
** Matt Leininger, June 1996 Rewrote for out-of-core 
**
** Arguments:
**    buffer4 = buffer to store b vectors 
**    buffer5 = buffer to store b vectors 
**    buf_size = size of buffers
**    extra_buf = size of the last extra buffer of a vector
**    num_buf = number of buffers of size buf_size in a vector length N
**    N = dimension (length) of b vectors 
**
*/
double *v_schmidt(double *buffer4, double *buffer5, int buf_size, int extra_buf,
                int num_buf, int N, int L, int b_file)
{
   double normval, tval, *dotval;
   int i, j, k, I, tridim;
   PSI_FPTR x_index = 0, jnk = 0, y_index, xwrit_index;

   /* determine array of dot products b.f */
   x_index = N*sizeof(double); 
   xwrit_index = x_index; 
   tridim = (L*(L+1))/2;
   dotval = init_array(tridim); 
 for (k=1; k<L; k++) { /* loop over b vectors */
   for (i=0; i<num_buf; i++) { /* loop over num_buf's in b vectors */
      wreadw(b_file, (char *) buffer4, sizeof(double)*buf_size, 
             x_index, &x_index); 
      y_index = i*buf_size*sizeof(double); 
      for (j=0; j<=k; j++) { /* loop over b vectors */
         wreadw(b_file, (char *) buffer5, sizeof(double)*buf_size, 
                y_index, &jnk); 
         dot_arr(buffer4, buffer5, buf_size, &tval);
         dotval[INDEX(k,j)] += tval; 
         y_index += N*sizeof(double); 
         } /* end j loop */
      } /* end i loop */
   if (extra_buf != 0) {
   y_index = num_buf*buf_size*sizeof(double); 
   wreadw(b_file, (char *) buffer4, sizeof(double)*extra_buf,
          x_index, &x_index); 
   for (i=0; i<=k; i++) {
      wreadw(b_file, (char *) buffer5, sizeof(double)*extra_buf, 
             y_index, &jnk);
      dot_arr(buffer4, buffer5, extra_buf, &tval);
      dotval[INDEX(i,k)] += tval; 
      y_index += N*sizeof(double);
      } /* end i loop */
     } /* end extra_buf */

   /* schmidt orthronormalize the set of b vectors */
   x_index -= N*sizeof(double); 
   normval = 0;
   for (i=0; i<num_buf; i++) {
      wreadw(b_file, (char *) buffer4, sizeof(double)*buf_size, 
             x_index, &x_index);
      y_index = i*buf_size*sizeof(double); 
      for (j=0; j<k; j++) { /* loop over b vectors */
         wreadw(b_file, (char *) buffer5, sizeof(double)*buf_size, 
                y_index, &jnk);
         y_index += N*sizeof(double);
         for (I=0; I<buf_size; I++) 
            buffer4[I] -= dotval[INDEX(j,k)]*buffer5[I]; 
         } /* end j loop */
      dot_arr(buffer4, buffer4, buf_size, &tval);
      normval += tval;
      wwritw(b_file, (char *) buffer4, sizeof(double)*buf_size, 
             xwrit_index, &xwrit_index); 
      } /* end i loop */
   if (extra_buf != 0) {
     y_index = num_buf*buf_size*sizeof(double);
     wreadw(b_file, (char *) buffer4, sizeof(double)*extra_buf, 
            x_index, &x_index); 
     for (i=0; i<k; i++) {
        wreadw(b_file, (char *) buffer5, sizeof(double)*extra_buf, 
               y_index, &jnk);
        y_index += N*sizeof(double);
        for (I=0; I<extra_buf; I++)
           buffer4[I] -= dotval[INDEX(i,k)]*buffer5[I];
        } /* end i loop */
     dot_arr(buffer4, buffer4, extra_buf, &tval);
     normval += tval;
     wwritw(b_file, (char *) buffer4, sizeof(double)*extra_buf, 
            xwrit_index, &xwrit_index); 
     } /* end extra_buf */

   /* Normalize vector */
   normval = sqrt(normval);
   x_index -= N*sizeof(double);
   for (i=0; i<num_buf; i++) {
      wreadw(b_file, (char *) buffer4, sizeof(double)*buf_size, 
             x_index, &jnk);
      for (I=0; I<buf_size; I++) 
         buffer4[I] /= normval;
      wwritw(b_file, (char *) buffer4, sizeof(double)*buf_size, 
             x_index, &x_index);
      } /* end i loop */
   if (extra_buf != 0) {
     wreadw(b_file, (char *) buffer4, sizeof(double)*extra_buf,
            x_index, &jnk);
     for (I=0; I<extra_buf; I++) 
         buffer4[I] /= normval;
      wwritw(b_file, (char *) buffer4, sizeof(double)*extra_buf, 
             x_index, &x_index);
     } /* end extra_buf */
   } /* end k loop */
 return(dotval);
}

void det2strings(BIGINT det, int *alp_code, int *alp_idx,
                 int *bet_code, int *bet_idx)
{
   int i;

   /* determine the CI block we're in */
   for (i=0; i<CIblks.num_blocks-1; i++) {
      if (CIblks.offset[i+1] > det) break;
      }
   *alp_code = CIblks.Ia_code[i];
   *bet_code = CIblks.Ib_code[i];

   *alp_idx = (det - CIblks.offset[i]) / CIblks.Ib_size[i];
   *bet_idx = (det - CIblks.offset[i]) % CIblks.Ib_size[i];

}

BIGINT strings2det(int alp_code, int alp_idx, int bet_code, int bet_idx)
{
   int blknum;
   BIGINT addr;

   blknum = CIblks.decode[alp_code][bet_code];
   if (blknum == -1) {
     fprintf(outfile, "CIvect::strings2det failed --- invalid block\n");
     exit(1);
   }

   addr = CIblks.offset[blknum];
   addr += alp_idx * CIblks.Ib_size[blknum] + bet_idx;

   return(addr);
}


/*
** unit_guess - uses a unit guess for the first b vector in the
**              davidson-liu algorithm
**
** 
*/
void unit_guess(int alp_code, int alp_idx, int bet_code, int bet_idx,
                int switch_buf3, double *buffer, int buf_size,
                int num_buf, int extra_buf, PSI_FPTR b_file,
                PSI_FPTR b_writ, int M, int N)
{
 int i;
 register int j;
  

   /* Form initial guess b vector */
   i = strings2det(CalcInfo.ref_alp_list, CalcInfo.ref_alp_rel,
                   CalcInfo.ref_bet_list, CalcInfo.ref_bet_rel);
   if (!switch_buf3) {
     buffer[i] = 1.0;
     wwritw(b_file, (char *) buffer, sizeof(double)*buf_size,
            b_writ, &b_writ);
     buffer[i] = 0.0;
     v_normalize(buffer, (b_writ-M*N*sizeof(double)), buf_size, 0,
                 1, b_file);
     }
   else if (switch_buf3) {
     for (j=0; j<num_buf; j++) {
        if ((i >= (j*buf_size)) && (i < ((j+1)*buf_size)))
          buffer[i-j*buf_size] = 1.0;
        wwritw(b_file, (char *) buffer, sizeof(double)*buf_size,
               b_writ, &b_writ);
        if ((i >= (j*buf_size)) && (i < ((j+1)*buf_size)))
          buffer[i-j*buf_size] = 0.0;
        }
     if (extra_buf != 0) {
       if ((i >= num_buf*buf_size) && (i < (extra_buf+num_buf*buf_size)))
         buffer[i-num_buf*buf_size] = 1.0;
       wwritw(b_file, (char *) buffer, sizeof(double)*extra_buf,
              b_writ, &b_writ);
       if ((i > num_buf*buf_size) && (i < (extra_buf+num_buf*buf_size)))
         buffer[i-num_buf*buf_size] = 0.0;
       }
     v_normalize(buffer, (b_writ-M*N*sizeof(double)), buf_size, extra_buf,
                 num_buf, b_file);
     }

   /* b_writ -= N*sizeof(double);
        wreadw(b_file, (char *) buffer, sizeof(double)*buf_size,
               b_writ, &b_writ);
        for (I=0; I<buf_size; I++)
           if (buffer[I] != 0.0) fprintf(outfile,"b[0][%d] = %lf\n\n",
               I, buffer[I]);
        printf("Done reading in b trial vector for check.\n"); */

}


/*
** max_element  - determines the maximum value in a list of diagonal
**               matrix elements.
** buffer       - array of matrix elements
** num_elements - number of elements in array buffer
** max          - pointer to max value
** max_num      - the element number to the maximum value
*/
void max_element(double *buffer, int num_elements, double *max, int *max_num)
{
 register int i;

 (*max) = buffer[0];
 (*max_num) = 0; 
 for (i=1; i<num_elements; i++) {
    if (buffer[i] > (*max)) {
      (*max) = buffer[i];
      (*max_num) = i;
      }
    } 

}


/*
** min_element  - determines the minimum value in a list of diagonal
**               matrix elements.
** buffer       - array of matrix elements
** num_elements - number of elements in array buffer
** min          - pointer to max value
** min_num      - the element number to the minimum value
*/
void min_element(double *buffer, int num_elements, double *min, int *min_num)
{
 register int i;

 (*min) = buffer[0];
 (*min_num) = 0;
 for (i=1; i<num_elements; i++) {
    if (buffer[i] < (*min)) {
      (*min) = buffer[i];
      (*min_num) = i;
      }
    }

}


/*
**  read_c() - reads in the CI vector for a restart of the calculation
**
*/
void read_c(int switch_buf3, double *buffer, int buf_size, int num_buf,
            int extra_buf, int b_file, PSI_FPTR b_writ, int c_file, 
            PSI_FPTR c_index)
{
 register int i;

 if (!switch_buf3) {
   wreadw(c_file, (char *) buffer, sizeof(double)*buf_size,
          c_index, &c_index);
   wwritw(b_file, (char *) buffer, sizeof(double)*buf_size,
          b_writ, &b_writ);
   }
 else if (switch_buf3) {
   for (i=0; i<num_buf; i++) {
      wreadw(c_file, (char *) buffer, sizeof(double)*buf_size,
          c_index, &c_index);
      wwritw(b_file, (char *) buffer, sizeof(double)*buf_size,
          b_writ, &b_writ);
      }
   if (extra_buf != 0) {
      wreadw(c_file, (char *) buffer, sizeof(double)*extra_buf,
          c_index, &c_index);
      wwritw(b_file, (char *) buffer, sizeof(double)*extra_buf,
          b_writ, &b_writ);
     }
   }

}

}} // namespace psi::detci

