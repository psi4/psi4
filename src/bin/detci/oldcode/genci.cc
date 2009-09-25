/* 

   GENCI.CC
  
   Part of the General Configuration Interaction (GENCI) package

   Matt Leininger
   Center for Computational Quantum Chemistry
   1996

   This code uses the string lists from the main DETCI code to explicitly
   compute H and its action on c to produce sigma, using the Slater
   Determinant class.  The idea was that completely arbitrary determinants
   could be specified (through an interface that was not ever implemented).
   Since the code as it was only used the same (RAS-selected) determinants
   as the main code, there was no extra functionality (other than perhaps
   for debugging purposes) and of course this approach is much, much slower.
   When the old wreadw/wwritw I/O calls were obsoleted in favor of the
   new libpsio library, it was not deemed useful to translate this genci
   code to the new library, as it did not add functionality.  However, it
   is retained here as a possible starting point for future work in general 
   CI expansions. --- C. David Sherrill, January 2008
  
*/ 

/*** DEFINES ***/
#define EXTERN

/*** INCLUDES ***/
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <libipv1/ip_lib.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include "structs.h"
#include "globals.h"
#include "genci.h"
#include "slaterd.h"

namespace psi { namespace detci {

#define MAX_B_ROWS 200
#define MIN_F_DENOM 1.0E-3
#define INDEX(i,j) ( (i>j) ? (ioff[(i)] + (j)): (ioff[(j)] + (i)) )

extern struct stringwr **alplist;  
extern struct stringwr **betlist; 

/*** FUNCTION PROTOTYPES ***/
void sem_genci(struct stringwr **alplist, struct stringwr **betlist,
   int N, int M, double *evals, double conv_e, double conv_rms, 
   int maxiter, double nucrep, double efzc, int maxnvect, 
   FILE *outfile, int print_lvl, int h0blocksize, int guess_vector,
   int restart, int restart_vecs);
void calc_h0_blk(double *h0, int *detref, double *detval, int h0blocksize, 
   int tridim, struct stringwr **alplist, struct stringwr **betlist);
void gen_mat(unsigned long i, unsigned long j, int buf_size, double *A,
   struct stringwr **alplist, struct stringwr **betlist);
void gendi_a(unsigned long i, unsigned long buf_size, double *A,
   struct stringwr **alplist, struct stringwr **betlist);
PSI_FPTR sigma(struct stringwr **alplist, struct stringwr **betlist,
   double **buffer, int L, int switch_buf3, int buf_size, int num_buf, 
   int num_new_vec, int extra_buf3, int b_file, int bA_file, 
   int buf_val3, int N, int extra_buf, PSI_FPTR byte, int print_lvl);
void h0_guess(int alp_code, int alp_idx, int bet_code, int bet_idx,
   int switch_buf3, double *buffer, int buf_size,
   int num_buf, int extra_buf, int b_file,
   PSI_FPTR b_writ, int M, int N, int h0blocksize,
   double nucrep, double efzc, struct stringwr **alplist,
   struct stringwr **betlist);


/*** GLOBAL VARIABLES THIS MODULE ***/
int print_lvl = 1; /* print flag */


void diag_h_genci(struct stringwr **alplist, struct stringwr **betlist)
{
  int nroots, i, j, size ;
  double conv_rms, conv_e, *evals, **evecs, nucrep, efzc ;
  int *tptr;
  double *cbuf;
  int restart, restart_vecs;

  nroots = Parameters.num_roots ;
  conv_rms = pow(10.0, -(Parameters.convergence));
  conv_e = pow(10.0, -(Parameters.energy_convergence));
  size = CIblks.vectlen;
  if (Parameters.nprint > size) Parameters.nprint = size;
  nucrep = CalcInfo.enuc;
  efzc = CalcInfo.efzc;
  restart = Parameters.restart;
  restart_vecs = Parameters.restart_vecs;
   
  /* Davidson-Liu SEM for GENCI */ 
  if (Parameters.diag_method == 3) {
     
    if (Parameters.print_lvl) {
      fprintf(outfile, 
        "\nFind the roots by the Simultaneous Expansion Method \n");
      fprintf(outfile, "Energy convergence = %3g\n", conv_e);
      fprintf(outfile, "RMS CI vector convergence = %3g\n\n", conv_rms);
      fflush(outfile);
    }
     
    evals = init_array(nroots);
    sem_genci(alplist, betlist, size, nroots, evals, conv_e, 
              conv_rms, Parameters.maxiter, nucrep, efzc, 
              Parameters.maxnvect, outfile, Parameters.print_lvl,
              Parameters.h0blocksize, Parameters.guess_vector,
              restart, restart_vecs);
     
  } /* end the Davidson-Liu section for GENCI */
   
  else { 
    fprintf(outfile, "\n Diagonalization Method not available for GENCI\n"); 
    return; 
  }
   
  /* write the CI energy to file30: later fix this to loop over roots */
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_etot(evals[Parameters.root]); /* need nucrep or efzc also? */
  chkpt_close();

}



/*
** sem_genci
**
** Based on the Simultaneous Expansion Method (SEM), a variation on the
** Davidson algorithm, by Bowen Liu.  Adapted for the GENCI module.
**
** Matt Leininger
** Center for Computational Quantum Chemistry, UGA
** June 1996
**
** Arguments:
**   N        =  size of matrix
**   M        =  number of eigenvalues to solve for
**   evals    =  array for eigenvalues
**   conv_e   =  convergence tolerance.  The lowest energy eigenvalue must
**                 be converged to within this range.  It is interesting
**                 that higher roots _may_ converge faster.
**   conv_rms =  the required tolerance for convergence of the CI correction
**                 vector
**   maxiter  =  max number of iterations allowed
**   nucrep   =  nucrep to add to eigenvalues in printing (e.g. enuc)
**   efzv     =  frozen core energy to add to eigenvalues in printing
**   vu       =  pointer to int to hold how many expansion vectors used
**
** Returns: none
*/
void sem_genci(struct stringwr **alplist, struct stringwr **betlist,
               int N, int M, double *evals, double conv_e, double conv_rms, 
               int maxiter, double nucrep, double efzc, int maxnvect, 
               FILE *outfile, int print_lvl, int h0blocksize, int guess_vector,
               int restart, int restart_vecs)
{
   double *G, *Gprime, *lambda;   /* arrays for G, Gprime, and lambda */ 
   double **alpha;                /* matrices for eigenvecs, */   
   double *dotval;                /* array for dot products of b vecs */
   double tval = 0.0;             /* value of a dot product */ 
   double **buffer;               /* matrix of buffers */ 
   double *buffer4;               /* tmp buffer */ 
   double *buffer5;               /* tmp buffer */ 
   double *lastroot;              /* stores previous eigenvalue array */
   double tmp = 0.0;              /* tmp component of a dot product */ 
   double tmp_val = 0.0;          /* used when forming G */
   double memory;                 /* amount of memory to be used in MB */ 
   register int i = 0, j;         /* various loop variables */ 
   register int k, I, m, l;
   int imdet;                     /* number of important determinants */
   int L;                         /* number of total b vectors */
   int converged = 0, iter = 1;   /* convergence flag and iter counter */
   int buf_size = 0;              /* the size of a buffer */ 
   int rownum = 0;                /* row number variable for H matrix */ 
   int num_new_vec = 0;           /* number of new b vectors in an iter */ 
   int num_buf, core_memory;      /* number of buffers and memory(bytes) */ 
   int extra_buf;                 /* number of extra buffers necessary */
   int A_col_index = 0;           /* column variable */
   int tridim, errcod;            /* dimension of low tri and error code */
   int old_vec = 0;               /* represents the last old b vector */ 
   int b_file = 96;               /* file number for b vectors */
   int bA_file = 97;              /* file number for sigma vectors */
   int d_file = 98;               /* file number for d corr. vectors */
   int c_file = 77;               /* file number for final ci vector */
   int buf_val2 = 0;              /* num bufs if large_core for 2 buf met */
   int buf_val3 = 0;              /* num bufs if large_core for 3 buf met */
   int switch_buf2 = 1;           /* if small core then switch buf sizes */  
   int switch_buf3 = 1;           /* if small core then switch buf sizes */
   PSI_FPTR bT_index = 0;/* index for b transpose during bA . bT */
   PSI_FPTR buf_size2;   /* buf size for two buf method */
   PSI_FPTR extra_buf2;  /* extra buf size for two buf method */
   PSI_FPTR num_buf2;    /* num bufs for two buf method */ 
   PSI_FPTR buf_size3 =0;/* buf size for three buf method */
   PSI_FPTR extra_buf3;  /* extra buf size for three buf method */
   PSI_FPTR num_buf3;    /* num bufs for three buf method */ 
   PSI_FPTR index = 0;   /* index for extra_buf of d corr. vector */
   PSI_FPTR d_index = 0;          /* index for d correction vector */
   PSI_FPTR dummy = 0;            /* dummy index */ 
   PSI_FPTR byte = 0;             /* position for writing bA to disk */
   PSI_FPTR jnk = 0;              /* another dummy index */
   PSI_FPTR c_index = 0;          /* position fow writting of CI vector */
   PSI_FPTR b_index2 = 0;         /* b index for convergence check */
   PSI_FPTR b_writ = 0;           /* b index for initial guess */
   PSI_FPTR bA_index = 0;/* index for bA during bA . bT */

   tridim = (maxnvect*(maxnvect+1))/2;

   memory = 16.0005;                 /* default memory */
   errcod = ip_data("MEMORY","%lf", &memory,0);

   L = M;                        /* set number of b vec to num roots */
   num_new_vec = L;              /* number of new b vectors */ 

   core_memory = (int) (((memory)*1E6)/sizeof(double));
   if (core_memory > (3*N)) {
     buf_size3 = N;
     num_buf3 = (core_memory/N)-2;
     if (num_buf3 > N) num_buf3 = N-1;
     extra_buf3 = (N%(num_buf3));
     switch_buf3 = 0;
     fprintf(outfile,"Using himem\n");
     fflush(outfile);
     }
   if (core_memory > (2*N)) {
     buf_size2 = N;
     /* buf_size2 = core_memory/2; */
     num_buf2 = (core_memory/N)-1;
     /* num_buf2 = (N/buf_size2 - (N%buf_size2)/buf_size2); */
     if (num_buf2 > N) num_buf2 = N; 
     extra_buf2 = (N%(num_buf2));
     switch_buf2 = 0; 
     fflush(outfile);
     if (!buf_size3) {
       fprintf(outfile,"Using core_memory > (2*N)\n");
       fflush(outfile);
       buf_size3 = core_memory/3;
       /* num_buf3 = N/buf_size3; */
       num_buf3 = (N/buf_size3 - (N%buf_size3)/buf_size3); 
       extra_buf3 = N - num_buf3*buf_size3;
       }
     /* buffer = block_matrix(num_buf2+1, buf_size2); */
     buffer = block_matrix(num_buf2+1, buf_size2);
     }
   else {
     buf_size3 = core_memory/3;
     extra_buf3 = N%buf_size3;
     num_buf3 = (N/buf_size3 - (N%buf_size3)/buf_size3);
     buf_val3 = 1;
     buf_size2 = core_memory/2;
     extra_buf2 = N%buf_size2;
     num_buf2 = (N/buf_size2 - (N%buf_size2)/buf_size2);
     buf_val2 = 1;
     buffer = block_matrix(2, buf_size2);
     fprintf(outfile,"Using lomem\n");
     fflush(outfile);
     }

   buf_size = buf_size3;
   extra_buf = extra_buf3;
   num_buf = num_buf3;

   fprintf(outfile,"\n");
   fprintf(outfile,"memory       = %lf\n", memory);
   fprintf(outfile,"num_roots    = %d\n", M);
   fprintf(outfile,"maxiter      = %d\n", maxiter);
   fprintf(outfile,"core_memeory = %d\n", core_memory);
   fprintf(outfile,"buf_size2    = %d\n", buf_size2);
   fprintf(outfile,"num_buf2     = %d\n", num_buf2);
   fprintf(outfile,"extra_buf2   = %d\n", extra_buf2);
   fprintf(outfile,"switch_buf2  = %d\n", switch_buf2);
   fprintf(outfile,"buf_size3    = %d\n", buf_size3);
   fprintf(outfile,"num_buf3     = %d\n", num_buf3);
   fprintf(outfile,"extra_buf3   = %d\n", extra_buf3);
   fprintf(outfile,"switch_buf3  = %d\n", switch_buf3);
   fprintf(outfile,"nucrep       = %lf\n", nucrep);
   fprintf(outfile,"efzc         = %lf\n", efzc);
   fprintf(outfile,"\n");
   fprintf(outfile,"\n");
   fflush(outfile);

   if (switch_buf3) {
     buffer4 = buffer[1];
     buffer5 = buffer[2];
     buffer[1] = buffer[0]+buf_size3;
     buffer[2] = buffer[0]+2*buf_size3;
     }

   /* allocate other arrays with ~fixed dimensions during iteration */
   lambda = init_array(maxnvect);
   alpha = init_matrix(maxnvect, maxnvect);
   lastroot = init_array(maxnvect); 
   evals = init_array(maxnvect); 
   G = init_array(tridim);  
   /* Gprime = init_array(tridim); */

   /* open up storage files for b vectors and product of bA */
     rfile(b_file);  /* file store b vectors */
     rfile(bA_file); /* file stores the bA matrix */
     rfile(d_file);  /* file store the d part of the correction vectors */
     rfile(c_file);  /* file to write CI vector to */

 fprintf(outfile,"Determining initial guess vector(s)....");
 if (restart == 1) 
   read_c(switch_buf3, buffer[0], buf_size, num_buf, extra_buf, b_file,
          b_writ, c_file, c_index);
 else if (guess_vector == 0 && M==1) 
   unit_guess(CalcInfo.ref_alp_list, CalcInfo.ref_alp_rel,
              CalcInfo.ref_bet_list, CalcInfo.ref_bet_rel,
              switch_buf3, buffer[0], buf_size, num_buf,
              extra_buf, b_file, b_writ, M, N);
 else if (guess_vector == 1 || M!=1)
   h0_guess(CalcInfo.ref_alp_list, CalcInfo.ref_alp_rel,
            CalcInfo.ref_bet_list, CalcInfo.ref_bet_rel,
            switch_buf3, buffer[0], buf_size, num_buf,
            extra_buf, b_file, b_writ, M, N, h0blocksize,
            nucrep, efzc, alplist, betlist);
 
 while (!converged && iter-1 <= maxiter) {

   buf_size = buf_size3;
   extra_buf = extra_buf3;
   num_buf = num_buf3;

   zero_arr(buffer[0], buf_size); 
   zero_arr(buffer[1], buf_size); 

   byte = sigma(alplist, betlist, buffer, L, switch_buf3, buf_size, 
                num_buf, num_new_vec, extra_buf3, b_file, bA_file, 
                buf_val3, N, extra_buf, byte, print_lvl); 

   /* form G matrix - bA.b*/
   buf_size = buf_size2;
   extra_buf = extra_buf2;
   num_buf = num_buf2;
   if (switch_buf3) {
     buffer[1] = buffer4;
     buffer[2] = buffer5;
     } 
   if (!switch_buf3) {
     extra_buf = 0;
     num_buf = 1;
     }
   bA_index = 0 ; 
   old_vec = (L-num_new_vec)*N*sizeof(double);
   i = 0;
   while (bA_index < (L*N*sizeof(double))) {
        if ((bA_index < old_vec) || ((bA_index == 0) && (iter == 1))){
          /* mult by all new b vectors */
             for (j=0; j<num_buf; j++) { 
                /* loop over # buffers in bA and bT */
                bT_index = j*buf_size*sizeof(double)+old_vec;
                           /* index for bT during bA . bT */
                wreadw(bA_file, (char *) buffer[0], sizeof(double)*buf_size, 
                       bA_index, &bA_index);
                for (k=(L-num_new_vec); k<L; k++) { /* loop over b vectors */
                   wreadw(b_file, (char *) buffer[1], sizeof(double)*buf_size, 
                          bT_index, &jnk);
                   dot_arr(buffer[0], buffer[1], buf_size, &tmp);
                   G[INDEX(i,k)] += tmp;
                   bT_index += N*sizeof(double);
                   } /* end k loop */
                } /* end j loop */
              if (extra_buf != 0) {
                bT_index = num_buf*buf_size*sizeof(double)+old_vec;
                wreadw(bA_file,(char *) buffer[0], sizeof(double)*extra_buf, 
                       bA_index, &bA_index);
                for (k=(L-num_new_vec); k<L; k++) { /* # columns of bA */
                   wreadw(b_file, (char *) buffer[1], sizeof(double)*extra_buf, 
                          bT_index, &jnk);
                   dot_arr(buffer[0], buffer[1], extra_buf, &tmp);
                   G[INDEX(i,k)] += tmp;
                   bT_index += N*sizeof(double);
                   } /* end k loop */
                } /* end if extra_buf */
          } /* end if */
         else {
             /* bA(i) is new and is mult by all b(i) vectors where i->L */
                   dummy = bA_index; 
                   for (j=0; j<num_buf; j++) { 
                      /* loop over # buffers in bA and bT */
                      bT_index = j*buf_size*sizeof(double)+dummy;
                                 /* index for bT during bA . bT */
                      wreadw(bA_file, (char *) buffer[0], 
                             sizeof(double)*buf_size, bA_index, &bA_index);
                      for (k=i; k<L; k++) { /* loop over b vectors */
                         wreadw(b_file, (char *) buffer[1], 
                                sizeof(double)*buf_size, bT_index, &jnk);
                         dot_arr(buffer[0], buffer[1], buf_size, &tmp);
                         G[INDEX(i,k)] += tmp;
                         bT_index += N*sizeof(double);
                         } /* end k loop */
                      } /* end j loop */
                    if (extra_buf != 0) {
                      bT_index = num_buf*buf_size*sizeof(double)+dummy;
                      wreadw(bA_file,(char *) buffer[0], 
                             sizeof(double)*extra_buf, bA_index, &bA_index);
                      for (k=i; k<L; k++) { /* # columns of bA */
                         wreadw(b_file, (char *) buffer[1], 
                                sizeof(double)*extra_buf, bT_index, &jnk);
                         dot_arr(buffer[0], buffer[1], extra_buf, &tmp);
                         G[INDEX(i,k)] += tmp;
                         bT_index += N*sizeof(double);
                         } /* end k loop */
                      } /* end if extra_buf */
                } /* end else */
             i++;
          } /* end while loop */
   /* end form G matrix */

   if (print_lvl > 4) {
     fprintf(outfile,"-----------------------------------------------------\n");
     print_array(G, L,outfile);
     fprintf(outfile,"-----------------------------------------------------\n");
     } 
 
  /* solve the L x L eigenvalue problem G a = lambda a for M roots */
    rsp(L, L, tridim, G, lambda, 1, alpha, 1E-14);

  /* form the d part of the correction vector */
      if (!switch_buf2) {
        num_buf = 1;
        extra_buf = 0;
        }
      for (k=0; k<M; k++) { /* loop over number of roots */
      d_index = (PSI_FPTR) k*N*sizeof(double); /* for iter and multiple roots*/ 
         for (j=0; j<num_buf; j++) { /* for each buffer of d */ 
            /* loop over # buffers for length of bA or b */
            zero_arr(buffer[1], buf_size); 
            for (i=0; i<L; i++) {
               tmp_val = alpha[i][k] * (-lambda[k]); 
               index = (i*N+j*buf_size)*sizeof(double);
                  /* will have to change with iterations */
                  wreadw(bA_file, (char *) buffer[0], sizeof(double)*buf_size, 
                         index, &jnk);
                  for (I=0; I<buf_size; I++) 
                     buffer[1][I] += buffer[0][I]*alpha[i][k];
                  wreadw(b_file,(char *)buffer[0], sizeof(double)*buf_size, 
                         index, &jnk);
                  for (I=0; I<buf_size; I++) 
                     buffer[1][I] += buffer[0][I] * tmp_val;
                } /* end i loop */
                wwritw(d_file,(char *)buffer[1],sizeof(double)*buf_size,
                       d_index,&d_index);
             } /* end j loop */
         if (extra_buf != 0) {
           zero_arr(buffer[1], extra_buf); 
           for (i=0; i<L; i++) {
              tmp_val = alpha[i][k] * (-lambda[k]); 
              index = (i*N+num_buf*buf_size)*sizeof(double);
              wreadw(bA_file, (char *) buffer[0], sizeof(double)*extra_buf, 
                     index, &jnk);
              for (I=0; I<extra_buf; I++)
                 buffer[1][I] += buffer[0][I] * alpha[i][k];
              wreadw(b_file, (char *) buffer[0], sizeof(double)*extra_buf, 
                     index, &jnk);
              for (I=0; I<extra_buf; I++) 
                 buffer[1][I] += buffer[0][I] * tmp_val;
              } /* end i loop */
              wwritw(d_file, (char *)buffer[1], sizeof(double)*extra_buf, 
                     d_index, &d_index);
             } /* end if extra_buf */
         } /* end k loop */

  /* check for convergence */
  for (i=0; i<M; i++) { /* loop over number of roots */
     tval = 0.0; /* necessary for iter and multi-roots */ 
     d_index = (PSI_FPTR) i*N*sizeof(double); 
     for (j=0; j<num_buf; j++) {
        wreadw(d_file, (char *) buffer[0], sizeof(double)*buf_size, 
               d_index, &d_index); 
        dot_arr(buffer[0], buffer[0], buf_size, &tmp);        
        tval += tmp;    
        } /* end j loop */
     if (extra_buf != 0) {
       wreadw(d_file, (char *) buffer[0], sizeof(double)*extra_buf, 
              d_index, &d_index); 
       dot_arr(buffer[0], buffer[0], extra_buf, &tmp);
       tval += tmp;
       } /* end if extra_buf */ 
     tval = sqrt(tval); 
     fprintf(outfile, "Iter %3d  Root %d = %13.9lf", iter-1, i+1, 
             (lambda[i] + nucrep + efzc));
     fprintf(outfile, "    Delta_E %.3E   Delta_C %.3E\n",  
             lambda[i] - lastroot[i], tval);
     fflush(outfile); 
     } /* end i loop */

  if (M > 1) fprintf(outfile,"\n");
  if (fabs(lambda[0]-lastroot[0]) <= conv_e && tval <= conv_rms) {
    converged = 1;
    for (i=0; i<M; i++) {
       evals[i] = lambda[i];
       fprintf(outfile, "\nROOT %d ECI = %17.13lf\n\n\n", i+1, 
               evals[i] + nucrep + efzc); 
       for (k=0; k<num_buf; k++) { /* for b */
          b_index2 = (PSI_FPTR) k*buf_size*sizeof(double); 
          zero_arr(buffer[1], buf_size);  
          for (j=0; j<L; j++) {
             tval = alpha[j][i];
             wreadw(b_file, (char *) buffer[0], sizeof(double)*buf_size, 
                    b_index2, &jnk);
             for (I=0; I<buf_size; I++)
                buffer[1][I] += tval*buffer[0][I]; 
             b_index2 += N*sizeof(double);
             } /* end j loop */
          wwritw(c_file, (char *) buffer[1], sizeof(double)*buf_size, 
                 c_index, &c_index);
          imdet = (N > 20) ? 20 : N; 
          if (k == -1) {
            fprintf(outfile,"The %d most important determinants\n",imdet);
            for (I=0; I<imdet; I++)
               fprintf(outfile,"     %d   %7.6f\n", I, buffer[1][I]); 
            }
          } /* end k loop */
       if (extra_buf != 0) {
         b_index2 = (PSI_FPTR) num_buf*buf_size*sizeof(double); 
         zero_arr(buffer[1], extra_buf); 
         for (j=0; j<L; j++) {
            tval = alpha[j][i]; 
            wreadw(b_file, (char *) buffer[0], sizeof(double)*extra_buf, 
                   b_index2, &jnk); 
            for (I=0; I<extra_buf; I++)
               buffer[1][I] += tval*buffer[0][I];
            b_index2 += N*sizeof(double);
            } /* end j loop */
         wwritw(c_file, (char *) buffer[1], sizeof(double)*extra_buf, 
                c_index, &c_index);
         } /* end if extra_buf */
       } /* end i loop */
     break;  
    } /* end if converged */
  else {
      for (i=0; i<M; i++) lastroot[i] = lambda[i]; 
      } /* end else */        

  /* form the correction vector(s) f from the d vector(s) and normalize */
  d_index = 0; 
  for (k=0; k<M; k++) {
     for (j=0; j<num_buf; j++) {
        wreadw(d_file, (char *) buffer[0], sizeof(double)*buf_size, 
               d_index, &jnk); 
        gendi_a((j*buf_size), buf_size, buffer[1], alplist, betlist);  
        for (I=0; I<buf_size; I++) {
           tval = lambda[k] - buffer[1][I]; 
           if (fabs(tval) < 1.0E-8) 
             buffer[0][I] = 0.0 ; /* the way guga does it */ 
           else 
                 buffer[0][I] /= tval; 
           } /* end I loop */
        wwritw(d_file, (char *) buffer[0], sizeof(double)*buf_size, 
               d_index, &d_index); 
        } /* end j loop */
     if (extra_buf != 0) {
       wreadw(d_file, (char *) buffer[0], sizeof(double)*extra_buf, 
              d_index, &jnk); 
       gendi_a((num_buf*buf_size), extra_buf, buffer[1], alplist, betlist);
       for (I=0; I<extra_buf; I++) {
          tval = lambda[k] - buffer[1][I];
          if (fabs(tval) < 1.0E-8) 
            buffer[0][I] = 0.0; /* the way guga does it */
          else  
              buffer[0][I] /= tval;
          } /* end I loop */
       wwritw(d_file, (char *) buffer[0], sizeof(double)*extra_buf,
              d_index, &d_index);
        } /* end extra_buf */
     } /* end k loop */

  d_index = (PSI_FPTR) 0; 
  for (i=0; i<num_new_vec; i++) {
     v_normalize(buffer[0], d_index, buf_size, extra_buf, num_buf, d_file);  
     d_index += N*sizeof(double); 
     } /* end i loop */

  /* Schmidt orthog and append f's to b */
  num_new_vec = L; 
  for (i=0; i<M; i++) {
     d_index = (PSI_FPTR) i*N*sizeof(double); 
     if(schmidt_addoc(buffer[0], buffer[1], buf_size, extra_buf, num_buf, 
                    d_index, N, L, b_file, d_file)) L++; 
     } /* end i loop */
  num_new_vec = L - num_new_vec;  

  if (L>maxnvect) {
    fprintf(outfile, "(test_sem): L(%2d) > maxnvect(%2d)!",L,maxnvect);
    fprintf(outfile, " Aborting!\n");
    exit(0);
    } 

  /* Again Schmidt orthog b's (minimize numerical error) */
     dotval = v_schmidt(buffer[0], buffer[1], buf_size, extra_buf, 
                       num_buf, N, L, b_file);  

  if (switch_buf3) { /* set up three buffer for b.A */
    buffer4 = buffer[1];
    buffer5 = buffer[2];
    buffer[1] = buffer[0]+buf_size3;
    buffer[2] = buffer[0]+2*buf_size3;
    }

  iter++; 
  free(dotval);
 } /* iterate till convergence (while loop) */ 

 if (converged == 0) {
   for (i=0; i<M; i++) {
      evals[i] = lambda[i];
      for (k=0; k<num_buf; k++) { /* for b */
         b_index2 = (PSI_FPTR) k*buf_size*sizeof(double);
         zero_arr(buffer[1], buf_size);
         for (j=0; j<L; j++) {
            tval = alpha[j][i];
            wreadw(b_file, (char *) buffer[0], sizeof(double)*buf_size,
                   b_index2, &jnk);
            for (I=0; I<buf_size; I++)
               buffer[1][I] += tval*buffer[0][I];
            b_index2 += N*sizeof(double);
            } /* end j loop */
         wwritw(c_file, (char *) buffer[1], sizeof(double)*buf_size,
                c_index, &c_index);
         } /* end k loop */
      if (extra_buf != 0) {
        b_index2 = (PSI_FPTR) num_buf*buf_size*sizeof(double);
        zero_arr(buffer[1], extra_buf);
        for (j=0; j<L; j++) {
           tval = alpha[j][i];
           wreadw(b_file, (char *) buffer[0], sizeof(double)*extra_buf,
                  b_index2, &jnk);
           for (I=0; I<extra_buf; I++)
              buffer[1][I] += tval*buffer[0][I];
           b_index2 += N*sizeof(double);
           } /* end j loop */
        wwritw(c_file, (char *) buffer[1], sizeof(double)*extra_buf,
               c_index, &c_index);
        } /* end if extra_buf */
      fprintf(outfile, "\nROOT %d ECI = %17.13lf\n\n\n", i+1,
              evals[i] + nucrep + efzc);
      } /* end i loop */
   fprintf(outfile,"Maximum number of iterations exceded.\n");
   }

/* free up arrays */
free(buffer); 
free(lastroot);
free(lambda); 
free(alpha); 
free(evals);
free(G);  
/* free(Gprime); */

/* close all files */
rclose(b_file,3); 
rclose(bA_file,3);
rclose(d_file,3); 
rclose(c_file,3); 

}


/*********************************************/
/* gendi_a generates the diagonal matrix     */ 
/* elements of A in batches                  */
/*********************************************/
void gendi_a(unsigned long i, unsigned long buf_size, double *A,
             struct stringwr **alplist, struct stringwr **betlist)  
{
 int k ; 
 int Ialist, Iarel, Iblist, Ibrel;
 static SlaterDeterminant I;

   for (k=i; k<(i+buf_size); k++) {
      det2strings(k, &Ialist, &Iarel, &Iblist, &Ibrel);
      I.set(CalcInfo.num_alp_expl, 
            alplist[Ialist][Iarel].occs, CalcInfo.num_bet_expl,
            betlist[Iblist][Ibrel].occs);
      A[k-i] = matrix_element(&I, &I);
      }
}


/*********************************************/
/* gen_mat generates the matrix elements of  */
/* A in batches using the slater class.      */
/*********************************************/
void gen_mat(unsigned long i, unsigned long j, int buf_size, double *A,
             struct stringwr **alplist, struct stringwr **betlist)  
{
 int k ; 
 int Ialist, Iarel, Iblist, Ibrel;
 static SlaterDeterminant I, J; 
 
   det2strings(j, &Ialist, &Iarel, &Iblist, &Ibrel);
   J.set(CalcInfo.num_alp_expl,
            alplist[Ialist][Iarel].occs, CalcInfo.num_bet_expl,
            betlist[Iblist][Ibrel].occs);
   for (k=i; k<(i+buf_size); k++) {
      det2strings(k, &Ialist, &Iarel, &Iblist, &Ibrel);
      I.set(CalcInfo.num_alp_expl,
            alplist[Ialist][Iarel].occs, CalcInfo.num_bet_expl,
            betlist[Iblist][Ibrel].occs);
      A[k-i] = matrix_element(&I, &J);
      }  
}


/*
** calc_h0_blk  - calculates h0 from the h0blocksize number of determinants
**
** h0           - matrix representation of H in the basis of h0blocksize
**                determinants.
** detref       - array of determinant numbers 
** detval       - array of diagonal determinant values  
** h0blocksize  - size of h0 matrix
** tridim       - size of lower diagonal of matrix h0
*/
void calc_h0_blk(double *h0, int *detref, double *detval, int h0blocksize, 
             int tridim, struct stringwr **alplist, 
             struct stringwr **betlist) 
{
 int i,j;
 int tmp, Ialist, Iarel, Iblist, Ibrel;
 static SlaterDeterminant I, J;

/* for (i=0; i<h0blocksize; i++) {
    tmp = INDEX(i,i);
    h0[tmp] = detval[i];    
    } */

 for (i=0; i<h0blocksize; i++) {
    det2strings(detref[i], &Ialist, &Iarel, &Iblist, &Ibrel);
    I.set(CalcInfo.num_alp_expl,
          alplist[Ialist][Iarel].occs, CalcInfo.num_bet_expl,
          betlist[Iblist][Ibrel].occs);
    for (j=i; j<h0blocksize; j++) {
       tmp = INDEX(i,j);
       det2strings(detref[j], &Ialist, &Iarel, &Iblist, &Ibrel);
       J.set(CalcInfo.num_alp_expl,
             alplist[Ialist][Iarel].occs, CalcInfo.num_bet_expl,
             betlist[Iblist][Ibrel].occs);
       h0[tmp] = matrix_element(&I, &J);
       }
    }
}

/*
**
** sigma: generates the sigma matrix (b . A) for use in the davidson
**        liu algorithm.  
*/
PSI_FPTR sigma(struct stringwr **alplist, struct stringwr **betlist,
                        double **buffer, int L, int switch_buf3,
                        int buf_size, int num_buf, int num_new_vec, 
                        int extra_buf3, int b_file, int bA_file, 
                        int buf_val3, int N, int extra_buf,
                        PSI_FPTR byte, int print_lvl)
{
 int rownum;
 register int i, j, k, I, m;
 PSI_FPTR b_index = 0;
 double tmp = 0;
 int A_col_index;

 if (!switch_buf3) { /* use large buffer method */
   rownum = 0;
   for (m=0; m<num_new_vec; m++) { /* loop of number of roots */
      b_index = (PSI_FPTR) sizeof(double)*N*(L-num_new_vec+m);
      wreadw(b_file, (char *) buffer[0], sizeof(double)*buf_size,
             b_index, &b_index);
      buf_val3 = (int) (N/num_buf);
      for (i=0; i<buf_val3; i++) {
         for (j=0; j<num_buf; j++)
            gen_mat(rownum, j+i*num_buf, buf_size, buffer[j+2],
                    alplist, betlist);
         for (j=0; j<num_buf; j++) {
            dot_arr(buffer[0], buffer[j+2], buf_size, &tmp);
            buffer[1][i*num_buf+j] += tmp;
            }
         }
      if (extra_buf3) {
        for (i=0; i<extra_buf3; i++) {
           gen_mat(rownum, buf_val3*(num_buf)+i, buf_size, buffer[i+2],
                   alplist, betlist);
           dot_arr(buffer[0], buffer[i+2], buf_size, &tmp);
           buffer[1][i+num_buf*buf_val3] += tmp;
           }
        }
      if (print_lvl > 5) 
        for (I=0; I<buf_size; I++) 
           fprintf(outfile,"sigma[%d] = %15.7f\n", I, buffer[1][I]); 

      wwritw(bA_file, (char *) buffer[1], sizeof(double)*buf_size, 
             byte, &byte);
      zero_arr(buffer[1], buf_size);
      }
   } /* end if switch_buf3 */

 else { /* use small buffer method */
   for (m=0; m<num_new_vec; m++) { /* loop of number of roots */
     b_index = (PSI_FPTR) sizeof(double)*N*(L-num_new_vec+m);
     for (k=0; k<num_buf; k++) { /* num_buf section across A */
        rownum = 0; /* row number variable */
        A_col_index = k*buf_size; /* column index for A */
        for (j=0; j<num_buf; j++) { /* for length of b_vectors */
           wreadw(b_file, (char *) buffer[0], sizeof(double)*buf_size,
                  b_index, &b_index);
           for (i=0; i<buf_size; i++) {
              /* run over columns of A buf_size at a time */
                 gen_mat(rownum, (i+A_col_index), buf_size, buffer[2],
                         alplist, betlist);
                 dot_arr(buffer[0], buffer[2], buf_size, &tmp);
                 buffer[1][i] += tmp;
              } /* end i loop */
            rownum += buf_size;
           } /* end j loop */

        if (extra_buf != 0) { /* for extra_buf of bottom rows of A */
          wreadw(b_file, (char *) buffer[0], sizeof(double)*extra_buf,
                 b_index, &b_index);
          for (i=0; i<buf_size; i++) { /* columns of A */
             gen_mat(rownum, (A_col_index+i), extra_buf, buffer[2],
                     alplist, betlist);
             dot_arr(buffer[0], buffer[2], extra_buf, &tmp);
             buffer[1][i] += tmp;
             } /* end i loop */
           } /* end if extra_buf */
       if (print_lvl > 5) { 
         for (I=0; I<buf_size; I++) 
            fprintf(outfile,"sigma[%d] = %15.7f\n", I+k*buf_size, buffer[1][I]);
         }
       wwritw(bA_file, (char *) buffer[1], sizeof(double)*buf_size, 
              byte, &byte);
       zero_arr(buffer[1], buf_size);
       b_index -= N*sizeof(double);
         } /* end k loop */

     if (extra_buf != 0) {
       /* extra columns of A that do not make up a full buf_size */
       rownum = 0; /* row number variable */
       A_col_index = num_buf*buf_size;
       for (j=0; j<num_buf; j++) { /* for length of b_vectors */
          wreadw(b_file, (char *) buffer[0], sizeof(double)*buf_size,
                 b_index, &b_index);
          for (i=0; i<extra_buf; i++) {
             /* run over columns of A of size extra_buf */
             gen_mat(rownum, (A_col_index+i), buf_size, buffer[2],
                     alplist, betlist);
             dot_arr(buffer[0], buffer[2], buf_size, &tmp);
             buffer[1][i] += tmp;
             } /* end i loop */
          rownum += buf_size;
          } /* end j loop */
      wreadw(b_file, (char *) buffer[0], sizeof(double)*extra_buf,
             b_index, &b_index);
      for (i=0; i<extra_buf; i++) { /* loop over low rows of A in extra_buf */
         gen_mat(rownum, (A_col_index+i), extra_buf, buffer[2],
                 alplist, betlist);
         dot_arr(buffer[0], buffer[2], extra_buf, &tmp);
         buffer[1][i] += tmp;
         } /* end i loop */
       } /* end if extra_buf */
 
     if (print_lvl > 5) {
       for (I=0; I<extra_buf; I++) 
          fprintf(outfile,"sigma[%d] = %15.7f\n", I+num_buf*buf_size,
             buffer[1][I]);
       }
     wwritw(bA_file, (char *) buffer[1], sizeof(double)*extra_buf, byte, &byte);
     zero_arr(buffer[1], buf_size);
     } /* end m loop */
   } /* end else */
   /* end form bA matrix */

 return(byte);
}


/*
** h0_guess - uses a unit guess for the first b vector in the
**              davidson-liu algorithm
**
**
*/
void h0_guess(int alp_code, int alp_idx, int bet_code, int bet_idx,
              int switch_buf3, double *buffer, int buf_size,
              int num_buf, int extra_buf, int b_file,
              PSI_FPTR b_writ, int M, int N, int h0blocksize,
              double nucrep, double efzc, struct stringwr **alplist,
              struct stringwr **betlist)
{
 int i;
 register int j, k, I;
 double *detval, *eigval, *h0, *tmpval;
 int *detref, *tmpref;
 double  **eigvec;
 int tridim, tmp2, num_h0blocksize, extra_h0blocksize;
 int max_num, min_num, offst; 
 double max, min, tmp1;

    /* Form blocksize of diagonal elements of H */
   tridim = h0blocksize*(h0blocksize+1)/2;
   detval = init_array(h0blocksize);
   detref = init_int_array(h0blocksize);
   tmpval = init_array(h0blocksize);
   tmpref = init_int_array(h0blocksize);
   eigval = init_array(h0blocksize);
   h0 = init_array(tridim);
   eigvec = block_matrix(h0blocksize, h0blocksize);
   i = strings2det(CalcInfo.ref_alp_list, CalcInfo.ref_alp_rel,
                   CalcInfo.ref_bet_list, CalcInfo.ref_bet_rel);
   fprintf(outfile,"reference det = %d\n", i);
   fflush(outfile);
   if (N < h0blocksize) {
     fprintf(outfile, "blocksize larger than number of determinants.\n");
     exit(0);
     } 

   num_h0blocksize = N/h0blocksize;
   extra_h0blocksize = N%h0blocksize;

   fprintf(outfile,"num_h0blocksize = %d\n", num_h0blocksize);
   fprintf(outfile,"extra_h0blocksize = %d\n", extra_h0blocksize);
   fflush(outfile);

   for (j=0; j<num_h0blocksize; j++) {
      offst = j*h0blocksize;
      gendi_a(offst, h0blocksize, tmpval, alplist, betlist);
       if (j==0) {
         for (I=0; I<h0blocksize; I++) {
            detref[I] = tmpref[I] = I;
            detval[I] = tmpval[I];
            }
         max_element(detval, h0blocksize, &max, &max_num);
         }
      else {
        for (I=0; I<h0blocksize; I++) { 
           tmpref[I] = I+offst;
           }
        min_element(tmpval, h0blocksize, &min, &min_num);
        while (max > min) {
          tmp1 = max;
          tmp2 = detref[max_num];
          detval[max_num] = min;
          detref[max_num] = tmpref[min_num];
          tmpval[min_num] = tmp1;
          tmpref[min_num] = tmp2;
          max_element(detval, h0blocksize, &max, &max_num);
          min_element(tmpval, h0blocksize, &min, &min_num);
          } 
        }    
     }
   
   if (extra_h0blocksize != 0) {
     offst = h0blocksize*num_h0blocksize;
     gendi_a(offst, extra_h0blocksize, tmpval, alplist, betlist);
     for (I=0; I<extra_h0blocksize; I++) { 
        tmpref[I] = I+offst;
        }
     min_element(tmpval, extra_h0blocksize, &min, &min_num);
     if (max > min) {
          offst = h0blocksize*num_h0blocksize;
          tmp1 = max;
          tmp2 = detref[max_num];
          detval[max_num] = min;
          detref[max_num] = tmpref[min_num];
          tmpval[min_num] = tmp1;
          tmpref[min_num] = tmp2;
          }
      }
 
   /* for (I=0; I<h0blocksize; I++) {
      fprintf(outfile,"detval[%d] = %lf\t detref[%d] = %d\n",
              I, detval[I], I, detref[I]); 
      } */
   calc_h0_blk(h0, detref, detval, h0blocksize, tridim, alplist, betlist);
   rsp(h0blocksize, h0blocksize, tridim, h0, eigval, 1, eigvec, 1.0e-14);
   
   fprintf(outfile,"eigval = %20.8f\n", eigval[0]);
   fprintf(outfile,"efzc   = %20.8f\n", efzc);
   fprintf(outfile,"nucrep = %20.8f\n", nucrep);
   fprintf(outfile,"*** H0 Block Eigenvalue = %20.8f\n\n",
           (eigval[0]+efzc+nucrep));
   fflush(outfile);
    /* Form initial guess b vector */
   if (!switch_buf3) {
     for (k=0; k<M; k++) {
        for (j=0; j<h0blocksize; j++) {
           tmp2 = detref[j];
           buffer[tmp2] = eigvec[j][k];
           }
        wwritw(b_file, (char *) buffer, sizeof(double)*buf_size,
               b_writ, &b_writ);
        } 
        v_normalize(buffer, (b_writ-M*N*sizeof(double)), buf_size, 0,
                    1, b_file);
     }
   else if (switch_buf3) {
     zero_arr(buffer, buf_size);
     for (k=0; k<M; k++) {
        for (j=0; j<num_buf; j++) {
           for (I=0; I<h0blocksize; I++) {
              tmp2 = detref[I];
              if ((tmp2 >= (j*buf_size)) && (tmp1 < ((j+1)*buf_size)))  
                buffer[tmp2-j*buf_size] = eigvec[I][k];
              }
           wwritw(b_file, (char *) buffer, sizeof(double)*buf_size,
                  b_writ, &b_writ);
           zero_arr(buffer, buf_size);
           }
        if (extra_buf != 0) {
          for (I=0; I<h0blocksize; I++) {
             tmp1 = detref[I];
             if ((tmp2 >= num_buf*buf_size) && (tmp2 < (extra_buf+num_buf*buf_size)))
               buffer[tmp2-num_buf*buf_size] = eigvec[I][k];
             }
          wwritw(b_file, (char *) buffer, sizeof(double)*extra_buf,
                 b_writ, &b_writ);
          zero_arr(buffer, extra_buf);
          }
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

   free(detval);
   free(tmpval);
   free(tmpref);
   free(h0);
   free(eigvec);
}

}} // namespace psi::detci

