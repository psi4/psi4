/*! \file
**  \ingroup DETCI
**  \brief Routines to compute sigma = H * c
**
** Here collect the stuff to calculate the sigma vector within the
** framework of the CI vector class.
** Rewrote lots of stuff to handle the three out-of-core cases better.
** 
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia
** February 1996
**
*/


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#define EXTERN
#include "globals.h"
#include "civect.h"

namespace psi { namespace detci {

extern void transp_sigma(double **a, int rows, int cols, int phase);
extern void H0block_gather(double **mat, int al, int bl, int cscode, 
   int mscode, int phase);
extern void b2brepl(unsigned char **occs, int *Jcnt, int **Jij, int **Joij,
   int **Jridx, signed char **Jsgn, struct olsen_graph *Graph,
   int Ilist, int Jlist, int len);
extern void b2brepl_test(unsigned char ***occs, int *Jcnt, int **Jij, 
   int **Joij, int **Jridx, signed char **Jsgn, struct olsen_graph *Graph);
extern void s3_block_bz(int Ialist, int Iblist, int Jalist, 
   int Jblist, int nas, int nbs, int cnas, 
   double *tei, double **C, double **S, double **Cprime, double **Sprime);
extern void set_row_ptrs(int rows, int cols, double **matrix);

#ifdef OLD_CDS_ALG
extern void s1_block_fci(struct stringwr **alplist, 
   struct stringwr **betlist, 
   double **C, double **S, double *oei, double *tei, double *F, 
   int nlists, int nas, int nbs, int sbc, int cbc, int cnbs);
extern void s2_block_fci(struct stringwr **alplist, 
   struct stringwr **betlist, 
   double **C, double **S, double *oei, double *tei, double *F,
   int nlists, int nas, int nbs, int Ia_list, int Ja_list, 
   int Ja_list_nas);
extern void s1_block_ras(struct stringwr **alplist, 
   struct stringwr **betlist, 
   double **C, double **S, double *oei, double *tei, double *F, 
   int nlists, int nas, int nbs, int sbc, int cbc, int cnbs);
extern void s1_block_ras_rotf(int *Cnt[2], int **Ij[2], int **Oij[2],
   int **Ridx[2], signed char **Sgn[2], unsigned char **Toccs,
   double **C, double **S,
   double *oei, double *tei, double *F, int nlists, int nas, int nbs,
   int Ib_list, int Jb_list, int Jb_list_nbs);
extern void s2_block_ras(struct stringwr **alplist, 
   struct stringwr **betlist, 
   double **C, double **S, double *oei, double *tei, double *F,
   int nlists, int nas, int nbs, int sac, int cac, int cnas);
extern void s2_block_ras_rotf(int *Cnt[2], int **Ij[2], int **Oij[2],
   int **Ridx[2], signed char **Sgn[2], unsigned char **Toccs,
   double **C, double **S,
   double *oei, double *tei, double *F, int nlists, int nas, int nbs,
   int Ia_list, int Ja_list, int Ja_list_nbs);
extern void s3_block(struct stringwr *alplist, struct stringwr *betlist,
   double **C, double **S, double *tei, int nas, int nbs,
   int Ja_list, int Jb_list);
extern void s3_block_diag(struct stringwr *alplist,struct stringwr *betlist,
   double **C, double **S, double *tei, int nas, int nbs,
   int Ja_list, int Jb_list);
extern void s3_block_diag_rotf(int *Cnt[2], int **Ij[2], 
   int **Ridx[2], signed char **Sgn[2], double **C, double **S,
   double *tei, int nas, int nbs);
extern void s3_block_rotf(int *Cnt[2], int **Ij[2], 
   int **Ridx[2], signed char **Sgn[2], double **C, double **S,
   double *tei, int nas, int nbs);
#else
extern void s1_block_vfci_thread(struct stringwr **alplist, 
   struct stringwr **betlist,
   double **C, double **S, double *oei, double *tei, double *F,
   int nlists, int nas, int nbs, int Ib_list, int Jb_list, 
   int Jb_list_nbs);
extern void s1_block_vfci(struct stringwr **alplist, 
   struct stringwr **betlist,
   double **C, double **S, double *oei, double *tei, double *F,
   int nlists, int nas, int nbs, int Ib_list, int Jb_list, 
   int Jb_list_nbs);
extern void s1_block_vras(struct stringwr **alplist, 
   struct stringwr **betlist, 
   double **C, double **S, double *oei, double *tei, double *F, 
   int nlists, int nas, int nbs, int sbc, int cbc, int cnbs);
extern void s1_block_vras_thread(struct stringwr **alplist, 
   struct stringwr **betlist, 
   double **C, double **S, double *oei, double *tei, double *F, 
   int nlists, int nas, int nbs, int sbc, int cbc, int cnbs);
extern void s1_block_vras_rotf(int *Cnt[2], int **Ij[2], int **Oij[2],
   int **Ridx[2], signed char **Sgn[2], unsigned char **Toccs,
   double **C, double **S,
   double *oei, double *tei, double *F, int nlists, int nas, int nbs,
   int Ib_list, int Jb_list, int Jb_list_nbs);
extern void s2_block_vfci_thread(struct stringwr **alplist, 
   struct stringwr **betlist, 
   double **C, double **S, double *oei, double *tei, double *F,
   int nlists, int nas, int nbs, int Ia_list, int Ja_list, 
   int Ja_list_nas);
extern void s2_block_vfci(struct stringwr **alplist, 
   struct stringwr **betlist, 
   double **C, double **S, double *oei, double *tei, double *F,
   int nlists, int nas, int nbs, int Ia_list, int Ja_list, 
   int Ja_list_nas);
extern void s2_block_vras(struct stringwr **alplist, 
   struct stringwr **betlist, 
   double **C, double **S, double *oei, double *tei, double *F,
   int nlists, int nas, int nbs, int sac, int cac, int cnas);
extern void s2_block_vras_thread(struct stringwr **alplist, 
   struct stringwr **betlist, 
   double **C, double **S, double *oei, double *tei, double *F,
   int nlists, int nas, int nbs, int sac, int cac, int cnas);
extern void s2_block_vras_rotf(int *Cnt[2], int **Ij[2], int **Oij[2],
   int **Ridx[2], signed char **Sgn[2], unsigned char **Toccs,
   double **C, double **S,
   double *oei, double *tei, double *F, int nlists, int nas, int nbs,
   int Ia_list, int Ja_list, int Ja_list_nbs);
extern void s3_block_vdiag(struct stringwr *alplist,
   struct stringwr *betlist,
   double **C, double **S, double *tei, int nas, int nbs, int cnas,
   int Ib_list, int Ja_list, int Jb_list, int Ib_sym, int Jb_sym,
   double **Cprime, double *F, double *V, double *Sgn, int *L, int *R);
extern void s3_block_v(struct stringwr *alplist,struct stringwr *betlist,
   double **C, double **S, double *tei, int nas, int nbs, int cnas,
   int Ib_list, int Ja_list, int Jb_list, int Ib_sym, int Jb_sym,
   double **Cprime, double *F, double *V, double *Sgn, int *L, int *R);
extern void s3_block_vrotf(int *Cnt[2], int **Ij[2], int **Ridx[2],
   signed char **Sn[2], double **C, double **S, 
   double *tei, int nas, int nbs, int cnas,
   int Ib_list, int Ja_list, int Jb_list, int Ib_sym, int Jb_sym,
   double **Cprime, double *F, double *V, double *Sgn, int *L, int *R);
extern void s3_block_vdiag_rotf(int *Cnt[2], int **Ij[2], int **Ridx[2],
   signed char **Sn[2], double **C, double **S, 
   double *tei, int nas, int nbs, int cnas,
   int Ib_list, int Ja_list, int Jb_list, int Ib_sym, int Jb_sym,
   double **Cprime, double *F, double *V, double *Sgn, int *L, int *R);
#endif

extern unsigned char ***Occs;
extern struct olsen_graph *AlphaG;
extern struct olsen_graph *BetaG;

extern int cc_reqd_sblocks[CI_BLK_MAX];

/* FUNCTION PROTOS THIS MODULE */

void sigma_block(struct stringwr **alplist, struct stringwr **betlist,
      double **cmat, double **smat, double *oei, double *tei, int fci, 
      int cblock, int sblock, int nas, int nbs, int sac, int sbc, 
      int cac, int cbc, int cnas, int cnbs, int cnac, int cnbc, 
      int sbirr, int cbirr, int Ms0);



/* GLOBALS THIS MODULE */

double *F;
int **Jij[2];
int **Joij[2];
int **Jridx[2];
signed char **Jsgn[2];
int *Jcnt[2];
unsigned char **Toccs;
double **transp_tmp = NULL;
double **cprime = NULL, **sprime = NULL;

#ifndef OLD_CDS_ALG
   double *V, *Sgn;
   int *L, *R;
#endif


/*
** sigma_init()
**
** This function initializes all the globals associated with calculating
** the sigma vector.
**
*/
void sigma_init(CIvect& C, CIvect &S, struct stringwr **alplist, 
      struct stringwr **betlist)
{
   int i,j;
   int maxcols=0, maxrows=0;
   int nsingles, max_dim=0;
   unsigned long int bufsz=0;

   if (CalcInfo.sigma_initialized) {
      printf("(sigma_init): sigma_initialized already set to 1\n");
      return;
      }

   for (i=0; i<C.num_blocks; i++) {
      if (C.Ib_size[i] > max_dim) max_dim = C.Ib_size[i];
      if (C.Ia_size[i] > max_dim) max_dim = C.Ia_size[i];
      }
   F = init_array(max_dim);

   #ifndef OLD_CDS_ALG
   Sgn = init_array(max_dim);
   V = init_array(max_dim);
   L = init_int_array(max_dim);
   R = init_int_array(max_dim);
   #endif

   if (Parameters.repl_otf) {
      max_dim += AlphaG->num_el_expl;
      nsingles = AlphaG->num_el_expl * AlphaG->num_orb;
      for (i=0; i<2; i++) {
         Jcnt[i] = init_int_array(max_dim);
         Jij[i] = init_int_matrix(max_dim, nsingles);
         Joij[i] = init_int_matrix(max_dim, nsingles);
         Jridx[i] = init_int_matrix(max_dim, nsingles);
         Jsgn[i] = (signed char **) malloc (max_dim * sizeof(signed char *));
         for (j=0; j<max_dim; j++) { 
            Jsgn[i][j] = (signed char *) malloc (nsingles * 
               sizeof(signed char));
            }
         }

      Toccs = (unsigned char **) malloc (sizeof(unsigned char *) * nsingles);

      /* test out the on-the-fly replacement routines */
      /*
      b2brepl_test(Occs,Jcnt[0],Jij[0],Joij[0],Jridx[0],Jsgn[0],AlphaG);
      */
      }

   /* figure out which C blocks contribute to s */
   s1_contrib = init_int_matrix(S.num_blocks, C.num_blocks);
   s2_contrib = init_int_matrix(S.num_blocks, C.num_blocks);
   s3_contrib = init_int_matrix(S.num_blocks, C.num_blocks);
   if (Parameters.repl_otf)
      sigma_get_contrib_rotf(C, S, s1_contrib, s2_contrib, s3_contrib,
         Jcnt, Jij, Joij, Jridx, Jsgn, Toccs);
   else  
      sigma_get_contrib(alplist, betlist, C, S, s1_contrib, s2_contrib,
         s3_contrib);

   if ((C.icore==2 && C.Ms0 && CalcInfo.ref_sym != 0) || (C.icore==0 &&
         C.Ms0)) {
     for (i=0, maxrows=0, maxcols=0; i<C.num_blocks; i++) {
       if (C.Ia_size[i] > maxrows) maxrows = C.Ia_size[i];
       if (C.Ib_size[i] > maxcols) maxcols = C.Ib_size[i];
     }
     if (maxcols > maxrows) maxrows = maxcols;
     transp_tmp = (double **) malloc (maxrows * sizeof(double *));
     if (transp_tmp == NULL) {
       printf("(sigma_init): Trouble with malloc'ing transp_tmp\n");
     }
     bufsz = C.get_max_blk_size();
     transp_tmp[0] = init_array(bufsz);
     if (transp_tmp[0] == NULL) {
       printf("(sigma_init): Trouble with malloc'ing transp_tmp[0]\n");
     }
   }

   /* make room for cprime and sprime if necessary */
   for (i=0, maxrows=0; i<C.num_blocks; i++) {
     if (C.Ia_size[i] > maxrows) maxrows = C.Ia_size[i];
     if (C.Ib_size[i] > maxcols) maxcols = C.Ib_size[i];
   }
   if ((C.icore==2 && C.Ms0 && CalcInfo.ref_sym != 0) || (C.icore==0 &&
         C.Ms0)) {
     if (maxcols > maxrows) maxrows = maxcols;
   }
   bufsz = C.get_max_blk_size();

   #ifndef OLD_CDS_ALG
   cprime = (double **) malloc (maxrows * sizeof(double *));
   if (cprime == NULL) {
      printf("(sigma_init): Trouble with malloc'ing cprime\n");
   }
   if (C.icore==0 && C.Ms0 && transp_tmp != NULL && transp_tmp[0] != NULL) 
     cprime[0] = transp_tmp[0];
   else 
     cprime[0] = init_array(bufsz);

   if (cprime[0] == NULL) {
     printf("(sigma_init): Trouble with malloc'ing cprime[0]\n");
   }

   if (Parameters.bendazzoli) {
     sprime = (double **) malloc (maxrows * sizeof(double *));
     if (sprime == NULL) {
       printf("(sigma_init): Trouble with malloc'ing sprime\n");
     }
     sprime[0] = init_array(bufsz);
     if (sprime[0] == NULL) {
       printf("(sigma_init): Trouble with malloc'ing sprime[0]\n");
     }
   }

   #else
   if (Parameters.bendazzoli) {
     cprime = (double **) malloc (maxrows * sizeof(double *));
     if (cprime == NULL) {
       printf("(sigma_init): Trouble with malloc'ing cprime\n");
     }
     if (C.icore==0 && C.Ms0 && transp_tmp != NULL && transp_tmp[0] != NULL) 
       cprime[0] = transp_tmp[0];
     else 
       cprime[0] = init_array(bufsz);

     if (cprime[0] == NULL) {
       printf("(sigma_init): Trouble with malloc'ing cprime[0]\n");
     }
     sprime = (double **) malloc (maxrows * sizeof(double *));
     if (sprime == NULL) {
       printf("(sigma_init): Trouble with malloc'ing sprime\n");
     }
     sprime[0] = init_array(bufsz);
     if (sprime[0] == NULL) {
       printf("(sigma_init): Trouble with malloc'ing sprime[0]\n");
     }
   }
   #endif
 
   CalcInfo.sigma_initialized = 1;
}


/*
** sigma()
**
** Routine to get the sigma vector using the CI vector class
**
** Changed into a master function which calls the appropriate subfunction
**
*/
void sigma(struct stringwr **alplist, struct stringwr **betlist,
      CIvect& C, CIvect& S, double *oei, double *tei, int fci, int ivec)
{
   if (!CalcInfo.sigma_initialized) sigma_init(C, S, alplist, betlist);

   switch (C.icore) {
      case 0: 
         sigma_a(alplist, betlist, C, S, oei, tei, fci, ivec);
         break;
      case 1:
         sigma_b(alplist, betlist, C, S, oei, tei, fci, ivec);
         break;
      case 2:
         sigma_c(alplist, betlist, C, S, oei, tei, fci, ivec);
         break;
      default:
         fprintf(stderr, "(sigma): Error, invalid icore option\n");
         break;
      } 

}


      
/*
** sigma_a(): This function computes the sigma vector for a given C vector
**    using the CIvector class.   Is somewhat intelligent about constructing
**    the sigma vector blockwise; will attempt to reduce I/O for out-of-core
**    cases.  This version is for icore==0.
**
** Parameters:
**    alplist  = list of alpha strings with replacements
**    betlist  = same for beta strings
**    C        = current ci vector
**    S        = sigma vector to be computed 
**    oei      = one-electron integrals
**    tei      = two-electron integrals
**    fci      = full-ci flag (helps determine which sigma1 routine called)
**    ivec     = sigma vector number (for write call)
**    
** Notes: assumes M_s = 0 for now
*/
void sigma_a(struct stringwr **alplist, struct stringwr **betlist,
      CIvect& C, CIvect& S, double *oei, double *tei, int fci, int ivec)
{

   int buf, cbuf;
   int sblock, cblock, cblock2;  /* id of sigma and C blocks */
   int i,j,k;
   int sac, sbc, nas, nbs;
   int cac, cbc, cnas, cnbs;
   int do_cblock, do_cblock2;
   int cairr, cbirr, sbirr;
   int did_sblock = 0;
   int phase;

   if (!Parameters.Ms0) phase = 1;
   else phase = ((int) Parameters.S % 2) ? -1 : 1;

   /* this does a sigma subblock at a time: icore==0 */
   for (buf=0; buf<S.buf_per_vect; buf++) {
      S.zero(); 
      did_sblock = 0;
      sblock = S.buf2blk[buf];
      sac = S.Ia_code[sblock];
      sbc = S.Ib_code[sblock];
      nas = S.Ia_size[sblock];
      nbs = S.Ib_size[sblock];
      sbirr = sbc / BetaG->subgr_per_irrep;
      if (sprime != NULL) set_row_ptrs(nas, nbs, sprime);

      for (cbuf=0; cbuf<C.buf_per_vect; cbuf++) {
         do_cblock=0; do_cblock2=0;
         cblock=C.buf2blk[cbuf];
         cblock2 = -1;
         cac = C.Ia_code[cblock];
         cbc = C.Ib_code[cblock];
         cbirr = cbc / BetaG->subgr_per_irrep;
         cairr = cac / AlphaG->subgr_per_irrep;
         if (C.Ms0) cblock2 = C.decode[cbc][cac];
         cnas = C.Ia_size[cblock];
         cnbs = C.Ib_size[cblock];
         if (s1_contrib[sblock][cblock] || s2_contrib[sblock][cblock] ||
             s3_contrib[sblock][cblock]) do_cblock = 1;
         if (C.buf_offdiag[cbuf] && (s1_contrib[sblock][cblock2] || 
             s2_contrib[sblock][cblock2] || s3_contrib[sblock][cblock2])) 
            do_cblock2 = 1;
         if (C.check_zero_block(cblock)) do_cblock = 0;
         if (cblock2 >= 0 && C.check_zero_block(cblock2)) do_cblock2 = 0;
         if (!do_cblock && !do_cblock2) continue;

         C.read(C.cur_vect, cbuf);

         if (do_cblock) {
            if (cprime != NULL) set_row_ptrs(cnas, cnbs, cprime);
            sigma_block(alplist, betlist, C.blocks[cblock], S.blocks[sblock], 
               oei, tei, fci, cblock, sblock, nas, nbs, sac, sbc, cac, cbc, 
               cnas, cnbs, C.num_alpcodes, C.num_betcodes, sbirr, cbirr,
               S.Ms0);
            did_sblock = 1;
            }

         /* what's with this bcopy stuff?  what's going on?  -DS 6/11/96 */
         /* I think I should copy to cblock2 not cblock */
         if (do_cblock2) {
            C.transp_block(cblock, transp_tmp);
//          bcopy((char *) transp_tmp[0], (char *) C.blocks[cblock][0], 
//            cnas * cnbs * sizeof(double));
//          bcopy is non-ANSI.  memcpy reverses the arguments.
            memcpy((void *) C.blocks[cblock][0], (void *) transp_tmp[0],
              cnas * cnbs * sizeof(double));
            /* set_row_ptrs(cnbs, cnas, C.blocks[cblock]); */
            if (cprime != NULL) set_row_ptrs(cnbs, cnas, cprime);
            sigma_block(alplist, betlist, C.blocks[cblock2], S.blocks[sblock], 
               oei, tei, fci, cblock2, sblock, nas, nbs, sac, sbc, 
               cbc, cac, cnbs, cnas, C.num_alpcodes, C.num_betcodes, sbirr,
               cairr, S.Ms0);
            did_sblock = 1;
            }

         } /* end loop over c buffers */

      if (did_sblock) {
         S.set_zero_block(sblock, 0);
         if (S.Ms0) S.set_zero_block(S.decode[sbc][sac], 0);
         }

      if (S.Ms0 && (sac==sbc)) 
         transp_sigma(S.blocks[sblock], nas, nbs, phase);

      H0block_gather(S.blocks[sblock], sac, sbc, 1, Parameters.Ms0, phase);

      if (S.Ms0) {
         if ((int) Parameters.S % 2) S.symmetrize(-1.0, sblock);
         else S.symmetrize(1.0, sblock);
         }
      S.write(ivec, buf);

      } /* end loop over sigma buffers */

}



/*
** sigma_b(): This function computes the sigma vector for a given C vector
**    using the CIvector class.   Is somewhat intelligent about constructing
**    the sigma vector blockwise; will attempt to reduce I/O for out-of-core
**    cases.  This version is for icore=1 (whole vector in-core)
**
** Parameters:
**    alplist  = list of alpha strings with replacements
**    betlist  = same for beta strings
**    C        = current ci vector
**    S        = sigma vector to be computed 
**    oei      = one-electron integrals
**    tei      = two-electron integrals
**    fci      = full-ci flag (helps determine which sigma1 routine called)
**    ivec     = sigma vector number (for write call)
**    
** Notes: I think I removed the M_s = 0 assumption from this one
*/
void sigma_b(struct stringwr **alplist, struct stringwr **betlist,
      CIvect& C, CIvect& S, double *oei, double *tei, int fci, int ivec)
{

   int sblock, cblock;  /* id of sigma and C blocks */
   int sac, sbc, nas, nbs;
   int cac, cbc, cnas, cnbs;
   int sbirr, cbirr;
   int did_sblock = 0;
   int phase;

   if (!Parameters.Ms0) phase = 1;
   else phase = ((int) Parameters.S % 2) ? -1 : 1;

   S.zero(); 
   C.read(C.cur_vect, 0);

   /* loop over unique sigma subblocks */ 
   for (sblock=0; sblock<S.num_blocks; sblock++) {
      if (Parameters.cc && !cc_reqd_sblocks[sblock]) continue;
      did_sblock = 0;
      sac = S.Ia_code[sblock];
      sbc = S.Ib_code[sblock];
      nas = S.Ia_size[sblock];
      nbs = S.Ib_size[sblock];
      if (nas==0 || nbs==0) continue;
      if (S.Ms0 && sbc > sac) continue;
      sbirr = sbc / BetaG->subgr_per_irrep;
      if (sprime != NULL) set_row_ptrs(nas, nbs, sprime);

      for (cblock=0; cblock<C.num_blocks; cblock++) {
         if (C.check_zero_block(cblock)) continue;
         cac = C.Ia_code[cblock];
         cbc = C.Ib_code[cblock];
         cnas = C.Ia_size[cblock];
         cnbs = C.Ib_size[cblock];
         cbirr = cbc / BetaG->subgr_per_irrep;
         if (s1_contrib[sblock][cblock] || s2_contrib[sblock][cblock] ||
             s3_contrib[sblock][cblock]) {
            if (cprime != NULL) set_row_ptrs(cnas, cnbs, cprime);
            sigma_block(alplist, betlist, C.blocks[cblock], S.blocks[sblock], 
               oei, tei, fci, cblock, sblock, nas, nbs, sac, sbc, 
               cac, cbc, cnas, cnbs, C.num_alpcodes, C.num_betcodes, sbirr,
               cbirr, S.Ms0);
            did_sblock = 1;
            }
         } /* end loop over c blocks */

      if (did_sblock) S.set_zero_block(sblock, 0);

      if (S.Ms0 && (sac==sbc)) 
         transp_sigma(S.blocks[sblock], nas, nbs, phase);
      H0block_gather(S.blocks[sblock], sac, sbc, 1, Parameters.Ms0, 
         phase);
      } /* end loop over sigma blocks */

      if (S.Ms0) {
         if ((int) Parameters.S % 2) S.symmetrize(-1.0, 0);
         else S.symmetrize(1.0, 0);
         }

   S.write(ivec, 0);

}


/*
** sigma_c(): This function computes the sigma vector for a given C vector
**    using the CIvector class.   Is somewhat intelligent about constructing
**    the sigma vector blockwise; will attempt to reduce I/O for out-of-core
**    cases.  This version is for icore=2 (irrep at a time)
**
** Parameters:
**    alplist  = list of alpha strings with replacements
**    betlist  = same for beta strings
**    C        = current ci vector
**    S        = sigma vector to be computed 
**    oei      = one-electron integrals
**    tei      = two-electron integrals
**    fci      = full-ci flag (helps determine which sigma1 routine called)
**    ivec     = sigma vector number (for write call)
**    
** Notes: tried to remove Ms=0 assumption
*/
void sigma_c(struct stringwr **alplist, struct stringwr **betlist,
      CIvect& C, CIvect& S, double *oei, double *tei, int fci, int ivec)
{

   int buf, cbuf;
   int sblock, cblock, cblock2;  /* id of sigma and C blocks */
   int sairr;                    /* irrep of alpha string for sigma block */
   int cairr;                    /* irrep of alpha string for C block */
   int sbirr, cbirr;
   int i,j,k;
   int sac, sbc, nas, nbs;
   int cac, cbc, cnas, cnbs;
   int did_sblock = 0;
   int phase;

   if (!Parameters.Ms0) phase = 1;
   else phase = ((int) Parameters.S % 2) ? -1 : 1;


   for (buf=0; buf<S.buf_per_vect; buf++) {
      sairr = S.buf2blk[buf];
      sbirr = sairr ^ CalcInfo.ref_sym;
      S.zero();
      for (cbuf=0; cbuf<C.buf_per_vect; cbuf++) {
         C.read(C.cur_vect, cbuf); /* go ahead and assume it will contrib */
         cairr = C.buf2blk[cbuf];
         cbirr = cairr ^ CalcInfo.ref_sym;

         for (sblock=S.first_ablk[sairr];sblock<=S.last_ablk[sairr];sblock++){
            sac = S.Ia_code[sblock];
            sbc = S.Ib_code[sblock];
            nas = S.Ia_size[sblock];
            nbs = S.Ib_size[sblock];
            did_sblock = 0;

            if (S.Ms0 && (sac < sbc)) continue;
            if (sprime != NULL) set_row_ptrs(nas, nbs, sprime);

            for (cblock=C.first_ablk[cairr]; cblock <= C.last_ablk[cairr];
                  cblock++) {

               cac = C.Ia_code[cblock];
               cbc = C.Ib_code[cblock];
               cnas = C.Ia_size[cblock];
               cnbs = C.Ib_size[cblock];

               if ((s1_contrib[sblock][cblock] || s2_contrib[sblock][cblock] ||
                    s3_contrib[sblock][cblock]) && 
                    !C.check_zero_block(cblock)) {
		  if (cprime != NULL) set_row_ptrs(cnas, cnbs, cprime);
                  sigma_block(alplist, betlist, C.blocks[cblock], 
                     S.blocks[sblock], oei, tei, fci, cblock,
                     sblock, nas, nbs, sac, sbc, cac, cbc, cnas, cnbs, 
                     C.num_alpcodes, C.num_betcodes, sbirr, cbirr, S.Ms0);
                  did_sblock = 1;
                  }

               if (C.buf_offdiag[cbuf]) {
                  cblock2 = C.decode[cbc][cac];
                  if ((s1_contrib[sblock][cblock2] || 
                       s2_contrib[sblock][cblock2] ||
                       s3_contrib[sblock][cblock2]) &&
                      !C.check_zero_block(cblock2)) {
                     C.transp_block(cblock, transp_tmp);
		     if (cprime != NULL) set_row_ptrs(cnbs, cnas, cprime);
                     sigma_block(alplist, betlist, transp_tmp,S.blocks[sblock],
                        oei, tei, fci, cblock2, sblock, nas, nbs, sac, sbc, 
                        cbc, cac, cnbs, cnas, C.num_alpcodes, C.num_betcodes,
                        sbirr, cairr, S.Ms0);
                     did_sblock = 1;
                     }
                  }
               } /* end loop over C blocks in this irrep */

            if (did_sblock) S.set_zero_block(sblock, 0); 
            } /* end loop over sblock */

         } /* end loop over cbuf */

      /* transpose the diagonal sigma subblocks in this irrep */
      for (sblock=S.first_ablk[sairr];sblock<=S.last_ablk[sairr];sblock++){
         sac = S.Ia_code[sblock];
         sbc = S.Ib_code[sblock];
         nas = S.Ia_size[sblock];
         nbs = S.Ib_size[sblock];
         if (S.Ms0 && (sac==sbc)) transp_sigma(S.blocks[sblock], nas, nbs, 
            phase);

         /* also gather the contributions from sigma to the H0block */
         if (!S.Ms0 || sac >= sbc) {
            H0block_gather(S.blocks[sblock], sac, sbc, 1, Parameters.Ms0,
               phase);
            }
         }

      if (S.Ms0) {
         if ((int) Parameters.S % 2) S.symmetrize(-1.0, sairr);
         else S.symmetrize(1.0, sairr);
         }
     S.write(ivec, buf);

     } /* end loop over sigma irrep */
}


/* 
** sigma_block()
**
** Calculate the contribution to sigma block sblock from C block cblock
**
*/
void sigma_block(struct stringwr **alplist, struct stringwr **betlist,
      double **cmat, double **smat, double *oei, double *tei, int fci, 
      int cblock, int sblock, int nas, int nbs, int sac, int sbc, 
      int cac, int cbc, int cnas, int cnbs, int cnac, int cnbc, 
      int sbirr, int cbirr, int Ms0)
{

   /* SIGMA2 CONTRIBUTION */
  if (s2_contrib[sblock][cblock]) {

    detci_time.s2_before_time = wall_time_new();

#ifndef OLD_CDS_ALG
      if (fci) {
          if (Parameters.nthreads > 1)
              s2_block_vfci_thread(alplist, betlist, cmat, smat, oei, tei, F, 
                 cnac, nas, nbs, sac, cac, cnas);
          else
              s2_block_vfci(alplist, betlist, cmat, smat, oei, tei, F, cnac,
                            nas, nbs, sac, cac, cnas);
        }
      else {
          if (Parameters.repl_otf) {
              s2_block_vras_rotf(Jcnt, Jij, Joij, Jridx, Jsgn,
                                 Toccs, cmat, smat, oei, tei, F, cnac, 
                                 nas, nbs, sac, cac, cnas);
            }
          else if (Parameters.nthreads > 1) {
            s2_block_vras_thread(alplist, betlist, cmat, smat, 
                            oei, tei, F, cnac, nas, nbs, sac, cac, cnas);  
            }
          else {
              s2_block_vras(alplist, betlist, cmat, smat, 
                            oei, tei, F, cnac, nas, nbs, sac, cac, cnas);
            }
        }
#else
      if (fci) {
          s2_block_fci(alplist, betlist, cmat, smat, oei, tei, F, cnac,
                       nas, nbs, sac, cac, cnas);
        }
      else {
          if (Parameters.repl_otf) {
              s2_block_ras_rotf(Jcnt, Jij, Joij, Jridx, Jsgn,
                                Toccs, cmat, smat, oei, tei, F, cnac, 
                                nas, nbs, sac, cac, cnas);
            }
          else {
              s2_block_ras(alplist, betlist, cmat, smat, 
                           oei, tei, F, cnac, nas, nbs, sac, cac, cnas);
            }
        }
#endif
    detci_time.s2_after_time = wall_time_new();
    detci_time.s2_total_time += detci_time.s2_after_time - detci_time.s2_before_time;

    } /* end sigma2 */

   
   #ifdef DEBUG
   fprintf(outfile, "s2: Contribution to sblock=%d from cblock=%d\n",
      sblock, cblock);
   print_mat(smat, nas, nbs, outfile);
   #endif

   /* SIGMA1 CONTRIBUTION */
   if (!Ms0 || (sac != sbc)) {
    detci_time.s1_before_time = wall_time_new();

      if (s1_contrib[sblock][cblock]) {
         #ifndef OLD_CDS_ALG
          if (fci) { 
              if (Parameters.nthreads > 1)
                  s1_block_vfci_thread(alplist, betlist, cmat, smat, oei, tei, F, cnbc,
                                       nas, nbs, sbc, cbc, cnbs);
              else
                  s1_block_vfci(alplist, betlist, cmat, smat, 
                                oei, tei, F, cnbc, nas, nbs, sbc, cbc, cnbs);
            } 
         else {
            if (Parameters.repl_otf) {
               s1_block_vras_rotf(Jcnt, Jij, Joij, Jridx, Jsgn,
                  Toccs, cmat, smat, oei, tei, F, cnbc, nas, nbs,
                  sbc, cbc, cnbs);
               }
            else if (Parameters.nthreads > 1) {
               s1_block_vras_thread(alplist, betlist, cmat, smat, oei, tei, F, cnbc, 
                  nas, nbs, sbc, cbc, cnbs);
              }
            else {
               s1_block_vras(alplist, betlist, cmat, smat, oei, tei, F, cnbc, 
                  nas, nbs, sbc, cbc, cnbs);
               }
            } 
         #else
         if (fci) { 
            s1_block_fci(alplist, betlist, cmat, smat, 
               oei, tei, F, cnbc, nas, nbs, sbc, cbc, cnbs);
            } 
         else {
            if (Parameters.repl_otf) {
               s1_block_ras_rotf(Jcnt, Jij, Joij, Jridx, Jsgn,
                  Toccs, cmat, smat, oei, tei, F, cnbc, nas, nbs,
                  sbc, cbc, cnbs);
               }
            else {
               s1_block_ras(alplist, betlist, cmat, smat, oei, tei, F, cnbc, 
                  nas, nbs, sbc, cbc, cnbs);
               }
            } 
         #endif
         }
          detci_time.s1_after_time = wall_time_new();
          detci_time.s1_total_time += detci_time.s1_after_time - detci_time.s1_before_time;

      } /* end sigma1 */

   #ifdef DEBUG
   fprintf(outfile, "s1: Contribution to sblock=%d from cblock=%d\n",
      sblock, cblock);
   print_mat(smat, nas, nbs, outfile);
   #endif

   /* SIGMA3 CONTRIBUTION */
   if (s3_contrib[sblock][cblock]) {
      detci_time.s3_before_time = wall_time_new();

      /* zero_mat(smat, nas, nbs); */

      if (!Ms0 || (sac != sbc)) {
         if (Parameters.repl_otf) {
            b2brepl(Occs[sac], Jcnt[0], Jij[0], Joij[0], Jridx[0],
               Jsgn[0], AlphaG, sac, cac, nas);  
            b2brepl(Occs[sbc], Jcnt[1], Jij[1], Joij[1], Jridx[1],
               Jsgn[1], BetaG, sbc, cbc, nbs);  
            #ifndef OLD_CDS_ALG
            s3_block_vrotf(Jcnt, Jij, Jridx, Jsgn, cmat, smat, tei, nas, nbs,
               cnas, sbc, cac, cbc, sbirr, cbirr, cprime, F, V, Sgn, L, R);
            #else
            s3_block_rotf(Jcnt, Jij, Jridx, Jsgn, cmat, smat, tei, nas, nbs);
            #endif
            }      
         else {
         #ifndef OLD_CDS_ALG
            s3_block_v(alplist[sac], betlist[sbc], cmat, smat, tei,
               nas, nbs, cnas, sbc, cac, cbc, sbirr, cbirr, 
               cprime, F, V, Sgn, L, R);
         #else
            s3_block(alplist[sac], betlist[sbc], cmat, smat, 
               tei, nas, nbs, cac, cbc);
         #endif
            }
         }

      else if (Parameters.bendazzoli) {
         s3_block_bz(sac, sbc, cac, cbc, nas, nbs, cnas, tei, cmat, smat, 
            cprime, sprime);
         }

      else {
         if (Parameters.repl_otf) {
            b2brepl(Occs[sac], Jcnt[0], Jij[0], Joij[0], Jridx[0],
               Jsgn[0], AlphaG, sac, cac, nas);  
            b2brepl(Occs[sbc], Jcnt[1], Jij[1], Joij[1], Jridx[1],
               Jsgn[1], BetaG, sbc, cbc, nbs);  
            #ifndef OLD_CDS_ALG
            s3_block_vdiag_rotf(Jcnt, Jij, Jridx, Jsgn, cmat, smat, tei, 
               nas, nbs, cnas, sbc, cac, cbc, sbirr, cbirr, cprime, F, V,
               Sgn, L, R);
            #else
            s3_block_diag_rotf(Jcnt, Jij, Jridx, Jsgn,
               cmat, smat, tei, nas, nbs);
            #endif
            }      
         else {
            #ifndef OLD_CDS_ALG
            s3_block_vdiag(alplist[sac], betlist[sbc], cmat, smat, tei,
               nas, nbs, cnas, sbc, cac, cbc, sbirr, cbirr, 
               cprime, F, V, Sgn, L, R);
            #else
            s3_block_diag(alplist[sac], betlist[sbc], cmat, smat, tei, 
               nas, nbs, cac, cbc);
            #endif
            }
         }

      #ifdef DEBUG
      fprintf(outfile, "s3: Contribution to sblock=%d from cblock=%d\n",
         sblock, cblock);
      print_mat(smat, nas, nbs, outfile);
      #endif
      
      detci_time.s3_after_time = wall_time_new();
      detci_time.s3_total_time += 
         detci_time.s3_after_time - detci_time.s3_before_time;

      } /* end sigma3 */
}
 

}} // namespace psi::detci

