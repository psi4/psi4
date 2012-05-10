/*! \file
**  \ingroup DETCI
**  \brief Code for the CI vector class
** 
** David Sherrill, 15 June 1995
** Center for Comptuational Quantum Chemistry
**
** Update Feb 1996 to ensure everything works symmetry block at a time
** Rewrite blk_xxx members to buf_xxx to avoid confusion with RAS blocks
** 
** Current working assumption for Ms=0: try to fix it so that no disk
** space is required for redundant buffers, but once in memory, assume
** that a buffer has been transposed to give redundant information.  That
** is, assume complete core storage for whole buffer, but don't write 
** redundant buffers to disk.
** Modification: actually, don't store redundant _buffers_, but store
** redundant blocks if needed.
**
*/


#define EXTERN

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libpsio/psio.h>
#include "structs.h"
#include "globals.h"
#include "ci_tol.h"
#include "civect.h"

namespace psi { namespace detci {

extern void calc_hd_block(struct stringwr *alplist, 
   struct stringwr *betlist,
   double **H0, double *oei, double *tei, double efzc,
   int nas, int nbs, int na, int nb, int nbf);
extern void calc_hd_block_ave(struct stringwr *alplist, 
   struct stringwr *betlist, double **H0, double *tf_oei, 
   double *tei, double efzc, int nas, int nbs, int na, int nb, int nbf);
extern void calc_hd_block_z_ave(struct stringwr *alplist,
   struct stringwr *betlist, double **H0, double pert_param,
   double *tei, double efzc, int nas, int nbs, int na, int nb, int nbf);
extern void calc_hd_block_orbenergy(struct stringwr *alplist, 
   struct stringwr *betlist, double **H0, double *oei, 
   double *tei, double efzc, int nas, int nbs, int na, 
   int nb, int nbf);
extern void calc_hd_block_mll(struct stringwr *alplist, 
   struct stringwr *betlist, double **H0, double *oei, 
   double *tei, double efzc, int nas, int nbs, int na, 
   int nb, int nbf);
extern void calc_hd_block_evangelisti(struct stringwr *alplist, 
   struct stringwr *betlist, double **H0, double *tf_oei, 
   double *tei, double efzc, int nas, int nbs, int na, 
   int nb, int nbf);
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
   int nlists, int nas, int nbs, int Ia_list, int Ja_list, 
   int Ja_list_nas);
extern void s2_block_ras_rotf(int *Cnt[2], int **Ij[2], int **Oij[2],
   int **Ridx[2], signed char **Sgn[2], unsigned char **Toccs,
   double **C, double **S,
   double *oei, double *tei, double *F, int nlists, int nas, int nbs,
   int Ia_list, int Ja_list, int Ja_list_nbs);
extern void s3_block(struct stringwr *alplist, struct stringwr *betlist,
   double **C, double **S, double *tei, int nas, int nbs,
   int Ja_list, int Jb_list);
extern void s3_block_diag(struct stringwr *alplist, struct stringwr *betlist,
   double **C, double **S, double *tei, int nas, int nbs,
   int Ja_list, int Jb_list);
extern void s3_block_diag_rotf(int *Cnt[2], int **Ij[2], 
   int **Ridx[2], signed char **Sgn[2], double **C, double **S,
   double *tei, int nas, int nbs);
extern void s3_block_rotf(int *Cnt[2], int **Ij[2], 
   int **Ridx[2], signed char **Sgn[2], double **C, double **S,
   double *tei, int nas, int nbs);
extern void transp_sigma(double **a, int rows, int cols, int phase);
extern void H0block_gather(double **mat, int al, int bl, int cscode, 
   int mscode, int phase);
extern double buf_xy1(double *c, double *hd, double E, int len);
extern void b2brepl(unsigned char **occs, int *Jcnt, int **Jij, int **Joij, 
   int **Jridx, signed char **Jsgn, struct olsen_graph *Graph,
   int Ilist, int Jlist, int len);
extern void b2brepl_test(unsigned char ***occs, int *Jcnt, int **Jij, 
   int **Joij, int **Jridx, signed char **Jsgn, struct olsen_graph *Graph);
extern void xey(double *x, double *y, int size);
extern void xeay(double *x, double a, double *y, int size);
extern void xpeay(double *x, double a, double *y, int size);
extern void xpey(double *x, double *y, int size);
extern void xeax(double *x, double a, int size);
extern void xexmy(double *x, double *y, int size);
extern void calc_d(double *target, double alpha, double *sigma,
   double lambda, double *c, int size);
extern double calc_d2(double *target, double lambda, double *Hd, 
   int size, int precon);
extern double calc_mpn_vec(double *target, double energy, double *Hd, 
   int size, double sign1, double sign2, int precon);
extern void xeaxmy(double *x, double *y, double a, int size);
extern void xeaxpby(double *x, double *y, double a, double b, int size);
extern void xexy(double *x, double *y, int size);
extern void buf_ols_denom(double *a, double *hd, double E, int len);
extern void buf_ols_updt(double *a, double *c, double *norm, double *ovrlap, 
   double *c1norm, int len, FILE *outfile);
extern int H0block_calc(double E);
extern double ssq(struct stringwr *alplist, struct stringwr *betlist,
  double **CL, double **CR, int nas, int nbs, int Ja_list, int Jb_list);

extern unsigned char ***Occs;
extern struct olsen_graph *AlphaG;
extern struct olsen_graph *BetaG;

extern void H0block_coupling_calc(double E, struct stringwr **alplist,
   struct stringwr **betlist);


CIvect::CIvect() // Default constructor
{
   vectlen = 0;
   num_blocks = 0;
   icore = 1;
   Ms0 = 0;
   Ia_code = NULL;
   Ib_code = NULL;
   Ia_size = NULL;
   Ib_size = NULL;
   offset = NULL;
   num_alpcodes = 0;
   num_betcodes = 0;
   nirreps = 0;
   codes_per_irrep = 0;
   buf_per_vect = 0;
   buf_total = 0;
   new_first_buf = 0;
   maxvect = 0;
   nvect = 0;
   nunits = 0;
   cur_vect = -1;
   cur_buf = -1;
   buffer_size = 0;
   units = NULL; 
   file_number = NULL;
   buf_size = NULL;
   buf2blk = NULL;
   buf_offdiag = NULL;
   first_ablk = NULL;
   last_ablk = NULL;
   decode = NULL;
   blocks = NULL;
   zero_blocks = NULL;
   buf_locked = 0;
   buffer = NULL;
   in_file = 0;
   extras = 0;
   units_used = 0;
   cur_unit = 0;
   cur_size = 0;
   first_unit = 0;
}


CIvect::CIvect(BIGINT vl, int nb, int incor, int ms0, int *iac,
         int *ibc, int *ias, int *ibs, BIGINT *offs, int nac, int nbc, 
         int nirr, int cdpirr, int mxv, int nu, 
         int funit, int *fablk, int *lablk, int **dc)
{

   vectlen = 0;
   num_blocks = 0;
   icore = 1;
   Ms0 = 0;
   Ia_code = NULL;
   Ib_code = NULL;
   Ia_size = NULL;
   Ib_size = NULL;
   offset = NULL;
   num_alpcodes = 0;
   num_betcodes = 0;
   nirreps = 0;
   codes_per_irrep = 0;
   buf_per_vect = 0;
   buf_total = 0;
   maxvect = 0;
   nvect = 0;
   nunits = 0;
   cur_vect = -1;
   cur_buf = -1;
   buffer_size = 0;
   units = NULL; 
   file_number = NULL;
   buf_size = NULL;
   buf2blk = NULL;
   buf_offdiag = NULL;
   first_ablk = NULL;
   last_ablk = NULL;
   decode = NULL;
   blocks = NULL;
   zero_blocks = NULL;
   buf_locked = 0;
   buffer = NULL;
   in_file = 0;
   extras = 0;
   units_used = 0;
   cur_unit = 0;
   cur_size = 0;

   set(vl, nb, incor, ms0, iac, ibc, ias, ibs, offs, nac, nbc, 
         nirr, cdpirr, mxv, nu, funit, fablk, lablk, dc);

   buffer = buf_malloc();
   blocks[0][0] = buffer;
   buf_lock(buffer);

}



void CIvect::set(BIGINT vl, int nb, int incor, int ms0, int *iac,
         int *ibc, int *ias, int *ibs, BIGINT *offs, int nac, int nbc, 
         int nirr, int cdpirr, int mxv, int nu, int fu, int *fablk, 
         int *lablk, int **dc)
{
   int i,j,ij,k,l;
   int maxrows = 0, maxcols = 0;
   unsigned long bufsize, maxbufsize;
   unsigned long size, cur_offset;
   static int first=1;
   /* int in_file, extras, units_used, cur_unit; */
   
   vectlen = vl;
   num_blocks = nb;
   icore = incor;
   Ms0 = ms0;
   nirreps = nirr;
   codes_per_irrep = cdpirr;
   maxvect = mxv;
   nvect = 1;
   nunits = nu;
   if (nunits) units = init_int_array(nunits);
   first_unit = fu;
   for (i=0; i<nunits; i++) units[i] = fu + i;

   Ia_code = init_int_array(nb);
   Ib_code = init_int_array(nb);
   Ia_size = init_int_array(nb);
   Ib_size = init_int_array(nb);
   offset = (BIGINT *) malloc (nb * sizeof(BIGINT));
 
   for (i=0; i<nb; i++) {
      Ia_code[i] = iac[i];
      Ib_code[i] = ibc[i];
      Ia_size[i] = ias[i];
      Ib_size[i] = ibs[i];
      offset[i] = offs[i];
      }

   num_alpcodes = nac;
   num_betcodes = nbc;

   first_ablk = init_int_array(nirr);
   last_ablk = init_int_array(nirr);
   for (i=0; i<nirr; i++) {
      first_ablk[i] = fablk[i];
      last_ablk[i] = lablk[i];
      }

   decode = init_int_matrix(nac, nbc);
   for (i=0; i<nac; i++) {
      for (j=0; j<nbc; j++) {
         decode[i][j] = dc[i][j];
         }
      }

   if (icore == 1) { /* whole vector in-core */
      buf_per_vect = 1;
      buf_total = maxvect;
      buf_size = (unsigned long *) malloc(buf_per_vect * sizeof(unsigned long));
      buf2blk = init_int_array(buf_per_vect);
      buf_offdiag = init_int_array(buf_per_vect);
      for (i=0; i<buf_per_vect; i++) buf_size[i] = 0;
      if (maxvect < nunits) nunits = maxvect;
      if (nvect > maxvect) nvect = maxvect;
      size = vectlen;  /* may want to change for Ms=0 later */
      for (i=0; i<buf_per_vect; i++) {
         buf2blk[i] = -1;
         buf_size[i] = size;
         }
      }

   if (icore == 2) { /* whole symmetry block in-core */
      /* figure out how many buffers per vector */
      for (i=0,buf_per_vect=0; i<nirreps; i++) {
         j = first_ablk[i];
         if (j < 0) continue;
         if (!Ms0 || CalcInfo.ref_sym==0) buf_per_vect++;
         else if (Ia_code[j]/codes_per_irrep > Ib_code[j]/codes_per_irrep) 
            buf_per_vect++;
         }

      buf2blk = init_int_array(buf_per_vect);
      buf_offdiag = init_int_array(buf_per_vect);

      for (i=0,j=0; i<nirreps; i++) {
         k = first_ablk[i];
         if (k < 0) continue;
         if (!Ms0 || CalcInfo.ref_sym==0) {
            buf2blk[j] = i;
            j++;
            } 
         else if (Ia_code[k]/codes_per_irrep > Ib_code[k]/codes_per_irrep) {
            buf_offdiag[j] = 1;
            buf2blk[j] = i;
            j++;
            }
         }
           
      buf_total = maxvect * buf_per_vect; 
      buf_size = (unsigned long *) malloc(buf_per_vect * sizeof(unsigned long));
      for (i=0; i<buf_per_vect; i++) {
         buf_size[i] = 0; 
         j = buf2blk[i];
         for (k=first_ablk[j]; k<=last_ablk[j]; k++) {
            buf_size[i] += (unsigned long) Ia_size[k] * 
                           (unsigned long) Ib_size[k];
            }
         }
      } /* end icore==2 */
            
   if (icore == 0) { /* one subblock in-core */

      /* figure out how many blocks are stored */
      buf_per_vect = 0;
      if (Ms0) {
         for (i=0; i<num_blocks; i++) {
            if (Ia_code[i] >= Ib_code[i] && Ia_size[i]>0 && Ib_size[i]>0)
               buf_per_vect++;
            }
         }
      else {
         for (i=0; i<num_blocks; i++) {
            if (Ia_size[i] > 0 && Ib_size[i] > 0) buf_per_vect++;
            }
         }          

      buf_total = buf_per_vect * maxvect; 
      buf2blk = init_int_array(buf_per_vect);
      buf_offdiag = init_int_array(buf_per_vect);
      buf_size = (unsigned long *) malloc(buf_per_vect * sizeof(unsigned long));

      if (Ms0) {
         for (i=0,j=0; i<num_blocks; i++) { 
            if (Ia_code[i] >= Ib_code[i] && Ia_size[i]>0 && Ib_size[i]>0) {
               buf2blk[j] = i;
               buf_size[j] = (unsigned long) Ia_size[i] * 
                             (unsigned long) Ib_size[i];
               if (Ia_code[i] != Ib_code[i]) buf_offdiag[j] = 1;
               j++;
               }
            }
         }
      else {
         for (i=0,j=0; i<num_blocks; i++) { 
            if (Ia_size[i]>0 && Ib_size[i]>0) {
               buf2blk[j] = i;
               buf_size[j] = (unsigned long) Ia_size[i] * 
                             (unsigned long) Ib_size[i];
               j++;
               }
            }
         }

      } /* end icore==0 */

   file_number = init_int_array(buf_total);

   if (nunits) {
      in_file = 0;
      extras = buf_total % nunits;
      units_used = 0;
      cur_unit = units[0];

      for (i=0; i<buf_total; i++) {
  
         if (in_file + 1 <= buf_total / nunits) {
            file_number[i] = cur_unit;
            in_file++;
            }
      
         else if ((in_file == buf_total / nunits) && extras) {
            file_number[i] = cur_unit;
            extras--;
            in_file++;
            }
  
         else {
            units_used++;
            cur_unit = units[units_used];
            file_number[i] = cur_unit;
            in_file = 1;
            } 
         } /* end loop over buffers */
     }
                   
    // do next step separately now to control OPEN_NEW vs OPEN_OLD
    // init_io_files();  

/*
   fprintf(outfile,"num_blocks = %d\n", num_blocks);
   for (i=0; i<buf_total; i++)
      fprintf(outfile,"file_offset[%d] = %lu\n ", i, file_offset[i]);
   for (i=0; i<maxvect; i++)
      fprintf(outfile,"zero_block_offset[%d] = %lu\n ", i,zero_block_offset[i]);
   for (i=0; i<buf_per_vect; i++) 
      fprintf(outfile,"buf_size[%d] = %lu\n ", i, buf_size[i]); 
*/

   // Set up the flags for zero blocks...at first, all blocks are all 0's
   // but we will put '0' indicating that the program should assume
   // nonzero.  It is hard to put this stuff in correctly at this point,
   // and dangerous to throw out any nonzero data, so we will only set
   // the blocks to all zero when we know we want them to be treated as
   // all 0's for certain.
   zero_blocks = init_int_array(num_blocks);

   // Figure out the buffer size and allocate some pointers (but not storage)
   blocks = (double ***) malloc (num_blocks * sizeof(double **));

   if (icore == 1) {   /* everything is in-core */
      for (i=0; i<num_blocks; i++) {
         if (Ia_size[i])
            blocks[i] = (double **) malloc (Ia_size[i] * sizeof(double *));
         else
            blocks[i] = (double **) malloc (sizeof(double *));
         }
      buffer_size = vectlen;
      } /* end icore==1 */
   else if (icore == 2) { /* one symmetry block is held at a time */
     /* figure out which symmetry block is largest */ 
      for (i=0, maxbufsize=0; i<nirreps; i++) {
         for (j=first_ablk[i],bufsize=0; j<=last_ablk[i]; j++) {
            bufsize += (unsigned long) Ia_size[j] * (unsigned long) Ib_size[j];
            }
         if (bufsize > maxbufsize) maxbufsize = bufsize;
         } 
      for (i=0; i<num_blocks; i++) {
         if (Ia_size[i])
            blocks[i] = (double **) malloc (Ia_size[i] * sizeof(double *));
         else
            blocks[i] = (double **) malloc (sizeof(double *));
         }
      buffer_size = maxbufsize;
      } /* end icore==2 */
   else if (icore == 0) { /* only one subblock in core at once */
      for (i=0, maxbufsize=0; i<num_blocks; i++) {
         if (Ia_size[i] > maxrows) maxrows = Ia_size[i];
         if (Ib_size[i] > maxcols) maxcols = Ib_size[i];
         if (Ia_size[i])
            blocks[i] = (double **) malloc (Ia_size[i] * sizeof(double *));
         else
            blocks[i] = (double **) malloc (sizeof(double *));
         bufsize = (unsigned long) Ia_size[i] * (unsigned long) Ib_size[i];
         if (bufsize > maxbufsize) maxbufsize = bufsize;
         }
      // CDS 11/5/97: Revise buffer_size, the size of the biggest buffer 
      // buffer_size = maxrows * maxcols;   Made buffers too large
      buffer_size = maxbufsize;  // Why didn't I do it this way before?
      } /* end icore==0 */
   else {
     printf("CIvect::set(): unrecognized option for icore = %d\n", icore);
     return;
     }
   
   // MLL 5/7/98: Want to know the subblock length of a vector //
   if (first) {
     if (Parameters.print_lvl) {
        fprintf(outfile,"\n CI vector/subblock length = %ld\n", buffer_size);
        fflush(outfile);
        }
     first=0;
     }


}


CIvect::~CIvect()
{
   int i;

   if (num_blocks) {
      if (buf_locked) free(buffer);
      for (i=0; i<num_blocks; i++) {
         free(blocks[i]);
         }

      /* this causes problems if more than one CIvect point to same file 
      for (i=0; i<nunits; i++) {
         rclose(units[i], 3);
         }
      */

      free(blocks);
      free(zero_blocks);
      free(Ia_code);
      free(Ib_code);
      free(Ia_size);
      free(Ib_size);
      free(units);
      free(file_number);
      free(buf_size);
      free(buf2blk);
      free(buf_offdiag);
      free(first_ablk);
      free(last_ablk);
      free_int_matrix(decode);
      free(offset);
      }

}


/* 
** CIvect::buf_malloc
** 
** Function malloc's memory for a buffer of the CI vector.  A buffer
** is the appropriate size to hold either one entire CI vector, 
** a symmetry block of the CI vector, or a subblock of the vector,
** depending on the value of icore.
** Each subblock is a matrix, but the matrices are allocated as
** contiguous blocks of memory.
**
** Parameters:
**    none
**
** Returns:
**    pointer to memory buffer (double *)
*/
double * CIvect::buf_malloc(void)
{
   double *tmp;

   tmp = init_array(buffer_size);

   return(tmp);
}   
   

/*
** CIvect::set_nvect
**
** Sets the value of nvect.  Usually set in CIvect::set().  Only need
** to call if it must be changed in an unusual way (i.e. resetting
** the Davidson subspace).
*/
void CIvect::set_nvect(int i)
{
  nvect = i;
}


/* 
** CIvect::operator *
**
** Function returns the scalar product of two CI vectors.
** Assumes that diagonal blocks are full for Ms=0 cases
*/
double CIvect::operator*(CIvect &b)
{
   double dotprod=0.0, tval;
   int i, len, buf;

   if (Ms0) {
      for (buf=0; buf<buf_per_vect; buf++) {
         read(cur_vect, buf);
         b.read(b.cur_vect, buf);
         dot_arr(buffer, b.buffer, buf_size[buf], &tval);
         if (buf_offdiag[buf]) tval *= 2.0;
         dotprod += tval;
         } 
      }

   else {
      for (buf=0; buf<buf_per_vect; buf++) {
         read(cur_vect, buf);
         b.read(b.cur_vect, buf);
         dot_arr(buffer, b.buffer, buf_size[buf], &tval);
         dotprod += tval;
         } 
      }

   return(dotprod);
}


void CIvect::setarray(const double *a, int len)
{
   double *aptr;
   int i;

   if (len > vectlen) len = vectlen;

   if (icore == 1) {
      aptr = buffer;
      for (i=0; i<len; i++) {
         aptr[i] = a[i];
         }
      }

   else {
      fprintf(stdout, "(CIvect::setarray): Invalid icore option!\n");
      fprintf(stdout, "   use only for icore=1\n");
      } 
}


void CIvect::max_abs_vals(int nval, int *iac, int *ibc, int *iaidx, int *ibidx,
      double *coeff, int neg_only)
{
   int i,buf,irrep;
   double minval=0.0;

 
   if (icore==1) {
      for (i=0; i<num_blocks; i++) {
         minval = blk_max_abs_vals(i, 0, nval, iac, ibc, iaidx, ibidx, coeff,
                    minval, neg_only);
         } 
      } /* end case icore==1 */

   if (icore==2) { /* symmetry block at a time */
      for (buf=0; buf<buf_per_vect; buf++) {
         if (!read(cur_vect, buf)) continue;
         irrep = buf2blk[buf];
         for (i=first_ablk[irrep]; i<=last_ablk[irrep]; i++) {
            minval = blk_max_abs_vals(i, buf_offdiag[buf], nval, iac, ibc, 
               iaidx, ibidx, coeff, minval, neg_only);
            }
         } 
      } /* end case icore==2 */

   if (icore==0) { /* RAS block at a time */
      for (buf=0; buf<buf_per_vect; buf++) {
         if (!read(cur_vect, buf)) continue;
         i = buf2blk[buf];
         minval = blk_max_abs_vals(i, buf_offdiag[buf], nval, iac, ibc, 
            iaidx, ibidx, coeff, minval, neg_only);
         }
      } /* end case icore==0 */

}


double CIvect::blk_max_abs_vals(int i, int offdiag, int nval, int *iac,
       int *ibc, int *iaidx, int *ibidx, double *coeff, double minval,
        int neg_only)
{
   int j,k,m,n;
   double value, abs_value;
   int iacode, ibcode;

   iacode = Ia_code[i];
   ibcode = Ib_code[i];
   for (j=0; j<Ia_size[i]; j++) {
      for (k=0; k<Ib_size[i]; k++) {
         value = blocks[i][j][k];
         if ((value > 0.0) && (neg_only)) continue;
         abs_value = fabs(value);
         if (abs_value >= fabs(minval)) {
            for (m=0; m<nval; m++) {
               if (abs_value > fabs(coeff[m])) {
                  for (n=nval-1; n>m; n--) {
                     coeff[n] = coeff[n-1];
                     iac[n] = iac[n-1];
                     ibc[n] = ibc[n-1];
                     iaidx[n] = iaidx[n-1];
                     ibidx[n] = ibidx[n-1];
                     } 
                  coeff[n] = value;
                  iac[n] = iacode;
                  ibc[n] = ibcode;
                  iaidx[n] = j;
                  ibidx[n] = k;
                  break;
                  } 
               }
            H0block.spin_cp_vals = minval;
            minval = coeff[nval-1];
            }
         if (offdiag) {
            if (Parameters.Ms0 && ((int) Parameters.S % 2) && 
                (!neg_only)) value -= value;
            if (abs_value >= minval) {
               for (m=0; m<nval; m++) {
                  if (abs_value > fabs(coeff[m])) {
                     for (n=nval-1; n>m; n--) {
                        coeff[n] = coeff[n-1];
                        iac[n] = iac[n-1];
                        ibc[n] = ibc[n-1];
                        iaidx[n] = iaidx[n-1];
                        ibidx[n] = ibidx[n-1];
                        } 
                     coeff[n] = value;
                     iac[n] = ibcode;
                     ibc[n] = iacode;
                     iaidx[n] = k;
                     ibidx[n] = j;
                     break;
                     } 
                  }
               H0block.spin_cp_vals = minval;
               minval = coeff[nval-1];
               }
            } 
         }
      }

 /*
   for (i=0; i<H0block.size+H0block.coupling_size; i++)
      fprintf(outfile,"H0block.H00[%d] = %lf\n",i,H0block.H00[i]);
   fprintf(outfile,"H0block.spin_cp_vals = %lf\n",H0block.spin_cp_vals);
   fprintf(outfile,"minval = %lf\n",minval);
   fprintf(outfile,"printed in civect::blk_max_abs_vals\n");
 */
   return(minval);
}


void CIvect::det2strings(BIGINT det, int *alp_code, int *alp_idx,
         int *bet_code, int *bet_idx)
{
   int i;

   /* determine the CI block we're in */
   for (i=0; i<num_blocks-1; i++) {
      if (offset[i+1] > det) break;
      }
   *alp_code = Ia_code[i];
   *bet_code = Ib_code[i];

   *alp_idx = (int) ((det - offset[i]) / (BIGINT) Ib_size[i]);
   *bet_idx = ((det - offset[i]) % (BIGINT) Ib_size[i]);

}

BIGINT CIvect::strings2det(int alp_code, int alp_idx,
      int bet_code, int bet_idx)
{
   int blknum;
   BIGINT addr;

   blknum = decode[alp_code][bet_code];
   addr = offset[blknum];
   addr += alp_idx * Ib_size[blknum] + bet_idx;

   return(addr);
}


void CIvect::diag_mat_els(struct stringwr **alplist, struct stringwr
      **betlist, double *oei, double *tei, double efzc, int na, int nb, 
      int nbf, int method)
{

   int block, buf, iac, ibc, ias, ibs, irrep;
   double minval=0.0;

   if (icore == 1) { /* whole vector in-core */
      for (block=0; block<num_blocks; block++) {
         iac = Ia_code[block];
         ibc = Ib_code[block];
         ias = Ia_size[block];
         ibs = Ib_size[block];
         if (method == HD_KAVE) 
           calc_hd_block_ave(alplist[iac], betlist[ibc], blocks[block], 
              oei, tei, efzc, ias, ibs, na, nb, nbf);
         else if (method == ORB_ENER) 
           calc_hd_block_orbenergy(alplist[iac], betlist[ibc], blocks[block], 
              oei, tei, efzc, ias, ibs, na, nb, nbf);
         else if (method == EVANGELISTI) 
           calc_hd_block_evangelisti(alplist[iac], betlist[ibc], blocks[block], 
              oei, tei, efzc, ias, ibs, na, nb, nbf);
         else if (method == LEININGER) 
           calc_hd_block_mll(alplist[iac], betlist[ibc], blocks[block], 
              oei, tei, efzc, ias, ibs, na, nb, nbf);
         else if (method == HD_EXACT) 
           calc_hd_block(alplist[iac], betlist[ibc], blocks[block],
               oei, tei, efzc, ias, ibs, na, nb, nbf);
         else if (method == Z_HD_KAVE) 
           calc_hd_block_z_ave(alplist[iac], betlist[ibc], blocks[block], 
              Parameters.perturbation_parameter, tei, efzc, ias, ibs, na, 
              nb, nbf);
         else {
           fprintf(outfile," hd_ave option not recognized.\n");
           exit(0);
           }

         if (Parameters.hd_otf && H0block.size) {
            minval = blk_max_abs_vals(block, 0, 
                    (H0block.size+H0block.coupling_size), H0block.alplist, 
                    H0block.betlist, H0block.alpidx, H0block.betidx, 
                    H0block.H00, minval, Parameters.neg_only);      
           }

         }
      if (!Parameters.hd_otf) write(0,0);

      /*
      fprintf(outfile,"Diagonal matrix elements\n");
      print(outfile); 
      */
      } /* end icore==1 */

   else if (icore == 2) { /* whole symmetry block at a time */
      for (buf=0; buf<buf_per_vect; buf++) {
         irrep = buf2blk[buf];
         for (block=first_ablk[irrep]; block<=last_ablk[irrep]; block++) {
            iac = Ia_code[block];
            ibc = Ib_code[block];
            ias = Ia_size[block];
            ibs = Ib_size[block];
            if (method == HD_KAVE)
              calc_hd_block_ave(alplist[iac], betlist[ibc], blocks[block],
                 oei, tei, efzc, ias, ibs, na, nb, nbf);
            else if (method == ORB_ENER)
              calc_hd_block_orbenergy(alplist[iac], betlist[ibc], 
                 blocks[block], oei, tei, efzc, ias, ibs, na, nb, nbf);
            else if (method == EVANGELISTI)
              calc_hd_block_evangelisti(alplist[iac], betlist[ibc], 
                 blocks[block], oei, tei, efzc, ias, ibs, na, nb, nbf);
            else if (method == LEININGER)
              calc_hd_block_mll(alplist[iac], betlist[ibc], 
                 blocks[block], oei, tei, efzc, ias, ibs, na, nb, nbf);
            else if (method == HD_EXACT)
              calc_hd_block(alplist[iac], betlist[ibc], blocks[block], 
                 oei, tei, efzc, ias, ibs, na, nb, nbf);
            else if (method == Z_HD_KAVE)
              calc_hd_block_z_ave(alplist[iac], betlist[ibc], blocks[block],
              Parameters.perturbation_parameter, tei, efzc, ias, ibs, na,
              nb, nbf);
            else {
              fprintf(outfile," hd_ave option not recognized.\n");
              exit(0);
              }

            if (Parameters.hd_otf && H0block.size) {
              minval = blk_max_abs_vals(block, buf_offdiag[buf], 
                    (H0block.size+H0block.coupling_size), 
                    H0block.alplist, H0block.betlist, 
                    H0block.alpidx, H0block.betidx, H0block.H00, 
                    minval, Parameters.neg_only);
              }

            }
         if (!Parameters.hd_otf) write(0,buf);
         }

      } /* end icore==2 */

   else if (icore == 0) { /* one subblock at a time */
      for (buf=0; buf<buf_per_vect; buf++) {
         block = buf2blk[buf];
         iac = Ia_code[block];
         ibc = Ib_code[block];
         ias = Ia_size[block];
         ibs = Ib_size[block];
        if (method == HD_KAVE)
          calc_hd_block_ave(alplist[iac], betlist[ibc], blocks[block],
             oei, tei, efzc, ias, ibs, na, nb, nbf);
        else if (method == ORB_ENER)
          calc_hd_block_orbenergy(alplist[iac], betlist[ibc], blocks[block],
             oei, tei, efzc, ias, ibs, na, nb, nbf);
        else if (method == EVANGELISTI)
          calc_hd_block_evangelisti(alplist[iac], betlist[ibc], blocks[block],
             oei, tei, efzc, ias, ibs, na, nb, nbf);
        else if (method == LEININGER)
          calc_hd_block_mll(alplist[iac], betlist[ibc], blocks[block],
             oei, tei, efzc, ias, ibs, na, nb, nbf);
         else if (method == HD_EXACT)
           calc_hd_block(alplist[iac], betlist[ibc], blocks[block], 
              oei, tei, efzc, ias, ibs, na, nb, nbf);
         else if (method == Z_HD_KAVE)
           calc_hd_block_z_ave(alplist[iac], betlist[ibc], blocks[block],
              Parameters.perturbation_parameter, tei, efzc, ias, ibs, na,
              nb, nbf);
         else {
           fprintf(outfile," hd_ave option not recognized.\n");
           exit(0);
           }
         if (Parameters.hd_otf && H0block.size) {
            minval = blk_max_abs_vals(block, buf_offdiag[buf], 
                    (H0block.size+H0block.coupling_size), 
                    H0block.alplist, H0block.betlist, H0block.alpidx, 
                    H0block.betidx, H0block.H00, minval, Parameters.neg_only);
           }
         if (!Parameters.hd_otf) write(0,buf);
         }
      } /* end icore==0 */

   else {
      printf("(diag_mat_els): Unrecognized icore option!\n");
      }
}


void CIvect::diag_mat_els_otf(struct stringwr **alplist, struct stringwr
      **betlist, double *oei, double *tei, double efzc, int na, int nb,
      int nbf, int buf, int method)
{

   int block, iac, ibc, ias, ibs, irrep;

   if (icore == 1) { /* whole vector in-core */
      for (block=0; block<num_blocks; block++) {
         iac = Ia_code[block];
         ibc = Ib_code[block];
         ias = Ia_size[block];
         ibs = Ib_size[block];
         if (method == HD_KAVE)
           calc_hd_block_ave(alplist[iac], betlist[ibc], blocks[block],
              oei, tei, efzc, ias, ibs, na, nb, nbf);
         else if (method == ORB_ENER)
           calc_hd_block_orbenergy(alplist[iac], betlist[ibc], blocks[block],
              oei, tei, efzc, ias, ibs, na, nb, nbf);
         else if (method == EVANGELISTI)
           calc_hd_block_evangelisti(alplist[iac], betlist[ibc], blocks[block],
              oei, tei, efzc, ias, ibs, na, nb, nbf);
         else if (method == LEININGER)
           calc_hd_block_mll(alplist[iac], betlist[ibc], blocks[block],
              oei, tei, efzc, ias, ibs, na, nb, nbf);
         else if (method == HD_EXACT)
           calc_hd_block(alplist[iac], betlist[ibc], blocks[block],
               oei, tei, efzc, ias, ibs, na, nb, nbf);
         else if (method == Z_HD_KAVE)
           calc_hd_block_z_ave(alplist[iac], betlist[ibc], blocks[block],
              Parameters.perturbation_parameter, tei, efzc, ias, ibs, na,
              nb, nbf);
         else {
           fprintf(outfile," hd_ave option not recognized.\n");
           exit(0);
           }
         }
      } /* end icore==1 */

   else if (icore == 2) { /* whole symmetry block at a time */
         irrep = buf2blk[buf];
         for (block=first_ablk[irrep]; block<=last_ablk[irrep]; block++) {
            iac = Ia_code[block];
            ibc = Ib_code[block];
            ias = Ia_size[block];
            ibs = Ib_size[block];
            if (method == HD_KAVE)
              calc_hd_block_ave(alplist[iac], betlist[ibc], blocks[block],
                 oei, tei, efzc, ias, ibs, na, nb, nbf);
            else if (method == ORB_ENER)
              calc_hd_block_orbenergy(alplist[iac], betlist[ibc],
                 blocks[block], oei, tei, efzc, ias, ibs, na, nb, nbf);
            else if (method == EVANGELISTI)
              calc_hd_block_evangelisti(alplist[iac], betlist[ibc],
                 blocks[block], oei, tei, efzc, ias, ibs, na, nb, nbf);
            else if (method == LEININGER)
              calc_hd_block_mll(alplist[iac], betlist[ibc],
                 blocks[block], oei, tei, efzc, ias, ibs, na, nb, nbf);
            else if (method == HD_EXACT)
              calc_hd_block(alplist[iac], betlist[ibc], blocks[block],
                 oei, tei, efzc, ias, ibs, na, nb, nbf);
            else if (method == Z_HD_KAVE)
              calc_hd_block_z_ave(alplist[iac], betlist[ibc], blocks[block],
              Parameters.perturbation_parameter, tei, efzc, ias, ibs, na,
              nb, nbf);
            else {
              fprintf(outfile," hd_ave option not recognized.\n");
              exit(0);
              }
            }
      } /* end icore==2 */

   else if (icore == 0) { /* one subblock at a time */
         block = buf2blk[buf];
         iac = Ia_code[block];
         ibc = Ib_code[block];
         ias = Ia_size[block];
         ibs = Ib_size[block];
        if (method == HD_KAVE)
          calc_hd_block_ave(alplist[iac], betlist[ibc], blocks[block],
             oei, tei, efzc, ias, ibs, na, nb, nbf);
        else if (method == ORB_ENER)
          calc_hd_block_orbenergy(alplist[iac], betlist[ibc], blocks[block],
             oei, tei, efzc, ias, ibs, na, nb, nbf);
        else if (method == EVANGELISTI)
          calc_hd_block_evangelisti(alplist[iac], betlist[ibc], blocks[block],
             oei, tei, efzc, ias, ibs, na, nb, nbf);
        else if (method == LEININGER)
          calc_hd_block_mll(alplist[iac], betlist[ibc], blocks[block],
             oei, tei, efzc, ias, ibs, na, nb, nbf);
         else if (method == HD_EXACT)
           calc_hd_block(alplist[iac], betlist[ibc], blocks[block],
              oei, tei, efzc, ias, ibs, na, nb, nbf);
         else if (method == Z_HD_KAVE)
           calc_hd_block_z_ave(alplist[iac], betlist[ibc], blocks[block],
              Parameters.perturbation_parameter, tei, efzc, ias, ibs, na,
              nb, nbf);
         else {
           fprintf(outfile," hd_ave option not recognized.\n");
           exit(0);
           }
      } /* end icore==0 */

   else {
      printf("(diag_mat_els): Unrecognized icore option!\n");
      }
}

void CIvect::print(FILE *outfile)
{

   int block, buf, irrep;

   if (cur_vect < 0 || cur_buf < 0) {
      printf("(CIvect::print): Warning...printing unlocked vector\n");
      fprintf(outfile, "[Can't print unlocked vector]\n");
      }

   if (vectlen > 100000) {
      fprintf(outfile, "Not printing long (>100000) vector...\n");
      return;
      }

   if (icore == 0) {
      for (buf=0; buf<buf_per_vect; buf++) {
         read(cur_vect, buf);
         block = buf2blk[buf];
         fprintf(outfile, "\nBlock %2d, codes = (%2d,%2d)\n", block, 
            Ia_code[block], Ib_code[block]);
         print_mat(blocks[block], Ia_size[block], Ib_size[block], outfile);
         }
      }

   else if (icore == 1) {
      for (block=0; block<num_blocks; block++) {
         fprintf(outfile, "\nBlock %2d, codes = (%2d,%2d)\n", block, 
            Ia_code[block], Ib_code[block]);
         print_mat(blocks[block], Ia_size[block], Ib_size[block], outfile);
         }
      }

   else if (icore == 2) {
      for (buf=0; buf<buf_per_vect; buf++) {
         read(cur_vect, buf);
         irrep = buf2blk[buf];
         for (block=first_ablk[irrep]; block<=last_ablk[irrep]; block++) {
            fprintf(outfile, "\nBlock %2d, codes = (%2d,%2d)\n", block, 
               Ia_code[block], Ib_code[block]);
            print_mat(blocks[block], Ia_size[block], Ib_size[block], outfile);
            }
         }
      }
 
   else {
      fprintf(outfile, "(CIvect::print): unrecognized icore option\n");
      }
}


void CIvect::init_vals(int ivect, int nvals, int *alplist, int *alpidx, 
      int *betlist, int *betidx, int *blknums, double *value) 
{
   int i, j, buf, irrep, blk, ai, bi;

   //ok here it seems safe to set zero blocks
   for (i=0; i<num_blocks; i++) zero_blocks[i] = 1;

   /* this used to read >= PARM_GUESS_VEC_H0_BLOCK... but these
      are now gathered from a symnorm so I'll comment this out
      CDS 8/03
   if (Parameters.guess_vector == PARM_GUESS_VEC_H0_BLOCK) {
     for (i=0; i<nvals; i++) 
        H0block.c0b[i] = value[i];
     }
   */

   if (icore == 1) { /* whole vector in-core */
      zero();
      for (i=0; i<nvals; i++) {
         blk = blknums[i];
         ai = alpidx[i];
         bi = betidx[i];
         blocks[blk][ai][bi] = value[i];
         zero_blocks[blk] = 0;
         }
      write(ivect, 0); 
      } /* end icore=1 */

   if (icore == 2) { /* whole symmetry block in core */
      for (buf=0; buf<buf_per_vect; buf++) {
         irrep = buf2blk[buf];
         if (first_ablk[irrep] < 0) continue;
         zero();
         for (blk=first_ablk[irrep]; blk<=last_ablk[irrep]; blk++) {
            for (j=0; j<nvals; j++) {
               if (blknums[j] == blk) {
                  ai = alpidx[j];
                  bi = betidx[j];
                  blocks[blk][ai][bi] = value[j];
                  zero_blocks[blk] = 0;
                  }
               }
            } /* end loop over blocks */
         write(ivect, buf);
         } /* end loop over irreps/bufs */

      } /* end icore=2 */

   if (icore == 0) { /* one subblock at a time */
      for (buf=0; buf<buf_per_vect; buf++) {
         zero();
         for (i=0; i<nvals; i++) {
            blk = blknums[i];
            if (blk == buf2blk[buf]) {
               ai = alpidx[i];
               bi = betidx[i];
               j = ai * Ib_size[blk] + bi;
               buffer[j] = value[i];
               zero_blocks[blk] = 0;
               if (Ms0) 
                  zero_blocks[decode[Ib_code[blk]][Ia_code[blk]]] = 0;
               }
            }
         write(ivect, buf);
         } /* end loop over buf */
      } /* end icore=0 */

}

void CIvect::set_vals(int ivect, int nvals, int *alplist, int *alpidx, 
      int *betlist, int *betidx, int *blknums, double *value) 
{
   int i, j, buf, irrep, blk, ai, bi, vec_modified;
   double tval;

   tval = value[0];

   /* this used to read >= PARM_GUESS_VEC_H0_BLOCK... but these
      are now gathered from a symnorm so I'll comment this out
      CDS 8/03
   if (Parameters.guess_vector == PARM_GUESS_VEC_H0_BLOCK) {
     for (i=0; i<nvals; i++) 
        H0block.c0b[i] = value[i];
     }
   */

   if (icore == 1) { /* whole vector in-core */
      read(ivect, 0);
      vec_modified = 0;
      for (i=0; i<nvals; i++) {
         blk = blknums[i];
         ai = alpidx[i];
         bi = betidx[i];
         blocks[blk][ai][bi] = value[i];
         zero_blocks[blk] = 0;
         vec_modified++;
         }
      if (vec_modified) write(ivect, 0); 
      } /* end icore=1 */

   if (icore == 2) { /* whole symmetry block in core */
      for (buf=0; buf<buf_per_vect; buf++) {
         vec_modified = 0;
         read(ivect, buf); 
         irrep = buf2blk[buf];
         if (first_ablk[irrep] < 0) continue;
         for (blk=first_ablk[irrep]; blk<=last_ablk[irrep]; blk++) {
            for (j=0; j<nvals; j++) {
               if (blknums[j] == blk) {
                  ai = alpidx[j];
                  bi = betidx[j];
                  blocks[blk][ai][bi] = value[j];
                  zero_blocks[blk] = 0;
                  vec_modified++;
                  }
               }
            } /* end loop over blocks */

         if (vec_modified) write(ivect, buf);
         } /* end loop over irreps/bufs */

      } /* end icore=2 */


   if (icore == 0) { /* one subblock at a time */
      for (buf=0; buf<buf_per_vect; buf++) {
         vec_modified = 0;
         read(ivect, buf); 
         for (i=0; i<nvals; i++) {
            blk = blknums[i];
            if (blk == buf2blk[buf]) {
               ai = alpidx[i];
               bi = betidx[i];
               j = ai * Ib_size[blk] + bi;
               buffer[j] = value[i];
               zero_blocks[blk] = 0;
               vec_modified++;
               }
            }
         if (vec_modified) write(ivect, buf);
         } /* end loop over buf */
      } /* end icore=0 */

}


void CIvect::extract_vals(int ivect, int nvals, int *alplist, int *alpidx,
      int *betlist, int *betidx, int *blknums, double *value)
{
   int i, j, buf, irrep, blk, ai, bi, vec_modified;
   double tval;

   tval = value[0];

   if (Parameters.guess_vector == PARM_GUESS_VEC_H0_BLOCK) {
     for (i=0; i<nvals; i++)
        H0block.c0b[i] = value[i];
     }

   if (icore == 1) { /* whole vector in-core */
      read(ivect, 0);
      vec_modified = 0;
      for (i=0; i<nvals; i++) {
         blk = blknums[i];
         ai = alpidx[i];
         bi = betidx[i];
         value[i] = blocks[blk][ai][bi];
         zero_blocks[blk] = 0;
         vec_modified++;
         }
      if (vec_modified) write(ivect, 0);
      } /* end icore=1 */

   if (icore == 2) { /* whole symmetry block in core */
      for (buf=0; buf<buf_per_vect; buf++) {
         vec_modified = 0;
         read(ivect, buf);
         irrep = buf2blk[buf];
         if (first_ablk[irrep] < 0) continue;
         for (blk=first_ablk[irrep]; blk<=last_ablk[irrep]; blk++) {
            for (j=0; j<nvals; j++) {
               if (blknums[j] == blk) {
                  ai = alpidx[j];
                  bi = betidx[j];
                  value[j] = blocks[blk][ai][bi];
                  zero_blocks[blk] = 0;
                  vec_modified++;
                  }
               }
            } /* end loop over blocks */

         if (vec_modified) write(ivect, buf);
         } /* end loop over irreps/bufs */

      } /* end icore=2 */


   if (icore == 0) { /* one subblock at a time */
      for (buf=0; buf<buf_per_vect; buf++) {
         vec_modified = 0;
         read(ivect, buf);
         for (i=0; i<nvals; i++) {
            blk = blknums[i];
            if (blk == buf2blk[buf]) {
               ai = alpidx[i];
               bi = betidx[i];
               j = ai * Ib_size[blk] + bi;
               value[i] = buffer[j];
               zero_blocks[blk] = 0;
               vec_modified++;
               }
            }
         if (vec_modified) write(ivect, buf);
         } /* end loop over buf */
      } /* end icore=0 */

}


/*
** CIvect::symnorm()
**
** This function simultaneously symmetrizes and normalizes a CI vector.
** The symmetrization is according to the rule C(Ia,Ib) = (-1)^S C(Ib,Ia)
** and is valid only for Ms=0.  
**
** Previously assumed that if Ms!=0, this function is not called.
** Now, re-route such calls to new function CIvect::scale().  Due to
** the old assumption, some required calls to symnorm may be missing
** in the iteration routines.
**
** a = norm
** vecode = 0 for C vector and 1 for Sigma vector
** gather_vec = 0 no gather and 1 for gather
**
*/
void CIvect::symnorm(double a, int vecode, int gather_vec)
{
   int i,j;
   int blk,buf,irrep,ac,bc,len,upper;
   double **mat,*arr, phase, tval;

   if (!Ms0) {
      scale(a,vecode,gather_vec);
      return;
      }

   if (!Parameters.Ms0) phase = 1.0;
   else phase = ((int) Parameters.S % 2) ? -1.0 : 1.0;

   if (icore == 1) {

      read(cur_vect, 0);
      for (blk=0; blk<num_blocks; blk++) {
         ac = Ia_code[blk];
         bc = Ib_code[blk];
         mat = blocks[blk];
         if (ac == bc) { /* diagonal block */
            for (i=0; i<Ia_size[blk]; i++) {
               mat[i][i] *= a;
               for (j=0; j<i; j++) {
                  mat[i][j] *= a;
                  mat[j][i] = mat[i][j] * phase;
                  }
               }
            }
         if (ac > bc) { /* off-diagonal block */
            xeax(blocks[blk][0], a, Ia_size[blk] * Ib_size[blk]); 
            upper = decode[bc][ac];
            if (upper >= 0) {
               zero_blocks[upper] = zero_blocks[blk];
               for (i=0; i<Ia_size[blk]; i++) {
                  for (j=0; j<Ib_size[blk]; j++) {
                     blocks[upper][j][i] = mat[i][j] * phase;
                     }
                  }
               }
            }
         } /* end loop over blocks */
        
      if (gather_vec) h0block_gather_vec(vecode); 
      write(cur_vect, 0); 
    
      } /* end icore == 1 */


   else if (icore == 2) { /* irrep at a time */

      for (buf=0; buf<buf_per_vect; buf++) {
         read(cur_vect, buf);
         irrep = buf2blk[buf];
         if (buf_offdiag[buf]) { /* normalize only. other part never stored */
            for (blk=first_ablk[irrep]; blk<=last_ablk[irrep]; blk++) {
                  xeax(blocks[blk][0], a, Ia_size[blk] * Ib_size[blk]); 
               }
            }
         else { /* diagonal irrep, symmetrize and normalize */
            for (blk=first_ablk[irrep]; blk<=last_ablk[irrep]; blk++) {
               ac = Ia_code[blk];
               bc = Ib_code[blk];
               mat = blocks[blk];
               if (ac == bc) { /* diagonal block */
                  for (i=0; i<Ia_size[blk]; i++) {
                     mat[i][i] *= a;
                     for (j=0; j<i; j++) {
                        mat[i][j] *= a;
                        mat[j][i] = mat[i][j] * phase;
                        }
                     }
                  }
               if (ac > bc) { /* off-diagonal block in lower triangle */
                  xeax(blocks[blk][0], a, Ia_size[blk] * Ib_size[blk]); 
                  upper = decode[bc][ac];
                  if (upper >= 0) {
                     zero_blocks[upper] = zero_blocks[blk];
                     for (i=0; i<Ia_size[blk]; i++) {
                        for (j=0; j<Ib_size[blk]; j++) {
                           blocks[upper][j][i] = mat[i][j] * phase;
                           }
                        }
                     }
                  }
               } /* end loop over blocks */
            } /* end diagonal irrep case */
         if (gather_vec) h0block_gather_vec(vecode); 
         write(cur_vect, buf);
         } /* end loop over buffers */
      } /* end icore==2 */

   else if (icore==0) { /* one RAS block at a time */
      for (buf=0; buf<buf_per_vect; buf++) {
         blk = buf2blk[buf];
         read(cur_vect, buf);
         ac = Ia_code[blk];
         bc = Ib_code[blk];
         mat = blocks[blk];
         if (ac == bc) { /* diagonal block */
            for (i=0; i<Ia_size[blk]; i++) {
               mat[i][i] *= a;
               for (j=0; j<i; j++) {
                  mat[i][j] *= a;
                  mat[j][i] = mat[i][j] * phase;
                  }
               }
            }
         else { /* off-diagonal block in lower triangle */
            xeax(blocks[blk][0], a, Ia_size[blk] * Ib_size[blk]); 
            }

         if (gather_vec) h0block_gather_vec(vecode); 
         write(cur_vect, buf);
         } /* end loop over buffers */

      } /* end case icore==0 */

   else {
      printf("(CIvect::symnorm): Unrecognized icore option\n");
      return;
      }

}

/* 
   the following subroutine isn't really effective in keeping out some
   lower-energy states that might creep in due to numerical contamination,
   at least not if it is only used after computing the correction vector...
   I added it trying to keep lower-lying sigmas out of delta state
   computations for C2, but it wasn't effective.  ---CDS October 2004
*/

/*
** CIvect::zero_det()
**
** Zero out a specified determinant
** Implement for icore==1 for now... easy to extend.
*/
double CIvect::zero_det(int iac, int ia, int ibc, int ib)
{
  int blk;
  double tval;

  if (icore != 1) {
    fprintf(outfile, "CIvect::zero_det: Implemented for icore==1 only\n");
    return (0.0);
  }

  blk = decode[iac][ibc];
  tval = blocks[blk][ia][ib];
  printf("zero_det reports coefficient %12.6lf\n", tval);
  tval = tval*tval;
  blocks[blk][ia][ib] = 0.0;

  return(tval); 
}


/*
** CIvect::scale()
**
** This function scales a CI vector by a given factor.
** Does not concern itself with any possible symmetries or redundancies
** in the CI vector.
**
** Parameters:
**   a = the scale factor
**   vecode = 0 if C vector 1 if Sigma vector
**   gather_vec = 1 if gather 0 otherwise
**
** Returns:
**   none
*/
void CIvect::scale(double a, int vecode, int gather_vec)
{
   int buf;

   for (buf=0; buf<buf_per_vect; buf++) {
      read(cur_vect, buf);
      xeax(buffer, a, buf_size[buf]);
      if (gather_vec) h0block_gather_vec(vecode);
      write(cur_vect, buf);
      }
}



/*
** CIvect::buf_lock()
**
** This function "locks in" a memory buffer for the use of a CIvector.
** The appropriate flag is set and pointers are made to point to the
** right regions of the buffer.
**
** Parameters: 
**    a  = array to use for the buffer
** Returns: none
*/
void CIvect::buf_lock(double *a)
{
   int i,j,k;

   if (buf_locked) {
      printf("Warning (CIvect::buf_lock): CIvector is already locked!\n");
      }

   if (icore == 1) { /* whole vector in-core */
      blocks[0][0] = a;
      for (j=1; j<Ia_size[0]; j++) {
         blocks[0][j] = blocks[0][0] + Ib_size[0] * j;
         }
      for (i=1; i<num_blocks; i++) {
         blocks[i][0] = blocks[i-1][0] + Ia_size[i-1] * Ib_size[i-1];
         for (j=1; j<Ia_size[i]; j++) {
            blocks[i][j] = blocks[i][0] + Ib_size[i] * j;
            }
         }
      } /* end icore==1 option */

   if (icore == 2) { /* one symmetry block is held at a time */
      blocks[0][0] = a;
      for (i=0; i<nirreps; i++) {
         for (j=first_ablk[i]; j<=last_ablk[i]; j++) {
            if (j==first_ablk[i]) 
               blocks[j][0] = a;
            else
               blocks[j][0] = blocks[j-1][0] + Ia_size[j-1] * Ib_size[j-1];
            for (k=1; k<Ia_size[j]; k++) 
               blocks[j][k] = blocks[j][0] + Ib_size[j] * k;
            }        
         }
      } /* end icore==2 option */

   if (icore == 0) { /* one subblock at a time */
      for (i=0; i<num_blocks; i++) {
         blocks[i][0] = a;
         for (j=1; j<Ia_size[i]; j++) {
            blocks[i][j] = blocks[i][0] + Ib_size[i] * j;
            }
         }
      } /* end icore==0 option */

   buffer = a;
   buf_locked = 1;
   /* zero(); * commented out 3/13/96: need to eliminate */
} 


/*
** CIvect::buf_unlock()
**
** This function "unlocks" a memory buffer from the use of a CIvector.
**
** Parameters: none
** Returns: none
*/
void CIvect::buf_unlock(void)
{
   buf_locked = 0;
   blocks[0][0] = NULL;
   buffer = NULL;
   cur_vect = -1;
   cur_buf = -1;
}



/*
** CIvect::symmetrize(): This function symmetrizes the CI vector
**    to maintain the appropriate spin symmetry.  Symmetrizes only
**    the in-core portion of the vector.
**    Assume that this function is called only if Ms=0
**
** Parameters:
**    phase = the exponent S in C(Ia,Ib) = (-1)^S * C(Ib,Ia)
**
*/
void CIvect::symmetrize(double phase, int iblock)
{
   int i,j;
   int blk,irrep,ac,bc,len,upper;
   double **mat,*arr;

   if (icore == 1) {

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
               zero_blocks[upper] = zero_blocks[blk];
               for (i=0; i<Ia_size[blk]; i++) {
                  for (j=0; j<Ib_size[blk]; j++) {
                     blocks[upper][j][i] = mat[i][j] * phase;
                     }
                  }
               }
            }
         } /* end loop over blocks */

      } /* end icore == 1 */


   else if (icore == 2) { /* irrep at a time */

      irrep = iblock;

      /* do only for diagonal irrep blocks */
      if (CalcInfo.ref_sym != 0) return;

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
               zero_blocks[upper] = zero_blocks[blk];
               for (i=0; i<Ia_size[blk]; i++) {
                  for (j=0; j<Ib_size[blk]; j++) {
                     blocks[upper][j][i] = mat[i][j] * phase;
                     }
                  }
               }
            }
         } /* end loop over blocks */
      } /* end icore==2 */

   else if (icore==0) { /* one RAS block at a time */
      ac = Ia_code[iblock];
      bc = Ib_code[iblock];
      mat = blocks[iblock];
      if (ac == bc) { /* diagonal block */
         for (i=0; i<Ia_size[iblock]; i++) {
            for (j=0; j<i; j++) {
               mat[j][i] = mat[i][j] * phase;
               }
            }
         }

      } /* end case icore==0 */


   else {
      printf("(CIvect::symmetrize): Unrecognized icore option\n");
      return;
      }

}



/*
** CIvect::blockptr(): Return a pointer to a block so that it can be 
**    accessed outside (yeah, this destroys data protection, but
**    it's better to do this here).
**
*/
double ** CIvect::blockptr(int blknum)
{
   return blocks[blknum];
}


/* 
** CIvect::init_io_files()
**
** Parameters: none
** 
** Returns: none
*/
void CIvect::init_io_files(bool open_old)
{
   int i;

   for (i=0; i<nunits; i++) {
     if (!psio_open_check((ULI) units[i])) {
       if (open_old) psio_open((ULI) units[i], PSIO_OPEN_OLD); 
       else psio_open((ULI) units[i], PSIO_OPEN_NEW);
     }
   }
}



/*
** CIvect::close_io_files()
**
** Parameters: 
**    keep = 1 to keep files, else 0 to delete
**
** Returns: none 
*/
void CIvect::close_io_files(int keep)
{
   int i;

   for (i=0; i<nunits; i++) {
     // rclose(units[i], keep ? 3 : 4); // old way 
     psio_close(units[i], keep); // new way   
   }
}



/*
** CIvect::read(): Read in a section of a CI vector from external storage.
**
** Parameters:
**    ivect  = vector number
**    ibuf   = buffer number (ibuf can specify an irrep or a subblock
**                within an irrep, depending on the value of icore.  If
**                icore = 1, then ibuf is ignored.)
**
** Returns: 1 for success, 0 for failure
*/
int CIvect::read(int ivect, int ibuf)
{
   int unit, buf, k, i;
   unsigned long int size;
   int blk;
   char key[20];

   detci_time.read_before_time = wall_time_new();

   if (nunits < 1) {
      cur_vect = ivect;
      cur_buf = ibuf;
      return(1);
      }

   if (ivect < 0 || ibuf < 0) {
      printf("(CIvect::read): Called with negative argument\n");
      return(0);
      }

   if (icore == 1) ibuf = 0;
   buf = ivect * buf_per_vect + ibuf;

   size = buf_size[ibuf] * (unsigned long int) sizeof(double);   

   /* translate buffer number in case we renumbered after collapse * */
   buf += new_first_buf;
   if (buf >= buf_total) buf -= buf_total;
   sprintf(key, "buffer %d", buf);
   unit = file_number[buf];

   psio_read_entry((ULI) unit, key, (char *) buffer, size);  

   cur_vect = ivect;
   cur_buf = ibuf;

   detci_time.read_after_time = wall_time_new();
   detci_time.read_total_time += detci_time.read_after_time - 
     detci_time.read_before_time;

   return(1);
}  


/*
** CIvect::write(): Write a section of a CI vector to external storage.
**
** Parameters:
**    ivect  = vector number
**    ibuf   = buffer number (ibuf can specify an irrep or a subblock
**                within an irrep, depending on the value of icore.  If
**                icore = 1, then ibuf is ignored.)
**
** Returns: 1 for success, 0 for failure
*/
int CIvect::write(int ivect, int ibuf)
{
   int unit, buf, i;
   unsigned long int size;
   int blk;
   char key[20];

   detci_time.write_before_time = wall_time_new();

   if (nunits < 1) return(1);

   if (ivect >= maxvect) {
      fprintf(outfile, "(CIvect::write): ivect >= maxvect\n");
      return(0);
      }

   if (ivect > nvect) { 
      fprintf(outfile, "(CIvect::write): ivect > nvect\n");
      return(0);
      }
   
   if (icore == 1) ibuf = 0;
   buf = ivect * buf_per_vect + ibuf;
   size = buf_size[ibuf] * (unsigned long int) sizeof(double);   

   /* translate buffer number in case we renumbered after collapse * */
   buf += new_first_buf;
   if (buf >= buf_total) buf -= buf_total;
   sprintf(key, "buffer %d", buf);
   unit = file_number[buf];
  
   psio_write_entry((ULI) unit, key, (char *) buffer, size);

   if (ivect >= nvect) nvect = ivect + 1;
   cur_vect = ivect;
   cur_buf = ibuf;
   
   detci_time.write_after_time = wall_time_new();
   detci_time.write_total_time += detci_time.write_after_time - 
     detci_time.write_before_time;

   return(1);
}  


/*
** CIvect::schmidt_add()
**
** This function Gram-Schmidt orthogonalizes a new vector d and adds it to
** the list of vectors in the CIvector c, which must contain room
** for the new vector (i.e. after the new vector is added, nvect <= maxvect).
** Don't add orthogonalized d' if norm(d') < SA_NORM_TOL.
**
** Parameters:
**    L   = number of vectors in CIvect to consider
**
** Returns: 1 if a vector is added, 0 otherwise
**
** Notes: Assumes vectors c,d are same size.  Should account for Ms0 now.
*/
int CIvect::schmidt_add(CIvect &c, int L)
{
   double tval, norm, *dotval;
   int buf, cvect;

   norm = 0.0;

   dotval = init_array(L);

   for (buf=0; buf<buf_per_vect; buf++) {
      read(cur_vect, buf);
      for (cvect=0; cvect<L; cvect++) {
         c.read(cvect, buf);
         dot_arr(buffer, c.buffer, buf_size[buf], &tval);
         if (buf_offdiag[buf]) tval *= 2.0;
         dotval[cvect] += tval;
         }
      }

   for (buf=0; buf<buf_per_vect; buf++) {
      read(cur_vect, buf);
      for (cvect=0; cvect<L; cvect++) {
         c.read(cvect, buf);
       /* 
         fprintf(outfile,"dotval[%d] = %2.15lf\n",cvect,dotval[cvect]); 
       */
         xpeay(buffer, -dotval[cvect], c.buffer, buf_size[buf]);
         }
      dot_arr(buffer, buffer, buf_size[buf], &tval);
      if (buf_offdiag[buf]) tval *= 2.0;
      norm += tval;
      write(cur_vect, buf);
      }

   free(dotval);

   norm = sqrt(norm);
   if (norm < SA_NORM_TOL) return(0);
   norm = 1.0 / norm;

   if (c.nvect > c.maxvect) {
      fprintf(stderr, "(CIvect::schmidt_add): no more room to add vectors!\n");
      fprintf(stderr, "   c.nvect = %d, c.maxvect = %d\n", c.nvect, c.maxvect);
      return(0);
      }
   else { /* add to c */
      c.cur_vect = c.nvect;
      c.nvect++;
      for (buf=0; buf<buf_per_vect; buf++) {
         read(cur_vect, buf);
         xeay(c.buffer, norm, buffer, buf_size[buf]);
         c.write(c.cur_vect, buf);
         }
      return(1);
      }

}



/*
** CIvect::schmidt_add2()
**
** This function Gram-Schmidt orthogonalizes a new vector d and adds it to
** the list of vectors in the CIvector c, which must contain room
** for the new vector (i.e. after the new vector is added, nvect <= maxvect).
** Don't add orthogonalized d' if norm(d') < SA_NORM_TOL.  This version
** differs from CIvect::schmidt_add() to allow the user somewhat finer
** control over the numbering of vectors, etc, and allows the return of
** dot products and the normalization factor.
**
** Parameters:
**    L          = number of vectors in CIvect to consider
**    first_vec  = first vec num in C to orthogonalize new vector against
**    last_vec   = last vec num in C to orthogonalize new vector against
**    source_vec = vector number in D file to read 
**    target_vec = vector number to write new vector to
**    dotval     = array of dot products of new vector with old vectors
**    nrm        = normalization constant of new vector after orthogonalization 
**    ovlpmax    = maximum overlap of current SO vectors to previous vectors
**
** Returns: 1 if a vector is added, 0 otherwise
**
** Notes: Assumes vectors c,d are same size.  Should account for Ms0 now.
*/
int CIvect::schmidt_add2(CIvect &c, int first_vec, int last_vec, 
   int source_vec, int target_vec, double *dotval, double *nrm, double *ovlpmax)
{
   double tval, norm, *dotchk, dotck, tmp_norm;
   int buf, cvect, i;

   norm = 0.0;
   tmp_norm = 0.0;
   dotchk = init_array(100);
   *ovlpmax = 0.0;

   for (buf=0; buf<buf_per_vect; buf++) {
      read(source_vec, buf);
      for (cvect=first_vec; cvect<=last_vec; cvect++) {
         c.read(cvect, buf);
         dot_arr(buffer, c.buffer, buf_size[buf], &tval);
         if (buf_offdiag[buf]) tval *= 2.0;
         dotval[cvect] += tval;
         }
      }

   for (i=first_vec; i<=last_vec; i++) {
      tval = fabs(dotval[i]);
      if (tval>*ovlpmax) *ovlpmax = tval;
      }

   /* Schmidt orthogonalize and double check orthogonalization */
   for (buf=0; buf<buf_per_vect; buf++) {
      read(cur_vect, buf);
      for (cvect=first_vec; cvect<=last_vec; cvect++) {
         c.read(cvect, buf);
         xpeay(buffer, -dotval[cvect], c.buffer, buf_size[buf]);
         }
      dot_arr(buffer, buffer, buf_size[buf], &tval); 
      if (buf_offdiag[buf]) tval *= 2.0;
      norm += tval;
      write(cur_vect, buf);
      }

   /* fprintf(outfile,"Norm of %d vec = %20.15f\n",target_vec,norm); */
   norm = sqrt(norm); 
   /* fprintf(outfile,"sqrt Norm of %d vec = %20.15f\n",target_vec,norm); */
   if (Parameters.mpn_schmidt) 
     if (norm < MPn_NORM_TOL) return(0);
   else if (norm < SA_NORM_TOL) return(0);
 /*  
   if (norm < SA_NORM_TOL && !Parameters.mpn) return(0);
 */
   norm = 1.0 / norm;
   /* fprintf(outfile,"1.0/sqrt(norm) of %d vec = %20.15f\n",target_vec,norm); */
   *nrm = norm;

   if (c.nvect > c.maxvect) {
      fprintf(stderr, "(CIvect::schmidt_add2): no more room to add vectors!\n");
      fprintf(stderr, "   c.nvect = %d, c.maxvect = %d\n", c.nvect, c.maxvect);
      return(0);
      }
   else { /* add to c */
      c.cur_vect = target_vec;
      if (c.cur_vect > c.nvect) c.nvect++;
      zero_arr(dotchk,100);

      for (buf=0; buf<buf_per_vect; buf++) {
         read(cur_vect, buf); 
         xeay(c.buffer, norm, buffer, buf_size[buf]);
         c.write(c.cur_vect, buf);
         }
     /*
      fprintf(outfile, "c.cur_vect = %d\n",c.cur_vect);
      fprintf(outfile, "dot product of normalized vector in SA2 = %20.10f\n",
              tmp_norm);
      c.print(outfile);
     */

      if (Parameters.mpn) {
        zero_arr(dotchk,100);
        for (buf=0; buf<buf_per_vect; buf++) {
           read(source_vec, buf);
           for (cvect=first_vec; cvect<=last_vec; cvect++) {
              c.read(cvect, buf);
              dot_arr(buffer, c.buffer, buf_size[buf], &tval);
              if (buf_offdiag[buf]) tval *= 2.0;
              dotchk[cvect] += tval;
              }
           }
        for (i=first_vec; i<=last_vec; i++)
           if (dotchk[i] > *ovlpmax) *ovlpmax = dotchk[i];
        }
      return(1);
      }

}



/*
** CIvect::zero()
**
** Zero out the current memory buffer for a CI vector.
**
** Parameters: none
** Returns: none
**/
void CIvect::zero(void)
{
   zero_arr(buffer, (int) buffer_size);
}



/*
** CIvect::sigma_renorm()
**
** Function calculates the numerator part of the Davidson correction vector d
**
** Parameters:
**    nr       = number of roots (=number of d vectors to calculate)
**    L        = number of previous vectors in CI subspace
**    alpha    = subspace CI eigenvector matrix
**    lambda   = array of subspace eigenvalues
**    norm_arr = norm array (hold norm of for each d vector) 
**    C        = CIvect for subspace vectors
**    printflag= 1 to print d vector(s), else 0
**    outfile  = where to put any output
**
** Returns: none
*/
void CIvect::sigma_renorm(int nr, int L, double renorm_C, CIvect &S, 
             double *buf1, int printflag, FILE *outfile)
{
   int buf, ivect, root;
   double tval;

      for (buf=0; buf<buf_per_vect; buf++) {
         for (ivect=0; ivect<L; ivect++) {
            S.buf_lock(buf1);
            S.read(ivect, buf);
            xeay(S.buffer, renorm_C, S.buffer, buf_size[buf]);
            S.buf_unlock();
            } /* end loop over ivect */

         write(nr, buf);
         if (printflag) {
            fprintf(outfile, "\nSigma renormalized matrix\n");
            print_buf(outfile);
            }
         } /* loop over buffers */
}


/*
** CIvect::dcalc()
**
** Function calculates the numerator part of the Davidson correction vector d
**
** Parameters:
**    nr       = number of roots (=number of d vectors to calculate)
**    L        = number of previous vectors in CI subspace
**    alpha    = subspace CI eigenvector matrix
**    lambda   = array of subspace eigenvalues
**    norm_arr = norm array (hold norm of for each d vector) 
**    C        = CIvect for subspace vectors
**    S        = CIvect for sigma vectors
**    root_converged = 1 if root has converged, else 0
**    printflag= 1 to print d vector(s), else 0
**    outfile  = where to put any output
**    E_est     = Intermediate is OLSEN update
**
** Returns: none
*/
void CIvect::dcalc(int nr, int L, double **alpha, double *lambda, 
      double *norm_arr, CIvect &C, CIvect &S, double *buf1, double *buf2, 
      int *root_converged, int printflag, FILE *outfile, double *E_est)
{
   int buf, ivect, root, tmproot, converged=0, i;
   double tval;


   buf_lock(buf2);

   /* Calculate the d vector for each converged root but do
   ** not form the complete correction (f) vector
   */

   /* 
   for (root=0; root<nr; root++) 
      if (root_converged[root]) {
        converged = root;
        break;
        }

   if (converged) nr++; 
   */

   for (root=0; root<nr; root++) {
      norm_arr[root] = 0.0;
      /* if (converged && root==nr-1 && L<nr) break; */
      for (buf=0; buf<buf_per_vect; buf++) {
         zero();
         if (Parameters.update==UPDATE_OLSEN) {
           read(root,buf);
           xeax(buffer, -E_est[root], buf_size[buf]); 
           /* buffer is know E_est*C^k */ 
           }
         for (ivect=0; ivect<L; ivect++) {
            if (Parameters.update == UPDATE_DAVIDSON) { /* DAVIDSON update formula */
              C.buf_lock(buf1);
              C.read(ivect, buf);
              tval = -alpha[ivect][root] * lambda[root];
              xpeay(buffer, tval, C.buffer, buf_size[buf]);
              C.buf_unlock();
              }
            S.buf_lock(buf1);
            S.read(ivect, buf);
            xpeay(buffer, alpha[ivect][root], S.buffer, buf_size[buf]);
            S.buf_unlock();
            } /* end loop over ivect */
         dot_arr(buffer, buffer, buf_size[buf], &tval);
         if (buf_offdiag[buf]) tval *= 2.0;
         norm_arr[root] += tval;
         write(root, buf);
         /* 
         if (root==nr-1 && converged) write(converged, buf);
         else write(root, buf);
         */

         if (printflag) {
            fprintf(outfile, "\nfirst D matrix\n");
            print_buf(outfile);
            }
         } /* loop over buffers */

      norm_arr[root] = sqrt(norm_arr[root]);

      } /* end loop over roots */

   buf_unlock(); 
}



/*
** 
** CIvect::dcalc2()
**
** Function calculates the denominator part of the Davidson correction 
** vector 'd'.
**
** Parameters:
**    rootnum   = number of current root
**    lambda    = current iteration's energy eigenvalue for the current root
**    Hd        = CIvector for the Hamiltonian diagonal
**
** Returns: sum of squares of coefficients in d.
*/
double CIvect::dcalc2(int rootnum, double lambda, CIvect &Hd, 
        int precon, struct stringwr **alplist, struct stringwr **betlist)
{
   int buf, errcod, i;
   double tval, norm = 0.0;

   for (buf=0; buf<buf_per_vect; buf++) {
      read(rootnum, buf);
      if (Parameters.hd_otf == FALSE) Hd.read(0, buf);
      else if (Parameters.hd_otf == TRUE) {
        if (Parameters.mpn) 
          Hd.diag_mat_els_otf(alplist, betlist, CalcInfo.onel_ints, 
             CalcInfo.twoel_ints, CalcInfo.e0_fzc, CalcInfo.num_alp_expl, 
             CalcInfo.num_bet_expl, CalcInfo.nmo, buf, Parameters.hd_ave);
        else 
          Hd.diag_mat_els_otf(alplist, betlist, CalcInfo.onel_ints, 
             CalcInfo.twoel_ints, CalcInfo.efzc, CalcInfo.num_alp_expl, 
             CalcInfo.num_bet_expl, CalcInfo.nmo, buf, Parameters.hd_ave);
        }

      if (Parameters.mpn) norm = calc_mpn_vec(buffer, lambda, Hd.buffer,
                           buf_size[buf], 1.0, -1.0, DIV);
      else {
        if (Parameters.precon >= PRECON_GEN_DAVIDSON) 
          h0block_gather_vec(CI_VEC);
        tval = calc_d2(buffer, lambda, Hd.buffer, buf_size[buf], precon);
       }

      if (buf_offdiag[buf]) tval *= 2.0;
      norm += tval;
      write(rootnum, buf);
      } 
   if (!Parameters.mpn) errcod = H0block_calc(lambda); /* MLL */
   return(norm);
}

/*
**
** CIvect::construct_kth_order_wf()
**
** Function constructs the kth order wavefunction from all 
** other nth order wavefunction (n<k) and the current Sigma vector
** Hc_k-1.  Uses only the two buffers.
**
** Parameters:
**    Hd         = CIvect for H0
**    S          = CIvect for Sigma vector
**    C          = CIvect for Cvec vector
**    alplist    = alpha string list
**    betlist    = beta  string list
**    buf1       = first buffer for a CIvect
**    buf2       = second buffer for a CIvect
**    
** Returns: none
*/
void CIvect::construct_kth_order_wf(CIvect &Hd, CIvect &S, CIvect &C,  
       struct stringwr **alplist, struct stringwr **betlist, double *buf1, 
       double *buf2, int k, double *mp_energy, double **cvec_coeff, 
       double *cvec_norm)
{
   int i, j, r, order, buf, block;
   double tval, norm;

   //fprintf(outfile,"\nCVEC_COEFF and CVEC_NORMS in CONSTRUCT\n");
   //print_mat(cvec_coeff, k-2, k-2, outfile);
   //for (i=0; i<k-2; i++) 
   //   fprintf(outfile,"cvec_norm[%d] = %lf\n",i, cvec_norm[i]);
   //fflush(outfile);

     for (buf=0; buf<buf_per_vect; buf++) {
        Hd.buf_lock(buf2);
        Hd.diag_mat_els_otf(alplist, betlist, CalcInfo.onel_ints,
             CalcInfo.twoel_ints, CalcInfo.e0_fzc, CalcInfo.num_alp_expl,
             CalcInfo.num_bet_expl, CalcInfo.nmo, buf, Parameters.hd_ave);
        read(k-1, buf);
        norm = calc_mpn_vec(buffer, (mp_energy[1]-CalcInfo.efzc), 
                Hd.buffer, buf_size[buf], 1.0, 1.0, MULT);
        Hd.buf_unlock();

        C.buf_lock(buf2);
        if (Parameters.mpn_schmidt) {
          for (i=0; i<=k-2; i++) {
             C.read(i, buf);
             tval = 0.0;
             for (r=2; r<=k; r++) { 
               if ((k-r)==i) tval+= mp_energy[r]*cvec_coeff[k-r][i]
                                    *(1.0/cvec_norm[k-r]);
               else tval+= mp_energy[r]*cvec_coeff[k-r][i];
               }
             xpeay(buffer, tval, C.buffer, buf_size[buf]); 
             }
          }
        else {
          for (i=2; i<=k; i++) {
             C.read(k-i, buf);
             xpeay(buffer, mp_energy[i], C.buffer, buf_size[buf]); 
             }
          }
        C.buf_unlock();

        S.buf_lock(buf2);
        S.read(0, buf); 
        xeaxmy(buffer, S.buffer, 1.0, S.buf_size[buf]);   
        S.buf_unlock();

        Hd.buf_lock(buf2);
        Hd.diag_mat_els_otf(alplist, betlist, CalcInfo.onel_ints,
             CalcInfo.twoel_ints, CalcInfo.e0_fzc, CalcInfo.num_alp_expl,
             CalcInfo.num_bet_expl, CalcInfo.nmo, buf, Parameters.hd_ave);
        norm = calc_mpn_vec(buffer, CalcInfo.e0, Hd.buffer, buf_size[buf],
                -1.0, 1.0, DIV);

        if (Ms0) {
          block = buf2blk[buf];  
          if ((int) Parameters.S % 2) symmetrize(-1.0, block);
          else symmetrize(1.0, block);
         }
        copy_zero_blocks(S);
        write(k, buf);
        Hd.buf_unlock();
        } 

}

/*
**
** CIvect::wigner_E2k_formula()
**
** Uses the kth order wavefunction and the Wigner formulas
**   to compute the 2k and 2k+1 th order energies
**
** Parameters:
**    Hd         = CIvect for H0
**    S          = CIvect for Sigma vector
**    C          = CIvect for Cvec vector
**    alplist    = alpha string list
**    betlist    = beta  string list
**    buf1       = first buffer for a CIvect
**    buf2       = second buffer for a CIvect
**    k          = kth order
**    mp2k_energy= array for storing the energy values 
**
** Returns: none
*/
void CIvect::wigner_E2k_formula(CIvect &Hd, CIvect &S, CIvect &C,
       struct stringwr **alplist, struct stringwr **betlist, double *buf1,
       double *buf2, int k, double *mp2k_energy, double **wfn_overlap,
       double **cvec_coeff, double *cvec_norm, int kvec_offset)
{
   int i, j, buf, I, J;
   double tval, E2k, E2kp1, tval2;

   E2k = E2kp1 = 0.0;

   /* First determine the overlap of kth order wavefunction with 
   ** all previous order wavefunctions
   */
   if (Parameters.mpn_schmidt) {
     zero_mat(wfn_overlap, Parameters.maxnvect+1, Parameters.maxnvect+1);
     /* for (i=0; i<k-1; i++) wfn_overlap[i][i] = 1.0; */
     }

   C.buf_lock(buf2);
   for (buf=0; buf<buf_per_vect; buf++) { 
      if (Parameters.mpn_schmidt) {
        for (i=0; i<=(k-kvec_offset); i++) {
           read(i, buf);
           for (j=i; j<=(k-kvec_offset); j++) {
              C.read(j,buf);
              dot_arr(buffer, C.buffer, C.buf_size[buf], &tval);
              if (buf_offdiag[buf]) tval *= 2.0;
              wfn_overlap[i+kvec_offset][j+kvec_offset] += tval;
              if (i!=j) wfn_overlap[j+kvec_offset][i+kvec_offset] += tval;
              }
           }
      /*
        read(k-1, buf);
        for (i=1; i<=k-1; i++) { 
           C.read(i, buf);          
           dot_arr(buffer, C.buffer, C.buf_size[buf], &tval);
           if (buf_offdiag[buf]) tval *= 2.0;
           wfn_overlap[k-1][i] += tval;
           if (i!=(k-1)) wfn_overlap[i][k-1] += tval;
           }
       */
        }
      else {
        read(k, buf);
        for (i=(1-kvec_offset); i<=(k-kvec_offset); i++) { 
           C.read(i, buf);          
           dot_arr(buffer, C.buffer, C.buf_size[buf], &tval);
           if (buf_offdiag[buf]) tval *= 2.0;
           wfn_overlap[k][i+kvec_offset] += tval;
           if ((i+kvec_offset)!=k) wfn_overlap[i+kvec_offset][k] += tval;
           }
        }
      }
   C.buf_unlock();

   
   if (Parameters.print_lvl > 3) {
     fprintf(outfile,"\nwfn_overlap = \n");
     print_mat(wfn_overlap, k+1, k+1, outfile);
     fprintf(outfile,"\t\t\t\t");
     }

   /* Compute E_2k and E_2k+1 */

   for (buf=0; buf<buf_per_vect; buf++) {
      S.buf_lock(buf2);
      S.read(0, buf);
      read(k-1-kvec_offset, buf);
      dot_arr(buffer, S.buffer, buf_size[buf], &tval);
      if (buf_offdiag[buf]) tval *= 2.0;
      E2k += tval;
      read(k-kvec_offset, buf);
      dot_arr(buffer, S.buffer, buf_size[buf], &tval);
      if (buf_offdiag[buf]) tval *= 2.0;
      E2kp1 += tval;
      S.buf_unlock();
      Hd.buf_lock(buf2);
      Hd.diag_mat_els_otf(alplist, betlist, CalcInfo.onel_ints,
           CalcInfo.twoel_ints, CalcInfo.e0_fzc, CalcInfo.num_alp_expl,
           CalcInfo.num_bet_expl, CalcInfo.nmo, buf, Parameters.hd_ave);
      xexy(Hd.buffer, buffer, buf_size[buf]);
      dot_arr(buffer, Hd.buffer, buf_size[buf], &tval);
      if (buf_offdiag[buf]) tval *= 2.0;
      E2kp1 -= tval;
      read(k-1-kvec_offset, buf);
      dot_arr(buffer, Hd.buffer, buf_size[buf], &tval);
      if (buf_offdiag[buf]) tval *= 2.0;
      E2k -= tval;      
      Hd.buf_unlock();
      }

   if (Parameters.mpn_schmidt) {
   /*
     C.buf_lock(buf2);
     for (i=1; i<=k-2; i++) {
        zero_arr(buffer, buf_size[0]);
        for (I=1; I<=k-2; I++) {
           C.read(I,0); 
           fprintf(outfile, " prescaled Bvec %d = \n", I);
           C.print(outfile); 
           if (i==I) tval = cvec_coeff[i][I]*(1.0/cvec_norm[i]);
           else tval = cvec_coeff[i][I];
           xpeay(buffer, tval, C.buffer, buf_size[0]); 
           }
        if (i==1) {
          fprintf(outfile, " Cvec %d = \n", i);
          print(outfile);
          }
        }
     C.buf_unlock();
    */        
     for (i=1-kvec_offset; i<=k-2-kvec_offset; i++) {
        for (j=1-kvec_offset; j<=k-2-kvec_offset; j++) {
           tval = 0.0;
           for (I=1-kvec_offset; I<=k-2-kvec_offset; I++) {
              if (I==j && I==i) tval += cvec_coeff[i+kvec_offset][I+kvec_offset]*
                                        cvec_coeff[j+kvec_offset][I+kvec_offset]*
                                        (1.0/cvec_norm[i+kvec_offset])*
                                        (1.0/cvec_norm[j+kvec_offset]);
              else if (I==i) tval += cvec_coeff[i+kvec_offset][I+kvec_offset]*
                                     cvec_coeff[j+kvec_offset][I+kvec_offset]*
                                     (1.0/cvec_norm[i+kvec_offset]); 
              else if (I==j) tval += cvec_coeff[i+kvec_offset][I+kvec_offset]*
                                     cvec_coeff[j+kvec_offset][I+kvec_offset]*
                                     (1.0/cvec_norm[j+kvec_offset]); 
              else tval += cvec_coeff[i+kvec_offset][I+kvec_offset]*
                           cvec_coeff[j+kvec_offset][I+kvec_offset];
              }
           E2k -= tval*mp2k_energy[2*k+kvec_offset-i+kvec_offset-j+kvec_offset];
           E2kp1 -= tval*mp2k_energy[2*k+1+kvec_offset-i+kvec_offset-j+kvec_offset];
           }
        tval = tval2 = 0.0;
        for (I=1-kvec_offset; I<=k-kvec_offset-2; I++) { 
           if (I==i) {
             tval += wfn_overlap[k][I]*cvec_coeff[i][I]*(1.0/cvec_norm[i]); 
             tval2 += wfn_overlap[k-1][I]*cvec_coeff[i][I]*(1.0/cvec_norm[i]); 
             }
           else { 
             tval += wfn_overlap[k][I]*cvec_coeff[i][I];
             tval2 += wfn_overlap[k-1][I]*cvec_coeff[i][I]; 
            }
           }
        E2kp1 -= tval*2.0*mp2k_energy[k+1-i]; 
        E2kp1 -= tval2*2.0*mp2k_energy[k+2-i]; 
        E2k -= tval*mp2k_energy[k-i];
        E2k -= tval2*2.0*mp2k_energy[k+1-i];
       /*
        fprintf(outfile, "E2kp1 -> - tval*2.0*mp2k_energy[k+1-i] = %20.10f\n", 
                tval*2.0*mp2k_energy[k+1-i]);
        fprintf(outfile, "E2kp1 -> - tval2*2.0*mp2k_energy[k+2-i] = %20.10f\n", 
                tval*2.0*mp2k_energy[k+2-i]);
        fprintf(outfile, "E2k -> - tval*mp2k_energy[k-i] = %20.10f\n", 
                tval*mp2k_energy[k-i]);
        fprintf(outfile, "E2k -> - tval*2.0*mp2k_energy[k+1-i] = %20.10f\n", 
                tval*2.0*mp2k_energy[k+1-i]);
        */
        }
     E2kp1 += (CalcInfo.efzc-mp2k_energy[1])*wfn_overlap[k][k];
     E2kp1 -= 2.0*mp2k_energy[2]*wfn_overlap[k-1][k];
     E2kp1 -= mp2k_energy[3]*wfn_overlap[k-1][k-1];
     E2k += (CalcInfo.efzc-mp2k_energy[1])*wfn_overlap[k][k-1]; 
     E2k -= mp2k_energy[2]*wfn_overlap[k-1][k-1];
     }

  /*
   else {
     C.buf_lock(buf2);
     for (I=1; I<=k; I++) {
        C.read(I,0);
        if (I==1) {
          fprintf(outfile,"Cvec %d = \n", I);
          C.print(outfile);
          }
        }
     C.buf_unlock();
     for (i=1; i<=k-2; i++) {
        E2kp1 -= 2.0 * mp2k_energy[k+1-i] * wfn_overlap[i][k];
        E2kp1 -= 2.0 * mp2k_energy[k+2-i] * wfn_overlap[i][k-1]; 
        for (j=1; j<=k-2; j++) 
           E2kp1 -= mp2k_energy[2*k+1-i-j] * wfn_overlap[i][j];
        }

     for (i=1; i<=k-2; i++) {
        E2k -= mp2k_energy[k-i] * wfn_overlap[k][i];
        E2k -= 2.0 * mp2k_energy[k+1-i] * wfn_overlap[k-1][i];
        for (j=1; j<=k-2; j++) 
           E2k -= mp2k_energy[2*k-i-j] * wfn_overlap[i][j];
        }
     E2k += (CalcInfo.efzc-mp2k_energy[1]) * wfn_overlap[k][k-1];
     E2k -= mp2k_energy[2] * wfn_overlap[k-1][k-1];
     E2kp1 += (CalcInfo.efzc-mp2k_energy[1])*wfn_overlap[k][k];
     E2kp1 -= mp2k_energy[3]*wfn_overlap[k-1][k-1];
     E2kp1 -= 2.0 * mp2k_energy[2] * wfn_overlap[k-1][k];
     }
  */

   else {
     for (i=1; i<=k; i++) 
        for (j=1; j<=k; j++) {
           E2kp1 -= mp2k_energy[2*k+1-i-j] * wfn_overlap[i][j];
           if ((i==k) && (j==k)) E2kp1 += CalcInfo.efzc * wfn_overlap[k][k];
           }

     for (i=1; i<=k; i++)
        for (j=1; j<k; j++) {
           E2k -= mp2k_energy[2*k-i-j] * wfn_overlap[i][j];
           if ((i==k) && (j==k-1)) E2k += CalcInfo.efzc * wfn_overlap[k][k-1];
           }
     }

  /*  fprintf(outfile, "final E2k = %lf\n", E2k);
   fprintf(outfile, "final E2kp1 = %lf\n", E2kp1);
  */

   mp2k_energy[2*k] = E2k;
   mp2k_energy[2*k+1] = E2kp1;

}


/*
** CIvect::print_buf()
**
** Function prints the in-core section of a CI vector
**
** Parameters: 
**    outfile  = file pointer for output file
**
** Returns: none
*/
void CIvect::print_buf(FILE *outfile)
{
   int blk;
   int irrep; 

   if (icore == 1) {
      for (blk=0; blk<num_blocks; blk++) {
         fprintf(outfile, "\nBlock %2d, codes = (%2d,%2d)\n", blk, 
            Ia_code[blk], Ib_code[blk]);
         print_mat(blocks[blk], Ia_size[blk], Ib_size[blk], outfile);
         }
      }

   if (icore == 2) { /* symmetry block in-core */
      irrep = buf2blk[cur_buf];
      if (first_ablk[irrep] < 0) {
         fprintf(outfile, "(CIvect::print_blk): No blks for irrep %d\n",irrep);
         return;
         }
      else {
         for (blk=first_ablk[irrep]; blk <= last_ablk[irrep]; blk++) {
            fprintf(outfile, "\nBlock %2d, codes = (%2d,%2d)\n", blk,
               Ia_code[blk], Ib_code[blk]);
            print_mat(blocks[blk], Ia_size[blk], Ib_size[blk], outfile);
            }
         }
      }

   if (icore == 0) { /* one subblock in-core */
      blk = buf2blk[cur_buf];

      fprintf(outfile, "\nBlock %2d, codes = (%2d,%2d)\n", blk,
         Ia_code[blk], Ib_code[blk]);
      print_mat(blocks[blk], Ia_size[blk], Ib_size[blk], outfile);
      }
}



/*
** CIvect::civ_xeay()
**
** Function does the operation X = a * Y for two CI vectors X and Y and
** some constant a.  
**
** Parameters:
**    a     =  constant in X = a * Y
**    Y     =  vector in X = a * Y
**    xvect = vector number in X  
**    yvect = vector number in Y  
** 
** Returns: none
*/
void CIvect::civ_xeay(double a, CIvect &Y, int xvect, int yvect)
{
   int buf;

   for (buf=0; buf<buf_per_vect; buf++) {
      Y.read(yvect, buf);
      xeay(buffer, a, Y.buffer, buf_size[buf]);
      write(xvect, buf);
      }
}



/*
** CIvect::civ_xpeay()
**
** Function does the operation X += a * Y for two CI vectors X and Y and
** some constant a.
**
** Parameters:
**    a     =  constant in X = a * Y
**    Y     =  vector in X = a * Y
**    xvect = vector number in X
**    yvect = vector number in Y
**
** Returns: none
*/
void CIvect::civ_xpeay(double a, CIvect &Y, int xvect, int yvect)
{
   int buf;

   for (buf=0; buf<buf_per_vect; buf++) {
      Y.read(yvect, buf);
      read(xvect, buf);
      xpeay(buffer, a, Y.buffer, buf_size[buf]);
      write(xvect, buf);
      }
}


/*
** CIvect::transp_block()
**
** Transpose a CI vector (or a piece of a CI vector)
**
** Parameters:
**    iblock = RAS subblock number
**    tmparr = scratch array to use as intermediate
**
** Returns: none
*/
void CIvect::transp_block(int iblock, double **tmparr)
{
   int i,j,nrows,ncols;
   double **src, *dest;

   src = blocks[iblock];
   dest = tmparr[0];

   /* bind pointers to subbblock topology */
   nrows = Ib_size[iblock];
   ncols = Ia_size[iblock];
   for (i=1; i<nrows; i++) {
      tmparr[i] = dest + i * ncols;
      }

   /* copy data */
   for (i=0; i<Ib_size[iblock]; i++) {
      for (j=0; j<Ia_size[iblock]; j++) {
         *dest++ = src[j][i];
         }
      }

}


/*
** CIvect::get_max_blk_size()
**
** Return the maximum RAS subblock size as a long unsigned integer
**
*/
unsigned long CIvect::get_max_blk_size(void)
{
   int i;
   unsigned long blksize, maxblksize=0;

   for (i=0; i<num_blocks; i++) {
      blksize = (unsigned long) Ia_size[i] * (unsigned long) Ib_size[i];
      if (blksize > maxblksize) maxblksize = blksize;
      } 

   return(maxblksize);
}


/*
** CIvect::checknorm()
**
** Check the norm of a CI vector
*/
double CIvect::checknorm(void)
{
   double tval, dotprod = 0.0;
   int buf;

   for (buf=0; buf<buf_per_vect; buf++) {
      read(cur_vect, buf);
      dot_arr(buffer, buffer, buf_size[buf], &tval);
      if (buf_offdiag[buf]) tval *= 2.0;
      dotprod += tval;
      }

   return(dotprod);
}


/*
** CIvect::copy()
** 
** This copies one CI vector to another
**
*/
void CIvect::copy(CIvect &Src, int targetvec, int srcvec)
{
   int buf, blk;

   for (buf=0; buf<buf_per_vect; buf++) {
      Src.read(srcvec, buf);
      xey(buffer, Src.buffer, buf_size[buf]);
      blk = buf2blk[buf];
      if ( (zero_blocks[blk]==0) || (Src.zero_blocks[blk]==0) )
         zero_blocks[blk] = 0;
      write(targetvec, buf);
      } 
}


/*
** CIvect::restart_gather()
**
** This function takes a linear combination of previous 'b' vectors to
** make a current 'c' vector during a Davidson iteration procedure.
** The coefficients are given by array alpha, the current vector number
** is given by ivec, and the number of previous vectors is nvec.
**
*/
void CIvect::restart_gather(int ivec, int nvec, int nroot, double **alpha,
      double *buffer1, double *buffer2)
{
   int buf, oldvec;

   for (buf=0; buf<buf_per_vect; buf++) {
      zero_arr(buffer2, buf_size[buf]);
      buf_lock(buffer1);
      for (oldvec=0; oldvec<nvec; oldvec++) {
         read(oldvec, buf);
         xpeay(buffer2, alpha[oldvec][nroot], buffer1, buf_size[buf]);
         }
      buf_unlock();
      buf_lock(buffer2);
      write(ivec, buf);
      buf_unlock();
      } 
}



/*
** CIvect::gather()
**
** This function takes a linear combination of previous 'b' vectors to
** make a current 'c' vector during a Davidson iteration procedure.
** The coefficients are given by array alpha, the current vector number
** is given by ivec, and the number of previous vectors is nvec.
** Similar to CIvect::restart_gather() except that we no longer assume
** two different CIvects.  Assume buffer locks are already done.
**
*/
void CIvect::gather(int ivec, int nvec, int nroot, double **alpha,
      CIvect &C)
{
   int buf, oldvec;

   /* fprintf(outfile,"In CIvect::gather\n"); */ 
   for (buf=0; buf<buf_per_vect; buf++) {
      zero_arr(buffer, buf_size[buf]);
      for (oldvec=0; oldvec<nvec; oldvec++) {
         C.read(oldvec, buf);
         xpeay(buffer, alpha[oldvec][nroot], C.buffer, buf_size[buf]);
         /* fprintf(outfile,"coef[%d][%d] = %10.7f\n",oldvec,nroot,alpha[oldvec][nroot]); */ 
         }
      write(ivec, buf);
      } 
}


/*
** CIvect::restart_reord_fp()
**
** This function reorders the file pointers to do a restart of the iterative
** diagonalization method.  The idea is that, for nroot > 1, it is not 
** possible to restart and make CI vectors 0...nroot equal to the restarted 
** approximate eigenvectors, since this overwrites info needed to construct
** them.  Thus, the new CI vectors must occupy the LAST nroot positions.
** However, it is nevertheless useful to index them as 0...nroot.  This
** is most easily accomplished by a remapping (rotation) of the file
** pointer info.  The one parameter is L, the new "0" vector.
**
** Actually, it's slightly more complex.  For multiple restarts in a given
** calc, the 0 position rotates around.  This routine should still work.
**
** In the latest version, I am phasing out the "offset" array in favor
** of a new_first_buf array which basically gives the buffer number of
** the new "0" vector.  This is more natural for the libpsio implementation.
*/
void CIvect::restart_reord_fp(int L)
{
   int buf, newbuf;
   int *tmp_file_number;

   new_first_buf = L*buf_per_vect + new_first_buf;
   if (new_first_buf >= buf_total) new_first_buf -= buf_total;

   /*
   tmp_file_offset = (unsigned long *) malloc (buf_total * 
                      sizeof(unsigned long));
   tmp_file_number = init_int_array(buf_total);


   for (buf=L*buf_per_vect,newbuf=0; buf<buf_total; buf++,newbuf++) {
     tmp_file_offset[newbuf] = file_offset[buf];
     tmp_file_number[newbuf] = file_number[buf];
   }

   for (buf=0; buf<L*buf_per_vect; buf++,newbuf++) { 
     tmp_file_offset[newbuf] = file_offset[buf];
     tmp_file_number[newbuf] = file_number[buf];
   }

   for (buf=0; buf<buf_total; buf++) {
     file_offset[buf] = tmp_file_offset[buf];
     file_number[buf] = tmp_file_number[buf];
   }

   free(tmp_file_offset);
   free(tmp_file_number);
   */
}


void CIvect::print_fptrs()
{
   int buf;

   fprintf(outfile,"Printing file pointer information\n");
   for (buf=0; buf<buf_total; buf++)  
      fprintf(outfile,"%d %d\n",buf, file_number[buf]);
      
    
}
void sigma_get_contrib(struct stringwr **alplist, struct stringwr **betlist, 
      CIvect &C, CIvect &S, int **s1_contrib, int **s2_contrib,
      int **s3_contrib)
{

   int sblock,cblock;
   int sac, sbc, cac, cbc;
   int nas, nbs;
   struct stringwr *Ib, *Ia, *Kb, *Ka;
   unsigned int Ibidx, Iaidx, Kbidx, Kaidx, Ib_ex, Ia_ex;
   unsigned int Ibcnt, Iacnt, *Ibridx, *Iaridx; 
   int Kb_list, Ka_list;
   int found,i,j;
    
   for (sblock=0; sblock<S.num_blocks; sblock++) {
      sac = S.Ia_code[sblock];
      sbc = S.Ib_code[sblock];
      nas = S.Ia_size[sblock];
      nbs = S.Ib_size[sblock];
      for (cblock=0; cblock<C.num_blocks; cblock++) {
         cac = C.Ia_code[cblock];
         cbc = C.Ib_code[cblock];

     
         /* does this c block contribute to sigma1? */
         if (sac == cac) {
            for (Ib=betlist[sbc], Ibidx=0, found=0; Ibidx < nbs && !found; 
               Ibidx++, Ib++) {
               /* loop over excitations E^b_{kl} from |B(I_b)> */ 
               for (Kb_list=0; Kb_list < S.num_betcodes && !found; Kb_list++) {
                  Ibcnt = Ib->cnt[Kb_list];
                  Ibridx = Ib->ridx[Kb_list];
                  for (Ib_ex=0; Ib_ex < Ibcnt; Ib_ex++) {
                     Kbidx = *Ibridx++;
                     Kb = betlist[Kb_list] + Kbidx;
                     if (Kb->cnt[cbc]) { found=1;  break; }
                     }
                  }
               }
            if (found) s1_contrib[sblock][cblock] = 1;
            }       

         /* does this c block contribute to sigma2? */
         if (sbc == cbc) {
         for (Ia=alplist[sac], Iaidx=0, found=0; Iaidx < nas && !found;
               Iaidx++, Ia++) {
               /* loop over excitations E^a_{kl} from |A(I_a)> */
               for (Ka_list=0; Ka_list < S.num_alpcodes && !found; Ka_list++) {
                  Iacnt = Ia->cnt[Ka_list];
                  Iaridx = Ia->ridx[Ka_list];
                  for (Ia_ex=0; Ia_ex < Iacnt; Ia_ex++) {
                     Kaidx = *Iaridx++;
                     Ka = alplist[Ka_list] + Kaidx;
                     if (Ka->cnt[cac]) { found=1;  break; }
                     }
                  }
               }
            if (found) s2_contrib[sblock][cblock] = 1;
            }

         /* does this c block contribute to sigma3? */
         for (Iaidx=0,found=0; Iaidx<S.Ia_size[sblock]; Iaidx++) {
            if (alplist[sac][Iaidx].cnt[cac]) found=1;
            } 
         if (found) { /* see if beta is ok */
            found=0;
            for (Ibidx=0; Ibidx<S.Ib_size[sblock]; Ibidx++) {
               if (betlist[sbc][Ibidx].cnt[cbc]) found=1;
               }
            if (found)
               s3_contrib[sblock][cblock] = 1;
            }        

         } /* end loop over c blocks */
      } /* end loop over sigma blocks */

   #ifdef DEBUG
   printf("\nSigma 1:\n");
   for (i=0; i<S.num_blocks; i++) {
     fprintf(outfile, "Contributions to sigma block %d\n", i);
     for (j=0; j<C.num_blocks; j++) {
       if (s1_contrib[i][j]) fprintf(outfile, "%3d ", j);
     }
     fprintf(outfile, "\n");
   }

   printf("\n\nSigma 2:\n");
   for (i=0; i<S.num_blocks; i++) {
     fprintf(outfile, "Contributions to sigma block %d\n", i);
     for (j=0; j<C.num_blocks; j++) {
       if (s2_contrib[i][j]) fprintf(outfile, "%3d ", j);
     }
     fprintf(outfile, "\n");
   }

   printf("\n\nSigma 3:\n");
   for (i=0; i<S.num_blocks; i++) {
     fprintf(outfile, "Contributions to sigma block %d\n", i);
     for (j=0; j<C.num_blocks; j++) {
       if (s3_contrib[i][j]) fprintf(outfile, "%3d ", j);
     }
     fprintf(outfile, "\n");
   }
   #endif

}


/*
** olsen_iter_xy(): This function evaluates the quantities x and y defined by
**    x = C^(i) * (Hd - E)^-1 * C^(i)
**    y = C^(i) * (Hd - E)^-1 * sigma^(i)
**    The diagonal elements of H (Hd) are used in this routine.
**    These quantities are subtracted out and substituted in later
**    routines if the gen_davidson or h0block_inv preconditioners
**    are used. A very intelligent approach originally implimented 
**    by a fellow Jedi.
**
** Parameters:
**   C       =  reference of current iteration's ci vector
**   S       =  reference of current iteration's sigma vector
**   Hd      =  reference of vector of diagonal elements of H
**   x       =  pointer to double to hold x
**   y       =  pointer to double to hold y 
**   buffer1 = pointer to first I/O buffer
**   buffer2 = pointer to second I/O buffer
**   E       =  current iteration's energy
**   curvect = current vector
**   L       = number of vectors in b or sigma files
**   alplist = alpha string list for use with diag energies on the fly
**   betlist = beta string list for use with diag energies on the fly
**
** Returns: none
*/
void olsen_iter_xy(CIvect &C, CIvect &S, CIvect &Hd, double *x, double *y,
      double *buffer1, double *buffer2, double E, int curvect, int L, 
      double **alpha, struct stringwr **alplist, struct stringwr **betlist)
{

   int buf, i,j;
   double tx = 0.0, ty = 0.0, tmpy = 0.0;
   double *sigma0b1, *sigma0b2;
   *x = 0.0;
   *y = 0.0;
 
   Hd.buf_lock(buffer2);
   if (Parameters.diag_method==METHOD_DAVIDSON_LIU_SEM) {
     sigma0b1 = init_array(H0block.size);
     sigma0b2 = init_array(H0block.size);
     }
   for (buf=0; buf<C.buf_per_vect; buf++) {
      tx = ty = 0.0;
      C.buf_lock(buffer1);
      C.read(curvect,buf);
      if (Parameters.diag_method==METHOD_DAVIDSON_LIU_SEM) 
        C.h0block_gather_vec(CI_VEC); 
      if (Parameters.hd_otf == FALSE) Hd.read(0,buf);
      else Hd.diag_mat_els_otf(alplist, betlist, CalcInfo.onel_ints,
           CalcInfo.twoel_ints, CalcInfo.efzc, CalcInfo.num_alp_expl,
           CalcInfo.num_bet_expl, CalcInfo.nmo, buf, Parameters.hd_ave);
      tx = buf_xy1(buffer1, buffer2, E, Hd.buf_size[buf]);
      /* buffer2 = Hd * Ci */
      C.buf_unlock();
      S.buf_lock(buffer1);
      if (Parameters.diag_method <= METHOD_MITRUSHENKOV) { 
        /* Olsen and Mitrushenkov iterators */
        S.read(curvect,buf);
        dot_arr(buffer1, buffer2, C.buf_size[buf], &ty);
        }
      else { /* Dot buffer2 with all Sigma vectors on disk */
        for (i=0; i<L; i++) {
           S.read(i,buf);
           dot_arr(buffer1, buffer2, C.buf_size[buf], &tmpy); 
           ty += tmpy * alpha[i][curvect]; 
           zero_arr(sigma0b1,H0block.size); 
           S.h0block_gather_multivec(sigma0b1);
           for (j=0; j<H0block.size; j++)
              sigma0b2[j] += alpha[i][curvect] * sigma0b1[j]; 
           }
        
       } 
      if (C.buf_offdiag[buf]) {
         *x += 2.0 * tx;
         *y += 2.0 * ty;
         }
      else {
         *x += tx;
         *y += ty;
         }
      S.buf_unlock();
      }

   Hd.buf_unlock();
   if (Parameters.diag_method==METHOD_DAVIDSON_LIU_SEM) {
     for (j=0; j<H0block.size; j++)
        H0block.s0b[j] = sigma0b2[j]; 
     free(sigma0b1);
     free(sigma0b2);
     }
           
}


/*
** olsen_update()
**
** Update the CI vector according to the Olsen method.
** Try to assume storage of only one C, one S, and Hd on disk.  Let
** the mitrush_update() function assume two C and two S.  No, can't 
** assume that or we'll never get Mitrush started.  Just fix it so that
** Cn can be identical to Cc.
** The update formula is c0' = c0 + (Ho-Eo)^-1(E_est*Co - sigma_o) 
** this is the same formula as the davidson procedure except lambda_i 
** has been substituted by E_est 
*/
void olsen_update(CIvect &C, CIvect &S, CIvect &Hd, double E, double E_est, 
      double *norm, double *c1norm, double *ovrlap, double *buffer1, 
      double *buffer2, 
      int curr, int next, FILE *outfile, int iter, struct stringwr **alplist,
      struct stringwr **betlist)
{

   int buf;
   double nx=0.0, ox=0.0, tmp1, tmp2, normc1=0.0, tmpnorm=0.0;
   double rnorm=0.0, rnormtmp=0.0;

   for (buf=0; buf<C.buf_per_vect; buf++) {
      tmp1 = 0.0;
      tmp2 = 0.0;
      C.buf_lock(buffer1);
      S.buf_lock(buffer2);
      C.read(curr, buf);
      S.read(curr, buf);
      /* C = E_est * C - S, C is buffer1*/
      xeaxmy(buffer1, buffer2, E_est, C.buf_size[buf]);
      C.buf_unlock();
      S.buf_unlock();
      Hd.buf_lock(buffer2);
      if (Parameters.hd_otf == FALSE) Hd.read(0,buf);
      else Hd.diag_mat_els_otf(alplist, betlist, CalcInfo.onel_ints,
           CalcInfo.twoel_ints, CalcInfo.efzc, CalcInfo.num_alp_expl,
           CalcInfo.num_bet_expl, CalcInfo.nmo, buf, Parameters.hd_ave);
      /* Check norm of residual vector i.e. before preconditioning */
      dot_arr(buffer1, buffer1, C.buf_size[buf], &rnormtmp);
      /* C = C/(Hd - E) */
      buf_ols_denom(buffer1, buffer2, E, S.buf_size[buf]); 
      /* buffer1 is now equal to C^1, i.e. the correction to C_i 
      ** without the H0block correction
      */ 
      Hd.buf_unlock();
      /* C_new = C_i + C^1 */
      C.buf_lock(buffer2);
      C.read(curr, buf);
      buf_ols_updt(buffer1,buffer2,&tmp1,&tmp2,&tmpnorm,C.buf_size[buf],outfile);
      if (Parameters.precon >= PRECON_GEN_DAVIDSON) 
        C.h0block_buf_ols(&tmp1,&tmp2,&tmpnorm,E_est); 
      if (C.buf_offdiag[buf]) {
         tmp1 *= 2.0;
         tmp2 *= 2.0;
         tmpnorm *= 2.0;
         rnormtmp *= 2.0;
         }
      normc1 += tmpnorm;
      nx += tmp1;
      ox += tmp2;
      rnorm += rnormtmp;
      C.write(next, buf);
      C.buf_unlock();
 } 

   *norm = nx;
   /* fprintf(outfile,"\n ovrlap(ox) = %20.16f\n", ox); */
   *ovrlap = ox; 
   if (normc1 <= 1.0E-13) { 
     fprintf(outfile,"Norm of correction vector = %5.4e\n", normc1);
     fprintf(outfile,"This may cause numerical errors which would" \
         " deteriorate the diagonalization procedure.\n");
     }
   *c1norm = sqrt(rnorm); 
   normc1 = sqrt(normc1); 
}


/*
** CIvect::h0block_buf_init()
**
** Initialize H0block stuff pertaining to buffers
**
*/
void CIvect::h0block_buf_init(void)
{
   int i, cnt, irrep, buf, blk;
   int *tmparr;

   H0block.nbuf = buf_per_vect;
   H0block.buf_num = init_int_array(buf_per_vect);
   if (H0block.size < 1) return;

   tmparr = init_int_array(H0block.size+H0block.coupling_size);

   if (icore == 1) {
      H0block.buf_member = 
        init_int_matrix(1, H0block.size+H0block.coupling_size);
      for (i=0; i<(H0block.size+H0block.coupling_size); i++) {
         H0block.buf_member[0][i] = i;
         } 
      H0block.buf_num[0] = H0block.size+H0block.coupling_size;
      }
   else if (icore == 2) {
      H0block.buf_member = (int **) malloc (buf_per_vect * sizeof(int *));
      for (buf=0; buf<buf_per_vect; buf++) {
         cnt = 0;
         irrep = buf2blk[buf];
         for (blk=first_ablk[irrep]; blk<=last_ablk[irrep]; blk++) {
            for (i=0; i<H0block.size+H0block.coupling_size; i++) {
               if (H0block.blknum[i] == blk) {
                  tmparr[cnt++] = i;
                  }
               }
            } 
         H0block.buf_num[buf] = cnt;
         if (cnt) H0block.buf_member[buf] = init_int_array(cnt);
         for (i=0; i<cnt; i++) {
            H0block.buf_member[buf][i] = tmparr[i];
            }
         } /* end loop over bufs */
      } /* end icore==2 */
   else { 
      H0block.buf_member = (int **) malloc (buf_per_vect * sizeof(int *));
      for (buf=0; buf<buf_per_vect; buf++) {
         cnt = 0; 
         blk = buf2blk[buf];
         for (i=0; i<H0block.size+H0block.coupling_size; i++) {
            if (H0block.blknum[i] == blk) {
               tmparr[cnt++] = i;
               }
            }
         H0block.buf_num[buf] = cnt;
         if (cnt) H0block.buf_member[buf] = init_int_array(cnt);
         for (i=0; i<cnt; i++) {
            H0block.buf_member[buf][i] = tmparr[i];
            }
         } /* end loop over bufs */
      } /* end icore==0 */

   free(tmparr);
}



/*
** CIvect::h0block_buf_ols()
**
** Parameters:
**    norm   = block's norm 
**    ovrlap = blocks' overlap 
*/
void CIvect::h0block_buf_ols(double *nx, double *ox, double *c1norm,double E_est)
{
   int i, j, k, blk, al, bl;
   double c, cn, tval, c1;

   for (i=0; i<H0block.buf_num[cur_buf]; i++) {
      j = H0block.buf_member[cur_buf][i];
      blk = H0block.blknum[j];
      al = H0block.alpidx[j];
      bl = H0block.betidx[j];
      c = H0block.c0b[j];
      cn = blocks[blk][al][bl]; 
      c1 = cn - c;
      *nx -= cn * cn;
      *ox -= cn * c;
      *c1norm -= c1 * c1;
      tval = c + E_est * H0block.c0bp[j];
      tval -= H0block.s0bp[j];
      blocks[blk][al][bl] = tval;
      /* H0block.c0b[j] = tval; */ /* this should gather all of c0b. Norm later */
                              /* What about the effect of symmetrization? */
                              /* symmetrization was the bug MLL */
     /*
      if (buf_offdiag[cur_buf]) {
         k = H0block.pair[j];
         if (k >= 0 && k != j) {
            H0block.c0b[k] = tval * phase;
            }
         }
      */
      *nx += tval * tval;
      *ox += tval * c;
      *c1norm += (tval - c) * (tval - c);
      }
}
       
/*
** CIvect::h0block_gather_vec(int vecode)
**
** Parameters:
**    curr = current vector number
**    vecode = 0 for C vector and 1 for Sigma vector
*/
void CIvect::h0block_gather_vec(int vecode)
{
   int buf, i, j, k, blk, al, bl;
   double c, cn, tval, phase, norm = 0.0;

   if (!Parameters.Ms0) phase = 1.0;
   else phase = ((int) Parameters.S % 2) ? -1.0 : 1.0;

   for (i=0; i<H0block.buf_num[cur_buf]; i++) {
      j = H0block.buf_member[cur_buf][i];
      blk = H0block.blknum[j];
      al = H0block.alpidx[j];
      bl = H0block.betidx[j];
      tval = blocks[blk][al][bl];
      if (vecode) H0block.s0b[j] = tval;
      else H0block.c0b[j] = tval;
      if (buf_offdiag[cur_buf]) {
        k = H0block.pair[j];
        if (k >= 0 && k != j) {
        /* if (k >= 0 && k != j && Parameters.Ms0)  */
          if (vecode) H0block.s0b[k] = tval * phase;
          else H0block.c0b[k] = tval * phase;
           }
        }
      }
    /*
     if (!vecode) {
       fprintf(outfile,"c0b in h0block_gather_vec = \n");
       print_mat(&(H0block.c0b), 1, H0block.size, outfile);
       }
    */
}


/*
** CIvect::h0block_gather_multivec(double *vec)
**
** Parameters:
**    curr = current vector number
**    vecode = 0 for C vector and 1 for Sigma vector
*/
void CIvect::h0block_gather_multivec(double *vec)
{
   int buf, i, j, k, blk, al, bl;
   double c, cn, tval, phase, norm = 0.0;

   if (!Parameters.Ms0) phase = 1.0;
   else phase = ((int) Parameters.S % 2) ? -1.0 : 1.0;

   for (i=0; i<H0block.buf_num[cur_buf]; i++) {
      j = H0block.buf_member[cur_buf][i];
      blk = H0block.blknum[j];
      al = H0block.alpidx[j];
      bl = H0block.betidx[j];
      tval = blocks[blk][al][bl];
      vec[j] = tval;
      if (buf_offdiag[cur_buf]) {
        k = H0block.pair[j];
        if (k >= 0 && k != j) {
        /* if (k >= 0 && k != j && Parameters.Ms0)  */
          vec[k] = tval * phase;
           }
        }
      }
}


/*
** CIvect::h0block_buf_precon(double *nx, int root)
**
** Routine used by sem
** Parameters:
**    norm   = block's norm
*/
void CIvect::h0block_buf_precon(double *nx, int root)
{
   int i, j, k, blk, al, bl, buf;
   double c, cn, tval, phase, norm = 0.0;

   if (!Parameters.Ms0) phase = 1.0;
   else phase = ((int) Parameters.S % 2) ? -1.0 : 1.0;

   for (buf=0; buf<buf_per_vect; buf++) {
      read(root,buf);
      for (i=0; i<H0block.buf_num[buf]; i++) {
         j = H0block.buf_member[buf][i];
         blk = H0block.blknum[j];
         al = H0block.alpidx[j];
         bl = H0block.betidx[j];
         tval = blocks[blk][al][bl] * blocks[blk][al][bl];
         *nx -= tval;
         if (buf_offdiag[buf]) {
           k = H0block.pair[j];
           if (k >= 0 && k!=j) *nx -= tval * phase;
           }
         tval = H0block.c0bp[j] * H0block.c0bp[j];
         *nx += tval;
         if (buf_offdiag[buf]) {
           k = H0block.pair[j];
           if (k>= 0 && k!=j) *nx += tval * phase; 
           }
         blocks[blk][al][bl] = -H0block.c0bp[j];
         }
      write(root,buf);
      }
}

/*
** mitrush_update()
** Perform the Mitrushenkov update.  New version 3/96
**
*/
void mitrush_update(CIvect &C, CIvect &S, double norm, double acur, 
   double alast, double *buffer1, double *buffer2, int curr, int next)
{
   int i, j, k, buf, blk, al, bl;
   double phase, tval; 

   if (!Parameters.Ms0) phase = 1.0;
   else phase = ((int) Parameters.S % 2) ? -1.0 : 1.0;

   for (buf=0; buf<C.buf_per_vect; buf++) {
      C.buf_lock(buffer1);
      C.read(curr, buf);
      C.buf_unlock();
      C.buf_lock(buffer2);
      C.read(next, buf);
      xeaxpby(buffer2, buffer1, alast, acur, C.buf_size[buf]);
      /* buffer2 is the new C(i) vector constructed from 
      ** alpha(i-1)*C(i-1) + alpha(i)*C(i).
      ** However, this vector is not normalized or symmetrized yet.
      ** Still need to construct the olsen_update of this new
      ** C(i) vector.
      */ 
      C.write(curr, buf);
      C.buf_unlock();
      }
   C.buf_lock(buffer1);
   C.read(curr,0);
   C.symnorm(norm,0,1); 
   C.buf_unlock();

   for (buf=0; buf<S.buf_per_vect; buf++) {
      S.buf_lock(buffer1);
      S.read(curr, buf);
      S.buf_unlock();
      S.buf_lock(buffer2);
      S.read(next, buf);
      xeaxpby(buffer2, buffer1, alast, acur, S.buf_size[buf]);
      S.write(curr, buf);
      S.buf_unlock();
      }
   S.buf_lock(buffer1);
   S.read(curr,0);
   S.symnorm(norm,1,1);
   S.buf_unlock();
}

void sigma_get_contrib_rotf(CIvect &C, CIvect &S, 
      int **s1_contrib, int **s2_contrib, int **s3_contrib, 
      int *Cnt[2], int **Ij[2], int **Oij[2], int **Ridx[2],
      signed char **Sgn[2], unsigned char **Toccs)
{

   int sblock,cblock;
   int sac, sbc, cac, cbc;
   int nas, nbs;
   int Ibidx, Iaidx, Ib_ex, Ia_ex;
   int Ibcnt, Iacnt; 
   int Kb_list, Ka_list;
   int found,i,j;
    
   for (sblock=0; sblock<S.num_blocks; sblock++) {
      sac = S.Ia_code[sblock];
      sbc = S.Ib_code[sblock];
      nas = S.Ia_size[sblock];
      nbs = S.Ib_size[sblock];
      for (cblock=0; cblock<C.num_blocks; cblock++) {
         cac = C.Ia_code[cblock];
         cbc = C.Ib_code[cblock];

     
         /* does this c block contribute to sigma1? */
         if (sac == cac) {
            found = 0;
            for (Kb_list=0; Kb_list < S.num_betcodes && !found; Kb_list++) {
               b2brepl(Occs[sbc], Cnt[0], Ij[0], Oij[0], Ridx[0],
                  Sgn[0], BetaG, sbc, Kb_list, nbs);
               for (Ibidx=0; Ibidx < nbs && !found; Ibidx++) {
                  Ibcnt = Cnt[0][Ibidx];
                  if (Ibcnt) {
                     for (i=0; i<Ibcnt; i++) {
                        j = Ridx[0][Ibidx][i];
                        Toccs[i] = Occs[Kb_list][j];
                        }
                     b2brepl(Toccs, Cnt[1], Ij[1], Oij[1], Ridx[1], Sgn[1],
                        BetaG, Kb_list, cbc, Ibcnt);
                     for (Ib_ex=0; Ib_ex < Ibcnt; Ib_ex++) {
                        if (Cnt[1][Ib_ex]) { found=1; break; }
                        } 
                     }
                  }
               }
            if (found) s1_contrib[sblock][cblock] = 1;
            }

         /* does this c block contribute to sigma2? */
         if (sbc == cbc) {
            found = 0;
            for (Ka_list=0; Ka_list < S.num_alpcodes && !found; Ka_list++) {
               b2brepl(Occs[sac], Cnt[0], Ij[0], Oij[0], Ridx[0],
                  Sgn[0], AlphaG, sac, Ka_list, nas);
               for (Iaidx=0; Iaidx < nas && !found; Iaidx++) {
                  Iacnt = Cnt[0][Iaidx];
                  if (Iacnt) {
                     for (i=0; i<Iacnt; i++) {
                        j = Ridx[0][Iaidx][i]; 
                        Toccs[i] = Occs[Ka_list][j];
                        }
                     b2brepl(Toccs, Cnt[1], Ij[1], Oij[1], Ridx[1], Sgn[1],
                        AlphaG, Ka_list, cac, Iacnt);
                     for (Ia_ex=0; Ia_ex < Iacnt; Ia_ex++) {
                        if (Cnt[1][Ia_ex]) { found=1; break; }
                        } 
                     }
                  }
               }
            if (found) s2_contrib[sblock][cblock] = 1;
            }

         /* does this c block contribute to sigma3? */
         b2brepl(Occs[sac], Cnt[0], Ij[0], Oij[0], Ridx[0],
            Sgn[0], AlphaG, sac, cac, nas);
         for (Iaidx=0,found=0; Iaidx<S.Ia_size[sblock]; Iaidx++) {
            if (Cnt[0][Iaidx]) found=1;
            }
         if (found) { /* see if beta is ok */
            found=0;
            b2brepl(Occs[sbc], Cnt[0], Ij[0], Oij[0], Ridx[0],
               Sgn[0], BetaG, sbc, cbc, nbs);
            for (Ibidx=0; Ibidx<S.Ib_size[sblock]; Ibidx++) {
               if (Cnt[0][Ibidx]) found=1;
               }
            if (found) s3_contrib[sblock][cblock] = 1;
            }
         } /* end loop over c blocks */
      } /* end loop over sigma blocks */

   #ifdef DEBUG
   printf("\nSigma 1:\n");
   for (i=0; i<S.num_blocks; i++) {
     fprintf(outfile, "Contributions to sigma block %d\n", i);
     for (j=0; j<C.num_blocks; j++) {
       if (s1_contrib[i][j]) fprintf(outfile, "%3d ", j);
     }
     fprintf(outfile, "\n");
   }

   printf("\n\nSigma 2:\n");
   for (i=0; i<S.num_blocks; i++) {
     fprintf(outfile, "Contributions to sigma block %d\n", i);
     for (j=0; j<C.num_blocks; j++) {
       if (s2_contrib[i][j]) fprintf(outfile, "%3d ", j);
     }
     fprintf(outfile, "\n");
   }

   printf("\n\nSigma 3:\n");
   for (i=0; i<S.num_blocks; i++) {
     fprintf(outfile, "Contributions to sigma block %d\n", i);
     for (j=0; j<C.num_blocks; j++) {
       if (s3_contrib[i][j]) fprintf(outfile, "%3d ", j);
     }
     fprintf(outfile, "\n");
   }
   #endif

}



/*
** CALC_SSQ: Computes S^2
**
*/

double CIvect::calc_ssq(double *buffer1, double *buffer2, 
             struct stringwr **alplist, struct stringwr **betlist, int vec_num)
{
   int bra_block, ket_block;  
   int ket_ac, ket_bc, ket_nas, ket_nbs;
   int bra_ac, bra_bc, bra_nas, bra_nbs;
   int ket_birr, bra_birr;
   double tval = 0.0;
   double tval2 = 0.0;
   double S2, Ms; 
   int i; 

   buf_lock(buffer1);
   read(vec_num, 0);

   #ifdef DEBUG
   for (i=0; i<num_blocks; i++) {
      ket_nas = Ia_size[i];
      ket_nbs = Ib_size[i];
      if (ket_nas==0 || ket_nbs==0) continue;
      print_mat(blocks[i], ket_nas, ket_nbs, outfile); 
      }
   #endif
       
   /* loop over ket blocks of c */ 
   for (ket_block=0; ket_block<num_blocks; ket_block++) {
      ket_ac = Ia_code[ket_block];
      ket_bc = Ib_code[ket_block];
      ket_nas = Ia_size[ket_block];
      ket_nbs = Ib_size[ket_block];
      if (ket_nas==0 || ket_nbs==0) continue;
      ket_birr = ket_block / BetaG->subgr_per_irrep;

      for (bra_block=0; bra_block<num_blocks; bra_block++) {
         bra_ac = Ia_code[bra_block];
         bra_bc = Ib_code[bra_block];
         bra_nas = Ia_size[bra_block];
         bra_nbs = Ib_size[bra_block];
         if (bra_nas==0 || bra_nbs==0) continue;
         bra_birr = bra_bc / BetaG->subgr_per_irrep;
         tval2 = ssq(alplist[ket_ac], betlist[ket_bc], blocks[bra_block], 
                   blocks[ket_block], ket_nas, ket_nbs, bra_ac, bra_bc);
         tval += tval2;
         #ifdef DEBUG
         fprintf(outfile,"\nbra_block = %d\n",bra_block);
         fprintf(outfile,"ket_block = %d\n",ket_block);
         fprintf(outfile,"Contribution to <S_S+> = %lf\n",tval2);
         #endif
         } /* end loop over bra_blocks */
    
    } /* end loop over ket_block */

    Ms = 0.5 * (CalcInfo.num_alp_expl - CalcInfo.num_bet_expl);
    #ifdef DEBUG
    fprintf(outfile,"\n\n<S_z> = %lf\n", Ms);
    fprintf(outfile,"<S_z>^2 = %lf\n", Ms*Ms);
    fprintf(outfile,"<S_S+> = %lf\n", tval);
    #endif
    S2 = CalcInfo.num_bet_expl + tval + Ms + Ms*Ms;
    
    fprintf(outfile,"Computed <S^2> vector %d = %20.15f\n\n", vec_num, S2);

  buf_unlock();
  return(S2);
} 

int CIvect::check_zero_block(int blocknum)
{
   if (blocknum < 0 || blocknum > num_blocks) {
      fprintf(stderr, "CIvect::check_zero_block(): Block %d out of range\n",
              blocknum);
   }

   return(zero_blocks[blocknum]);
}

void CIvect::set_zero_block(int blocknum, int value)
{
   if (blocknum < 0 || blocknum > num_blocks) {
      fprintf(stderr, "CIvect::set_zero_block(): Block %d out of range\n",
              blocknum);
   }

   if (value != 0 && value != 1) {
      fprintf(stderr, "CIvect::set_zero_block(): Value %d out of range\n",
              value);
   }

   zero_blocks[blocknum] = value;
}

void CIvect::set_zero_blocks_all(void)
{
   int i;

   for (i=0; i<num_blocks; i++) zero_blocks[i] = 1;

}

void CIvect::copy_zero_blocks(CIvect &src)
{
   int i;

   for (i=0; i<num_blocks; i++) {
      zero_blocks[i] = src.zero_blocks[i];
      // fprintf(outfile, "zero_block[%d] = %d\n", i, zero_blocks[i]);
   }
}

void CIvect::print_zero_blocks(void)
{
   int i;

   for (i=0; i<num_blocks; i++) {
      fprintf(outfile, "zero_block[%d] = %d\n", i, zero_blocks[i]);
   }

}


/*
**
** scale_sigma()
**
** Modifies a sigma vector to account for scaling in
** perturbed Hamiltonian. 
**
** Parameters:
**    Hd         = CIvect for H0
**    S          = CIvect for Sigma vector
**    C          = CIvect for Cvec vector
**    alplist    = alpha string list
**    betlist    = beta  string list
**    buf1       = first buffer for a CIvect
**    buf2       = second buffer for a CIvect
**
** Returns: none
*/
void CIvect::scale_sigma(CIvect &Hd, CIvect &C,
       struct stringwr **alplist, struct stringwr **betlist, int i, 
       double *buf1, double *buf2)
{
   int buf;

   for (buf=0; buf<buf_per_vect; buf++) {
      /* fprintf(outfile," i = %d\n", i);
      fprintf(outfile,"In scale_sigma\n"); */
      Hd.buf_lock(buf1);
      Hd.diag_mat_els_otf(alplist, betlist, CalcInfo.onel_ints,
         CalcInfo.twoel_ints, CalcInfo.e0_fzc, CalcInfo.num_alp_expl,
         CalcInfo.num_bet_expl, CalcInfo.nmo, buf, ORB_ENER);
      C.buf_lock(buf2);
      C.read(i, buf);
      xexy(buf1, buf2, C.buf_size[buf]);
      C.buf_unlock();
      buf_lock(buf2);
      read(i, buf);
      xexmy(buf2, buf1, buf_size[buf]); 
      xpeay(buf1, Parameters.perturbation_parameter, buf2, buf_size[buf]); 
      buf_unlock();
      Hd.buf_unlock();
      buf_lock(buf1);
      write(i, buf);
      buf_unlock();
      }
     
} 

/*
**
** CIvect::dcalc_evangelisti()
**
** Function calculates the denominator part of the Davidson correction
** vector 'd'.
**
** Parameters:
**    rootnum   = number of current root
**    lambda    = current iteration's energy eigenvalue for the current root
**    Hd        = CIvector for the Hamiltonian diagonal
**    C         = CI vector 
**
** Returns: sum of squares of coefficients in d.
*/
double CIvect::dcalc_evangelisti(int rootnum, int num_vecs, double lambda, 
        CIvect &Hd, CIvect &C, double *buf1, double *buf2, int precon, int L, 
        struct stringwr **alplist, struct stringwr **betlist, double **alpha)
{
   int buf, errcod, i;
   double tval, norm = 0.0;

   for (buf=0; buf<buf_per_vect; buf++) {
      Hd.buf_unlock();

      buf_unlock();
      zero_arr(buf1, buf_size[buf]);
      C.buf_lock(buf2);
      for (i=0; i<L; i++) {
         C.read(i, buf); 
         xpeay(buf1, alpha[rootnum][i], buf2, C.buf_size[buf]);         
         }
      C.buf_unlock();
      buf_lock(buf2);
      read(rootnum, buf);

      xexy(buf2, buf1, buf_size[buf]); /* r_I*c_I */
      xeax(buf2, -2.0, buf_size[buf]); /* -2*r_I*c_I */
      xexy(buf1, buf1, buf_size[buf]); /* c_I*c_I */
      xpey(buf1, buf2, buf_size[buf]); /* -2*r_I*c_I + c_I*c_I */
      buf_unlock();
      Hd.buf_lock(buf2);
      if (Parameters.hd_otf == FALSE) Hd.read(0, buf);
      else if (Parameters.hd_otf == TRUE) {
          Hd.diag_mat_els_otf(alplist, betlist, CalcInfo.onel_ints,
             CalcInfo.twoel_ints, CalcInfo.efzc, CalcInfo.num_alp_expl,
             CalcInfo.num_bet_expl, CalcInfo.nmo, buf, Parameters.hd_ave);
        }
      xpey(buf2, buf1, buf_size[buf]); /* Hd -2*r_I*c_I + c_I*c_I */
      buf_lock(buf1);
      read(rootnum, buf);
      tval = calc_d2(buf1, lambda, buf2, buf_size[buf], precon);
      if (buf_offdiag[buf]) tval *= 2.0;
      norm += tval;
      write(rootnum, buf);
      }

  return(norm);
}

/*
** Write the number of the new first buffer to disk.
** The new first buffer is the buffer which is renumbered as "zero" after
** a collapse of the subspace.  The labels on disk are not actually 
** changed so that it is a little easier to deal with two logical CIvectors
** which point to the same physical CIvector (as happens if nodfile).
*/
void CIvect::write_new_first_buf(void)
{
  int unit;

  unit = first_unit;
  psio_write_entry((ULI) unit, "New First Buffer", (char *) &new_first_buf,
    sizeof(int));
}

/*
** Read the number of the new first buffer from disk.
** The new first buffer is the buffer which is renumbered as "zero" after
** a collapse of the subspace.  The labels on disk are not actually 
** changed so that it is a little easier to deal with two logical CIvectors
** which point to the same physical CIvector (as happens if nodfile).
** Return -1 if "New First Buffer" is not stored in the file yet.
*/
int CIvect::read_new_first_buf(void)
{
  int unit;
  int nfb;

  unit = first_unit;
  if (psio_tocscan((ULI) unit, "New First Buffer") == NULL) return(-1);
  psio_read_entry((ULI) unit, "New First Buffer", (char *) &nfb, 
    sizeof(int));
  return(nfb);

}

/*
** Set the number of the new first buffer.
** The new first buffer is the buffer which is renumbered as "zero" after
** a collapse of the subspace.  The labels on disk are not actually 
** changed so that it is a little easier to deal with two logical CIvectors
** which point to the same physical CIvector (as happens if nodfile).
*/
void CIvect::set_new_first_buf(int nfb)
{
  new_first_buf = nfb;
}


/*
** Read the number of valid vectors in this object.  That will be stored
** in the first unit.
*/
int CIvect::read_num_vecs(void)
{
  int unit;
  int nv;

  unit = first_unit;
  if (psio_tocscan((ULI) unit, "Num Vectors") == NULL) return(-1);
  psio_read_entry((ULI) unit, "Num Vectors", (char *) &nv, sizeof(int));
  return(nv);
}


/*
** Write the number of valid vectors in this object.  That will be stored
** in the first unit.
*/
void CIvect::write_num_vecs(int nv)
{
  int unit;

  unit = first_unit;
  psio_write_entry((ULI) unit, "Num Vectors", (char *) &nv, sizeof(int));
  write_toc();
  //civect_psio_debug();
}


/*
** Write the libpsio table of contents to disk in case we crash before
** we're done.  The TOC is written to the end of the file.  If we aren't
** done filling it up, it will be wiped out by the next write but written
** again at the new end of file next time we call this function.
*/
void CIvect::write_toc(void)
{
  int i,unit;

  for (i=0; i<nunits; i++) { 
    psio_tocwrite(units[i]);
  }

}


/*
** Print libpsio debug info
*/
void CIvect::civect_psio_debug(void)
{
  int i, unit;

  /* psio_tocprint not available right now in PSI4; re-enable it if you
     need this functionality.
  */
  /*
  for (i=0; i<nunits; i++)
    psio_tocprint(units[i], outfile);
  */
  fprintf(outfile, "Number of vectors = %d\n", read_num_vecs());
  fprintf(outfile, "New first buffer = %d\n", read_new_first_buf());
  fprintf(outfile, "Internal new first buffer = %d\n", new_first_buf);
}


/*
** perturbation theory correction
** CDS 2/04
*/
void CIvect::pt_correction(struct stringwr **alplist, struct stringwr
      **betlist)
{
  int block, iac, ibc;
  int nas, nbs;

  if (icore==1) { /* whole vector at once */
    for (block=0; block<num_blocks; block++) {
      iac = Ia_code[block];  nas = Ia_size[block];
      ibc = Ib_code[block];  nbs = Ib_size[block];
      // calc_pt_block(alplist[iac], betlist[ibc], blocks[block], nas, nbs);
    }
  }
  else {
    fprintf(outfile, "only icore=1 works for now\n");
  }

}

/*
** CIvect::compute_follow_overlap
**
** Computes the overlap with some user-supplied vector
** Only works for icore==1 (whole vector) for now
**
** Returns: the overlap
*/
double CIvect::compute_follow_overlap(int troot, int ncoef, double *coef, 
  int *Iac, int *Iaridx, int *Ibc, int *Ibridx)
{
  int i, a, b, blk;
  double tval;

  if (icore != 1) {
    fprintf(outfile, "CIvect::compute_follow_overlap: can't use icore != 1\n");
    return(0.0);
  }

  read(troot,0);

  tval = 0.0;

  for (i=0; i<ncoef; i++) {
    blk = decode[Iac[i]][Ibc[i]];
    a = Iaridx[i];
    b = Ibridx[i];
    tval += blocks[blk][a][b] * coef[i];
  }

  tval = fabs(tval); 
  return(tval);

}

}} // namespace psi::detci

