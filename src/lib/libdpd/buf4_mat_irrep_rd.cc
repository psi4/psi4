/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libqt/qt.h>
#include "dpd.h"

namespace psi {
	
/* dpd_buf4_mat_irrep_rd(): Reads an entire irrep from disk into a dpd
** four-index buffer using the "rules" specified when the buffer was
** initialized by dpd_buf4_init().
**
** Arguments:
**   dpdbuf4 *Buf: A pointer to the dpdbuf4 where the data will
**                 be stored.
**   int irrep: The irrep number to be read.
**
** Tested methods: 11,12,22,31,44. (12/3/96)
**
** There are definitely some problems with the read routines here.
** Mainly in the case of unpacking bra indices, there is a danger that
** the input row buffer won't be filled with new values or zeros.  See,
** e.g., method 21.
**
** T. Daniel Crawford
** December 1996
**
** Minor modifications for newest dpd version.
** TDC
** October 1997
**
** Converted to latest version.
** TDC
** September 1999
*/

int dpd_buf4_mat_irrep_rd(dpdbuf4 *Buf, int irrep)
{
  int method, filerow, all_buf_irrep;
  int pq, rs;  /* dpdbuf row and column indices */
  int p, q, r, s;  /* orbital indices */
  int filepq, filers, filesr;  /* Input dpdfile row and column indices */
  int rowtot, coltot;  /* dpdbuf row and column dimensions */
  int b_perm_pq, b_perm_rs, b_peq, b_res;
  int f_perm_pq, f_perm_rs, f_peq, f_res;
  int pq_permute, permute;
  double value;
  long int size;

#ifdef DPD_TIMER
  timer_on("buf_rd");
#endif

  all_buf_irrep = Buf->file.my_irrep;

  rowtot = Buf->params->rowtot[irrep];
  coltot = Buf->params->coltot[irrep^all_buf_irrep];
  size = ((long) rowtot) * ((long) coltot);

  b_perm_pq = Buf->params->perm_pq; b_perm_rs = Buf->params->perm_rs;
  f_perm_pq = Buf->file.params->perm_pq; f_perm_rs = Buf->file.params->perm_rs;
  b_peq = Buf->params->peq; b_res = Buf->params->res;
  f_peq = Buf->file.params->peq; f_res = Buf->file.params->res;

  if((b_perm_pq == f_perm_pq) && (b_perm_rs == f_perm_rs) &&
     (b_peq == f_peq) && (b_res == f_res)) {
      if(Buf->anti) method = 11;
      else method = 12;
      }
  else if((b_perm_pq != f_perm_pq) && (b_perm_rs == f_perm_rs) &&
	  (b_res == f_res)) {
      if(f_perm_pq && !b_perm_pq) {
	  if(Buf->anti) {
	      printf("\n\tUnpack pq and antisymmetrize?\n");
	      exit(PSI_RETURN_FAILURE);
	    }
	  method = 21;
	}
      else if(!f_perm_pq && b_perm_pq) {
	  if(Buf->anti) method = 22;
	  else method = 23;
	}
      else {
	  printf("\n\tInvalid second-level method!\n");
	  exit(PSI_RETURN_FAILURE);
	}
    }
  else if((b_perm_pq == f_perm_pq) && (b_perm_rs != f_perm_rs) &&
	  (b_peq == f_peq)) {
      if(f_perm_rs && !b_perm_rs) {
	  if(Buf->anti) {
	      printf("\n\tUnpack rs and antisymmetrize?\n");
	      exit(PSI_RETURN_FAILURE);
	    }
	  method = 31;
	}
      else if(!f_perm_rs && b_perm_rs) {
	  if(Buf->anti) method = 32;
	  else method = 33;
	}
      else {
	  printf("\n\tInvalid third-level method!\n");
	  exit(PSI_RETURN_FAILURE);
	}
    }
  else if((b_perm_pq != f_perm_pq) && (b_perm_rs != f_perm_rs)) {
      if(f_perm_pq && !b_perm_pq) {
	  if(f_perm_rs && !b_perm_rs) {
	      if(Buf->anti) {
		  printf("\n\tUnpack pq and rs and antisymmetrize?\n");
		  exit(PSI_RETURN_FAILURE);
		}
	      else method = 41;
	    }
	  else if(!f_perm_rs && b_perm_rs) {
	      if(Buf->anti) {
		  printf("\n\tUnpack pq and antisymmetrize?\n");
		  exit(PSI_RETURN_FAILURE);
		}
	      else method = 42;
	    }
	}
      else if(!f_perm_pq && b_perm_pq) {
	  if(f_perm_rs && !b_perm_rs) {
	      if(Buf->anti) {
		  printf("\n\tUnpack rs and antisymmetrize?\n");
		  exit(PSI_RETURN_FAILURE);
		}
	      else method = 43;
	    }
	  else if(!f_perm_rs && b_perm_rs) {
	      if(Buf->anti) method = 44;
	      else method = 45;
	    }
	}
      else {
	  printf("\n\tInvalid fourth-level method!\n");
	  exit(PSI_RETURN_FAILURE);
	}
    }
  else {
      printf("\n\tInvalid method in dpd_buf_mat_irrep_rd!\n");
      exit(PSI_RETURN_FAILURE);
    }


  switch(method) {
  case 11: /* No change in pq or rs; antisymmetrize */

#ifdef DPD_TIMER
      timer_on("buf_rd_11");
#endif

      /* Prepare the input buffer from the input file */
      dpd_file4_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf/dpdfile */
      for(pq=0; pq < rowtot; pq++) {

	  /* Fill the buffer */
	  dpd_file4_mat_irrep_row_rd(&(Buf->file), irrep, pq);

	  filerow = Buf->file.incore ? pq : 0;
	  
	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      r = Buf->params->colorb[irrep^all_buf_irrep][rs][0];
	      s = Buf->params->colorb[irrep^all_buf_irrep][rs][1];

	      /* Column indices in the dpdfile */
	      filers = rs;
	      filesr = Buf->file.params->colidx[s][r];

	      value = Buf->file.matrix[irrep][filerow][filers];

	      value -= Buf->file.matrix[irrep][filerow][filesr];

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = value;
	    }
	}

      /* Close the input buffer */
      dpd_file4_mat_irrep_row_close(&(Buf->file), irrep);

#ifdef DPD_TIMER
      timer_off("buf_rd_11");
#endif

      break;
  
  case 12: /* No change in pq or rs */

#ifdef DPD_TIMER
      timer_on("buf_rd_12");
#endif

      if(Buf->file.incore && size) {
          
          /* We shouldn't actually have to do anything here since the
             pointer to the data should already have been copied in
             buf4_mat_irrep_init(). */
          1;
        }
      else {
	  Buf->file.matrix[irrep] = Buf->matrix[irrep];
	  dpd_file4_mat_irrep_rd(&(Buf->file), irrep);
	}

#ifdef DPD_TIMER
      timer_off("buf_rd_12");
#endif
      
      break;
  case 21: /* Unpack pq; no change in rs */

#ifdef DPD_TIMER
      timer_on("buf_rd_21");
#endif
      
      /* Prepare the input buffer from the input file */
      dpd_file4_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf */
      for(pq=0; pq < rowtot; pq++) {
	  p = Buf->params->roworb[irrep][pq][0];
	  q = Buf->params->roworb[irrep][pq][1];
	  filepq = Buf->file.params->rowidx[p][q];

	  filerow = Buf->file.incore ? filepq : 0;

	  /* Set the permutation operator's value */
	  permute = ((p < q) && (f_perm_pq < 0) ? -1 : 1);

	  /* Fill the buffer */
	  if(filepq >= 0)
	      dpd_file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);
	  else
	      dpd_file4_mat_irrep_row_zero(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      filers = rs;

	      if(filepq >= 0)
		  value = Buf->file.matrix[irrep][filerow][filers];
	      else
		  value = 0;

	      /* Assign the value, keeping track of the sign */
	      Buf->matrix[irrep][pq][rs] = permute*value;
	    }
	}

      /* Close the input buffer */
      dpd_file4_mat_irrep_row_close(&(Buf->file), irrep);

#ifdef DPD_TIMER
      timer_off("buf_rd_21");
#endif

      break;
  case 22: /* Pack pq; no change in rs; antisymmetrize */

#ifdef DPD_TIMER
      timer_on("buf_rd_22");
#endif
      
      /* Prepare the input buffer from the input file */
      dpd_file4_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf */
      for(pq=0; pq < rowtot; pq++) {
	  p = Buf->params->roworb[irrep][pq][0];
	  q = Buf->params->roworb[irrep][pq][1];
	  filepq = Buf->file.params->rowidx[p][q];

	  filerow = Buf->file.incore ? filepq : 0;

	  /* Fill the buffer */
	  dpd_file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      r = Buf->params->colorb[irrep^all_buf_irrep][rs][0];
	      s = Buf->params->colorb[irrep^all_buf_irrep][rs][1];

	      /* Column indices in the dpdfile */
	      filers = rs;
	      filesr = Buf->file.params->colidx[s][r];

	      value = Buf->file.matrix[irrep][filerow][filers];

	      value -= Buf->file.matrix[irrep][filerow][filesr];

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = value;
	    }
	}

      /* Close the input buffer */
      dpd_file4_mat_irrep_row_close(&(Buf->file), irrep);

#ifdef DPD_TIMER
      timer_off("buf_rd_22");
#endif

      break;
  case 23: /* Pack pq; no change in rs */

#ifdef DPD_TIMER
      timer_on("buf_rd_23");
#endif
      
      /* Prepare the input buffer from the input file */
      dpd_file4_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf */
      for(pq=0; pq < rowtot; pq++) {
	  p = Buf->params->roworb[irrep][pq][0];
	  q = Buf->params->roworb[irrep][pq][1];
	  filepq = Buf->file.params->rowidx[p][q];

	  filerow = Buf->file.incore ? filepq : 0;

	  dpd_file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      filers = rs;

	      value = Buf->file.matrix[irrep][filerow][filers];

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = value;
	    }
	}

      /* Close the input buffer */
      dpd_file4_mat_irrep_row_close(&(Buf->file), irrep);

#ifdef DPD_TIMER
      timer_off("buf_rd_23");
#endif

      break;
  case 31: /* No change in pq; unpack rs */

#ifdef DPD_TIMER
      timer_on("buf_rd_31");
#endif
      
      /* Prepare the input buffer from the input file */
      dpd_file4_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf/dpdfile */
      for(pq=0; pq < rowtot; pq++) {
	  filepq = pq;

	  filerow = Buf->file.incore ? filepq : 0;

	  /* Fill the buffer */
	  dpd_file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      r = Buf->params->colorb[irrep^all_buf_irrep][rs][0];
	      s = Buf->params->colorb[irrep^all_buf_irrep][rs][1];
	      filers = Buf->file.params->colidx[r][s];

	      /* rs permutation operator */
	      permute = ((r < s) && (f_perm_rs < 0) ? -1 : 1);

	      /* Is this fast enough? */
	      value = ((filers < 0) ? 0 :
		       Buf->file.matrix[irrep][filerow][filers]);

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = permute*value;
	    }
	}

      /* Close the input buffer */
      dpd_file4_mat_irrep_row_close(&(Buf->file), irrep);

#ifdef DPD_TIMER
      timer_off("buf_rd_31");
#endif

      break;
  case 32: /* No change in pq; pack rs; antisymmetrize */

#ifdef DPD_TIMER
      timer_on("buf_rd_32");
#endif
      
      /* Prepare the input buffer from the input file */
      dpd_file4_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf/dpdfile */
      for(pq=0; pq < rowtot; pq++) {
	  filepq = pq;

	  filerow = Buf->file.incore ? filepq : 0;

	  /* Fill the buffer */
	  dpd_file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      r = Buf->params->colorb[irrep^all_buf_irrep][rs][0];
	      s = Buf->params->colorb[irrep^all_buf_irrep][rs][1];

	      /* Column indices in the dpdfile */
	      filers = Buf->file.params->colidx[r][s];
	      filesr = Buf->file.params->colidx[s][r];

	      value = Buf->file.matrix[irrep][filerow][filers];
	      value -= Buf->file.matrix[irrep][filerow][filesr];

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = value;
	    }
	}

      /* Close the input buffer */
      dpd_file4_mat_irrep_row_close(&(Buf->file), irrep);

#ifdef DPD_TIMER
      timer_off("buf_rd_32");
#endif

      break;
  case 33: /* No change in pq; pack rs */

#ifdef DPD_TIMER
      timer_on("buf_rd_33");
#endif
      
      /* Prepare the input buffer from the input file */
      dpd_file4_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf/dpdfile */
      for(pq=0; pq < rowtot; pq++) {
	  filepq = pq;

	  filerow = Buf->file.incore ? filepq : 0;

	  /* Fill the buffer */
	  dpd_file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      r = Buf->params->colorb[irrep^all_buf_irrep][rs][0];
	      s = Buf->params->colorb[irrep^all_buf_irrep][rs][1];
	      filers = Buf->file.params->colidx[r][s];

	      value = Buf->file.matrix[irrep][filerow][filers];

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = value;
	    }
	}

      /* Close the input buffer */
      dpd_file4_mat_irrep_row_close(&(Buf->file), irrep);

#ifdef DPD_TIMER
      timer_off("buf_rd_33");
#endif

      break;
  case 41: /* Unpack pq and rs */

#ifdef DPD_TIMER
      timer_on("buf_rd_41");
#endif
      
      /* Prepare the input buffer from the input file */
      dpd_file4_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf */
      for(pq=0; pq < rowtot; pq++) {
	  p = Buf->params->roworb[irrep][pq][0];
	  q = Buf->params->roworb[irrep][pq][1];
	  filepq = Buf->file.params->rowidx[p][q];

	  filerow = Buf->file.incore ? filepq : 0;

	  /* Set the value of the pq permutation operator */
	  pq_permute = ((p < q) && (f_perm_pq < 0) ? -1 : 1);

	  /* Fill the buffer */
	  if(filepq >= 0)
	      dpd_file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);
	  else
	      dpd_file4_mat_irrep_row_zero(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      r = Buf->params->colorb[irrep^all_buf_irrep][rs][0];
	      s = Buf->params->colorb[irrep^all_buf_irrep][rs][1];
	      filers = Buf->file.params->colidx[r][s];

	      /* Set the value of the pqrs permutation operator */
	      permute = ((r < s) && (f_perm_rs < 0) ? -1 : 1)*pq_permute;

              value = 0;

	      if(filers >= 0 && filepq >= 0)
		  value = Buf->file.matrix[irrep][filerow][filers];

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = permute*value;
	    }
	}

      /* Close the input buffer */
      dpd_file4_mat_irrep_row_close(&(Buf->file), irrep);

#ifdef DPD_TIMER
      timer_off("buf_rd_41");
#endif
      break;
  case 42: /* Pack pq; unpack rs */
      printf("\n\tHaven't programmed method 42 yet!\n");
      exit(PSI_RETURN_FAILURE);

      break;
  case 43: /* Unpack pq; pack rs */
      printf("\n\tHaven't programmed method 43 yet!\n");
      exit(PSI_RETURN_FAILURE);

      break;
  case 44: /* Pack pq; pack rs; antisymmetrize */

#ifdef DPD_TIMER
      timer_on("buf_rd_44");
#endif
      
      /* Prepare the input buffer from the input file */
      dpd_file4_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf */
      for(pq=0; pq < rowtot; pq++) {
	  p = Buf->params->roworb[irrep][pq][0];
	  q = Buf->params->roworb[irrep][pq][1];
	  filepq = Buf->file.params->rowidx[p][q];

	  filerow = Buf->file.incore ? filepq : 0;

	  /* Fill the buffer */
	  dpd_file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      r = Buf->params->colorb[irrep^all_buf_irrep][rs][0];
	      s = Buf->params->colorb[irrep^all_buf_irrep][rs][1];

	      /* Column indices in the dpdfile */
	      filers = Buf->file.params->colidx[r][s];
	      filesr = Buf->file.params->colidx[s][r];

	      value = Buf->file.matrix[irrep][filerow][filers];
	      value -= Buf->file.matrix[irrep][filerow][filesr];

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = value;
	    }
	}

      /* Close the input buffer */
      dpd_file4_mat_irrep_row_close(&(Buf->file), irrep);

#ifdef DPD_TIMER
      timer_off("buf_rd_44");
#endif

      break;
  case 45: /* Pack pq and rs */

#ifdef DPD_TIMER
      timer_on("buf_rd_45");
#endif
      
      /* Prepare the input buffer from the input file */
      dpd_file4_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf */
      for(pq=0; pq < rowtot; pq++) {
	  p = Buf->params->roworb[irrep][pq][0];
	  q = Buf->params->roworb[irrep][pq][1];
	  filepq = Buf->file.params->rowidx[p][q];

	  filerow = Buf->file.incore ? filepq : 0;

	  dpd_file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      r = Buf->params->colorb[irrep^all_buf_irrep][rs][0];
	      s = Buf->params->colorb[irrep^all_buf_irrep][rs][1];
	      filers = Buf->file.params->colidx[r][s];

	      if(filers < 0) {
		  printf("\n\tNegative colidx in method 44?\n");
		  exit(PSI_RETURN_FAILURE);
		}

	      value = Buf->file.matrix[irrep][filerow][filers];

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = value;
	    }
	}

      /* Close the input buffer */
      dpd_file4_mat_irrep_row_close(&(Buf->file), irrep);

#ifdef DPD_TIMER
      timer_off("buf_rd_45");
#endif

      break;
  default:  /* Error trapping */
      printf("\n\tInvalid switch case in dpd_buf_mat_irrep_rd!\n");
      exit(PSI_RETURN_FAILURE);
      break;
    }

#ifdef DPD_TIMER
  timer_off("buf_rd");
#endif
  
  return 0;

}

}

