/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <stdlib.h>
#include "dpd.h"

extern "C" {

int dpd_buf4_mat_irrep_row_wrt(dpdbuf4 *Buf, int irrep, int pq)
{
  int method, filerow, all_buf_irrep;
  int rs;  /* dpdfile row and column indices */
  int p, q, r, s;  /* orbital indices */
  int bufpq, bufrs;  /* Input dpdbuf row and column indices */
  int filepq;
  int rowtot, coltot;  /* dpdfile row and column dimensions */
  int b_perm_pq, b_perm_rs, b_peq, b_res;
  int f_perm_pq, f_perm_rs, f_peq, f_res;
  int permute;
  double value; 

  all_buf_irrep = Buf->file.my_irrep;
  /* Row and column dimensions in the DPD file */
  rowtot = Buf->file.params->rowtot[irrep];
  coltot = Buf->file.params->coltot[irrep^all_buf_irrep];

  /* Index packing information */
  b_perm_pq = Buf->params->perm_pq; b_perm_rs = Buf->params->perm_rs;
  f_perm_pq = Buf->file.params->perm_pq; f_perm_rs = Buf->file.params->perm_rs;
  b_peq = Buf->params->peq; b_res = Buf->params->res;
  f_peq = Buf->file.params->peq; f_res = Buf->file.params->res;

  /* Exit if buffer is antisymmetrized */
  if(Buf->anti) {
      fprintf(stderr, "\n\tCannot write antisymmetrized buffer\n");
      fprintf(stderr,   "\tback to original DPD file!\n");
      exit(PSI_RETURN_FAILURE);
    }

  if((b_perm_pq == f_perm_pq) && (b_perm_rs == f_perm_rs) &&
     (b_peq == f_peq) && (b_res == f_res))   method = 12;
  else if((b_perm_pq != f_perm_pq) && (b_perm_rs == f_perm_rs) &&
	  (b_res == f_res)) {
      if(f_perm_pq && !b_perm_pq) method = 21;
      else if(!f_perm_pq && b_perm_pq) method = 23;
      else {
	  fprintf(stderr, "\n\tInvalid second-level method!\n");
	  exit(PSI_RETURN_FAILURE);
	}
    }
  else if((b_perm_pq == f_perm_pq) && (b_perm_rs != f_perm_rs) &&
	  (b_peq == f_peq)) {
      if(f_perm_rs && !b_perm_rs) method = 31;
      else if(!f_perm_rs && b_perm_rs) method = 33;
      else {
	  fprintf(stderr, "\n\tInvalid third-level method!\n");
	  exit(PSI_RETURN_FAILURE);
	}
    }
  else if((b_perm_pq != f_perm_pq) && (b_perm_rs != f_perm_rs)) {
      if(f_perm_pq && !b_perm_pq) {
	  if(f_perm_rs && !b_perm_rs) method = 41;
	  else if(!f_perm_rs && b_perm_rs) method = 42;
	}
      else if(!f_perm_pq && b_perm_pq) {
	  if(f_perm_rs && !b_perm_rs) method = 43;
	  else if(!f_perm_rs && b_perm_rs) method = 45;
	}
      else {
	  fprintf(stderr, "\n\tInvalid fourth-level method!\n");
	  exit(PSI_RETURN_FAILURE);
	}
    }
  else {
      fprintf(stderr, "\n\tInvalid method in dpd_buf_mat_irrep_rd!\n");
      exit(PSI_RETURN_FAILURE);
    }


  switch(method) {
  case 12: /* No change in pq or rs */

      if(Buf->file.incore) {
	  for(rs=0; rs < rowtot; rs++)
	      Buf->file.matrix[irrep][pq][rs] = Buf->matrix[irrep][0][rs];
          dpd_file4_cache_dirty(&(Buf->file));
        }
      else {
	  Buf->file.matrix[irrep] = Buf->matrix[irrep];
	  dpd_file4_mat_irrep_row_wrt(&(Buf->file), irrep, pq);
	}
      
      break;
  case 21: /* Pack pq; no change in rs */
      /* Prepare the output buffer for the output DPD file */
      dpd_file4_mat_irrep_row_init(&(Buf->file), irrep);

      p = Buf->file.params->roworb[irrep][pq][0];
      q = Buf->file.params->roworb[irrep][pq][1];
      filepq = Buf->file.params->rowidx[p][q];

      filerow = Buf->file.incore ? filepq : 0;

      /* Loop over the columns in the dpdbuf */
      for(rs=0; rs < coltot; rs++) {
	  bufrs = rs;

	  value = Buf->matrix[irrep][0][bufrs];

	  /* Assign the value */
	  Buf->file.matrix[irrep][filerow][rs] = value;
	}

      /* Write out the row */
      dpd_file4_mat_irrep_row_wrt(&(Buf->file), irrep, filepq);

      /* Close the input buffer */
      dpd_file4_mat_irrep_row_close(&(Buf->file), irrep);

      break;
  case 23: /* Unpack pq; no change in rs */
      /* I don't know if I'll ever use this, so I'll avoid it for now */
      fprintf(stderr, "\n\tShould you be using method %d?\n", method);
      exit(PSI_RETURN_FAILURE);

      break;
  case 31: /* No change in pq; pack rs */
      /* Prepare the output buffer for the output DPD file */
      dpd_file4_mat_irrep_row_init(&(Buf->file), irrep);

      filerow = Buf->file.incore ? pq : 0;

      /* Loop over the columns in the dpdfile */
      for(rs=0; rs < coltot; rs++) {
	  r = Buf->file.params->colorb[irrep^all_buf_irrep][rs][0];
	  s = Buf->file.params->colorb[irrep^all_buf_irrep][rs][1];
	  bufrs = Buf->params->colidx[r][s];

	  value = Buf->matrix[irrep][0][bufrs];

	  /* Assign the value */
	  Buf->file.matrix[irrep][filerow][rs] = value;
	}

      /* Write out the row */
      dpd_file4_mat_irrep_row_wrt(&(Buf->file), irrep, pq);

      /* Close the input buffer */
      dpd_file4_mat_irrep_row_close(&(Buf->file), irrep);

      break;
  case 33: /* No change in pq; unpack rs */
      /* I'm not sure if I'll ever need this, so I'm removing it for now */
      fprintf(stderr, "\n\tShould you be using method %d?\n", method);
      exit(PSI_RETURN_FAILURE);

      break;
  case 41: /* Pack pq and rs */
      /* Prepare the output buffer for the output DPD file */
      dpd_file4_mat_irrep_row_init(&(Buf->file), irrep);

      p = Buf->file.params->roworb[irrep][pq][0];
      q = Buf->file.params->roworb[irrep][pq][1];
      filepq = Buf->file.params->rowidx[p][q];

      filerow = Buf->file.incore ? filepq : 0;


      /* Loop over the columns in the dpdfile */
      for(rs=0; rs < coltot; rs++) {
	  r = Buf->file.params->colorb[irrep^all_buf_irrep][rs][0];
	  s = Buf->file.params->colorb[irrep^all_buf_irrep][rs][1];
	  bufrs = Buf->params->colidx[r][s];

	  value = Buf->matrix[irrep][0][bufrs];

	  /* Assign the value */
	  Buf->file.matrix[irrep][filerow][rs] = value;
	}

      /* Write out the row */
      dpd_file4_mat_irrep_row_wrt(&(Buf->file), irrep, filepq);

      /* Close the input buffer */
      dpd_file4_mat_irrep_row_close(&(Buf->file), irrep);

      break;
  case 42: /* Pack pq; unpack rs */
      fprintf(stderr, "\n\tHaven't programmed method 42 yet!\n");
      exit(PSI_RETURN_FAILURE);

      break;
  case 43: /* Unpack pq; pack rs */
      fprintf(stderr, "\n\tHaven't programmed method 43 yet!\n");
      exit(PSI_RETURN_FAILURE);

      break;
  case 45: /* Unpack pq and rs */
      /* I'm not sure if I'll ever need this, so I'm removing it for now */
      fprintf(stderr, "\n\tShould you be using method %d?\n", method);
      exit(PSI_RETURN_FAILURE);

      break;
  default:  /* Error trapping */
      fprintf(stderr, "\n\tInvalid switch case in dpd_buf_mat_irrep_rd!\n");
      exit(PSI_RETURN_FAILURE);
      break;
    }
  
  return 0;

}

} /* extern "C" */
