/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include <libqt/qt.h>
#include <libpsio/psio.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

namespace psi {

/* dpd_contract444(): Contracts a pair of four-index quantities to
** give a product four-index quantity.
**
** Arguments:
**   dpdbuf4 *X: A pointer to the leftmost dpd four-index
**               buffer in the product.
**   dpdbuf4 *Y: A pointer to the rightmost dpd four-index
**               buffer in the product.
**   int target_X: Indicates which pair of indices (0 = bra, 1 =
**                 ket) of X is the target pair.
**   int target_Y: Indicates which pair of indices (0 = bra, 1 =
**                 ket) of Y is the target pair.
**   double alpha: A prefactor for the product alpha * X * Y.
**   double beta: A prefactor for the target beta * Z.
*/

int dpd_contract444(dpdbuf4 *X, dpdbuf4 *Y, dpdbuf4 *Z,
		    int target_X, int target_Y, double alpha,
		    double beta)
{
  int n, Hx, Hy, Hz, GX, GY, GZ, nirreps, Xtrans, Ytrans, *numlinks, symlink;
  long int size_Y, size_Z, size_file_X_row;
  int incore, nbuckets;
  long int memoryd, core, rows_per_bucket, rows_left, memtotal;
#if DPD_DEBUG
  int *xrow, *xcol, *yrow, *ycol, *zrow, *zcol;
  double byte_conv;
#endif

  nirreps = X->params->nirreps;
  GX = X->file.my_irrep;
  GY = Y->file.my_irrep;
  GZ = Z->file.my_irrep;

  if(target_X == 0) { Xtrans = 0; numlinks = X->params->coltot; symlink=GX; }
  else if(target_X == 1) { Xtrans = 1; numlinks = X->params->rowtot; symlink=0; }

  if(target_Y == 0) Ytrans = 1;
  else if(target_Y == 1) Ytrans = 0;

#ifdef DPD_DEBUG
  if(Xtrans) { xrow = X->params->coltot; xcol = X->params->rowtot; }
  else { xrow = X->params->rowtot; xcol = X->params->coltot; }

  if(Ytrans) { yrow = Y->params->coltot; ycol = Y->params->rowtot; }
  else { yrow = Y->params->rowtot; ycol = Y->params->coltot; }

  zrow = Z->params->rowtot; zcol = Z->params->coltot;
  
  if((zrow != xrow) || (zcol != ycol) || (xcol != yrow)) {
    fprintf(stderr, "** Alignment error in contract444 **\n");
    dpd_error("dpd_contract444",stderr);
  }

#endif
  

  for(Hx=0; Hx < nirreps; Hx++) {

    if      ((!Xtrans)&&(!Ytrans))  {Hy = Hx^GX;    Hz = Hx;    }
    else if ((!Xtrans)&&( Ytrans))  {Hy = Hx^GX^GY; Hz = Hx;    }
    else if (( Xtrans)&&(!Ytrans))  {Hy = Hx;       Hz = Hx^GX; }
    else /* (( Xtrans)&&( Ytrans))*/{Hy = Hx^GY;    Hz = Hx^GX; }

    size_Y = ((long) Y->params->rowtot[Hy]) * ((long) Y->params->coltot[Hy^GY]);
    size_Z = ((long) Z->params->rowtot[Hz]) * ((long) Z->params->coltot[Hz^GZ]);
    size_file_X_row = ((long) X->file.params->coltot[0]); /* need room for a row of the X->file */
	
    memoryd = dpd_memfree() - (size_Y + size_Z + size_file_X_row);

    if(X->params->rowtot[Hx] && X->params->coltot[Hx^GX]) {

      if(X->params->coltot[Hx^GX])
	rows_per_bucket = memoryd/X->params->coltot[Hx^GX];
      else rows_per_bucket = -1;

      if(rows_per_bucket > X->params->rowtot[Hx])
	rows_per_bucket = X->params->rowtot[Hx];

      if(!rows_per_bucket)
	dpd_error("contract444: Not enough memory for one row", stderr);

      nbuckets = (int) ceil((double) X->params->rowtot[Hx]/
		      (double) rows_per_bucket);

      rows_left = X->params->rowtot[Hx] % rows_per_bucket;
      
      incore = 1;
      if(nbuckets > 1) incore = 0;
    }
    else incore = 1;

#if DPD_DEBUG
    if(!incore) {
      fprintf(stderr, "Contract444: memory information.\n");
      fprintf(stderr, "Contract444: h = %d, row = %d, col = %d, tot = %d\n", 
	      Hx, X->params->rowtot[Hx], X->params->coltot[Hx^GX],
	      X->params->rowtot[Hx] * X->params->coltot[Hx^GX]);

      fprintf(stderr, "Contract444: nbuckets = %d\n", nbuckets);
      fprintf(stderr, "Contract444: rows_per_bucket = %d\n",rows_per_bucket);
      fprintf(stderr, "Contract444: rows_left = %d\n",rows_left);
      memtotal = X->params->rowtot[Hx] * X->params->coltot[Hx^GX];
      byte_conv = ((double) sizeof(double))/1e6;
      fprintf(stderr, "Contract444: out of core algorithm used.\n");
      fprintf(stderr, "Contract444: memtotal = %d.\n", memtotal);
      fprintf(stderr, "Contract444: Need %5.2f MB to run in memory.\n",
	      ((double) memtotal)*byte_conv);
      dpd_file4_cache_print(stderr);
      fflush(stderr);
    }
#endif

    /*
      if(!incore && Xtrans) {
      dpd_file4_cache_print(stderr);
      dpd_error("out-of-core contract444 Xtrans=1 not coded", stderr);
      }
    */

    if(incore) {
      dpd_buf4_mat_irrep_init(X, Hx);
      dpd_buf4_mat_irrep_rd(X, Hx);

      dpd_buf4_mat_irrep_init(Y, Hy);
      dpd_buf4_mat_irrep_rd(Y, Hy);
      dpd_buf4_mat_irrep_init(Z, Hz);
      if(fabs(beta) > 0.0) dpd_buf4_mat_irrep_rd(Z, Hz);

      if(Z->params->rowtot[Hz] &&
	 Z->params->coltot[Hz^GZ] && 
	 numlinks[Hx^symlink]) {
	C_DGEMM(Xtrans?'t':'n', Ytrans?'t':'n', 
		Z->params->rowtot[Hz], Z->params->coltot[Hz^GZ],
		numlinks[Hx^symlink], alpha, 
		&(X->matrix[Hx][0][0]), X->params->coltot[Hx^GX], 
		&(Y->matrix[Hy][0][0]), Y->params->coltot[Hy^GY], beta, 
		&(Z->matrix[Hz][0][0]), Z->params->coltot[Hz^GZ]);
      }

      dpd_buf4_mat_irrep_close(X, Hx);

      dpd_buf4_mat_irrep_wrt(Z, Hz);
      dpd_buf4_mat_irrep_close(Y, Hy);
      dpd_buf4_mat_irrep_close(Z, Hz);
    }
    else {

      /* out-of-core algorithm coded only for NT and TN arrangements, not NN or TT */
      if(!Ytrans && !Xtrans || Ytrans && Xtrans) { 
	fprintf(stderr, "Out-of-core algorithm not yet coded for NN or TT DGEMM.\n");
	dpd_error("contract444", stderr);
      }

      dpd_buf4_mat_irrep_init_block(X, Hx, rows_per_bucket);

      dpd_buf4_mat_irrep_init(Y, Hy);
      dpd_buf4_mat_irrep_rd(Y, Hy);
      dpd_buf4_mat_irrep_init(Z, Hz);
      if(fabs(beta) > 0.0) dpd_buf4_mat_irrep_rd(Z, Hz);

      for(n=0; n < (rows_left ? nbuckets-1 : nbuckets); n++) {

	dpd_buf4_mat_irrep_rd_block(X, Hx, n*rows_per_bucket, rows_per_bucket);

	if(!Xtrans && Ytrans) {
	  if(Z->params->coltot[Hz^GZ] && rows_per_bucket && numlinks[Hx^symlink])
	    C_DGEMM('n', 't', rows_per_bucket, Z->params->coltot[Hz^GZ],
		    numlinks[Hx^symlink], alpha, &(X->matrix[Hx][0][0]), numlinks[Hx^symlink],
		    &(Y->matrix[Hy][0][0]), numlinks[Hx^symlink], beta,
		    &(Z->matrix[Hz][n*rows_per_bucket][0]), Z->params->coltot[Hz^GZ]);
	}
	else if(Xtrans && !Ytrans) {
	  /* CAUTION: We need to accumulate the results of DGEMM for
             each bucket in this case.  So, we set beta="user value"
             on the first bucket, but beta=1 for every bucket
             thereafter. */
	  if(Z->params->coltot[Hz^GZ] && Z->params->rowtot[Hz] && rows_per_bucket)
	    C_DGEMM('t', 'n', Z->params->rowtot[Hz], Z->params->coltot[Hz^GZ],
		    rows_per_bucket, alpha, &(X->matrix[Hx][0][0]), X->params->coltot[Hx^GX],
		    &(Y->matrix[Hy][n*rows_per_bucket][0]), Y->params->coltot[Hy^GY], (n==0 ? beta : 1.0),
		    &(Z->matrix[Hz][0][0]), Z->params->coltot[Hz^GZ]);
	}

      }

      if(rows_left) {

	dpd_buf4_mat_irrep_rd_block(X, Hx, n*rows_per_bucket, rows_left);

	if(!Xtrans && Ytrans) {
	  if(Z->params->coltot[Hz^GZ] && rows_left && numlinks[Hx^symlink])
	    C_DGEMM('n', 't', rows_left, Z->params->coltot[Hz^GZ],
		    numlinks[Hx^symlink], alpha, &(X->matrix[Hx][0][0]), numlinks[Hx^symlink],
		    &(Y->matrix[Hy][0][0]), numlinks[Hx^symlink], beta,
		    &(Z->matrix[Hz][n*rows_per_bucket][0]), Z->params->coltot[Hz^GZ]);
	}
	else if(Xtrans && !Ytrans) {
	  if(Z->params->coltot[Hz^GZ] && Z->params->rowtot[Hz] && rows_left)
	    C_DGEMM('t', 'n', Z->params->rowtot[Hz], Z->params->coltot[Hz^GZ],
		    rows_left, alpha, &(X->matrix[Hx][0][0]), X->params->coltot[Hx^GX],
		    &(Y->matrix[Hy][n*rows_per_bucket][0]), Y->params->coltot[Hy^GY], 1.0,
		    &(Z->matrix[Hz][0][0]), Z->params->coltot[Hz^GZ]);
	}

      }

      dpd_buf4_mat_irrep_close_block(X, Hx, rows_per_bucket);

      dpd_buf4_mat_irrep_close(Y, Hy);
      dpd_buf4_mat_irrep_wrt(Z, Hz);
      dpd_buf4_mat_irrep_close(Z, Hz);
    }
  }

  return 0;
}

}
