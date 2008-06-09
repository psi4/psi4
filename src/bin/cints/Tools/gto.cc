/*! \file gto.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<psiconfig.h>
#include<libciomr/libciomr.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

namespace psi { namespace CINTS {

static const int use_cca_integrals_standard = (PSI_INTEGRALS_STANDARD == 1);

/*-------------------------------
  Explicit function declarations
 -------------------------------*/
static double **init_bf_norm(int);
static double ***init_cart2pureang(int);
static double ****init_cc2pp(int);
static void init_sparse_cc2pp(int);
static double xyz2lm_Coeff(int, int, int, int, int);

void init_gto()
{
  GTOs.bf_norm = init_bf_norm(BasisSet.max_am);
  if (BasisSet.puream && UserOptions.symm_ints) {
     GTOs.cart2pureang = init_cart2pureang(BasisSet.max_am);
#if USE_MM && !SPARSE_C2P
     GTOs.cc2pp = init_cc2pp(BasisSet.max_am);
#endif
#if USE_MM && SPARSE_C2P
     init_sparse_cc2pp(BasisSet.max_am);
#endif
  }

  return;
}


void cleanup_gto()
{
  int l;

  free_matrix(GTOs.bf_norm,BasisSet.max_am);
  
  if (BasisSet.puream && UserOptions.symm_ints) {
     for(l=0;l<BasisSet.max_am;l++)
       free_block(GTOs.cart2pureang[l]);
     free(GTOs.cart2pureang);
   }
}

/*!---------------------------------------------------------
  Computes normalization constants for cartesian Gaussians
 ---------------------------------------------------------*/
double **init_bf_norm(int max_am)
{
  double **bf_norm;
  int am,bf,i,j,l1,m1,n1;

  bf_norm = (double **) malloc(sizeof(double *)*max_am);
  for(am=0; am<max_am; am++) {
    bf = 0;
    bf_norm[am] = init_array(ioff[am+1]);
    for(i=0; i<=am; i++) {
      l1 = am - i;
      for(j=0; j<=i; j++) {
	m1 = i-j;
	n1 = j;
	if (use_cca_integrals_standard)
	  bf_norm[am][bf++] = 1.0;
	else
	  bf_norm[am][bf++] = sqrt(df[2*am]/(df[2*l1]*df[2*m1]*df[2*n1]));
      }
    }
  }

  return bf_norm;
}


double ***init_cart2pureang(int max_am)
{
  int i,j,l,m,bf;
  double ***cart2pureang;

  cart2pureang = (double ***) malloc(sizeof(double **)*(max_am));
  for(l=0;l<max_am;l++)
    cart2pureang[l] = block_matrix(2*l+1,ioff[l+1]);
  for(l=0;l<max_am;l++) {
    for(m=-l;m<=l;m++) {
      bf = 0;
      for(i=0;i<=l;i++)
	for(j=0;j<=i;j++)
	  cart2pureang[l][m+l][bf++] = xyz2lm_Coeff(l,m,l-i,i-j,j);
    }
    if (UserOptions.print_lvl >= PRINT_DEBUG) {
      fprintf(outfile,"  -cart->sph.harm matrix (l = %d):\n",l);
      print_mat(cart2pureang[l],2*l+1,ioff[l+1],outfile);
    }
  }

  return cart2pureang;
}


double ****init_cc2pp(int max_am)
{
  int i1,j1,l1,m1,bf1,i2,j2,l2,m2,bf2;
  int nc2,np2;
  double ****cc2pp;

  cc2pp = (double ****) malloc(sizeof(double ***)*max_am);
  for(l1=0;l1<max_am;l1++)
    cc2pp[l1] = (double ***) malloc(sizeof(double **)*max_am);
  for(l1=0;l1<max_am;l1++)
    for(l2=0;l2<max_am;l2++) {
      nc2 = ioff[l2+1];
      np2 = 2*l2+1;
      cc2pp[l1][l2] = block_matrix((2*l1+1)*np2,ioff[l1+1]*nc2);
      for(m1=-l1;m1<=l1;m1++)
	for(m2=-l2;m2<=l2;m2++) {
	  bf1 = 0;
	  for(i1=0;i1<=l1;i1++)
	    for(j1=0;j1<=i1;j1++,bf1++) {
	      bf2 = 0;
	      for(i2=0;i2<=l2;i2++)
		for(j2=0;j2<=i2;j2++,bf2++)
		  cc2pp[l1][l2][(m1+l1)*np2+m2+l2][bf1*nc2+bf2] = xyz2lm_Coeff(l1,m1,l1-i1,i1-j1,j1)*
								  xyz2lm_Coeff(l2,m2,l2-i2,i2-j2,j2);
	    }
	}
    }

  return cc2pp;
}


void init_sparse_cc2pp(int max_am)
{
  /*--- Used for convenience ---*/
  mat_elem ****cc2pp_sparse;
  int ***cc2pp_rowlength;
  mat_elem ****pp2cc_sparse;
  int ***pp2cc_rowlength;
  mat_elem *scratch;

  int i1,j1,l1,m1,bf1,i2,j2,l2,m2,bf2;
  int nc2,np2;
  int row_index, col_index;
  int col_count, total_count, length;
  int increment;
  double value;

  /*---------------------------------------------------------------------
    cc2pp matrices are computed for all l1 and l2, not only for l1 >= l2
    (1/13/2004) EFV
   ---------------------------------------------------------------------*/
  cc2pp_sparse = (mat_elem ****) malloc(sizeof(mat_elem ***)*max_am);
  cc2pp_rowlength = (int ***) malloc(sizeof(int **)*max_am);
  for(l1=0;l1<max_am;l1++) {
    cc2pp_sparse[l1] = (mat_elem ***) malloc(sizeof(mat_elem **)*max_am);
    cc2pp_rowlength[l1] = (int **) malloc(sizeof(int *)*max_am);
  }
  for(l1=0;l1<max_am;l1++)
    for(l2=0;l2<max_am;l2++) {
      nc2 = ioff[l2+1];
      np2 = 2*l2+1;
      cc2pp_sparse[l1][l2] = (mat_elem **) malloc(sizeof(mat_elem *)*ioff[l1+1]*nc2);
      cc2pp_rowlength[l1][l2] = (int *) malloc(sizeof(int)*ioff[l1+1]*nc2);
      /* allocate some space for now, realloc later */
      increment = (2*l1+1)*np2*ioff[l1+1]*nc2;
      scratch = (mat_elem *) malloc(sizeof(mat_elem)*increment);
      length = increment;
      bf1 = 0;
      total_count = 0;
      for(i1=0;i1<=l1;i1++)
	for(j1=0;j1<=i1;j1++,bf1++) {
	  bf2 = 0;
	  for(i2=0;i2<=l2;i2++)
	    for(j2=0;j2<=i2;j2++,bf2++) {
	      row_index = bf1*nc2+bf2;
	      cc2pp_sparse[l1][l2][row_index] = &(scratch[total_count]);
	      col_count = 0;
	      for(m1=-l1;m1<=l1;m1++)
		for(m2=-l2;m2<=l2;m2++) {
		  cc2pp_sparse[l1][l2][row_index][col_count].column = (m1+l1)*np2+m2+l2;
		  cc2pp_sparse[l1][l2][row_index][col_count].row = row_index;
		  value = xyz2lm_Coeff(l1,m1,l1-i1,i1-j1,j1) * xyz2lm_Coeff(l2,m2,l2-i2,i2-j2,j2);
		  if (fabs(value) > ZERO)
		    cc2pp_sparse[l1][l2][row_index][col_count].value = value;
		  else
		    continue;
		  col_count++;
		  total_count++;
		  if (total_count == length) {
		    scratch = (mat_elem *) realloc(scratch, sizeof(mat_elem)*increment);
		    length += increment;
		  }
		}
	      cc2pp_rowlength[l1][l2][row_index] = col_count;
	    }
	}
    }

  pp2cc_sparse = (mat_elem ****) malloc(sizeof(mat_elem ***)*max_am);
  pp2cc_rowlength = (int ***) malloc(sizeof(int **)*max_am);
  for(l1=0;l1<max_am;l1++) {
    pp2cc_sparse[l1] = (mat_elem ***) malloc(sizeof(mat_elem **)*max_am);
    pp2cc_rowlength[l1] = (int **) malloc(sizeof(int *)*max_am);
  }
  for(l1=0;l1<max_am;l1++)
    for(l2=0;l2<max_am;l2++) {
      nc2 = ioff[l2+1];
      np2 = 2*l2+1;
      pp2cc_sparse[l1][l2] = (mat_elem **) malloc(sizeof(mat_elem *)*(2*l1+1)*np2);
      pp2cc_rowlength[l1][l2] = (int *) malloc(sizeof(int)*(2*l1+1)*np2);
      /* allocate some space for now, realloc later */
      increment = (2*l1+1)*np2*ioff[l1+1]*nc2;
      scratch = (mat_elem *) malloc(sizeof(mat_elem)*increment);
      length = increment;
      total_count = 0;
      for(m1=-l1;m1<=l1;m1++)
	for(m2=-l2;m2<=l2;m2++) {
	  row_index = (m1+l1)*np2+m2+l2;
	  pp2cc_sparse[l1][l2][row_index] = &(scratch[total_count]);
	  col_count = 0;
	  bf1 = 0;
	  for(i1=0;i1<=l1;i1++)
	    for(j1=0;j1<=i1;j1++,bf1++) {
	      bf2 = 0;
	      for(i2=0;i2<=l2;i2++)
		for(j2=0;j2<=i2;j2++,bf2++) {
		  pp2cc_sparse[l1][l2][row_index][col_count].row = row_index;
		  pp2cc_sparse[l1][l2][row_index][col_count].column = bf1*nc2+bf2;
		  value = xyz2lm_Coeff(l1,m1,l1-i1,i1-j1,j1) * xyz2lm_Coeff(l2,m2,l2-i2,i2-j2,j2);
		  if (fabs(value) > ZERO)
		    pp2cc_sparse[l1][l2][row_index][col_count].value = value;
		  else
		    continue;
		  col_count++;
		  total_count++;
		  if (total_count == length) {
		    scratch = (mat_elem *) realloc(scratch, sizeof(mat_elem)*increment);
		    length += increment;
		  }
		}
	      pp2cc_rowlength[l1][l2][row_index] = col_count;
	    }
	}
    }

  GTOs.cc2pp_sparse = cc2pp_sparse;
  GTOs.cc2pp_rowlength = cc2pp_rowlength;
  GTOs.pp2cc_sparse = pp2cc_sparse;
  GTOs.pp2cc_rowlength = pp2cc_rowlength;

  return;
}

  

#define parity(m) ((m)%2 ? -1 : 1)  /*Returns (-1)^m */
/*!---------------------------------------------------------------------------------------------
  Computes transformation coefficients from cartesian to real pure angular momentum functions.
  See IJQC 54, 83 (1995), eqn (15). If m is negative, imaginary part is computed, whereas
  a positive m indicates that the real part of spherical harmonic Ylm is requested.
 ---------------------------------------------------------------------------------------------*/
double xyz2lm_Coeff(int l, int m, int lx, int ly, int lz)
{
  int i,j,k,i_max;
  int k_min, k_max;
  int abs_m;
  int comp;
  int q;
  double pfac, pfac1, sum, sum1;

  abs_m = abs(m);
  if ((lx + ly - abs(m))%2)
    return 0.0;
  else
    j = (lx + ly - abs(m))/2;

  if (j < 0)
    return 0.0;

  /*----------------------------------------------------------------------------------------
    Checking whether the cartesian polynomial contributes to the requested component of Ylm
   ----------------------------------------------------------------------------------------*/
  comp = (m >= 0) ? 1 : -1;
/*  if (comp != ((abs_m-lx)%2 ? -1 : 1))*/
  i = abs_m-lx;
  if (comp != parity(abs(i)))
    return 0.0;

  pfac = sqrt(fac[2*lx]*fac[2*ly]*fac[2*lz]*fac[l-abs_m]/(fac[2*l]*fac[l]*fac[lx]*fac[ly]*fac[lz]*fac[l+abs_m]));
/*  pfac = sqrt(fac[l-abs_m]/(fac[l]*fac[l]*fac[l+abs_m]));*/
  pfac /= (1 << l);
  if (m < 0)
    pfac *= parity((i-1)/2);
  else
    pfac *= parity(i/2);

  i_max = (l-abs_m)/2;
  sum = 0.0;
  for(i=0;i<=i_max;i++) {
    pfac1 = bc[l][i]*bc[i][j];
    if (pfac1 == 0.0)
      continue;
    else
      pfac1 *= (parity(i)*fac[2*(l-i)]/fac[l-abs_m-2*i]);
    sum1 = 0.0;
    k_min = MAX((lx-abs_m)/2,0);
    k_max = MIN(j,lx/2);
    for(k=k_min;k<=k_max;k++)
      sum1 += bc[j][k]*bc[abs_m][lx-2*k]*parity(k);
    sum += pfac1*sum1;
  }

  if (use_cca_integrals_standard)
    sum *= sqrt(df[2*l]/(df[2*lx]*df[2*ly]*df[2*lz]));

  if (m == 0)
    return pfac*sum;
  else
    return M_SQRT2*pfac*sum;
}

};};
