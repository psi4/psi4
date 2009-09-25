/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <psiconfig.h>
#include "includes.h"
#include "prototypes.h"
#include "globals.h"

namespace psi { namespace oeprop {

/* This function generates three matrices (xpow_bf, etc.) that
   determine exponents of x, y, and z in every basis function
   arising from the shell of a given angular momentum l in addition to
   computing normalization constants. Normalization constants are
   defined relative to the normalization constant of x^ly^0z^0 which
   is assumed to be unity. Each matrix has lmax rows and every lth row
   has (l+1)*(l+2)/2 entries. */

void init_xyz()
{
  int l, i, j, bf;
  static int use_cca_integrals_standard = (PSI_INTEGRALS_STANDARD == 1);
   
  /* allocate matrices, and generate the content */
  xpow_bf = init_int_matrix(lmax+1,(lmax+1)*(lmax+2)/2);
  ypow_bf = init_int_matrix(lmax+1,(lmax+1)*(lmax+2)/2);
  zpow_bf = init_int_matrix(lmax+1,(lmax+1)*(lmax+2)/2);
  norm_bf = init_matrix(lmax+1,(lmax+1)*(lmax+2)/2);
  for(l=0;l<=lmax;l++) {
    bf = 0;
    for(i = 0; i <= l; i++){
      for(j = 0; j <= i; j++){
	xpow_bf[l][bf] = l - i;
        ypow_bf[l][bf] = i - j;
	zpow_bf[l][bf] = j;
	if (use_cca_integrals_standard)
	  norm_bf[l][bf] = 1.0;
	else
	  norm_bf[l][bf] = sqrt(df[2*l]/(df[2*(l-i)]*df[2*(i-j)]*df[2*j]));
	bf++;
      }
    }
  }

  return;
}


double ***init_box(int a, int b, int c)
{
  int i,j,k;
  double ***box;

  box = (double ***) malloc(sizeof(double **)*a);
  for(i=0;i<a;i++)
    box[i] = (double **) malloc(sizeof(double *)*b);
  for(i=0;i<a;i++)
    for(j=0;j<b;j++) {
      box[i][j] = (double *) malloc(sizeof(double)*c);
      //bzero((char *) box[i][j],sizeof(double)*c);
      memset(box[i][j], '\0', sizeof(double)*c);
    }

  return box;
}


void free_box(double ***box, int a, int b)
{
  int i,j;

  for(i=0;i<a;i++)
    for(j=0;j<b;j++)
      free(box[i][j]);

  for(i=0;i<a;i++)
    free(box[i]);

  free(box);
}

}} // namespace psi::oeprop
