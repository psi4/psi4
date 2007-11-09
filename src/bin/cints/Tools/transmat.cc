/*! \file transmat.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include<stdlib.h>
#include<libciomr/libciomr.h>
#include<cmath>
#include<symmetry.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

namespace psi { namespace CINTS {

/*!------------------------------------------------------------------------------------
  This function builds matrices of transformation coefficients for basis functions of
  each angular momentum type present in the basis
 ------------------------------------------------------------------------------------*/

double ***build_transmat(int *sym_oper, int nirreps, int max_am)
{
  int i,j,l;
  int irr,coeff,ao,ao_max;
  int oper,symop;
  int **xexp_ao, **yexp_ao, **zexp_ao;
  int symX[8], symY[8], symZ[8];
  int x,y,z;
  double ***ao_type_transmat;

    
    /*-------------------------------------------------------------------------------------
      Initialize global arrays of x, y, and z exponents for particular basis function type
     -------------------------------------------------------------------------------------*/
    xexp_ao = init_int_matrix(max_am,ioff[max_am]);
    yexp_ao = init_int_matrix(max_am,ioff[max_am]);
    zexp_ao = init_int_matrix(max_am,ioff[max_am]);
    for(l=0;l<max_am;l++) {
      ao = 0;
      for(i=0; i<=l; i++) {
        x = l - i;
        for(j=0; j<=i; j++) {
          y = i-j;
          z = j;
	  xexp_ao[l][ao] = x;
	  yexp_ao[l][ao] = y;
	  zexp_ao[l][ao] = z;
	  ao++;
	}
      }
    }

    symX[EFLAG]     =  1;  symY[EFLAG]     =  1;  symZ[EFLAG]     =  1;  
    symX[C2XFLAG]   =  1;  symY[C2XFLAG]   = -1;  symZ[C2XFLAG]   = -1;
    symX[C2YFLAG]   = -1;  symY[C2YFLAG]   =  1;  symZ[C2YFLAG]   = -1;  
    symX[C2ZFLAG]   = -1;  symY[C2ZFLAG]   = -1;  symZ[C2ZFLAG]   =  1;  
    symX[IFLAG]     = -1;  symY[IFLAG]     = -1;  symZ[IFLAG]     = -1;  
    symX[SIGXYFLAG] =  1;  symY[SIGXYFLAG] =  1;  symZ[SIGXYFLAG] = -1;  
    symX[SIGXZFLAG] =  1;  symY[SIGXZFLAG] = -1;  symZ[SIGXZFLAG] =  1;  
    symX[SIGYZFLAG] = -1;  symY[SIGYZFLAG] =  1;  symZ[SIGYZFLAG] =  1;

    /*Initialize global arrays*/
    ao_type_transmat = (double ***) malloc(max_am*sizeof(double **));
    for(i=0;i<max_am;i++)
      ao_type_transmat[i] = init_matrix(nirreps,ioff[i+1]);

    for(l=0;l<max_am;l++) {
      for(j=0;j<ioff[l+1];j++)
	ao_type_transmat[l][EFLAG][j] = 1.0;
      switch(l) {
        case 0:
	  for(i=0;i<nirreps;i++)
	    ao_type_transmat[l][i][0] = 1.0;
	  break;

        default:
	  for(i=1;i<nirreps;i++) {
	    oper = sym_oper[i];
	    for(j=0;j<ioff[l+1];j++) {
	      ao_type_transmat[l][i][j] = pow((double)symX[oper],(double)xexp_ao[l][j])*
	                                  pow((double)symY[oper],(double)yexp_ao[l][j])*
	                                  pow((double)symZ[oper],(double)zexp_ao[l][j]);
	    }
	  }
      }
    }

/*    printf("  Transformation matrices :\n\n");
    for(l=0;l<max_am;l++) {
      printf("   l = %d\n",l);
      for(i=0;i<nirreps;i++) {
	printf("Symmetry Operation %d\n",i);
	for(j=0;j<ioff[l+1];j++)
	  printf(" %d  %lf\n",j+1,ao_type_transmat[l][i][j]);
	printf("\n");
      }
      printf("\n");
    }*/
    
    return ao_type_transmat;
}


};};
