/*! \file 
    \ingroup (INPUT)
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <math.h>
#include <symmetry.h>
#include "input.h"
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

/*-----------------------------------------------------------------------------------------------------------------
  This function builds global transformation matrices for each angular momentum type present in the basis
 -----------------------------------------------------------------------------------------------------------------*/

void build_transmat()
{
  int i,j,l;
  int irr,coeff,ao,ao_max;
  int oper,symop;
  int *xexp, *yexp, *zexp;
  int symX[8], symY[8], symZ[8];
  int max_l;
  int x,y,z;

    
    /*-------------------------------------------------------------------------------------
      Initialize global arrays of x, y, and z exponents for particular basis function type
     -------------------------------------------------------------------------------------*/
    xexp_ao = (int **) malloc(MAXANGMOM*sizeof(int *));
    yexp_ao = (int **) malloc(MAXANGMOM*sizeof(int *));
    zexp_ao = (int **) malloc(MAXANGMOM*sizeof(int *));
    for(l=0;l<MAXANGMOM;l++) {
      xexp_ao[l] = init_int_array(ioff[l+1]);
      yexp_ao[l] = init_int_array(ioff[l+1]);
      zexp_ao[l] = init_int_array(ioff[l+1]);
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
    max_l = MAX(max_angmom+1,MAXANGMOM);
    ao_type_transmat = (double ***) malloc(max_l*sizeof(double **));
    for(i=0;i<max_l;i++)
      ao_type_transmat[i] = init_matrix(nirreps,ioff[i+1]);
    num_cart_so = init_int_matrix(max_l,nirreps);
    ao_type_irr = (int **) malloc(sizeof(int *)*max_l);
    for(l=0;l<max_l;l++)
      ao_type_irr[l] = init_int_array(ioff[l+1]);


    for(l=0;l<max_l;l++) {
      xexp = xexp_ao[l];
      yexp = yexp_ao[l];
      zexp = zexp_ao[l];
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
	    for(j=0;j<ioff[l+1];j++)
	      ao_type_transmat[l][i][j] = pow((double)symX[oper],(double)xexp[j])*
	                                  pow((double)symY[oper],(double)yexp[j])*
	                                  pow((double)symZ[oper],(double)zexp[j]);
	  }
	  break;
      }
    }


    /*--------------------------------------------------------------------------------------------
      Determine the number of cartesian SOs in each symmetry block for each angular momentum type
      and irreps to which every basis function type belongs
     --------------------------------------------------------------------------------------------*/

    num_cart_so[0][0] = 1;
    for(l=1;l<max_l;l++) {
      ao_max = ioff[l+1];
      for(ao=0;ao<ao_max;ao++)
	for(irr=0;irr<nirreps;irr++) {
	  coeff = 0;
	  for(symop=0;symop<nirreps;symop++)
	    coeff += ao_type_transmat[l][symop][ao]*irr_char[irr][symop];
	  if (coeff != 0) {
	    ao_type_irr[l][ao] = irr;
	    num_cart_so[l][irr]++;
	  }
	}
    }

    /*
    fprintf(outfile,"  Transformation matrices :\n\n");
    for(l=0;l<max_l;l++) {
      fprintf(outfile,"   l = %d\n",l);
      for(i=0;i<nirreps;i++) {
	fprintf(outfile,"Symmetry Operation %d\n",i);
	for(j=0;j<ioff[l+1];j++)
	  fprintf(outfile," %d  %lf\n",j+1,ao_type_transmat[l][i][j]);
	fprintf(outfile,"\n");
      }
      fprintf(outfile,"\n");
    }
    */
    
    return;
}

}} // namespace psi::input
