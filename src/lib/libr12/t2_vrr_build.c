/*! \file
    \ingroup R12
    \brief Enter brief description of file here 
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <libint/libint.h>
#include "libr12.h"

extern void punt(char *);
static int hash(int a[2][3], int b[2]);

REALTYPE *t2_vrr_build_xxxx(int am_in[2], prim_data *Data, contr_data *ShellData, REALTYPE *vp, REALTYPE *i0,
			 const REALTYPE *i1, const REALTYPE *i2, const REALTYPE *i3, const REALTYPE *i4)
{
  int i, j, k, l;
  int am[2][3];
  int t1, t2, t3, t4;
  int xyz;
  int la, lc;
  REALTYPE AC[3], U1[3], U0, lzdon, r12t2int;
  static int io[] = {0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153,171,190,210,231,253,276,300,325,351,378,406,435,465};

  la = am_in[0];
  lc = am_in[1];
  lzdon = Data->twozeta_d*Data->oo2n;
  U0 = (Data->twozeta_c - Data->twozeta_d*(Data->twozeta_c*ShellData->CDdotCA + lc + 1))*Data->oo2n;
  AC[0] = ShellData->AC[0];
  AC[1] = ShellData->AC[1];
  AC[2] = ShellData->AC[2];
  U1[0] = ShellData->CD[0]*(Data->twozeta_c*Data->twozeta_d*Data->oo2n);
  U1[1] = ShellData->CD[1]*(Data->twozeta_c*Data->twozeta_d*Data->oo2n);
  U1[2] = ShellData->CD[2]*(Data->twozeta_c*Data->twozeta_d*Data->oo2n);

  t1 = t2 = t3 = t4 = 0;

  for(i = 0; i <= la; i++){
    am[0][0] = la - i;
    for(j = 0; j <= i; j++){
      am[0][1] = i - j;
      am[0][2] = j;

      for(k = 0; k <= lc; k++){
	am[1][0] = lc - k;
	for(l = 0; l <= k; l++){
	  am[1][1] = k - l;
	  am[1][2] = l;
	  
	  /*--- (a0|c0) ---*/
	  r12t2int = U0*(*i0);

	  /*--- add (a+1,0|c0) and (a0|c+1,0) ---*/
	  for(xyz=0;xyz<3;xyz++) {
	    am[1][xyz] = am[1][xyz] + 1;
	    am_in[1] = am_in[1] + 1;
	    t1 = hash(am,am_in);
	    am[1][xyz] = am[1][xyz] - 1;
	    am_in[1] = am_in[1] - 1;
	    
	    am[0][xyz] = am[0][xyz] + 1;
	    am_in[0] = am_in[0] + 1;
	    t2 = hash(am,am_in);
	    am[0][xyz] = am[0][xyz] - 1;
	    am_in[0] = am_in[0] - 1;
	    
	    r12t2int -= U1[xyz]*(i1[t1]-i2[t2]);
	  }
	    
	  /*--- (a-1,0|c+1,0) and (a-1,0|c0) if possible ---*/
	  for(xyz=0;xyz<3;xyz++) {
	    if(am[1][xyz]){
	      am[0][xyz] = am[0][xyz] + 1;
	      am_in[0] = am_in[0] + 1;
	      am[1][xyz] = am[1][xyz] - 1;
	      am_in[1] = am_in[1] - 1;
	      t3 = hash(am,am_in);
	      am[0][xyz] = am[0][xyz] - 1;
	      am_in[0] = am_in[0] - 1;
	      t4 = hash(am,am_in);
	      am[1][xyz] = am[1][xyz] + 1;
	      am_in[1] = am_in[1] + 1;

	      r12t2int += am[1][xyz]*lzdon*(i3[t3] + AC[xyz]*i4[t4]);
	    }
	  }
	  *vp = r12t2int;
	  vp++;
	  i0++;
	}
      }
    }
  }
  
  return vp;
}


int hash(a, b)
  int a[2][3];
  int b[2];
{
  int c[2] = {0,0};
  int i;
  static int io[] = {0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153,171,190,210,231,253,276,300,325,351,378,406,435,465};

  if(b[0]){
    i=b[0]-a[0][0];
    c[0]=i+io[i]-a[0][1];
    }
  if(b[1]){
    i=b[1]-a[1][0];
    c[1]=i+io[i]-a[1][1];
    }

  return c[0]*io[b[1]+1]+c[1];
}
