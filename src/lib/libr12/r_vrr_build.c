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

REALTYPE *r_vrr_build_xxxx(int am_in[2], prim_data *Data, REALTYPE *vp, const REALTYPE *i0, const REALTYPE *i1, REALTYPE *i2,
    const REALTYPE *i3, const REALTYPE *i4, const REALTYPE *i5)
{
  int i, j, k, l;
  int a;
  int flag = 0;
  int am[2][3];
  int t1, t2, t3, t4, t2max;
  int xyz;
  int la, lc;
  REALTYPE PA[3], U1[3], loo4zn, loo2z, loo2p, r12int;
  static int io[] = {0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153,171,190,210,231,253,276,300,325,351,378,406,435,465};

  la = am_in[0];
  lc = am_in[1];
  loo4zn = Data->oo2z*Data->oo2n;
  loo2p = Data->oo2p;
  if (la == 0) { /*--- Decrement on C ---*/
    a = 1;
    t2max = io[lc];
    PA[0] = Data->U[2][0];
    PA[1] = Data->U[2][1];
    PA[2] = Data->U[2][2];
    U1[0] = Data->U[2][0]*Data->oo2z + Data->U[3][0]*Data->oo2n;
    U1[1] = Data->U[2][1]*Data->oo2z + Data->U[3][1]*Data->oo2n;
    U1[2] = Data->U[2][2]*Data->oo2z + Data->U[3][2]*Data->oo2n;
    loo2z = Data->oo2n;
  }
  else {         /*--- Decrement on A ---*/
    a = 0;
    t2max = io[la]*io[lc+1];
    PA[0] = Data->U[0][0];
    PA[1] = Data->U[0][1];
    PA[2] = Data->U[0][2];
    U1[0] = Data->U[0][0]*Data->oo2n + Data->U[1][0]*Data->oo2z;
    U1[1] = Data->U[0][1]*Data->oo2n + Data->U[1][1]*Data->oo2z;
    U1[2] = Data->U[0][2]*Data->oo2n + Data->U[1][2]*Data->oo2z;
    loo2z = Data->oo2z;
  }

  t2 = t3 = t4 = 0;

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
	  
	  if(am[a][2]) xyz = 2;
	  if(am[a][1]) xyz = 1;
	  if(am[a][0]) xyz = 0;
	  
	  if (t2 == t2max) {
	    /*--- reset indices (read Justin Fermann's thesis, pp 36-41 ---*/
	    /*--- (a-1,0|c0) ---*/
	    am[a][xyz] = am[a][xyz] - 1;
	    am_in[a] = am_in[a] - 1;
	    t2 = hash(am,am_in);
	    
	    /*--- (a-2,0|c0) ---*/
	    if (am_in[a]) {
	      am[a][xyz] = am[a][xyz] - 1;
	      am_in[a] = am_in[a] - 1;
	      t3 = hash(am,am_in);
	      am[a][xyz] = am[a][xyz] + 1;
	      am_in[a] = am_in[a] + 1;
	    }
	    
	    /*--- (a-1,0|c-1,0) ---*/
	    if (am_in[a^1]) {
	      am[a^1][0] = am[a^1][0] - 1;
	      am_in[a^1] = am_in[a^1] - 1;
	      t4 = hash(am,am_in);
	      am[a^1][0] = am[a^1][0] + 1;
	      am_in[a^1] = am_in[a^1] + 1;
	    }

	    am[a][xyz] = am[a][xyz] + 1;
	    am_in[a] = am_in[a] + 1;
	  }
	  r12int = PA[xyz]*i0[t2] - U1[xyz]*i3[t2] + loo2p*(*i2);
	  t2++;
	  if(am[a][xyz] > 1){
	    r12int += (am[a][xyz]-1)*(loo2z*i1[t3] - loo4zn*i4[t3]);
	    t3++;
	  }
	  if(am[a^1][xyz] > 0){
	    r12int -= am[a^1][xyz]*loo4zn*i5[t4];
	    t4++;
	  }
	  *vp = r12int;
	  vp++;
	  i2++;
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
