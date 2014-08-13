
#include "ga.h"
#include "cscc.h"  
#include <math.h>   
#include <stdio.h> 

void input(void)
{
//.......................................................................
// Input configuration from an XYZ format file call be.inpt and set up
// initial data structures. Atomic numbers from XYZ file are ignored and
// an atomic number of 4 (Beryllium) is used instead.
//.......................................................................

  int i, j;
  static double r;
  static int ifcnt;
  double ax_test;

  //   initialize variables

  natom = 0;
  for (i = 0; i < maxatom; i++) {
    ax[i] = 0.0;
    ay[i] = 0.0;
    az[i] = 0.0;
  }

  if (GA_Nodeid() == 0) {
    FILE *fp;

    fp = fopen("be.inpt","r");
    fscanf(fp,"%d",&natom);

    /* Read in coordinates */

    for (i = 0; i < natom; i++) {
      fscanf(fp,"%d",&j);
      fscanf(fp,"%lf",&ax[i]);
      fscanf(fp,"%lf",&ay[i]);
      fscanf(fp,"%lf",&az[i]);
    }

    fclose(fp);
  }

  GA_Igop(&natom, 1, "+");
  GA_Dgop(&ax[0], natom, "+");
  GA_Dgop(&ay[0], natom, "+");
  GA_Dgop(&az[0], natom, "+");

  // Set up s-function centers and nuclear charges

  //ifcnt = 1;
  ifcnt = 0; //Fortran Array Index 1===> C Array Index 0===>
  for (i = 0; i < natom; i++) {
    q[i] = 4.0;

    expnt[ifcnt] = 1741.0;
    expnt[ifcnt + 1] = 262.1;
    expnt[ifcnt + 2] = 60.33;
    expnt[ifcnt + 3] = 17.62;
    expnt[ifcnt + 4] = 5.933;
    expnt[ifcnt + 5] = 2.185;
    expnt[ifcnt + 6] = 0.859;
    expnt[ifcnt + 7] = 0.1806;
    expnt[ifcnt + 8] = 0.05835;
    expnt[ifcnt + 9] = 0.3;
    expnt[ifcnt + 10] = 0.3;
    expnt[ifcnt + 11] = 0.3;
    expnt[ifcnt + 12] = 0.3;
    expnt[ifcnt + 13] = 0.3;
    expnt[ifcnt + 14] = 0.3;

    for (j = 0; j < 15; j++) {
      x[ifcnt] = ax[i];
      y[ifcnt] = ay[i];
      z[ifcnt] = az[i];
      if (j == 9) //10-->9 index tranlation 
	x[ifcnt] += 1.6;
      if (j == 10) 
	x[ifcnt] -= 1.6;
      if (j == 11) 
	y[ifcnt] += 1.6;
      if (j == 12) 
	y[ifcnt] -= 1.6;
      if (j == 13) 
	z[ifcnt] += 1.6;
      if (j == 14) 
	z[ifcnt] -= 1.6;
      ifcnt++;
    }
  }

  //  evaluate repulsion energy

  enrep = 0.0;
  for (i = 0; i < natom; i++) {
    for (j = i + 1; j < natom; j++) {
      r = sqrt((ax[i] - ax[j]) * (ax[i] - ax[j]) + (ay[i] - ay[j]) * (ay[i] - ay[j]) +
	       (az[i] - az[j]) * (az[i] - az[j]));
      enrep += q[i] * q[j] / r;
    }
  }

  nocc = natom << 1;
  nbfn = natom * 15;
  nnbfn = nbfn * (nbfn + 1) / 2;
    
  return;
}
