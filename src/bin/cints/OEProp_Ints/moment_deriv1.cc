/*! \file moment_deriv1.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<cmath>
#include<cstdio>
#include<stdlib.h>
#include<cstring>
#include<libipv1/ip_lib.h>
#include<libiwl/iwl.h>
#include<libciomr/libciomr.h>
#include<libint/libint.h>
#include<psifiles.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

namespace psi { namespace CINTS {

static double overlap_int(double a1, int l1, int m1, int n1, double norm1,
			  double a2, int l2, int m2, int n2, double norm2,
			  struct coordinates AB,
			  struct coordinates PA,
			  struct coordinates PB);
static double f_n(int k, int l1, int l2, double A, double B);
static double int_pow(double a, int p);

void moment_deriv1(void)
{
  struct coordinates PA, PB, AB, A, B;
  struct shell_pair *sp;
  register int i, j, k, l, ii, jj, kk, ll;
  int si, sj;
  int ni, nj;
  int l1, m1, n1, l2, m2, n2;
  int am_i, am_j;
  int ai, aj, I, J;
  int ioffset, joffset;
  int atom1, atom2;
  int ax, ay, az, bx, by, bz; /* nuclear coordinate indices */
  int coord, ntri, ij;
  double a1, a2;
  double ab2, oog, gam;
  double inorm, jnorm, over_pf;
  double ***mxderiv, ***myderiv, ***mzderiv; /* Form: mideriv[nuc_coord][ao][ao] */
  double value1, value2, value3, value4;
  double norm1, norm12;
  double *ptr1, *ptr2;
  double *scratch;
  char *label;

  mxderiv = (double ***) malloc(Molecule.num_atoms * 3 * sizeof(double **));
  myderiv = (double ***) malloc(Molecule.num_atoms * 3 * sizeof(double **));
  mzderiv = (double ***) malloc(Molecule.num_atoms * 3 * sizeof(double **));
  for(coord=0; coord < Molecule.num_atoms * 3; coord++) {
    mxderiv[coord] = block_matrix(BasisSet.num_ao, BasisSet.num_ao);
    myderiv[coord] = block_matrix(BasisSet.num_ao, BasisSet.num_ao);
    mzderiv[coord] = block_matrix(BasisSet.num_ao, BasisSet.num_ao);
  }

  for (si=0; si<BasisSet.num_shells; si++){
    am_i = BasisSet.shells[si].am-1;
    ni = ioff[BasisSet.shells[si].am];
    A = Molecule.centers[BasisSet.shells[si].center-1];
    ioffset = BasisSet.shells[si].fao - 1;
    atom1 = BasisSet.shells[si].center-1;

    ax = atom1 * 3;
    ay = ax + 1;
    az = ax + 2;

    for (sj=0; sj<=si; sj++){
      nj = ioff[BasisSet.shells[sj].am];
      am_j = BasisSet.shells[sj].am-1;
      B = Molecule.centers[BasisSet.shells[sj].center-1];
      joffset = BasisSet.shells[sj].fao - 1;
      atom2 = BasisSet.shells[sj].center-1;

      bx = atom2 * 3;
      by = bx + 1;
      bz = bx + 2;

      sp = &(BasisSet.shell_pairs[si][sj]);
      AB.x = sp->AB[0];
      AB.y = sp->AB[1];
      AB.z = sp->AB[2];
      ab2 = AB.x * AB.x;
      ab2 += AB.y * AB.y;
      ab2 += AB.z * AB.z;

      /*--- contract by primitives here ---*/
      for (i = 0; i < BasisSet.shells[si].n_prims; i++) {
	a1 = sp->a1[i];
	inorm = sp->inorm[i];
	for (j = 0; j < BasisSet.shells[sj].n_prims; j++) {
	  a2 = sp->a2[j];
	  gam = sp->gamma[i][j];
	  jnorm = sp->jnorm[j];
	  PA.x = sp->PA[i][j][0];
	  PA.y = sp->PA[i][j][1];
	  PA.z = sp->PA[i][j][2];
	  PB.x = sp->PB[i][j][0];
	  PB.y = sp->PB[i][j][1];
	  PB.z = sp->PB[i][j][2];
	  oog = 1.0/gam;
	  over_pf = exp(-a1*a2*ab2*oog)*sqrt(M_PI*oog)*M_PI*oog*inorm*jnorm;

	  /*--- create all am components of si ---*/
	  ai = 0;
	  for(ii = 0; ii <= am_i; ii++){
	    l1 = am_i - ii;
	    for(jj = 0; jj <= ii; jj++){
	      m1 = ii - jj;
	      n1 = jj;
	      I = ioffset + ai;

	      /*--- create all am components of sj ---*/
	      aj = 0;
	      for(kk = 0; kk <= am_j; kk++){
		l2 = am_j - kk;
		for(ll = 0; ll <= kk; ll++){
		  m2 = kk - ll;
		  n2 = ll;
		  J = joffset + aj;

		  /***************** mu-x derivatives ****************/

		  /**** A derivatives ****/

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_x|b+1_x) */
		  value1 = overlap_int(a1, l1+1, m1, n1, inorm, a2, l2+1, m2, n2, jnorm, AB, PA, PB);
		  /* (a+1_x|b) */
		  value2 = overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  if(l1) {
		    /* (a-1_x|b+1_x) */
		    value3 = overlap_int(a1, l1-1, m1, n1, inorm, a2, l2+1, m2, n2, jnorm, AB, PA, PB);
		    /* (a-1_x|b) */
		    value4 = overlap_int(a1, l1-1, m1, n1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  }
		  mxderiv[ax][I][J] -= 2.0*a1*(value1 + B.x*value2) - l1 * (value3 + B.x*value4);

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_y|b+1_x) */
		  value1 = overlap_int(a1, l1, m1+1, n1, inorm, a2, l2+1, m2, n2, jnorm, AB, PA, PB);
		  /* (a+1_y|b) */
		  value2 = overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  if(m1) {
		    /* (a-1_y|b+1_x) */
		    value3 = overlap_int(a1, l1, m1-1, n1, inorm, a2, l2+1, m2, n2, jnorm, AB, PA, PB);
		    /* (a-1_y|b) */
		    value4 = overlap_int(a1, l1, m1-1, n1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  }
		  mxderiv[ay][I][J] -= 2.0*a1*(value1 + B.x*value2) - m1 * (value3 + B.x*value4);

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_z|b+1_x) */
		  value1 = overlap_int(a1, l1, m1, n1+1, inorm, a2, l2+1, m2, n2, jnorm, AB, PA, PB);
		  /* (a+1_z|b) */
		  value2 = overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  if(n1) {
		    /* (a-1_z|b+1_x) */
		    value3 = overlap_int(a1, l1, m1, n1-1, inorm, a2, l2+1, m2, n2, jnorm, AB, PA, PB);
		    /* (a-1_z|b) */
		    value4 = overlap_int(a1, l1, m1, n1-1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  }
		  mxderiv[az][I][J] -= 2.0*a1*(value1 + B.x*value2) - n1 * (value3 + B.x*value4);


		  /**** B derivatives ****/

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_x|b+1_x) */
		  value1 = overlap_int(a1, l1+1, m1, n1, inorm, a2, l2+1, m2, n2, jnorm, AB, PA, PB);
		  /* (a|b+1_x) */
		  value2 = overlap_int(a1, l1, m1, n1, inorm, a2, l2+1, m2, n2, jnorm, AB, PA, PB);
		  if(l2) {
		    /* (a+1_x|b-1_x) */
		    value3 = overlap_int(a1, l1+1, m1, n1, inorm, a2, l2-1, m2, n2, jnorm, AB, PA, PB);
		    /* (a|b-1_x) */
		    value4 = overlap_int(a1, l1, m1, n1, inorm, a2, l2-1, m2, n2, jnorm, AB, PA, PB);
		  }
		  mxderiv[bx][I][J] -= 2.0*a2*(value1 + A.x*value2) - l2 * (value3 + A.x*value4);

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_x|b+1_y) */
		  value1 = overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2+1, n2, jnorm, AB, PA, PB);
		  /* (a|b+1_y) */
		  value2 = overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2+1, n2, jnorm, AB, PA, PB);
		  if(m2) {
		    /* (a+1_x|b-1_y) */
		    value3 = overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2-1, n2, jnorm, AB, PA, PB);
		    /* (a|b-1_y) */
		    value4 = overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2-1, n2, jnorm, AB, PA, PB);
		  }
		  mxderiv[by][I][J] -= 2.0*a2*(value1 + A.x*value2) - m2 * (value3 + A.x*value4);

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_x|b+1_z) */
		  value1 = overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2, n2+1, jnorm, AB, PA, PB);
		  /* (a|b+1_z) */
		  value2 = overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2, n2+1, jnorm, AB, PA, PB);
		  if(n2) {
		    /* (a+1_x|b-1_z) */
		    value3 = overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2, n2-1, jnorm, AB, PA, PB);
		    /* (a|b-1_z) */
		    value4 = overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2, n2-1, jnorm, AB, PA, PB);
		  }
		  mxderiv[bz][I][J] -= 2.0*a2*(value1 + A.x*value2) - n2 * (value3 + A.x*value4);

		  /***************** mu-y derivatives ****************/

		  /**** A derivatives ****/

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_x|b+1_y) */
		  value1 = overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2+1, n2, jnorm, AB, PA, PB);
		  /* (a+1_x|b) */
		  value2 = overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  if(l1) {
		    /* (a-1_x|b+1_y) */
		    value3 = overlap_int(a1, l1-1, m1, n1, inorm, a2, l2, m2+1, n2, jnorm, AB, PA, PB);
		    /* (a-1_x|b) */
		    value4 = overlap_int(a1, l1-1, m1, n1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  }
		  myderiv[ax][I][J] -= 2.0*a1*(value1 + B.y*value2) - l1 * (value3 + B.y*value4);

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_y|b+1_y) */
		  value1 = overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2+1, n2, jnorm, AB, PA, PB);
		  /* (a+1_y|b) */
		  value2 = overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  if(m1) {
		    /* (a-1_y|b+1_y) */
		    value3 = overlap_int(a1, l1, m1-1, n1, inorm, a2, l2, m2+1, n2, jnorm, AB, PA, PB);
		    /* (a-1_y|b) */
		    value4 = overlap_int(a1, l1, m1-1, n1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  }
		  myderiv[ay][I][J] -= 2.0*a1*(value1 + B.y*value2) - m1 * (value3 + B.y*value4);

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_z|b+1_y) */
		  value1 = overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2+1, n2, jnorm, AB, PA, PB);
		  /* (a+1_z|b) */
		  value2 = overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  if(n1) {
		    /* (a-1_z|b+1_y) */
		    value3 = overlap_int(a1, l1, m1, n1-1, inorm, a2, l2, m2+1, n2, jnorm, AB, PA, PB);
		    /* (a-1_z|b) */
		    value4 = overlap_int(a1, l1, m1, n1-1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  }
		  myderiv[az][I][J] -= 2.0*a1*(value1 + B.y*value2) - n1 * (value3 + B.y*value4);


		  /**** B derivatives ****/

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_y|b+1_x) */
		  value1 = overlap_int(a1, l1, m1+1, n1, inorm, a2, l2+1, m2, n2, jnorm, AB, PA, PB);
		  /* (a|b+1_x) */
		  value2 = overlap_int(a1, l1, m1, n1, inorm, a2, l2+1, m2, n2, jnorm, AB, PA, PB);
		  if(l2) {
		    /* (a+1_y|b-1_x) */
		    value3 = overlap_int(a1, l1, m1+1, n1, inorm, a2, l2-1, m2, n2, jnorm, AB, PA, PB);
		    /* (a|b-1_x) */
		    value4 = overlap_int(a1, l1, m1, n1, inorm, a2, l2-1, m2, n2, jnorm, AB, PA, PB);
		  }
		  myderiv[bx][I][J] -= 2.0*a2*(value1 + A.y*value2) - l2 * (value3 + A.y*value4);

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_y|b+1_y) */
		  value1 = overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2+1, n2, jnorm, AB, PA, PB);
		  /* (a|b+1_y) */
		  value2 = overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2+1, n2, jnorm, AB, PA, PB);
		  if(m2) {
		    /* (a+1_y|b-1_y) */
		    value3 = overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2-1, n2, jnorm, AB, PA, PB);
		    /* (a|b-1_y) */
		    value4 = overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2-1, n2, jnorm, AB, PA, PB);
		  }
		  myderiv[by][I][J] -= 2.0*a2*(value1 + A.y*value2) - m2 * (value3 + A.y*value4);

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_y|b+1_z) */
		  value1 = overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2, n2+1, jnorm, AB, PA, PB);
		  /* (a|b+1_z) */
		  value2 = overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2, n2+1, jnorm, AB, PA, PB);
		  if(n2) {
		    /* (a+1_y|b-1_z) */
		    value3 = overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2, n2-1, jnorm, AB, PA, PB);
		    /* (a|b-1_z) */
		    value4 = overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2, n2-1, jnorm, AB, PA, PB);
		  }
		  myderiv[bz][I][J] -= 2.0*a2*(value1 + A.y*value2) - n2 * (value3 + A.y*value4);

		  /***************** mu-z derivatives ****************/

		  /**** A derivatives ****/

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_x|b+1_z) */
		  value1 = overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2, n2+1, jnorm, AB, PA, PB);
		  /* (a+1_x|b) */
		  value2 = overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  if(l1) {
		    /* (a-1_x|b+1_z) */
		    value3 = overlap_int(a1, l1-1, m1, n1, inorm, a2, l2, m2, n2+1, jnorm, AB, PA, PB);
		    /* (a-1_x|b) */
		    value4 = overlap_int(a1, l1-1, m1, n1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  }
		  mzderiv[ax][I][J] -= 2.0*a1*(value1 + B.z*value2) - l1 * (value3 + B.z*value4);

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_y|b+1_z) */
		  value1 = overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2, n2+1, jnorm, AB, PA, PB);
		  /* (a+1_y|b) */
		  value2 = overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  if(m1) {
		    /* (a-1_y|b+1_z) */
		    value3 = overlap_int(a1, l1, m1-1, n1, inorm, a2, l2, m2, n2+1, jnorm, AB, PA, PB);
		    /* (a-1_y|b) */
		    value4 = overlap_int(a1, l1, m1-1, n1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  }
		  mzderiv[ay][I][J] -= 2.0*a1*(value1 + B.z*value2) - m1 * (value3 + B.z*value4);

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_z|b+1_z) */
		  value1 = overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2, n2+1, jnorm, AB, PA, PB);
		  /* (a+1_z|b) */
		  value2 = overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  if(n1) {
		    /* (a-1_z|b+1_z) */
		    value3 = overlap_int(a1, l1, m1, n1-1, inorm, a2, l2, m2, n2+1, jnorm, AB, PA, PB);
		    /* (a-1_z|b) */
		    value4 = overlap_int(a1, l1, m1, n1-1, inorm, a2, l2, m2, n2, jnorm, AB, PA, PB);
		  }
		  mzderiv[az][I][J] -= 2.0*a1*(value1 + B.z*value2) - n1 * (value3 + B.z*value4);


		  /**** B derivatives ****/

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_z|b+1_x) */
		  value1 = overlap_int(a1, l1, m1, n1+1, inorm, a2, l2+1, m2, n2, jnorm, AB, PA, PB);
		  /* (a|b+1_x) */
		  value2 = overlap_int(a1, l1, m1, n1, inorm, a2, l2+1, m2, n2, jnorm, AB, PA, PB);
		  if(l2) {
		    /* (a+1_z|b-1_x) */
		    value3 = overlap_int(a1, l1, m1, n1+1, inorm, a2, l2-1, m2, n2, jnorm, AB, PA, PB);
		    /* (a|b-1_x) */
		    value4 = overlap_int(a1, l1, m1, n1, inorm, a2, l2-1, m2, n2, jnorm, AB, PA, PB);
		  }
		  mzderiv[bx][I][J] -= 2.0*a2*(value1 + A.z*value2) - l2 * (value3 + A.z*value4);

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_z|b+1_y) */
		  value1 = overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2+1, n2, jnorm, AB, PA, PB);
		  /* (a|b+1_y) */
		  value2 = overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2+1, n2, jnorm, AB, PA, PB);
		  if(m2) {
		    /* (a+1_z|b-1_y) */
		    value3 = overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2-1, n2, jnorm, AB, PA, PB);
		    /* (a|b-1_y) */
		    value4 = overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2-1, n2, jnorm, AB, PA, PB);
		  }
		  mzderiv[by][I][J] -= 2.0*a2*(value1 + A.z*value2) - m2 * (value3 + A.z*value4);

		  value1 = value2 = value3 = value4 = 0.0;

		  /* (a+1_z|b+1_z) */
		  value1 = overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2, n2+1, jnorm, AB, PA, PB);
		  /* (a|b+1_z) */
		  value2 = overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2, n2+1, jnorm, AB, PA, PB);
		  if(n2) {
		    /* (a+1_z|b-1_z) */
		    value3 = overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2, n2-1, jnorm, AB, PA, PB);
		    /* (a|b-1_z) */
		    value4 = overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2, n2-1, jnorm, AB, PA, PB);
		  }
		  mzderiv[bz][I][J] -= 2.0*a2*(value1 + A.z*value2) - n2 * (value3 + A.z*value4);

		  aj++;
		}
	      }
	      ai++;
	    }
	  } /* done with primitives */

	}
      } /* done with contractions */

      /* normalize the contracted integrals */
      ptr1 = GTOs.bf_norm[am_i];
      ptr2 = GTOs.bf_norm[am_j];
      for(i=0; i<ni; i++) {
	norm1 = ptr1[i];
	I = ioffset + i;
	for(j=0; j<nj; j++) {
	  norm12 = norm1*ptr2[j];
	  J = joffset + j;
	  for(coord=0; coord < Molecule.num_atoms*3; coord++) {
	    mxderiv[coord][I][J] *= norm12;
	    myderiv[coord][I][J] *= norm12;
	    mzderiv[coord][I][J] *= norm12;
	  }
	}
      }
    }
  } /* done with this shell pair */

  /* dump the derivative integrals to disk */
  ntri = BasisSet.num_ao * (BasisSet.num_ao + 1)/2;
  scratch = init_array(ntri);
  label = (char *) malloc(PSIO_KEYLEN * sizeof(char));
  for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';
  for(coord=0; coord < Molecule.num_atoms*3; coord++) {

    for(i=0,ij=0; i < BasisSet.num_ao; i++)
      for(j=0; j <= i; j++,ij++)
	scratch[ij] = mxderiv[coord][i][j];
    sprintf(label, "AO-basis MUX Derivs (%d)", coord);
    iwl_wrtone(PSIF_OEI, label, ntri, scratch);
    for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';

    for(i=0,ij=0; i < BasisSet.num_ao; i++)
      for(j=0; j <= i; j++,ij++)
	scratch[ij] = myderiv[coord][i][j];
    sprintf(label, "AO-basis MUY Derivs (%d)", coord);
    iwl_wrtone(PSIF_OEI, label, ntri, scratch);
    for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';

    for(i=0,ij=0; i < BasisSet.num_ao; i++)
      for(j=0; j <= i; j++,ij++)
	scratch[ij] = mzderiv[coord][i][j];
    sprintf(label, "AO-basis MUZ Derivs (%d)", coord);
    iwl_wrtone(PSIF_OEI, label, ntri, scratch);
    for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';

    free_block(mxderiv[coord]);
    free_block(myderiv[coord]);
    free_block(mzderiv[coord]);
  }
  free(mxderiv);
  free(myderiv);
  free(mzderiv);
  free(scratch);
  free(label);

}

/*-----------------------------------------------
  This computes overlap of 2 primitive gaussians
 -----------------------------------------------*/
double overlap_int(double a1, int l1, int m1, int n1, double norm1,
		   double a2, int l2, int m2, int n2, double norm2,
		   struct coordinates AB,
		   struct coordinates PA,
		   struct coordinates PB)
{
  int i, j, k, l;
  int imax, jmax, kmax;
  int print = 0;
  double Ix, Iy, Iz;
  double I;
  double gam;
  double AB2;
  double tval, tval1, tval2 ;
  double norm_fact ;

  
  AB2 = AB.x*AB.x+AB.y*AB.y+AB.z*AB.z;
  gam = a1+a2;
  norm_fact = norm1*norm2; 

  tval1 = 2*gam;
  imax = (l1+l2)/2;
  Ix = 0.0;
  for (i = 0; i<= imax; i++){
    tval = f_n(i*2, l1, l2, PA.x, PB.x);
    tval2 = int_pow(tval1, i);
    Ix += tval*(num_ser[i])/(tval2);
    }
  jmax = (m1+m2)/2;
  Iy = 0.0;
  for (j = 0; j<= jmax; j++){
    tval = f_n(j*2, m1, m2, PA.y, PB.y);
    tval2 = int_pow(tval1, j);
    Iy += tval*num_ser[j]/(tval2);
    }
  kmax = (n1+n2)/2;
  Iz = 0.0;
  for (k = 0; k<= kmax; k++){
    tval = f_n(k*2, n1, n2, PA.z, PB.z);
    tval2 = int_pow(tval1, k);
    Iz += tval*num_ser[k]/(tval2);
    }
 
  I=exp(-1*a1*a2*AB2/gam)*Ix*Iy*Iz*sqrt(M_PI/gam)*(M_PI/gam);

  return I*norm_fact;
}

double f_n(int k, int l1, int l2, double A, double B)
{
  int i, j;
  double sum = 0.0;
  double tval;
  double tmp1, tmp2, tmp3, tmp4 ;

  for(i=0; i<=l1; i++){
    for(j=0; j<=l2; j++){
      tmp1 = tmp2 = tmp3 = tmp4 = 0.0 ;
      if((i+j) == k){
        tmp1 = int_pow(A, (l1-i));
        tmp2 = int_pow(B, (l2-j));
        tmp3 = bc[l1][i];
        tmp4 = bc[l2][j];
        tval = tmp1 * tmp2 * tmp3 * tmp4 ;
        sum += tval ;
        }
      }
    }
  return sum;
}

double int_pow(double a, int p)
{
  register int i;
  double b = 1.0;
  
  for(i=0; i<p; i++) b = b*a;
  return b;
}
};};
