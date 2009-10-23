/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<string>
#include<libipv1/ip_lib.h>
#include<libiwl/iwl.h>
#include<libciomr/libciomr.h>
#include<libint/libint.h>
#include<psifiles.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

namespace psi { namespace cints {

static double overlap_int(double a1, int l1, int m1, int n1, double norm1,
			  double a2, int l2, int m2, int n2, double norm2,
			  struct coordinates AB,
			  struct coordinates PA,
			  struct coordinates PB);
static double f_n(int k, int l1, int l2, double A, double B);
static double int_pow(double a, int p);

void angmom_ints(void)
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
  double **Lx, **Ly, **Lz;
  struct coordinates C;
  double Sxy1, Sxy2, Sxz1, Sxz2;
  double Syx1, Syx2, Syz1, Syz2;
  double Szx1, Szx2, Szy1, Szy2;
  double S0x1, S0x2, S0y1, S0y2, S0z1, S0z2;
  double muxy1, muxy2, muxz1, muxz2;
  double muyx1, muyx2, muyz1, muyz2;
  double muzx1, muzx2, muzy1, muzy2;
  double norm1, norm12, *ptr1, *ptr2;
  double *scratch_x, *scratch_y, *scratch_z;

  Lx = block_matrix(BasisSet.num_ao, BasisSet.num_ao);
  Ly = block_matrix(BasisSet.num_ao, BasisSet.num_ao);
  Lz = block_matrix(BasisSet.num_ao, BasisSet.num_ao);
  scratch_x = init_array(ioff[BasisSet.num_ao]);
  scratch_y = init_array(ioff[BasisSet.num_ao]);
  scratch_z = init_array(ioff[BasisSet.num_ao]);

  C = UserOptions.origin;
  if(UserOptions.print_lvl >= PRINT_OEI)
    fprintf(outfile, "    Reference point for ang. mom. ints. = (%5.3f, %5.3f, %5.3f)\n", C.x, C.y, C.z);

  for (si=0; si<BasisSet.num_shells; si++){
    am_i = BasisSet.shells[si].am-1;
    ni = ioff[BasisSet.shells[si].am];
    A = Molecule.centers[BasisSet.shells[si].center-1];
    ioffset = BasisSet.shells[si].fao - 1;
    atom1 = BasisSet.shells[si].center-1;

    ax = atom1 * 3;
    ay = ax + 1;
    az = ax + 2;

    for (sj=0; sj<BasisSet.num_shells; sj++){
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

		  /*** overlap integrals ***/

		  Sxy1 = Sxy2 = Sxz1 = Sxz2 = 0.0;

		  /* (a+1x|b+1y) */
		  Sxy1 = overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2+1, n2, jnorm, AB, PA, PB);
		  /* (a+1x|b-1y) */
		  if(m2)
		    Sxy2 = overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2-1, n2, jnorm, AB, PA, PB);
		  /* (a+1x|b+1z) */
		  Sxz1 = overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2, n2+1, jnorm, AB, PA, PB);
		  /* (a+1x|b-1z) */
		  if(n2)
		    Sxz2 = overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2, n2-1, jnorm, AB, PA, PB);

		  Syx1 = Syx2 = Syz1 = Syz2 = 0.0;

		  /* (a+1y|b+1x) */
		  Syx1 = overlap_int(a1, l1, m1+1, n1, inorm, a2, l2+1, m2, n2, jnorm, AB, PA, PB);
		  /* (a+1y|b-1x) */
		  if(l2)
		    Syx2 = overlap_int(a1, l1, m1+1, n1, inorm, a2, l2-1, m2, n2, jnorm, AB, PA, PB);
		  /* (a+1y|b+1z) */
		  Syz1 = overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2, n2+1, jnorm, AB, PA, PB);
		  /* (a+1y|b-1z) */
		  if(n2)
		    Syz2 = overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2, n2-1, jnorm, AB, PA, PB);

		  Szx1 = Szx2 = Szy1 = Szy2 = 0.0;

		  /* (a+1z|b+1x) */
		  Szx1 = overlap_int(a1, l1, m1, n1+1, inorm, a2, l2+1, m2, n2, jnorm, AB, PA, PB);
		  /* (a+1z|b-1x) */
		  if(l2)
		    Szx2 = overlap_int(a1, l1, m1, n1+1, inorm, a2, l2-1, m2, n2, jnorm, AB, PA, PB);
		  /* (a+1z|b+1y) */
		  Szy1 = overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2+1, n2, jnorm, AB, PA, PB);
		  /* (a+1z|b-1y) */
		  if(m2)
		    Szy2 = overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2-1, n2, jnorm, AB, PA, PB);

		  S0x1 = S0x2 = S0y1 = S0y2 = S0z1 = S0z2 = 0.0;

		  /* (a|b+1x) */
		  S0x1 = overlap_int(a1, l1, m1, n1, inorm, a2, l2+1, m2, n2, jnorm, AB, PA, PB);
		  /* (a|b-1x) */
		  if(l2)
		    S0x2 = overlap_int(a1, l1, m1, n1, inorm, a2, l2-1, m2, n2, jnorm, AB, PA, PB);
		  /* (a|b+1y) */
		  S0y1 = overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2+1, n2, jnorm, AB, PA, PB);
		  /* (a|b-1y) */
		  if(m2)
		    S0y2 = overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2-1, n2, jnorm, AB, PA, PB);
		  /* (a|b+1z) */
		  S0z1 = overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2, n2+1, jnorm, AB, PA, PB);
		  /* (a|b-1z) */
		  if(n2)
		    S0z2 = overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2, n2-1, jnorm, AB, PA, PB);

		  /*** moment integrals ***/

		  muxy1 = muxy2 = muxz1 = muxz2 = 0.0;

		  /* (a|(x-Cx)|b+1y) */
		  muxy1 = Sxy1 + (A.x - C.x) * S0y1;
		  /* (a|(x-Cx)|b+1z) */
		  muxz1 = Sxz1 + (A.x - C.x) * S0z1;
		  /* (a|(x-Cx)|b-1y) */
		  muxy2 = Sxy2 + (A.x - C.x) * S0y2;
		  /* (a|(x-Cx)|b-1z) */
		  muxz2 = Sxz2 + (A.x - C.x) * S0z2;

		  muyx1 = muyx2 = muyz1 = muyz2 = 0.0;

		  /* (a|(y-Cy)|b+1x) */
		  muyx1 = Syx1 + (A.y - C.y) * S0x1;
		  /* (a|(y-Cy)|b+1z) */
		  muyz1 = Syz1 + (A.y - C.y) * S0z1;
		  /* (a|(y-Cy)|b-1x) */
		  muyx2 = Syx2 + (A.y - C.y) * S0x2;
		  /* (a|(y-Cy)|b-1z) */
		  muyz2 = Syz2 + (A.y - C.y) * S0z2;

		  muzx1 = muzx2 = muzy1 = muzy2 = 0.0;

		  /* (a|(z-Cz)|b+1x) */
		  muzx1 = Szx1 + (A.z - C.z) * S0x1;
		  /* (a|(z-Cz)|b+1y) */
		  muzy1 = Szy1 + (A.z - C.z) * S0y1;
		  /* (a|(z-Cz)|b+1x) */
		  muzx2 = Szx2 + (A.z - C.z) * S0x2;
		  /* (a|(z-Cz)|b+1y) */
		  muzy2 = Szy2 + (A.z - C.z) * S0y2;

		  /********** Lx integrals ***********/

		  /* (a|Lx|b) = 2 a2 * (a|(y-Cy)|b+1z) - B.z * (a|(y-Cy)|b-1z)
		     - 2 a2 * (a|(z-Cz)|b+1y) + B.y * (a|(z-Cz)|b-1y) */

		  Lx[I][J] += 2.0*a2*muyz1 - n2*muyz2 - 2.0*a2*muzy1 + m2*muzy2;
                  
		  /********** Ly integrals ***********/

		  /* (a|Ly|b) = 2 a2 * (a|(z-Cz)|b+1x) - B.x * (a|(z-Cz)|b-1x)
		     - 2 a2 * (a|(x-Cx)|b+1z) + B.z * (a|(x-Cx)|b-1z) */

		  Ly[I][J] += 2.0*a2*muzx1 - l2*muzx2 - 2.0*a2*muxz1 + n2*muxz2;
                  
		  /********** Lz integrals ***********/

		  /* (a|Lz|b) = 2 a2 * (a|(x-Cx)|b+1y) - B.y * (a|(x-Cx)|b-1y)
		     - 2 a2 * (a|(y-Cy)|b+1x) + B.x * (a|(y-Cy)|b-1x) */

		  Lz[I][J] += 2.0*a2*muxy1 - m2*muxy2 - 2.0*a2*muyx1 + l2*muyx2;

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

	  Lx[I][J] *= norm12;
	  Ly[I][J] *= norm12;
	  Lz[I][J] *= norm12;
	}
      }
    }
  } /* done with this shell pair */

  for(i=0,ij=0; i < BasisSet.num_ao; i++)
    for(j=0; j <= i; j++,ij++) {
      scratch_x[ij] = Lx[i][j];
      scratch_y[ij] = Ly[i][j];
      scratch_z[ij] = Lz[i][j];
    }

  if (UserOptions.print_lvl >= PRINT_OEI) {
    fprintf(outfile, "  -AO-Basis LX AngMom Integrals:\n");
    print_mat(Lx, BasisSet.num_ao, BasisSet.num_ao, outfile);
    fprintf(outfile, "  -AO-Basis LY AngMom Integrals:\n");
    print_mat(Ly, BasisSet.num_ao, BasisSet.num_ao, outfile);
    fprintf(outfile, "  -AO-Basis LZ AngMom Integrals:\n");
    print_mat(Lz, BasisSet.num_ao, BasisSet.num_ao, outfile);
    fprintf(outfile,"\n");
  }

  /* dump the integrals to disk here */
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_LX, ioff[BasisSet.num_ao], scratch_x);
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_LY, ioff[BasisSet.num_ao], scratch_y);
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_LZ, ioff[BasisSet.num_ao], scratch_z);

  free(scratch_x);
  free(scratch_y);
  free(scratch_z);
  free_block(Lx);
  free_block(Ly);
  free_block(Lz);
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
}}
