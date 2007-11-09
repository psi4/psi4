/*! \file oe_deriv2.cc
    \ingroup (CINTS)
    \brief One-electron integral second derivatives
*/
#include <cstring>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <libciomr/libciomr.h>
#include <libint/libint.h>
#include <libderiv/libderiv.h>
#include <Tools/prints.h>
#include "defines.h"
#define EXTERN
#include "global.h"
#include "oe_osrr.h"
#include "oe_deriv2_osrr.h"
#include"symmetrize.h"
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#endif
#include "small_fns.h"

#define INCLUDE_OVERLAP 1
#define INCLUDE_KINETIC 1
#define INCLUDE_POTENTIAL 1
#define INCLUDE_VAA 1
#define INCLUDE_VCA 1
#define INCLUDE_VCC 1

namespace psi { namespace CINTS {
/*--- These frequently used numbers are to avoid costs of passing parameters ---*/
static double oo2g, oog, gam;

/*-------------------------------
  Explicit function declarations
 -------------------------------*/
static double overlap_int(double a1, int l1, int m1, int n1, double norm1,
			  double a2, int l2, int m2, int n2, double norm2,
			  struct coordinates AB,
			  struct coordinates PA,
			  struct coordinates PB);
static double ke_int(double a1, int l1, int m1, int n1, double norm1,
		     double a2, int l2, int m2, int n2, double norm2,
		     struct coordinates AB,
		     struct coordinates PA,
		     struct coordinates PB);
static double f_n(int k, int l1, int l2, double A, double B);
static double int_pow(double a, int p);

/*!-------------------------------------------------------------
  This function computes derivatives of one-electron integrals
 -------------------------------------------------------------*/
void oe_deriv2()
{
   /* Computing up to second-order derivatives here */
   const int deriv_lvl = 2;
   const int deriv0_lvl = 0;
   const int deriv1_lvl = 1;
   const int deriv2_lvl = 2;

   double **hess_ov, **hess_oe;

   struct coordinates PA, PB, AB, PC;
   struct shell_pair *sp;
   int i, j, k, l, ii, jj, kk, ll;
   int count;
   int si, sj;
   int np_i, np_j;
   int si_fao, sj_fao;
   int I, J;
   int sz;
   int l1, l2, m1, m2, n1, n2;
   int ioffset, joffset ;
   int ij;
   int ijpack;
   int h1;
   int am;
   int dimension ;
   int ni,li,nj,lj,ai,aj;
   int am_i, am_j;
   int ixm, iym, izm, jxm, jym, jzm;
   int ixm1, iym1, izm1, jxm1, jym1, jzm1;
   int indmax, iind,jind;
   int atom, atom1, atom2;
   int coord_ax, coord_ay, coord_az,
     coord_bx, coord_by, coord_bz,
     coord_cx, coord_cy, coord_cz;
   int bf;
   int iimax, jjmax;
   double a1, a2;
   double ab2;
   double tmp;
   double inorm, jnorm, over_pf;
   double norm_pf, normover_pf, dens_pf, wdens_pf, zdens_pf, znormover_pf, hds_norm_pf;
   double twozeta_a, twozeta_b, upuppfac, updownpfac, s_int, t_int, v_int;
   double *ptr1, *ptr2, norm1, norm12;
   double ***AI0;
   double ***AIX, ***AIY, ***AIZ;
   double ***AIXX, ***AIXY, ***AIXZ,
     ***AIYY, ***AIYZ, ***AIZZ;

#ifdef USE_TAYLOR_FM
  /*--- +4*deriv_lvl because of the way we invoke AI_Deriv2_OSrecurs ---*/
  init_Taylor_Fm_Eval(BasisSet.max_am*4-4+4*deriv_lvl,UserOptions.cutoff);
#endif  

  hess_ov = block_matrix(Molecule.num_atoms*3,Molecule.num_atoms*3);
  hess_oe = block_matrix(Molecule.num_atoms*3,Molecule.num_atoms*3);

  indmax = (BasisSet.max_am+deriv_lvl-1)*(BasisSet.max_am+deriv_lvl)*(BasisSet.max_am+deriv_lvl)+1;
  AI0 = init_box(indmax,indmax,2*(BasisSet.max_am+deriv2_lvl)+1);
  AIX = init_box(indmax,indmax,2*(BasisSet.max_am+deriv1_lvl)+1);
  AIY = init_box(indmax,indmax,2*(BasisSet.max_am+deriv1_lvl)+1);
  AIZ = init_box(indmax,indmax,2*(BasisSet.max_am+deriv1_lvl)+1);
  AIXX = init_box(indmax,indmax,2*(BasisSet.max_am+deriv0_lvl)+1);
  AIXY = init_box(indmax,indmax,2*(BasisSet.max_am+deriv0_lvl)+1);
  AIXZ = init_box(indmax,indmax,2*(BasisSet.max_am+deriv0_lvl)+1);
  AIYY = init_box(indmax,indmax,2*(BasisSet.max_am+deriv0_lvl)+1);
  AIYZ = init_box(indmax,indmax,2*(BasisSet.max_am+deriv0_lvl)+1);
  AIZZ = init_box(indmax,indmax,2*(BasisSet.max_am+deriv0_lvl)+1);
  
  for (si=0; si<BasisSet.num_shells; si++){
    ni = ioff[BasisSet.shells[si].am];
    am_i = BasisSet.shells[si].am-1;
    izm = 1;
    iym = am_i+1;
    ixm = iym*iym;
    izm1 = 1;
    iym1 = am_i+deriv2_lvl+1;
    ixm1 = iym1*iym1;
    atom1 = BasisSet.shells[si].center-1;
    si_fao = BasisSet.shells[si].fao-1;

    coord_ax = atom1*3;
    coord_ay = coord_ax+1;
    coord_az = coord_ax+2;

    for (sj=0; sj<=si; sj++){
      nj = ioff[BasisSet.shells[sj].am];
      am_j = BasisSet.shells[sj].am-1;
      jzm = 1;
      jym = am_j+1;
      jxm = jym*jym;
      jzm1 = 1;
      jym1 = am_j+deriv2_lvl+1;
      jxm1 = jym1*jym1;
      atom2 = BasisSet.shells[sj].center-1;
      sj_fao = BasisSet.shells[sj].fao - 1;
	
      coord_bx = atom2*3;
      coord_by = coord_bx+1;
      coord_bz = coord_bx+2;

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
	  twozeta_a = 2.0 * a1;
	  inorm = sp->inorm[i];
	  for (j = 0; j < BasisSet.shells[sj].n_prims; j++) {
	    a2 = sp->a2[j];
	    twozeta_b = 2.0 * a2;
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
		n1 = jj ;
		I = si_fao + ai;

		/*--- create all am components of sj ---*/
		aj = 0;
		for(kk = 0; kk <= am_j; kk++){
		  l2 = am_j - kk;
		  for(ll = 0; ll <= kk; ll++){
		    m2 = kk - ll;
		    n2 = ll;
		    J = sj_fao + aj;

		    if (si == sj && ai < aj)
		      break;

		    norm_pf = (GTOs.bf_norm[am_i][ai] * GTOs.bf_norm[am_j][aj]);
		    dens_pf = Dens[I][J];
		    dens_pf *= norm_pf;
		    wdens_pf = (-1.0)*Lagr[I][J];
		    wdens_pf *= norm_pf;
		    hds_norm_pf = norm_pf;
		    if (I != J) {
		      dens_pf *= 2.0;
		      wdens_pf *= 2.0;
		      norm_pf *= 2.0;
		    }
		    /*--- A factor of 1/2 for correlated Lagrangians ---*/
		    if (strcmp(UserOptions.wfn,"SCF"))
		      wdens_pf *= 0.5;

		    /*----------------------------------
		      Half-derivative overap integrals
		     ----------------------------------*/
		    /*--- d/dAx ---*/
		    tmp = 2.0*a1*overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2, 
					     n2, jnorm, AB, PA, PB);
		    if (l1)
		      tmp -= l1*overlap_int(a1, l1-1, m1, n1, inorm, a2, l2, m2, 
					    n2, jnorm, AB, PA, PB);
		    HDS[coord_ax][I][J] += tmp*hds_norm_pf;

		    /*--- d/dAy ---*/
		    tmp = 2.0*a1*overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2, 
					     n2, jnorm, AB, PA, PB);
		    if (m1)
		      tmp -= m1*overlap_int(a1, l1, m1-1, n1, inorm, a2, l2, m2, 
					    n2, jnorm, AB, PA, PB);
		    HDS[coord_ay][I][J] += tmp*hds_norm_pf;

		    /*--- d/dAz ---*/
		    tmp = 2.0*a1*overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2, 
					     n2, jnorm, AB, PA, PB);
		    if (n1)
		      tmp -= n1*overlap_int(a1, l1, m1, n1-1, inorm, a2, l2, m2, 
					    n2, jnorm, AB, PA, PB);
		    HDS[coord_az][I][J] += tmp*hds_norm_pf;

		    if(I != J) {
                      /*--- d/dAx ---*/
                      tmp = 2.0*a2*overlap_int(a2, l2+1, m2, n2, jnorm, a1, l1, m1,
                                               n1, inorm, AB, PB, PA);
                      if (l2)
                        tmp -= l2*overlap_int(a2, l2-1, m2, n2, jnorm, a1, l1, m1,
                                              n1, inorm, AB, PB, PA);
                      HDS[coord_bx][J][I] += tmp*hds_norm_pf;
  
                      /*--- d/dAy ---*/
                      tmp = 2.0*a2*overlap_int(a2, l2, m2+1, n2, jnorm, a1, l1, m1,
                                               n1, inorm, AB, PB, PA);
                      if (m2)
                        tmp -= m2*overlap_int(a2, l2, m2-1, n2, jnorm, a1, l1, m1,
                                              n1, inorm, AB, PB, PA);
                      HDS[coord_by][J][I] += tmp*hds_norm_pf;
  
                      /*--- d/dAz ---*/
                      tmp = 2.0*a2*overlap_int(a2, l2, m2, n2+1, jnorm, a1, l1, m1,
                                               n1, inorm, AB, PB, PA);
                      if (n2)
                        tmp -= n2*overlap_int(a2, l2, m2, n2-1, jnorm, a1, l1, m1,
                                              n1, inorm, AB, PB, PA);
                      HDS[coord_bz][J][I] += tmp*hds_norm_pf;
		    }

		    /*----------------------------------
		      First derivative overap integrals
		     ----------------------------------*/
		    /*--- d/dAx ---*/
		    tmp = 2.0*a1*overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2, 
					     n2, jnorm, AB, PA, PB);
		    if (l1)
		      tmp -= l1*overlap_int(a1, l1-1, m1, n1, inorm, a2, l2, m2, 
					    n2, jnorm, AB, PA, PB);
		    S[coord_ax][I][J] += tmp*norm_pf;

		    /*--- d/dAy ---*/
		    tmp = 2.0*a1*overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2, 
					     n2, jnorm, AB, PA, PB);
		    if (m1)
		      tmp -= m1*overlap_int(a1, l1, m1-1, n1, inorm, a2, l2, m2, 
					    n2, jnorm, AB, PA, PB);
		    S[coord_ay][I][J] += tmp*norm_pf;

		    /*--- d/dAz ---*/
		    tmp = 2.0*a1*overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2, 
					     n2, jnorm, AB, PA, PB);
		    if (n1)
		      tmp -= n1*overlap_int(a1, l1, m1, n1-1, inorm, a2, l2, m2, 
					    n2, jnorm, AB, PA, PB);
		    S[coord_az][I][J] += tmp*norm_pf;

		    /*--- d/dBx ---*/
		    tmp = 2.0*a2*overlap_int(a1, l1, m1, n1, inorm, a2, l2+1, m2, 
					     n2, jnorm, AB, PA, PB);
		    if (l2)
		      tmp -= l2*overlap_int(a1, l1, m1, n1, inorm, a2, l2-1, m2, 
					    n2, jnorm, AB, PA, PB);
		    S[coord_bx][I][J] += tmp*norm_pf;

		    /*--- d/dBy ---*/
		    tmp = 2.0*a2*overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2+1, 
					     n2, jnorm, AB, PA, PB);
		    if (m2)
		      tmp -= m2*overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2-1, 
					    n2, jnorm, AB, PA, PB);
		    S[coord_by][I][J] += tmp*norm_pf;

		    /*--- d/dBz ---*/
		    tmp = 2.0*a2*overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2, 
					     n2+1, jnorm, AB, PA, PB);
		    if (n2)
		      tmp -= n2*overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2, 
					    n2-1, jnorm, AB, PA, PB);
		    S[coord_bz][I][J] += tmp*norm_pf;


		    /*------------------------------------------
		      First derivative kinetic energy integrals
		     ------------------------------------------*/
		    /*--- d/dAx ---*/
		    tmp = 2.0*a1*ke_int(a1, l1+1, m1, n1, inorm, a2, l2, m2, 
					n2, jnorm, AB, PA, PB);
		    if (l1)
		      tmp -= l1*ke_int(a1, l1-1, m1, n1, inorm, a2, l2, m2, 
				       n2, jnorm, AB, PA, PB);
		    F[coord_ax][I][J] += tmp*norm_pf;

		    /*--- d/dAy ---*/
		    tmp = 2.0*a1*ke_int(a1, l1, m1+1, n1, inorm, a2, l2, m2, 
					n2, jnorm, AB, PA, PB);
		    if (m1)
		      tmp -= m1*ke_int(a1, l1, m1-1, n1, inorm, a2, l2, m2, 
				       n2, jnorm, AB, PA, PB);
		    F[coord_ay][I][J] += tmp*norm_pf;

		    /*--- d/dAz ---*/
		    tmp = 2.0*a1*ke_int(a1, l1, m1, n1+1, inorm, a2, l2, m2, 
					n2, jnorm, AB, PA, PB);
		    if (n1)
		      tmp -= n1*ke_int(a1, l1, m1, n1-1, inorm, a2, l2, m2, 
				       n2, jnorm, AB, PA, PB);
		    F[coord_az][I][J] += tmp*norm_pf;

		    /*--- d/dBx ---*/
		    tmp = 2.0*a2*ke_int(a1, l1, m1, n1, inorm, a2, l2+1, m2, 
					n2, jnorm, AB, PA, PB);
		    if (l2)
		      tmp -= l2*ke_int(a1, l1, m1, n1, inorm, a2, l2-1, m2, 
				       n2, jnorm, AB, PA, PB);
		    F[coord_bx][I][J] += tmp*norm_pf;

		    /*--- d/dBy ---*/
		    tmp = 2.0*a2*ke_int(a1, l1, m1, n1, inorm, a2, l2, m2+1, 
					n2, jnorm, AB, PA, PB);
		    if (m2)
		      tmp -= m2*ke_int(a1, l1, m1, n1, inorm, a2, l2, m2-1, 
				       n2, jnorm, AB, PA, PB);
		    F[coord_by][I][J] += tmp*norm_pf;

		    /*--- d/dBz ---*/
		    tmp = 2.0*a2*ke_int(a1, l1, m1, n1, inorm, a2, l2, m2, 
					n2+1, jnorm, AB, PA, PB);
		    if (n2)
		      tmp -= n2*ke_int(a1, l1, m1, n1, inorm, a2, l2, m2, 
				       n2-1, jnorm, AB, PA, PB);
		    F[coord_bz][I][J] += tmp*norm_pf;


		    /*------------------------------------
		      Second derivative overlap integrals
		     ------------------------------------*/
#if INCLUDE_OVERLAP
		    s_int = overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2, 
					n2, jnorm, AB, PA, PB);

		    upuppfac = twozeta_a*twozeta_a;
		    /*--- d2/dAx2 ---*/
		    tmp = upuppfac*overlap_int(a1, l1+2, m1, n1, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    updownpfac = l1 + 1; if (l1) updownpfac += l1;
		    tmp -= twozeta_a*updownpfac*s_int;
		    if (l1 >= 2)
		      tmp += l1*(l1-1)*overlap_int(a1, l1-2, m1, n1, inorm, a2, l2, m2, 
						   n2, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_ax][coord_ax] += tmp*wdens_pf;

		    /*--- d2/dAy2 ---*/
		    tmp = upuppfac*overlap_int(a1, l1, m1+2, n1, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    updownpfac = m1 + 1; if (m1) updownpfac += m1;
		    tmp -= twozeta_a*updownpfac*s_int;
		    if (m1 >= 2)
		      tmp += m1*(m1-1)*overlap_int(a1, l1, m1-2, n1, inorm, a2, l2, m2, 
						   n2, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_ay][coord_ay] += tmp*wdens_pf;

		    /*--- d2/dAz2 ---*/
		    tmp = upuppfac*overlap_int(a1, l1, m1, n1+2, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    updownpfac = n1 + 1; if (n1) updownpfac += n1;
		    tmp -= twozeta_a*updownpfac*s_int;
		    if (n1 >= 2)
		      tmp += n1*(n1-1)*overlap_int(a1, l1, m1, n1-2, inorm, a2, l2, m2, 
						   n2, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_az][coord_az] += tmp*wdens_pf;

		    /*--- d2/dAxdAy ---*/
		    tmp = upuppfac*overlap_int(a1, l1+1, m1+1, n1, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    if (l1)
		      tmp -= twozeta_a*l1*overlap_int(a1, l1-1, m1+1, n1, inorm, a2, l2, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (m1)
		      tmp -= twozeta_a*m1*overlap_int(a1, l1+1, m1-1, n1, inorm, a2, l2, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (l1 > 0 && m1 > 0)
		      tmp += l1*m1*overlap_int(a1, l1-1, m1-1, n1, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_ax][coord_ay] += tmp*wdens_pf;

		    /*--- d2/dAxdAz ---*/
		    tmp = upuppfac*overlap_int(a1, l1+1, m1, n1+1, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    if (l1)
		      tmp -= twozeta_a*l1*overlap_int(a1, l1-1, m1, n1+1, inorm, a2, l2, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (n1)
		      tmp -= twozeta_a*n1*overlap_int(a1, l1+1, m1, n1-1, inorm, a2, l2, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (l1 > 0 && n1 > 0)
		      tmp += l1*n1*overlap_int(a1, l1-1, m1, n1-1, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_ax][coord_az] += tmp*wdens_pf;

		    /*--- d2/dAydAz ---*/
		    tmp = upuppfac*overlap_int(a1, l1, m1+1, n1+1, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    if (m1)
		      tmp -= twozeta_a*m1*overlap_int(a1, l1, m1-1, n1+1, inorm, a2, l2, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (n1)
		      tmp -= twozeta_a*n1*overlap_int(a1, l1, m1+1, n1-1, inorm, a2, l2, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (m1 > 0 && n1 > 0)
		      tmp += m1*n1*overlap_int(a1, l1, m1-1, n1-1, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_ay][coord_az] += tmp*wdens_pf;


		    upuppfac = twozeta_b*twozeta_b;
		    /*--- d2/dBx2 ---*/
		    tmp = upuppfac*overlap_int(a1, l1, m1, n1, inorm, a2, l2+2, m2, 
					       n2, jnorm, AB, PA, PB);
		    updownpfac = l2 + 1; if (l2) updownpfac += l2;
		    tmp -= twozeta_b*updownpfac*s_int;
		    if (l2 >= 2)
		      tmp += l2*(l2-1)*overlap_int(a1, l1, m1, n1, inorm, a2, l2-2, m2, 
						   n2, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_bx][coord_bx] += tmp*wdens_pf;

		    /*--- d2/dBy2 ---*/
		    tmp = upuppfac*overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2+2, 
					       n2, jnorm, AB, PA, PB);
		    updownpfac = m2 + 1; if (m2) updownpfac += m2;
		    tmp -= twozeta_b*updownpfac*s_int;
		    if (m2 >= 2)
		      tmp += m2*(m2-1)*overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2-2, 
						   n2, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_by][coord_by] += tmp*wdens_pf;

		    /*--- d2/dBz2 ---*/
		    tmp = upuppfac*overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2, 
					       n2+2, jnorm, AB, PA, PB);
		    updownpfac = n2 + 1; if (n2) updownpfac += n2;
		    tmp -= twozeta_b*updownpfac*s_int;
		    if (n2 >= 2)
		      tmp += n2*(n2-1)*overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2, 
						   n2-2, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_bz][coord_bz] += tmp*wdens_pf;

		    /*--- d2/dBxdBy ---*/
		    tmp = upuppfac*overlap_int(a1, l1, m1, n1, inorm, a2, l2+1, m2+1, 
					       n2, jnorm, AB, PA, PB);
		    if (l2)
		      tmp -= twozeta_b*l2*overlap_int(a1, l1, m1, n1, inorm, a2, l2-1, m2+1, 
						      n2, jnorm, AB, PA, PB);
		    if (m2)
		      tmp -= twozeta_b*m2*overlap_int(a1, l1, m1, n1, inorm, a2, l2+1, m2-1, 
						      n2, jnorm, AB, PA, PB);
		    if (l2 > 0 && m2 > 0)
		      tmp += l2*m2*overlap_int(a1, l1, m1, n1, inorm, a2, l2-1, m2-1, 
					       n2, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_bx][coord_by] += tmp*wdens_pf;

		    /*--- d2/dBxdBz ---*/
		    tmp = upuppfac*overlap_int(a1, l1, m1, n1, inorm, a2, l2+1, m2, 
					       n2+1, jnorm, AB, PA, PB);
		    if (l2)
		      tmp -= twozeta_b*l2*overlap_int(a1, l1, m1, n1, inorm, a2, l2-1, m2, 
						      n2+1, jnorm, AB, PA, PB);
		    if (n2)
		      tmp -= twozeta_b*n2*overlap_int(a1, l1, m1, n1, inorm, a2, l2+1, m2, 
						      n2-1, jnorm, AB, PA, PB);
		    if (l2 > 0 && n2 > 0)
		      tmp += l2*n2*overlap_int(a1, l1, m1, n1, inorm, a2, l2-1, m2, 
					       n2-1, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_bx][coord_bz] += tmp*wdens_pf;

		    /*--- d2/dBydBz ---*/
		    tmp = upuppfac*overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2+1, 
					       n2+1, jnorm, AB, PA, PB);
		    if (m2)
		      tmp -= twozeta_b*m2*overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2-1, 
						      n2+1, jnorm, AB, PA, PB);
		    if (n2)
		      tmp -= twozeta_b*n2*overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2+1, 
						      n2-1, jnorm, AB, PA, PB);
		    if (m2 > 0 && n2 > 0)
		      tmp += m2*n2*overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2-1, 
					       n2-1, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_by][coord_bz] += tmp*wdens_pf;


		    upuppfac = twozeta_a*twozeta_b;
		    /*--- d2/dAxdBx ---*/
		    tmp = upuppfac*overlap_int(a1, l1+1, m1, n1, inorm, a2, l2+1, m2, 
					       n2, jnorm, AB, PA, PB);
		    if (l1)
		      tmp -= twozeta_b*l1*overlap_int(a1, l1-1, m1, n1, inorm, a2, l2+1, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (l2)
		      tmp -= twozeta_a*l2*overlap_int(a1, l1+1, m1, n1, inorm, a2, l2-1, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (l1 > 0 && l2 > 0)
		      tmp += l1*l2*overlap_int(a1, l1-1, m1, n1, inorm, a2, l2-1, m2, 
					       n2, jnorm, AB, PA, PB);
		    if (coord_ax == coord_bx)
		      tmp *= 2.0;
		    hess_ov[coord_ax][coord_bx] += tmp*wdens_pf;

		    /*--- d2/dAxdBy ---*/
		    tmp = upuppfac*overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2+1, 
					       n2, jnorm, AB, PA, PB);
		    if (l1)
		      tmp -= twozeta_b*l1*overlap_int(a1, l1-1, m1, n1, inorm, a2, l2, m2+1, 
						      n2, jnorm, AB, PA, PB);
		    if (m2)
		      tmp -= twozeta_a*m2*overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2-1, 
						      n2, jnorm, AB, PA, PB);
		    if (l1 > 0 && m2 > 0)
		      tmp += l1*m2*overlap_int(a1, l1-1, m1, n1, inorm, a2, l2, m2-1, 
					       n2, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_ax][coord_by] += tmp*wdens_pf;

		    /*--- d2/dAxdBz ---*/
		    tmp = upuppfac*overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2, 
					       n2+1, jnorm, AB, PA, PB);
		    if (l1)
		      tmp -= twozeta_b*l1*overlap_int(a1, l1-1, m1, n1, inorm, a2, l2, m2, 
						      n2+1, jnorm, AB, PA, PB);
		    if (n2)
		      tmp -= twozeta_a*n2*overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2, 
						      n2-1, jnorm, AB, PA, PB);
		    if (l1 > 0 && n2 > 0)
		      tmp += l1*n2*overlap_int(a1, l1-1, m1, n1, inorm, a2, l2, m2, 
					       n2-1, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_ax][coord_bz] += tmp*wdens_pf;

		    /*--- d2/dAydBx ---*/
		    tmp = upuppfac*overlap_int(a1, l1, m1+1, n1, inorm, a2, l2+1, m2, 
					       n2, jnorm, AB, PA, PB);
		    if (m1)
		      tmp -= twozeta_b*m1*overlap_int(a1, l1, m1-1, n1, inorm, a2, l2+1, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (l2)
		      tmp -= twozeta_a*l2*overlap_int(a1, l1, m1+1, n1, inorm, a2, l2-1, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (m1 > 0 && l2 > 0)
		      tmp += m1*l2*overlap_int(a1, l1, m1-1, n1, inorm, a2, l2-1, m2, 
					       n2, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_ay][coord_bx] += tmp*wdens_pf;

		    /*--- d2/dAydBy ---*/
		    tmp = upuppfac*overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2+1, 
					       n2, jnorm, AB, PA, PB);
		    if (m1)
		      tmp -= twozeta_b*m1*overlap_int(a1, l1, m1-1, n1, inorm, a2, l2, m2+1, 
						      n2, jnorm, AB, PA, PB);
		    if (m2)
		      tmp -= twozeta_a*m2*overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2-1, 
						      n2, jnorm, AB, PA, PB);
		    if (m1 > 0 && m2 > 0)
		      tmp += m1*m2*overlap_int(a1, l1, m1-1, n1, inorm, a2, l2, m2-1, 
					       n2, jnorm, AB, PA, PB);
		    if (coord_ay == coord_by)
		      tmp *= 2.0;
		    hess_ov[coord_ay][coord_by] += tmp*wdens_pf;

		    /*--- d2/dAydBz ---*/
		    tmp = upuppfac*overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2, 
					       n2+1, jnorm, AB, PA, PB);
		    if (m1)
		      tmp -= twozeta_b*m1*overlap_int(a1, l1, m1-1, n1, inorm, a2, l2, m2, 
						      n2+1, jnorm, AB, PA, PB);
		    if (n2)
		      tmp -= twozeta_a*n2*overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2, 
						      n2-1, jnorm, AB, PA, PB);
		    if (m1 > 0 && n2 > 0)
		      tmp += m1*n2*overlap_int(a1, l1, m1-1, n1, inorm, a2, l2, m2, 
					       n2-1, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_ay][coord_bz] += tmp*wdens_pf;

		    /*--- d2/dAzdBx ---*/
		    tmp = upuppfac*overlap_int(a1, l1, m1, n1+1, inorm, a2, l2+1, m2, 
					       n2, jnorm, AB, PA, PB);
		    if (n1)
		      tmp -= twozeta_b*n1*overlap_int(a1, l1, m1, n1-1, inorm, a2, l2+1, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (l2)
		      tmp -= twozeta_a*l2*overlap_int(a1, l1, m1, n1+1, inorm, a2, l2-1, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (n1 > 0 && l2 > 0)
		      tmp += n1*l2*overlap_int(a1, l1, m1, n1-1, inorm, a2, l2-1, m2, 
					       n2, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_az][coord_bx] += tmp*wdens_pf;

		    /*--- d2/dAzdBy ---*/
		    tmp = upuppfac*overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2+1, 
					       n2, jnorm, AB, PA, PB);
		    if (n1)
		      tmp -= twozeta_b*n1*overlap_int(a1, l1, m1, n1-1, inorm, a2, l2, m2+1, 
						      n2, jnorm, AB, PA, PB);
		    if (m2)
		      tmp -= twozeta_a*m2*overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2-1, 
						      n2, jnorm, AB, PA, PB);
		    if (n1 > 0 && m2 > 0)
		      tmp += n1*m2*overlap_int(a1, l1, m1, n1-1, inorm, a2, l2, m2-1, 
					       n2, jnorm, AB, PA, PB);
		    
		    hess_ov[coord_az][coord_by] += tmp*wdens_pf;

		    /*--- d2/dAzdBz ---*/
		    tmp = upuppfac*overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2, 
					       n2+1, jnorm, AB, PA, PB);
		    if (n1)
		      tmp -= twozeta_b*n1*overlap_int(a1, l1, m1, n1-1, inorm, a2, l2, m2, 
						      n2+1, jnorm, AB, PA, PB);
		    if (n2)
		      tmp -= twozeta_a*n2*overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2, 
						      n2-1, jnorm, AB, PA, PB);
		    if (n1 > 0 && n2 > 0)
		      tmp += n1*n2*overlap_int(a1, l1, m1, n1-1, inorm, a2, l2, m2, 
					       n2-1, jnorm, AB, PA, PB);
		    if (coord_az == coord_bz)
		      tmp *= 2.0;
		    hess_ov[coord_az][coord_bz] += tmp*wdens_pf;
#endif


		    /*-------------------------------------------
		      Second derivative kinetic energy integrals
		     -------------------------------------------*/
#if INCLUDE_KINETIC
		    t_int = ke_int(a1, l1, m1, n1, inorm, a2, l2, m2, 
				   n2, jnorm, AB, PA, PB);

		    upuppfac = twozeta_a*twozeta_a;
		    /*--- d2/dAx2 ---*/
		    tmp = upuppfac*ke_int(a1, l1+2, m1, n1, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    updownpfac = l1 + 1; if (l1) updownpfac += l1;
		    tmp -= twozeta_a*updownpfac*t_int;
		    if (l1 >= 2)
		      tmp += l1*(l1-1)*ke_int(a1, l1-2, m1, n1, inorm, a2, l2, m2, 
						   n2, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_ax][coord_ax] += tmp*dens_pf;

		    /*--- d2/dAy2 ---*/
		    tmp = upuppfac*ke_int(a1, l1, m1+2, n1, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    updownpfac = m1 + 1; if (m1) updownpfac += m1;
		    tmp -= twozeta_a*updownpfac*t_int;
		    if (m1 >= 2)
		      tmp += m1*(m1-1)*ke_int(a1, l1, m1-2, n1, inorm, a2, l2, m2, 
						   n2, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_ay][coord_ay] += tmp*dens_pf;

		    /*--- d2/dAz2 ---*/
		    tmp = upuppfac*ke_int(a1, l1, m1, n1+2, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    updownpfac = n1 + 1; if (n1) updownpfac += n1;
		    tmp -= twozeta_a*updownpfac*t_int;
		    if (n1 >= 2)
		      tmp += n1*(n1-1)*ke_int(a1, l1, m1, n1-2, inorm, a2, l2, m2, 
						   n2, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_az][coord_az] += tmp*dens_pf;

		    /*--- d2/dAxdAy ---*/
		    tmp = upuppfac*ke_int(a1, l1+1, m1+1, n1, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    if (l1)
		      tmp -= twozeta_a*l1*ke_int(a1, l1-1, m1+1, n1, inorm, a2, l2, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (m1)
		      tmp -= twozeta_a*m1*ke_int(a1, l1+1, m1-1, n1, inorm, a2, l2, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (l1 > 0 && m1 > 0)
		      tmp += l1*m1*ke_int(a1, l1-1, m1-1, n1, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_ax][coord_ay] += tmp*dens_pf;

		    /*--- d2/dAxdAz ---*/
		    tmp = upuppfac*ke_int(a1, l1+1, m1, n1+1, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    if (l1)
		      tmp -= twozeta_a*l1*ke_int(a1, l1-1, m1, n1+1, inorm, a2, l2, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (n1)
		      tmp -= twozeta_a*n1*ke_int(a1, l1+1, m1, n1-1, inorm, a2, l2, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (l1 > 0 && n1 > 0)
		      tmp += l1*n1*ke_int(a1, l1-1, m1, n1-1, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_ax][coord_az] += tmp*dens_pf;

		    /*--- d2/dAydAz ---*/
		    tmp = upuppfac*ke_int(a1, l1, m1+1, n1+1, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    if (m1)
		      tmp -= twozeta_a*m1*ke_int(a1, l1, m1-1, n1+1, inorm, a2, l2, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (n1)
		      tmp -= twozeta_a*n1*ke_int(a1, l1, m1+1, n1-1, inorm, a2, l2, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (m1 > 0 && n1 > 0)
		      tmp += m1*n1*ke_int(a1, l1, m1-1, n1-1, inorm, a2, l2, m2, 
					       n2, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_ay][coord_az] += tmp*dens_pf;

		    upuppfac = twozeta_b*twozeta_b;
		    /*--- d2/dBx2 ---*/
		    tmp = upuppfac*ke_int(a1, l1, m1, n1, inorm, a2, l2+2, m2, 
					       n2, jnorm, AB, PA, PB);
		    updownpfac = l2 + 1; if (l2) updownpfac += l2;
		    tmp -= twozeta_b*updownpfac*t_int;
		    if (l2 >= 2)
		      tmp += l2*(l2-1)*ke_int(a1, l1, m1, n1, inorm, a2, l2-2, m2, 
						   n2, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_bx][coord_bx] += tmp*dens_pf;

		    /*--- d2/dBy2 ---*/
		    tmp = upuppfac*ke_int(a1, l1, m1, n1, inorm, a2, l2, m2+2, 
					       n2, jnorm, AB, PA, PB);
		    updownpfac = m2 + 1; if (m2) updownpfac += m2;
		    tmp -= twozeta_b*updownpfac*t_int;
		    if (m2 >= 2)
		      tmp += m2*(m2-1)*ke_int(a1, l1, m1, n1, inorm, a2, l2, m2-2, 
						   n2, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_by][coord_by] += tmp*dens_pf;

		    /*--- d2/dBz2 ---*/
		    tmp = upuppfac*ke_int(a1, l1, m1, n1, inorm, a2, l2, m2, 
					       n2+2, jnorm, AB, PA, PB);
		    updownpfac = n2 + 1; if (n2) updownpfac += n2;
		    tmp -= twozeta_b*updownpfac*t_int;
		    if (n2 >= 2)
		      tmp += n2*(n2-1)*ke_int(a1, l1, m1, n1, inorm, a2, l2, m2, 
						   n2-2, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_bz][coord_bz] += tmp*dens_pf;

		    /*--- d2/dBxdBy ---*/
		    tmp = upuppfac*ke_int(a1, l1, m1, n1, inorm, a2, l2+1, m2+1, 
					       n2, jnorm, AB, PA, PB);
		    if (l2)
		      tmp -= twozeta_b*l2*ke_int(a1, l1, m1, n1, inorm, a2, l2-1, m2+1, 
						      n2, jnorm, AB, PA, PB);
		    if (m2)
		      tmp -= twozeta_b*m2*ke_int(a1, l1, m1, n1, inorm, a2, l2+1, m2-1, 
						      n2, jnorm, AB, PA, PB);
		    if (l2 > 0 && m2 > 0)
		      tmp += l2*m2*ke_int(a1, l1, m1, n1, inorm, a2, l2-1, m2-1, 
					       n2, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_bx][coord_by] += tmp*dens_pf;

		    /*--- d2/dBxdBz ---*/
		    tmp = upuppfac*ke_int(a1, l1, m1, n1, inorm, a2, l2+1, m2, 
					       n2+1, jnorm, AB, PA, PB);
		    if (l2)
		      tmp -= twozeta_b*l2*ke_int(a1, l1, m1, n1, inorm, a2, l2-1, m2, 
						      n2+1, jnorm, AB, PA, PB);
		    if (n2)
		      tmp -= twozeta_b*n2*ke_int(a1, l1, m1, n1, inorm, a2, l2+1, m2, 
						      n2-1, jnorm, AB, PA, PB);
		    if (l2 > 0 && n2 > 0)
		      tmp += l2*n2*ke_int(a1, l1, m1, n1, inorm, a2, l2-1, m2, 
					       n2-1, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_bx][coord_bz] += tmp*dens_pf;

		    /*--- d2/dBydBz ---*/
		    tmp = upuppfac*ke_int(a1, l1, m1, n1, inorm, a2, l2, m2+1, 
					       n2+1, jnorm, AB, PA, PB);
		    if (m2)
		      tmp -= twozeta_b*m2*ke_int(a1, l1, m1, n1, inorm, a2, l2, m2-1, 
						      n2+1, jnorm, AB, PA, PB);
		    if (n2)
		      tmp -= twozeta_b*n2*ke_int(a1, l1, m1, n1, inorm, a2, l2, m2+1, 
						      n2-1, jnorm, AB, PA, PB);
		    if (m2 > 0 && n2 > 0)
		      tmp += m2*n2*ke_int(a1, l1, m1, n1, inorm, a2, l2, m2-1, 
					       n2-1, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_by][coord_bz] += tmp*dens_pf;


		    upuppfac = twozeta_a*twozeta_b;
		    /*--- d2/dAxdBx ---*/
		    tmp = upuppfac*ke_int(a1, l1+1, m1, n1, inorm, a2, l2+1, m2, 
					       n2, jnorm, AB, PA, PB);
		    if (l1)
		      tmp -= twozeta_b*l1*ke_int(a1, l1-1, m1, n1, inorm, a2, l2+1, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (l2)
		      tmp -= twozeta_a*l2*ke_int(a1, l1+1, m1, n1, inorm, a2, l2-1, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (l1 > 0 && l2 > 0)
		      tmp += l1*l2*ke_int(a1, l1-1, m1, n1, inorm, a2, l2-1, m2, 
					       n2, jnorm, AB, PA, PB);
		    if (coord_ax == coord_bx)
		      tmp *= 2.0;
		    hess_oe[coord_ax][coord_bx] += tmp*dens_pf;

		    /*--- d2/dAxdBy ---*/
		    tmp = upuppfac*ke_int(a1, l1+1, m1, n1, inorm, a2, l2, m2+1, 
					       n2, jnorm, AB, PA, PB);
		    if (l1)
		      tmp -= twozeta_b*l1*ke_int(a1, l1-1, m1, n1, inorm, a2, l2, m2+1, 
						      n2, jnorm, AB, PA, PB);
		    if (m2)
		      tmp -= twozeta_a*m2*ke_int(a1, l1+1, m1, n1, inorm, a2, l2, m2-1, 
						      n2, jnorm, AB, PA, PB);
		    if (l1 > 0 && m2 > 0)
		      tmp += l1*m2*ke_int(a1, l1-1, m1, n1, inorm, a2, l2, m2-1, 
					       n2, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_ax][coord_by] += tmp*dens_pf;

		    /*--- d2/dAxdBz ---*/
		    tmp = upuppfac*ke_int(a1, l1+1, m1, n1, inorm, a2, l2, m2, 
					       n2+1, jnorm, AB, PA, PB);
		    if (l1)
		      tmp -= twozeta_b*l1*ke_int(a1, l1-1, m1, n1, inorm, a2, l2, m2, 
						      n2+1, jnorm, AB, PA, PB);
		    if (n2)
		      tmp -= twozeta_a*n2*ke_int(a1, l1+1, m1, n1, inorm, a2, l2, m2, 
						      n2-1, jnorm, AB, PA, PB);
		    if (l1 > 0 && n2 > 0)
		      tmp += l1*n2*ke_int(a1, l1-1, m1, n1, inorm, a2, l2, m2, 
					       n2-1, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_ax][coord_bz] += tmp*dens_pf;

		    /*--- d2/dAydBx ---*/
		    tmp = upuppfac*ke_int(a1, l1, m1+1, n1, inorm, a2, l2+1, m2, 
					       n2, jnorm, AB, PA, PB);
		    if (m1)
		      tmp -= twozeta_b*m1*ke_int(a1, l1, m1-1, n1, inorm, a2, l2+1, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (l2)
		      tmp -= twozeta_a*l2*ke_int(a1, l1, m1+1, n1, inorm, a2, l2-1, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (m1 > 0 && l2 > 0)
		      tmp += m1*l2*ke_int(a1, l1, m1-1, n1, inorm, a2, l2-1, m2, 
					       n2, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_ay][coord_bx] += tmp*dens_pf;

		    /*--- d2/dAydBy ---*/
		    tmp = upuppfac*ke_int(a1, l1, m1+1, n1, inorm, a2, l2, m2+1, 
					       n2, jnorm, AB, PA, PB);
		    if (m1)
		      tmp -= twozeta_b*m1*ke_int(a1, l1, m1-1, n1, inorm, a2, l2, m2+1, 
						      n2, jnorm, AB, PA, PB);
		    if (m2)
		      tmp -= twozeta_a*m2*ke_int(a1, l1, m1+1, n1, inorm, a2, l2, m2-1, 
						      n2, jnorm, AB, PA, PB);
		    if (m1 > 0 && m2 > 0)
		      tmp += m1*m2*ke_int(a1, l1, m1-1, n1, inorm, a2, l2, m2-1, 
					       n2, jnorm, AB, PA, PB);
		    if (coord_ay == coord_by)
		      tmp *= 2.0;
		    hess_oe[coord_ay][coord_by] += tmp*dens_pf;

		    /*--- d2/dAydBz ---*/
		    tmp = upuppfac*ke_int(a1, l1, m1+1, n1, inorm, a2, l2, m2, 
					       n2+1, jnorm, AB, PA, PB);
		    if (m1)
		      tmp -= twozeta_b*m1*ke_int(a1, l1, m1-1, n1, inorm, a2, l2, m2, 
						      n2+1, jnorm, AB, PA, PB);
		    if (n2)
		      tmp -= twozeta_a*n2*ke_int(a1, l1, m1+1, n1, inorm, a2, l2, m2, 
						      n2-1, jnorm, AB, PA, PB);
		    if (m1 > 0 && n2 > 0)
		      tmp += m1*n2*ke_int(a1, l1, m1-1, n1, inorm, a2, l2, m2, 
					       n2-1, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_ay][coord_bz] += tmp*dens_pf;

		    /*--- d2/dAzdBx ---*/
		    tmp = upuppfac*ke_int(a1, l1, m1, n1+1, inorm, a2, l2+1, m2, 
					       n2, jnorm, AB, PA, PB);
		    if (n1)
		      tmp -= twozeta_b*n1*ke_int(a1, l1, m1, n1-1, inorm, a2, l2+1, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (l2)
		      tmp -= twozeta_a*l2*ke_int(a1, l1, m1, n1+1, inorm, a2, l2-1, m2, 
						      n2, jnorm, AB, PA, PB);
		    if (n1 > 0 && l2 > 0)
		      tmp += n1*l2*ke_int(a1, l1, m1, n1-1, inorm, a2, l2-1, m2, 
					       n2, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_az][coord_bx] += tmp*dens_pf;

		    /*--- d2/dAzdBy ---*/
		    tmp = upuppfac*ke_int(a1, l1, m1, n1+1, inorm, a2, l2, m2+1, 
					       n2, jnorm, AB, PA, PB);
		    if (n1)
		      tmp -= twozeta_b*n1*ke_int(a1, l1, m1, n1-1, inorm, a2, l2, m2+1, 
						      n2, jnorm, AB, PA, PB);
		    if (m2)
		      tmp -= twozeta_a*m2*ke_int(a1, l1, m1, n1+1, inorm, a2, l2, m2-1, 
						      n2, jnorm, AB, PA, PB);
		    if (n1 > 0 && m2 > 0)
		      tmp += n1*m2*ke_int(a1, l1, m1, n1-1, inorm, a2, l2, m2-1, 
					       n2, jnorm, AB, PA, PB);
		    
		    hess_oe[coord_az][coord_by] += tmp*dens_pf;

		    /*--- d2/dAzdBz ---*/
		    tmp = upuppfac*ke_int(a1, l1, m1, n1+1, inorm, a2, l2, m2, 
					       n2+1, jnorm, AB, PA, PB);
		    if (n1)
		      tmp -= twozeta_b*n1*ke_int(a1, l1, m1, n1-1, inorm, a2, l2, m2, 
						      n2+1, jnorm, AB, PA, PB);
		    if (n2)
		      tmp -= twozeta_a*n2*ke_int(a1, l1, m1, n1+1, inorm, a2, l2, m2, 
						      n2-1, jnorm, AB, PA, PB);
		    if (n1 > 0 && n2 > 0)
		      tmp += n1*n2*ke_int(a1, l1, m1, n1-1, inorm, a2, l2, m2, 
					       n2-1, jnorm, AB, PA, PB);
		    if (coord_az == coord_bz)
		      tmp *= 2.0;
		    hess_oe[coord_az][coord_bz] += tmp*dens_pf;
#endif
		    
		    aj++;
		  }
		}  
		ai++;
	      }
	    } /*--- end cartesian components for (si,sj) with primitives (i,j) ---*/


	    /*---------------------------------------------------------
	      First and second derivative nuclear attraction integrals
	     ---------------------------------------------------------*/
	    for(atom=0;atom<Molecule.num_atoms;atom++) {
	      coord_cx = atom*3;
	      coord_cy = coord_cx + 1;
	      coord_cz = coord_cx + 2;
	      PC.x = sp->P[i][j][0] - Molecule.centers[atom].x;
	      PC.y = sp->P[i][j][1] - Molecule.centers[atom].y;
	      PC.z = sp->P[i][j][2] - Molecule.centers[atom].z;
	      AI_Deriv2_OSrecurs(AI0,AIX,AIY,AIZ,AIXX,AIXY,AIXZ,AIYY,AIYZ,AIZZ,PA,PB,PC,gam,am_i+deriv2_lvl,am_j+deriv2_lvl);

	      /*--- create all am components of si ---*/
	      ai = 0;
	      for(ii = 0; ii <= am_i; ii++){
		l1 = am_i - ii;
		for(jj = 0; jj <= ii; jj++){
		  m1 = ii - jj;
		  n1 = jj;
		  iind = n1*izm1 + m1*iym1 + l1*ixm1;
		  I = si_fao + ai;

		  /*--- create all am components of sj ---*/
		  aj = 0;
		  for(kk = 0; kk <= am_j; kk++){
		    l2 = am_j - kk;
		    for(ll = 0; ll <= kk; ll++){
		      m2 = kk - ll;
		      n2 = ll ;
		      jind = n2*jzm1 + m2*jym1 + l2*jxm1;
		      J = sj_fao + aj;

		      if (si == sj && ai < aj)
		      break;

		      norm_pf = (GTOs.bf_norm[am_i][ai] * GTOs.bf_norm[am_j][aj]);
		      normover_pf = norm_pf * over_pf;
		      dens_pf = Dens[I][J];
		      dens_pf *= norm_pf;

		      if (I != J) {
			dens_pf *= 2.0;
			normover_pf *= 2.0;
		      }

		      zdens_pf = over_pf * dens_pf * Molecule.centers[atom].Z_nuc;
		      znormover_pf = normover_pf * Molecule.centers[atom].Z_nuc;

		      tmp = 2.0*a1*AI0[iind+ixm1][jind][0];
		      if (l1)
			tmp -= l1*AI0[iind-ixm1][jind][0];
		      F[coord_ax][I][J] -= tmp * Molecule.centers[atom].Z_nuc * (normover_pf);

		      tmp = 2.0*a2*AI0[iind][jind+jxm1][0];
		      if (l2)
			tmp -= l2*AI0[iind][jind-jxm1][0];
		      F[coord_bx][I][J] -= tmp * Molecule.centers[atom].Z_nuc * (normover_pf);

		      F[coord_cx][I][J] -= AIX[iind][jind][0] * Molecule.centers[atom].Z_nuc * normover_pf;

		      tmp = 2.0*a1*AI0[iind+iym1][jind][0];
		      if (m1)
			tmp -= m1*AI0[iind-iym1][jind][0];
		      F[coord_ay][I][J] -= tmp * Molecule.centers[atom].Z_nuc * (normover_pf);

		      tmp = 2.0*a2*AI0[iind][jind+jym1][0];
		      if (m2)
			tmp -= m2*AI0[iind][jind-jym1][0];
		      F[coord_by][I][J] -= tmp * Molecule.centers[atom].Z_nuc * (normover_pf);

		      F[coord_cy][I][J] -= AIY[iind][jind][0] * Molecule.centers[atom].Z_nuc * normover_pf;

		      tmp = 2.0*a1*AI0[iind+izm1][jind][0];
		      if (n1)
			tmp -= n1*AI0[iind-izm1][jind][0];
		      F[coord_az][I][J] -= tmp * Molecule.centers[atom].Z_nuc * (normover_pf);

		      tmp = 2.0*a2*AI0[iind][jind+jzm1][0];
		      if (n2)
			tmp -= n2*AI0[iind][jind-jzm1][0];
		      F[coord_bz][I][J] -= tmp * Molecule.centers[atom].Z_nuc * (normover_pf);

		      F[coord_cz][I][J] -= AIZ[iind][jind][0] * Molecule.centers[atom].Z_nuc * normover_pf;

#if INCLUDE_POTENTIAL
		      v_int = AI0[iind][jind][0];
#if INCLUDE_VAA
		      upuppfac = twozeta_a*twozeta_a;
		      /*--- d2/dAx2 ---*/
		      tmp = upuppfac*AI0[iind+ixm1+ixm1][jind][0];
		      updownpfac = l1 + 1; if (l1) updownpfac += l1;
		      tmp -= twozeta_a*updownpfac*v_int;
		      if (l1 >= 2)
			tmp += l1*(l1-1)*AI0[iind-ixm1-ixm1][jind][0];
		      
		      hess_oe[coord_ax][coord_ax] -= tmp * zdens_pf;
		      
		      /*--- d2/dAy2 ---*/
		      tmp = upuppfac*AI0[iind+iym1+iym1][jind][0];
		      updownpfac = m1 + 1; if (m1) updownpfac += m1;
		      tmp -= twozeta_a*updownpfac*v_int;
		      if (m1 >= 2)
			tmp += m1*(m1-1)*AI0[iind-iym1-iym1][jind][0];
		      
		      hess_oe[coord_ay][coord_ay] -= tmp * zdens_pf;
		      
		      /*--- d2/dAz2 ---*/
		      tmp = upuppfac*AI0[iind+izm1+izm1][jind][0];
		      updownpfac = n1 + 1; if (n1) updownpfac += n1;
		      tmp -= twozeta_a*updownpfac*v_int;
		      if (n1 >= 2)
			tmp += n1*(n1-1)*AI0[iind-izm1-izm1][jind][0];
		      
		      hess_oe[coord_az][coord_az] -= tmp * zdens_pf;
		      
		      /*--- d2/dAxdAy ---*/
		      tmp = upuppfac*AI0[iind+ixm1+iym1][jind][0];
		      if (l1)
			tmp -= twozeta_a*l1*AI0[iind-ixm1+iym1][jind][0];
		      if (m1)
			tmp -= twozeta_a*m1*AI0[iind+ixm1-iym1][jind][0];
		      if (l1 > 0 && m1 > 0)
			tmp += l1*m1*AI0[iind-ixm1-iym1][jind][0];
		      
		      hess_oe[coord_ax][coord_ay] -= tmp*zdens_pf;

		      /*--- d2/dAxdAz ---*/
		      tmp = upuppfac*AI0[iind+ixm1+izm1][jind][0];
		      if (l1)
			tmp -= twozeta_a*l1*AI0[iind-ixm1+izm1][jind][0];
		      if (n1)
			tmp -= twozeta_a*n1*AI0[iind+ixm1-izm1][jind][0];
		      if (l1 > 0 && n1 > 0)
			tmp += l1*n1*AI0[iind-ixm1-izm1][jind][0];
		      
		      hess_oe[coord_ax][coord_az] -= tmp*zdens_pf;

		      /*--- d2/dAydAz ---*/
		      tmp = upuppfac*AI0[iind+iym1+izm1][jind][0];
		      if (m1)
			tmp -= twozeta_a*m1*AI0[iind-iym1+izm1][jind][0];
		      if (n1)
			tmp -= twozeta_a*n1*AI0[iind+iym1-izm1][jind][0];
		      if (m1 > 0 && n1 > 0)
			tmp += m1*n1*AI0[iind-iym1-izm1][jind][0];
		      
		      hess_oe[coord_ay][coord_az] -= tmp*zdens_pf;


		      upuppfac = twozeta_b*twozeta_b;
		      /*--- d2/dBx2 ---*/
		      tmp = upuppfac*AI0[iind][jind+jxm1+jxm1][0];
		      updownpfac = l2 + 1; if (l2) updownpfac += l2;
		      tmp -= twozeta_b*updownpfac*v_int;
		      if (l2 >= 2)
			tmp += l2*(l2-1)*AI0[iind][jind-jxm1-jxm1][0];
		      
		      hess_oe[coord_bx][coord_bx] -= tmp * zdens_pf;
		      
		      /*--- d2/dBy2 ---*/
		      tmp = upuppfac*AI0[iind][jind+jym1+jym1][0];
		      updownpfac = m2 + 1; if (m2) updownpfac += m2;
		      tmp -= twozeta_b*updownpfac*v_int;
		      if (m2 >= 2)
			tmp += m2*(m2-1)*AI0[iind][jind-jym1-jym1][0];
		      
		      hess_oe[coord_by][coord_by] -= tmp * zdens_pf;
		      
		      /*--- d2/dBz2 ---*/
		      tmp = upuppfac*AI0[iind][jind+jzm1+jzm1][0];
		      updownpfac = n2 + 1; if (n2) updownpfac += n2;
		      tmp -= twozeta_b*updownpfac*v_int;
		      if (n2 >= 2)
			tmp += n2*(n2-1)*AI0[iind][jind-jzm1-jzm1][0];
		      
		      hess_oe[coord_bz][coord_bz] -= tmp * zdens_pf;
		      
		      /*--- d2/dBxdBy ---*/
		      tmp = upuppfac*AI0[iind][jind+jxm1+jym1][0];
		      if (l2)
			tmp -= twozeta_b*l2*AI0[iind][jind-jxm1+jym1][0];
		      if (m2)
			tmp -= twozeta_b*m2*AI0[iind][jind+jxm1-jym1][0];
		      if (l2 > 0 && m2 > 0)
			tmp += l2*m2*AI0[iind][jind-jxm1-jym1][0];
		      
		      hess_oe[coord_bx][coord_by] -= tmp*zdens_pf;

		      /*--- d2/dBxdBz ---*/
		      tmp = upuppfac*AI0[iind][jind+jxm1+jzm1][0];
		      if (l2)
			tmp -= twozeta_b*l2*AI0[iind][jind-jxm1+jzm1][0];
		      if (n2)
			tmp -= twozeta_b*n2*AI0[iind][jind+jxm1-jzm1][0];
		      if (l2 > 0 && n2 > 0)
			tmp += l2*n2*AI0[iind][jind-jxm1-jzm1][0];
		      
		      hess_oe[coord_bx][coord_bz] -= tmp*zdens_pf;

		      /*--- d2/dBydBz ---*/
		      tmp = upuppfac*AI0[iind][jind+jym1+jzm1][0];
		      if (m2)
			tmp -= twozeta_b*m2*AI0[iind][jind-jym1+jzm1][0];
		      if (n2)
			tmp -= twozeta_b*n2*AI0[iind][jind+jym1-jzm1][0];
		      if (m2 > 0 && n2 > 0)
			tmp += m2*n2*AI0[iind][jind-jym1-jzm1][0];
		      
		      hess_oe[coord_by][coord_bz] -= tmp*zdens_pf;


		      upuppfac = twozeta_a*twozeta_b;
		      /*--- d2/dAxdBx ---*/
		      tmp = upuppfac*AI0[iind+ixm1][jind+jxm1][0];
		      if (l1)
			tmp -= twozeta_b*l1*AI0[iind-ixm1][jind+jxm1][0];
		      if (l2)
			tmp -= twozeta_a*l2*AI0[iind+ixm1][jind-jxm1][0];
		      if (l1 > 0 && l2 > 0)
			tmp += l1*l2*AI0[iind-ixm1][jind-jxm1][0];
		      if (coord_ax == coord_bx)
			tmp *= 2.0;
		      hess_oe[coord_ax][coord_bx] -= tmp*zdens_pf;

		      /*--- d2/dAxdBy ---*/
		      tmp = upuppfac*AI0[iind+ixm1][jind+jym1][0];
		      if (l1)
			tmp -= twozeta_b*l1*AI0[iind-ixm1][jind+jym1][0];
		      if (m2)
			tmp -= twozeta_a*m2*AI0[iind+ixm1][jind-jym1][0];
		      if (l1 > 0 && m2 > 0)
			tmp += l1*m2*AI0[iind-ixm1][jind-jym1][0];
		      hess_oe[coord_ax][coord_by] -= tmp*zdens_pf;

		      /*--- d2/dAxdBz ---*/
		      tmp = upuppfac*AI0[iind+ixm1][jind+jzm1][0];
		      if (l1)
			tmp -= twozeta_b*l1*AI0[iind-ixm1][jind+jzm1][0];
		      if (n2)
			tmp -= twozeta_a*n2*AI0[iind+ixm1][jind-jzm1][0];
		      if (l1 > 0 && n2 > 0)
			tmp += l1*n2*AI0[iind-ixm1][jind-jzm1][0];
		      hess_oe[coord_ax][coord_bz] -= tmp*zdens_pf;

		      /*--- d2/dAydBx ---*/
		      tmp = upuppfac*AI0[iind+iym1][jind+jxm1][0];
		      if (m1)
			tmp -= twozeta_b*m1*AI0[iind-iym1][jind+jxm1][0];
		      if (l2)
			tmp -= twozeta_a*l2*AI0[iind+iym1][jind-jxm1][0];
		      if (m1 > 0 && l2 > 0)
			tmp += m1*l2*AI0[iind-iym1][jind-jxm1][0];
		      hess_oe[coord_ay][coord_bx] -= tmp*zdens_pf;

		      /*--- d2/dAydBy ---*/
		      tmp = upuppfac*AI0[iind+iym1][jind+jym1][0];
		      if (m1)
			tmp -= twozeta_b*m1*AI0[iind-iym1][jind+jym1][0];
		      if (m2)
			tmp -= twozeta_a*m2*AI0[iind+iym1][jind-jym1][0];
		      if (m1 > 0 && m2 > 0)
			tmp += m1*m2*AI0[iind-iym1][jind-jym1][0];
		      if (coord_ay == coord_by)
			tmp *= 2.0;
		      hess_oe[coord_ay][coord_by] -= tmp*zdens_pf;

		      /*--- d2/dAydBz ---*/
		      tmp = upuppfac*AI0[iind+iym1][jind+jzm1][0];
		      if (m1)
			tmp -= twozeta_b*m1*AI0[iind-iym1][jind+jzm1][0];
		      if (n2)
			tmp -= twozeta_a*n2*AI0[iind+iym1][jind-jzm1][0];
		      if (m1 > 0 && n2 > 0)
			tmp += m1*n2*AI0[iind-iym1][jind-jzm1][0];
		      hess_oe[coord_ay][coord_bz] -= tmp*zdens_pf;

		      /*--- d2/dAzdBx ---*/
		      tmp = upuppfac*AI0[iind+izm1][jind+jxm1][0];
		      if (n1)
			tmp -= twozeta_b*n1*AI0[iind-izm1][jind+jxm1][0];
		      if (l2)
			tmp -= twozeta_a*l2*AI0[iind+izm1][jind-jxm1][0];
		      if (n1 > 0 && l2 > 0)
			tmp += n1*l2*AI0[iind-izm1][jind-jxm1][0];
		      hess_oe[coord_az][coord_bx] -= tmp*zdens_pf;

		      /*--- d2/dAzdBy ---*/
		      tmp = upuppfac*AI0[iind+izm1][jind+jym1][0];
		      if (n1)
			tmp -= twozeta_b*n1*AI0[iind-izm1][jind+jym1][0];
		      if (m2)
			tmp -= twozeta_a*m2*AI0[iind+izm1][jind-jym1][0];
		      if (n1 > 0 && m2 > 0)
			tmp += n1*m2*AI0[iind-izm1][jind-jym1][0];
		      hess_oe[coord_az][coord_by] -= tmp*zdens_pf;

		      /*--- d2/dAzdBz ---*/
		      tmp = upuppfac*AI0[iind+izm1][jind+jzm1][0];
		      if (n1)
			tmp -= twozeta_b*n1*AI0[iind-izm1][jind+jzm1][0];
		      if (n2)
			tmp -= twozeta_a*n2*AI0[iind+izm1][jind-jzm1][0];
		      if (n1 > 0 && n2 > 0)
			tmp += n1*n2*AI0[iind-izm1][jind-jzm1][0];
		      if (coord_az == coord_bz)
			tmp *= 2.0;
		      hess_oe[coord_az][coord_bz] -= tmp*zdens_pf;
#endif

#if INCLUDE_VCA
		      /*--- d2/dCxdAx ---*/
		      tmp = twozeta_a*AIX[iind+ixm1][jind][0];
		      if (l1)
			tmp -= l1*AIX[iind-ixm1][jind][0];
		      if (coord_ax == coord_cx) tmp *= 2.0;
		      hess_oe[coord_ax][coord_cx] -= tmp*zdens_pf;

		      /*--- d2/dCxdAy ---*/
		      tmp = twozeta_a*AIX[iind+iym1][jind][0];
		      if (m1)
			tmp -= m1*AIX[iind-iym1][jind][0];
		      hess_oe[coord_ay][coord_cx] -= tmp*zdens_pf;

		      /*--- d2/dCxdAz ---*/
		      tmp = twozeta_a*AIX[iind+izm1][jind][0];
		      if (n1)
			tmp -= n1*AIX[iind-izm1][jind][0];
		      hess_oe[coord_az][coord_cx] -= tmp*zdens_pf;

		      /*--- d2/dCydAx ---*/
		      tmp = twozeta_a*AIY[iind+ixm1][jind][0];
		      if (l1)
			tmp -= l1*AIY[iind-ixm1][jind][0];
		      hess_oe[coord_ax][coord_cy] -= tmp*zdens_pf;

		      /*--- d2/dCydAy ---*/
		      tmp = twozeta_a*AIY[iind+iym1][jind][0];
		      if (m1)
			tmp -= m1*AIY[iind-iym1][jind][0];
		      if (coord_ay == coord_cy) tmp *= 2.0;
		      hess_oe[coord_ay][coord_cy] -= tmp*zdens_pf;

		      /*--- d2/dCydAz ---*/
		      tmp = twozeta_a*AIY[iind+izm1][jind][0];
		      if (n1)
			tmp -= n1*AIY[iind-izm1][jind][0];
		      hess_oe[coord_az][coord_cy] -= tmp*zdens_pf;

		      /*--- d2/dCzdAx ---*/
		      tmp = twozeta_a*AIZ[iind+ixm1][jind][0];
		      if (l1)
			tmp -= l1*AIZ[iind-ixm1][jind][0];
		      hess_oe[coord_ax][coord_cz] -= tmp*zdens_pf;

		      /*--- d2/dCzdAy ---*/
		      tmp = twozeta_a*AIZ[iind+iym1][jind][0];
		      if (m1)
			tmp -= m1*AIZ[iind-iym1][jind][0];
		      hess_oe[coord_ay][coord_cz] -= tmp*zdens_pf;

		      /*--- d2/dCzdAz ---*/
		      tmp = twozeta_a*AIZ[iind+izm1][jind][0];
		      if (n1)
			tmp -= n1*AIZ[iind-izm1][jind][0];
		      if (coord_az == coord_cz) tmp *= 2.0;
		      hess_oe[coord_az][coord_cz] -= tmp*zdens_pf;


		      /*--- d2/dCxdBx ---*/
		      tmp = twozeta_b*AIX[iind][jind+jxm1][0];
		      if (l2)
			tmp -= l2*AIX[iind][jind-jxm1][0];
		      if (coord_bx == coord_cx) tmp *= 2.0;
		      hess_oe[coord_bx][coord_cx] -= tmp*zdens_pf;

		      /*--- d2/dCxdBy ---*/
		      tmp = twozeta_b*AIX[iind][jind+jym1][0];
		      if (m2)
			tmp -= m2*AIX[iind][jind-jym1][0];
		      hess_oe[coord_by][coord_cx] -= tmp*zdens_pf;

		      /*--- d2/dCxdBz ---*/
		      tmp = twozeta_b*AIX[iind][jind+jzm1][0];
		      if (n2)
			tmp -= n2*AIX[iind][jind-jzm1][0];
		      hess_oe[coord_bz][coord_cx] -= tmp*zdens_pf;

		      /*--- d2/dCydBx ---*/
		      tmp = twozeta_b*AIY[iind][jind+jxm1][0];
		      if (l2)
			tmp -= l2*AIY[iind][jind-jxm1][0];
		      hess_oe[coord_bx][coord_cy] -= tmp*zdens_pf;

		      /*--- d2/dCydBy ---*/
		      tmp = twozeta_b*AIY[iind][jind+jym1][0];
		      if (m2)
			tmp -= m2*AIY[iind][jind-jym1][0];
		      if (coord_by == coord_cy) tmp *= 2.0;
		      hess_oe[coord_by][coord_cy] -= tmp*zdens_pf;

		      /*--- d2/dCydBz ---*/
		      tmp = twozeta_b*AIY[iind][jind+jzm1][0];
		      if (n2)
			tmp -= n2*AIY[iind][jind-jzm1][0];
		      hess_oe[coord_bz][coord_cy] -= tmp*zdens_pf;

		      /*--- d2/dCzdBx ---*/
		      tmp = twozeta_b*AIZ[iind][jind+jxm1][0];
		      if (l2)
			tmp -= l2*AIZ[iind][jind-jxm1][0];
		      hess_oe[coord_bx][coord_cz] -= tmp*zdens_pf;

		      /*--- d2/dCzdBy ---*/
		      tmp = twozeta_b*AIZ[iind][jind+jym1][0];
		      if (m2)
			tmp -= m2*AIZ[iind][jind-jym1][0];
		      hess_oe[coord_by][coord_cz] -= tmp*zdens_pf;

		      /*--- d2/dCzdBz ---*/
		      tmp = twozeta_b*AIZ[iind][jind+jzm1][0];
		      if (n2)
			tmp -= n2*AIZ[iind][jind-jzm1][0];
		      if (coord_bz == coord_cz) tmp *= 2.0;
		      hess_oe[coord_bz][coord_cz] -= tmp*zdens_pf;
#endif

#if INCLUDE_VCC		      
		      /*---- d2/dCx2 ---*/
		      hess_oe[coord_cx][coord_cx] -= AIXX[iind][jind][0]*zdens_pf;
		      /*---- d2/dCy2 ---*/
		      hess_oe[coord_cy][coord_cy] -= AIYY[iind][jind][0]*zdens_pf;
		      /*---- d2/dCz2 ---*/
		      hess_oe[coord_cz][coord_cz] -= AIZZ[iind][jind][0]*zdens_pf;
		      /*---- d2/dCxdCy ---*/
		      hess_oe[coord_cx][coord_cy] -= AIXY[iind][jind][0]*zdens_pf;
		      /*---- d2/dCxdCz ---*/
		      hess_oe[coord_cx][coord_cz] -= AIXZ[iind][jind][0]*zdens_pf;
		      /*---- d2/dCydCz ---*/
		      hess_oe[coord_cy][coord_cz] -= AIYZ[iind][jind][0]*zdens_pf;
#endif
#endif
		      aj++;
		    }
		  }  
		  ai++;
		}
	      } /*--- end cartesian components for (si,sj) with primitives (i,j) ---*/
	    }
	  }
	} /*--- end primitive contraction ---*/
    

    }
  }

  symmetrize_hessian(hess_ov);
  symmetrize_hessian(hess_oe);
  if (UserOptions.print_lvl >= PRINT_OEDERIV) {
    print_atommat("Overlap component of the molecular Hessian (a.u.)",hess_ov);
    print_atommat("Core Hamiltonian component of the molecular Hessian (a.u.)",hess_oe);
  }

  add_mat(Hess,hess_ov,Hess,Molecule.num_atoms*3,Molecule.num_atoms*3);
  add_mat(Hess,hess_oe,Hess,Molecule.num_atoms*3,Molecule.num_atoms*3);

  /*--- flush it all away ---*/
  fflush(outfile);
  free_box(AI0,indmax,indmax);
  free_box(AIX,indmax,indmax);
  free_box(AIY,indmax,indmax);
  free_box(AIZ,indmax,indmax);
  free_box(AIXX,indmax,indmax);
  free_box(AIXY,indmax,indmax);
  free_box(AIXZ,indmax,indmax);
  free_box(AIYY,indmax,indmax);
  free_box(AIYZ,indmax,indmax);
  free_box(AIZZ,indmax,indmax);
  free_block(hess_ov);
  free_block(hess_oe);
}   


/*!-----------------------------------------------
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


/*!-----------------------------------------------------------------------------
  This computes matrix element of kinetic energy between 2 primitive gaussians
 -----------------------------------------------------------------------------*/
double ke_int(double a1, int l1, int m1, int n1, double norm1,
	      double a2, int l2, int m2, int n2, double norm2,
	      struct coordinates AB,
	      struct coordinates PA,
	      struct coordinates PB)
{
  int localprint = 0 ;
  double Ix, Iy, Iz;
  double I1 = 0.0, I2 = 0.0, I3 = 0.0, I4 = 0.0 ;
  double gam;
  double norm_fact;

  I2 = overlap_int(a1, l1+1, m1, n1, norm1, a2, l2+1, m2, n2, norm2, 
                    AB, PA, PB);
  I1 = overlap_int(a1, l1-1, m1, n1, norm1, a2, l2-1, m2, n2, norm2, 
                    AB, PA, PB);
  I3 = overlap_int(a1, l1+1, m1, n1, norm1, a2, l2-1, m2, n2, norm2, 
                    AB, PA, PB);
  I4 = overlap_int(a1, l1-1, m1, n1, norm1, a2, l2+1, m2, n2, norm2, 
                    AB, PA, PB);
 
  Ix = (l1*l2*I1/2.0 + 2*a1*a2*I2 - a1*l2*I3 - a2*l1*I4);

  I1 = 0.0;
  I2 = 0.0;
  I3 = 0.0;
  I4 = 0.0;

  I2 = overlap_int(a1, l1, m1+1, n1, norm1, a2, l2, m2+1, n2, norm2, 
                    AB, PA, PB);
  I1 = overlap_int(a1, l1, m1-1, n1, norm1, a2, l2, m2-1, n2, norm2, 
                    AB, PA, PB);
  I3 = overlap_int(a1, l1, m1+1, n1, norm1, a2, l2, m2-1, n2, norm2, 
                    AB, PA, PB);
  I4 = overlap_int(a1, l1, m1-1, n1, norm1, a2, l2, m2+1, n2, norm2, 
                    AB, PA, PB);
  
  Iy = (m1*m2*I1/2.0 + 2*a1*a2*I2 - a1*m2*I3 - a2*m1*I4);

  I1 = 0.0;
  I2 = 0.0;
  I3 = 0.0;
  I4 = 0.0;

  I2 = overlap_int(a1, l1, m1, n1+1, norm1, a2, l2, m2, n2+1, norm2, 
                    AB, PA, PB);
  I1 = overlap_int(a1, l1, m1, n1-1, norm1, a2, l2, m2, n2-1, norm2, 
                    AB, PA, PB);
  I3 = overlap_int(a1, l1, m1, n1+1, norm1, a2, l2, m2, n2-1, norm2, 
                    AB, PA, PB);
  I4 = overlap_int(a1, l1, m1, n1-1, norm1, a2, l2, m2, n2+1, norm2, 
                    AB, PA, PB);
  
  Iz = (n1*n2*I1/2.0 + 2*a1*a2*I2 - a1*n2*I3 - a2*n1*I4);

  return((Ix+Iy+Iz));

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
