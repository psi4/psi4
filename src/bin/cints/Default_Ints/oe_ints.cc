/*! \file oe_ints.cc
    \ingroup (CINTS)
    \brief One-electron integrals.
*/
#include<cmath>
#include<cstring>
#include<stdio.h>
#include<stdlib.h>
#include<libipv1/ip_lib.h>
#include<libiwl/iwl.h>
#include<libciomr/libciomr.h>
#include<libint/libint.h>
#include<psifiles.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"oe_osrr.h"
#include"small_fns.h"

#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#endif

namespace psi { namespace CINTS {
/*--- These frequently used numbers are to avoid costs of passing parameters ---*/
static double oo2g, oog, gam;

/*--------------------------------------------------------------------------
  This function computes one-electron integrals and writes them out to disk
 --------------------------------------------------------------------------*/
void oe_ints()
{
   struct coordinates PA, PB, AB, PC;
   struct shell_pair *sp;
   struct unique_shell_pair *usp;
   register int i, j, k, l, ii, jj, kk, ll;
   int count;
   int si, sj, sjj, usi, usj;
   int np_i, np_j;
   int so_i, so_j, so_ij, bf_i, bf_j;
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
   int indmax,iind,jind;
   int atom;
   int bf;
   int ud, dcr_ij, stab_i, stab_j;
   int *R_list;
   int R;
   int lambda_R;
   int *sj_arr;
   int iimax, jjmax;
   int num_unique_doublets;
   double a1, a2;
   double ab2;
   double x0, y0, z0;
   double tx, ty, tz;
   double ***Stmp, ***Ttmp, ***Vtmp;
   double **stemp, **ttemp, **vtemp, **temp, **tmp_dptr;
   double *S, *T, *V;
   double inorm, jnorm, over_pf;
   double *ptr1, *ptr2, norm1, norm12;
   double ***AI0;
   double **OIX, **OIY, **OIZ;

#ifdef USE_TAYLOR_FM
  init_Taylor_Fm_Eval(BasisSet.max_am*4-4,UserOptions.cutoff);
#endif  

  /*--- allocate room for the one-e matrices ---*/
  dimension = ioff[Symmetry.num_so];
  S = init_array(dimension);
  V = init_array(dimension);
  T = init_array(dimension);

  /*--- allocate storage for shell blocks of one electron integrals ---*/
  dimension = ioff[BasisSet.max_am];
  if (Symmetry.nirreps > 1) { /*--- Non-C1 symmetry case ---*/
    Stmp = (double ***) malloc(Symmetry.nirreps*sizeof(double **));
    Ttmp = (double ***) malloc(Symmetry.nirreps*sizeof(double **));
    Vtmp = (double ***) malloc(Symmetry.nirreps*sizeof(double **));
    for(i=0;i<Symmetry.nirreps;i++) {
      Stmp[i] = init_matrix(dimension,dimension);
      Ttmp[i] = init_matrix(dimension,dimension);
      Vtmp[i] = init_matrix(dimension,dimension);
    }
  }
  else { /*--- C1 symmetry case ---*/
    stemp = init_matrix(dimension,dimension);
    ttemp = init_matrix(dimension,dimension);
    vtemp = init_matrix(dimension,dimension);
  }
  sj_arr = (int *)malloc(sizeof(int)*Symmetry.nirreps);

  if (BasisSet.puream) /*--- allocate a scratch matrix for transfromation
		    into puream basis ---*/
    temp = init_matrix(dimension,dimension);

  indmax = (BasisSet.max_am-1)*BasisSet.max_am*BasisSet.max_am+1;
  AI0 = init_box(indmax,indmax,2*BasisSet.max_am+1);
  OIX = block_matrix(BasisSet.max_am+2,BasisSet.max_am+2);
  OIY = block_matrix(BasisSet.max_am+2,BasisSet.max_am+2);
  OIZ = block_matrix(BasisSet.max_am+2,BasisSet.max_am+2);
  
  for (usi=0; usi<Symmetry.num_unique_shells; usi++){
    si = Symmetry.us2s[usi];
    am_i = BasisSet.shells[si].am-1;
    izm = 1;
    iym = am_i+1;
    ixm = iym*iym;
    for (usj=0; usj<=usi; usj++){
      if (Symmetry.nirreps > 1) { /*--- Non-C1 symmetry case ---*/
        usp = &(Symmetry.us_pairs[usi][usj]); 
	sjj = Symmetry.us2s[usj];
	stab_i = Symmetry.atom_positions[BasisSet.shells[si].center-1];
	stab_j = Symmetry.atom_positions[BasisSet.shells[sjj].center-1];
	ni = ioff[BasisSet.shells[si].am];
	nj = ioff[BasisSet.shells[sjj].am];
	am_j = BasisSet.shells[sjj].am-1;
	jzm = 1;
	jym = am_j+1;
	jxm = jym*jym;
	R_list = Symmetry.dcr[stab_i][stab_j];
	lambda_R = Symmetry.nirreps/Symmetry.dcr_deg[stab_i][stab_j];
	memset(sj_arr,0,sizeof(int)*Symmetry.nirreps);
	count = 0;
	/*--- generate petite list ---*/
	for(dcr_ij=0;dcr_ij<Symmetry.dcr_dim[stab_i][stab_j];dcr_ij++){
	  R = R_list[dcr_ij];
	  sj = BasisSet.shells[sjj].trans_vec[R]-1;
	  sj_arr[count] = sj;
	  count++;
	}
	num_unique_doublets = count;
      }
      else { /*--- C1 symmetry case ---*/
	num_unique_doublets = 1;
	sj_arr[0] = usj;
	ni = ioff[BasisSet.shells[usi].am];
	nj = ioff[BasisSet.shells[usj].am];
	am_j = BasisSet.shells[usj].am-1;
	jzm = 1;
	jym = am_j+1;
	jxm = jym*jym;
	ioffset = BasisSet.shells[usi].fbf - 1;
	joffset = BasisSet.shells[usj].fbf - 1;
      }
	
      for(ud=0;ud<num_unique_doublets;ud++) {
	sj = sj_arr[ud];

	if (Symmetry.nirreps > 1) {
	  stemp = Stmp[ud];
	  ttemp = Ttmp[ud];
	  vtemp = Vtmp[ud];
	}

	sp = &(BasisSet.shell_pairs[si][sj]);
	AB.x = sp->AB[0];
	AB.y = sp->AB[1];
	AB.z = sp->AB[2];
	ab2 = AB.x * AB.x;
	ab2 += AB.y * AB.y;
	ab2 += AB.z * AB.z;
	
	/*--- zero the temporary storage for accumulating contractions ---*/
	for(i=0;i<dimension;i++) {
	  memset(stemp[i],0,sizeof(double)*dimension);
	  memset(ttemp[i],0,sizeof(double)*dimension);
	  memset(vtemp[i],0,sizeof(double)*dimension);
	}
      
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

	    OI_OSrecurs(OIX,OIY,OIZ,PA,PB,gam,am_i+2,am_j+2);

	    /*--- create all am components of si ---*/
	    ai = 0;
	    for(ii = 0; ii <= am_i; ii++){
	      l1 = am_i - ii;
	      for(jj = 0; jj <= ii; jj++){
		m1 = ii - jj;
		n1 = jj ;
		/*--- create all am components of sj ---*/
		aj = 0;
		for(kk = 0; kk <= am_j; kk++){
		  l2 = am_j - kk;
		  for(ll = 0; ll <= kk; ll++){
		    m2 = kk - ll;
		    n2 = ll ;

		    x0 = OIX[l1][l2];  y0 = OIY[m1][m2];  z0 = OIZ[n1][n2];
		    stemp[ai][aj] += over_pf*x0*y0*z0;
		    tx = a2*(2*l2+1)*OIX[l1][l2] - 2*a2*a2*OIX[l1][l2+2];
		    if (l2 >= 2)
		      tx -= 0.5*l2*(l2-1)*OIX[l1][l2-2];
		    ty = a2*(2*m2+1)*OIY[m1][m2] - 2*a2*a2*OIY[m1][m2+2];
		    if (m2 >= 2)
		      ty -= 0.5*m2*(m2-1)*OIY[m1][m2-2];
		    tz = a2*(2*n2+1)*OIZ[n1][n2] - 2*a2*a2*OIZ[n1][n2+2];
		    if (n2 >= 2)
		      tz -= 0.5*n2*(n2-1)*OIZ[n1][n2-2];
		    ttemp[ai][aj] += over_pf*(tx*y0*z0 + x0*ty*z0 + x0*y0*tz);

		    aj++;
		  }
		}  
		ai++;
	      }
	    } /*--- end cartesian components for (si,sj) with primitives (i,j) ---*/

      /*--- create all am components of si ---*/
	    for(atom=0;atom<Molecule.num_atoms;atom++) {
	      PC.x = sp->P[i][j][0] - Molecule.centers[atom].x;
	      PC.y = sp->P[i][j][1] - Molecule.centers[atom].y;
	      PC.z = sp->P[i][j][2] - Molecule.centers[atom].z;
	      AI_OSrecurs(AI0,PA,PB,PC,gam,am_i,am_j);
	      ai = 0;
	      for(ii = 0; ii <= am_i; ii++){
		l1 = am_i - ii;
		for(jj = 0; jj <= ii; jj++){
		  m1 = ii - jj;
		  n1 = jj ;
		  iind = n1*izm + m1*iym + l1*ixm;
		  /*--- create all am components of sj ---*/
		  aj = 0;
		  for(kk = 0; kk <= am_j; kk++){
		    l2 = am_j - kk;
		    for(ll = 0; ll <= kk; ll++){
		      m2 = kk - ll;
		      n2 = ll ;

		      jind = n2*jzm + m2*jym + l2*jxm;

		      vtemp[ai][aj] += -AI0[iind][jind][0] * 
		                       Molecule.centers[atom].Z_nuc * over_pf;

		      aj++;
		    }
		  }  
		  ai++;
		}
	      } /*--- end cartesian components for (si,sj) with primitives (i,j) ---*/
	    }
	  }
	} /*--- end primitive contraction ---*/
    
      /*--- Normalize the contracted integrals ---*/
      ptr1 = GTOs.bf_norm[am_i];
      ptr2 = GTOs.bf_norm[am_j];
      for(i=0; i<ni; i++) {
	norm1 = ptr1[i];
	for(j=0; j<nj; j++) {
	  norm12 = norm1*ptr2[j];
	  stemp[i][j] *= norm12;
	  ttemp[i][j] *= norm12;
	  vtemp[i][j] *= norm12;
	}
      }

      /*--- if puream - transform integrals into puream basis ---*/
      li = ni;
      if (BasisSet.puream && am_i > 0) {
	li = 2*am_i + 1;
	mmult(GTOs.cart2pureang[am_i],0,stemp,0,temp,0,li,ni,nj,0);
	for(i=0;i<li;i++)
	  for(j=0;j<nj;j++)
	    stemp[i][j] = temp[i][j];
	mmult(GTOs.cart2pureang[am_i],0,ttemp,0,temp,0,li,ni,nj,0);
	for(i=0;i<li;i++)
	  for(j=0;j<nj;j++)
	    ttemp[i][j] = temp[i][j];
	mmult(GTOs.cart2pureang[am_i],0,vtemp,0,temp,0,li,ni,nj,0);
	for(i=0;i<li;i++)
	  for(j=0;j<nj;j++)
	    vtemp[i][j] = temp[i][j];
      }
      lj = nj;
      if (BasisSet.puream && am_j > 0) {
	lj = 2*am_j + 1;
	mmult(stemp,0,GTOs.cart2pureang[am_j],1,temp,0,li,nj,lj,0);
	for(i=0;i<li;i++)
	  for(j=0;j<lj;j++)
	    stemp[i][j] = temp[i][j];
	mmult(ttemp,0,GTOs.cart2pureang[am_j],1,temp,0,li,nj,lj,0);
	for(i=0;i<li;i++)
	  for(j=0;j<lj;j++)
	    ttemp[i][j] = temp[i][j];
	mmult(vtemp,0,GTOs.cart2pureang[am_j],1,temp,0,li,nj,lj,0);
	for(i=0;i<li;i++)
	  for(j=0;j<lj;j++)
	    vtemp[i][j] = temp[i][j];
      }
      } /*--- All unique shell pairs are done ---*/

      ni = li;
      nj = lj;
      if (Symmetry.nirreps > 1) {
	ioffset = BasisSet.shells[si].fbf-1;
	for(ij=0;ij<usp->SOpair_npi[0];ij++) {
	  i = usp->SOpair_bf_i[0][ij];
	  j = usp->SOpair_bf_j[0][ij];
	  so_i = usp->SOpair_so_i[0][ij];
	  so_j = usp->SOpair_so_j[0][ij];
	  so_ij = INDEX(so_i,so_j);
	  S[so_ij] = T[so_ij] = 0.0;
	  for(ud=0;ud<num_unique_doublets;ud++){
	    sj = sj_arr[ud];
	    bf_i = ioffset + i;
	    bf_j = BasisSet.shells[sj].fbf-1+j;
	    S[so_ij] += (Symmetry.usotao[so_i][bf_i]*
	                 Symmetry.usotao[so_j][bf_j])*
	                 Stmp[ud][i][j]*lambda_R;
	    T[so_ij] += (Symmetry.usotao[so_i][bf_i]*
	                 Symmetry.usotao[so_j][bf_j])*
	                 Ttmp[ud][i][j]*lambda_R;
	    V[so_ij] += (Symmetry.usotao[so_i][bf_i]*
	                 Symmetry.usotao[so_j][bf_j])*
	                 Vtmp[ud][i][j]*lambda_R;
	  }
	}
      }
      else /*--- C1 symmetry case ---*/
	for(i=0;i<ni;i++)
	  for(j=0;j<nj;j++) {
	    ij = INDEX(ioffset + i, joffset + j);
	    S[ij] = stemp[i][j];
	    T[ij] = ttemp[i][j];
	    V[ij] = vtemp[i][j];
	  }
    }
  }  /*--- This unique shell pair is done ---*/

  /*--- flush it all away ---*/
  fflush(outfile);
  free_block(OIX);
  free_block(OIY);
  free_block(OIZ);
  free_box(AI0,indmax,indmax);
  if (Symmetry.nirreps > 1) {
    for(i=0;i<Symmetry.nirreps;i++) {
      free_matrix(Stmp[i],dimension);
      free_matrix(Ttmp[i],dimension);
      free_matrix(Vtmp[i],dimension);
    }
    free(Stmp);
    free(Ttmp);
    free(Vtmp);
  }
  if (BasisSet.puream)
    free_matrix(temp,dimension); 
  dimension=ioff[Symmetry.num_so];
  iwl_wrtone(IOUnits.itapS,PSIF_SO_S,dimension,S);
  iwl_wrtone(IOUnits.itapT,PSIF_SO_T,dimension,T);
  iwl_wrtone(IOUnits.itapV,PSIF_SO_V,dimension,V);
  if (UserOptions.print_lvl >= PRINT_OEI) {
    fprintf(outfile,"  -Overlap integrals:\n\n");
    print_array(S,Symmetry.num_so,outfile);
    fprintf(outfile,"\n  -Kinetic energy integrals:\n\n");
    print_array(T,Symmetry.num_so,outfile);
    fprintf(outfile,"\n  -Nuclear attraction energy integrals:\n\n");
    print_array(V,Symmetry.num_so,outfile);
    fprintf(outfile,"\n");
  }
  free(S);
  free(T);
  free(V);

#ifdef USE_TAYLOR_FM
  free_Taylor_Fm_Eval();
#endif

  return;
}   
};};
