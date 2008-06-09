/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>

#include "input.h"
#include "defines.h"
#define EXTERN
#include "global.h"

namespace psi { namespace input {

using namespace psi;

/*--- External function declaration ---*/
extern void OI_OSrecurs(double **OIX, double **OIY, double **OIZ, struct coordinates PA, struct coordinates PB,
			double gamma, int lmaxi, int lmaxj);

/*------------------------------------------------------
  This function computes the overlap matrix S12 between
  the new and old bases (in SO basis)
 ------------------------------------------------------*/
double **overlap_new_old()
{
  struct coordinates A, B, P, PA, PB;

  int ii, jj, kk, ll;
  int maxsz1, maxsz2;
  int x1, x2, y1, y2, z1, z2;
  int s1, s2, ao1, ao2, l1, l2, n1, n2, p1, p2,
    prim1, prim2, fprim1, fprim2, nprim1, nprim2, atom1, atom2, bf1, bf2;
  double exp1, exp2, cc1, cc2, norm1;
  double ab2, pfac, gam, oog;
  double *ptr1, *ptr2;
  double **stemp, **S12_AO, **S12, **OIX, **OIY, **OIZ;
  double **tmpmat;

  S12_AO = block_matrix(num_ao, Oldcalc.num_ao);
  maxsz1 = ioff[max_angmom+1];
  maxsz2 = ioff[Oldcalc.max_angmom+1];
  stemp = block_matrix(maxsz1,maxsz2);

  OIX = block_matrix(max_angmom+1,Oldcalc.max_angmom+1);
  OIY = block_matrix(max_angmom+1,Oldcalc.max_angmom+1);
  OIZ = block_matrix(max_angmom+1,Oldcalc.max_angmom+1);
  
  /*--- Loop over shells ---*/
  for(s1=0;s1<num_shells;s1++) {
    ao1 = first_ao_shell[s1];
    nprim1 = nprim_in_shell[s1];
    fprim1 = first_prim_shell[s1];
    bf1 = first_ao_shell[s1];
    l1 = shell_ang_mom[s1];
    n1 = ioff[l1+1];
    atom1 = shell_nucleus[s1];
    A.x = geometry[atom1][0];
    A.y = geometry[atom1][1];
    A.z = geometry[atom1][2];
    for(s2=0;s2<Oldcalc.num_shells;s2++) {
      ao2 = Oldcalc.first_ao_shell[s2];
      nprim2 = Oldcalc.nprim_in_shell[s2];
      fprim2 = Oldcalc.first_prim_shell[s2];
      bf2 = Oldcalc.first_ao_shell[s2];
      l2 = Oldcalc.shell_ang_mom[s2];
      n2 = ioff[l2+1];
      atom2 = Oldcalc.shell_nucleus[s2];
      B.x = Oldcalc.geometry[atom2][0];
      B.y = Oldcalc.geometry[atom2][1];
      B.z = Oldcalc.geometry[atom2][2];
	  
      ab2  = (A.x-B.x)*(A.x-B.x);
      ab2 += (A.y-B.y)*(A.y-B.y);
      ab2 += (A.z-B.z)*(A.z-B.z);

      /*--- zero the temporary storage for accumulating contractions ---*/
      memset(stemp[0],0,sizeof(double)*maxsz1*maxsz2);
	  
      /*--- Loop over primitives here ---*/
      for(p1=0,prim1=fprim1;p1<nprim1;p1++,prim1++) {
	exp1 = exponents[prim1];
	cc1 = contr_coeff[prim1];
	for(p2=0,prim2=fprim2;p2<nprim2;p2++,prim2++) {
	  exp2 = Oldcalc.exponents[prim2];
	  cc2 = Oldcalc.contr_coeff[prim2][l2];

	  gam = exp1 + exp2;
	  oog = 1.0/gam;
	  P.x = (exp1*A.x + exp2*B.x)*oog;
	  P.y = (exp1*A.y + exp2*B.y)*oog;
	  P.z = (exp1*A.z + exp2*B.z)*oog;
	  PA.x = P.x - A.x;
	  PA.y = P.y - A.y;
	  PA.z = P.z - A.z;
	  PB.x = P.x - B.x;
	  PB.y = P.y - B.y;
	  PB.z = P.z - B.z;

	  pfac = exp(-exp1*exp2*ab2*oog)*sqrt(M_PI*oog)*M_PI*oog*cc1*cc2;

	  OI_OSrecurs(OIX,OIY,OIZ,PA,PB,gam,l1,l2);

	  /*--- create all am components of si ---*/
	  ao1 = 0;
	  for(ii = 0; ii <= l1; ii++){
	    x1 = l1 - ii;
	    for(jj = 0; jj <= ii; jj++){
	      y1 = ii - jj;
	      z1 = jj;
	      /*--- create all am components of sj ---*/
	      ao2 = 0;
	      for(kk = 0; kk <= l2; kk++){
		x2 = l2 - kk;
		for(ll = 0; ll <= kk; ll++){
		  y2 = kk - ll;
		  z2 = ll;
				  
		  stemp[ao1][ao2] += pfac*OIX[x1][x2]*OIY[y1][y2]*OIZ[z1][z2];
				  
		  ao2++;
		}
	      }
	      ao1++;
	    }
	  } /*--- end cartesian components for (si,sj) with primitives (i,j) ---*/
	}
      } /*--- end primitive contraction ---*/

      /*--- Normalize the contracted integrals and put them into S12 ---*/
      ptr1 = GTOs.bf_norm[l1];
      ptr2 = GTOs.bf_norm[l2];
      for(ii=0; ii<n1; ii++) {
	norm1 = ptr1[ii];
	for(jj=0; jj<n2; jj++) {
	  S12_AO[bf1+ii][bf2+jj] = stemp[ii][jj]*norm1*ptr2[jj];
	}
      }
	  
    }
  }   /*--- This shell pair is done ---*/

  free_block(OIX);
  free_block(OIY);
  free_block(OIZ);
  free_block(stemp);

  /*--- transform to SO basis ---*/
  tmpmat = block_matrix(num_so,Oldcalc.num_ao);
  mmult(usotao,0,S12_AO,0,tmpmat,0,num_so,num_ao,Oldcalc.num_ao,0);
  free_block(S12_AO);
  S12 = block_matrix(num_so,Oldcalc.num_so);
  mmult(tmpmat,0,Oldcalc.usotao,1,S12,0,num_so,Oldcalc.num_ao,Oldcalc.num_so,0);
  free_block(tmpmat);

  if (print_lvl >= DEBUGPRINT) {
    fprintf(outfile,"  -Overlap matrix between the new and old basis (in SO basis):\n");
    print_mat(S12,num_so,Oldcalc.num_so,outfile);
  }
  
  return S12;
}


/*--------------------------------------------------
  This function computes the overlap matrix S11 for
  the new basis set (in SO basis)
 --------------------------------------------------*/
double **overlap()
{
  struct coordinates A, B, P, PA, PB;

  int ii, jj, kk, ll;
  int maxsz1, maxsz2;
  int x1, x2, y1, y2, z1, z2;
  int s1, s2, ao1, ao2, l1, l2, n1, n2, p1, p2,
      prim1, prim2, fprim1, fprim2, nprim1, nprim2, atom1, atom2, bf1, bf2;
  double exp1, exp2, cc1, cc2, norm1;
  double ab2, pfac, gam, oog;
  double *ptr1, *ptr2;
  double **stemp, **S11_AO, **S11, **OIX, **OIY, **OIZ;
  double **tmpmat;

  S11_AO = block_matrix(num_ao, num_ao);
  maxsz1 = ioff[max_angmom+1];
  maxsz2 = ioff[max_angmom+1];
  stemp = block_matrix(maxsz1,maxsz2);

  OIX = block_matrix(max_angmom+1,max_angmom+1);
  OIY = block_matrix(max_angmom+1,max_angmom+1);
  OIZ = block_matrix(max_angmom+1,max_angmom+1);
  
  /*--- Loop over shells ---*/
  for(s1=0;s1<num_shells;s1++) {
      ao1 = first_ao_shell[s1];
      nprim1 = nprim_in_shell[s1];
      fprim1 = first_prim_shell[s1];
      bf1 = first_ao_shell[s1];
      l1 = shell_ang_mom[s1];
      n1 = ioff[l1+1];
      atom1 = shell_nucleus[s1];
      A.x = geometry[atom1][0];
      A.y = geometry[atom1][1];
      A.z = geometry[atom1][2];
      for(s2=0;s2<num_shells;s2++) {
	  ao2 = first_ao_shell[s2];
	  nprim2 = nprim_in_shell[s2];
	  fprim2 = first_prim_shell[s2];
	  bf2 = first_ao_shell[s2];
	  l2 = shell_ang_mom[s2];
	  n2 = ioff[l2+1];
	  atom2 = shell_nucleus[s2];
	  B.x = geometry[atom2][0];
	  B.y = geometry[atom2][1];
	  B.z = geometry[atom2][2];
	  
	  ab2  = (A.x-B.x)*(A.x-B.x);
	  ab2 += (A.y-B.y)*(A.y-B.y);
	  ab2 += (A.z-B.z)*(A.z-B.z);

	  /*--- zero the temporary storage for accumulating contractions ---*/
	  memset(stemp[0],0,sizeof(double)*maxsz1*maxsz2);

	  /*--- Loop over primitives here ---*/
	  for(p1=0,prim1=fprim1;p1<nprim1;p1++,prim1++) {
	      exp1 = exponents[prim1];
	      cc1 = contr_coeff[prim1];
	      for(p2=0,prim2=fprim2;p2<nprim2;p2++,prim2++) {
		  exp2 = exponents[prim2];
		  cc2 = contr_coeff[prim2];

		  gam = exp1 + exp2;
		  oog = 1.0/gam;
		  P.x = (exp1*A.x + exp2*B.x)*oog;
		  P.y = (exp1*A.y + exp2*B.y)*oog;
		  P.z = (exp1*A.z + exp2*B.z)*oog;
		  PA.x = P.x - A.x;
		  PA.y = P.y - A.y;
		  PA.z = P.z - A.z;
		  PB.x = P.x - B.x;
		  PB.y = P.y - B.y;
		  PB.z = P.z - B.z;

		  pfac = exp(-exp1*exp2*ab2*oog)*sqrt(M_PI*oog)*M_PI*oog*cc1*cc2;

		  OI_OSrecurs(OIX,OIY,OIZ,PA,PB,gam,l1,l2);

		  /*--- create all am components of si ---*/
		  ao1 = 0;
		  for(ii = 0; ii <= l1; ii++){
		      x1 = l1 - ii;
		      for(jj = 0; jj <= ii; jj++){
			  y1 = ii - jj;
			  z1 = jj;
			  /*--- create all am components of sj ---*/
			  ao2 = 0;
			  for(kk = 0; kk <= l2; kk++){
			      x2 = l2 - kk;
			      for(ll = 0; ll <= kk; ll++){
				  y2 = kk - ll;
				  z2 = ll;
				  
				  stemp[ao1][ao2] += pfac*OIX[x1][x2]*OIY[y1][y2]*OIZ[z1][z2];
				  
				  ao2++;
			      }
			  }
			  ao1++;
		      }
		  } /*--- end cartesian components for (si,sj) with primitives (i,j) ---*/
	      }
	  } /*--- end primitive contraction ---*/

	  /*--- Normalize the contracted integrals and put them into S11 ---*/
	  ptr1 = GTOs.bf_norm[l1];
	  ptr2 = GTOs.bf_norm[l2];
	  for(ii=0; ii<n1; ii++) {
	      norm1 = ptr1[ii];
	      for(jj=0; jj<n2; jj++) {
		  S11_AO[bf1+ii][bf2+jj] = stemp[ii][jj]*norm1*ptr2[jj];
	      }
	  }
      }
  }   /*--- This shell pair is done ---*/

  free_block(OIX);
  free_block(OIY);
  free_block(OIZ);
  free_block(stemp);

  /*--- transform to SO basis ---*/
  tmpmat = block_matrix(num_so,num_ao);
  mmult(usotao,0,S11_AO,0,tmpmat,0,num_so,num_ao,num_ao,0);
  free_block(S11_AO);
  S11 = block_matrix(num_so,num_so);
  mmult(tmpmat,0,usotao,1,S11,0,num_so,num_ao,num_so,0);
  free_block(tmpmat);

  if (print_lvl >= DEBUGPRINT) {
      fprintf(outfile,"  -Overlap matrix in the new basis (in SO basis):\n");
      print_mat(S11,num_so,num_so,outfile);
  }
  
  return S11;
}

}} // namespace psi::input
