/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cmath>
#include<cstdio>
#include<cstdlib>
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
#include"oe_osrr.h"

namespace psi {
  namespace cints {

/*!---------------------------------------------------------------------------
  This function computes AO overlap, dipole and quadrupole moment integrals,
  and nabla integrals and writes them out to disk
 ---------------------------------------------------------------------------*/
void moment_ints()
{
   struct coordinates PA, PB, AB, A, B;
   struct shell_pair *sp;
   struct unique_shell_pair *usp;
   register int i, j, k, l, ii, jj, kk, ll;
   int count;
   int si, sj;
   int np_i, np_j;
   int sz;
   int l1, l2, m1, m2, n1, n2;
   int ioffset, joffset ;
   int ij;
   int h1;
   int am;
   int dimension ;
   int ni,nj,ai,aj;
   int am_i, am_j;
   double a1, a2;
   double ab2, oog, gam;
   double x00, y00, z00;
   double x10, y10, z10;
   double x01, y01, z01;
   double x11, y11, z11;
   double nx, ny, nz;
   double **stemp, **mxtemp, **mytemp, **mztemp, **qxxtemp, **qxytemp, **qxztemp,
     **qyytemp, **qyztemp, **qzztemp, **nxtemp, **nytemp, **nztemp,
     **temp, **tmp_dptr;
   double *S, *MX, *MY, *MZ, *QXX, *QXY, *QXZ, *QYY, *QYZ, *QZZ, *NX, *NY, *NZ;
   double inorm, jnorm, over_pf;
   double *ptr1, *ptr2, norm1, norm12;
   double **OIX, **OIY, **OIZ;
   struct coordinates C;


  /*--- allocate room for the one-e matrices ---*/
  dimension = ioff[BasisSet.num_ao];
  S = init_array(dimension);
  MX = init_array(dimension);
  MY = init_array(dimension);
  MZ = init_array(dimension);
  QXX = init_array(dimension);
  QXY = init_array(dimension);
  QXZ = init_array(dimension);
  QYY = init_array(dimension);
  QYZ = init_array(dimension);
  QZZ = init_array(dimension);
  NX = init_array(dimension);
  NY = init_array(dimension);
  NZ = init_array(dimension);

  /*--- allocate storage for shell blocks of one electron integrals ---*/
  dimension = ioff[BasisSet.max_am];
  stemp = block_matrix(dimension,dimension);
  mxtemp = block_matrix(dimension,dimension);
  mytemp = block_matrix(dimension,dimension);
  mztemp = block_matrix(dimension,dimension);
  qxxtemp = block_matrix(dimension,dimension);
  qxytemp = block_matrix(dimension,dimension);
  qxztemp = block_matrix(dimension,dimension);
  qyytemp = block_matrix(dimension,dimension);
  qyztemp = block_matrix(dimension,dimension);
  qzztemp = block_matrix(dimension,dimension);
  nxtemp = block_matrix(dimension,dimension);
  nytemp = block_matrix(dimension,dimension);
  nztemp = block_matrix(dimension,dimension);

  /* Note the "+1" -- we are computing dipole and quadrupole moment integrals here */
  OIX = block_matrix(BasisSet.max_am+1,BasisSet.max_am+1);
  OIY = block_matrix(BasisSet.max_am+1,BasisSet.max_am+1);
  OIZ = block_matrix(BasisSet.max_am+1,BasisSet.max_am+1);
  
  C = UserOptions.origin;
  if(UserOptions.print_lvl >= PRINT_OEI) 
    fprintf(outfile, "    Reference point for elec. mom. ints. = (%5.3f, %5.3f, %5.3f)\n", C.x, C.y, C.z);

  for (si=0; si<BasisSet.num_shells; si++){
    am_i = BasisSet.shells[si].am-1;
    ni = ioff[BasisSet.shells[si].am];
    A = Molecule.centers[BasisSet.shells[si].center-1];
    ioffset = BasisSet.shells[si].fao - 1;
    for (sj=0; sj<=si; sj++){
      nj = ioff[BasisSet.shells[sj].am];
      am_j = BasisSet.shells[sj].am-1;
      B = Molecule.centers[BasisSet.shells[sj].center-1];
      joffset = BasisSet.shells[sj].fao - 1;

      sp = &(BasisSet.shell_pairs[si][sj]);
      AB.x = sp->AB[0];
      AB.y = sp->AB[1];
      AB.z = sp->AB[2];
      ab2 = AB.x * AB.x;
      ab2 += AB.y * AB.y;
      ab2 += AB.z * AB.z;
	
      /*--- zero the temporary storage for accumulating contractions ---*/
      memset(stemp[0],0,sizeof(double)*dimension*dimension);
      memset(mxtemp[0],0,sizeof(double)*dimension*dimension);
      memset(mytemp[0],0,sizeof(double)*dimension*dimension);
      memset(mztemp[0],0,sizeof(double)*dimension*dimension);
      memset(qxxtemp[0],0,sizeof(double)*dimension*dimension);
      memset(qxytemp[0],0,sizeof(double)*dimension*dimension);
      memset(qxztemp[0],0,sizeof(double)*dimension*dimension);
      memset(qyytemp[0],0,sizeof(double)*dimension*dimension);
      memset(qyztemp[0],0,sizeof(double)*dimension*dimension);
      memset(qzztemp[0],0,sizeof(double)*dimension*dimension);
      memset(nxtemp[0],0,sizeof(double)*dimension*dimension);
      memset(nytemp[0],0,sizeof(double)*dimension*dimension);
      memset(nztemp[0],0,sizeof(double)*dimension*dimension);
      
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

	  OI_OSrecurs(OIX,OIY,OIZ,PA,PB,gam,am_i+1,am_j+1);

	  /*--- create all am components of si ---*/
	  ai = 0;
	  for(ii = 0; ii <= am_i; ii++){
	    l1 = am_i - ii;
	    for(jj = 0; jj <= ii; jj++){
	      m1 = ii - jj;
	      n1 = jj;
	      /*--- create all am components of sj ---*/
	      aj = 0;
	      for(kk = 0; kk <= am_j; kk++){
		l2 = am_j - kk;
		for(ll = 0; ll <= kk; ll++){
		  m2 = kk - ll;
		  n2 = ll;

		  x00 = OIX[l1][l2];   y00 = OIY[m1][m2];    z00 = OIZ[n1][n2];
		  x01 = OIX[l1][l2+1]; y01 = OIY[m1][m2+1];  z01 = OIZ[n1][n2+1];
		  x10 = OIX[l1+1][l2]; y10 = OIY[m1+1][m2];  z10 = OIZ[n1+1][n2];
		  x11 = OIX[l1+1][l2+1]; y11 = OIY[m1+1][m2+1];  z11 = OIZ[n1+1][n2+1];
		  stemp[ai][aj] += over_pf*x00*y00*z00;
		  /*--- electrons have negative charge ---*/
		  mxtemp[ai][aj] -= over_pf*(x01+x00*(B.x-C.x))*y00*z00;
		  mytemp[ai][aj] -= over_pf*x00*(y01+y00*(B.y-C.y))*z00;
		  mztemp[ai][aj] -= over_pf*x00*y00*(z01+z00*(B.z-C.z));
		  qxxtemp[ai][aj] -= over_pf*(x11 + x10*(B.x-C.x) + x01*(A.x-C.x)
                                              + x00*(A.x-C.x)*(B.x-C.x))*y00*z00;
		  qyytemp[ai][aj] -= over_pf*(y11 + y10*(B.y-C.y) + y01*(A.y-C.y)
                                              + y00*(A.y-C.y)*(B.y-C.y))*x00*z00;
		  qzztemp[ai][aj] -= over_pf*(z11 + z10*(B.z-C.z) + z01*(A.z-C.z)
                                              + z00*(A.z-C.z)*(B.z-C.z))*x00*y00;
		  qxytemp[ai][aj] -= over_pf*(x01+x00*(B.x-C.x))*(y01+y00*(B.y-C.y))*z00;
		  qxztemp[ai][aj] -= over_pf*(x01+x00*(B.x-C.x))*y00*(z01+z00*(B.z-C.z));
		  qyztemp[ai][aj] -= over_pf*x00*(y01+y00*(B.y-C.y))*(z01+z00*(B.z-C.z));

		  nx = -2.0*a2*x01;
		  if (l2 >= 1)
		    nx += l2*OIX[l1][l2-1];
		  nxtemp[ai][aj] += nx*y00*z00*over_pf;
		  ny = -2.0*a2*y01;
		  if (m2 >= 1)
		    ny += m2*OIY[m1][m2-1];
		  nytemp[ai][aj] += x00*ny*z00*over_pf;
		  nz = -2.0*a2*z01;
		  if (n2 >= 1)
		    nz += n2*OIZ[n1][n2-1];
		  nztemp[ai][aj] += x00*y00*nz*over_pf;
		  
		  aj++;
		}
	      }
	      ai++;
	    }
	  } /*--- end cartesian components for (si,sj) with primitives (i,j) ---*/
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
	  mxtemp[i][j] *= norm12;
	  mytemp[i][j] *= norm12;
	  mztemp[i][j] *= norm12;
	  qxxtemp[i][j] *= norm12;
	  qxytemp[i][j] *= norm12;
	  qxztemp[i][j] *= norm12;
	  qyytemp[i][j] *= norm12;
	  qyztemp[i][j] *= norm12;
	  qzztemp[i][j] *= norm12;
	  nxtemp[i][j] *= norm12;
	  nytemp[i][j] *= norm12;
	  nztemp[i][j] *= norm12;
	}
      }

      for(i=0;i<ni;i++)
	for(j=0;j<nj;j++) {
	  ij = INDEX(ioffset + i, joffset + j);
	  S[ij] = stemp[i][j];
	  MX[ij] = mxtemp[i][j];
	  MY[ij] = mytemp[i][j];
	  MZ[ij] = mztemp[i][j];
	  QXX[ij] = qxxtemp[i][j];
	  QXY[ij] = qxytemp[i][j];
	  QXZ[ij] = qxztemp[i][j];
	  QYY[ij] = qyytemp[i][j];
	  QYZ[ij] = qyztemp[i][j];
	  QZZ[ij] = qzztemp[i][j];
	  NX[ij] = nxtemp[i][j];
	  NY[ij] = nytemp[i][j];
	  NZ[ij] = nztemp[i][j];
	}
    }
  }  /*--- This shell pair is done ---*/

  /*--- flush it all away ---*/
  fflush(outfile);
  free_block(OIX);
  free_block(OIY);
  free_block(OIZ);
  dimension=ioff[BasisSet.num_ao];
  iwl_wrtone(IOUnits.itapS_AO, PSIF_AO_S, dimension,S);
  iwl_wrtone(IOUnits.itapMX_AO,PSIF_AO_MX,dimension,MX);
  iwl_wrtone(IOUnits.itapMY_AO,PSIF_AO_MY,dimension,MY);
  iwl_wrtone(IOUnits.itapMZ_AO,PSIF_AO_MZ,dimension,MZ);
  iwl_wrtone(IOUnits.itapQXX_AO,PSIF_AO_QXX,dimension,QXX);
  iwl_wrtone(IOUnits.itapQXY_AO,PSIF_AO_QXY,dimension,QXY);
  iwl_wrtone(IOUnits.itapQXZ_AO,PSIF_AO_QXZ,dimension,QXZ);
  iwl_wrtone(IOUnits.itapQYY_AO,PSIF_AO_QYY,dimension,QYY);
  iwl_wrtone(IOUnits.itapQYZ_AO,PSIF_AO_QYZ,dimension,QYZ);
  iwl_wrtone(IOUnits.itapQZZ_AO,PSIF_AO_QZZ,dimension,QZZ);
  iwl_wrtone(IOUnits.itapNablaX_AO,PSIF_AO_NablaX,dimension,NX);
  iwl_wrtone(IOUnits.itapNablaY_AO,PSIF_AO_NablaY,dimension,NY);
  iwl_wrtone(IOUnits.itapNablaZ_AO,PSIF_AO_NablaZ,dimension,NZ);
  if (UserOptions.print_lvl >= PRINT_OEI) {
    fprintf(outfile,"  -Overlap AO integrals:\n\n");
    print_array(S,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -mu(x) AO integrals:\n\n");
    print_array(MX,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -mu(y) AO integrals:\n\n");
    print_array(MY,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -mu(z) AO integrals:\n\n");
    print_array(MZ,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -q(xx) AO integrals:\n\n");
    print_array(QXX,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -q(xy) AO integrals:\n\n");
    print_array(QXY,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -q(xz) AO integrals:\n\n");
    print_array(QXZ,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -q(yy) AO integrals:\n\n");
    print_array(QYY,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -q(yz) AO integrals:\n\n");
    print_array(QYZ,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -q(zz) AO integrals:\n\n");
    print_array(QZZ,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -Nabla_x AO integrals:\n\n");
    print_array(NX,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -Nabla_y AO integrals:\n\n");
    print_array(NY,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -Nabla_z AO integrals:\n\n");
    print_array(NZ,BasisSet.num_ao,outfile);
    fprintf(outfile,"\n");
  }

  free_block(stemp);
  free_block(mxtemp);
  free_block(mytemp);
  free_block(mztemp);
  free_block(qxxtemp);
  free_block(qxytemp);
  free_block(qxztemp);
  free_block(qyytemp);
  free_block(qyztemp);
  free_block(qzztemp);
  free_block(nxtemp);
  free_block(nytemp);
  free_block(nztemp);

  free(S);
  free(MX);
  free(MY);
  free(MZ);
  free(QXX);
  free(QXY);
  free(QXZ);
  free(QYY);
  free(QYZ);
  free(QZZ);
  free(NX);
  free(NY);
  free(NZ);

  return;
}   
  }
}
