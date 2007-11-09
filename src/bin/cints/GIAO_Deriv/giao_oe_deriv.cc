/*! \file giao_oe_deriv.cc
    \ingroup (CINTS)
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
#include"taylor_fm_eval.h"
#include"oe_osrr.h"
#include"small_fns.h"


namespace psi { namespace CINTS {

/*!-----------------------------------------------------------------------------------------
  This function computes some derivatives of integrals over GIAO Gaussians with respect to
  E and B fields (at E=0, B=0, i.e. in absence of external fields) and writes them out.
  
                                   t
  d2h/dBdE = -0.5 i <mu| (AB x r) r  |nu>
  
  d h/dB   =  0.5   <mu| L   + i (AB x r) h |nu>
                          B

  d h/dE   =  - <mu| r |nu>   -- simply electric dipole integral, already computed, hence
                                 not computed here
 -----------------------------------------------------------------------------------------*/
void giao_oe_deriv()
{

#ifdef USE_TAYLOR_FM
  init_Taylor_Fm_Eval(BasisSet.max_am*4+1,UserOptions.cutoff);
#endif  

  /*--- allocate room for the one-e matrices ---*/
  double** D2HDBXDEX = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  double** D2HDBXDEY = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  double** D2HDBXDEZ = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  double** D2HDBYDEX = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  double** D2HDBYDEY = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  double** D2HDBYDEZ = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  double** D2HDBZDEX = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  double** D2HDBZDEY = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  double** D2HDBZDEZ = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  double** DHDBX = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  double** DHDBY = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  double** DHDBZ = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  double** DSDBX = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  double** DSDBY = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  double** DSDBZ = block_matrix(BasisSet.num_ao,BasisSet.num_ao);

  /*--- allocate storage for shell blocks of one electron integrals ---*/
  int dimension = ioff[BasisSet.max_am];
  double** d2HdBxdExtemp = block_matrix(dimension,dimension);
  double** d2HdBxdEytemp = block_matrix(dimension,dimension);
  double** d2HdBxdEztemp = block_matrix(dimension,dimension);
  double** d2HdBydExtemp = block_matrix(dimension,dimension);
  double** d2HdBydEytemp = block_matrix(dimension,dimension);
  double** d2HdBydEztemp = block_matrix(dimension,dimension);
  double** d2HdBzdExtemp = block_matrix(dimension,dimension);
  double** d2HdBzdEytemp = block_matrix(dimension,dimension);
  double** d2HdBzdEztemp = block_matrix(dimension,dimension);
  double** dHdBxtemp = block_matrix(dimension,dimension);
  double** dHdBytemp = block_matrix(dimension,dimension);
  double** dHdBztemp = block_matrix(dimension,dimension);
  double** dSdBxtemp = block_matrix(dimension,dimension);
  double** dSdBytemp = block_matrix(dimension,dimension);
  double** dSdBztemp = block_matrix(dimension,dimension);

  /* Note the "+1" on the row dimension -- we are computing kinetic energy
     integrals with one extra quantum in the bra function here */
  double** OIX = block_matrix(BasisSet.max_am+1,BasisSet.max_am+2);
  double** OIY = block_matrix(BasisSet.max_am+1,BasisSet.max_am+2);
  double** OIZ = block_matrix(BasisSet.max_am+1,BasisSet.max_am+2);
  /* same here -- need to compute integrals with extra quantum in the bra */
  int indmax_0 = (BasisSet.max_am-1)*BasisSet.max_am*BasisSet.max_am+1;
  int indmax_p1 = (BasisSet.max_am)*(BasisSet.max_am+1)*(BasisSet.max_am+1)+1;
  double*** AI0 = init_box(indmax_p1,indmax_0,2*BasisSet.max_am+3);
  
  for (int si=0; si<BasisSet.num_shells; si++){
    int am_i = BasisSet.shells[si].am-1;
    int ni = ioff[BasisSet.shells[si].am];
    struct coordinates A = Molecule.centers[BasisSet.shells[si].center-1];
    int ioffset = BasisSet.shells[si].fao - 1;
    int izm = 1;
    int iym = am_i+2;     // Note +2 instead of +1 here -- we'll be computing nuclear attraction
                          // integrals with extra quantum in bra!
    int ixm = iym*iym;

    for (int sj=0; sj<BasisSet.num_shells; sj++){
      int nj = ioff[BasisSet.shells[sj].am];
      int am_j = BasisSet.shells[sj].am-1;
      struct coordinates B = Molecule.centers[BasisSet.shells[sj].center-1];
      int joffset = BasisSet.shells[sj].fao - 1;
      int jzm = 1;
      int jym = am_j+1;
      int jxm = jym*jym;

      struct shell_pair* sp = &(BasisSet.shell_pairs[si][sj]);
      struct coordinates AB;
      AB.x = sp->AB[0];
      AB.y = sp->AB[1];
      AB.z = sp->AB[2];
      double ab2 = AB.x * AB.x;
      ab2 += AB.y * AB.y;
      ab2 += AB.z * AB.z;

      /*--- zero the temporary storage for accumulating contractions ---*/
      memset(d2HdBxdExtemp[0],0,sizeof(double)*dimension*dimension);
      memset(d2HdBxdEytemp[0],0,sizeof(double)*dimension*dimension);
      memset(d2HdBxdEztemp[0],0,sizeof(double)*dimension*dimension);
      memset(d2HdBydExtemp[0],0,sizeof(double)*dimension*dimension);
      memset(d2HdBydEytemp[0],0,sizeof(double)*dimension*dimension);
      memset(d2HdBydEztemp[0],0,sizeof(double)*dimension*dimension);
      memset(d2HdBzdExtemp[0],0,sizeof(double)*dimension*dimension);
      memset(d2HdBzdEytemp[0],0,sizeof(double)*dimension*dimension);
      memset(d2HdBzdEztemp[0],0,sizeof(double)*dimension*dimension);
      memset(dHdBxtemp[0],0,sizeof(double)*dimension*dimension);
      memset(dHdBytemp[0],0,sizeof(double)*dimension*dimension);
      memset(dHdBztemp[0],0,sizeof(double)*dimension*dimension);
      memset(dSdBxtemp[0],0,sizeof(double)*dimension*dimension);
      memset(dSdBytemp[0],0,sizeof(double)*dimension*dimension);
      memset(dSdBztemp[0],0,sizeof(double)*dimension*dimension);
      
      /*--- contract by primitives here ---*/
      for (int i = 0; i < BasisSet.shells[si].n_prims; i++) {
	double a1 = sp->a1[i];
	double inorm = sp->inorm[i];
	for (int j = 0; j < BasisSet.shells[sj].n_prims; j++) {
	  double a2 = sp->a2[j];
	  double gam = sp->gamma[i][j];
	  double jnorm = sp->jnorm[j];
          struct coordinates PA, PB;
	  PA.x = sp->PA[i][j][0];
	  PA.y = sp->PA[i][j][1];
	  PA.z = sp->PA[i][j][2];
	  PB.x = sp->PB[i][j][0];
	  PB.y = sp->PB[i][j][1];
	  PB.z = sp->PB[i][j][2];
	  double oog = 1.0/gam;
	  double over_pf = exp(-a1*a2*ab2*oog)*sqrt(M_PI*oog)*M_PI*oog*inorm*jnorm;

	  OI_OSrecurs(OIX,OIY,OIZ,PA,PB,gam,am_i+1,am_j+2);

	  /*--- create all am components of si ---*/
	  int ai = 0;
	  for(int ii = 0; ii <= am_i; ii++){
	    int l1 = am_i - ii;
	    for(int jj = 0; jj <= ii; jj++){
	      int m1 = ii - jj;
	      int n1 = jj;
	      /*--- create all am components of sj ---*/
	      int aj = 0;
	      for(int kk = 0; kk <= am_j; kk++){
		int l2 = am_j - kk;
		for(int ll = 0; ll <= kk; ll++){
		  int m2 = kk - ll;
		  int n2 = ll;

		  double x00 = OIX[l1][l2];
                  double y00 = OIY[m1][m2];
                  double z00 = OIZ[n1][n2];
                  
		  double x01 = OIX[l1][l2+1];
                  double y01 = OIY[m1][m2+1];
                  double z01 = OIZ[n1][n2+1];
                  
		  double x10 = OIX[l1+1][l2];
                  double y10 = OIY[m1+1][m2];
                  double z10 = OIZ[n1+1][n2];
                  
		  double x11 = OIX[l1+1][l2+1];
                  double y11 = OIY[m1+1][m2+1];
                  double z11 = OIZ[n1+1][n2+1];

                  double x0 = x00;
                  double y0 = y00;
                  double z0 = z00;
                  double x1 = x01 + x00*B.x;
                  double y1 = y01 + y00*B.y;
                  double z1 = z01 + z00*B.z;
                  double x2 = x11 + x10*B.x + x01*A.x + x00*A.x*B.x;
                  double y2 = y11 + y10*B.y + y01*A.y + y00*A.y*B.y;
                  double z2 = z11 + z10*B.z + z01*A.z + z00*A.z*B.z;

                  double dSdB_pfac = 0.5*over_pf;
                  double x = x1*y0*z0;
                  double y = x0*y1*z0;
                  double z = x0*y0*z1;
                  dSdBxtemp[ai][aj] += dSdB_pfac * (AB.y * z - AB.z * y);
                  dSdBytemp[ai][aj] += dSdB_pfac * (AB.z * x - AB.x * z);
                  dSdBztemp[ai][aj] += dSdB_pfac * (AB.x * y - AB.y * x);

/*                  dSdBxtemp[ai][aj] += 2.0*dSdB_pfac * x;
                  dSdBytemp[ai][aj] += 2.0*dSdB_pfac * y;
                  dSdBztemp[ai][aj] += 2.0*dSdB_pfac * z;*/
                  
                  double d2BE_pfac = -0.5*over_pf;
                  double qxx = x2*y0*z0;
                  double qxy = x1*y1*z0;
                  double qxz = x1*y0*z1;
                  double qyx = qxy;
                  double qyy = x0*y2*z0;
                  double qyz = x0*y1*z1;
                  double qzx = qxz;
                  double qzy = qyz;
                  double qzz = x0*y0*z2;
                  d2HdBxdExtemp[ai][aj] += d2BE_pfac * ( AB.y * qzx - AB.z * qyx);
                  d2HdBxdEytemp[ai][aj] += d2BE_pfac * ( AB.y * qzy - AB.z * qyy);
                  d2HdBxdEztemp[ai][aj] += d2BE_pfac * ( AB.y * qzz - AB.z * qyz);
                  d2HdBydExtemp[ai][aj] += d2BE_pfac * ( AB.z * qxx - AB.x * qzx);
                  d2HdBydEytemp[ai][aj] += d2BE_pfac * ( AB.z * qxy - AB.x * qzy);
                  d2HdBydEztemp[ai][aj] += d2BE_pfac * ( AB.z * qxz - AB.x * qzz);
                  d2HdBzdExtemp[ai][aj] += d2BE_pfac * ( AB.x * qyx - AB.y * qxx);
                  d2HdBzdEytemp[ai][aj] += d2BE_pfac * ( AB.x * qyy - AB.y * qxy);
                  d2HdBzdEztemp[ai][aj] += d2BE_pfac * ( AB.x * qyz - AB.y * qxz);
                  
		  double nx = -2.0*a2*x01;
		  if (l2 >= 1)
		    nx += l2*OIX[l1][l2-1];
		  double ny = -2.0*a2*y01;
		  if (m2 >= 1)
		    ny += m2*OIY[m1][m2-1];
		  double nz = -2.0*a2*z01;
		  if (n2 >= 1)
		    nz += n2*OIZ[n1][n2-1];

                  double dB_pfac = -0.5*over_pf;
		  dHdBxtemp[ai][aj] += dB_pfac * x00 * ( y01 * nz - z01 * ny);
		  dHdBytemp[ai][aj] += dB_pfac * y00 * ( z01 * nx - x01 * nz);
		  dHdBztemp[ai][aj] += dB_pfac * z00 * ( x01 * ny - y01 * nx);

                  double tx00 = a2*(2*l2+1)*OIX[l1][l2] - 2*a2*a2*OIX[l1][l2+2];
                  if (l2 >= 2)
                    tx00 -= 0.5*l2*(l2-1)*OIX[l1][l2-2];
                  double ty00 = a2*(2*m2+1)*OIY[m1][m2] - 2*a2*a2*OIY[m1][m2+2];
                  if (m2 >= 2)
                    ty00 -= 0.5*m2*(m2-1)*OIY[m1][m2-2];
                  double tz00 = a2*(2*n2+1)*OIZ[n1][n2] - 2*a2*a2*OIZ[n1][n2+2];
                  if (n2 >= 2)
                    tz00 -= 0.5*n2*(n2-1)*OIZ[n1][n2-2];

                  double tx10 = a2*(2*l2+1)*(OIX[l1+1][l2] + A.x*OIX[l1][l2]) - 2*a2*a2*(OIX[l1+1][l2+2] + A.x*OIX[l1][l2+2]);
                  if (l2 >= 2)
                    tx10 -= 0.5*l2*(l2-1)*(OIX[l1+1][l2-2] + A.x*OIX[l1][l2-2]);
                  double ty10 = a2*(2*m2+1)*(OIY[m1+1][m2] + A.y*OIY[m1][m2]) - 2*a2*a2*(OIY[m1+1][m2+2] + A.y*OIY[m1][m2+2]);
                  if (m2 >= 2)
                    ty10 -= 0.5*m2*(m2-1)*(OIY[m1+1][m2-2] + A.y*OIY[m1][m2-2]);
                  double tz10 = a2*(2*n2+1)*(OIZ[n1+1][n2] + A.z*OIZ[n1][n2]) - 2*a2*a2*(OIZ[n1+1][n2+2] + A.z*OIZ[n1][n2+2]);
                  if (n2 >= 2)
                    tz10 -= 0.5*n2*(n2-1)*(OIZ[n1+1][n2-2] + A.z*OIZ[n1][n2-2]);

                  double xT = tx10*y0*z0 + x1*ty00*z0 + x1*y0*tz00;
                  double yT = tx00*y1*z0 + x0*ty10*z0 + x0*y1*tz00;
                  double zT = tx00*y0*z1 + x0*ty00*z1 + x0*y0*tz10;

                  dB_pfac = 0.5*over_pf;
                  dHdBxtemp[ai][aj] += dB_pfac * ( AB.y * zT - AB.z * yT);
                  dHdBytemp[ai][aj] += dB_pfac * ( AB.z * xT - AB.x * zT);
                  dHdBztemp[ai][aj] += dB_pfac * ( AB.x * yT - AB.y * xT);

/*                  dHdBxtemp[ai][aj] += 2.0*dB_pfac * (tx00*y0*z0 + x0*ty00*z0+x0*y0*tz00);
                  dHdBztemp[ai][aj] += 2.0*dB_pfac * (tx00*y0*z0 + x0*ty00*z0+x0*y0*tz00);*/
                  
                  /*double T = tx00*y0*z0 + x0*ty00*z0 + x0*y0*tz00;
                  dHdBytemp[ai][aj] += 2.0*dB_pfac*T;*/
		  
		  aj++;
		}
	      }
	      ai++;
	    }
	  } /*--- end cartesian components for (si,sj) with primitives (i,j) ---*/
          
          for(int atom=0;atom<Molecule.num_atoms;atom++) {
            struct coordinates PC;
            PC.x = sp->P[i][j][0] - Molecule.centers[atom].x;
            PC.y = sp->P[i][j][1] - Molecule.centers[atom].y;
            PC.z = sp->P[i][j][2] - Molecule.centers[atom].z;
            AI_OSrecurs(AI0,PA,PB,PC,gam,am_i+1,am_j);

            double dB_pfac = 0.5 * over_pf * (-Molecule.centers[atom].Z_nuc);

            int ai = 0;
            for(int ii = 0; ii <= am_i; ii++){
              int l1 = am_i - ii;
              for(int jj = 0; jj <= ii; jj++){
                int m1 = ii - jj;
                int n1 = jj ;
                int iind_000 = n1*izm + m1*iym + l1*ixm;
                int iind_100 = n1*izm + m1*iym + (l1+1)*ixm;
                int iind_010 = n1*izm + (m1+1)*iym + l1*ixm;
                int iind_001 = (n1+1)*izm + m1*iym + l1*ixm;
                /*--- create all am components of sj ---*/
                int aj = 0;
                for(int kk = 0; kk <= am_j; kk++){
                  int l2 = am_j - kk;
                  for(int ll = 0; ll <= kk; ll++){
                    int m2 = kk - ll;
                    int n2 = ll ;

                    int jind = n2*jzm + m2*jym + l2*jxm;

                    double V = AI0[iind_000][jind][0];
                    double xV = AI0[iind_100][jind][0] + A.x * V;
                    double yV = AI0[iind_010][jind][0] + A.y * V;
                    double zV = AI0[iind_001][jind][0] + A.z * V;

                    dHdBxtemp[ai][aj] += dB_pfac * ( AB.y * zV - AB.z * yV);
                    dHdBytemp[ai][aj] += dB_pfac * ( AB.z * xV - AB.x * zV);
                    dHdBztemp[ai][aj] += dB_pfac * ( AB.x * yV - AB.y * xV);

/*                    dHdBxtemp[ai][aj] += dB_pfac * xV;
                    dHdBytemp[ai][aj] += dB_pfac * yV;
                    dHdBztemp[ai][aj] += dB_pfac * zV;*/
                  
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
      double* ptr1 = GTOs.bf_norm[am_i];
      double* ptr2 = GTOs.bf_norm[am_j];
      for(int i=0; i<ni; i++) {
	double norm1 = ptr1[i];
	for(int j=0; j<nj; j++) {
	  double norm12 = norm1*ptr2[j];
	  d2HdBxdExtemp[i][j] *= norm12;
	  d2HdBxdEytemp[i][j] *= norm12;
	  d2HdBxdEztemp[i][j] *= norm12;
	  d2HdBydExtemp[i][j] *= norm12;
	  d2HdBydEytemp[i][j] *= norm12;
	  d2HdBydEztemp[i][j] *= norm12;
	  d2HdBzdExtemp[i][j] *= norm12;
	  d2HdBzdEytemp[i][j] *= norm12;
	  d2HdBzdEztemp[i][j] *= norm12;
	  dHdBxtemp[i][j] *= norm12;
	  dHdBytemp[i][j] *= norm12;
	  dHdBztemp[i][j] *= norm12;
	  dSdBxtemp[i][j] *= norm12;
	  dSdBytemp[i][j] *= norm12;
	  dSdBztemp[i][j] *= norm12;
	}
      }

      for(int i=0;i<ni;i++)
	for(int j=0;j<nj;j++) {
          int ii = ioffset + i;
          int jj = joffset + j;
	  D2HDBXDEX[ii][jj] = d2HdBxdExtemp[i][j];
	  D2HDBXDEY[ii][jj] = d2HdBxdEytemp[i][j];
	  D2HDBXDEZ[ii][jj] = d2HdBxdEztemp[i][j];
	  D2HDBYDEX[ii][jj] = d2HdBydExtemp[i][j];
	  D2HDBYDEY[ii][jj] = d2HdBydEytemp[i][j];
	  D2HDBYDEZ[ii][jj] = d2HdBydEztemp[i][j];
	  D2HDBZDEX[ii][jj] = d2HdBzdExtemp[i][j];
	  D2HDBZDEY[ii][jj] = d2HdBzdEytemp[i][j];
	  D2HDBZDEZ[ii][jj] = d2HdBzdEztemp[i][j];
	  DHDBX[ii][jj] = dHdBxtemp[i][j];
	  DHDBY[ii][jj] = dHdBytemp[i][j];
	  DHDBZ[ii][jj] = dHdBztemp[i][j];
	  DSDBX[ii][jj] = dSdBxtemp[i][j];
	  DSDBY[ii][jj] = dSdBytemp[i][j];
	  DSDBZ[ii][jj] = dSdBztemp[i][j];
	}
    }
  }  /*--- This shell pair is done ---*/

  /*--- flush it all away ---*/
  fflush(outfile);
  free_block(OIX);
  free_block(OIY);
  free_block(OIZ);
  free_box(AI0,indmax_p1,indmax_0);
  dimension=BasisSet.num_ao*BasisSet.num_ao;
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_D2HDBDE_XX, dimension,D2HDBXDEX[0]);
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_D2HDBDE_XY, dimension,D2HDBXDEY[0]);
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_D2HDBDE_XZ, dimension,D2HDBXDEZ[0]);
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_D2HDBDE_YX, dimension,D2HDBYDEX[0]);
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_D2HDBDE_YY, dimension,D2HDBYDEY[0]);
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_D2HDBDE_YZ, dimension,D2HDBYDEZ[0]);
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_D2HDBDE_ZX, dimension,D2HDBZDEX[0]);
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_D2HDBDE_ZY, dimension,D2HDBZDEY[0]);
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_D2HDBDE_ZZ, dimension,D2HDBZDEZ[0]);
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_DHDB_X, dimension,DHDBX[0]);
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_DHDB_Y, dimension,DHDBY[0]);
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_DHDB_Z, dimension,DHDBZ[0]);
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_DSDB_X, dimension,DSDBX[0]);
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_DSDB_Y, dimension,DSDBY[0]);
  iwl_wrtone(IOUnits.itapOEInt_Misc, PSIF_AO_DSDB_Z, dimension,DSDBZ[0]);
  if (UserOptions.print_lvl >= PRINT_OEI) {
    fprintf(outfile,"  -dS/dB_x AO integrals:\n\n");
    print_mat(DSDBX,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -dS/dB_y AO integrals:\n\n");
    print_mat(DSDBY,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -dS/dB_z AO integrals:\n\n");
    print_mat(DSDBZ,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -dh/dB_x AO integrals:\n\n");
    print_mat(DHDBX,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -dh/dB_y AO integrals:\n\n");
    print_mat(DHDBY,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -dh/dB_z AO integrals:\n\n");
    print_mat(DHDBZ,BasisSet.num_ao,BasisSet.num_ao,outfile);
/*    fprintf(outfile,"  -d2h/dB_x dE_x AO integrals:\n\n");
    print_mat(D2HDBXDEX,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -d2h/dB_x dE_y AO integrals:\n\n");
    print_mat(D2HDBXDEY,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -d2h/dB_x dE_z AO integrals:\n\n");
    print_mat(D2HDBXDEZ,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -d2h/dB_y dE_x AO integrals:\n\n");
    print_mat(D2HDBYDEX,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -d2h/dB_y dE_y AO integrals:\n\n");
    print_mat(D2HDBYDEY,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -d2h/dB_y dE_z AO integrals:\n\n");
    print_mat(D2HDBYDEZ,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -d2h/dB_z dE_x AO integrals:\n\n");
    print_mat(D2HDBZDEX,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -d2h/dB_z dE_y AO integrals:\n\n");
    print_mat(D2HDBZDEY,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -d2h/dB_z dE_z AO integrals:\n\n");
    print_mat(D2HDBZDEZ,BasisSet.num_ao,BasisSet.num_ao,outfile);*/
    fprintf(outfile,"\n");
  }

  free_block(d2HdBxdExtemp);
  free_block(d2HdBxdEytemp);
  free_block(d2HdBxdEztemp);
  free_block(d2HdBydExtemp);
  free_block(d2HdBydEytemp);
  free_block(d2HdBydEztemp);
  free_block(d2HdBzdExtemp);
  free_block(d2HdBzdEytemp);
  free_block(d2HdBzdEztemp);
  free_block(dHdBxtemp);
  free_block(dHdBytemp);
  free_block(dHdBztemp);
  free_block(dSdBxtemp);
  free_block(dSdBytemp);
  free_block(dSdBztemp);

  free_block(D2HDBXDEX);
  free_block(D2HDBXDEY);
  free_block(D2HDBXDEZ);
  free_block(D2HDBYDEX);
  free_block(D2HDBYDEY);
  free_block(D2HDBYDEZ);
  free_block(D2HDBZDEX);
  free_block(D2HDBZDEY);
  free_block(D2HDBZDEZ);
  free_block(DHDBX);
  free_block(DHDBY);
  free_block(DHDBZ);
  free_block(DSDBX);
  free_block(DSDBY);
  free_block(DSDBZ);

#ifdef USE_TAYLOR_FM
  free_Taylor_Fm_Eval();
#endif

  return;
}   
};};
