/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "includes.h"
#include "globals.h"
#include "prototypes.h"

namespace psi { namespace oeprop {

void compute_grid_oeprops()
{
  int i,j,k,l,ig,jg,ibf,jbf,ib,jb,jlim,kk,ll;
  int iatom, jatom, iang, jang, i_1stbf, j_1stbf, nbfi, nbfj;
  int igmin, igmax, jgmin, jgmax;
  int ixm, iym, izm, jxm, jym, jzm, iind, jind;
  int eq_shell, eq_atoms, indmax;
  int ix, iy, iz;
  double ax, ay, az, bx, by, bz, xab, yab, zab, rab2;
  double ai, bj, gamma, over_pf, norm_pf, prod_pf, dens_pf, sdens_pf, zvec_pf;
  double px, py, pz, pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz;
  double cax, cax2, cay, cay2, caz, caz2, cbx, cbx2, cby, cby2, cbz, cbz2, rca2, rcb2;
  double cax_l, cax_l_1, cax_l_2, cay_l, cay_l_1, cay_l_2, caz_l, caz_l_1, caz_l_2;
  double cbx_l, cbx_l_1, cbx_l_2, cby_l, cby_l_1, cby_l_2, cbz_l, cbz_l_1, cbz_l_2;
  int lx1, ly1, lz1, lx2, ly2, lz2;
  double x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;
  double tx, ty, tz, t2x, t2y, t2z;
  double x,y,z;
  double r,r2;
  double tmp, sum1, sum2, sum3;
  double *Density;
  double **temp_mat;

  
		/* Initialize some intermediates */
  grid_pts = init_matrix(nix+1,niy+1);
  
  /* A comment on calculation of electrostatic potentials, 
     electric field and gradient integrals - since these ints are 
     not separable into products of 3 integrals we have to pack 3 indices n,l, and 
     m (powers of x,y, and z in preexponential of the primitive) into one.
     It's done by transforming the number written as nlm in the system with base 
     lmax+1 into decimal number iind (or jind). Maximum value of these numbers is 
     lmax*(lmax+1)^2, therefore AI$ matrices are boxes 
     (lmax*(lmax+1) + 1) x (lmax*(lmax+1) + 1) x (2*lmax+1); */

  indmax = lmax*(lmax+1)*(lmax+1)+1;
  if (grid == 1) {
    AI0 = init_box(indmax,indmax,2*lmax+3);
    AIX = init_box(indmax,indmax,2*lmax+2);
    AIY = init_box(indmax,indmax,2*lmax+2);
    AIZ = init_box(indmax,indmax,2*lmax+2);
    AIXX = init_box(indmax,indmax,2*lmax+1);
    AIYY = init_box(indmax,indmax,2*lmax+1);
    AIZZ = init_box(indmax,indmax,2*lmax+1);
    AIXY = init_box(indmax,indmax,2*lmax+1);
    AIXZ = init_box(indmax,indmax,2*lmax+1);
    AIYZ = init_box(indmax,indmax,2*lmax+1);
  }
  
	/* I-shell loop */

  for(i=0;i<nshell;i++) {
    iatom = snuc[i] - 1;
    ax = geom[iatom][0];
    ay = geom[iatom][1];
    az = geom[iatom][2];
    iang = stype[i]-1;
    izm = 1;
    iym = iang+1;
    ixm = iym*iym;
    i_1stbf = sloc[i];
    nbfi = (iang+2)*(iang+1)/2;
    igmin = sprim[i] - 1;
    igmax = igmin + snumg[i] - 1;
    
    
    	/* J-shell loop */
    
    for(j=0;j<=i;j++) {
      jatom = snuc[j] - 1;
      bx = geom[jatom][0];
      by = geom[jatom][1];
      bz = geom[jatom][2];
      jang = stype[j]-1;
      jzm = 1;
      jym = jang+1;
      jxm = jym*jym;
      j_1stbf = sloc[j];
      nbfj = (jang+2)*(jang+1)/2;
      jgmin = sprim[j] - 1;
      jgmax = jgmin + snumg[j] - 1;
      
      eq_shell = (i == j);
      eq_atoms = (iatom == jatom);
      
      xab = ax - bx;
      yab = ay - by;
      zab = az - bz;
      rab2 = xab*xab + yab*yab + zab*zab;
      
      for(ig=igmin;ig<=igmax;ig++) {
        ai = exps[ig];
        if (eq_shell)
          jgmax = ig;
        for(jg=jgmin;jg<=jgmax;jg++) {
          bj = exps[jg];
          gamma = ai + bj;
          tmp = ai*bj*rab2/gamma;
          
          over_pf = contr[ig]*contr[jg]*exp(-tmp)*pow(M_PI/gamma,1.5);
          if (eq_shell && ig != jg)
            over_pf *= 2.0;
          
          /* Compute P and intermediates */
          
          px = (ax*ai + bx*bj)/gamma;
          py = (ay*ai + by*bj)/gamma;
          pz = (az*ai + bz*bj)/gamma;
          if (eq_atoms)
            pax = pay = paz = pbx = pby = pbz = 0.0;
          else {
            pax = px - ax;  pay = py - ay;  paz = pz - az;
            pbx = px - bx;  pby = py - by;  pbz = pz - bz;
          }

	  
	  /*****************************************************
	  *      Computing various properties over a grid      *
	  *****************************************************/
	  
   	    switch(grid) {

	        /* Electrostatic potential */
	      case 1:
	      for(ix=0;ix<=nix;ix++) {
		pcx = px - (grid_origin[0] + grid_step_x[0]*ix) + grid_step_y[0];
		pcy = py - (grid_origin[1] + grid_step_x[1]*ix) + grid_step_y[1];
		pcz = pz - (grid_origin[2] + grid_step_x[2]*ix) + grid_step_y[2];
	      for(iy=0;iy<=niy;iy++) {
		pcx -= grid_step_y[0];
		pcy -= grid_step_y[1];
                pcz -= grid_step_y[2];
                AI_OSrecurs(pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz,gamma,iang,jang);
                for(ibf=0;ibf<nbfi;ibf++) {
                  if (eq_shell)
                   jlim = ibf + 1;
                  else
                    jlim = nbfj;
                  for(jbf=0;jbf<jlim;jbf++) {
                    ib = i_1stbf + ibf - 1; jb = j_1stbf + jbf - 1;
		    norm_pf = norm_bf[iang][ibf]*norm_bf[jang][jbf];
		    /*---
		      
		      ATTENTION: Normalization prefactor is included in dens_pf
		      
		     ---*/
                    dens_pf = Ptot[ioff[ib]+jb] * ((ib == jb) ? 1.0 : 2.0) * norm_pf;
                    if (fabs(dens_pf) < PFACCUTOFF) continue;
                    iind = zpow_bf[iang][ibf]*izm + ypow_bf[iang][ibf]*iym +
                           xpow_bf[iang][ibf]*ixm;
                    jind = zpow_bf[jang][jbf]*jzm + ypow_bf[jang][jbf]*jym +
                           xpow_bf[jang][jbf]*jxm;
                    
                    grid_pts[ix][iy] += -over_pf*AI0[iind][jind][0] * dens_pf;
                  }
                }
	      }
	      }
	      break;

	    }
	  /* End of loops over primitives (ig,jg)*/
        }
      }
      /* End of loops over shells (i,j)*/
    }
  }


	/* Computing nuclear contribution to the electrostatic potential over a grid */

  if (grid == 1)
    for(ix=0;ix<=nix;ix++) {
      x = grid_origin[0] + grid_step_x[0]*ix - grid_step_y[0]*iy;
      y = grid_origin[1] + grid_step_x[1]*ix - grid_step_y[1]*iy;
      z = grid_origin[2] + grid_step_x[2]*ix - grid_step_y[2]*iy;
      for(iy=0;iy<=niy;iy++) {
        x += grid_step_y[0]*iy;
        y += grid_step_y[1]*iy;
        z += grid_step_y[2]*iy;
        for(k=0;k<natom;k++) {
          r = sqrt((x-geom[k][0])*(x-geom[k][0]) + 
                   (y-geom[k][1])*(y-geom[k][1]) +
                   (z-geom[k][2])*(z-geom[k][2]));
          if (r > 1.0E-9)
            grid_pts[ix][iy] += zvals[k]/r;
        }
      }
    }

  

	/* Cleaning up */
  if (grid == 1) {
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
  }    

}

}} // namespace psi::oeprop
