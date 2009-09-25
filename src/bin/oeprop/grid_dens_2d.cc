/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "includes.h"
#include "globals.h"
#include "prototypes.h"

namespace psi { namespace oeprop {

void compute_grid_dens_2d()
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
  if (grid == 3) {
      grid_vecX = init_matrix(nix+1,niy+1);
      grid_vecY = init_matrix(nix+1,niy+1);
      grid_vecZ = init_matrix(nix+1,niy+1);
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

		/* Electron density */
	      case 2:
	      if (spin_prop)    /* If spin_prop is set, compute spin density... */
		Density = Pspin;
	      else              /* .. otherwise compute electron density */
		Density = Ptot;
	      for(ix=0;ix<=nix;ix++) {
		pcx = grid_origin[0] + grid_step_x[0]*ix - grid_step_y[0];
		pcy = grid_origin[1] + grid_step_x[1]*ix - grid_step_y[1];
		pcz = grid_origin[2] + grid_step_x[2]*ix - grid_step_y[2];
	      for(iy=0;iy<=niy;iy++) {
		pcx += grid_step_y[0];
		pcy += grid_step_y[1];
		pcz += grid_step_y[2];
		cax = pcx-ax; cbx = pcx-bx;
		cax2 = cax*cax; cbx2 = cbx*cbx;
		cay = pcy-ay; cby = pcy-by;
		cay2 = cay*cay; cby2 = cby*cby;
		caz = pcz-az; cbz = pcz-bz;
		caz2 = caz*caz; cbz2 = cbz*cbz;
		prod_pf = contr[ig]*contr[jg]*exp(-ai*(cax2 + cay2 + caz2))*exp(-bj*(cbx2 + cby2 + cbz2));
		sum1 = 0.0;
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
                    dens_pf = Density[ioff[ib]+jb] * ((ib == jb) ? 1.0 : 2.0) * norm_pf;
		    if (fabs(dens_pf) < PFACCUTOFF) continue;
		    if (eq_shell && ig != jg)
                      dens_pf *= 2.0;
		    /* For **efficiency** reasons I'm trying to avoid calling pow() for up to
		       g-functions through switch() constructs */
                    lx1 = xpow_bf[iang][ibf];
		    switch (lx1) {
		      case 0:
		      cax_l = 1.0;
		      break;
		      case 1:
		      cax_l = cax;
		      break;
		      case 2:
		      cax_l = cax2;
		      break;
		      case 3:
		      cax_l = cax2*cax;
		      break;
		      case 4:
		      cax_l = cax2*cax2;
		      break;
		      default:
		      cax_l = pow(cax,lx1);
		    }
		    ly1 = ypow_bf[iang][ibf];
		    switch (ly1) {
		      case 0:
		      cay_l = 1.0;
		      break;
		      case 1:
		      cay_l = cay;
		      break;
		      case 2:
		      cay_l = cay2;
		      break;
		      case 3:
		      cay_l = cay2*cay;
		      break;
		      case 4:
		      cay_l = cay2*cay2;
		      break;
		      default:
		      cay_l = pow(cay,ly1);
		    }
		    lz1 = zpow_bf[iang][ibf];
		    switch (lz1) {
		      case 0:
		      caz_l = 1.0;
		      break;
		      case 1:
		      caz_l = caz;
		      break;
		      case 2:
		      caz_l = caz2;
		      break;
		      case 3:
		      caz_l = caz2*caz;
		      break;
		      case 4:
		      caz_l = caz2*caz2;
		      break;
		      default:
		      caz_l = pow(caz,lz1);
		    }
		    lx2 = xpow_bf[jang][jbf];
		    switch (lx2) {
		      case 0:
		      cbx_l = 1.0;
		      break;
		      case 1:
		      cbx_l = cbx;
		      break;
		      case 2:
		      cbx_l = cbx2;
		      break;
		      case 3:
		      cbx_l = cbx2*cbx;
		      break;
		      case 4:
		      cbx_l = cbx2*cbx2;
		      break;
		      default:
		      cbx_l = pow(cbx,lx2);
		    }
		    ly2 = ypow_bf[jang][jbf];
		    switch (ly2) {
		      case 0:
		      cby_l = 1.0;
		      break;
		      case 1:
		      cby_l = cby;
		      break;
		      case 2:
		      cby_l = cby2;
		      break;
		      case 3:
		      cby_l = cby2*cby;
		      break;
		      case 4:
		      cby_l = cby2*cby2;
		      break;
		      default:
		      cby_l = pow(cby,ly2);
		    }
		    lz2 = zpow_bf[jang][jbf];
		    switch (lz2) {
		      case 0:
		      cbz_l = 1.0;
		      break;
		      case 1:
		      cbz_l = cbz;
		      break;
		      case 2:
		      cbz_l = cbz2;
		      break;
		      case 3:
		      cbz_l = cbz2*cbz;
		      break;
		      case 4:
		      cbz_l = cbz2*cbz2;
		      break;
		      default:
		      cbz_l = pow(cbz,lz2);
		    }
		    tmp = dens_pf*cax_l*cay_l*caz_l*cbx_l*cby_l*cbz_l;
		    sum1 += tmp;
		  }
		}
		grid_pts[ix][iy] += sum1*prod_pf;
	      }
	      }
	      break;

		/* Electron density gradient vector */
	      case 3:
	      if (spin_prop)    /* If spin_prop is set, compute gradient of spin density... */
		Density = Pspin;
	      else              /* .. otherwise compute gradient of the electron density */
		Density = Ptot;
	      for(ix=0;ix<=nix;ix++) {
	        pcx = grid_origin[0] + grid_step_x[0]*ix - grid_step_y[0];
		pcy = grid_origin[1] + grid_step_x[1]*ix - grid_step_y[1];
		pcz = grid_origin[2] + grid_step_x[2]*ix - grid_step_y[2];
	      for(iy=0;iy<=niy;iy++) {
	        pcx += grid_step_y[0];
		pcy += grid_step_y[1];
		pcz += grid_step_y[2];
		cax = pcx-ax; cbx = pcx-bx;
		cax2 = cax*cax; cbx2 = cbx*cbx;
		cay = pcy-ay; cby = pcy-by;
		cay2 = cay*cay; cby2 = cby*cby;
		caz = pcz-az; cbz = pcz-bz;
		caz2 = caz*caz; cbz2 = cbz*cbz;
		prod_pf = contr[ig]*contr[jg]*exp(-ai*(cax2 + cay2 + caz2))*exp(-bj*(cbx2 + cby2 + cbz2));
		sum1 = sum2 = sum3 = 0.0;
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
                    dens_pf = Density[ioff[ib]+jb] * ((ib == jb) ? 1.0 : 2.0) * norm_pf;
		    if (fabs(dens_pf) < PFACCUTOFF) continue;
		    if (eq_shell && ig != jg)
                      dens_pf *= 2.0;
		    /* For **efficiency** reasons I'm trying to avoid calling pow() for up to
		       g-functions through switch() constructs */
		    lx1 = xpow_bf[iang][ibf];
		    switch (lx1) {
		      case 0:
		      cax_l_1 = cax_l = 1.0;
		      break;
		      case 1:
		      cax_l_1 = 1.0;
		      cax_l = cax;
		      break;
		      case 2:
		      cax_l_1 = cax;
		      cax_l = cax2;
		      break;
		      case 3:
		      cax_l_1 = cax2;
		      cax_l = cax2*cax;
		      break;
		      case 4:
		      cax_l_1 = cax2*cax;
		      cax_l = cax2*cax2;
		      break;
		      default:
		      cax_l_1 = pow(cax,lx1-1);
		      cax_l = cax_l_1*cax;
		    }
		    ly1 = ypow_bf[iang][ibf];
		    switch (ly1) {
		      case 0:
		      cay_l_1 = cay_l = 1.0;
		      break;
		      case 1:
		      cay_l_1 = 1.0;
		      cay_l = cay;
		      break;
		      case 2:
		      cay_l_1 = cay;
		      cay_l = cay2;
		      break;
		      case 3:
		      cay_l_1 = cay2;
		      cay_l = cay2*cay;
		      break;
		      case 4:
		      cay_l_1 = cay2*cay;
		      cay_l = cay2*cay2;
		      break;
		      default:
		      cay_l_1 = pow(cay,ly1-1);
		      cay_l = cay_l_1*cay;
		    }
		    lz1 = zpow_bf[iang][ibf];
		    switch (lz1) {
		      case 0:
		      caz_l_1 = caz_l = 1.0;
		      break;
		      case 1:
		      caz_l_1 = 1.0;
		      caz_l = caz;
		      break;
		      case 2:
		      caz_l_1 = caz;
		      caz_l = caz2;
		      break;
		      case 3:
		      caz_l_1 = caz2;
		      caz_l = caz2*caz;
		      break;
		      case 4:
		      caz_l_1 = caz2*caz;
		      caz_l = caz2*caz2;
		      break;
		      default:
		      caz_l_1 = pow(caz,lz1-1);
		      caz_l = caz_l_1*caz;
		    }
		    lx2 = xpow_bf[jang][jbf];
		    switch (lx2) {
		      case 0:
		      cbx_l_1 = cbx_l = 1.0;
		      break;
		      case 1:
		      cbx_l_1 = 1.0;
		      cbx_l = cbx;
		      break;
		      case 2:
		      cbx_l_1 = cbx;
		      cbx_l = cbx2;
		      break;
		      case 3:
		      cbx_l_1 = cbx2;
		      cbx_l = cbx2*cbx;
		      break;
		      case 4:
		      cbx_l_1 = cbx2*cbx;
		      cbx_l = cbx2*cbx2;
		      break;
		      default:
		      cbx_l_1 = pow(cbx,lx2-1);
		      cbx_l = cbx_l_1*cbx;
		    }
		    ly2 = ypow_bf[jang][jbf];
		    switch (ly2) {
		      case 0:
		      cby_l_1 = cby_l = 1.0;
		      break;
		      case 1:
		      cby_l_1 = 1.0;
		      cby_l = cby;
		      break;
		      case 2:
		      cby_l_1 = cby;
		      cby_l = cby2;
		      break;
		      case 3:
		      cby_l_1 = cby2;
		      cby_l = cby2*cby;
		      break;
		      case 4:
		      cby_l_1 = cby2*cby;
		      cby_l = cby2*cby2;
		      break;
		      default:
		      cby_l_1 = pow(cby,ly2-1);
		      cby_l = cby_l_1*cby;
		    }
		    lz2 = zpow_bf[jang][jbf];
		    switch (lz2) {
		      case 0:
		      cbz_l_1 = cbz_l = 1.0;
		      break;
		      case 1:
		      cbz_l_1 = 1.0;
		      cbz_l = cbz;
		      break;
		      case 2:
		      cbz_l_1 = cbz;
		      cbz_l = cbz2;
		      break;
		      case 3:
		      cbz_l_1 = cbz2;
		      cbz_l = cbz2*cbz;
		      break;
		      case 4:
		      cbz_l_1 = cbz2*cbz;
		      cbz_l = cbz2*cbz2;
		      break;
		      default:
		      cbz_l_1 = pow(cbz,lz2-1);
		      cbz_l = cbz_l_1*cbz;
		    }
		    x2 = cax_l*cbx_l;
		    y2 = cay_l*cby_l;
		    z2 = caz_l*cbz_l;
		      /* X-component of the gradient */
		    sum1 += dens_pf*y2*z2*((lx1*cax_l_1 - 2*ai*cax_l*cax)*cbx_l +
			                   (lx2*cbx_l_1 - 2*bj*cbx_l*cbx)*cax_l );
		      /* Y-component of the gradient */
		    sum2 += dens_pf*x2*z2*((ly1*cay_l_1 - 2*ai*cay_l*cay)*cby_l +
			                   (ly2*cby_l_1 - 2*bj*cby_l*cby)*cay_l );
		      /* Z-component of the gradient */
		    sum3 += dens_pf*x2*y2*((lz1*caz_l_1 - 2*ai*caz_l*caz)*cbz_l +
			                   (lz2*cbz_l_1 - 2*bj*cbz_l*cbz)*caz_l );
		  }
		}
		grid_vecX[ix][iy] += sum1*prod_pf;
		grid_vecY[ix][iy] += sum2*prod_pf;
		grid_vecZ[ix][iy] += sum3*prod_pf;
	      }
	      }
	      break;


	        /* Laplacian of electron density */
	      case 4:
	      if (spin_prop)    /* If spin_prop is set, compute the Laplacian of the spin density... */
		Density = Pspin;
	      else              /* .. otherwise compute the Laplacian of the electron density */
		Density = Ptot;
	      for(ix=0;ix<=nix;ix++) {
	        pcx = grid_origin[0] + grid_step_x[0]*ix - grid_step_y[0];
		pcy = grid_origin[1] + grid_step_x[1]*ix - grid_step_y[1];
		pcz = grid_origin[2] + grid_step_x[2]*ix - grid_step_y[2];
	      for(iy=0;iy<=niy;iy++) {
		pcx += grid_step_y[0];
		pcy += grid_step_y[1];
		pcz += grid_step_y[2];
		cax = pcx-ax; cbx = pcx-bx;
		cax2 = cax*cax; cbx2 = cbx*cbx;
		cay = pcy-ay; cby = pcy-by;
		cay2 = cay*cay; cby2 = cby*cby;
		caz = pcz-az; cbz = pcz-bz;
		caz2 = caz*caz; cbz2 = cbz*cbz;
		prod_pf = contr[ig]*contr[jg]*exp(-ai*(cax2 + cay2 + caz2))*exp(-bj*(cbx2 + cby2 + cbz2));
		sum1 = 0.0;
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
                    dens_pf = Density[ioff[ib]+jb] * ((ib == jb) ? 1.0 : 2.0) * norm_pf;
		    if (fabs(dens_pf) < PFACCUTOFF) continue;
		    if (eq_shell && ig != jg)
                      dens_pf *= 2.0;
		    lx1 = xpow_bf[iang][ibf];
		    ly1 = ypow_bf[iang][ibf];
		    lz1 = zpow_bf[iang][ibf];
		    /* For **efficiency** reasons I'm trying to avoid calling pow() for up to
			 g-functions through switch() constructs */
		    switch (lx1) {
		      case 0:
		      cax_l_2 = cax_l_1 = cax_l = 1.0;
		      break;
		      case 1:
		      cax_l_2 = cax_l_1 = 1.0;
		      cax_l = cax;
		      break;
		      case 2:
		      cax_l_2 = 1.0;
		      cax_l_1 = cax;
		      cax_l = cax2;
		      break;
		      case 3:
		      cax_l_2 = cax;
		      cax_l_1 = cax2;
		      cax_l = cax2*cax;
		      break;
		      case 4:
		      cax_l_2 = cax2;
		      cax_l_1 = cax2*cax;
		      cax_l = cax2*cax2;
		      break;
		      default:
		      cax_l_2 = pow(cax,lx1-2);
		      cax_l_1 = cax_l_2*cax;
		      cax_l = cax_l_1*cax;
		    }
		    switch (ly1) {
		      case 0:
		      cay_l_2 = cay_l_1 = cay_l = 1.0;
		      break;
		      case 1:
		      cay_l_2 = cay_l_1 = 1.0;
		      cay_l = cay;
		      break;
		      case 2:
		      cay_l_2 = 1.0;
		      cay_l_1 = cay;
		      cay_l = cay2;
		      break;
		      case 3:
		      cay_l_2 = cay;
		      cay_l_1 = cay2;
		      cay_l = cay2*cay;
		      break;
		      case 4:
		      cay_l_2 = cay2;
		      cay_l_1 = cay2*cay;
		      cay_l = cay2*cay2;
		      break;
		      default:
		      cay_l_2 = pow(cay,ly1-2);
		      cay_l_1 = cay_l_2*cay;
		      cay_l = cay_l_1*cay;
		    }
		    switch (lz1) {
		      case 0:
		      caz_l_2 = caz_l_1 = caz_l = 1.0;
		      break;
		      case 1:
		      caz_l_2 = caz_l_1 = 1.0;
		      caz_l = caz;
		      break;
		      case 2:
		      caz_l_2 = 1.0;
		      caz_l_1 = caz;
		      caz_l = caz2;
		      break;
		      case 3:
		      caz_l_2 = caz;
		      caz_l_1 = caz2;
		      caz_l = caz2*caz;
		      break;
		      case 4:
		      caz_l_2 = caz2;
		      caz_l_1 = caz2*caz;
		      caz_l = caz2*caz2;
		      break;
		      default:
		      caz_l_2 = pow(caz,lz1-2);
		      caz_l_1 = caz_l_2*caz;
		      caz_l = caz_l_1*caz;
		    }
		    lx2 = xpow_bf[jang][jbf];
		    ly2 = ypow_bf[jang][jbf];
		    lz2 = zpow_bf[jang][jbf];
		    switch (lx2) {
		      case 0:
		      cbx_l_2 = cbx_l_1 = cbx_l = 1.0;
		      break;
		      case 1:
		      cbx_l_2 = cbx_l_1 = 1.0;
		      cbx_l = cbx;
		      break;
		      case 2:
		      cbx_l_2 = 1.0;
		      cbx_l_1 = cbx;
		      cbx_l = cbx2;
		      break;
		      case 3:
		      cbx_l_2 = cbx;
		      cbx_l_1 = cbx2;
		      cbx_l = cbx2*cbx;
		      break;
		      case 4:
		      cbx_l_2 = cbx2;
		      cbx_l_1 = cbx2*cbx;
		      cbx_l = cbx2*cbx2;
		      break;
		      default:
		      cbx_l_2 = pow(cbx,lx2-2);
		      cbx_l_1 = cbx_l_2*cbx;
		      cbx_l = cbx_l_1*cbx;
		    }
		    switch (ly2) {
		      case 0:
		      cby_l_2 = cby_l_1 = cby_l = 1.0;
		      break;
		      case 1:
		      cby_l_2 = cby_l_1 = 1.0;
		      cby_l = cby;
		      break;
		      case 2:
		      cby_l_2 = 1.0;
		      cby_l_1 = cby;
		      cby_l = cby2;
		      break;
		      case 3:
		      cby_l_2 = cby;
		      cby_l_1 = cby2;
		      cby_l = cby2*cby;
		      break;
		      case 4:
		      cby_l_2 = cby2;
		      cby_l_1 = cby2*cby;
		      cby_l = cby2*cby2;
		      break;
		      default:
		      cby_l_2 = pow(cby,ly2-2);
		      cby_l_1 = cby_l_2*cby;
		      cby_l = cby_l_1*cby;
		    }
		    switch (lz2) {
		      case 0:
		      cbz_l_2 = cbz_l_1 = cbz_l = 1.0;
		      break;
		      case 1:
		      cbz_l_2 = cbz_l_1 = 1.0;
		      cbz_l = cbz;
		      break;
		      case 2:
		      cbz_l_2 = 1.0;
		      cbz_l_1 = cbz;
		      cbz_l = cbz2;
		      break;
		      case 3:
		      cbz_l_2 = cbz;
		      cbz_l_1 = cbz2;
		      cbz_l = cbz2*cbz;
		      break;
		      case 4:
		      cbz_l_2 = cbz2;
		      cbz_l_1 = cbz2*cbz;
		      cbz_l = cbz2*cbz2;
		      break;
		      default:
		      cbz_l_2 = pow(cbz,lz2-2);
		      cbz_l_1 = cbz_l_2*cbz;
		      cbz_l = cbz_l_1*cbz;
		    }
		    x2 = cax_l*cbx_l;
		    y2 = cay_l*cby_l;
		    z2 = caz_l*cbz_l;
		      /* d2/dx2 part */
		    sum2 = -y2*z2*
			 ( (lx1*(lx1-1)*cax_l_2 + 4*ai*ai*cax_l*cax2 - 2*ai*(2*lx1+1)*cax_l)*cbx_l +
			   (lx2*(lx2-1)*cbx_l_2 + 4*bj*bj*cbx_l*cbx2 - 2*bj*(2*lx2+1)*cbx_l)*cax_l +
			   2*(lx1*cax_l_1 - 2*ai*cax*cax_l)*(lx2*cbx_l_1 - 2*bj*cbx*cbx_l) );
		      /* d2/dy2 part */
		    sum2 -= x2*z2*
			 ( (ly1*(ly1-1)*cay_l_2 + 4*ai*ai*cay_l*cay2 - 2*ai*(2*ly1+1)*cay_l)*cby_l +
			   (ly2*(ly2-1)*cby_l_2 + 4*bj*bj*cby_l*cby2 - 2*bj*(2*ly2+1)*cby_l)*cay_l +
			   2*(ly1*cay_l_1 - 2*ai*cay*cay_l)*(ly2*cby_l_1 - 2*bj*cby*cby_l) );
		      /* d2/dz2 part */
		    sum2 -= x2*y2*
			 ( (lz1*(lz1-1)*caz_l_2 + 4*ai*ai*caz_l*caz2 - 2*ai*(2*lz1+1)*caz_l)*cbz_l +
			   (lz2*(lz2-1)*cbz_l_2 + 4*bj*bj*cbz_l*cbz2 - 2*bj*(2*lz2+1)*cbz_l)*caz_l +
			   2*(lz1*caz_l_1 - 2*ai*caz*caz_l)*(lz2*cbz_l_1 - 2*bj*cbz*cbz_l) );
		      /* Contract with the density */
		    sum1 += dens_pf*sum2;
		  }
		}
	      grid_pts[ix][iy] += sum1*prod_pf;
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


        /* Projecting density gradient onto the plane of the grid and computing the magnitude */

  if (grid == 3)
    for(ix=0;ix<=nix;ix++)
      for(iy=0;iy<=niy;iy++) {
	x = grid_vecX[ix][iy]*grid_unit_x[0] + grid_vecY[ix][iy]*grid_unit_x[1] + grid_vecZ[ix][iy]*grid_unit_x[2];
        y = grid_vecX[ix][iy]*grid_unit_y[0] + grid_vecY[ix][iy]*grid_unit_y[1] + grid_vecZ[ix][iy]*grid_unit_y[2];
	grid_pts[ix][iy] = sqrt(x*x + y*y);
	if (edgrad_logscale) {
	  tmp = (double)edgrad_logscale*log(1.0 + grid_pts[ix][iy]);
	  x /= tmp;
	  y /= tmp;
	  grid_pts[ix][iy] /= tmp;
	}
	grid_vecX[ix][iy] = x;
	grid_vecY[ix][iy] = y;
      }


    

	/* Cleaning up */

  return;
}

}} // namespace psi::oeprop
