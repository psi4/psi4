/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "includes.h"
#include "globals.h"
#include "prototypes.h"

namespace psi { namespace oeprop {

void compute_grid_mos()
{
  int i,j,k,l,ig,jg,ibf,jbf,ib,jb,jlim,kk,ll;
  int iatom, jatom, iang, jang, i_1stbf, j_1stbf, nbfi, nbfj;
  int igmin, igmax;
  int ixm, iym, izm, iind;
  int ix, iy, iz;
  int lx, ly, lz;
  int mo, mo_num;
  double ax, ay, az, xa, ya, za, ra2;
  double ai, norm_pf, ang_pf, exponent;
  double x,y,z;
  double r,r2,tmp;
  double *bf_values, *values;


  /* Initialize some intermediates */
  grid3d_pts = (double****) malloc(num_mos_to_plot*sizeof(double***));
  for(mo=0;mo<num_mos_to_plot;mo++)
    grid3d_pts[mo] = init_box(nix+1,niy+1,niz+1);
  bf_values = init_array(nbfao);

  for(ix=0;ix<=nix;ix++) {
      for(iy=0;iy<=niy;iy++) {
	x = grid_origin[0] + grid_step_x[0]*ix + grid_step_y[0]*iy;
	y = grid_origin[1] + grid_step_x[1]*ix + grid_step_y[1]*iy;
	z = grid_origin[2] + grid_step_x[2]*ix + grid_step_y[2]*iy;

	for(iz=0;iz<=niz;iz++) {
	  /* Shell loop */
	  for(i=0;i<nshell;i++) {
	      iatom = snuc[i] - 1;
	      ax = geom[iatom][0];
	      ay = geom[iatom][1];
	      az = geom[iatom][2];
	      iang = stype[i]-1;
	      izm = 1;
	      iym = iang+1;
	      ixm = iym*iym;
	      i_1stbf = sloc[i] - 1;
	      nbfi = (iang+2)*(iang+1)/2;
	      igmin = sprim[i] - 1;
	      igmax = igmin + snumg[i] - 1;
    
	      xa = x - ax;
	      ya = y - ay;
	      za = z - az;
	      ra2 = xa*xa + ya*ya + za*za;

	      /*--- Compute exponentional factor - loop over primitives ---*/
	      exponent = 0.0;
	      for(ig=igmin;ig<=igmax;ig++) {
		  ai = exps[ig];
		  exponent += contr[ig]*exp(-ai*ra2);
	      }

	      /*--- Loop over basis functions in the shell ---*/
	      values = &(bf_values[i_1stbf]);
	      for(ibf=0;ibf<nbfi;ibf++,values++) {
		  norm_pf = norm_bf[iang][ibf];
		  lx = xpow_bf[iang][ibf];
		  ly = ypow_bf[iang][ibf];
		  lz = zpow_bf[iang][ibf];
		  tmp = 1.0;
		  switch (lx) {
		  case 0:
		      break;
		  case 1:
		      tmp *= xa;
		      break;
		  case 2:
		      tmp *= xa*xa;
		      break;
		  case 3:
		      tmp *= xa*xa*xa;
		      break;
		  case 4:
		      tmp *= (xa*xa)*(xa*xa);
		      break;
		  default:
		      tmp *= pow(xa,lx);
		  }
		  ang_pf = tmp;
		  tmp = 1.0;
		  switch (ly) {
		  case 0:
		      break;
		  case 1:
		      tmp *= ya;
		      break;
		  case 2:
		      tmp *= ya*ya;
		      break;
		  case 3:
		      tmp *= ya*ya*ya;
		      break;
		  case 4:
		      tmp *= (ya*ya)*(ya*ya);
		      break;
		  default:
		      tmp *= pow(ya,ly);
		  }
		  ang_pf *= tmp;
		  tmp = 1.0;
		  switch (lz) {
		  case 0:
		      break;
		  case 1:
		      tmp *= za;
		      break;
		  case 2:
		      tmp *= za*za;
		      break;
		  case 3:
		      tmp *= za*za*za;
		      break;
		  case 4:
		      tmp *= (za*za)*(za*za);
		      break;
		  default:
		      tmp *= pow(za,lz);
		  }
		  ang_pf *= tmp;
		  
		  *values = norm_pf * ang_pf * exponent;
	      }
	  } /*--- end of shell loop ---*/
	  for(mo=0;mo<num_mos_to_plot;mo++) {
	    tmp = 0;
	    mo_num = mos_to_plot[mo];
	    for(i=0;i<nbfao;i++)
	      tmp += bf_values[i]*scf_evec_ao[i][mo_num];
	    grid3d_pts[mo][ix][iy][iz] = tmp;
	  }
          x += grid_step_z[0];
	  y += grid_step_z[1];
	  z += grid_step_z[2];
	}
      }
  }
}

}} // namespace psi::oeprop
