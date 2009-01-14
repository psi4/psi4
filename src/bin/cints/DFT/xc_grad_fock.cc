/*! \file xc_grad_fock.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cmath>
#include<cstring>
#include<cstdio>
#include<cstdlib>
#include<libipv1/ip_lib.h>
#include<libciomr/libciomr.h>
#include<libpsio/psio.h>
#include<libint/libint.h>
#include<pthread.h>
#include<libqt/qt.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

#include"dft_init.h"
#include"weighting.h"
#include"calc_den_fast.h"
#include"calc_den.h"
#include"calc_grad_fast.h"
#include"functional.h"
#include"physconst.h"
#include"grid_init.h"
#include"free_grid_structs.h"
#include"dcr.h"
#include"init_close_shell_info.h"
#include"calc_close_basis.h"

namespace psi { namespace CINTS {

void xc_grad_fock(void){
  int i,j,k,l,m,n,q,s,t,u;
  int ua, atom, ua_deg;
  int rpoints;
  int ang_points;
  int num_ao;
  int point_count=0;
  int dum;
  int num;
  int moff,noff,mtmp,ntmp;
  int am2shell1,am2shell2;
  int shell_type1,shell_type2;
  int ioff1,ioff2;
  int chek = 1;
  int close_shells;
  int close_aos;
  int tmp;
  
  int nstri;
  double temp;
  double *temp_arr;
  double *omega_arr;
  double **tmpmat1;
  double *Gtri;  /* Total and open-shell 
			     G matrices and lower 
			     triagonal form in SO basis */  
  double r;
  double rind;
  double rpoints_double;
  double ua_deg_d;
  double bragg;
  double qr;
  double drdq;
  double jacobian;
  double xa,ya,za;
  double Becke_weight;
  double ang_quad;
  double four_pi_div_by_rps;
  double den_val=0.0;
  double exch_pval=0.0;
  double exch_gval=0.0;
  double corr_pval=0.0;
  double corr_gval=0.0;
  double exch_eval=0.0;
  double corr_eval=0.0;
  double pval = 0.0;
  double gval = 0.0;
  double eval = 0.0;
  double bas1 = 0.0;
  double bas2 = 0.0;
  double vvalbas = 0.0;
  double dpx,dpy,dpz;
  
  struct coordinates geom;
  struct den_info_s den_info;
  struct xc_info_s xc_info;
  
  struct atomic_grid_s *atm_grd;
  struct leb_chunk_s *chnk;
  leb_sphere_t *sphr;
  leb_point_t *pnt;
  
  
  num_ao = BasisSet.num_ao;
  
  DFT_options.basis = (double *)malloc(sizeof(double)*num_ao);
  DFT_options.gradx = (double *)malloc(sizeof(double)*num_ao);
  DFT_options.grady = (double *)malloc(sizeof(double)*num_ao);
  DFT_options.gradz = (double *)malloc(sizeof(double)*num_ao);
  DFT_options.gamma_basis = (double *)malloc(sizeof(double)*num_ao);
  omega_arr = (double *)malloc(sizeof(double)*num_ao);
  
  G = init_matrix(num_ao,num_ao);
  
  /* ----Initialize Close shell data structure ---- */
  
  DFT_options.close_shell_info = init_close_shell_info();

  timer_init();
  timer_on("DFT");
  
  grid_init();
  
/*-------------------------------------------------------
    Loop over symmetry-unique atoms only since integration
    domains around symm.-unique atoms are equivalent
    We are NOT employing the symmetry of angular grids
    about atoms in special positions (on symm. elements)
    like Handy & co. do
   -------------------------------------------------------*/
  for(ua=0;ua<Symmetry.num_unique_atoms;ua++){
      atm_grd = &(DFT_options.grid.atomic_grid[ua]);
      
      atom = atm_grd->atom_num;
      
/*--- Cheap trick to get degeneracies of each unique atom ---*/
      
      ua_deg = atm_grd->atom_degen;
      ua_deg_d = (double) ua_deg;

      xa = atm_grd->atom_center.x;
      ya = atm_grd->atom_center.y;
      za = atm_grd->atom_center.z;
      
      bragg = atm_grd->Bragg_radii;
      
      for(j=0;j<atm_grd->chunk_num;j++){
	  chnk = &(atm_grd->leb_chunk[j]);
	  timer_on("close basis");
	  
	  calc_close_basis(ua,j);
	  
	  close_shells = DFT_options.close_shell_info.num_close_shells;
	  close_aos = DFT_options.close_shell_info.num_close_aos;
	  timer_off("close basis");
	  for(k=0;k<chnk->size;k++){
	      
	      sphr = &(chnk->spheres[k]);
	      
	      r = sphr->r*bragg;
	      /*r = 0.87895;*/
	      
	      
	      drdq = sphr->drdq*bragg*bragg*bragg;
	      
	      for(l=0;l<sphr->n_ang_points;l++){
		  pnt = &(sphr->points[l]);
		  /* ----------------------------------
		     Calculate the cartesian points of the point
		     relative to the center of mass
		     -------------------------------------*/
      		  geom.x = bragg*pnt->p_cart.x+xa;
		  geom.y = bragg*pnt->p_cart.y+ya;
		  geom.z = bragg*pnt->p_cart.z+za;
		  
		  /*-----------------------------------
		    Calculate the weighting funtion 
		    ----------------------------------*/
		  Becke_weight = weight_calc(atom,geom,3);
		  
		  if(Becke_weight> WEIGHT_CUTOFF){
		      /*-----------------------------------
			Get the density information for this 
			point
			----------------------------------*/
		      timer_on("DEN1");
		      /*den_info = DFT_options.den_calc(geom,ua);*/
		      den_info = calc_grad_fast(geom,ua);
		      
		      timer_off("DEN1");
		      
		      if(den_info.den > DEN_CUTOFF){
			  /*-------------------------------------
			    Weight from Lebedev
			    -----------------------------------*/
		      
			  ang_quad = 4.0*_pi*pnt->ang_weight;
		      
			  /*------------------------------------
			    Calculate the potential functional
			    and energy functional at this
			    point
			    -----------------------------------*/
			  /*fprintf(outfile,"\nua_deg = %10.10lf",ua_deg_d);*/
			  den_val += 2.0*ua_deg_d*drdq*ang_quad
			      *Becke_weight*den_info.den;
			  
			  xc_info.exch_info = DFT_options.
			      exchange_func(den_info);
			  
			  xc_info.corr_info = DFT_options.
			      correlation_func(den_info);
			 
			  exch_pval = ua_deg_d*drdq*ang_quad
			      *Becke_weight*xc_info.exch_info.dpval;
			  exch_gval = ua_deg_d*drdq*ang_quad
			      *Becke_weight*xc_info.exch_info.dgval;
	 
			  corr_pval = ua_deg_d*drdq*ang_quad
			      *Becke_weight*xc_info.corr_info.dpval;
			  corr_gval = ua_deg_d*drdq*ang_quad
			      *Becke_weight*xc_info.corr_info.dgval;

			  pval = exch_pval+corr_pval;
			  gval = exch_gval+corr_gval;
			  
			  exch_eval += ua_deg_d*drdq*ang_quad
			      *Becke_weight*xc_info.exch_info.eval;
			  corr_eval += ua_deg_d*drdq*ang_quad
			      *Becke_weight*xc_info.corr_info.eval;
			  eval += ua_deg_d*drdq*ang_quad
			      *Becke_weight*(xc_info.exch_info.eval+
					     xc_info.corr_info.eval);
			  
			  /*------------------------------------
			    Update the G matrix
			    -----------------------------------*/
			  timer_on("FOCK");
			  t=0;
			  
			  /* Form omega array <  p+dphiu  > */
			  
			  dpx = den_info.gradx;
			  dpy = den_info.grady;
			  dpz = den_info.gradz;
			  
			  /* Omega is the df/dgamma * delp . delp(basu*basv) */
			  /* See Benny Johnson's Thesis page 109 */
			  
			  /*for(m=0;m<close_aos;m++){
			      omega_arr[m] = 2.0*gval*(dpx*DFT_options.gradx[m]
						       +dpy*DFT_options.grady[m]
						       +dpz*DFT_options.gradz[m]);
			  }
			 
		      for(m=0;m<close_aos;m++){
			      bas1 = pval*DFT_options.basis[m]
				  + omega_arr[m];
			      moff = DFT_options.close_shell_info.
				  aos_close_to_chunk[m];
			      for(n=m;n<close_aos;n++){
				  bas2 = bas1*DFT_options.basis[n]
				      +DFT_options.basis[m]*omega_arr[n];
				  noff = DFT_options.close_shell_info.
				      aos_close_to_chunk[n];
				  if(noff > moff)
				      G[moff][noff] += bas2;
				  else
				      G[noff][moff] += bas2;
			      }
			  }*/
		      /*fprintf(outfile,"\ngval = %10.10lf",gval);*/
		      for(m=0;m<close_aos;m++){
			  for(n=0;n<close_aos;n++){
			      moff = DFT_options.close_shell_info.
				  aos_close_to_chunk[m];
			      noff = DFT_options.close_shell_info.
				      aos_close_to_chunk[n];
			      bas1 = pval*DFT_options.basis[m]
				  *DFT_options.basis[n];
			      bas2 = 2.0*gval*((dpx*DFT_options.gradx[m]+
					      dpy*DFT_options.grady[m]+
					      dpz*DFT_options.gradz[m])*
					     DFT_options.basis[n]
					     +(dpx*DFT_options.gradx[n]+
					      dpy*DFT_options.grady[n]+
					      dpz*DFT_options.gradz[n])*
					     DFT_options.basis[m]);
			      /*fprintf(outfile,"\nbas1 = %10.10lf bas2 = %10.10lf\n\n\n",bas1,bas2);*/
			      G[moff][noff] += pval*DFT_options.basis[m]
				  *DFT_options.basis[n]
				  +2.0*gval*((dpx*DFT_options.gradx[m]+
					      dpy*DFT_options.grady[m]+
					      dpz*DFT_options.gradz[m])*
					     DFT_options.basis[n]
					     +(dpx*DFT_options.gradx[n]+
					      dpy*DFT_options.grady[n]+
					      dpz*DFT_options.gradz[n])*
					     DFT_options.basis[m]);
			  }
		      }
		      timer_off("FOCK");      
		      
		      }
		  }
	      }  
	  }
      }
  }
  
  free_close_shell_info(DFT_options.close_shell_info);
  print_mat(G,num_ao,num_ao,outfile);
  for(m=0;m<num_ao;m++)
      for(n=m;n<num_ao;n++)
	  G[n][m]=G[m][n];
  
  timer_off("DFT");
  timer_done();
  
  
  /*----------------------
    Unsort the Fock matrix
    back to shell ordering
    ---------------------*/
  
  
/*----------------------
  Transform to SO basis
  ----------------------*/
  if (Symmetry.nirreps > 1 || BasisSet.puream) {
      tmpmat1 = block_matrix(Symmetry.num_so,BasisSet.num_ao);
      mmult(Symmetry.usotao,0,G,0,tmpmat1,0,Symmetry.num_so,BasisSet.num_ao,BasisSet.num_ao,0);
      mmult(tmpmat1,0,Symmetry.usotao,1,G,0,Symmetry.num_so,BasisSet.num_ao,Symmetry.num_so,0);
      if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
	  mmult(Symmetry.usotao,0,Go,0,tmpmat1,0,Symmetry.num_so,BasisSet.num_ao,BasisSet.num_ao,0);
	  mmult(tmpmat1,0,Symmetry.usotao,1,Go,0,Symmetry.num_so,BasisSet.num_ao,Symmetry.num_so,0);
      }
      free_block(tmpmat1);
  }
  
  /*-------------------------
    Write G-matrices to disk
    -------------------------*/
  
  nstri = ioff[Symmetry.num_so];
  Gtri = init_array(nstri);
  sq_to_tri(G,Gtri,Symmetry.num_so);
  free_block(G);
  fprintf(outfile,"\nDFT_energy = %10.10lf",eval);
  fprintf(outfile,"\nX-Energy = %10.10lf",exch_eval);
  fprintf(outfile,"\nC-Energy = %10.10lf",corr_eval);
  fprintf(outfile,"\ntrace of density = %10.10lf\n",den_val);
  
  psio_open(IOUnits.itapDSCF, PSIO_OPEN_OLD);
  psio_write_entry(IOUnits.itapDSCF,"DFT X-energy",
		   (char *) &(exch_eval), sizeof(double));
  psio_write_entry(IOUnits.itapDSCF,"DFT C-energy",
		   (char *) &(corr_eval), sizeof(double));
  psio_write_entry(IOUnits.itapDSCF,"DFT XC-energy",
		   (char *) &(eval), sizeof(double));
  psio_write_entry(IOUnits.itapDSCF,"DFT Den",
		   (char *) &(den_val), sizeof(double));
  
  psio_write_entry(IOUnits.itapDSCF, "Total XC G-matrix"
		   , (char *) Gtri, sizeof(double)*nstri);
  free(Gtri);
 
  /*-- Cleanup the DFT stuff and close files --*/
  cleanup_grid_type(DFT_options.grid);
  psio_close(IOUnits.itapDSCF, 1);
  return;
}
}}
