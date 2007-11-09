/*! \file xc_fock_u.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<cmath>
#include<cstring>
#include<stdio.h>
#include<stdlib.h>
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
#include"calc_den_u.h"
#include"functional_u.h"
#include"physconst.h"
#include"grid_init.h"
#include"free_grid_structs.h"
#include"dcr.h"
#include"init_close_shell_info.h"
#include"calc_close_basis_u.h"
#include"pade.h"

namespace psi { namespace CINTS {
 
void xc_fock_u(void){
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
  double **tmpmat1;
  double *Gtria, *Gtrib;  /* Total and open-shell 
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
  double exch_vval_a=0.0;
  double exch_vval_b=0.0;
  double corr_vval_a=0.0;
  double corr_vval_b=0.0;
  double exch_eval=0.0;
  double corr_eval=0.0;
  double vval_a = 0.0;
  double vval_b = 0.0;
  double eval = 0.0;
  double bas1 = 0.0;
  double bas1a = 0.0;
  double bas1b = 0.0;
  double bas2 = 0.0;
  double bas2a = 0.0;
  double bas2b = 0.0;
  
  struct coordinates geom;
  struct den_info_s den_info;
  struct xc_info_s xc_info;
  
  struct atomic_grid_s *atm_grd;
  struct leb_chunk_s *chnk;
  leb_sphere_t *sphr;
  leb_point_t *pnt;
 
  /*fprintf(outfile,"\nPade = %10.10lf",Pade_int(1.0,-0.00475840
					       ,1.13107,13.0045,
					       -0.0337737,7.123108918));*/
  num_ao = BasisSet.num_ao;
  DFT_options.basis = (double *)malloc(sizeof(double)*num_ao);
  
  Ga = block_matrix(num_ao,num_ao);
  Gb = block_matrix(num_ao,num_ao);
  
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
	  calc_close_basis_u(ua,j);
	  close_shells = DFT_options.close_shell_info.num_close_shells;
	  close_aos = DFT_options.close_shell_info.num_close_aos;
	  timer_off("close basis");
	  for(k=0;k<chnk->size;k++){
	      
	      sphr = &(chnk->spheres[k]);
	      
	      r = sphr->r*bragg;
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
                      den_info = calc_density_u(geom,ua);
		      timer_off("DEN1");
		      if(den_info.dena > DEN_CUTOFF || 
			 den_info.denb > DEN_CUTOFF){
			  /*-------------------------------------
			    Weight from Lebedev
			    -----------------------------------*/
		      
			  ang_quad = pnt->ang_weight;
		      
			  /*------------------------------------
			    Calculate the potential functional
			    and energy functional at this
			    point
			    -----------------------------------*/
			 
			  den_val += ua_deg_d*drdq*ang_quad
			      *Becke_weight*(den_info.dena + den_info.denb);
			 
			  xc_info.exch_info = DFT_options.
			      exchange_func(den_info);
			  
			  xc_info.corr_info = DFT_options.
			      correlation_func(den_info);
			  
			  exch_vval_a = ua_deg_d*drdq*ang_quad
			      *Becke_weight*xc_info.exch_info.dvala;
			  exch_vval_b = ua_deg_d*drdq*ang_quad
			      *Becke_weight*xc_info.exch_info.dvalb;
			  
			  corr_vval_a = ua_deg_d*drdq*ang_quad
			      *Becke_weight*xc_info.corr_info.dvala;
			  corr_vval_b = ua_deg_d*drdq*ang_quad
			      *Becke_weight*xc_info.corr_info.dvalb;
			  
			  vval_a = exch_vval_a+corr_vval_a;
			  vval_b = exch_vval_b+corr_vval_b;
			  
			  exch_eval += ua_deg_d*drdq*ang_quad
			      *Becke_weight*xc_info.exch_info.eval;
			  corr_eval += ua_deg_d*drdq*ang_quad
			      *Becke_weight*xc_info.corr_info.eval;;
			  eval += ua_deg_d*drdq*ang_quad
			      *Becke_weight*(xc_info.exch_info.eval+
					     xc_info.corr_info.eval);
			  
			  
			  
			  /*------------------------------------
			    Update the G matrix
			    -----------------------------------*/
			  timer_on("FOCK");
			  t=0;
			  
			  for(m=0;m<close_aos;m++){
			      bas1 = DFT_options.basis[m];
			      bas1a = vval_a*bas1;
			      bas1b = vval_b*bas1;
			      moff = DFT_options.close_shell_info.
				  aos_close_to_chunk[m];
			      for(n=m;n<close_aos;n++){
				  bas2 = DFT_options.basis[n];
				  bas2a = bas1a*bas2;
				  bas2b = bas1b*bas2;
				  noff = DFT_options.close_shell_info.
				      aos_close_to_chunk[n];
				  if(noff > moff){
				      Ga[moff][noff] += bas2a;
				      Gb[moff][noff] += bas2b;
				  }				  
				  else{
				      Ga[noff][moff] += bas2a;
				      Gb[noff][moff] += bas2b;
				  }
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
  
  for(m=0;m<num_ao;m++)
      for(n=m;n<num_ao;n++)
	  Ga[n][m]=Ga[m][n];
  for(m=0;m<num_ao;m++)
      for(n=m;n<num_ao;n++)
	  Gb[n][m]=Gb[m][n];
  /*print_mat(Ga,num_ao,num_ao,outfile);
    print_mat(Gb,num_ao,num_ao,outfile);*/
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
      mmult(Symmetry.usotao,0,Ga,0,tmpmat1,0,Symmetry.num_so,
	    BasisSet.num_ao,BasisSet.num_ao,0);
      mmult(tmpmat1,0,Symmetry.usotao,1,Ga,0,Symmetry.num_so,
	    BasisSet.num_ao,Symmetry.num_so,0);
      mmult(Symmetry.usotao,0,Gb,0,tmpmat1,0,Symmetry.num_so,
	    BasisSet.num_ao,BasisSet.num_ao,0);
      mmult(tmpmat1,0,Symmetry.usotao,1,Gb,0,Symmetry.num_so,
	    BasisSet.num_ao,Symmetry.num_so,0);
      
      free_block(tmpmat1);
  }
  
  /*-------------------------
    Write G-matrices to disk
    -------------------------*/  
 

  
  /*fprintf(outfile,"\nDFT_energy = %10.10lf",eval);
  fprintf(outfile,"\nX-Energy = %10.10lf",exch_eval);
  fprintf(outfile,"\nC-Energy = %10.10lf",corr_eval);
  fprintf(outfile,"\ntrace of density = %10.10lf\n",den_val);*/
  psio_open(IOUnits.itapDSCF, PSIO_OPEN_OLD);
  psio_write_entry(IOUnits.itapDSCF,"DFT X-energy",
		   (char *) &(exch_eval), sizeof(double));
  psio_write_entry(IOUnits.itapDSCF,"DFT C-energy",
		   (char *) &(corr_eval), sizeof(double));
  psio_write_entry(IOUnits.itapDSCF,"DFT XC-energy",
		   (char *) &(eval), sizeof(double));
  psio_write_entry(IOUnits.itapDSCF,"DFT Den",
		   (char *) &(den_val), sizeof(double));
  nstri = ioff[Symmetry.num_so];
  Gtria = init_array(nstri);
  sq_to_tri(Ga,Gtria,Symmetry.num_so);
  free_block(Ga);
  Gtrib = init_array(nstri);
  sq_to_tri(Gb,Gtrib,Symmetry.num_so);
  free_block(Gb);
  
  psio_write_entry(IOUnits.itapDSCF, "Alpha XC G-matrix"
		   , (char *) Gtria, sizeof(double)*nstri);
  psio_write_entry(IOUnits.itapDSCF, "Beta XC G-matrix"
		   , (char *) Gtrib, sizeof(double)*nstri);
  free(Gtria);
  free(Gtrib);

  /*-- Cleanup the DFT stuff and close files --*/
  cleanup_grid_type(DFT_options.grid);
  psio_close(IOUnits.itapDSCF, 1);
  return;
}
};}
