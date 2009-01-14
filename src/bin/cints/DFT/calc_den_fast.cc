/*! \file calc_den_fast.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
    \author Shawn Brown
    
    The code contains functions that will retrieve 
    and process the density at a given x,y,z coord
    in space
    
    --------------------------------------------------*/
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"bas_comp_functions.h"

namespace psi { namespace CINTS {
  
  struct den_info_s calc_density_fast(struct coordinates geom, int atom_num){
    
    int i,j,k,l,m;
    int am2shell;
    int shell_type;
    int shell_start;
    int shell_end;
    int n_shells;
    int num_ao,ndocc;
    int shell_center;
    double x,y,z;
    double xa,ya,za;
    double rr;
    double rrtmp;
    double bastmp;
    double den_sum;
    double coeff;
    double expon;
    double *norm_ptr;
    double *dist_atom;
    double *temp_arr;
    double *coord_array[14];
    double *expon_arr;
    double norm_ptr0,norm_ptr1,norm_ptr2;
    double norm_ptr3,norm_ptr4,norm_ptr5;
    double norm_ptr6,norm_ptr7,norm_ptr8;
    double norm_ptr9,norm_ptr10,norm_ptr11;
    double norm_ptr12,norm_ptr13,norm_ptr14;
    double xx,yy,zz,xy,xz,yz;
    double bxx,byy,bzz,bxy,bxz,byz;

    
    struct coordinates *dist_coord;
    struct leb_chunk_s *chunk;
    struct den_info_s den_info;
    struct close_shell_info_s *close;
    
    x = geom.x;
    y = geom.y;
    z = geom.z;
    
    num_ao = BasisSet.num_ao;
    ndocc = MOInfo.ndocc;
    temp_arr = init_array(ndocc);
    close = &(DFT_options.close_shell_info);
    /* ---------------------------------
       Compute distances from atom that 
       the basis function is centered on
       to the grid point
       --------------------------------*/
    zero_arr(DFT_options.basis,num_ao);
    dist_atom = init_array(Molecule.num_atoms);
    dist_coord = (struct coordinates *)
      malloc(sizeof(struct coordinates)*Molecule.num_atoms);
    timer_on("distance");
    for(i=0;i<Molecule.num_atoms;i++){
      dist_coord[i].x = x-Molecule.centers[i].x;
      dist_coord[i].y = y-Molecule.centers[i].y;
      dist_coord[i].z = z-Molecule.centers[i].z;
      dist_atom[i] = dist_coord[i].x*dist_coord[i].x
	+dist_coord[i].y*dist_coord[i].y
	+dist_coord[i].z*dist_coord[i].z;
    }
    
    n_shells = BasisSet.num_shells;
    timer_off("distance");
    timer_on("basis");
    
    l=0;
    m=0;
    for(i=0;i<BasisSet.max_am;i++){
	
      norm_ptr = GTOs.bf_norm[i];
      shell_type = i;
      switch(i){
	
	/* ------------------------------------------------------*/
	/* S-type functions */
	
      case 0:
	norm_ptr0 = norm_ptr[0];
	for(j=0;j<close->close_shells_per_am[i];j++){
	  
	  am2shell = close->shells_close_to_chunk[m];
	  shell_center = BasisSet.shells[am2shell].center-1;
	  xa = dist_coord[shell_center].x;
	  ya = dist_coord[shell_center].y;
	  za = dist_coord[shell_center].z;
	  rr = dist_atom[shell_center];
	  
	  shell_start = BasisSet.shells[am2shell].fprim-1;
	  shell_end = shell_start
	    +BasisSet.shells[am2shell].n_prims;
	  
	  
	  timer_on("exponent");
	  bastmp = 0.0;
	  
		/*bastmp = calc_exp_basis(am2shell,rr);*/
	  for(k=shell_start;k<shell_end;k++){
	    expon = -BasisSet.cgtos[k].exp;
	    coeff = BasisSet.cgtos[k].ccoeff[shell_type];
	    bastmp += coeff*exp(expon*rr);
	    
	  }
	  
	  timer_off("exponent");
	  
	  DFT_options.basis[l] = norm_ptr0*bastmp;
	  l++;
	  
	  m++;
	}
	
	break;
	/* ------------------------------------------------------*/
	
	
	/* ------------------------------------------------------*/
	/* P-type functions */
	
      case 1:
	norm_ptr0 = norm_ptr[0];
	norm_ptr1 = norm_ptr[1];
	norm_ptr2 = norm_ptr[2];
	
	for(j=0;j<close->close_shells_per_am[i];j++){
	  
	  am2shell = close->shells_close_to_chunk[m];
	  shell_center = BasisSet.shells[am2shell].center-1;
	  
	  xa = dist_coord[shell_center].x;
	  ya = dist_coord[shell_center].y;
	  za = dist_coord[shell_center].z;
	  rr = dist_atom[shell_center];
	  
	  shell_start = BasisSet.shells[am2shell].fprim-1;
	  shell_end = shell_start
	    +BasisSet.shells[am2shell].n_prims;
	  
	  
	  timer_on("exponent");
	  bastmp = 0.0;
	  for(k=shell_start;k<shell_end;k++){
	    expon = -BasisSet.cgtos[k].exp*rr;
	    coeff = BasisSet.cgtos[k].ccoeff[shell_type];
	    bastmp += coeff*exp(expon);
	    
	  }
	  
	  timer_off("exponent");
	  
	  DFT_options.basis[l] = norm_ptr0*bastmp*xa;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr1*bastmp*ya;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr2*bastmp*za;
	  l++;
	  m++;
	}
	break;
	/* ------------------------------------------------------*/
	
	/* ------------------------------------------------------*/
	/* D-type functions */ 
	
      case 2:
	norm_ptr0 = norm_ptr[0];
	norm_ptr1 = norm_ptr[1];
	norm_ptr2 = norm_ptr[2];
	norm_ptr3 = norm_ptr[3];
	norm_ptr4 = norm_ptr[4];
	norm_ptr5 = norm_ptr[5];
	for(j=0;j<close->close_shells_per_am[i];j++){
	  am2shell = close->shells_close_to_chunk[m];
	  shell_center = BasisSet.shells[am2shell].center-1;
	  xa = dist_coord[shell_center].x;
	  ya = dist_coord[shell_center].y;
	  za = dist_coord[shell_center].z;
	  rr = dist_atom[shell_center];
	  
	  shell_start = BasisSet.shells[am2shell].fprim-1;
	  shell_end = shell_start
	    +BasisSet.shells[am2shell].n_prims;
	  
	  
	  timer_on("exponent");
	  bastmp = 0.0;
	  for(k=shell_start;k<shell_end;k++){
	    expon = -BasisSet.cgtos[k].exp*rr;
	    coeff = BasisSet.cgtos[k].ccoeff[shell_type];
	    bastmp += coeff*exp(expon);
	    
	  }
	  timer_off("exponent");
	  
	  DFT_options.basis[l] = norm_ptr0*bastmp*xa*xa;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr1*bastmp*xa*ya;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr2*bastmp*xa*za;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr3*bastmp*ya*ya;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr4*bastmp*ya*za;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr5*bastmp*za*za;
	  l++;
	  m++;
	}
	break;
	/* ------------------------------------------------------*/
	
	/* ------------------------------------------------------*/
	/* F-type functions */
	
      case 3:
	norm_ptr0 = norm_ptr[0];
	norm_ptr1 = norm_ptr[1];
	norm_ptr2 = norm_ptr[2];
	norm_ptr3 = norm_ptr[3];
	norm_ptr4 = norm_ptr[4];
	norm_ptr5 = norm_ptr[5];
	norm_ptr6 = norm_ptr[6];
	norm_ptr7 = norm_ptr[7];
	norm_ptr8 = norm_ptr[8];
	norm_ptr9 = norm_ptr[9];
	
	for(j=0;j<close->close_shells_per_am[i];j++){
	  am2shell = close->shells_close_to_chunk[m];
	  shell_center = BasisSet.shells[am2shell].center-1;
	  xa = dist_coord[shell_center].x;
	  ya = dist_coord[shell_center].y;
	  za = dist_coord[shell_center].z;
	  rr = dist_atom[shell_center];
	  
	  xx = xa*xa;
	  xy = xa*ya;
	  yy = ya*ya;
	  zz = za*za;
	  
	  bxx = bastmp*xx;
	  byy = bastmp*yy;
	  bxy = bastmp*xy;
	  bzz = bastmp*zz;
	  
	  shell_start = BasisSet.shells[am2shell].fprim-1;
	  shell_end = shell_start
	    +BasisSet.shells[am2shell].n_prims;
	  
	  
	  timer_on("exponent");
	  bastmp = 0.0;
	  for(k=shell_start;k<shell_end;k++){
	    expon = -BasisSet.cgtos[k].exp*rr;
	    coeff = BasisSet.cgtos[k].ccoeff[shell_type];
	    bastmp += coeff*exp(expon);
	  }
	  timer_off("exponent");
	  
	  DFT_options.basis[l] = norm_ptr0*bxx*xa;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr1*bxx*ya;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr2*bxx*za;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr3*byy*xa;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr4*bxy*za;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr5*bzz*xa;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr6*byy*ya;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr7*byy*za;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr8*bzz*ya;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr9*bzz*za;
	  l++;
	  m++;
	}
	
	break;
	
	/* ------------------------------------------------------*/
	
	/* ------------------------------------------------------*/
	/* G-type functions */
      case 4:
	norm_ptr0 = norm_ptr[0];
	norm_ptr1 = norm_ptr[1];
	norm_ptr2 = norm_ptr[2];
	norm_ptr3 = norm_ptr[3];
	norm_ptr4 = norm_ptr[4];
	norm_ptr5 = norm_ptr[5];
	norm_ptr6 = norm_ptr[6];
	norm_ptr7 = norm_ptr[7];
	norm_ptr8 = norm_ptr[8];
	norm_ptr9 = norm_ptr[9];
	norm_ptr10 = norm_ptr[10];
	norm_ptr11 = norm_ptr[11];
	norm_ptr12 = norm_ptr[12];
	norm_ptr13 = norm_ptr[13];
	norm_ptr14 = norm_ptr[14];
	for(j=0;j<close->close_shells_per_am[i];j++){
	  am2shell = close->shells_close_to_chunk[m];
	  shell_center = BasisSet.shells[am2shell].center-1;
	  xa = dist_coord[shell_center].x;
	  ya = dist_coord[shell_center].y;
	  za = dist_coord[shell_center].z;
	  rr = dist_atom[shell_center];
	  
	  shell_start = BasisSet.shells[am2shell].fprim-1;
	  shell_end = shell_start
	    +BasisSet.shells[am2shell].n_prims;
	  
	  
	  timer_on("exponent");
	  bastmp = 0.0;
	  for(k=shell_start;k<shell_end;k++){
	    expon = -BasisSet.cgtos[k].exp;
	    coeff = BasisSet.cgtos[k].ccoeff[shell_type];
	    bastmp += coeff*exp(expon*rr);
	  }
	  timer_off("exponent");
	  DFT_options.basis[l] = norm_ptr0*bastmp*xa*xa*xa*xa;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr1*bastmp*xa*xa*xa*ya;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr2*bastmp*xa*xa*xa*za;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr3*bastmp*xa*xa*ya*ya;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr4*bastmp*xa*xa*ya*za;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr5*bastmp*xa*xa*za*za;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr6*bastmp*xa*ya*ya*ya;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr7*bastmp*xa*ya*ya*za;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr8*bastmp*xa*ya*za*za;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr9*bastmp*xa*za*za*za;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr10*bastmp*ya*ya*ya*ya;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr11*bastmp*ya*ya*ya*za;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr12*bastmp*ya*ya*za*za;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr13*bastmp*ya*za*za*za;
	  l++;
	  
	  DFT_options.basis[l] = norm_ptr14*bastmp*za*za*za*za;
	  l++;
	  m++;
	}
	
	break;
	
      default:
	fprintf(stderr,"Basis Functions of Angular Momentum \n \
%d not implemented in get_density function",shell_type);
	fprintf(outfile,"Basis Functions of Angular Momentum \n \
%d not implemented in get_density function",shell_type);
	throw std::domain_error("");
      }
      
    }
    timer_off("basis"); 
    /* Now contract the basis functions with the AO density matrix elements */
    timer_on("density"); 
    /*for(i=0;i<close->num_close_aos;i++)
      fprintf(outfile,"\nBasis[%d] = %10.10lf",i,DFT_options.basis[i]);*/
    
    den_sum = 0.0; 
#if USE_BLAS
    C_DGEMV('t',close->num_close_aos,ndocc,1.0,close->close_COCC[0],ndocc,
	    DFT_options.basis,1,0.0,temp_arr,1);
    den_sum = C_DDOT(ndocc,temp_arr,1,temp_arr,1);
#else
    for(i=0;i<ndocc;i++){
      for(j=0;j<close->num_close_aos;j++){
	temp_arr[i] += close->close_COCC[j][i]*DFT_options.basis[j];
      }
      /*fprintf(outfile,"\ntemp_arr[%d] = %10.10lf",i,temp_arr[i]);*/
    }
    dot_arr(temp_arr,temp_arr,MOInfo.ndocc,&den_sum);
#endif
    den_info.den = den_sum;
    /*fprintf(outfile,"\nden = %10.10lf",den_info.den);*/
    
    
    free(temp_arr);
    timer_off("density");
    free(dist_coord);
    free(dist_atom);
    return den_info;
  }
}}
