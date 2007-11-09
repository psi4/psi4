/*! \file calc_den.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
    \author  Shawn Brown
    
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
  
  struct den_info_s calc_density(struct coordinates geom){
    
    int i,j,k;
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
    double bastmp1;
    double den_sum;
    double coeff;
    double expon;
    double *norm_ptr;
    double *dist_atom;
    double *temp_arr;
    
    struct coordinates *dist_coord;
    struct den_info_s den_info;
    struct shell_pair *sp;
    
    x = geom.x;
    y = geom.y;
    z = geom.z;
    
    num_ao = BasisSet.num_ao;
    ndocc = MOInfo.ndocc;
    temp_arr = init_array(num_ao);
    dist_atom = init_array(Molecule.num_atoms);
    dist_coord = (struct coordinates *)malloc(sizeof(struct coordinates)*Molecule.num_atoms);
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
    for(i=k=0;i<n_shells;i++){
      
      shell_type = BasisSet.shells[i].am;
      shell_center = BasisSet.shells[i].center-1;
      xa = dist_coord[shell_center].x;
      ya = dist_coord[shell_center].y;
      za = dist_coord[shell_center].z;
      rr = dist_atom[shell_center];
      
      
      shell_start = BasisSet.shells[i].fprim-1;
      shell_end = shell_start+BasisSet.shells[i].n_prims;
      
      norm_ptr = GTOs.bf_norm[shell_type-1];
      timer_on("exponent");
      /*bastmp1 = calc_exp_basis(i,rr);*/
      
      bastmp = 0;
      for(j=shell_start;j<shell_end;j++){
	expon = -BasisSet.cgtos[j].exp;
	coeff = BasisSet.cgtos[j].ccoeff[shell_type-1];
	bastmp += coeff*exp(expon*rr);
      }
      /*if(bastmp != bastmp1)
	fprintf(outfile,"\nbastmp1 = %10.15lf bastmp2 = %10.15lf for shell %d at rr = %e",bastmp1,bastmp,i,rr);*/
      timer_off("exponent");	
      /*----------------------------------
	Compute values of basis functions
	
	NOTE: using Psi 3 ordering of
	functions within shells
	----------------------------------*/
      switch (shell_type) {
	
      case 1:
	DFT_options.basis[k] = norm_ptr[0]*bastmp;
	k++;
	break;
	
      case 2:
	DFT_options.basis[k] = norm_ptr[0]*bastmp*xa;
	k++;
	
	DFT_options.basis[k] = norm_ptr[1]*bastmp*ya;
	k++;
	
	DFT_options.basis[k] = norm_ptr[2]*bastmp*za;
	k++;
	break;
      case 3:
	DFT_options.basis[k] = norm_ptr[0]*bastmp*xa*xa;
	k++;
	
	DFT_options.basis[k] = norm_ptr[1]*bastmp*xa*ya;
	k++;
	
	DFT_options.basis[k] = norm_ptr[2]*bastmp*xa*za;
	k++;
	
	DFT_options.basis[k] = norm_ptr[3]*bastmp*ya*ya;
	k++;
	
	DFT_options.basis[k] = norm_ptr[4]*bastmp*ya*za;
	k++;
	
	DFT_options.basis[k] = norm_ptr[5]*bastmp*za*za;
	k++;
	break;
      case 4:
	DFT_options.basis[k] = norm_ptr[0]*bastmp*xa*xa*xa;
	k++;
	
	DFT_options.basis[k] = norm_ptr[1]*bastmp*xa*xa*ya;
	k++;
	
	DFT_options.basis[k] = norm_ptr[2]*bastmp*xa*xa*za;
	k++;
	
	DFT_options.basis[k] = norm_ptr[3]*bastmp*xa*ya*ya;
	k++;
	
	DFT_options.basis[k] = norm_ptr[4]*bastmp*xa*ya*za;
	k++;
	
	DFT_options.basis[k] = norm_ptr[5]*bastmp*xa*za*za;
	k++;
	
	DFT_options.basis[k] = norm_ptr[6]*bastmp*ya*ya*ya;
	k++;
	
	DFT_options.basis[k] = norm_ptr[7]*bastmp*ya*ya*za;
	k++;
	
	DFT_options.basis[k] = norm_ptr[8]*bastmp*ya*za*za;
	k++;
	
	DFT_options.basis[k] = norm_ptr[9]*bastmp*za*za*za;
	k++;
	break;
      case 5:
	DFT_options.basis[k] = norm_ptr[0]*bastmp*xa*xa*xa*xa;
	k++;
	
	DFT_options.basis[k] = norm_ptr[1]*bastmp*xa*xa*xa*ya;
	k++;
	    
	DFT_options.basis[k] = norm_ptr[2]*bastmp*xa*xa*xa*za;
	k++;
	
	DFT_options.basis[k] = norm_ptr[3]*bastmp*xa*xa*ya*ya;
	k++;
	
	DFT_options.basis[k] = norm_ptr[4]*bastmp*xa*xa*ya*za;
	k++;
	
	DFT_options.basis[k] = norm_ptr[5]*bastmp*xa*xa*za*za;
	k++;
	
	DFT_options.basis[k] = norm_ptr[6]*bastmp*xa*ya*ya*ya;
	k++;
	    
	DFT_options.basis[k] = norm_ptr[7]*bastmp*xa*ya*ya*za;
	k++;
	
	DFT_options.basis[k] = norm_ptr[8]*bastmp*xa*ya*za*za;
	k++;
	
	DFT_options.basis[k] = norm_ptr[9]*bastmp*xa*za*za*za;
	k++;
	
	DFT_options.basis[k] = norm_ptr[10]*bastmp*ya*ya*ya*ya;
	k++;
	
	DFT_options.basis[k] = norm_ptr[11]*bastmp*ya*ya*ya*za;
	k++;
	
	DFT_options.basis[k] = norm_ptr[12]*bastmp*ya*ya*za*za;
	k++;
	
	DFT_options.basis[k] = norm_ptr[13]*bastmp*ya*za*za*za;
	k++;
	
	DFT_options.basis[k] = norm_ptr[14]*bastmp*za*za*za*za;
	k++;
	break;
      default:
	fprintf(stderr,"Basis Functions of Angular Momentum %d not implemented in get_density function",shell_type);
	fprintf(outfile,"Basis Functions of Angular Momentum %d not implemented in get_density function",shell_type);
	throw std::domain_error("");
      }
    }
    
    timer_off("basis");
    /* Now contract the basis functions with the AO density matrix elements */
    timer_on("density"); 
    
    if(UserOptions.reftype == rhf){
      den_sum = 0.0;
#if USE_BLAS
      C_DGEMV('t',num_ao,ndocc,1.0,Cocc[0],ndocc,
	      DFT_options.basis,1,0.0,temp_arr,1);
      for(i=0;i<ndocc;i++)
	den_sum = C_DDOT(ndocc,temp_arr,1,temp_arr,1);
#else
      for(i=0;i<ndocc;i++){
	for(j=0;j<num_ao;j++){
	  temp_arr[i] += Cocc[j][i]*DFT_options.basis[j];
	}
	
      }
      dot_arr(temp_arr,temp_arr,MOInfo.ndocc,&den_sum);
#endif
      den_info.den = den_sum;
      
    }
    free(temp_arr);
    timer_off("density");
    free(dist_coord);
    free(dist_atom);
    return den_info;
  }
};
  
};	
