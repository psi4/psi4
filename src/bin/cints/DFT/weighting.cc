/*! \file weighting.cc
    \ingroup CINTS
    \author Shawn Brown
   
    This code contains some functions used for 
    the calculation of the weighting functions
    for the Becke scheme.
    
    Ref.  Becke, J. Chem. Phys., Vol. 88, pg. 2547, 1988.
    
    ----------------------------------------------*/
#include <cmath>
#include <cstdio>
#include <memory.h>
#include <cstdlib>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"small_fns.h"

/* Declare functions in this code */

namespace psi { namespace cints {
double u_calc(int i, int j, struct coordinates geom); 
double v_calc(int i, int j, double uij);
double f_u(double vij);
double s_u(double f);

double weight_calc(int atomn,struct coordinates  geom,int k_order){
    int i,j,k,l;
    int natoms;
    double utemp,vtemp,ftemp;
    double **s_mat;
    double *ptemp;
    double sum=0.0;
    double weight;
    
    natoms = Molecule.num_atoms;
    /*natoms = Symmetry.num_unique_atoms;*/
    
    ptemp = init_array(natoms);
    s_mat = block_matrix(natoms,natoms);
    
    /* calculate all possible s's*/
    for(i=0;i<natoms;i++){
	for(j=0;j<natoms;j++){
	    if(i!=j){
		utemp = u_calc(i,j,geom);
		vtemp = v_calc(i,j,utemp);
		ftemp = f_u(utemp);
		for(k=1;k<k_order;k++){
		    ftemp = f_u(ftemp);
		}
		s_mat[i][j] = s_u(ftemp);
	    }
	}
    }
    	
    /* calculate all of the cell functions into an array */
    
    for(i=0;i<natoms;i++){
	ptemp[i] = 1.0;
	for(j=0;j<natoms;j++){
	    if(i!=j)
		ptemp[i] *= s_mat[i][j];
	}
	sum += ptemp[i];
    }
    
    weight=ptemp[atomn]/sum;
    free(ptemp);
    free_block(s_mat);
    
    return weight;
 
}
		
double u_calc(int i, int j, struct coordinates geom){
    
    double ri;
    double rj;
    double rij;
    /* calculate the distance of the grid point from
       each atom */
    
    ri = distance_calc(Molecule.centers[i],geom);
    rj = distance_calc(Molecule.centers[j],geom);
    rij = distance_calc(Molecule.centers[i],Molecule.centers[j]);
    return (ri-rj)/rij;
}


double f_u(double vij){
    
    return 1.5*vij-0.5*(vij*vij*vij);
}

double s_u(double f){
    
    return 0.5*(1.0-f);
}
    
double v_calc(int atomi, int atomj,double uij){
    
    double aij;
    double utmp;
    double rattmp;
    double tmp1,tmp2,tmp3;

    rattmp = DFT_options.grid.atomic_grid[atomi].Bragg_radii/
	DFT_options.grid.atomic_grid[atomj].Bragg_radii;
    utmp = (rattmp-1)/(rattmp+1);
    aij=utmp/((utmp*utmp)-1);
    if(aij>0.5) aij = 0.5;
    if(aij<-0.5) aij = -0.5;
    tmp1 = uij*uij;
    tmp2 = aij*(1-tmp1);
    tmp3 = uij+tmp2;

    return uij+tmp2;
}               
}}
