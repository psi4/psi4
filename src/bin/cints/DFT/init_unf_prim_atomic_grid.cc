/*! \file init_unf_prim_atomic_grid.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"lebedev_init.h"
#include"physconst.h"

namespace psi { namespace CINTS {

prim_atomic_grid_t init_uniform_prim_atomic_grid(int n_rpoints,int n_angpoints,int num_chunks){
    int i,j,k;
    int start,end;
    int chunk_size;
    int n_rpoints_plus_two;
    double qr;
    double r;
    double rind;
    double n_rpoints_d;
    double four_pi_div_by_rps;
    double x,y,z;
    double drdq;
    
    prim_atomic_grid_t prim_atomic_grid;
    leb_sphere_t unit_sphere;
    leb_sphere_t *sph;
    
    /* Constants */
    
    n_rpoints_plus_two = n_rpoints+1.0;
    n_rpoints_d = (double) n_rpoints+1.0;
    four_pi_div_by_rps = 4.0*_pi/n_rpoints_d;
    
    prim_atomic_grid.chunk_num = num_chunks;

/*-------------------------
  Initialize the unit sphere, 
  there is only one here 
  ---------------------------*/
    
    unit_sphere = lebedev_init(n_angpoints);
    
    /* ------------------------
       Set up primitive chunks
       -----------------------*/
    
    chunk_size = n_rpoints/num_chunks;
    
    prim_atomic_grid.leb_chunk = (prim_leb_chunk_t *)
	malloc(sizeof(prim_leb_chunk_t)*num_chunks);
    
    for(i=0;i<num_chunks;i++){
	
	/* ----- Set up radial offsets for each chunk ------*/
	
	start = i*chunk_size+1;
	end = start+chunk_size;
	
	
	if(i == num_chunks-1){
	    end = n_rpoints+1;
	    chunk_size = end-start; 
	}
	
	
	
	/*------------------------------
	  Here I am actually going to 
	  calculate the r values as
	  if the Bragg radii was 1.0.
	  This way the primitive 
	  atomic grid will be self contained
	  ------------------------------*/
	
	prim_atomic_grid.leb_chunk[i].spheres = (leb_sphere_t *)
	    malloc(sizeof(leb_sphere_t)*chunk_size);
	
	for(j=0;j<chunk_size;j++){
	    sph = &(prim_atomic_grid.leb_chunk[i].spheres[j]);
	    
	    rind = (double) j + (double) start;
	    
	    qr = rind/n_rpoints_d;
	    
	    /* -------------------------------
	       Straight from the Murray, Handy, Laming paper
	       for mr = 2 
	       ----------------------------------*/
		
	    r = rind*rind/((n_rpoints_d  - rind)
			   *(n_rpoints_d  - rind));
	    
	    /*drdq = four_pi_div_by_rps*r*r*2.0*qr/((1-qr)*(1-qr)*(1-qr));*/
	    drdq = 2.0*pow(rind,5)*(n_rpoints_d)*pow(n_rpoints_d-rind,-7.0);
	    sph->points = (leb_point_t *)malloc(sizeof(leb_point_t)*n_angpoints);
	    
	    for(k=0;k<n_angpoints;k++){
		
		sph->points[k].p_cart.x = 
		    unit_sphere.points[k].p_cart.x*r;
		sph->points[k].p_cart.y =
		    unit_sphere.points[k].p_cart.y*r;
		sph->points[k].p_cart.z =
		    unit_sphere.points[k].p_cart.z*r;
		sph->points[k].ang_weight =
		    unit_sphere.points[k].ang_weight;
	    }

	    sph->r = r;
	    
	    sph->drdq = drdq;
	    sph->n_ang_points = unit_sphere.n_ang_points;
	}
	
	prim_atomic_grid.leb_chunk[i].radial_start = start;
	prim_atomic_grid.leb_chunk[i].radial_end = end;
	prim_atomic_grid.leb_chunk[i].size = chunk_size;    
    }
    return prim_atomic_grid;
}	   
};}
