/*! \file init_prun_prim_atomic_grid.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"lebedev_init.h"
#include"physconst.h"

namespace psi {
namespace CINTS {
prim_atomic_grid_t *init_pruned_prim_atomic_grid(int n_rpoints, int num_chunk, struct pruned_info_s pruned_info){
    int i,j,k,l;
    int start,end;
    int chunk_size;
    int cutoff_index;
    int angpoints;
    double qr;
    double r;
    double rind;
    double rext;
    double n_rpoints_d;
    double four_pi_div_by_rps;
    double x,y,z;
    double drdq;
    
    prim_atomic_grid_t *prim_atomic_grid;
    leb_sphere_t *unit_sphere;
    leb_sphere_t *sph;
    
    /* Constants */
    
    n_rpoints_d = (double) n_rpoints+1.0;
    four_pi_div_by_rps = 4.0*_pi/n_rpoints_d;
       
    /*-------------------------
      Initialize the unit sphere, 
      one for each type of angular grid used 
      ---------------------------*/
    
    unit_sphere = (leb_sphere_t *)malloc(sizeof(leb_sphere_t)*pruned_info.n_tot_ang_grids);

    for(i=0;i<pruned_info.n_tot_ang_grids;i++)
	unit_sphere[i] = lebedev_init(pruned_info.param_set[0].angpoints[i]);    
    
    /* ------------------------------
       In case there is more than one
       primitive atomic grid 
       -----------------------------*/
    
    prim_atomic_grid = (prim_atomic_grid_t *)
	malloc(sizeof(prim_atomic_grid_t)*pruned_info.n_param_sets);
    
/*------------------------
  Set up primitive chunks
  -----------------------*/
    
    chunk_size = n_rpoints/num_chunk;
    
    for(i=0;i<pruned_info.n_param_sets;i++){
	cutoff_index = 0;
	prim_atomic_grid[i].chunk_num = num_chunk;
	
	prim_atomic_grid[i].leb_chunk = (prim_leb_chunk_t *)
	    malloc(sizeof(prim_leb_chunk_t)*num_chunk);
	
	for(j=0;j<num_chunk;j++){
	    
	    /* ----- Set up radial offsets for each chunk ------*/
	    
	    start = j*chunk_size+1;
	    end = start+chunk_size;
	    
	    
	    if(j == num_chunk-1){
		end = n_rpoints;
		chunk_size = end-start; 
	    }
	    
	    /*------------------------------
	      Here I am actually going to 
	      calculate the r values as
	      if the Bragg radii was 1.0.
	      This way the primitive 
	      atomic grid will be self contained
	      ------------------------------*/
	    
	    prim_atomic_grid[i].leb_chunk[j].spheres = (leb_sphere_t *)
		malloc(sizeof(leb_sphere_t)*chunk_size);
	    
	    for(k=0;k<chunk_size;k++){
		sph = &(prim_atomic_grid[i].leb_chunk[j].spheres[k]);
		
		rind = (double) k + (double) start;
		
		qr = rind/n_rpoints_d;
		
		/* -------------------------------
		   Straight from the Murray, Handy, Laming paper
		   for mr = 2 
		   ----------------------------------*/
		
		
		r = rind*rind/((n_rpoints_d  - rind)
			       *(n_rpoints_d  - rind));
		
		/*drdq = four_pi_div_by_rps*r*r*2.0*qr/((1-qr)*(1-qr)*(1-qr));*/
		
		/* ---------------------------------
		   If I have moved into another 
		   pruned shell, increment the 
		   parameters
		   ---------------------------------*/
		
		if(r >= pruned_info.param_set[i].alpha[cutoff_index] 
		   && cutoff_index<pruned_info.param_set[i].n_ang_grids-1)
		    cutoff_index++;
		
		angpoints = pruned_info.param_set[i].angpoints[cutoff_index];
		
		
		drdq = 2.0*pow(rind,5)*(n_rpoints_d)*pow(n_rpoints_d-rind,-7.0);
		sph->points = (leb_point_t *)malloc(sizeof(leb_point_t)*angpoints);
		
		for(l=0;l<angpoints;l++){
		    
		    sph->points[l].p_cart.x = 
			unit_sphere[cutoff_index].points[l].p_cart.x*r;
		    sph->points[l].p_cart.y =
			unit_sphere[cutoff_index].points[l].p_cart.y*r;
		    sph->points[l].p_cart.z =
			unit_sphere[cutoff_index].points[l].p_cart.z*r;
		    sph->points[l].ang_weight =
			unit_sphere[cutoff_index].points[l].ang_weight;
		}
		
		sph->r = r;
		sph->drdq = drdq;
		sph->n_ang_points = unit_sphere[cutoff_index].n_ang_points;
	    }
	    
	    prim_atomic_grid[i].leb_chunk[j].radial_start = start;
	    prim_atomic_grid[i].leb_chunk[j].radial_end = end;
	    prim_atomic_grid[i].leb_chunk[j].size = chunk_size;    
	}
    }
    return prim_atomic_grid;
}	      
}
}
