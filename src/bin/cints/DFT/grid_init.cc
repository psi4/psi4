/*! \file grid_init.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include <cstdlib>
#include <cstring>
#include<cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"init_unf_prim_atomic_grid.h"
#include"init_prun_prim_atomic_grid.h"
#include"init_unf_conc_grid.h"
#include"init_prun_conc_grid.h"
#include"free_grid_structs.h"
#include"bragg.h"
#include"SG1.h"
#include"physconst.h"

namespace psi { namespace CINTS {
  void grid_init(){
    int errcod,i,j;
    int depth;
    int num_ang_grids;
    int rpointstmp;
    int angpoints;
    int chunktmp;
    int start,end;
    int chunk_size;
    int u_atom_num;
    int Z_nuc;
    
    char *gridstring;
    
    /* Contcrete Object */
    grid_t *grid;    
    
    /* ---------------------------
       Do all of the input parsing
       and then initialize primitive
       atomic grid
       --------------------------*/
   
    grid = &(DFT_options.grid);
    
    depth = 0;
    
    /* read the number of Abrams chunks */
    chunktmp = 4;
    errcod = ip_data("CHUNK_NUM","%d",&chunktmp,0);
    
    /* The SG-1 grid will be the default for now */
    errcod = ip_count("GRID",&depth,0);
    
    /* If depth is 0 then it is a special gridtype */
    if(depth==0){
	gridstring = "SG1";
	errcod = ip_string("GRID",&gridstring,0);
	
	/* This is the initialization fo the SG1 grid */

	if(!strcmp(gridstring,"SG1") || !strcmp(gridstring,"SG-1")){
	    
	    rpointstmp = 50;
	    
	    grid->pruned_info.n_param_sets = 3;
	    grid->pruned_info.n_diff_ang_grids = 1;
	    grid->pruned_info.n_tot_ang_grids = 5;
	    grid->pruned_info.a2param = 
		(int *)init_array(Symmetry.num_unique_atoms);
	
	    grid->pruned_info.param_set = 
		(struct param_set_s *)malloc(sizeof(struct param_set_s)*
					     grid->pruned_info.n_param_sets);
	
	    for(i=0;i<grid->pruned_info.n_param_sets;i++){
		grid->pruned_info.param_set[i].n_ang_grids = 5;
		
		grid->pruned_info.param_set[i].alpha = 
		    init_array(grid->pruned_info.param_set[i].n_ang_grids-1);
	    
		grid->pruned_info.param_set[i].angpoints = 
		    (int *) init_array(grid->pruned_info.param_set[i].n_ang_grids);
	   
		grid->pruned_info.param_set[i].alpha[0] = SG1alpha1[i];
		grid->pruned_info.param_set[i].alpha[1] = SG1alpha2[i];
		grid->pruned_info.param_set[i].alpha[2] = SG1alpha3[i];
		grid->pruned_info.param_set[i].alpha[3] = SG1alpha4[i];
		for(j=0;j<5;j++){
		    grid->pruned_info.param_set[i].angpoints[j] = SG1angular[j];
		}
	    }
	    
	    for(i=0;i<Symmetry.num_unique_atoms;i++){
		Z_nuc = (int) Molecule.centers[Symmetry.ua2a[i]].Z_nuc;
		grid->pruned_info.a2param[i] = SG1a2param[Z_nuc];
	    }
	    
	    grid->prim_pruned_atomic_grids = 
		init_pruned_prim_atomic_grid(rpointstmp,
					     chunktmp,
					     grid->pruned_info);
	    init_pruned_con_grid();
	    
	
	    grid->label = "SG-1";
	    grid->n_rad_points = 51;
	    
	}
	else
	    throw std::domain_error("No Special Grids have been implemented in this code yet");
	
    }

    /* if the depth is 2 then it is a EML grid that is uniform */

    else if(depth == 2){
	errcod = ip_data("GRID","%d",&rpointstmp,1,0);
	if(errcod != IPE_OK)
	    throw std::domain_error("Problem with Grid specification: the radial points");
	
	errcod = ip_data("GRID","%d",&angpoints,1,1);
	if(errcod != IPE_OK)
	    throw std::domain_error("Problem with Grid specification: the angular points");
       
	/* -------------------------
	   Construct Primitive atomic grid data 
	   ------------------------*/
	
	grid->prim_atomic_grid = 
	    init_uniform_prim_atomic_grid(rpointstmp,angpoints,chunktmp);
	init_uniform_con_grid();
       
	grid->label = "Euler-Mclaren / Lebedev Spheres";
	grid->n_rad_points = rpointstmp+1;
    }
    else
	throw std::domain_error("Problem with Grid specification: Wrong number of elements for keyword Grid");
    
    
    for(i=0;i<Symmetry.num_unique_atoms;i++){
	grid->atomic_grid[i].Bragg_radii = Bragg_radii[(int) Molecule.centers[Symmetry.ua2a[i]].Z_nuc];
    }
   
    return;
  }
}}
