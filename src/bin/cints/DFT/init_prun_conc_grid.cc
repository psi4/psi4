/*! \file init_prun_conc_grid.cc
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
#include"physconst.h"

namespace psi { namespace CINTS {
void init_pruned_con_grid(void){
    int i,j;
    int u_atom_num;
        
    struct leb_chunk_s *chnk;
    prim_leb_chunk_t *prim_chnk;
    
    grid_t *grid;
    
    grid = &(DFT_options.grid);
    
    grid->pruned_flag = 1;
    
    grid->atomic_grid = (struct atomic_grid_s *)
	malloc(sizeof(struct atomic_grid_s)*Symmetry.num_unique_atoms);
    
    for(i=0;i<Symmetry.num_unique_atoms;i++){
	
	/* ------- Must take care of symmetry here ------- */
	
	u_atom_num = Symmetry.ua2a[i];
	grid->atomic_grid[i].atom_num = u_atom_num;
	grid->atomic_grid[i].atom_center = Molecule.centers[u_atom_num];
	
	/* --- Calculate the degeneracy of symmetry unique atoms ---*/
	
	if(Symmetry.nirreps > 1)
	    grid->atomic_grid[i].atom_degen = Symmetry.nirreps/
		Symmetry.dcr_deg[Symmetry.atom_positions[u_atom_num]]
		[Symmetry.atom_positions[u_atom_num]];
	else
	    grid->atomic_grid[i].atom_degen = 1;
	
	/* --- These may not exactly be the Bragg-Slater Radii ---*/
	/* ------------ See Bragg.h for the reference ------------*/
	/*
	grid->atomic_grid[i].Bragg_radii = 
	    Bragg_radii[(int) Molecule.centers[Symmetry.ua2a[i]].Z_nuc];
	    */
	/* ---- Set up chunk information ----*/
	/* ---- This is the only part that depends on whether you it is a pruned grid or not */
	
	grid->atomic_grid[i].chunk_num = 
	    grid->prim_pruned_atomic_grids[grid->pruned_info.a2param[i]].chunk_num;
	
	grid->atomic_grid[i].leb_chunk = (struct leb_chunk_s *)
	    malloc(sizeof(struct leb_chunk_s)*grid->atomic_grid[i].chunk_num);
	
	for(j=0;j<grid->atomic_grid[i].chunk_num;j++){
	    chnk = &(grid->atomic_grid[i].leb_chunk[j]);
	    prim_chnk = 
		&(grid->prim_pruned_atomic_grids[grid->pruned_info.a2param[i]].leb_chunk[j]);
	    
	    chnk->radial_start = prim_chnk->radial_start;
	    chnk->radial_end = prim_chnk->radial_end;
	    chnk->size = prim_chnk->size;
	    chnk->spheres = prim_chnk->spheres;
	    
	    chnk->shells_close_to_chunk = 
		(int *)malloc(sizeof(int)*BasisSet.num_ao);
	    chnk->close_shells_per_am = 
		(int *)malloc(sizeof(int)*BasisSet.max_am);
	}
    }
    return;
}   
};};
